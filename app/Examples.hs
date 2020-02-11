{-# LANGUAGE CPP                                      #-}
{-# LANGUAGE DataKinds                                #-}
{-# LANGUAGE DeriveFoldable                           #-}
{-# LANGUAGE DeriveFunctor                            #-}
{-# LANGUAGE FlexibleContexts                         #-}
{-# LANGUAGE FlexibleInstances                        #-}
{-# LANGUAGE GADTs                                    #-}
{-# LANGUAGE LambdaCase                               #-}
{-# LANGUAGE OverloadedStrings                        #-}
{-# LANGUAGE PatternSynonyms                          #-}
{-# LANGUAGE RecordWildCards                          #-}
{-# LANGUAGE ScopedTypeVariables                      #-}
{-# LANGUAGE StandaloneDeriving                       #-}
{-# LANGUAGE TupleSections                            #-}
{-# LANGUAGE TypeApplications                         #-}
{-# LANGUAGE TypeOperators                            #-}
{-# LANGUAGE TypeSynonymInstances                     #-}
{-# LANGUAGE ViewPatterns                             #-}
{-# OPTIONS_GHC -fno-warn-orphans                     #-}
{-# OPTIONS_GHC -fplugin GHC.TypeLits.KnownNat.Solver #-}

-- | Hamilton example suite
--
-- See: https://github.com/mstksg/hamilton#example-app-runner
--
-- Or just run with:
--
-- > $ hamilton-examples --help
-- > $ hamilton-examples [EXAMPLE] --help
--

import           Control.Concurrent
import           Control.Monad
import           Data.Bifunctor
import           Data.Finite
import           Data.Foldable
import           Data.IORef
import           Data.List
import           Data.Semigroup                      ((<>))
import           GHC.TypeLits
import           Graphics.Vty hiding                 (Config, (<|>))
import           Numeric.Hamilton
import           Numeric.LinearAlgebra.Static hiding (dim, (<>))
import           Numeric.LinearAlgebra.Static.Vector
import           Options.Applicative
import           System.Exit
import           Text.Printf
import           Text.Read
import qualified Data.List.NonEmpty                  as NE
import qualified Data.Map.Strict                     as M
import qualified Data.Vector                         as VV
import qualified Data.Vector.Sized                   as V
import qualified Data.Vector.Storable.Sized          as VS
import qualified Text.PrettyPrint.ANSI.Leijen        as PP

data SysExample where
    SE :: (KnownNat m, KnownNat n)
       => { seName   :: String
          , seCoords :: V.Vector n String
          , seSystem :: System m n
          , seDraw   :: R m -> [V2 Double]
          , seInit   :: Phase n
          }
       -> SysExample

pendulum :: Double -> Double -> SysExample
pendulum θ0 ω0 = SE "Single pendulum" (V1 "θ") s f (toPhase s c0)
  where
    s :: System 2 1
    s = mkSystem' (vec2 1 1                             )     -- masses
                  (\(V1 θ)   -> V2 (sin θ) (0.5 - cos θ))     -- coordinates
                  (\(V2 _ y) -> y                       )     -- potential
    f :: R 2 -> [V2 Double]
    f xs = [grVec xs]
    c0 :: Config 1
    c0 = Cfg (konst θ0 :: R 1) (konst ω0 :: R 1)

doublePendulum :: Double -> Double -> SysExample
doublePendulum m1 m2 = SE "Double pendulum" (V2 "θ1" "θ2") s f (toPhase s c0)
  where
    s :: System 4 2
    s = mkSystem' (vec4 m1 m1 m2 m2)     -- masses
                  (\(V2 θ1 θ2)     -> V4 (sin θ1)            (1 - cos θ1)
                                         (sin θ1 + sin θ2/2) (1 - cos θ1 - cos θ2/2)
                  )                      -- coordinates
                  (\(V4 _ y1 _ y2) -> 5 * (realToFrac m1 * y1 + realToFrac m2 * y2))
                                         -- potential
    f :: R 4 -> [V2 Double]
    f (split->(xs,ys))= grVec <$> [xs, ys]
    c0 :: Config 2
    c0 = Cfg (vec2 (pi/2) 0) (vec2 0 0)

room :: Double -> SysExample
room θ = SE "Room" (V2 "x" "y") s f (toPhase s c0)
  where
    s :: System 2 2
    s = mkSystem (vec2 1 1)         -- masses
                 id                 -- coordinates
                 (\(V2 x y) -> sum [ 2 * y                      -- gravity
                                   , 1 - logistic (-1) 10 0.1 y  -- bottom wall
                                   ,     logistic 1 10 0.1 y     -- top wall
                                   , 1 - logistic (-2) 10 0.1 x  -- left wall
                                   ,     logistic 2 10 0.1 x     -- right wall
                                   ]
                 )                  -- potential
    f :: R 2 -> [V2 Double]
    f xs = [grVec xs]
    c0 :: Config 2
    c0 = Cfg (vec2 (-1) 0.25) (vec2 (cos θ) (sin θ))

twoBody :: Double -> Double -> Double -> SysExample
twoBody m1 m2 ω0 = SE "Two-Body" (V2 "r" "θ") s f (toPhase s c0)
  where
    mT :: Double
    mT = m1 + m2
    s :: System 4 2
    s = mkSystem (vec4 m1 m1 m2 m2) -- masses
                 -- positions are calculated assuming (0,0) is the center
                 -- of mass
                 (\(V2 r θ) -> let r1 = r * realToFrac (-m2 / mT)
                                   r2 = r * realToFrac (m1 / mT)
                               in  V4 (r1 * cos θ) (r1 * sin θ)
                                      (r2 * cos θ) (r2 * sin θ)
                 )                 -- coordinates
                 (\(V2 r _) -> - realToFrac (m1 * m2) / r)  -- potential
    f :: R 4 -> [V2 Double]
    f (split->(xs,ys))= grVec <$> [xs, ys]
    c0 :: Config 2
    c0 = Cfg (vec2 2 0) (vec2 0 ω0)

spring
    :: Double -> Double -> Double -> Double -> SysExample
spring mB mW k x0 = SE "Spring hanging from block" (V3 "r" "x" "θ") s f (toPhase s c0)
  where
    s :: System 3 3
    s = mkSystem (vec3 mB mW mW)                                                  -- masses
                 (\(V3 r x θ)  -> V3 r (r + (1 + x) * sin θ) ((1 + x) * (-cos θ))) -- coordinates
                 (\(V3 r x θ) -> realToFrac k * x**2 / 2        -- spring
                              + (1 - logistic (-1.5) 25 0.1 r)  -- left rail wall
                              + (    logistic   1.5  25 0.1 r)  -- right rail wall
                              + realToFrac mB * ((1 + x) * (-cos θ))  -- gravity
                 )
    f :: R 3 -> [V2 Double]
    f (headTail->(b,w)) = [V2 b 1, V2 0 1 + grVec w]
    c0 :: Config 3
    c0 = Cfg (vec3 0 x0 0) (vec3 1 0 (-0.5))

bezier
    :: forall n. KnownNat (1 + n)
    => V.Vector (1 + n) (V2 Double)
    -> SysExample
bezier ps = SE "Bezier" (V1 "t") s f (toPhase s c0)
  where
    s :: System 2 1
    s = mkSystem (vec2 1 1)                                             -- masses
                 (\(V1 t) -> bezierCurve (fmap realToFrac <$> ps) t)    -- coordinates
                 (\(V1 t) -> (1 - logistic 0 5 0.05 t)           -- left wall
                           +      logistic 1 5 0.05 t            -- right wall
                 )
    f :: R 2 -> [V2 Double]
    f xs = [grVec xs]
    c0 :: Config 1
    c0 = Cfg (0.5 :: R 1) (0.25 :: R 1)


data ExampleOpts = EO { eoChoice :: SysExampleChoice }

data SysExampleChoice =
        SECDoublePend Double Double
      | SECPend Double Double
      | SECRoom Double
      | SECTwoBody Double Double Double
      | SECSpring Double Double Double Double
      | SECBezier (NE.NonEmpty (V2 Double))

parseEO :: Parser ExampleOpts
parseEO = EO <$> (parseSEC <|> pure (SECDoublePend 1 1))

parseSEC :: Parser SysExampleChoice
parseSEC = subparser . mconcat $
    [ command "doublepend" $
        info (helper <*> parseDoublePend)
             (progDesc "Double pendulum (default)")
    , command "pend"       $
        info (helper <*> parsePend      )
             (progDesc "Single pendulum")
    , command "room"       $
        info (helper <*> parseRoom      )
        (progDesc "Ball in room, bouncing off of walls")
    , command "twobody"    $
        info (helper <*> parseTwoBody    )
        (progDesc "Two-body graviational simulation.  Note that bodies will only orbit if H < 0.")
    , command "spring"    $
        info (helper <*> parseSpring    )
        (progDesc "A spring hanging from a block on a rail, holding up a mass.  Block is constrained to bounce between -1.5 and 1.5.")
    , command "bezier"     $
        info (helper <*> parseBezier    )
        (progDesc "Particle moving along a parameterized bezier curve")
    , metavar "EXAMPLE"
    ]
  where
    parsePend
      = SECPend       <$> option auto ( long "angle"
                                     <> short 'a'
                                     <> metavar "ANGLE"
                                     <> help "Intitial rightward angle (in degrees) of bob"
                                     <> value 0
                                     <> showDefault
                                      )
                      <*> option auto ( long "vel"
                                     <> short 'v'
                                     <> metavar "VELOCITY"
                                     <> help "Initial rightward angular velocity of bob"
                                     <> value 1
                                     <> showDefault
                                      )
    parseDoublePend
      = SECDoublePend <$> option auto ( long "m1"
                                     <> metavar "MASS"
                                     <> help "Mass of first bob"
                                     <> value 1
                                     <> showDefault
                                      )
                      <*> option auto ( long "m2"
                                     <> metavar "MASS"
                                     <> help "Mass of second bob"
                                     <> value 1
                                     <> showDefault
                                      )
    parseRoom
      = SECRoom    <$> option auto ( long "angle"
                                  <> short 'a'
                                  <> metavar "ANGLE"
                                  <> help "Initial upward launch angle (in degrees) of object"
                                  <> value 45
                                  <> showDefault
                                   )
    parseTwoBody
      = SECTwoBody <$> option auto ( long "m1"
                                  <> metavar "MASS"
                                  <> help "Mass of first body"
                                  <> value 5
                                  <> showDefault
                                   )
                   <*> option auto ( long "m2"
                                  <> metavar "MASS"
                                  <> help "Mass of second body"
                                  <> value 0.5
                                  <> showDefault
                                   )
                   <*> option auto ( long "vel"
                                  <> short 'v'
                                  <> metavar "VELOCITY"
                                  <> help "Initial angular velocity of system"
                                  <> value 0.5
                                  <> showDefault
                                   )
    parseSpring
      = SECSpring <$> option auto ( long "block"
                                 <> short 'b'
                                 <> metavar "MASS"
                                 <> help "Mass of block on rail"
                                 <> value 2
                                 <> showDefault
                                  )
                  <*> option auto ( long "weight"
                                 <> short 'w'
                                 <> metavar "MASS"
                                 <> help "Mass of weight hanging from spring"
                                 <> value 1
                                 <> showDefault
                                  )
                  <*> option auto ( short 'k'
                                 <> metavar "NUM"
                                 <> help "Spring constant / stiffness of spring"
                                 <> value 10
                                 <> showDefault
                                  )
                  <*> option auto ( short 'x'
                                 <> metavar "DIST"
                                 <> help "Initial displacement of spring"
                                 <> value 0.1
                                 <> showDefault
                                  )
    parseBezier
      = SECBezier <$> option f ( long "points"
                              <> short 'p'
                              <> metavar "POINTS"
                              <> help "List of control points (at least one), as tuples"
                              <> value (V2 (-1) (-1) NE.:| [V2 (-2) 1, V2 0 1, V2 1 (-1), V2 2 1])
                              <> showDefaultWith (show . map (\(V2 x y) -> (x, y)) . toList)
                               )
      where f = eitherReader $ \s -> do
              ps  <- maybe (Left "Bad parse") Right
                  $ readMaybe s
              maybe (Left "At least one control point required") Right
                  $ NE.nonEmpty (uncurry V2 <$> ps)

data SimOpts = SO { soZoom :: Double
                  , soRate :: Double
                  , soHist :: Int
                  }
  deriving (Show)

data SimEvt = SEQuit
            | SEZoom Double
            | SERate Double
            | SEHist Int

main :: IO ()
main = do
    EO{..} <- execParser $ info (helper <*> parseEO)
        ( fullDesc
       <> header "hamilton-examples - hamilton library example suite"
       <> progDescDoc (Just descr)
        )

    vty <- mkVty =<< standardIOConfig

    opts <- newIORef $ SO 0.5 1 25

    t <- forkIO . loop vty opts $ case eoChoice of
      SECDoublePend m1 m2        -> doublePendulum m1 m2
      SECPend       d0 ω0        -> pendulum (d0 / 180 * pi) ω0
      SECRoom       d0           -> room (d0 / 180 * pi)
      SECTwoBody    m1 m2 ω0     -> twoBody m1 m2 ω0
      SECSpring     mB mW k x0   -> spring mB mW k x0
      SECBezier     (p NE.:| ps) -> V.withSized (VV.fromList ps)
                                      (bezier . V.cons p)


    forever $ do
      e <- nextEvent vty
      forM_ (processEvt e) $ \case
        SEQuit -> do
          killThread t
          shutdown vty
          exitSuccess
        SEZoom s ->
          modifyIORef opts $ \o -> o { soZoom = soZoom o * s }
        SERate r ->
          modifyIORef opts $ \o -> o { soRate = soRate o * r }
        SEHist h ->
          modifyIORef opts $ \o -> o { soHist = soHist o + h }
  where
    fps :: Double
    fps = 12
    screenRatio :: Double
    screenRatio = 2.1
    ptAttrs :: [(Char, Color)]
    ptAttrs  = ptChars `zip` ptColors
      where
        ptColors = cycle [white,yellow,blue,red,green]
        ptChars  = cycle "o*+~"
    loop :: Vty -> IORef SimOpts -> SysExample -> IO ()
    loop vty oRef SE{..} = go M.empty seInit
      where
        qVec = intercalate "," . V.toList $ seCoords
        go hists p = do
          SO{..} <- readIORef oRef
          let p'   = stepHam (soRate / fps) seSystem p  -- progress the simulation
              xb   = (- recip soZoom, recip soZoom)
              infobox = vertCat . map (string defAttr) $
                  [ printf "[ %s ]" seName
                  , printf " <%s>   : <%s>" qVec . intercalate ", "
                     . map (printf "%.4f") . VS.toList . rVec . phsPositions $ p
                  , printf "d<%s>/dt: <%s>" qVec . intercalate ", "
                     . map (printf "%.4f") . VS.toList . rVec . velocities seSystem $ p
                  , printf "KE: %.4f" . keP seSystem           $ p
                  , printf "PE: %.4f" . pe seSystem . phsPositions $ p
                  , printf "H : %.4f" . hamiltonian seSystem   $ p
                  , " "
                  , printf "rate: x%.2f <>" $ soRate
                  , printf "hist: % 5d []" $ soHist
                  , printf "zoom: x%.2f -+" $ soZoom
                  ]
              pts  = (`zip` ptAttrs) . seDraw . underlyingPos seSystem . phsPositions
                   $ p
              hists' = foldl' (\h (r, a) -> M.insertWith (addHist soHist) a [r] h) hists pts
          dr <- displayBounds $ outputIface vty
          update vty . picForLayers . (infobox:) . plot dr (PX xb (RR 0.5 screenRatio)) $
               ((second . second) (defAttr `withForeColor`) <$> pts)
            ++ (map (\((_,c),r) -> (r, ('.', defAttr `withForeColor` c)))
                  . concatMap sequence
                  . M.toList
                  $ hists'
               )
          threadDelay (round (1000000 / fps))
          go hists' p'
    addHist hl new old = take hl (new ++ old)
    descr :: PP.Doc
    descr = PP.vcat
      [ "Run examples from the hamilton library example suite."
      , "Use with [EXAMPLE] --help for more per-example options."
      , ""
      , "To adjust rate/history/zoom, use keys <>/[]/-+, respectively."
      , ""
      , "See: https://github.com/mstksg/hamilton#example-app-runner"
      ]

processEvt
    :: Event -> Maybe SimEvt
processEvt = \case
    EvKey KEsc        []      -> Just SEQuit
    EvKey (KChar 'c') [MCtrl] -> Just SEQuit
    EvKey (KChar 'q') []      -> Just SEQuit
    EvKey (KChar '+') []      -> Just $ SEZoom (sqrt 2)
    EvKey (KChar '-') []      -> Just $ SEZoom (sqrt 0.5)
    EvKey (KChar '>') []      -> Just $ SERate (sqrt 2)
    EvKey (KChar '<') []      -> Just $ SERate (sqrt (1/2))
    EvKey (KChar ']') []      -> Just $ SEHist 5
    EvKey (KChar '[') []      -> Just $ SEHist (-5)
    _                         -> Nothing

data RangeRatio = RR { -- | Where on the screen (0 to 1) to place the other axis
                       rrZero  :: Double
                       -- | Ratio of height of a terminal character to width
                     , rrRatio :: Double
                     }
                deriving (Show)

data PlotRange = PXY (Double, Double) (Double, Double)
               | PX  (Double, Double) RangeRatio
               | PY  RangeRatio       (Double, Double)

plot
    :: (Int, Int)               -- ^ display bounds
    -> PlotRange
    -> [(V2 Double, (Char, Attr))]   -- ^ points to plot
    -> [Image]
plot (wd,ht) pr = map (crop wd ht)
                . (++ bgs)
                . map (\(p, (c, a)) -> place EQ EQ p $ char a c)
  where
    wd' = fromIntegral wd
    ht' = fromIntegral ht
    ((xmin, xmax), (ymin, ymax)) = mkRange (wd', ht') pr
    origin = place EQ EQ (V2 0 0) $ char defAttr '+'
    xaxis  = place EQ EQ (V2 0 0) $ charFill defAttr '-' wd 1
    yaxis  = place EQ EQ (V2 0 0) $ charFill defAttr '|' 1 ht
    xrange = xmax - xmin
    yrange = ymax - ymin
    bg     = backgroundFill wd ht
    scale (V2 pX pY) = V2 x y
      where
        x = round $ (pX - xmin) * (wd' / xrange)
        y = round $ (pY - ymin) * (ht' / yrange)
    place aX aY p i = case scale p of
      V2 pX pY -> translate (fAlign aX (imageWidth  i))
                            (fAlign aY (imageHeight i))
                . translate pX pY
                $ i
    labels = [ place LT EQ (V2 xmin 0) . string defAttr $ printf "%.2f" xmin
             , place GT EQ (V2 xmax 0) . string defAttr $ printf "%.2f" xmax
             , place EQ LT (V2 0 ymin) . string defAttr $ printf "%.2f" ymin
             , place EQ GT (V2 0 ymax) . string defAttr $ printf "%.2f" ymax
             ]
    bgs    = labels ++ [origin, xaxis, yaxis, bg]
    fAlign = \case
      LT -> const 0
      EQ -> negate . (`div` 2)
      GT -> negate

mkRange
    :: (Double, Double)
    -> PlotRange
    -> ((Double, Double), (Double, Double))
mkRange (wd, ht) = \case
    PXY xb     yb     -> (xb, yb)
    PX  xb     RR{..} ->
      let yr = (uncurry (-) xb) * ht / wd * rrRatio
          y0 = (rrZero - 1) * yr
      in  (xb, (y0, y0 + yr))
    PY  RR{..} yb ->
      let xr = (uncurry (-) yb) * wd / ht / rrRatio
          x0 = (rrZero - 1) * xr
      in  ((x0, x0 + xr), yb)

pattern V1 :: a -> V.Vector 1 a
pattern V1 x <- (V.head->x)
  where
    V1 x = V.singleton x
#if __GLASGOW_HASKELL__ >= 802
{-# COMPLETE V1 #-}
#endif

type V2 = V.Vector 2
pattern V2 :: a -> a -> V2 a
pattern V2 x y <- (V.toList->[x,y])
  where
    V2 x y = V.fromTuple (x, y)
#if __GLASGOW_HASKELL__ >= 802
{-# COMPLETE V2 #-}
#endif

pattern V3 :: a -> a -> a -> V.Vector 3 a
pattern V3 x y z <- (V.toList->[x,y,z])
  where
    V3 x y z = V.fromTuple (x, y, z)
#if __GLASGOW_HASKELL__ >= 802
{-# COMPLETE V3 #-}
#endif

pattern V4 :: a -> a -> a -> a -> V.Vector 4 a
pattern V4 x y z a <- (V.toList->[x,y,z,a])
  where
    V4 x y z a = V.fromTuple (x, y, z, a)
#if __GLASGOW_HASKELL__ >= 802
{-# COMPLETE V4 #-}
#endif

logistic
    :: Floating a => a -> a -> a -> a -> a
logistic pos ht width = \x -> ht / (1 + exp (- beta * (x - pos)))
  where
    beta = log (0.9 / (1 - 0.9)) / width


bezierCurve
    :: forall n f a. (KnownNat (1 + n), Applicative f, Num a)
    => V.Vector (1 + n) (f a)
    -> a
    -> f a
bezierCurve ps t =
      foldl' (liftA2 (+)) (pure 0)
    . V.imap (\i -> let i' = fromIntegral i
                    in  fmap (* (fromIntegral (n' `choose` i') * (1 - t)^(n' - i') * t^i))
             )
    $ ps
  where
    n' :: Int
    n' = fromIntegral (maxBound :: Finite (1 + n))
    choose :: Int -> Int -> Int
    n `choose` k = factorial n `div` (factorial (n - k) * factorial k)
    factorial :: Int -> Int
    factorial m = product [1..m]

deriving instance Ord Color

