{-# LANGUAGE DataKinds            #-}
{-# LANGUAGE DeriveFoldable       #-}
{-# LANGUAGE DeriveFunctor        #-}
{-# LANGUAGE GADTs                #-}
{-# LANGUAGE LambdaCase           #-}
{-# LANGUAGE PatternSynonyms      #-}
{-# LANGUAGE RecordWildCards      #-}
{-# LANGUAGE StandaloneDeriving   #-}
{-# LANGUAGE TupleSections        #-}
{-# LANGUAGE ViewPatterns         #-}
{-# OPTIONS_GHC -fno-warn-orphans #-}

import           Control.Concurrent
import           Control.Monad
import           Data.Bifunctor
import           Data.IORef
import           Data.List
import           Data.Maybe
import           GHC.TypeLits
import           Graphics.Vty hiding                 (Config)
import           Numeric.Hamilton
import           Numeric.LinearAlgebra.Static hiding (dim)
import           System.Exit
import           Text.Printf
import qualified Data.Map.Strict                     as M
import qualified Data.Vector.Sized                   as V
import qualified Data.Vector.Storable                as VS

deriving instance Ord Color

data SysExample where
    SE :: (KnownNat m, KnownNat n)
       => { seName   :: String
          , seCoords :: V.Vector n String
          , seSystem :: System m n
          , seDraw   :: R m -> [Point Double]
          , seInit   :: Phase n
          }
       -> SysExample

pendulum :: SysExample
pendulum = SE "Single pendulum" (V1 "θ") s f (toPhase s c0)
  where
    s :: System 2 1
    s = mkSystem (vec2 1 1                        )
                 (\(V1 θ)   -> V2 (sin θ) (-cos θ))
                 (\(V2 _ y) -> y                  )
    f :: R 2 -> [Point Double]
    f = (:[]) . r2pt
    c0 :: Config 1
    c0 = Cfg (0 :: R 1) (1 :: R 1)

doublePendulum :: Double -> Double -> SysExample
doublePendulum m1 m2 = SE "Double pendulum" (V2 "θ1" "θ2") s f (toPhase s c0)
  where
    s :: System 4 2
    s = mkSystem (vec4 m1 m1 m2 m2)
                 (\(V2 θ1 θ2)     -> V4 (sin θ1)            (-cos θ1)
                                        (sin θ1 + sin θ2/2) (-cos θ1 - cos θ2/2)
                 )
                 (\(V4 _ y1 _ y2) -> 5 * (realToFrac m1 * y1 + realToFrac m2 * y2))
    f :: R 4 -> [Point Double]
    f = (\(p1,p2) -> [p1,p2]) . bimap r2pt r2pt . split
    c0 :: Config 2
    c0 = Cfg (vec2 (pi/2) 0) (vec2 0 0)

room :: SysExample
room = SE "Room" (V2 "x" "y") s f (toPhase s c0)
  where
    s :: System 2 2
    s = mkSystem (vec2 1 1)
                 id
                 (\(V2 x y) -> sum [ 2 * y                      -- gravity
                                   , 1 - logistic (-1) 5 0.1 y  -- bottom wall
                                   , logistic 1 5 0.1 y         -- top wall
                                   , 1 - logistic (-2) 5 0.1 x  -- left wall
                                   , logistic 2 5 0.1 x         -- right wall
                                   ]
                 )
    f :: R 2 -> [Point Double]
    f = (:[]) . r2pt
    c0 :: Config 2
    c0 = Cfg (vec2 (-0.5) 0.5) (vec2 1 0)


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
    cfg <- standardIOConfig
    vty <- mkVty cfg

    opts <- newIORef $ SO 0.5 1 25

    -- t <- forkIO $ loop vty opts pendulum
    -- t <- forkIO $ loop vty opts (doublePendulum 1 1)
    t <- forkIO $ loop vty opts room

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
        qVec = intercalate ", " . V.toList $ seCoords
        go hists p = do
          SO{..} <- readIORef oRef
          let p'   = stepHam (soRate / fps) seSystem p  -- progress the simulation
              xb   = (- recip soZoom, recip soZoom)
              info = vertCat . map (string defAttr) $
                       [ printf "[ %s ]" seName
                       , printf "<%s>: <%s>" qVec . intercalate ", "
                          . map (printf "%.5f") . r2list . phsPos   $ p
                       , printf "KE: %.5f" . keP seSystem           $ p
                       , printf "PE: %.5f" . pe seSystem . phsPos   $ p
                       , printf "H: %.5f" . hamiltonian seSystem    $ p
                       , " "
                       , printf "rate: x%.2f" $ soRate
                       , printf "hist: %d"    $ soHist
                       ]
              pts  = (`zip` ptAttrs) . seDraw . underlyingPos seSystem . phsPos
                   $ p
              hists' = foldl' (\h (r, a) -> M.insertWith (addHist soHist) a [r] h) hists pts
          dr <- displayBounds $ outputIface vty
          update vty . picForLayers . (info:) . plot dr (PX xb (RR 0.5 screenRatio)) $
               ((second . second) (defAttr `withForeColor`) <$> pts)
            ++ (map (\((_,c),r) -> (r, ('.', defAttr `withForeColor` c)))
                  . concatMap sequence
                  . M.toList
                  $ hists'
               )
          threadDelay (round (1000000 / fps))
          go hists' p'
    addHist hl new old = take hl (new ++ old)

processEvt
    :: Event -> Maybe SimEvt
processEvt = \case
    EvKey KEsc        []      -> Just SEQuit
    EvKey (KChar 'c') [MCtrl] -> Just SEQuit
    EvKey (KChar 'q') []      -> Just SEQuit
    EvKey (KChar '+') []      -> Just $ SEZoom 2
    EvKey (KChar '-') []      -> Just $ SEZoom 0.5
    EvKey (KChar '>') []      -> Just $ SERate 1.5
    EvKey (KChar '<') []      -> Just $ SERate (1/1.5)
    EvKey (KChar ']') []      -> Just $ SEHist 1
    EvKey (KChar '[') []      -> Just $ SEHist (-1)
    EvKey (KChar '}') []      -> Just $ SEHist 5
    EvKey (KChar '{') []      -> Just $ SEHist (-5)
    _                         -> Nothing

data Point a = Pt { pX :: a, pY :: a }
             deriving (Show, Functor)

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
    -> [(Point Double, (Char, Attr))]   -- ^ points to plot
    -> [Image]
plot (wd,ht) pr = map (crop wd ht)
                . (++ bgs)
                . map (\(p, (c, a)) -> place EQ EQ p $ char a c)
  where
    wd' = fromIntegral wd
    ht' = fromIntegral ht
    ((xmin, xmax), (ymin, ymax)) = mkRange (wd', ht') pr
    origin = place EQ EQ (Pt 0 0) $ char defAttr '+'
    xaxis  = place EQ EQ (Pt 0 0) $ charFill defAttr '-' wd 1
    yaxis  = place EQ EQ (Pt 0 0) $ charFill defAttr '|' 1 ht
    xrange = xmax - xmin
    yrange = ymax - ymin
    bg     = backgroundFill wd ht
    scale Pt{..} = Pt x y
      where
        x = round $ (pX - xmin) * (wd' / xrange)
        y = round $ (pY - ymin) * (ht' / yrange)
    place aX aY (scale->Pt{..}) i
        = translate (fAlign aX (imageWidth  i))
                    (fAlign aY (imageHeight i))
        . translate pX pY
        $ i
    labels = [ place LT EQ (Pt xmin 0) . string defAttr $ printf "%.2f" xmin
             , place GT EQ (Pt xmax 0) . string defAttr $ printf "%.2f" xmax
             , place EQ LT (Pt 0 ymin) . string defAttr $ printf "%.2f" ymin
             , place EQ GT (Pt 0 ymax) . string defAttr $ printf "%.2f" ymax
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

pattern V2 :: a -> a -> V.Vector 2 a
pattern V2 x y <- (V.toList->[x,y])
  where
    V2 x y = fromJust (V.fromList [x,y])

pattern V4 :: a -> a -> a -> a -> V.Vector 4 a
pattern V4 x y z a <- (V.toList->[x,y,z,a])
  where
    V4 x y z a = fromJust (V.fromList [x,y,z,a])

r2list
    :: KnownNat n
    => R n
    -> [Double]
r2list = VS.toList . extract

r2pt
    :: R 2
    -> Point Double
r2pt (r2list->[x,y]) = Pt x y

logistic
    :: Floating a => a -> a -> a -> a -> a
logistic pos ht width = \x -> ht / (1 + exp (- beta * (x - pos)))
  where
    beta = log (0.9 / (1 - 0.9)) / width

