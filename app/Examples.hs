{-# LANGUAGE DataKinds       #-}
{-# LANGUAGE DeriveFoldable  #-}
{-# LANGUAGE DeriveFunctor   #-}
{-# LANGUAGE GADTs           #-}
{-# LANGUAGE LambdaCase      #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE TupleSections   #-}
{-# LANGUAGE ViewPatterns    #-}

-- import           Data.Fixed
-- import           Data.Foldable
import           Control.Concurrent
import           Control.Monad
import           Data.Maybe
import           GHC.TypeLits
import           Graphics.Vty hiding          (Config)
import           Numeric.Hamilton
import           Numeric.LinearAlgebra.Static
import           System.Exit
import           Text.Printf
import qualified Data.Vector.Sized            as V
import qualified Data.Vector.Storable         as VS

fps :: Double
fps = 20

pendulum :: SysExample
pendulum = SE s f (toPhase s c0)
  where
    s :: System 2 1
    s = mkSystem (vec2 1 1) (fromJust . V.fromList . (\θ -> [cos θ, sin θ]) . V.head) (`V.unsafeIndex` 1)
    f :: R 2 -> Point Double
    f = (\[x,y] -> Pt x y) . VS.toList . extract
    c0 :: Config 1
    c0 = Cfg (vector [-pi/4]) (vector [0])

data SysExample where
    SE :: (KnownNat m, KnownNat n)
       => { seSystem :: System m n
          , seDraw   :: R m -> Point Double
          , seInit   :: Phase n
          }
       -> SysExample

main :: IO ()
main = do
    cfg <- standardIOConfig
    vty <- mkVty cfg

    SE{..} <- return pendulum
    -- let f r p = evolveHam' seSystem p [0, r] !! 1
    let f r = (!! 100) . iterate (stepHam (r / 100) seSystem)

    t <- forkIO . forM_ (iterate (f (1/fps)) seInit) $ \p -> do
      let pt  = seDraw . underlyingPos seSystem . phsPos $ p
      dr <- displayBounds $ outputIface vty
      let disp = vertCat [ string mempty . printf "%.5f" . head . VS.toList . extract . phsPos $ p
                         , string mempty . printf "%.5f" . keP seSystem $ p
                         , string mempty . printf "%.5f" . pe seSystem . phsPos $ p
                         , string mempty . printf "%.5f" . hamiltonian seSystem $ p
                         ]
      update vty . picForLayers . (disp:) $
        plot dr (PX (-2, 2) (RR 0.5 2)) [(pt, '*')]
      threadDelay (round (1000000 / fps))

    forever $ do
      e <- nextEvent vty
      when (quitter e) $ do
        killThread t
        shutdown vty
        exitSuccess

quitter :: Event -> Bool
quitter = \case
    EvKey KEsc        []      -> True
    EvKey (KChar 'c') [MCtrl] -> True
    EvKey (KChar 'q') []      -> True
    _                         -> False

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
    :: (Int, Int)         -- ^ display bounds
    -> PlotRange
    -> [(Point Double, Char)] -- ^ points to plot
    -> [Image]
plot (wd,ht) pr = map (crop wd ht)
                . (++ bgs)
                . map (\(p, c) -> place EQ EQ p $ char mempty c)
  where
    wd' = fromIntegral wd
    ht' = fromIntegral ht
    ((xmin, xmax), (ymin, ymax)) = mkRange (wd', ht') pr
    origin = place EQ EQ (Pt 0 0) $ char mempty '+'
    xaxis  = place EQ EQ (Pt 0 0) $ charFill mempty '-' wd 1
    yaxis  = place EQ EQ (Pt 0 0) $ charFill mempty '|' 1 ht
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
    labels = [ place LT EQ (Pt xmin 0) . string mempty $ printf "%.2f" xmin
             , place GT EQ (Pt xmax 0) . string mempty $ printf "%.2f" xmax
             , place EQ LT (Pt 0 ymin) . string mempty $ printf "%.2f" ymin
             , place EQ GT (Pt 0 ymax) . string mempty $ printf "%.2f" ymax
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
