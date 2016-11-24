{-# LANGUAGE DeriveFoldable  #-}
{-# LANGUAGE DeriveFunctor   #-}
{-# LANGUAGE LambdaCase      #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE ViewPatterns    #-}

import           Control.Concurrent
import           Control.Monad
import           Data.Fixed
import           Numeric
import           Graphics.Vty
import           System.Exit
import           Text.Printf

data Stream a = a :< Stream a
  deriving (Functor, Foldable)

fps :: Double
fps = 20

main :: IO ()
main = do
    cfg <- standardIOConfig
    vty <- mkVty cfg

    t <- forkIO . forM_ (stream succ 0) $ \i -> do
      dr <- displayBounds $ outputIface vty
      update vty . picForLayers $
        plot dr (PX (-5, 5) (RR 0.5 2)) [(Pt (-5 + (i / 20) `mod'` 10) 1, '*')]
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
plot (wd,ht) pr = (++ bgs) . map (\(p, c) -> place EQ EQ p $ char mempty c)
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


unfold
    :: (s -> (a, s))
    -> s
    -> Stream a
unfold f = go
  where
    go s = x :< go s'
      where
        (x, s') = f s

stream
    :: (a -> a)
    -> a
    -> Stream a
stream f = unfold (\s -> (s, f s))
