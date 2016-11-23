{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE KindSignatures      #-}
{-# LANGUAGE RankNTypes          #-}
{-# LANGUAGE RecordWildCards     #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeInType          #-}

module Numeric.Hamilton
  ( Config(..)
  , Phase(..)
  , System(..)
  , mkSystem
  , mkPhase
  , momenta
  , velocities
  , ke
  , keP
  , lagrangian
  , hamiltonian
  ) where

import           Data.Foldable
import           Data.Kind
import           Data.Maybe
import           GHC.TypeLits
import           Numeric.AD
import           Numeric.LinearAlgebra.Static
import qualified Control.Comonad              as C
import qualified Control.Comonad.Cofree       as C
import qualified Data.Vector.Generic.Sized    as VG
import qualified Data.Vector.Sized            as V

data Config :: Nat -> Type where
    Cfg :: { cfgPos :: !(R n)
           , cfgVel :: !(R n)
           }
        -> Config n

data Phase :: Nat -> Type where
    Phs :: { phsPos :: !(R n)
           , phsMom :: !(R n)
           }
        -> Phase n

-- TODO: cache inverses, JMJs?
data System :: Nat -> Nat -> Type where
    Sys :: { sysMass          :: R m
           , sysCoords        :: R n -> R m
           , sysJacobian      :: R n -> L m n
           , sysJacobian2     :: R n -> V.Vector m (Sym n)
           , sysPotential     :: R m -> Double
           , sysPotentialGrad :: R m -> R m
           }
        -> System m n

vec2r
    :: KnownNat n => V.Vector n Double -> R n
vec2r = fromJust . create . VG.fromSized . VG.convert

r2vec
    :: KnownNat n => R n -> V.Vector n Double
r2vec = VG.convert . fromJust . VG.toSized . extract

vec2l
    :: (KnownNat m, KnownNat n)
    => V.Vector m (V.Vector n Double)
    -> L m n
vec2l = fromJust . (\rs -> withRows rs exactDims) . toList . fmap vec2r

mkSystem
    :: forall m n. (KnownNat m, KnownNat n)
    => R m
    -> (forall a. RealFloat a => V.Vector n a -> V.Vector m a)
    -> (forall a. RealFloat a => V.Vector m a -> a)
    -> System m n
mkSystem m f u = Sys m
                     (vec2r . f . r2vec)
                     (vec2l . jacobian f . r2vec)
                     (fmap (sym . vec2l . j2 . C.hoistCofree VG.convert)
                        . VG.convert
                        . jacobians f
                        . r2vec
                        )
                     (u . r2vec)
                     (vec2r . grad u . r2vec)
  where
    j2  :: C.Cofree (V.Vector n) Double
        -> V.Vector n (V.Vector n Double)
    j2 = fmap (fmap C.extract . C.unwrap) . C.unwrap

momenta
    :: (KnownNat m, KnownNat n)
    => System m n
    -> Config n
    -> R n
momenta Sys{..} Cfg{..} = tr j #> diag sysMass #> j #> cfgVel
  where
    j = sysJacobian cfgPos

mkPhase
    :: (KnownNat m, KnownNat n)
    => System m n
    -> Config n
    -> Phase n
mkPhase s = Phs <$> cfgPos <*> momenta s

ke  :: (KnownNat m, KnownNat n)
    => System m n
    -> Config n
    -> Double
ke s = do
    vs <- cfgVel
    ps <- momenta s
    return $ (vs <.> ps) / 2

lagrangian
    :: (KnownNat m, KnownNat n)
    => System m n
    -> Config n
    -> Double
lagrangian s = do
    t <- ke s
    u <- sysPotential s . sysCoords s . cfgPos
    return (t - u)

velocities
    :: (KnownNat m, KnownNat n)
    => System m n
    -> Phase n
    -> R n
velocities Sys{..} Phs{..} = inv jmj #> phsMom
  where
    j   = sysJacobian phsPos
    jmj = tr j <> diag sysMass <> j

keP :: (KnownNat m, KnownNat n)
    => System m n
    -> Phase n
    -> Double
keP s = do
    ps <- phsMom
    vs <- velocities s
    return $ (vs <.> ps) / 2

hamiltonian
    :: (KnownNat m, KnownNat n)
    => System m n
    -> Phase n
    -> Double
hamiltonian s = do
    t <- keP s
    u <- sysPotential s . sysCoords s . phsPos
    return (t + u)
