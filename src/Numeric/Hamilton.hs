{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE DeriveGeneric       #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE KindSignatures      #-}
{-# LANGUAGE RankNTypes          #-}
{-# LANGUAGE RecordWildCards     #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE StandaloneDeriving  #-}
{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE TypeInType          #-}
{-# LANGUAGE TypeOperators       #-}

-- |
-- Module      : Numeric.Hamilton
-- Description : Hamiltonian dynamics for physical systems on generalized
--               coordinates using automatic differentiation
-- Copyright   : (c) Justin Le 2016
-- License     : BSD-3
-- Maintainer  : justin@jle.im
-- Stability   : unstable
-- Portability : portable
--
-- Simulate physical systems on generalized/arbitrary coordinates using
-- Hamiltonian mechanics and automatic differentiation!
--
-- See the <https://github.com/mstksg/hamilton#readme> for more
-- information on usage!
--

module Numeric.Hamilton
  ( -- * Systems and states
    -- ** Systems
    System
  , mkSystem
  , mkSystem'
  , underlyingPos
    -- ** States
  , Config(..)
  , Phase(..)
  , toPhase
  , fromPhase
    -- * State functions
  , momenta
  , velocities
  , keC
  , keP
  , pe
  , lagrangian
  , hamiltonian
  , hamEqs
    -- * Simulating hamiltonian dynamics
    -- ** Over phase space
  , stepHam
  , evolveHam
  , evolveHam'
    -- ** Over configuration space
    -- | Convenience wrappers over the normal phase-space
    -- steppers/simulators that allow you to provide input and expect
    -- output in configuration space instead of in phase space.  Note that
    -- the simulation itself still runs in phase space, so these all
    -- require conversions to and from phase space under the hood.
  , stepHamC
  , evolveHamC
  , evolveHamC'
  ) where

import           Control.Monad
import           Data.Bifunctor
import           Data.Foldable
import           Data.Kind
import           Data.Maybe
import           Data.Proxy
import           Data.Type.Equality
import           GHC.Generics                        (Generic)
import           GHC.TypeLits
import           GHC.TypeLits.Compare
import           Numeric.AD
import           Numeric.GSL.ODE
import           Numeric.LinearAlgebra.Static        as H
import           Numeric.LinearAlgebra.Static.Vector
import qualified Data.Vector.Generic.Sized           as VG
import qualified Data.Vector.Sized                   as V
import qualified Numeric.LinearAlgebra               as LA

-- | Represents the full state of a system of @n@ generalized coordinates
-- in configuration space (informally, "positions and velocities")
--
-- A configuration space representaiton is more directly "physically
-- meaningful" and intuitive/understandable to humans than a phase space
-- representation.  However, it's much less mathematically ideal to work
-- with because of the lack of some neat underlying symmetries.
--
-- You can convert a @'Config' n@ into a @'Phase' n@ (convert from
-- configuration space to phase space) for a given system with 'toPhase'.
-- This allows you to state your system in configuration space and then
-- convert it to phase space before handing it off to the hamiltonian
-- machinery.
data Config :: Nat -> Type where
    Cfg :: { -- | The current values ("positions") of each of the @n@
             -- generalized coordinates
             cfgPositions :: !(R n)
             -- | The current rate of changes ("velocities") of each of the
             -- @n@ generalized coordinates
           , cfgVelocities :: !(R n)
           }
        -> Config n
  deriving (Generic)

deriving instance KnownNat n => Show (Config n)

-- | Represents the full state of a system of @n@ generalized coordinates
-- in phase space (informally, "positions and momentums").
--
-- Phase space representations are much nicer to work with mathematically
-- because of some neat underlying symmetries.  For one, positions and
-- momentums are "interchangeable" in a system; if you swap every
-- coordinate's positions with their momentums, and also swap them in the
-- equations of motions, you get the same system back.  This isn't the case
-- with configuration space representations.
--
-- A hamiltonian simulation basically describes the trajectory of each
-- coordinate through phase space, so this is the /state/ of the
-- simulation.  However, configuration space representations are much more
-- understandable to humans, so it might be useful to give an initial state
-- in configuration space using 'Config', and then convert it to a 'Phase'
-- with 'toPhase'.
data Phase :: Nat -> Type where
    Phs :: { -- | The current values ("positions") of each of the @n@
             -- generalized coordinates.
             phsPositions :: !(R n)
             -- | The current conjugate momenta ("momentums") to each of
             -- the @n@ generalized coordinates
           , phsMomenta :: !(R n)
           }
        -> Phase n
  deriving (Generic)

deriving instance KnownNat n => Show (Phase n)

-- | Represents a physical system in which physics happens.  A @'System'
-- m n@ is a system whose state described using @n@ generalized coordinates
-- (an "@n@-dimensional" system), where the underlying cartesian coordinate
-- space is @m@-dimensional.
--
-- For the most part, you are supposed to be able to ignore @m@.  @m@ is
-- only provided because it's useful when plotting/drawing the system with
-- a given state back in rectangular coordinates. (The only function that
-- use the @m@ at the moment is 'underlyingPos')
--
-- A @'System' m n@'s state is described using a @'Config' n@ (which
-- describes the system in configuration space) or a @'Phase' n@ (which
-- describes the system in phase space).
data System :: Nat -> Nat -> Type where
    Sys :: { _sysInertia       :: R m
           , _sysCoords        :: R n -> R m
           , _sysJacobian      :: R n -> L m n
           , _sysHessian       :: R n -> V.Vector n (L m n)
           , _sysPotential     :: R n -> Double
           , _sysPotentialGrad :: R n -> R n
           }
        -> System m n

-- | Converts the position of generalized coordinates of a system to the
-- coordinates of the system's underlying cartesian coordinate system.
-- Useful for plotting/drawing the system in cartesian space.
underlyingPos
    :: System m n
    -> R n
    -> R m
underlyingPos = _sysCoords

-- | The potential energy of a system, given the position in the
-- generalized coordinates of the system.
pe  :: System m n
    -> R n
    -> Double
pe = _sysPotential

vec2l
    :: (KnownNat m, KnownNat n)
    => V.Vector m (V.Vector n Double)
    -> L m n
vec2l = rowsL . fmap (vecR . VG.convert)

-- | Create a system with @n@ generalized coordinates by describing its
-- coordinate space (by a function from the generalized coordinates to the
-- underlying cartesian coordinates), the inertia of each of those
-- underlying coordinates, and the pontential energy function.
--
-- The potential energy function is expressed in terms of the genearlized
-- coordinate space's positions.
mkSystem
    :: forall m n. (KnownNat m, KnownNat n)
    => R m      -- ^ The "inertia" of each of the @m@ coordinates
                -- in the underlying cartesian space of the system.  This
                -- should be mass for linear coordinates and rotational
                -- inertia for angular coordinates.
    -> (forall a. RealFloat a => V.Vector n a -> V.Vector m a)
                -- ^ Conversion function to convert points in the
                -- generalized coordinate space to the underlying cartesian
                -- space of the system.
    -> (forall a. RealFloat a => V.Vector n a -> a)
                -- ^ The potential energy of the system as a function of
                -- the generalized coordinate space's positions.
    -> System m n
mkSystem m f u = Sys
    { _sysInertia       =                     m
    , _sysCoords        = vecR . VG.convert . f           . VG.convert . rVec
    , _sysJacobian      = tr   . vec2l      . jacobianT f . VG.convert . rVec
    , _sysHessian       = tr2  . fmap vec2l . hessianF f  . VG.convert . rVec
    , _sysPotential     =                     u           . VG.convert . rVec
    , _sysPotentialGrad = vecR . VG.convert . grad u      . VG.convert . rVec
    }
  where
    tr2 :: forall o. (KnownNat n, KnownNat o)
        => V.Vector m (L n o)
        -> V.Vector n (L m o)
    tr2 = fmap rowsL . traverse lRows
    {-# INLINE tr2 #-}


-- | Convenience wrapper over 'mkSystem' that allows you to specify the
-- potential energy function in terms of the underlying cartesian
-- coordinate space.
mkSystem'
    :: forall m n. (KnownNat m, KnownNat n)
    => R m      -- ^ The "inertia" of each of the @m@ coordinates
                -- in the underlying cartesian space of the system.  This
                -- should be mass for linear coordinates and rotational
                -- inertia for angular coordinates.
    -> (forall a. RealFloat a => V.Vector n a -> V.Vector m a)
                -- ^ Conversion function to convert points in the
                -- generalized coordinate space to the underlying cartesian
                -- space of the system.
    -> (forall a. RealFloat a => V.Vector m a -> a)
                -- ^ The potential energy of the system as a function of
                -- the underlying cartesian coordinate space's positions.
    -> System m n
mkSystem' m f u = mkSystem m f (u . f)


-- | Compute the generalized momenta conjugate to each generalized
-- coordinate of a system by giving the configuration-space state of the
-- system.
--
-- Note that getting the momenta from a @'Phase' n@ involves just using
-- 'phsMomenta'.
momenta
    :: (KnownNat m, KnownNat n)
    => System m n
    -> Config n
    -> R n
momenta Sys{..} Cfg{..} = tr j #> diag _sysInertia #> j #> cfgVelocities
  where
    j = _sysJacobian cfgPositions

-- | Convert a configuration-space representaiton of the state of the
-- system to a phase-space representation.
--
-- Useful because the hamiltonian simulations use 'Phase' as its working
-- state, but 'Config' is a much more human-understandable and intuitive
-- representation.  This allows you to state your starting state in
-- configuration space and convert to phase space for your simulation to
-- use.
toPhase
    :: (KnownNat m, KnownNat n)
    => System m n
    -> Config n
    -> Phase n
toPhase s = Phs <$> cfgPositions <*> momenta s

-- | The kinetic energy of a system, given the system's state in
-- configuration space.
keC :: (KnownNat m, KnownNat n)
    => System m n
    -> Config n
    -> Double
keC s = do
    vs <- cfgVelocities
    ps <- momenta s
    return $ (vs <.> ps) / 2

-- | The Lagrangian of a system (the difference between the kinetic energy
-- and the potential energy), given the system's state in configuration
-- space.
lagrangian
    :: (KnownNat m, KnownNat n)
    => System m n
    -> Config n
    -> Double
lagrangian s = do
    t <- keC s
    u <- pe s . cfgPositions
    return (t - u)

-- | Compute the rate of change of each generalized coordinate by giving
-- the state of the system in phase space.
--
-- Note that getting the velocities from a @'Config' n@ involves just using
-- 'cfgVelocities'.
velocities
    :: (KnownNat m, KnownNat n)
    => System m n
    -> Phase n
    -> R n
velocities Sys{..} Phs{..} = inv jmj #> phsMomenta
  where
    j   = _sysJacobian phsPositions
    jmj = tr j H.<> diag _sysInertia H.<> j

-- | Invert 'toPhase' and convert a description of a system's state in
-- phase space to a description of the system's state in configuration
-- space.
--
-- Possibly useful for showing the phase space representation of a system's
-- state in a more human-readable/human-understandable way.
fromPhase
    :: (KnownNat m, KnownNat n)
    => System m n
    -> Phase n
    -> Config n
fromPhase s = Cfg <$> phsPositions <*> velocities s

-- | The kinetic energy of a system, given the system's state in
-- phase space.
keP :: (KnownNat m, KnownNat n)
    => System m n
    -> Phase n
    -> Double
keP s = do
    ps <- phsMomenta
    vs <- velocities s
    return $ (vs <.> ps) / 2

-- | The Hamiltonian of a system (the sum of kinetic energy and the
-- potential energy), given the system's state in phase space.
hamiltonian
    :: (KnownNat m, KnownNat n)
    => System m n
    -> Phase n
    -> Double
hamiltonian s = do
    t <- keP s
    u <- pe s . phsPositions
    return (t + u)

-- | The "hamiltonian equations" for a given system at a given state in
-- phase space.  Returns the rate of change of the positions and
-- conjugate momenta, which can be used to progress the simulation through
-- time.
--
-- Computed using the maths derived in
-- <https://blog.jle.im/entry/hamiltonian-dynamics-in-haskell.html>.
hamEqs
    :: (KnownNat m, KnownNat n)
    => System m n
    -> Phase n
    -> (R n, R n)
hamEqs Sys{..} Phs{..} = (dHdp, -dHdq)
  where
    mm   = diag _sysInertia
    j    = _sysJacobian phsPositions
    trj  = tr j
    jmj  = trj H.<> mm H.<> j
    ijmj = inv jmj
    dTdq = vecR . VG.convert
         . flip fmap (_sysHessian phsPositions) $ \djdq ->
             -phsMomenta <.> ijmj #> trj #> mm #> djdq #> ijmj #> phsMomenta
    dHdp = ijmj #> phsMomenta
    dHdq = dTdq + _sysPotentialGrad phsPositions

-- | Step a system through phase space over over a single timestep.
stepHam
    :: forall m n. (KnownNat m, KnownNat n)
    => Double           -- ^ timestep to step through
    -> System m n       -- ^ system to simulate
    -> Phase n          -- ^ initial state, in phase space
    -> Phase n
stepHam r s p = evolveHam @m @n @2 s p (fromJust $ V.fromList [0, r])
                  `V.unsafeIndex` 1

-- | Evolve a system using a hamiltonian stepper, with the given initial
-- phase space state.
--
-- Desired solution times provided as a list instead of a sized 'V.Vector'.
-- The output list should be the same length as the input list.
evolveHam'
    :: forall m n. (KnownNat m, KnownNat n)
    => System m n  -- ^ system to simulate
    -> Phase n     -- ^ initial state, in phase space
    -> [Double]    -- ^ desired solution times
    -> [Phase n]
evolveHam' _ _ [] = []
evolveHam' s p0 ts = V.withSizedList (toList ts') $ \(v :: V.Vector s Double) ->
                       case Proxy @2 %<=? Proxy @s of
                         LE Refl -> (if l1 then tail else id)
                                  . toList
                                  $ evolveHam s p0 v
                         NLE{}   -> error "evolveHam': Internal error"
  where
    (l1, ts') = case ts of
      [x] -> (True , [0,x])
      _   -> (False, ts   )

-- | Evolve a system using a hamiltonian stepper, with the given initial
-- phase space state.
evolveHam
    :: forall m n s. (KnownNat m, KnownNat n, KnownNat s, 2 <= s)
    => System m n           -- ^ system to simulate
    -> Phase n              -- ^ initial state, in phase space
    -> V.Vector s Double    -- ^ desired solution times
    -> V.Vector s (Phase n)
evolveHam s p0 ts = fmap toPs . fromJust . V.fromList . LA.toRows
                  $ odeSolveV RKf45 hi eps eps (const f) (fromPs p0) ts'
  where
    hi  = (V.unsafeIndex ts 1 - V.unsafeIndex ts 0) / 100
    eps = 1.49012e-08
    f :: LA.Vector Double -> LA.Vector Double
    f   = uncurry (\p m -> LA.vjoin [p,m])
        . join bimap extract . hamEqs s . toPs
    ts' = VG.fromSized . VG.convert $ ts
    n = fromInteger $ natVal (Proxy @n)
    fromPs :: Phase n -> LA.Vector Double
    fromPs p = LA.vjoin . map extract $ [phsPositions p, phsMomenta p]
    toPs :: LA.Vector Double -> Phase n
    toPs v = Phs pP pM
      where
        Just [pP, pM] = traverse create . LA.takesV [n, n] $ v

-- | A convenience wrapper for 'evolveHam'' that works on configuration
-- space states instead of phase space states.
--
-- Note that the simulation itself still runs in phase space; this function
-- just abstracts over converting to and from phase space for the inputs
-- and outputs.
evolveHamC'
    :: forall m n. (KnownNat m, KnownNat n)
    => System m n       -- ^ system to simulate
    -> Config n         -- ^ initial state, in configuration space
    -> [Double]         -- ^ desired solution times
    -> [Config n]
evolveHamC' s c0 = fmap (fromPhase s) . evolveHam' s (toPhase s c0)

-- | A convenience wrapper for 'evolveHam' that works on configuration
-- space states instead of phase space states.
--
-- Note that the simulation itself still runs in phase space; this function
-- just abstracts over converting to and from phase space for the inputs
-- and outputs.
evolveHamC
    :: forall m n s. (KnownNat m, KnownNat n, KnownNat s, 2 <= s)
    => System m n           -- ^ system to simulate
    -> Config n             -- ^ initial state, in configuration space
    -> V.Vector s Double    -- ^ desired solution times
    -> V.Vector s (Config n)
evolveHamC s c0 = fmap (fromPhase s) . evolveHam s (toPhase s c0)

-- | Step a system through configuration space over over a single timestep.
--
-- Note that the simulation itself still runs in phase space; this function
-- just abstracts over converting to and from phase space for the input
-- and output.
stepHamC
    :: forall m n. (KnownNat m, KnownNat n)
    => Double           -- ^ timestep to step through
    -> System m n       -- ^ system to simulate
    -> Config n         -- ^ initial state, in phase space
    -> Config n
stepHamC r s = fromPhase s . stepHam r s . toPhase s

