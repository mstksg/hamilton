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
{-# LANGUAGE ViewPatterns        #-}

module Numeric.Hamilton
  ( System
  , Config(..)
  , Phase(..)
  , mkSystem
  , underlyingPos
  , underlyingInertia
  , pe
  , toPhase
  , fromPhase
  , momenta
  , velocities
  , keC
  , keP
  , lagrangian
  , hamiltonian
  , hamEqs
  , stepHam
  , evolveHam
  , evolveHam'
  , evolveHamC
  , evolveHamC'
  ) where

import           Control.Monad
import           Data.Bifunctor
import           Data.Foldable
import           Data.Kind
import           Data.Maybe
import           Data.Proxy
import           GHC.Generics                 (Generic)
import           GHC.TypeLits
import           Numeric.AD
import           Numeric.GSL.ODE
import           Numeric.LinearAlgebra.Static
import qualified Control.Comonad              as C
import qualified Control.Comonad.Cofree       as C
import qualified Data.Vector.Generic.Sized    as VG
import qualified Data.Vector.Sized            as V
import qualified Numeric.LinearAlgebra        as LA

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
             cfgPos :: !(R n)
             -- | The current rate of changes ("velocities") of each of the
             -- @n@ generalized coordinates
           , cfgVel :: !(R n)
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
             phsPos :: !(R n)
             -- | The current conjugate momenta ("momentums") to each of
             -- the @n@ generalized coordinates
           , phsMom :: !(R n)
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
-- a given state back in rectangular coordinates.
--
-- A @'System' m n@'s state is described using a @'Config' n@ (which
-- describes the system in configuration space) or a @'Phase' n@ (which
-- describes the system in phase space).
data System :: Nat -> Nat -> Type where
    Sys :: { _sysInertia       :: R m
           , _sysCoords        :: R n -> R m
           , _sysJacobian      :: R n -> L m n
           , _sysJacobian2     :: R n -> V.Vector m (Sym n)
           , _sysPotential     :: R m -> Double
           , _sysPotentialGrad :: R m -> R m
           }
        -> System m n

underlyingPos
    :: System m n
    -> R n
    -> R m
underlyingPos = _sysCoords

underlyingInertia
    :: System m n
    -> R m
underlyingInertia = _sysInertia

pe  :: System m n
    -> R n
    -> Double
pe s = _sysPotential s . underlyingPos s

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

-- l2vec
--     :: (KnownNat m, KnownNat n)
--     => L m n
--     -> V.Vector m (V.Vector n Double)
-- l2vec = fromJust . V.fromList . map r2vec . toRows

-- | Create a system with @n@ generalized coordinates by providing the
-- underlying inertials, describing its coordinate space, and giving
-- pontential energy function.
--
-- Given by describing some fundamental properties of the underlying
-- cartesian coordinate space of the system (the inertia of each coordinate
-- and the potential energy function), along with a conversion function
-- from the generalized coordinate space to the underlying cartesian
-- coordinate space.
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
    -> (forall a. RealFloat a => V.Vector m a -> a)
                -- ^ The potential energy of the system as a function of
                -- points in the underlying cartesian space of the system.
    -> System m n
mkSystem m f u = Sys m
                     (vec2r . f . r2vec)
                     (tr . vec2l . jacobianT f . r2vec)
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

-- | Compute the generalized momenta conjugate to each generalized
-- coordinate of a system by giving the configuration-space state of the
-- system.
momenta
    :: (KnownNat m, KnownNat n)
    => System m n
    -> Config n
    -> R n
momenta Sys{..} Cfg{..} = tr j #> diag _sysInertia #> j #> cfgVel
  where
    j = _sysJacobian cfgPos

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
toPhase s = Phs <$> cfgPos <*> momenta s

-- | The kinetic energy of a system, given the system's state in
-- configuration space.
keC :: (KnownNat m, KnownNat n)
    => System m n
    -> Config n
    -> Double
keC s = do
    vs <- cfgVel
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
    u <- pe s . cfgPos
    return (t - u)

-- | Compute the rate of change of each generalized coordinate by giving
-- the state of the system in phase space.
velocities
    :: (KnownNat m, KnownNat n)
    => System m n
    -> Phase n
    -> R n
velocities Sys{..} Phs{..} = inv jmj #> phsMom
  where
    j   = _sysJacobian phsPos
    jmj = tr j <> diag _sysInertia <> j

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
fromPhase s = Cfg <$> phsPos <*> velocities s

-- | The kinetic energy of a system, given the system's state in
-- phase space.
keP :: (KnownNat m, KnownNat n)
    => System m n
    -> Phase n
    -> Double
keP s = do
    ps <- phsMom
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
    u <- pe s . phsPos
    return (t + u)

hamEqs
    :: (KnownNat m, KnownNat n)
    => System m n
    -> Phase n
    -> (R n, R n)
hamEqs Sys{..} Phs{..} = (dHdp, -dHdq)
  where
    mm   = diag _sysInertia
    j    = _sysJacobian phsPos
    trj  = tr j
    j'   = unSym <$> _sysJacobian2 phsPos
    jmj  = trj <> mm <> j
    ijmj = inv jmj
    dTdq = vec2r
         . flip fmap (tr2 j') $ \djdq ->
             (phsMom <.> tr ijmj #> tr djdq #> mm #> djdq #> ijmj #> phsMom) / 2
    dHdp = ijmj #> phsMom
    dHdq = dTdq + trj #> _sysPotentialGrad (_sysCoords phsPos)

stepHam
    :: (KnownNat m, KnownNat n)
    => Double
    -> System m n
    -> Phase n
    -> Phase n
stepHam r s p@Phs{..} = Phs (phsPos + konst r * dq) (phsMom + konst r * dp)
  where
    (dq, dp) = hamEqs s p

tr2
    :: (KnownNat m, KnownNat n, KnownNat o)
    => V.Vector m (L n o)
    -> V.Vector n (L m o)
tr2 = fmap (fromJust . (\rs -> withRows rs exactDims) . toList)
    . sequenceA
    . fmap (fromJust . V.fromList . toRows)

-- | Evolve a system using a hamiltonian stepper, with the given initial
-- phase space state.
--
-- Desired solution times provided as a list instead of a sized 'V.Vector'.
-- The output list should be the same length as the input list.
evolveHam'
    :: forall m n. (KnownNat m, KnownNat n)
    => System m n
    -> Phase n
    -> [Double]
    -> [Phase n]
evolveHam' s p0 ts = V.withSizedList ts (toList . evolveHam s p0)

-- | Evolve a system using a hamiltonian stepper, with the given initial
-- phase space state.
evolveHam
    :: forall m n s. (KnownNat m, KnownNat n, KnownNat s)
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
    fromPs p = LA.vjoin . map extract $ [phsPos p, phsMom p]
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
    => System m n           -- ^ system to simulate
    -> Config n             -- ^ initial state, in configuration space
    -> [Double]             -- ^ desired solution times
    -> [Config n]
evolveHamC' s c0 = fmap (fromPhase s) . evolveHam' s (toPhase s c0)

-- | A convenience wrapper for 'evolveHam' that works on configuration
-- space states instead of phase space states.
--
-- Note that the simulation itself still runs in phase space; this function
-- just abstracts over converting to and from phase space for the inputs
-- and outputs.
evolveHamC
    :: forall m n s. (KnownNat m, KnownNat n, KnownNat s)
    => System m n           -- ^ system to simulate
    -> Config n             -- ^ initial state, in configuration space
    -> V.Vector s Double    -- ^ desired solution times
    -> V.Vector s (Config n)
evolveHamC s c0 = fmap (fromPhase s) . evolveHam s (toPhase s c0)
