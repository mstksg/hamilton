{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE DeriveFoldable      #-}
{-# LANGUAGE DeriveFunctor       #-}
{-# LANGUAGE DeriveGeneric       #-}
{-# LANGUAGE DeriveTraversable   #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE KindSignatures      #-}
{-# LANGUAGE LambdaCase          #-}
{-# LANGUAGE RankNTypes          #-}
{-# LANGUAGE RecordWildCards     #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE StandaloneDeriving  #-}
{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE TypeInType          #-}
{-# LANGUAGE TypeOperators       #-}
{-# LANGUAGE ViewPatterns        #-}

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
import           Data.Functor.Compose
import           Data.Kind
import           Data.Maybe
import           Data.Proxy
import           Data.Type.Equality hiding           (sym)
import           GHC.Generics                        (Generic)
import           GHC.TypeLits
import           GHC.TypeLits.Compare
import           GHC.TypeLits.Witnesses
import           Numeric.AD
import           Numeric.GSL.ODE
import           Numeric.LinearAlgebra.Static hiding ((&))
import qualified Control.Comonad                     as C
import qualified Control.Comonad.Cofree              as C
import qualified Data.Vector.Generic.Sized           as VG
import qualified Data.Vector.Sized                   as V
import qualified Data.Vector.Storable                as VS
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
    Cfg :: { -- | The current time of the state of the system
             cfgTime       :: !Double
           , -- | The current values ("positions") of each of the @n@
             -- generalized coordinates
             cfgPositions  :: !(R n)
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
    Phs :: { -- | The current time of the state of the system
             phsTime      :: !Double
             -- | The current values ("positions") of each of the @n@
             -- generalized coordinates.
           , phsPositions :: !(R n)
             -- | The current conjugate momenta ("momentums") to each of
             -- the @n@ generalized coordinates
           , phsMomenta   :: !(R n)
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
    Sys :: { _sysInertia        :: R m
           , _sysCoords         :: Double -> R n -> R m
           , _sysJacobian       :: Double -> R n -> (R m, L m n)
           , _sysJacobian2      :: Double -> R n -> V.Vector m (Sym' n)
           , _sysPotential      :: Double -> R n -> Double
           , _sysPotentialGrad  :: Double -> R n -> (Double, R n)
           }
        -> System m n

-- | Helper data type to describe a symmetric matrix in block form, where
-- the last row/column and lower corner are stored separately.
data Sym' n = Sym' { sUpper  :: !(Sym n)
                   , sCross  :: !(R n)
                   , sCorner :: !(Double)
                   }

-- sysJacobian
--     :: (KnownNat m, KnownNat n, KnownNat (n + 1))
--     => System m n
--     -> Double
--     -> R n
--     -> L m (n + 1)
-- sysJacobian Sys{..} t x = jx ||| col jt
--   where
--     (jt, jx) = _sysJacobian t x


-- coordShift
--     :: (KnownNat m, KnownNat n, KnownNat o)
--     => (R o -> R n)
--     -> (R o -> L n o)
--     -> (R o -> V.Vector n (Sym o))
--     -> System m n
--     -> System m o
-- coordShift c j j2 = \case
--     Sys i c0 j0 j20 p g -> Sys i (c0 . c)
--                                ((<>) <$> j0 . c <*> j)
--                                ((\d -> fmap _) <$> j2 <*> j20 . c)
--                                p g

-- | Converts the position of generalized coordinates of a system to the
-- coordinates of the system's underlying cartesian coordinate system.
-- Useful for plotting/drawing the system in cartesian space.
underlyingPos
    :: System m n
    -> Double
    -> R n
    -> R m
underlyingPos = _sysCoords

-- | The potential energy of a system, given the position in the
-- generalized coordinates of the system.
pe  :: System m n
    -> Double
    -> R n
    -> Double
pe = _sysPotential

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

l2vec
    :: (KnownNat m, KnownNat n)
    => L m n
    -> V.Vector m (V.Vector n Double)
l2vec = fromJust . V.fromList . map r2vec . toRows

data Chron f a = Chron { chronTime :: !a
                       , chronVal  :: !(f a)
                       }
    deriving (Show, Eq, Ord, Functor, Foldable, Traversable)

hoistChron
    :: (f a -> g a)
    -> Chron f a
    -> Chron g a
hoistChron f (Chron x fx) = Chron x (f fx)

unChron
    :: Chron f a
    -> (a, f a)
unChron Chron{..} = (chronTime, chronVal)

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
    -> (forall a. RealFloat a => a -> V.Vector n a -> V.Vector m a)
                -- ^ Conversion function to convert points in the
                -- generalized coordinate space to the underlying cartesian
                -- space of the system.
    -> (forall a. RealFloat a => a -> V.Vector n a -> a)
                -- ^ The potential energy of the system as a function of
                -- the generalized coordinate space's positions.
    -> System m n
mkSystem m f u =
    Sys m
        (\t -> vec2r . f t . r2vec)
        (\t x -> let jt :: V.Vector m Double
                     jx :: V.Vector n (V.Vector m Double)
                     Chron jt jx = jacobianT f' (Chron t (r2vec x))
                 in  (vec2r jt, tr (vec2l jx))
        )
        (\t x -> fmap ( (\(t2, jt, jx) -> Sym' (sym (vec2l jx)) (vec2r jt) t2)
                      . j2
                      . C.hoistCofree (hoistChron VG.convert)
                      )
               . VG.convert
               $ jacobians f' (Chron t (r2vec x))
        )
        (\t -> u t . r2vec)
        (\t -> second vec2r . unChron . grad u' . Chron t . r2vec)
  where
    f'  :: RealFloat a => Chron (V.Vector n) a -> V.Vector m a
    f' (Chron t x) = f t x
    u'  :: RealFloat a => Chron (V.Vector n) a -> a
    u' (Chron t x) = u t x
    j2  :: C.Cofree (Chron (V.Vector n)) Double
        -> (Double, V.Vector n Double, V.Vector n (V.Vector n Double))
    j2 c = (t2, jt, chronVal <$> jx)
      where
        t2 :: Double
        jt :: V.Vector n Double
        jx :: V.Vector n (Chron (V.Vector n) Double)
        Chron (Chron t2 jt) jx = fmap (fmap C.extract . C.unwrap)
                               . C.unwrap
                               $ c

-- | Convenience wrapper over 'mkSystem' that allows you to specify the
-- potential energy function in terms of the underlying cartesian
-- coordinate space.
mkSystem'
    :: forall m n. (KnownNat m, KnownNat n)
    => R m      -- ^ The "inertia" of each of the @m@ coordinates
                -- in the underlying cartesian space of the system.  This
                -- should be mass for linear coordinates and rotational
                -- inertia for angular coordinates.
    -> (forall a. RealFloat a => a -> V.Vector n a -> V.Vector m a)
                -- ^ Conversion function to convert points in the
                -- generalized coordinate space to the underlying cartesian
                -- space of the system.
    -> (forall a. RealFloat a => a -> V.Vector m a -> a)
                -- ^ The potential energy of the system as a function of
                -- the underlying cartesian coordinate space's positions.
    -> System m n
mkSystem' m f u = mkSystem m f (\t -> u t . f t)


-- | Compute the generalized momenta conjugate to each generalized
-- coordinate of a system by giving the configuration-space state of the
-- system.
--
-- Note that getting the momenta from a @'Phase' n@ involves just using
-- 'phsMomenta'.
momenta
    :: forall m n. (KnownNat m, KnownNat n)
    => System m n
    -> Config n
    -> R n
momenta s Cfg{..} = case sysJMJ s cfgTime cfgPositions of
    Sym'{..} -> (unSym sUpper #> cfgPositions) + sCross

quadForm
    :: KnownNat n
    => R n
    -> L n n
    -> Double
quadForm x a = x <.> a #> x

sysJMJ
    :: (KnownNat m, KnownNat n)
    => System m n
    -> Double
    -> R n
    -> Sym' n
sysJMJ Sys{..} t x = Sym' (sym $ tr jx <> m <> jx) (tr jx #> mt) (jt <.> mt)
  where
    m = diag _sysInertia
    (jt, jx) = _sysJacobian t x
    mt = m #> jt

sysJMJ'
    :: (KnownNat m, KnownNat n)
    => System m n
    -> Double
    -> R n
    -> V.Vector n (Sym n, R n, Double)
sysJMJ' Sys{..} t x = flip fmap (chronVal (trJ jj)) $ \(jjx, jjt) ->
    let upper  = 2 * tr jjx <> m <> jjx
        -- why is there no <# ?
        cross  = tr (m <> jx) #> jjt + tr (m <> jjx) #> jt
        corner = 2 * jt <.> m #> jjt
    in  (sym upper, cross, corner)
  where
    m        = diag _sysInertia
    (jt, jx) = _sysJacobian t x
    jj       = _sysJacobian2 t x

trJ
    :: (KnownNat m, KnownNat n)
    => V.Vector m (Sym' n)
    -> Chron (V.Vector n) (L m n, R m)
trJ j0 = Chron jt jx
  where
    jx = fmap (bimap vec2l vec2r . V.unzip) . sequenceA
       . fmap (\Sym'{..} -> (,) <$> l2vec (unSym sUpper) <*> r2vec sCross)
       $ j0
    jt = bimap vec2l vec2r . V.unzip
       . fmap (\Sym'{..} -> (r2vec sCross, sCorner))
       $ j0

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
toPhase s = Phs <$> cfgTime <*> cfgPositions <*> momenta s

-- | The kinetic energy of a system, given the system's state in
-- configuration space.
keC :: (KnownNat m, KnownNat n)
    => System m n
    -> Config n
    -> Double
keC s = do
    vs <- cfgVelocities
    Sym'{..} <- sysJMJ s <$> cfgTime <*> cfgPositions
    return $ (quadForm vs (unSym sUpper) + 2 * vs <.> sCross + sCorner) / 2

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
    u <- pe s <$> cfgTime <*> cfgPositions
    return (t - u)

-- | Compute the rate of change of each generalized coordinate by giving
-- the state of the system in phase space.
--
-- Note that getting the velocities from a @'Config' n@ involves just using
-- 'cfgVelocities'.
velocities
    :: forall m n. (KnownNat m, KnownNat n)
    => System m n
    -> Phase n
    -> R n
velocities s Phs{..} = case sysJMJ s phsTime phsPositions of
    Sym'{..} -> inv (unSym sUpper) #> (phsMomenta - sCross)

-- momenta s Cfg{..} = case jmj s cfgTime cfgPositions of
--     Sym'{..} -> (unSym sUpper #> cfgPositions) + sCross
-- velocities s Phs{..} = withNatOp (%+) (Proxy @n) (Proxy @1) $
--     let j   = sysJacobian s phsTime phsPositions
--         jmj = tr j <> diag (_sysInertia s) <> j
--     in  initR $ inv jmj #> (phsMomenta & 1)

    -- let j = sysJacobian s cfgTime cfgPositions
    -- in  initR $ tr j #> diag (_sysInertia s) #> j #> (cfgVelocities & 1)

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
fromPhase s = Cfg <$> phsTime <*> phsPositions <*> velocities s

-- | The kinetic energy of a system, given the system's state in
-- phase space.
keP :: (KnownNat m, KnownNat n)
    => System m n
    -> Phase n
    -> Double
keP s = do
    ps <- phsMomenta
    Sym'{..} <- sysJMJ s <$> phsTime <*> phsPositions
    let jmji = inv $ unSym sUpper
    return $ (quadForm ps jmji - quadForm sCross jmji + sCorner) / 2

-- | The Hamiltonian of a system (the sum of kinetic energy and the
-- potential energy), given the system's state in phase space.
hamiltonian
    :: (KnownNat m, KnownNat n)
    => System m n
    -> Phase n
    -> Double
hamiltonian s = do
    cfg@Cfg{..} <- fromPhase s
    let lg = lagrangian s cfg
    ps <- phsMomenta
    return $ (cfgVelocities <.> ps) - lg

-- | The "hamiltonian equations" for a given system at a given state in
-- phase space.  Returns the rate of change of the positions and
-- conjugate momenta, which can be used to progress the simulation through
-- time.
hamEqs
    :: (KnownNat m, KnownNat n)
    => System m n
    -> Phase n
    -> (R n, R n)
hamEqs s p@Phs{..} = (dHdp, -dHdq)
  where
    jmj  = sysJMJ s phsTime phsPositions
    jmj' = sysJMJ' s phsTime phsPositions
    dTdq = vec2r
         . flip fmap jmj' $ \(upper', cross', corner') ->
             let iUpper  = inv (unSym (sUpper jmj))
                 iUpper' = - iUpper <> unSym upper' <> iUpper
             in  phsMomenta <.> iUpper' #> phsMomenta / 2
               - phsMomenta <.> iUpper' #> sCross jmj
               - phsMomenta <.> iUpper  #> cross'
               - sCross jmj <.> iUpper' #> sCross jmj / 2
               - cross'     <.> iUpper  #> sCross jmj
               - corner' / 2
    dHdp = velocities s p
    dHdq = dTdq + snd (_sysPotentialGrad s phsTime phsPositions)

tr2
    :: (KnownNat m, KnownNat n, KnownNat o)
    => V.Vector m (L n o)
    -> V.Vector n (L m o)
tr2 = fmap (fromJust . (\rs -> withRows rs exactDims) . toList)
    . sequenceA
    . fmap (fromJust . V.fromList . toRows)

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
                       case (Proxy %<=? Proxy) :: (2 :<=? s) of
                         LE Refl -> (if l1 then tail else id)
                                  . toList
                                  $ evolveHam s p0 v
                         NLE Refl -> error "evolveHam': Internal error"
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
    toPs v = Phs (VS.unsafeHead t) pP pM
      where
        t : vs = LA.takesV [1,n,n] v
        Just [pP, pM] = traverse create vs

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

