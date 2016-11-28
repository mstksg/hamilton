Hamilton
========

[![Build Status](https://travis-ci.org/mstksg/hamilton.svg?branch=master)](https://travis-ci.org/mstksg/hamilton)

Simulate physics on arbitrary coordinate systems using [automatic
differentiation][ad] and [Hamiltonian mechanics][].  State only an arbitrary
parameterization of your system and a potential energy function!

[ad]: http://hackage.haskell.org/package/ad
[Hamiltonian mechanics]: https://en.wikipedia.org/wiki/Hamiltonian_mechanics

For example, a simulating a [double pendulum system][dps] by simulating the
progression of the angles of each bob:

[dps]: https://en.wikipedia.org/wiki/Double_pendulum

[![My name is William Rowan Hamilton](http://i.imgur.com/Vaaa2EC.gif)][gifv]

[gifv]: http://i.imgur.com/Vaaa2EC.gifv

You only need:

1.  Your generalized coordinates (in this case, `θ1` and `θ2`), and equations
    to convert them to cartesian coordinates of your objects:

    ~~~haskell
    x1 = sin θ1
    y1 = -cos θ1
    x2 = sin θ1 + sin θ2 / 2      -- second pendulum is half-length
    y2 = -cos θ1 - cos θ2 / 2
    ~~~

2.  The masses/inertias of each of those cartesian coordinates (`m1` for `x1`
    and `y1`, `m2` for `x2` and `y2`)

3.  A potential energy function for your objects:

    ~~~haskell
    U = (m1 y1 + m2 y2) * g
    ~~~

And that's it! Hamiltonian mechanics steps your generalized coordinates (`θ1`
and `θ2`) through time, without needing to do any simulation involving
`x1`/`y1`/`x2`/`y2`!  And you don't need to worry about tension or any other
stuff like that.  All you need is a description of your coordinate system
itself, and the potential energy!

~~~haskell
doublePendulum :: System 4 2
doublePendulum =
    mkSystem' (vec4 m1 m1 m2 m2)            -- masses
              (\(V2 θ1 θ2)     -> V4 (sin θ1)            (-cos θ1)
                                     (sin θ1 + sin θ2/2) (-cos θ1 - cos θ2/2)
              )                             -- coordinates
              (\(V4 _ y1 _ y2) -> (m1 * y1 + m2 * y2) * g)
                                            -- potential
~~~

Thanks to [~~Alexander~~ William Rowan Hamilton][WRH], we can express our
system parameterized by arbitrary coordinates and get back equations of motions
as first-order differential equations.  This library solves those first-order
differential equations for you using automatic differentiation and some matrix
manipulation.

[WRH]: https://www.youtube.com/watch?v=SZXHoWwBcDc

See a [blog post][] I wrote on this, and also the [hackage documentation][] and the
[example runner user guide][] (and its [source][example runner]).

[blog post]: https://blog.jle.im/entry/introducing-the-hamilton-library.html
[documentation]: http://hackage.haskell.org/package/hamilton
[example runner]: https://github.com/mstksg/hamilton/blob/master/app/Examples.hs
[user guide]: https://github.com/mstksg/hamilton#example-app-runner

### Full Exmaple

Let's turn our double pendulum (with the second pendulum half as long) into an
actual running program.  Let's say that `g = 5`, `m1 = 1`, and `m2 = 2`.

First, the system:

~~~haskell
import           Numeric.LinearAlgebra.Static
import qualified Data.Vector.Sized as V


doublePendulum :: System 4 2
doublePendulum = mkSystem' masses coordinates potential
  where
    masses :: R 4
    masses = vec4 1 1 2 2
    coordinates
        :: Floating a
        => V.Vector 2 a
        -> V.Vector 4 a
    coordinates (V2 θ1 θ2) = V4 (sin θ1)            (-cos θ1)
                                (sin θ1 + sin θ2/2) (-cos θ1 - cos θ2/2)
    potential
        :: Num a
        => V.Vector 4 a
        -> a
    potential (V4 _ y1 _ y2) = (y1 + 2 * y2) * 5


-- some helper patterns to pattern match on sized vectors
pattern V2 :: a -> a -> V.Vector 2 a
pattern V2 x y <- (V.toList->[x,y])
  where
    V2 x y = fromJust (V.fromList [x,y])

pattern V4 :: a -> a -> a -> a -> V.Vector 4 a
pattern V4 x y z a <- (V.toList->[x,y,z,a])
  where
    V4 x y z a = fromJust (V.fromList [x,y,z,a])
~~~

Neat!  Easy, right?

Okay, now let's run it.  Let's pick a starting configuration (state of the
system) of `θ1` and `θ2`:

~~~haskell
config0 :: Config 2
config0 = Cfg (vec2 1 0  )  -- initial positions
              (vec2 0 0.5)  -- initial velocities
~~~

Configurations are nice, but Hamiltonian dynamics is all about motion through
phase space, so let's convert this configuration-space representation of the
state into a phase-space representation of the state:

~~~haskell
phase0 :: Phase 2
phase0 = toPhase doublePendulum config0
~~~

And now we can ask for the state of our system at any amount of points in time!

~~~haskell
ghci> evolveHam doublePendulum phase0 [0,0.1 .. 1]
-- result: state of the system at times 0, 0.1, 0.2, 0.3 ... etc.
~~~

Or, if you want to run the system step-by-step:


~~~haskell
evolution :: [Phase 2]
evolution = iterate (stepHam 0.1 doublePendulum) phase0
~~~

And you can get the position of the coordinates as:

~~~haskell
positions :: [R 2]
positions = phsPositions <$> evolution
~~~

And the position in the underlying cartesian space as:

~~~haskell
positions' :: [R 4]
positions' = underlyingPos doublePendulum <$> positions
~~~

Example App runner
------------------

*([Source][example runner])*

Installation:

~~~bash
$ git clone https://github.com/mstksg/hamilton
$ cd hamilton
$ stack install
~~~

Usage:

~~~bash
$ hamilton-examples [EXAMPLE] (options)
$ hamilton-examples --help
$ hamilton-examples [EXAMPLE] --help
~~~

The example runner is a command line application that plots the progression of
several example system through time.


| Example      | Description                                                | Coordinates                                                         | Options                                                       |
|--------------|------------------------------------------------------------|---------------------------------------------------------------------|---------------------------------------------------------------|
| `doublepend` | Double pendulum, described above                           | `θ1`, `θ2` (angles of bobs)                                         | Masses of each bob                                            |
| `pend`       | Single pendulum                                            | `θ` (angle of bob)                                                  | Initial angle and velocity of bob                             |
| `room`       | Object bounding around walled room                         | `x`, `y`                                                            | Initial launch angle of object                                |
| `twobody`    | Two gravitationally attracted bodies, described below      | `r`, `θ` (distance between bodies, angle of rotation)               | Masses of bodies and initial angular veocity                  |
| `spring`     | Spring hanging from a block on a rail, holding up a weight | `r`, `x`, `θ` (position of block, spring compression, spring angle) | Masses of block, weight, spring constant, initial compression |
| `bezier`     | Bead sliding at constant velocity along bezier curve       | `t` (Bezier time parameter)                                         | Control points for arbitrary bezier curve                     |

Call with `--help` (or `[EXAMPLE] --help`) for more information.

More examples
-------------

### Two-body system under gravity

[![The two-body solution](http://i.imgur.com/TDEHTcb.gif)][gifv2]

[gifv2]: http://i.imgur.com/TDEHTcb.gifv

1.  The generalized coordinates are just:

    *   `r`, the distance between the two bodies
    *   `θ`, the current angle of rotation

    ~~~haskell
    x1 =  m2/(m1+m2) * r * sin θ        -- assuming (0,0) is the center of mass
    y1 =  m2/(m1+m2) * r * cos θ
    x2 = -m1/(m1+m2) * r * sin θ
    y2 = -m1/(m1+m2) * r * cos θ
    ~~~

2.  The masses/inertias are again `m1` for `x1` and `y1`, and `m2` for `x2` and
    `y2`

3.  The potential energy function is the classic gravitational potential:

    ~~~haskell
    U = - m1 * m2 / r
    ~~~

And...that's all you need!

Here is the actual code for the two-body system, assuming `m1` is `100` and
`m2` is `1`:

~~~haskell
twoBody :: System 4 2
twoBody = mkSystem masses coordinates potential
  where
    masses :: R 4
    masses = vec4 100 100 1 1
    coordinates
        :: Floating a
        => V.Vector 2 a
        -> V.Vector 4 a
    coordinates (V2 r θ) = V4 (r1 * cos θ) (r1 * sin θ)
                              (r2 * cos θ) (r2 * sin θ)
      where
        r1 =   r *   1 / 101
        r2 = - r * 100 / 101
    potential
        :: Num a
        => V.Vector 4 a
        -> a
    potential (V2 r _) = - 100 / r
~~~

Potential improvements
----------------------

*   **Time-dependent systems**:  Shouldn't be an problem in theory/math; just
    add a time parameter before all of the functions.  This opens a lot of
    doors, like deriving inertial forces for free (like the famous Coriolis
    force and centrifugal force).

    The only thing is that it makes the API pretty inconvenient, because it'd
    require all of the functions to also take a time parameter.  Of course, the
    easy way out/ugly solution would be to just offer two versions of the same
    function (one for time-independent systems and one for time-dependent
    systems.  But this is un-ideal.

*   Velocity-dependent potentials:  Would give us the ability to model
    dissipative systems, like systems with friction (dependent on `signum v`)
    and linear & quadratic wind resistance.

    This issue is much harder, theoretically.  It involves inverting arbitrary
    functions `forall a. RealFloat a => V.Vector n a -> V.Vector m a`.  It
    might be possible with the help of some
    [bidirectionalization techniques][bff-pearl], but I can't get the [bff][]
    package to compile, and I'm not sure how to get [bff-mono][] to work with
    numeric functions.

    If anyone is familiar with bidirectionalization techniques and is willing
    to help out, please send me a message or open an issue! :)

[bff-pearl]: https://pdfs.semanticscholar.org/5f0d/ef02dbd96e102be9104d2ceb728d2a2a5beb.pdf
[bff]: http://hackage.haskell.org/package/bff
[bff-mono]: http://hackage.haskell.org/package/bff-mono

