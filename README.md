Hamilton
========

Simulate physics on arbitrary coordinate systems using [automatic
differentiation][ad] and [Hamiltonian mechanics][].

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

And Hamiltonian mechanics steps your generalized coordinates (`θ1` and `θ2`)
through time, without needing to do any simulation involving
`x1`/`y1`/`x2`/`y2`!  And you don't need to worry about tension or any other
stuff like that.  All you need is a description of your coordinate system
itself, and the potential energy!

Thanks to [~~Alexander~~ William Rowan Hamilton][WRH], we can express our
system parameterized by arbitrary coordinates and get back equations of motions
as first-order differential equations.  This library solves those first-order
differential equations for you using automatic differentiation and some matrix
manipulation.

[WRH]: https://www.youtube.com/watch?v=SZXHoWwBcDc

See [documentation][] and [example runner][].

[documentation]: https://mstksg.github.io/hamilton/
[example runner]: https://github.com/mstksg/hamilton/blob/master/app/Examples.hs

Example runner
--------------

Usage:

~~~bash
$ hamilton-examples [EXAMPLE] (options)
$ hamilton-examples --help
$ hamilton-examples [EXAMPLE] --help
~~~

The example runner is a command line application that plots the progression of
several example system through time.

| Example      | Description                                          | Coordinates                                           | Options                                      |
|--------------|------------------------------------------------------|-------------------------------------------------------|----------------------------------------------|
| `doublepend` | Double pendulum, described above                     | `θ1`, `θ2` (angles of bobs)                           | Masses of each bob                           |
| `pend`       | Single pendulum                                      | `θ` (angle of bob)                                    | Initial angle and velocity of bob            |
| `room`       | Object bounding around walled room                   | `x`, `y`                                              | Initial launch angle of object               |
| `twobody`    | Two gravitationally attracted bodies                 | `r`, `θ` (distance between bodies, angle of rotation) | Masses of bodies and initial angular veocity |
| `bezier`     | Bead sliding at constant velocity along bezier curve | `t` (Bezier time parameter)                           | Control points for arbitrary bezier curve    |

Call with `--help` (or `[EXAMPLE] --help`) for more information.
