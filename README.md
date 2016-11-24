Hamilton
========

Simulate physics on arbitrary coordinate systems using automatic
differentiation and Hamiltonian mechanics.

For example, a simulating a double pendulum system by simulating the
progression of the angles of each bob:

[![My name is William Rowan Hamilton](http://i.imgur.com/Vaaa2EC.gif)](http://i.imgur.com/Vaaa2EC.gifv)

You only need:

1.  Your generalized coordinates (in this case, `θ1` and `θ2`), and equations
    to convert them to cartesian coordinates of your objects:

    ~~~haskell
    x1 = sin θ1
    y1 = -cos θ1
    x2 = sin θ1 + sin θ2
    y2 = -cos θ1 - cos θ2
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


See [documentation](https://mstksg.github.io/hamilton/).

