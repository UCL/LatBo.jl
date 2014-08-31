Lattice-Boltzmann equations
===========================

Let $f_i(\vec{r}, t)$ be the probability of finding a particle with velocity
$\vec{c}_i$ at position $\vec{r}$ on the lattice at time $t$. Then
$f_i(\vec{r}, t + \Delta t)$ is given by:

$$
    f_i(\vec{r}, t+\Delta t) = f_i(\vec{r}, t) +
        \frac{1}{\tau}\left[
            f_i^{(eq)}(\vec{r}, t) - f_i(\vec{r}, t)
        \right]
$$

Where $\tau\Delta t$ is a relaxation time and $f_i^{(eq)}$ is the local
equilibrium population obtained at position $\vec{r}$ and time $t$:

$$
    f_i^{(eq)} = \omega_i\rho\left[
        1 + \frac{1}{c_s}\vec{u}\cdot\vec{c_i} 
        + \frac{1}{2c_s^2}
            \vec{u}\cdot\bar{\bar{Q_i}}\cdot\vec{u}
    \right]
$$
$$
    c_s^2 = \sum_i \omega_i c_i^2
$$
$$
    \vec{u} = \sum_i \frac{f_i \vec{c}_i}{\rho}
$$
$$
    \rho = \sum_i f_i
$$
$$
    Q_{i, \alpha\beta} = c_{i, \alpha}c_{i, \beta} - c_s^2\delta_{\alpha\beta}
$$

The $\omega_i$ are weights that depend on the discretization scheme used.
We will likely look at D2Q9 for two dimensional calculations, D3Q13 for three
dimensional calculations.


D2Q9
----

This a two dimensional discretization with on-site, nearest neighbors and next
nearest neighbors interactions on a *square* lattice. The weights are
defined as:

- on site: $w_0 = \frac{4}{9}$, null velocity
- nearest-neighbor: $w_{1-4}=\frac{1}{9}$,
    $c_i = (0, 0, 1)\frac{\Delta x}{\Delta t}$
- next nearest-neighbor: $w_{5-8} = \frac{1}{36}$,
    $c_i = (0, 1, 1)\frac{\Delta x}{\Delta t}$

The velocities are such that it takes time $t$ to go from $i$ to the relevant
neighboring site.


D3Q19
-----

Similarly, this is a three-dimensional discretization with on-site, nearest
neighbors and next nearest neighbors interactions on a *cubic* lattice.

- on site: $w_0 = \frac{1}{3}$, null velocity
- nearest-neighbor: $w_{1-6}=\frac{1}{18}$, velocities alongst lattice vectors
- next nearest-neighbor: $w_{7-18} = \frac{1}{18}$, velocities alongst short
  diagonals


Thermodynamic Quanitities
=========================

$$
    \vec{u} = \sum_i \frac{f_i \vec{c}_i}{\rho}
$$
$$
    \rho = \sum_i f_i
$$
$$
    \bar{\bar{\sigma}} = \left(1-\frac{1}{2\tau}\right)
        \sum_i (f_i - f_i^{eq})\vec{c}_i\vec{c}_i
$$

Where $\bar{\bar{\sigma}}$ is the deviatoric stress tensor.


Boundary Conditions
===================

Inlets and outlets
------------------

Bounce-back
-----------

Bounce-back is an implementation of the no-slip boundary condition, meaning
that at the wall particles have the same velocity as the wall (zero in our
case).


Initial Conditions
==================
