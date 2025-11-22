# Surfaces of Section


## Integrating orbits with Symplectic Integrators
For the potential

$$
\Phi(R, z) = \frac{1}{2}v_0^2 \ln\left(R^2 + \frac{z^2}{q^2}\right),
$$

we can numerically integrate orbits using a symplectic integrator; given the hamiltonian nature of the sytem, this represents the more accurate way to integrate orbits over long timescales [@Candy1991].

```python
>>> from galactic_dynamics_bt.chapter03.logarithmic_potential import (
... integrate_orbit,
... LogarithmicPotential,
... LogarithmicPotentialParams,
... )
>>> model = LogarithmicPotential(LogarithmicPotentialParams(q=0.9, v0=1.0), Lz=0.2)
>>> t, R, z, pR, pz = integrate_orbit(model, E=-0.8, R0=0.35, pR0=0.1, t_max=20.0, dt=0.01)
```

This snippet will generate an initial condition with energy $E = -0.8$, angular momentum $L_z = 0.2$, and integrate the orbit for 20 time units with a timestep of $\Delta t= 0.01$. Behind the scenes, the integrator is using a second-order scheme to perform the numerical evolution of the orbits.

The top two panels of *Figure 1* show the resulting orbit in the meridional plane $(R, z)$ using this strategy for $q = 0.9$ (left) and $q = 0.6$ (right).

For the same initial conditions on the top left panel, I also calculate the total angular momentum

$$
\begin{align}
L^2 &= \mathbf{L} \cdot \mathbf{L} \\
&= (z p_R - R p_z)^2 + (R^2 + z^2)\frac{L_z^2}{R^2},
\end{align}
$$

which is displayed in *Figure 2*, in two different windows of time. As expected, the total angular momentum is not exactly conserved due to the flattening of the potential, but it oscillates around a mean value.

![Orbits in a flatten logarithmic potential](assets/generated/orbital_properties.png)

*Figure: 1: Top panels: Orbits in the meridional plane $(R, z)$ for a flatten logarithmic potential with $q = 0.9$ (left) and $q = 0.6$ (right). Bottom panels: Corresponding Poincaré surfaces of section. The solid line indicates the zero-velocity curve, while the dashed line shows the orbit if the total angular momentum were exactly conserved.*

## Poincaré Surface of Section

To calculate the Poincaré surface of section, I record the coordinates along the orbit and look for the points where it crosses the plane $z = 0$: $z(t) < 0 < z(t + \Delta t)$ and $p_z(t) > 0$. The points are then subsequently integrated by changing the coordinates frp $(R(t), z(t), p_R(t), p_z(t))$ to $(R(z), t(z), p_R(z), p_z(z))$ [@Cvitanovic2020],

$$
\begin{align}
\frac{dt}{dz} &= \frac{1}{dz/dt} \\
\frac{dR}{dz} &= \frac{dR/dt}{dz/dt} \\
\frac{dp_R}{dz} &= \frac{dp_R/dt}{dz/dt} \\
\frac{dp_z}{dz} &= \frac{dp_z/dt}{dz/dt}.
\end{align}
$$

These equations are then integrated to $z = 0$ using a simple Runge-Kutta scheme. The resulting surface of section is shown in the bottom panels of *Figure 1* for $q = 0.9$ (left) and $q = 0.6$ (right) (dots).

To each surface of section I also add:

1. The zero-velocity curve (solid line), which is obtained by setting $p_z = 0$ and solving for $p_R(R)$ given the energy constraint.
2. The orbit if the total angular momentum were exactly conserved (dashed line). This is obtained by solving for $p_R(R)$ as above, but swapping $L_z$ with $L$.

![Total angular momentum](assets/generated/angular_momentum.png)

*Figure: 2: Total angular momentum for the orbit shown in the top left panel of Figure 1, calculated in two different time windows.*

<!-- ======================= -->
<!-- REFERENCES.             -->
<!-- ======================= -->

## References
\bibliography
