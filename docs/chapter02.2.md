<!-- ======================= -->
<!-- PROBLEM 2.11            -->
<!-- ======================= -->
## Problem 2.11

For an exponential disk with density profile $\Sigma(R, z) = \Sigma_0 e^{-R/R_d} \delta(z)$, the potential in the equatorial plane is

$$
\Phi(R, 0) = -\pi G \Sigma_0 R \left[ I_0\left(\frac{R}{2 R_d}\right) K_1\left(\frac{R}{2 R_d}\right) - I_1\left(\frac{R}{2 R_d}\right) K_0\left(\frac{R}{2 R_d}\right) \right].
$$

The potential energy is

$$
W = -\frac{1}{2} \int d^3\mathbf{x} \, \rho(\mathbf{x}) \Phi(\mathbf{x}) =
\pi \int_0^\infty dR \, R \, \Sigma(R) \Phi(R, 0) \approx -11.6274 \, G \Sigma_0^2 R_d^3.
$$

```python
>>> from galactic_dynamics_bt.chapter02.exponential_disk import potential_energy
>>> potential_energy()
-11.62735...
```

The angular momentum of a mass element $dm = d^2\mathbf{x} \, \rho(\mathbf{x})$ is $dJ = dm R v_c(R)$ -assuming circular orbits-

$$
J = 2\pi \int_0^\infty dR \, R^2 \, \Sigma(R) v_c(R) \approx 17.4654 \sqrt{G \Sigma_0^3 R_d^7},
$$

```python
>>> from galactic_dynamics_bt.chapter02.exponential_disk import angular_momentum
>>> angular_momentum()
17.46538...
```

Similarly for the kinetic energy

$$
K = \pi \int_0^\infty dR \, R \, \Sigma(R) v_c^2(R) \approx 5.8134 \, G \Sigma_0^2 R_d^3.
$$

```python
>>> from galactic_dynamics_bt.chapter02.exponential_disk import kinetic_energy
>>> kinetic_energy()
5.81367...
```

<!-- ======================= -->
<!-- PROBLEM 2.12            -->
<!-- ======================= -->
## Problem 2.12

$$
\Phi(f) = -\frac{GM}{r}e^{-\alpha r},
$$

$$
\frac{d\Phi}{dr} = GM \frac{e^{-\alpha r}}{r^2} (1 + \alpha r),
$$

and

$$
\frac{d}{dr}\left(r^2 \frac{d\Phi}{dr}\right) = -GM \alpha^2 r e^{-\alpha r}.
$$

Which means

$$
\nabla^2 \Phi = \frac{1}{r^2} \frac{d}{dr}\left(r^2 \frac{d\Phi}{dr}\right) = -GM \alpha^2 \frac{e^{-\alpha r}}{r}.
$$

<!-- ======================= -->
<!-- PROBLEM 2.13            -->
<!-- ======================= -->
## Problem 2.13

Any uniform density confocal spheroid can be written as a stack of infinitesimal confocal shells. By using the Homeoid Theorem, the exterior isopotentials of a thin confocal homoeoid are confocal with the shell, and the potential is constant on every exterior confocal spheroid. Hence, for any fixed exterior confocal surface $\Sigma_*$, the contribution from a given thin shell is a single constant over all of $\Sigma_*$. Summing over all shells, the total potential on $\Sigma_*$ is constant. And because the Newtonian potential is linear in density, that constant is proportional to the shell's mass.

Summing constant contributions over all shells shows that the potential on $\Sigma_*$ produced by the whole body is a constant that depends only on the total mass enclosed, not on how that mass is distributed among the interior shells.

Now, let's consider the uniqueness of the solution on the exterior region. Pick any $\Sigma_*$ enclosing both spheroids. Outside $\Sigma_*$, the potential satisfies Laplace's equation and decays as $-GM/r$ with the same $M$ (equal masses). Since both bodies induce the same constant Dirichlet value on $\Sigma_*$, the exterior solutions are identical by the uniqueness theorem. This holds for every choice of enclosing $\Sigma_*$, so the two exterior potentials (and fields) coincide everywhere outside the mass.

## References
\bibliography
