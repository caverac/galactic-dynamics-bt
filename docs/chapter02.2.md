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

<!-- ======================= -->
<!-- PROBLEM 2.14            -->
<!-- ======================= -->
## Problem 2.14

For the density profile

$$
\rho(m) = \frac{\rho_0}{(1 + m^2/a_1^2)^2},
$$

we have

$$
\psi(\infty) - \psi(m) = \int_{m^2}^\infty du \, \rho(u) = \rho_0 a_1^2  \frac{1}{1 + m^2/a_1^2}.
$$

We then need to compute

$$
\begin{align}
m(\mathbf{x}, t) & = a_1^2\sum_i\frac{x_i^2}{a_i^2 + \tau} \\
&= a_1^2 \frac{R^2}{a_1^2 + \tau} + a_1^2 \frac{z^2}{a_3^2 + \tau} \\
&= \frac{a_1^2}{a_1^2 + \tau}\Delta^2\cosh^2 u \sin^2 v + \frac{a_1^2}{a_3^2 + \tau}\Delta^2\sinh^2 u \cos^2 v \\
&= \frac{a_1^2 \Delta^2}{a_1^2 + \tau}(1 + \sinh^2 u - \cos^2 v - \sinh^2 u \cos^2 v)  + \frac{a_1^2 \Delta^2 }{a_3^2 + \tau}(\sinh^2 u \cos^2 v) \\
&= \frac{a_1^2}{a_1^2 + \tau}\left(\Delta^2 + \lambda + \mu + \frac{\lambda \mu}{\Delta^2}\right) - \frac{a_1^2}{a_3^2 + \tau}\frac{\lambda \mu}{\Delta^2} \\
&= \frac{a_1^2}{a_1^2 + \tau}(\Delta^2 + \lambda + \mu) + \frac{a_1^2 \lambda \mu}{(a_1^2 + \tau)(a_3^2 + \tau)} \\
&= \frac{a_1^2}{(a_1^2 + \tau)(a_3^2 + \tau)}\left[(a_3^2 + \tau)(\Delta^2 + \lambda + \mu) + \lambda \mu\right]
\end{align}
$$

Which means

$$
\begin{align}
1 + \frac{m^2}{a_1^2} &= \frac{1}{(a_1^2 + \tau)(a_3^2 + \tau)}\left[(a_1^2 + \tau)(a_3^2 + \tau) + (a_3^2 + \tau)(a_3^2 - a_1^2 + \lambda + \mu) + \lambda\mu\right] \\
&= \frac{(\tau + a_3^2 + \lambda)(\tau + a_3^2 + \mu)}{(a_1^2 + \tau)(a_3^2 + \tau)}.
\end{align}
$$

Thus,

$$
\begin{align}
\Phi(\mathbf{x}) &= -\pi G a_1 a_2 a_3 \int_0^\infty d\tau \frac{\psi(\infty) - \psi(m)}{\sqrt{(\tau + a_1^2)(\tau + a_1^2)(\tau + a_3^2)}} \\
&= -\pi G a_3 a_1^2\rho_0 \int_0^\infty d\tau \, \frac{(a_3^2 + \tau)^{1/2}}{(\tau + a_3^2 + \lambda)(\tau + a_3^2 + \mu)} \\
&= \frac{\pi G \rho_0 a_1^2a_3}{\lambda - \mu} \int_0^\infty d\tau \, (a_3^2 + \tau)^{1/2} \left(\frac{1}{\tau + a_3^2 + \mu} - \frac{1}{\tau + a_3^2 + \lambda}\right) \\
&= \frac{H(\lambda) - H(\mu)}{\lambda - \mu},
\end{align}
$$

with

$$
H(x) = \pi G \rho_0 a_1^2 a_3 \int_0^\infty d\tau \, \frac{(a_3^2 + \tau)^{1/2}}{\tau + a_3^2 + x}.
$$

<!-- ======================= -->
<!-- PROBLEM 2.15            -->
<!-- ======================= -->
## Problem 2.15

Start from

$$
\Sigma(R, \phi) = -\frac{1}{2\pi G} \sum_{m= -\infty}^\infty \int_0^\infty dk\, k S_m(k) J_m(kR) e^{im\phi}.
$$

Multiplying by $e^{-im'\phi}$ and averaging over $\phi$ gives

$$
\frac{1}{2\pi}\int_0^{2\pi} d\phi \, \Sigma(R, \phi) e^{-im'\phi} = -\frac{1}{2\pi G} \int_0^\infty dk \, kS_{m'}(k) J_{m'}(kR),
$$

where we used the orthogonality of the complex exponentials $\int_0^{2\pi} e^{i(m-m')\phi} d\phi = 2\pi \delta_{mm'}$. Similarly, multiplying by $R J_{m'}(k'R)$ and integrating over $R$ gives

$$
S_m(k) =-G \int_0^{2\pi}d\phi \, e^{-im\phi} \int_0^\infty dR \, R  \Sigma(R, \phi) J_m(kR),
$$

where we used the orthogonality of Bessel functions $\int_0^\infty dR \, R J_m(kR) J_m(k'R) = \delta(k-k')/k$. The potential can be written as

$$
\Phi(R, \phi, z) = \sum_{m=-\infty}^\infty \int_0^\infty dk \, e^{-k|z| + im\phi} J_m(kR) S_m(k).
$$

For an axisymmetric disk, only the $m=0$ term contributes,

$$
S_m(k) = -2\pi G \delta_{m0} \int_0^\infty dR \, R \, \Sigma(R) J_0(kR),
$$

and the potential at the origin ($R=0$, $z=0$) is

$$
\begin{align}
\Phi(0, 0) &= \int_0^\infty dx \, J_0(0)e^0 S_0(k) \\
&= -2\pi G \int_0^\infty dk \int_0^\infty dR \, R \, \Sigma(R) J_0(kR) \\
&= -2\pi G \int_0^\infty dR \, \Sigma(R) \int_0^\infty d(kR) \, J_0(kR) \\
&= -2\pi G \int_0^\infty dR \, \Sigma(R)
\end{align}
$$

<!-- ======================= -->
<!-- PROBLEM 2.16            -->
<!-- ======================= -->
## Problem 2.16

**Spherical System**

$$
\begin{align}
\frac{d\Phi(r)}{dr} &= -4\pi G\frac{d}{dr}\left[ \frac{1}{r}\int_0^r dr'\, r'^2\rho(r') + \int_r^\infty dr' \, r' \rho(r') \right] \\
&= -4\pi G \left[-\frac{1}{r^2}\int_0^r dr' \, r'^2 \rho(r') + r \rho(r) - r\rho(r) \right] \\
&= \frac{GM(r)}{r^2} \ge 0.
\end{align}
$$

The potential is thus a non-decreasing function of $r$.

**Axisymmetric System**

Consider the density profile

$$
\Sigma(R) = \frac{M}{2\pi a}\delta(R - a).
$$


We have

$$
S_0(k) = -\frac{GM}{a}\int_0^\infty dR \, R \, \delta(R - a) J_0(kR) = -GM J_0(ka).
$$

And

$$
\Phi(R, 0) = -GM \int_0^\infty dk \, J_0(kR) J_0(ka) = -\frac{2GM}{a\pi}
\begin{cases}K\left(R^2/a^2\right) & R < a \\ aK\left(a^2/R^2\right)/R & R > a \end{cases},
$$

which decreases from $R=0$ to $R=a$ and then increases for $R>a$, that is, the potential is not a monotonic function of $R$.

<!-- ======================= -->
<!-- PROBLEM 2.17            -->
<!-- ======================= -->
## Problem 2.17

With $R^2 = X^2 + z^2$

$$
\begin{align}
\mu(X) &= 2\int_0^\infty dz \, \Sigma(\sqrt{X^2 + z^2}) \\
&= 2 \int_X^\infty dR\, \frac{R}{\sqrt{R^2 - X^2}} \Sigma(R).
\end{align}
$$

Using Abel's inversion formula

$$
\Sigma(R) = -\frac{1}{\pi} \int_R^\infty dX \, \frac{1}{\sqrt{X^2 - R^2}}\frac{d\mu}{dX}.
$$

Therefore

$$
\begin{align}
\Phi(R, z) &= 4G\int_0^\infty da\, \sin^{-1}\left(\frac{2a}{\sqrt{+} + \sqrt{-}} \right)  \frac{d}{da} \int_a^\infty dR' \, \frac{R'}{\sqrt{R'^2 - a^2}} \Sigma(R') \\
&= 2G \int_0^\infty da \, \sin^{-1}\left(\frac{2a}{\sqrt{+} + \sqrt{-}} \right) \frac{d\mu}{da}, \quad\text{with } \sqrt{\pm} = \sqrt{z^2 + (R \pm a)^2}.
\end{align}
$$

<!-- ======================= -->
<!-- PROBLEM 2.18            -->
<!-- ======================= -->
## Problem 2.18

$$
\begin{align}
\Sigma(R) &= \frac{1}{2\pi G} \int_0^\infty dk \, k J_0(kR) \int_0^\infty dR' \, v_c^2(R') J_1(kR') \\
&= \frac{v_0^2}{2\pi G} \int_0^\infty dk \, k J_0(kR) \int_0^\infty dR' \, \frac{J_1(kR')}{1 + R'^2/a^2} \\
&= \frac{v_0^2}{2\pi G} \int_0^\infty dk \, J_0(k R) \left(1 - e^{-ak} \right) \\
&= \frac{v_0^2}{2\pi G R}\left[1 - \frac{1}{\sqrt{1 + a^2 / R^2}} \right].
\end{align}
$$

<!-- ======================= -->
<!-- PROBLEM 2.19            -->
<!-- ======================= -->
## Problem 2.19


## References
\bibliography
