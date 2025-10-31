# Chapter 1: Introduction

<!-- ======================= -->
<!-- PROBLEM 1.1             -->
<!-- ======================= -->
## Problem 1.1
Consider a constant density background $\rho$, a planet a distance $r$ from the sun, feels a gravitational force from the background mass inside a sphere of radius $r$ given by the solution of

$$
\nabla\Phi_\rho = 4\pi G \rho = \frac{1}{r^2}\frac{d}{dr}\left(r^2\frac{d\Phi_\rho}{dr}\right),
$$

which is

$$
\Phi_\rho = \frac{2\pi G \rho}{3} r^2 = \frac{1}{2}\omega_0^2 r^2 ~~\text{ where }~~ \omega_0^2 = \frac{4\pi G \rho}{3}.
$$

A planet orbiting the sun the feels both the force from the sun and the background density, that is

$$
\Phi(r) = -\frac{GM}{r} + \frac{1}{2}\omega_0^2 r^2.
$$

The circular-orbit (azimuthal) frequency and the epicyclic frequency are given by

$$
\begin{align}
\Omega_\phi^2 &= \frac{1}{r}\frac{d\Phi}{dr} = \frac{GM}{r^3} + \omega_0^2\\
\Omega_r^2 &= R\frac{d\Omega_\phi^2}{dR} + 4\Omega_\phi^2.
\end{align}
$$

At the perihelion $a$, these frequencies are


$$
\Omega_\phi^2(a) = \frac{GM}{a^3} + \omega_0^2 = n^2 + \omega_0^2,
\quad
\Omega_r^2(a) = \frac{GM}{a^3} + 4\omega_0^2 = n^2 + 4\omega_0^2,
$$

where $n^2 = GM/a^3$ is the mean motion of the planet. The rate of perihelion precession is then (for small $\omega_0^2/n^2$, i.e., small background density)

$$
\begin{align}
\delta\varpi &= \Omega_\phi - \Omega_r = \sqrt{n^2 + \omega_0^2} - \sqrt{n^2 + 4\omega_0^2} \\
&\approx n\left(1 + \frac{\omega_0^2}{2n^2}\right) - n\left(1 + \frac{4\omega_0^2}{2n^2}\right) \\
&= -\frac{3\omega_0^2}{2n} = -\frac{2\pi G \rho}{n}.
\end{align}
$$

Plugging in numbers

$$
\delta\varpi \approx -9.279 \times 10^{-10} \frac{\rmarcsec}{\rmyr} \left(\frac{\rho}{\Msun/\rmpc^3}\right)\left(\frac{a}{\rmAU}\right)^{3/2},
$$

!!! warning "Potential error"
    For Neptune, and $\rho = 0.1 \Msun/\rmpc^3$, we get $\delta\varpi = -1.5\times 10^{-8}\rmarcsec/\rmyr$ this differs from the value given in Binney & Tremaine (2008) by six orders of magnitude; unfortunately I cannot come up with a better derivation that can explain the discrepancy.

<!-- ======================= -->
<!-- PROBLEM 1.2             -->
<!-- ======================= -->
## Problem 1.2

### a. luminosity density to surface brightness
Using the same construct of Fig. 2.3 define $R, z$ such that $r^2 = R^2 + z^2$. The surface brightness at projected radius $R$ is given by

$$
\begin{align}
I(R) &= \int_{-\infty}^{+\infty} dz j(\sqrt{R^2 + z^2}) \\
&= 2\int_{0}^{+\infty} dz j(\sqrt{R^2 + z^2}) \\
&= 2\int_R^\infty dr\frac{r j(r)}{\sqrt{r^2 - R^2}}, \quad\text{with}\quad R = r\cos\theta \\
&= 2R\int_0^{\pi/2} d\theta j(\sec\theta)\sec^2\theta
\end{align}
$$

### b. Plummer sphere
For a Plummer sphere

$$
j(r) = \frac{j_0}{(1 + r^2/b^2)^{5/2}},
$$

we have

$$
\begin{align}
I(R) &= 2R j_0 \int_0^{\pi/2} d\theta \frac{\sec^2\theta}{(1 + R^2\sec^2\theta/b^2)^{5/2}} \quad\text{with}\quad u = \tan\theta \\
&= 2R j_0 \int_0^{\infty} du \frac{1}{(1 + R^2(1 + u^2)/b^2)^{5/2}} \\
&= 2R j_0 \int_0^{\infty} du \frac{1}{(A + B u^2)^{5/2}}
\end{align}
$$

with $A = 1 + R^2/b^2$ and $B = R^2/b^2$. The integral can be solved with

$$
\int_0^{\infty} du \frac{1}{(A + B u^2)^{5/2}} = \frac{2}{3A^2 B^{1/2}}
$$

Putting everything together we have

$$
I(R) = \frac{4}{3} \frac{j_0 b}{(1 + R^2/b^2)^2}.
$$

### c. Inversion formula

Using Eq. (B.72a) and Eq. (B.72b) with

$$
x = R^2, \quad t = r^2, \quad f(x) = I(\sqrt{x}), \quad g(t) = j(\sqrt{t}) \text{ and } \alpha = 1/2
$$

we have

$$
\begin{align}
f(R^2) &= \int_{R^2}^\infty dt \frac{g(t)}{\sqrt{t - R^2}} =
\int_R^\infty dr \frac{2r g(r^2)}{\sqrt{r^2 - R^2}} \\
&= 2\int_R^\infty dr \frac{r j(r)}{\sqrt{r^2 - R^2}} = I(R)
\end{align}
$$

And for the inverse we start with

$$
f'(x) = \frac{d}{dx}I(\sqrt{x}) = \left.\frac{1}{2\sqrt{x}}I'(\sqrt{x})\right|_{R=\sqrt{x}} = \frac{1}{2R}I'(R)
$$

so that

$$
j(r) = g(r^2) = -\frac{1}{\pi}\int_{r^2}^\infty dx \frac{f'(x)}{\sqrt{x - r^2}} = -\frac{1}{\pi}\int_r^\infty dR \frac{1}{\sqrt{R^2 - r^2}}\frac{dI}{dR}.
$$

### d. de Vaucouleurs profile

For the de Vaucouleurs profile we have ($m=4$)

$$
I_m(R) = I_e \exp\left[-b_m\left(\left(\frac{R}{R_e}\right)^{1/m} - 1\right)\right], \quad b_m \approx 2m - 0.324
$$

A numerical implementation of the inversion formula is implemented in the solutions module

```python
from galactic_dynamics_bt.chapter01.sersic_profile import plot_sersic_profile

plot_sersic_profile()
```


![FRW Model Phase Diagram](assets/generated/sersic_profile.png)

*Figure P1.2: Sersic luminosity density profile derived from the surface brightness profile using the inversion formula. For comparison different values of the parameter $m$ are shown.*

<!-- ======================= -->
<!-- PROBLEM 1.11            -->
<!-- ======================= -->
## Problem 1.11

From Eq.~(1.50) we know that

$$
\dot{a}^2 - \frac{8\pi G \rho}{3}a^2 = 2E.
$$

At present day $a_0 = 1$ and this becomes (Neglecting any contribution from radiation)

$$
2E = H_0^2(1 - \Omega_0) = (1 - \Omega_{\Lambda 0} - \Omega_{m0})H_0^2 = -\frac{kc^2}{x_u^2}
$$

The boundary between a bound and unbound universe is then given by $k=0$ or $\Omega_{\Lambda 0} + \Omega_{m0} = 1$.

Now, Eq.~(1.49) can be rewritten as

$$
\frac{\ddot{a}}{a} = \frac{H_0^2}{2}(2\Omega_{\Lambda 0} - \Omega_{m0}a^{-3}).
$$

The transition between a decelerating and accelerating universe occurs when $\ddot{a} = 0$, for all values of $a$, in particular at present day $a_0 = 1$, this happens when

$$
2\Omega_{\Lambda 0} - \Omega_{m0} = 0
$$

![FRW Model Phase Diagram](assets/generated/frw_model.png)

*Figure P1.11: Phase diagram of FRW cosmological models in the $\Omega_{m0}$-$\Omega_{\Lambda 0}$ plane. The solid line indicates the boundary between bound and unbound models, while the dashed line indicates the boundary between decelerating and accelerating models.*
