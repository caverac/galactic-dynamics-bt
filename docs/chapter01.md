# Chapter 1: Introduction

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

*Figure 1.1: Phase diagram of FRW cosmological models in the $\Omega_{m0}$-$\Omega_{\Lambda 0}$ plane. The solid line indicates the boundary between bound and unbound models, while the dashed line indicates the boundary between decelerating and accelerating models.*