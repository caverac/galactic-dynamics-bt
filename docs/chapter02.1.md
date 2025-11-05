<!-- ======================= -->
<!-- PROBLEM 2.1             -->
<!-- ======================= -->
## Problem 2.1

$$
\begin{align}
W &= -4\pi G\int_0^\infty dr \, r\rho(r)M(r) \\
&= -G \int_0^\infty dr \, \frac{1}{r}\left(M \frac{dM}{dr}\right) \\
&= -\frac{G}{2}\int_0^\infty dr \, \frac{d}{dr}\left(\frac{M^2}{r}\right)  - \frac{G}{2}\int_0^\infty dr \, \frac{M^2}{r^2}
\end{align}
$$

Since $M$ is finite, the term $M^2/r$ vanishes as $r\to\infty$. Now, for $r\to 0$, we have that $M^2/r \sim r^5 \rho^2$ which implies that this term also vanishes at $r=0$ provided that $\rho(r)$ diverges more slowly than $r^{-5/2}$. Thus, the first term in the last line above is zero and we obtain the desired result:

$$
W = -\frac{G}{2}\int_0^\infty dr \, \frac{M^2(r)}{r^2}
$$

<!-- ======================= -->
<!-- PROBLEM 2.2             -->
<!-- ======================= -->
## Problem 2.2

For a spherical potential

$$
\begin{align}
W_{jk} &= - \int d^3\mathbf{x} \, \rho(r) x_j \frac{\partial \Phi}{\partial x_k} \\
&= - \int d^3\mathbf{x} \, \rho(r) r \left(\frac{d\Phi}{dr} \right)\frac{x_j x_k}{r^2}  \\
&= -\int d^3\mathbf{x} \, \rho(r) r \left(\frac{GM(r)}{r^2} \right)\frac{x_j x_k}{r^2} \\
&= -G\int_0^\infty dr \, r \rho(r) M(r) \int d\Omega \, \frac{x_j x_k}{r^2} \\
&= -\frac{4\pi}{3} G \delta_{jk} \int_0^\infty dr \, r \rho(r) M(r) \\
&= \frac{1}{3} W \delta_{jk}
\end{align}
$$

One way to calculate the angular integral is to use Schur's lemma. Or we can represent the Cartesian coordinates in terms of spherical coordinates and perform the integral explicitly.

<!-- ======================= -->
<!-- PROBLEM 2.3             -->
<!-- ======================= -->
## Problem 2.3

### Using Gauss's theorem

Draw a cylindrical surface of radius $R$ and height $h$ aligned with the $z$-axis,

$$
2\frac{d\Phi}{dz} = 4\pi G \Sigma (\pi R^2) \implies \frac{d\Phi}{dz} = 2\pi G \Sigma,
$$

thus,

$$
\Phi(z) = 2\pi G \Sigma |z| + \mathrm{constant}.
$$

### Using Poisson's equation
In cylindrical coordinates, Poisson's equation reads

$$
\frac{1}{R}\frac{\partial}{\partial R}\left(R \frac{\partial \Phi}{\partial R}\right) + \frac{\partial^2 \Phi}{\partial z^2} = 4\pi G \Sigma \delta(z).
$$

Since the potential only depends on $z$, the first term vanishes and we have

$$
\frac{d^2 \Phi}{dz^2} = 4\pi G\Sigma \delta(z).
$$

which can be integrated twice to obtain the same result as before.

<!-- ======================= -->
<!-- PROBLEM 2.4             -->
<!-- ======================= -->
## Problem 2.4

### a. Laplace's equation

For $\Phi = \ln[r(1 + |\cos\theta|)]$ we have

1. The potential is not defined at $r=0$.
2. At $\theta = \pi/2$ the potential is continuous but not differentiable.

At any other location all derivatives are smooth and well defined

$$
\begin{align}
\nabla^2 \Phi &= \frac{1}{r^2}\frac{\partial}{\partial r}\left(r^2 \frac{\partial \Phi}{\partial r}\right) + \frac{1}{r^2 \sin\theta}\frac{\partial}{\partial \theta}\left(\sin\theta \frac{\partial \Phi}{\partial \theta}\right) \\
&= \frac{1}{r^2} + \frac{1}{r^2\sin\theta}\frac{d}{d\theta}\left(\sin\theta \frac{d}{d\theta}(\ln(1 + |\cos\theta|))\right) \\
\end{align}
$$

If $0 < \theta < \pi/2$

$$
\nabla^2 \Phi = \frac{1}{r^2} + \frac{1}{r^2\sin\theta}\frac{d}{d\theta}\left(\sin\theta \frac{d}{d\theta}(\ln(1 + \cos\theta))\right) = \frac{1}{r^2} - \frac{1}{r^2} = 0,
$$

and if $\pi/2 < \theta < \pi$

$$
\nabla^2 \Phi = \frac{1}{r^2} + \frac{1}{r^2\sin\theta}\frac{d}{d\theta}\left(\sin\theta \frac{d}{d\theta}(\ln(1 - \cos\theta))\right) = \frac{1}{r^2} - \frac{1}{r^2} = 0.
$$

In both cases, the Laplacian vanishes. Thus, the density is zero everywhere except at $r=0$ and $\theta = \pi/2$. Because of the absolute value, we conclude that a potential of the form $\Phi(R, z) = v_c^2\ln(R + |z|) + \mathrm{const}$ introduces a $\delta$-function at $\theta = \pi/2$, that is $\rho(R, z) = \Sigma(R)\delta(z)$.

### b. Disk density
We can integrate Poisson's equation across a thin pillbox of thickness $2h$ around the plane $z=0$.

$$
\begin{align}
\int_{-h}^{+h} dz \, \nabla^2 \Phi &= \left(\frac{\partial \Phi}{\partial z}\right)_{-h}^{h} = 4\pi G \Sigma(R) \\
&= \frac{v_c^2}{R + h} - \left(-\frac{v_c^2}{R + h}\right) = \frac{2v_c^2}{R + h} \\
&\to \frac{2v_c^2}{R} \quad\text{for } R \to \infty
\end{align}
$$

That is

$$
\Sigma(R) = \frac{v_c^2}{2\pi G R}.
$$

<!-- ======================= -->
<!-- PROBLEM 2.5             -->
<!-- ======================= -->
## Problem 2.5

### a. Constant circular speed

$$
M(r) = 2\pi A \int_0^r dr' = 2\pi A r \quad\text{for } r < R_0,
$$

and $M(r > R_0) = M(R_0) = 2\pi A R_0$. The cicular speed is then

$$
v_c^2(r) = 2\pi AG \begin{cases}
1 & r < R_0 \\
R_0/r & r \geq R_0
\end{cases}
$$

### b. Flattened model

Consider now the density profile

$$
\rho(m^2) = \frac{A}{2m^2} \quad\text{with}\quad m^2 = R^2 + \frac{z^2}{q^2}.
$$

The circular velocity is

$$
\begin{align}
v_c^2(R) &= 4\pi G \sqrt{1 - e^2}\int_0^R dm \, \frac{m^2\rho(m^2)}{\sqrt{R^2 - m^2 e^2}} \\
&= 2\pi AG \sqrt{1 - e^2}\int_0^R dm \, \frac{1}{\sqrt{R^2 - m^2 e^2}} \\
&= 2\pi AG q \frac{\arcsin e}{e}
\end{align}
$$

Note that in the limit $q\to 1$ (i.e. $e\to 0$) we have that $\arcsin e / e \to 1$ and we recover the result from part (a).

### c. Razor-thin limit

At fixed $R < R_0$ the body model extends in $z$ only for $z < q\sqrt{R_0^2 - R^2}$, and the density in the limit $q\to 0$ becomes

$$
\begin{align}
\Sigma(R) &= \lim_{q\to 0} \int_{-\infty}^{+\infty} dz \, \rho(m^2) = 2\int_{0}^{q\sqrt{R_0^2 - R^2}} dz \, \frac{A}{2(R^2 + z^2/q^2)} \\
&= Aq \int_0^{\sqrt{R_0^2 - R^2}} \frac{du}{R^2 + u^2} = \frac{Aq}{R}\arctan\left(\frac{\sqrt{R_0^2 - R^2}}{R}\right) \\
&= \frac{Aq}{R}\arccos\left(\frac{R}{R_0}\right).
\end{align}
$$

<!-- ======================= -->
<!-- PROBLEM 2.6             -->
<!-- ======================= -->
## Problem 2.6

$$
\begin{align}
R^2 + (a + |z|)^2 &= a^2\sinh^2u \sin^2v + (a + a \cosh u |\cos v|)^2 \\
&= a^2(\sinh^2u \sin^2v + 1 + 2\cosh u |\cos v| + \cosh^2u \cos^2v) \\
&= a^2(\cosh^2u + \cos^2v + 2\cosh u |\cos v|) \\
&= a^2(\cosh u + |\cos v|)^2
\end{align}
$$

The Kuzmin potential is then

$$
\begin{align}
\Phi_K &= -\frac{GM}{\sqrt{R^2 + (a + |z|)^2}} = -\frac{GM}{a(\cosh u + |\cos v|)} \\
&= -\frac{GM}{a}\frac{1}{\cosh u + |\cos v|}\frac{\cosh u - |\cos v|}{\cosh u - |\cos v|} \\
&= -\frac{GM}{a} \frac{\cosh u - |\cos v|}{\sinh^2u + \sin^2v}
\end{align}
$$

<!-- ======================= -->
<!-- PROBLEM 2.7             -->
<!-- ======================= -->
## Problem 2.7

The an answer is **No**

Outside a bounded mass, the Newtonian potential is determined only by the multipole moments of the interior density. Saying

$$
\Phi(r) = -\frac{GM}{r} \quad\text{for } r > R
$$

is equivalent to saying that all exterior multipoles with $l \geq 1$ vanish. This is only possible if the density is spherically symmetric inside $r < R$. That does not force the interior density to be spherically symmetric, only that the integrated multipole moments are zero.

As an example consider a spherically symmetric mas $M_0$ and add two thin concentric shells at radiia $0 < b < a < R$, whose surface densities are given by

$$
\Sigma_a(\Omega) = \Sigma_{0,a} + \epsilon Y_{20}(\Omega), \quad \Sigma_b(\Omega) = \Sigma_{0,b} - \epsilon (a / b)^2 Y_{20}(\Omega),
$$

For a thin shell at radius $r$, the exterior $l = 2$ coefficient scales as $r^2$ times the shell's $Y_{20}$ coefficient. TThus the two sheel give equal and opposite contributions to the total $l=2$ multipole moment for $r > a$. The next exterior potential for $r > R$ is exactly $-GM/r$ but the interior density is not spherically symmetric.

<!-- ======================= -->
<!-- PROBLEM 2.8             -->
<!-- ======================= -->
## Problem 2.8

For a density mass-to-light ratio $\Upsilon$ we know that $\rho(r) = \Upsilon j(r)$, and from (Problem 1.3)[./chapter01.3.md] the mass density is

$$
\rho(r) = -\frac{\Upsilon}{2\pi r}\frac{dS}{dr},
$$

therefore the potential energy is

$$
\begin{align}
\Phi(r) &= - 4\pi G \left[\frac{1}{r}\int_0^r dr' \, r'^2 \rho(r') + \int_r^{+\infty} dr' \, r' \rho(r')\right] \\
&= -4\pi G  \left[-\frac{\Upsilon}{2\pi r}\int_0^r dr' \, r' \frac{dS}{dr'} - \frac{\Upsilon}{2\pi}\int_r^{+\infty} dr' \, \frac{dS}{dr'}\right] \\
&= 2G\Upsilon \left[ \frac{1}{r}\left(rS(r) - \int_0^r dr'\, S(r')\right) - S(r)\right] \\
&= -\frac{2G\Upsilon}{r}\int_0^r dr' \, S(r')
\end{align}
$$

<!-- ======================= -->
<!-- PROBLEM 2.9             -->
<!-- ======================= -->
## Problem 2.9

$$
\begin{align}
\int_0^\infty dR\, I(R) &= 2\int_0^\infty dR \int_R^\infty dr \, \frac{r j(r)}{\sqrt{r^2 - R^2}} \\
&= 2\int_0^\infty dr \, r j(r) \int_0^r \frac{dR}{\sqrt{r^2 - R^2}} \\
&= \frac{2}{\Upsilon}\int_0^\infty dr \, r \rho(r) \cdot \frac{\pi}{2} \\
&= -\frac{\Phi(0)}{4G\Upsilon}
\end{align}
$$

<!-- ======================= -->
<!-- PROBLEM 2.10             -->
<!-- ======================= -->
## Problem 2.10

### a. Quadrupole expansion

Consider an axisymmetric density distribution $\rho(R, z)$ with total mass $M$. Since the density is axisymmetric, $\rho_{lm} = \delta_{m0} \rho_{l0}$, and

$$
\rho_{l0}(a) = \int d\Omega Y_{l0}^*(\Omega) \rho(a, \Omega),
$$

furthermore, since $\rho(R, z) = R(\rho, -z)$, all odd $l$ terms vanish. Thus, the multipole expansion reduces to

$$
\begin{align}
\Phi(r, \Omega) &= -4\pi G\sum_{l, \, \text{even}} \frac{Y_{l0}(\Omega)}{2l + 1}\left[\frac{1}{r^{l+1}}\int_0^r da \, a^{l+2} \rho_{l0}(a) + r^l \int_r^\infty da \, a^{1-l} \rho_{l0}(a)\right] \\
&= -4\pi G\sum_{l, \, \text{even}} \frac{Y_{l0}(\Omega)}{2l + 1}\left[\frac{1}{r^{l+1}}\int_0^r da \, a^{l+2} \int d\Omega' Y_{l0}^*(\Omega') \rho(a, \Omega') + r^l \int_r^\infty da \, a^{1-l} \int d\Omega' Y_{l0}^*(\Omega') \rho(a, \Omega')\right] \\
&= -4\pi G\sum_{l, \, \text{even}} \frac{Y_{l0}(\Omega)}{2l + 1}\left[\frac{1}{r^{l+1}}\int_0^r d^3\mathbf{x}' \, a^{l} Y_{l0}^*(\Omega') \rho(\mathbf{x}') + r^l \int_r^\infty d^3\mathbf{x}' \, a^{-l-1} Y_{l0}^*(\Omega') \rho(\mathbf{x}')\right] \\
&\equiv \sum_{l, \, \text{even}} \Phi_l(r, \Omega)
\end{align}
$$

**Monopole term ($l=0$)**

$$
Y_{00}(\Omega) = 1/\sqrt{4\pi}
$$

thus

$$
\begin{align}
\Phi_0(r, \Omega) &= -G \left[\frac{1}{r}\int_0^r d^3\mathbf{x}' \, \rho(\mathbf{x}') + \int_r^\infty d^3\mathbf{x}' \, \frac{1}{a'} \rho(\mathbf{x}')\right] \\
&= -\frac{GM}{r} - G\int_r^\infty d^3\mathbf{x}' \, \frac{1}{a'} \rho(\mathbf{x}'),
\end{align}
$$

where we have used the fact that $r > r_{\rm max}$.

**Quadrupole term ($l=2$)**

$$
Y_{20}(\Omega) = \sqrt{\frac{5}{4\pi}} \frac{1}{2}(3\cos^2\theta - 1)
$$

thus

$$
\begin{align}
\Phi_2(r, \Omega) &= -\frac{(3 \cos^2\theta - 1)G}{4r^3}\int_0^r d^3\mathbf{x}' \, a'^{2} (3\cos^2\theta' - 1) \rho(\mathbf{x}') \\
& - \frac{(3 \cos^2\theta - 1)Gr^2}{4}\int_r^\infty d^3\mathbf{x}' \, a'^{-3} (3\cos^2\theta' - 1) \rho(\mathbf{x}') \\
&= -\frac{3z^2 - r^2}{4r^5} G \int d^3\mathbf{x}' \, (3z'^2 - a'^2) \rho(\mathbf{x}') \\
&- \frac{3z^2 - r^2}{4} Gr^2 \int d^3\mathbf{x}' \, \frac{(3z'^2 - a'^2)}{a'^5} \rho(\mathbf{x}')
\end{align}
$$

Adding the two terms we obtain

$$
\Phi(R, z) \approx -\frac{GM}{r} - \frac{G}{4}\frac{R^2 - 2z^2}{r^5}\int d^3\mathbf{x}' \, (R'^2 - 2z'^2) \rho(\mathbf{x}')
$$

### b. Exponential disk

For $\rho(R, z) = \Sigma_0\exp(-R/R_d)\delta(z)$ we have

$$
\begin{align}
\Phi(R, z) &\approx -\frac{GM}{r} - \frac{G}{4}\frac{R^2 - 2z^2}{r^5}\int_0^\infty dR' \, R' \int_0^{2\pi} d\phi' \, (R'^2 \cos^2\phi' - 2\cdot 0^2) \Sigma_0 e^{-R'/R_d} \\
&= -\frac{GM}{r} - \frac{G\pi \Sigma_0}{2}\frac{R^2 - 2z^2}{r^5}\int_0^\infty dR' \, R'^3 e^{-R'/R_d} \\
&= -\frac{GM}{r} - 3\pi G \Sigma_0 R_d^4 \frac{R^2 - 2z^2}{r^5} \\
&= -\frac{GM}{r} - \frac{3\pi GM}{2} \frac{R^2 - 2z^2}{r^5}
\end{align}
$$

## References
\bibliography
