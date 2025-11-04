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

## References
\bibliography
