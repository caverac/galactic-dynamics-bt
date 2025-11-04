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


<!-- ======================= -->
<!-- REFERENCES.             -->
<!-- ======================= -->

## References
\bibliography
