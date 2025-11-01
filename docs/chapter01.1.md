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
    For Neptune, and $\rho = 0.1 \Msun/\rmpc^3$, we get $\delta\varpi = -1.5\times 10^{-8}\rmarcsec/\rmyr$ this differs from the value given in Binney & Tremaine (2008) [@binneytremaine2008] by six orders of magnitude; unfortunately I cannot come up with a better derivation that can explain the discrepancy.

<!-- ======================= -->
<!-- PROBLEM 1.2             -->
<!-- ======================= -->
## Problem 1.2

### a. Luminosity density to surface brightness
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
>>> from galactic_dynamics_bt.chapter01.sersic_profile import plot_sersic_profile
>>> plot_sersic_profile()
```


![FRW Model Phase Diagram](assets/generated/sersic_profile.png)

*Figure P1.2: Sersic luminosity density profile derived from the surface brightness profile using the inversion formula. For comparison different values of the parameter $m$ are shown.*

<!-- ======================= -->
<!-- PROBLEM 1.3            -->
<!-- ======================= -->
## Problem 1.3

### a. Strip brightness

$$
\begin{align}
S(x) &= \int_{-\infty}^{+\infty} dy I(\sqrt{x^2 + y^2}) \\
&= 2\int_{0}^{+\infty} dy I(\sqrt{x^2 + y^2}) \quad\text{with}\quad R^2 = x^2 + y^2 \\
&= 2\int_x^\infty dR \frac{R I(R)}{\sqrt{R^2 - x^2}}
\end{align}
$$

### b. From strip brightness to brightness profiles

Start with

$$
\begin{align}
S(x) &= 2\int_x^\infty dR \frac{R I(R)}{\sqrt{R^2 - x^2}} \\
&= 4\int_x^\infty dR\int_R^\infty dr \frac{r j(r)}{\sqrt{r^2 - R^2}\sqrt{R^2 - x^2}} \\
&= 4\int_x^\infty dr r j(r) \underbrace{\int_x^r dR \frac{R}{\sqrt{r^2 - R^2}\sqrt{R^2 - x^2}}}_{\pi/2} \\
&= 2\pi\int_x^\infty dr r j(r).
\end{align}
$$

Taking the derivative with respect to $x$ we have

$$
j(x) = -\frac{1}{2\pi x}\frac{dS}{dx}.
$$

The cumulative luminosity inside radius $r$ is given by


$$
L(r) = 4\pi \int_0^r dx x^2 j(x) = 4\pi \int_0^r dx x^2\left[-\frac{1}{2\pi x}\frac{dS}{dx}\right] = -2\int_0^r dx x \frac{dS}{dx}.
$$


<!-- ======================= -->
<!-- PROBLEM 1.4            -->
<!-- ======================= -->
## Problem 1.4

### a. Central surface brightness of an axisymmetric galaxy
We are assuming the axisymmetric density density distribution $j$ can be written as

$$
j = j(m) \quad\text{ with }\quad m^2 = R^2 + \frac{z^2}{q^2}.
$$

The projected surface brightness is simple the line-of-sight integral of the luminosity density

$$
I(R) = \int_{-\infty}^{+\infty} ds j(m).
$$

where $s$ measures distance along the line of sight. Let's consider two cases

**Face-on view** For a face-on line of sight, we have $R = 0$ and $m = z/q$, thus the central intensity is

$$
I_n = 2\int_0^\infty dz j(m) = 2q\int_0^\infty dm j(m).
$$

**Inclined view** The coodinate transformation for the central line of sight inclined by an angle $i$ is given by $R = s\sin i$ and $z = s\cos i$, so that the ellipsoidal radius along the line of sight is given by

$$
m^2 = s^2\sin^2 i + \frac{s^2\cos^2 i}{q^2} = s^2\left(\sin^2 i + \frac{\cos^2 i}{q^2}\right).
$$

Which means that

$$
\frac{dm}{ds} = \sqrt{\sin^2 i + \frac{\cos^2 i}{q^2}},
$$

and

$$
I_0(i) = 2\int_0^\infty ds j(m) = 2\int_0^\infty dm \frac{j(m)}{dm/ds} = \frac{2}{\sqrt{\sin^2 i + \cos^2 i/q^2}}\int_0^\infty dm j(m).
$$

The ratio is then

$$
\frac{I_0(i)}{I_n} = \frac{1/q}{\sqrt{\sin^2 i + \cos^2 i/q^2}}
$$

We now consider two cases depending on $q$

**Oblate case** ($q<1$) In this case $Q^2 = \cos^2 i + q^2\sin^2 i < 1$ and we can write

$$
I_0 = \frac{I_n}{Q}
$$

$Q$ decreases from $1$ (face-on) to $q$ (edge-on).

**Prolate case** ($q>1$) In this case $Q^2 = \sin^2 i + \cos^2 i/q^2 < 1$ and we can write

$$
I_0 = \frac{I_n}{\sqrt{q^2 + 1 - q^2Q^2}}
$$

$Q$ decreases from $1$ (viewed along the long axis) to $1/q$ (viewed perpendicular to it).

### b. Relation between apparent and intrinsic axis ratios

The relation is given by

$$
\begin{align}
Q^2 &= \cos^2 i + q^2\sin^2 i, \quad\text{oblate case}\\
Q^2 &= \sin^2 i + \frac{\cos^2 i}{q^2}, \quad\text{prolate case}
\end{align}
$$

### c. Probability distribution of apparent axis ratios

For random orientation, the inclination $u$ has a PDF $p(i) = \sin i$ for $i\in[0,\pi/2]$.
The fraction of galaxiess seen from a line of sight that lies within an angle $x$ of the symmetry axis is given by

$$
f_{\textrm{axis}} = \frac{\int_0^x di\sin i}{\int_0^{\pi/2} di\sin i} = 1 - \cos x.
$$

<!-- ======================= -->
<!-- PROBLEM 1.5            -->
<!-- ======================= -->
## Problem 1.5

### a. Mass-to-light ratio
The luminosity can be calculate from the observed flux $F$ as

$$
L = 4\pi d_L^2 F,
$$

where $d_L$ is the luminosity distance, which for small redshifts can be approximated as

$$
d_L \approx \frac{cz}{H_0}.
$$

That is $L \propto d_L^2 \propto H_0^{-2} \propto h_7^{-2}$, hence $\Upsilon_R \propto h_7^{-2}$.

### b. Correcting Zwicky's estimate

Assuming

$$
H_0 = 558 \rmkm/\rms^{-1}\rmMpc^{-1} = 70(7.97) \rmkm/\rms^{-1}\rmMpc^{-1}
$$

Since $\Gamma \propto h_7^2$, his estimated mass-to-light ratio should be corrected by a factor or $7.97^{-2}$, that is, his corrected estimadated value for the mass-to-light ratio is $~6$.


<!-- ======================= -->
<!-- PROBLEM 1.6             -->
<!-- ======================= -->
## Problem 1.6

In a flat universe completely dominated by $\Lambda$ we have

$$
a(t) \propto \exp\left(H_0t\right) = \exp\left[\left(\frac{8\pi G \rho_\Lambda}{3}\right)^{1/2}t\right].
$$

$(G\rho_\Lambda)^{-1/2}$ is then the e-folding time of expansion.

<!-- ======================= -->
<!-- PROBLEM 1.7             -->
<!-- ======================= -->
## Problem 1.7

Let's start from

$$
dr^2 = a^2(t)\left[\frac{dx^2}{1 - kx^2/x_u^2} + x^2(d\theta^2 + \sin^2\theta d\phi^2)\right].
$$

For $k=+1$ the spatial metric is a diagonal with components

$$
g_{xx} = \frac{a^2(t)}{1 - x^2/x_u^2}, \quad g_{\theta\theta} = a^2(t)x^2, \quad g_{\phi\phi} = a^2(t)x^2\sin^2\theta.
$$

So the determinant is given by

$$
\det g_3 = \frac{a^6(t)x^4\sin^2\theta}{(1 - x^2/x_u^2)},
$$

and the 3-volume is

$$
V/2 = \int dxd\theta d\phi\sqrt{\det g_3} = a^3(t) \int_0^{x_u} dx \frac{x^2}{\sqrt{1 - x^2/x_u^2}} \int_0^\pi d\theta \sin\theta \int_0^{2\pi} d\phi = \pi^2 a^3(t) x_u^3.
$$

Note that we integrated only over half the 3-sphere, hence the factor of $1/2$ in the left-hand side.
The 3-sphere of radius $x_u$ can be embedded in 4D Euclidean space, when we run the integral on the right-hand side, we are only integrating over half the 3-sphere, hence the factor of $1/2$ in the left-hand side.

The volume of the universe is then

$$
V = 2\pi^2 a^3(t) x_u^3.
$$

<!-- ======================= -->
<!-- PROBLEM 1.8             -->
<!-- ======================= -->
## Problem 1.8

### a. Static universe

Start from the acceleration equation (1.49)

$$
\frac{\ddot{a}}{a} = -\frac{4\pi G}{3}\left(\rho + \frac{3p}{c^2}\right).
$$

For a static universe $\ddot{a} = 0$, if $\rho_\gamma = 0$, then this equation implies that

$$
\begin{align}
0 &= -\frac{4\pi G}{3}\left((\rho_m + \rho_\Lambda) + \frac{3p}{c^2}(p_m + p_\Lambda)\right) \\
&= -\frac{4\pi G}{3}\left(\rho_m + \rho_\Lambda + \frac{3}{c^2}(0 -\rho_\Lambda c^2)\right) \\
&= -\frac{4\pi G}{3}\left(\rho_m + \rho_\Lambda - 3\rho_\Lambda\right) \\
\end{align}
$$

which means that

$$
\rho_m = 2\rho_\Lambda.
$$

### b. Unstable universe

Consider the perturbation $a(t) = a_0 + \delta a_1(t) = a_0[1 + \epsilon x(t)]$ around the static solution $a_0$, with $x(t) \equiv a_1(t) / a_0$ and $\epsilon \ll 1$. Each component of the density evolves as

| Component | Density evolution | Perturbation |
|-----------|-------------------|--------------|
| Matter    | $\rho_m(t) = \rho_{m0}\left(\frac{a(t)}{a_0}\right)^{-3} \approx \rho_{m0}(1 - 3\epsilon x(t))$ | $\delta \rho_m(t) \approx -3\epsilon \rho_{m0} x(t)$ |
| Dark energy | $\rho_\Lambda(t) = \rho_{\Lambda 0}$ | $\delta \rho_\Lambda(t) = 0$ |

The perturbed acceleration equation becomes

$$
\begin{align}
\delta(\ddot{a}/a) &= -\frac{4\pi G}{3}\left(\delta \rho_m + \delta \rho_\Lambda + \frac{3}{c^2}(\delta p_m + \delta p_\Lambda)\right) \\
&= -\frac{4\pi G}{3}(-3\epsilon \rho_{m0} x(t)) = 4\pi G \epsilon \rho_{m0} x(t)\\
\epsilon \ddot{x} &= 4\pi G \epsilon \rho_{m0} x(t).
\end{align}
$$

Whose solution is

$$
x(t) = A e^{\sqrt{4\pi G \rho_{m0}}t} + B e^{-\sqrt{4\pi G \rho_{m0}}t}.
$$

where $A$ and $B$ are constants determined by the initial conditions.
The perturbation grows exponentially with time, with characteristic timescale

$$
t_{\textrm{inst}} = \frac{1}{\sqrt{4\pi G \rho_{m0}}}.
$$


<!-- ======================= -->
<!-- PROBLEM 1.9            -->
<!-- ======================= -->
## Problem 1.9

For a flat universe,

$$
H^2 = \frac{\dot{a}^2}{a^2} = \frac{8\pi G}{3}\left( \rho_{m0}a^{-3} + \rho_{\Lambda 0}\right) =
H_0^2\left(\Omega_{m0}a^{-3} + \Omega_{\Lambda 0}\right).
$$

In the distant future, $a\gg 1$, the matter term becomes negligible and we have

$$
H_\infty^2 \approx H_0^2 \Omega_{\Lambda 0} < H_0^2.
$$

The timescale for the expansion is then $a(t)\approx e^{H_\infty t}$, which is accelerated expansion since $\ddot{a} = H_\infty^2 e^{H_\infty t} > 0$.

<!-- ======================= -->
<!-- PROBLEM 1.10            -->
<!-- ======================= -->
## Problem 1.10

The mean free path has the form

$$
\ell = \frac{1}{n\pi r^2} = \frac{4r \rho_{\rm int}}{3\Omega_a \rho_c},
$$

where $\rho_{\rm int}$ is the internal density of the asteroids, $\rho_c$ is the critical density of the universe, and $\Omega_a$ is the cosmological density parameter of asteroids. Requiring that the Universe is optically thin to distance quasars over a Hubble distance $D \sim c / H_0$ gives

$$
\tau \sim \frac{D}{\ell} \lesssim 1 \implies r \gtrsim \frac{3\Omega_a \rho_c}{4\rho_{\rm int}} \frac{c}{H_0}.
$$

Using $\rho_c = 1.88\times 10^{-29} h^2 \, \mathrm{g \, cm^{-3}}$, $c / H_0 \simeq 9.26\times 10^{27} h^{-1} \, \mathrm{cm}$, $\Omega_a = 0.1$ and $\rho_{\rm int} = 8 \, \mathrm{g \, cm^{-3}}$, we find

$$
r_{\min} \approx 1~\mu\mathrm{m}.
$$

<!-- ======================= -->
<!-- REFERENCES.             -->
<!-- ======================= -->

## References
\bibliography
