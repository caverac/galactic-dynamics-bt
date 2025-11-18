<!-- ======================= -->
<!-- PROBLEM 2.21            -->
<!-- ======================= -->
## Problem 2.21

Define

$$
\begin{align}
I(R, R', z) &\equiv \int_0^{2\pi}d\phi' \, \frac{1}{(R^2 + R'^2 - 2RR'\cos\phi' + z^2)^{1/2}} \\
&= \int_0^{2\pi}d\phi' \, \frac{1}{((R - R')^2 + z^2 + 4RR'\sin^2(\phi'/2))^{1/2}} \\
&= \int_0^{2\pi}d\phi' \, \frac{1}{((R + R')^2 + z^2 - 4RR'\cos^2(\phi'/2))^{1/2}} \\
&= \int_0^{2\pi}d\phi' \, \frac{1}{(a^2 - b^2\cos^2(\phi'/2))^{1/2}} \quad a^2 \equiv (R + R')^2 + z^2, \, b^2 \equiv 4RR' \\
&= 2\int_0^{\pi}d\theta \, \frac{1}{(a^2 - b^2\cos^2\theta)^{1/2}} \quad\text{with } \theta = \phi'/2 \\
&= 2\int_0^{\pi}d\theta \, \frac{1}{(a^2(1 - k^2) + a^2k^2\sin^2\theta)^{1/2}} \quad k^2 \equiv \frac{b^2}{a^2} = \frac{4RR'}{(R + R')^2 + z^2} \\
&= \frac{2}{a}\int_0^{\pi}d\theta \, \frac{1}{(1 - k^2\sin^2\theta)^{1/2}} \\
&= \frac{4}{\sqrt{(R + R')^2 + z^2}}K(k),
\end{align}
$$

where

$$
K(k) = \int_0^{\pi/2} \frac{d\theta}{(1 - k^2\sin^2\theta)^{1/2}}.
$$

For a thin axisymmetric disk with density $\rho(R, z) = \Sigma(R)\delta(z)$, the potential is thus

$$
\begin{align}
\Phi(R, z) &= -G\int_0^{\infty}R'dR'\int_0^{2\pi}d\phi' \, \frac{\Sigma(R')\delta(z')}{(R^2 + R'^2 - 2RR'\cos\phi' + (z - z')^2)^{1/2}} \\
&= -G\int_0^{\infty}R'dR' \,\Sigma(R') \, I(R, R', z) \\
&= -4G\int_0^{\infty}dR' \, \frac{R'\Sigma(R')}{\sqrt{(R + R')^2 + z^2}}K(k).
\end{align}
$$

Now we just need to massage this expression a bit. From the deinifition of $k$, we have

$$
\sqrt{(R + R')^2 + z^2} = \frac{2\sqrt{RR'}}{k},
$$

therefore

$$
\frac{R'}{((R + R')^2 + z^2)^{1/2}} = \frac{k}{2}\sqrt{\frac{R'}{R}},
$$

and the integral above becomes

$$
\Phi(R, z) = -\frac{2G}{\sqrt{R}}\int_0^{\infty}dR' \, \Sigma(R') k K(k)\sqrt{R'}.
$$

<!-- ======================= -->
<!-- PROBLEM 2.22            -->
<!-- ======================= -->
## Problem 2.22 🌶️

For a thin disk with surface density profile $\Sigma(\mathbf{R}')$ the potential on the plane is

$$
\Phi(\mathbf{R}, 0) = -G \int d^2\mathbf{R}' \, \frac{\Sigma(\mathbf{R}')}{|\mathbf{R} - \mathbf{R}'|}.
$$

### Laplacian on the plane

Define $\nabla^2_\perp = \partial^2_x + \partial^2_y$, then

$$
\nabla^2_\perp \Phi(\mathbf{R}, 0) = -G \int d^2\mathbf{R}'' \, \Sigma(\mathbf{R}'') \nabla^2_\perp \frac{1}{|\mathbf{R} - \mathbf{R}''|}.
$$

Define $I(\mathbf{R}')$ as

$$
\begin{align}
I(\mathbf{R}') &\equiv \int d^2\mathbf{R} \, \frac{1}{|\mathbf{R} - \mathbf{R}'|} \nabla^2_\perp \Phi(\mathbf{R}, 0) \\
&= \int d^2\mathbf{R} \, \frac{1}{|\mathbf{R} - \mathbf{R}'|} \left[ -G \int d^2\mathbf{R}'' \, \Sigma(\mathbf{R}'') \nabla^2_\perp \frac{1}{|\mathbf{R} - \mathbf{R}''|} \right] \\
&= -G \int d^2\mathbf{R}'' \, \Sigma(\mathbf{R}'') \int d^2\mathbf{R} \, \frac{1}{|\mathbf{R} - \mathbf{R}'|} \nabla^2_\perp \frac{1}{|\mathbf{R} - \mathbf{R}''|} \\
&= -G \int d^2\mathbf{R}'' \, \Sigma(\mathbf{R}'') \, K(\mathbf{R}', \mathbf{R}'')
\end{align}
$$

Where we have defined

$$
\begin{align}
K(\mathbf{R}', \mathbf{R}'') &= \int d^2\mathbf{R} \, \frac{1}{|\mathbf{R} - \mathbf{R}'|} \nabla^2_\perp \frac{1}{|\mathbf{R} - \mathbf{R}''|} \\
&= \int d^2\mathbf{u} \, \frac{1}{|\mathbf{u}|} \nabla^2_\perp \frac{1}{|\mathbf{u} + \mathbf{R}' - \mathbf{R}''|} \quad\text{with } \mathbf{u} = \mathbf{R} - \mathbf{R}'' \\
&= K(\mathbf{R}' - \mathbf{R}'') = K(\mathbf{\Delta}) \quad\text{with } \mathbf{\Delta} = \mathbf{R}' - \mathbf{R}''.
\end{align}
$$

### Interpreting the kernel

If we define $\varphi(\mathbf{R}) = 1 / |\mathbf{R}|$, then from our definition of $K$ we can see that it is nothing but the convolution of $\varphi$ with $\nabla^2_\perp \varphi$:

$$
K(\mathbf{\Delta}) = (\varphi * \nabla^2_\perp \varphi)(\mathbf{\Delta}).
$$

So its Fourier transform is simply the product of the Fourier transforms of each function:

$$
\begin{align}
\tilde{K}(\mathbf{k}) &= \tilde{\varphi}(\mathbf{k}) \cdot \mathcal{F}\{\nabla^2_\perp \varphi\}(\mathbf{k}) = \tilde{\varphi}(\mathbf{k}) \cdot (-k^2) \tilde{\varphi}(\mathbf{k}) = -k^2 |\tilde{\varphi}(\mathbf{k})|^2 \\
&= -k^2 \cdot\left( \frac{2\pi}{k} \right)^2 = -4\pi^2,
\end{align}
$$

Therefore, the inverse Fourier transform is

$$
K(\mathbf{\Delta}) = \mathcal{F}^{-1}\{\tilde{K}(\mathbf{k})\}(\mathbf{\Delta}) = -4\pi^2 \delta^{(2)}(\mathbf{\Delta}) = -4\pi^2 \delta^{(2)}(\mathbf{R}' - \mathbf{R}'').
$$

### Putting it all together

$$
\begin{align}
I(\mathbf{R}') &= -G \int d^2\mathbf{R}'' \, \Sigma(\mathbf{R}'') \, K(\mathbf{R}' - \mathbf{R}'') \\
&= 4\pi^2 G \int d^2\mathbf{R}'' \, \Sigma(\mathbf{R}'') \, \delta^{(2)}(\mathbf{R}' - \mathbf{R}'') \\
&= 4\pi^2 G \Sigma(\mathbf{R}') \\
&= \int d^2\mathbf{R} \, \frac{1}{|\mathbf{R} - \mathbf{R}'|} \nabla^2_\perp \Phi(\mathbf{R}, 0).
\end{align}
$$

That is

$$
\Sigma(x', y') = \frac{1}{4\pi^2 G} \iint dx\, dy \, \frac{1}{|\mathbf{R} - \mathbf{R}'|} \left( \frac{\partial^2\Phi}{\partial x^2} + \frac{\partial^2\Phi}{\partial y^2} \right).
$$

<!-- ======================= -->
<!-- PROBLEM 2.23            -->
<!-- ======================= -->
## Problem 2.23

With $\mathbf{r} = (l, m, n)$ and $\mathbf{k} = (k_x, k_y, k_z)$

$$
\begin{align}
(\nabla^2 \Phi)_{\mathbf{r}} &= \sum_{\mathbf{k}} \frac{\hat{\Phi}_{\mathbf{k}}}{\Delta^2}\left(e^{2\pi i (k_x(l+1) + k_ym + k_zn)/K} + e^{2\pi i (k_x(l-1) + k_ym + k_zn)/K} \right.\\
&\quad + \left.\cdots - 6e^{2\pi i (k_xl + k_ym + k_zn)/K}\right) \\
&= \sum_{\mathbf{k}} \frac{\hat{\Phi}_{\mathbf{k}}}{\Delta^2}\left( e^{2\pi i k_x/K} e^{2\pi i \mathbf{k}\cdot\mathbf{r}/K} + e^{-2\pi i k_x/K} e^{2\pi i \mathbf{k}\cdot\mathbf{r}/K} \right. \\
&\quad + \left. \cdots - 6e^{2\pi i \mathbf{k}\cdot\mathbf{r}/K} \right) \\
&= \sum_{\mathbf{k}} \frac{\hat{\Phi}_{\mathbf{k}}}{\Delta^2} \left(2\cos(2\pi k_x/K) + 2\cos(2\pi k_y/K) + 2\cos(2\pi k_z/K) - 6 \right) e^{2\pi i \mathbf{k}\cdot\mathbf{r}/K} \\
&= \sum_{\mathbf{k}} \frac{2\hat{\Phi}_{\mathbf{k}}}{\Delta^2} (\cos(2\pi k_x/K) + \cos(2\pi k_y/K) + \cos(2\pi k_z/K) - 3) e^{2\pi i \mathbf{k}\cdot\mathbf{r}/K} \\
&= 4\pi G \rho_{\mathbf{r}} \\
&= 4\pi G \sum_{\mathbf{k}} \hat{\rho}_{\mathbf{k}} e^{2\pi i \mathbf{k}\cdot\mathbf{r}/K}.
\end{align}
$$

From which,

$$
\hat{\Phi}_{\mathbf{k}} = \frac{2\pi G \Delta^2}{\cos(2\pi k_x/K) + \cos(2\pi k_y/K) + \cos(2\pi k_z/K) - 3} \hat{\rho}_{\mathbf{k}}.
$$

<!-- ======================= -->
<!-- PROBLEM 2.24            -->
<!-- ======================= -->
## Problem 2.24


### a. Force on a shell
The density for the system can be written as

$$
\rho(r) = \sum_{n=1}^N \frac{m_n}{4\pi r_n^2}\delta(r - r_n).
$$

The force at radius $r_{n - 1} < r < r_n$ is completely determined by the mass interior to $r$:

$$
\mathbf{F}(r) = -\frac{GM(r)}{r^2}\hat{e}_r = -\frac{G}{r^2}\left( \sum_{j=1}^{n - 1} m_j \right)\hat{e}_r = -\frac{GM_n}{r^2}\hat{e}_r.
$$

The for on the shell at $r_n$ is then

$$
\mathbf{F}_n = \frac{GM_n}{r_n^2}\hat{e}_r.
$$

### b. Computational complexity

The enclosed mass $M_n$ can be computed in $\mathcal{O}(N)$ operations by precomputing the cumulative sum of the masses: $M_n = M_{n-1} + m_n$.

After each $M_n$ is known, computing the force on each shell is $\mathcal{O}(1)$, so the total complexity for computing the forces on all shells is $\mathcal{O}(N)$.

The total computational complexity is therefore $\mathcal{O}(N)$.

## References
\bibliography
