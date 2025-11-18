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


## References
\bibliography
