# Multipole Expansion

For axially symmetric mass distributions, the density can be written as

$$
\rho(r, \Omega) = \rho(R, z) = \rho(r\sin\theta, r\cos\theta),
$$

that is, it does not depend on the azimuthal angle $\phi$. For this case Eq. (2.94) of [@BinneyTremaine2008] reduces to

$$
\rho_{lm}(a) = 2\pi \delta_{m0} \sqrt{\frac{2l + 1}{4\pi}} \int_{-1}^{1} du P_l(u) \rho(a\sqrt{1 - u^2}, au),
$$

which can numerically pre-computed for a given mass model. In the example below I use the Satoh model

$$
\rho(R, z) = \frac{ab^2M}{4\pi S^3(z^2 + b^2)}\left[ \frac{1}{\sqrt{z^2 + b^2}} + \frac{3}{a}\left(1 - \frac{R^2 + z^2}{S^2}\right)\right]
$$

The potential can be computed using the multipole expansion as

$$
\Phi(r \sin\theta, r\cos\theta) = -G\sum_{l=0}^{l_\mathrm{max}} P_l(\cos\theta)\sqrt{\frac{4\pi}{2l + 1}} \left[ \frac{1}{r^{l+1}} \int_0^r da\, a^{l+2} \rho_{l0}(a) + r^l \int_r^\infty da\, a^{1 - l} \rho_{l0}(a) \right],
$$

In code this looks like this

```python
>>> from galactic_dynamics_bt.chapter02.multipole_expansion import SatohModel, SatohModelParams
>>> model = SatohModel(SatohModelParams(q=0.6))
>>> expansion = MultipoleExpansion(model, 6)
```

This will create a multipole expansion up to $l_\mathrm{max} = 6$ for the Satoh model with flattening parameter $q = b/a = 0.6$. You can then evaluate the potential at any point $(R, z)$ using

```python
>>> R, z = 1.0, 0.5
>>> potential = expansion.potential(R, z)
```

![Multipole Expansion of Satoh Model](assets/generated/multipole_expansion_satoh.png)

*Figure: Multipole Expansion of Satoh Model, $b/a=0.6$. I included a flatter potential than the one in [@BinneyTremaine2008] to better illustrate the convergence of the expansion with $l_\mathrm{max}$. The solid black line shows the exact potential, other lines show the multipole expansion with different $l_\mathrm{max}$ values. $l_\mathrm{max} = 10$ is practically indistinguishable from the exact potential.*

<!-- ======================= -->
<!-- REFERENCES.             -->
<!-- ======================= -->

## References
\bibliography
