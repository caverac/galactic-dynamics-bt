<!-- ======================= -->
<!-- PROBLEM 1.11            -->
<!-- ======================= -->
## Problem 1.11

### a. Bound vs. unbound universes
From Eq. (1.50) we know that

$$
\dot{a}^2 - \frac{8\pi G \rho}{3}a^2 = 2E.
$$

At present day $a_0 = 1$ and this becomes (Neglecting any contribution from radiation)

$$
2E = H_0^2(1 - \Omega_0) = (1 - \Omega_{\Lambda 0} - \Omega_{m0})H_0^2 = -\frac{kc^2}{x_u^2}
$$

The boundary between a bound and unbound universe is then given by $k=0$ or $\Omega_{\Lambda 0} + \Omega_{m0} = 1$.


### b. Accelerating vs. decelerating universes
Now, Eq. (1.49) can be rewritten as

$$
\frac{\ddot{a}}{a} = \frac{H_0^2}{2}(2\Omega_{\Lambda 0} - \Omega_{m0}a^{-3}).
$$


The transition between a decelerating and accelerating universe occurs when $\ddot{a} = 0$, for all values of $a$, in particular at present day $a_0 = 1$, this happens when

$$
2\Omega_{\Lambda 0} - \Omega_{m0} = 0
$$

### c. Recollpapse vs. expand forever

Let's write the Friedmann equation as

$$
\frac{H^2}{H_0^2} = \Omega_{m0}a^{-3} + \Omega_{\Lambda 0} + (1 - \Omega_{m0} - \Omega_{\Lambda 0})a^{-2}.
$$

A condition for recollapse is that there exists a time $t_{\textrm{coll}}$ when $H(t_{\textrm{coll}}) = 0$, or equivalently a scale factor $a_{\textrm{coll}} > 1$ such that

$$
\Omega_{m0}a^{-3} + \Omega_{\Lambda 0} + (1 - \Omega_{m0} - \Omega_{\Lambda 0})a^{-2} = 0.
$$

Define the function

$$
f(a) = \Omega_{\Lambda 0}a^{3} + (1 - \Omega_{m0} - \Omega_{\Lambda 0})a + \Omega_{m0}.
$$

The separatrix in the parameter space between recollapsing and ever-expanding universes is given by the condition that $f(a)$ has a double root at some $a = a_{\textrm{coll}} > 1$. This requires that both $f(a_{\textrm{coll}}) = 0$ and $f'(a_{\textrm{coll}}) = 0$.

$$
0 = 3\Omega_{\Lambda 0}a^2 + (1 - \Omega_{m0} - \Omega_{\Lambda 0})
$$

A parametric solution of the separatrix is given by

$$
\Omega_{m0} = \frac{2 a^3}{1 - 3 a^2 + 2 a^3}, \quad
\Omega_{\Lambda 0} = \frac{1}{1 - 3 a^2 + 2 a^3}, \quad a > 1.
$$


A plot of the different regions in the $\Omega_{m0}$-$\Omega_{\Lambda 0}$ plane is shown below, generated with the following code:

```python
>>> from galactic_dynamics_bt.chapter01.frw_model import plot_frw_model
>>> plot_frw_model()
```

![FRW Model Phase Diagram](assets/generated/frw_models.png)

*Figure P1.11: Phase diagram of FRW cosmological models in the $\Omega_{m0}$-$\Omega_{\Lambda 0}$ plane. The solid line indicates the boundary between bound and unbound models, while the dashed line indicates the boundary between decelerating and accelerating models.*

<!-- ======================= -->
<!-- PROBLEM 1.12            -->
<!-- ======================= -->
## Problem 1.12

The age of the universe at redshift $z$ can be computed as

```python
>>> from galactic_dynamics_bt.chapter01.universe_age import find_universe_age
>>> z = 0
>>> h7 = 1.05
>>> find_universe_age(
...    z,
...    omega_m0=0.237,
...    omega_lambda0=0.763,
...    omega_gamma0=8.84e-5 / h7**2,
...    H0=70.0 * h7,
...)
14.26742...
```

For different redshifts we have

| Redshift | Age of the universe |
|----------|---------------------|
| 0        | 14.267 Gyr          |
| 1        | 6.3286 Gyr          |
| 1000     | 0.4605 Myr          |

The next figure shows the age of the universe as a function of redshift.

```python
>>> from galactic_dynamics_bt.chapter01.universe_age import plot_universe_age
>>> plot_universe_age()
```

![Age of the Universe vs Redshift](assets/generated/universe_age.png)
*Figure P1.12: Age of the universe as a function of redshift for flat cosmologies. Solid line represents the cosmology of Eq (1.73), dashed line are the results of the Planck 2018 cosmology [@planck2018].*


<!-- ======================= -->
<!-- PROBLEM 1.13            -->
<!-- ======================= -->
## Problem 1.13


### a. Scale factor at matter-radiation equality
Matter density scales as $\rho_m \propto a^{-3}$, while radiation density scales as $\rho_\gamma \propto a^{-4}$. The equality between matter and radiation densities occurs when

$$
\rho_m = \rho_\gamma \Rightarrow \Omega_{m0} a^{-3} = \Omega_{\gamma 0} a^{-4} \Rightarrow \Omega_{m0} a = \Omega_{\gamma 0}.
$$

This implies that the equality occurs at a scale factor

$$
1 + z_{\gamma m} = \frac{1}{a_{\gamma m}} = \frac{\Omega_{m0}}{\Omega_{\gamma 0}}
= 1.18\times 10^4 h_7^2 \Omega_{m0}.
$$


### b. Age of the universe at matter-radiation equality

Using the subsitution $u = \Omega_{\gamma 0} + a \Omega_{m 0}$ we can show that (For $\Omega_\Lambda = 0$)

$$
H_0t = \int_0^{\Omega_{\gamma 0}/\Omega_{m0}} da \frac{a}{\sqrt{\Omega_{m0} a + \Omega_{\gamma 0}}} =
2(2 - 2^{1/2})\frac{\Omega_{\gamma 0}^{1/2}(1 - \Omega_{m0})}{3\Omega_{m0}^2}.
$$

Since $\Omega_{m0} + \Omega_{\gamma 0} = 1$, we have that the age of the universe at matter-radiation equality is

$$
t_{\gamma m} = \frac{2(2 - 2^{1/2})}{3H_0}\frac{\Omega_{\gamma 0}^{3/2}}{\Omega_{m0}^2}.
$$

### c. Comoving horizon at matter-radiation equality

The comoving horizon at matter-radiation equality is given by

$$
x_{\gamma m} = c\int_0^{t_{\gamma m}} \frac{dt}{a(t)} = c\int_0^{a_{\gamma m}} \frac{da}{a^2 H(a)}
= c\int_0^{a_{\gamma m}} \frac{da}{H_0\sqrt{\Omega_{m0} a + \Omega_{\gamma 0}}}.
$$

Similarly as before, using the substitution $u = \Omega_{\gamma 0} + a \Omega_{m 0}$ we can show that

$$
x_{\gamma m} = (2^{1/2} - 1)\frac{c}{H_0} \frac{\Omega_{\gamma 0}^{1/2}}{\Omega_{m0}}.
$$

<!-- ======================= -->
<!-- PROBLEM 1.14            -->
<!-- ======================= -->
## Problem 1.14
The universe is not opaque for $z \lesssim 6$ even though it's reionized because the mean electron density is now extremely low due to cosmic expansion, so the Thomson scattering optical depth is tiny, and photons can propagate freely.

<!-- ======================= -->
<!-- REFERENCES.             -->
<!-- ======================= -->

## References
\bibliography
