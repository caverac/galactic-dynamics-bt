# Bifurcations of Orbit Families in a Flattened Logarithmic Potential


This is a small section to illustrate the claim about the bifurcation of orbit families in planar non-axisymmetric potentials made in Section 3.3. for this we will consider the potential

$$
\Phi(x, y) = \frac{1}{2}v_0^2\ln \left(x^2 + \frac{y^2}{q^2}\right), \quad 0 < q \leq 1.
$$

As we did in the [Surfaces of Section](chapter03-surfaces-of-section.md) section we numerically integrate orbits in the potential, using a symplectic integrator and calculate the corresponding Poincaré surfaces of section using the same procedure.

```python
>>> from galactic_dynamics_bt.chapter03.bifurcation import plot_loop_orbit_fraction
>>> plot_loop_orbit_fraction()
```


![Orbits in a flatten logarithmic potential](assets/generated/surfaces_of_section.png)

*Figure: 1: Surface of section for the flattened logarithmic potential with different flattening parameters.*

*Figure 1* shows the surfaces of section for different values of the flattening parameter $q$. As we decrease the value of $q$ from 1 (bottom-right panel) to 0.7 (top-right panel), we can see how the family of box orbits (the central region of the surface of section) starts to shrink, while the family of loop orbits (the outer region of the surface of section) starts to dominate. A quantitative measure of this effect is shown in *Figure 2*, where I plot the fraction of loop orbits as a function of the flattening parameter $q$.

```pythonpython
>>> from galactic_dynamics_bt.chapter03.bifurcation import  plot_loop_orbit_fraction()
>>> plot_loop_orbit_fraction()
```

![Fraction of loop orbits as a function of flattening](assets/generated/loop_orbit_fraction.png)

*Figure: 2: Fraction of loop orbits as a function of the flattening parameter $q$. As the potential becomes more flattened (smaller $q$), the fraction of loop orbits increases, indicating a bifurcation in the orbit families.*

This is certainly not a complete picture, since I only explored a very limitd grid of initial conditions, but it illustrates the main point: as the potential becomes rounder (larger $q$), the family of loop orbits becomes more dominant, while the family of box orbits shrinks, indicating a bifurcation in the orbit families.
