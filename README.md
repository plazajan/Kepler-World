Kepler's World

See abstract.pdf and coursePlanning.txt for info on using these course 
materials in the classroom.

A simple simulation (in Python) of movements of celestial bodies 
around the Sun - a computaitonal experiment for the purpose of 
studying Kepler's laws and how they emerge from Newton's law of gravity.
Designed for high school students or lower level undergraduate students.

Movements of the eight planets are simulated, as well as movements of
some celestial bodes made up for the purpose of experiments.
The simulation uses a discrete approximation of the real-world continuous effect
of the Newton's law of gravity (F = GMm/(r^2) on the moving body,
taking into account only that body and the Sun (the "two-body problem").
The calculations are done for the actual celestial body,
and the results are scaled down to be visualized in turtle graphics.
The simulation tests if the resulting orbits and bodies' movements obey
Kepler's laws.

Concepts of the methodology of science:
* Branches of natural science: 
  **physics, chemistry, biology, Earth and atmospheric (planetary) science.**
* **Natural science** vs. **mathematics and logic**
* **Empirical evidence.**
* Methods of science: 
  **experimental/observational, theoretical/mathematical, computational.**
* **Laws of science** vs **theorems of mathematics**
  Newton's laws and Kepler's laws, theorems about ellipses.
* The **scientific method.**
* **Units of measurement** and **dimensional analysis:** 
  units of the gravitational constant G,
  Vis Viva Equation: v^2 = G(M+m)(2/r - 1/a)
  Orbital period T = 2 pi sqrt( a^3 / (G(M+m) )
* **Discrete** simulation of a **continuous** process.
* **Approximation** and **approximation error**.
* **Absolute error** and **relative error**.
* **Rounding** numbers.
* **Floating-point arithmetic**
* **Recurrent/periodic process** - for instance planet's movement around the Sun.
* **Qualitative** vs. **quantitative** statements
  ("a planet moves fastest when close to the Sun" vs. Kepler's 2nd law).
* **Local** effects (of the law of gravity) vs.
  **global** properties (of orbits described by Kepler's laws).
* **Visualization** vs. **simulation**.
* Note: Newton gave a **(mathematical) proof**, using calculus
  deriving Kepler's laws from base laws of motion and gravity.

The visualization shows orbits to scale, 
but disregards some details not relevant to Kepler's laws.
Namely, in this simulation (unlike in reality):
* the orbits of all planets are in the same plane, and
* the foci of elliptical orbits of all planets are on the same straight line.

The period in which the turtle completes the orbit is proportional to the
actual orbital period of the planet, but the scaling factor is not the same
as for the orbit's size.

A slideshow is also available: Discovering the Mechanics of the Solar System.
