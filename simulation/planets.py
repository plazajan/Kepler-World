
#==================================================================================

from data import *

#=================================================================================
# PLANETS

class SimulatedPlanet(object):

    def __init__(self, name, mass, perihelionDistance, maxSpeed,
                 color="green"
                 # , timeStep=1000
                 ):
        """mass in kg, perihelionDistance in m, maxSpeed in m/s,
           simulation timeStep in s.
           Assumes that the Sun is at (0, 0) in a coordinate system.
           Puts the planet/self at coordinates (perihelionDistance, 0)
           giving it velocity (0, maxSpeed).
        """
        self._name = name
        self._mass = mass
        self._perihelion = perihelionDistance
        self._maxSpeed = maxSpeed
        self._color = color
        #self._timeStep = timeStep
        self.reset()

    def reset(self):
        """Return the planet to the perihelion giving it its maximum speed."""
        self._x = self._perihelion
        self._y = 0
        self._r2 = self._x*self._x + self._y*self._y # radius squared
        # we store radius squared to avoid loosing precision while recalculating it
        self._r = sqrt(self._r2) # radius = distance from Sun's center.
        self._vx = 0 # horizontal speed in m/s 
        self._vy = self._maxSpeed # vertical speed in m/s
        self._ax = -MU*self._x/(self._r2*self._r) # horizontal acceleration m/s^2
        self._ay = -MU*self._y/(self._r2*self._r) # vertical acceleration m/s^2
        mu = G * (MASS_SUN + self._mass) # gravitational parameter incl. planet      
        a = self._perihelion * mu / (
            2 * mu - self._perihelion * self._maxSpeed * self._maxSpeed) 
        t = 2 * pi * sqrt(a * a * a / mu) # predicted orbital period in s
        self._timeStep = round(sqrt(t/100_000)) # seconds

    def move(self):
        """Updates position, velocity, acceleration, radius
           to those after timeStep.
           timeStep is an optional parameter.
           If timeStep is is not provided, uses the default value
           with which the SimulatedPlanet object was created.
        """
        self._x += self._vx * self._timeStep
        self._y += self._vy * self._timeStep
        self._r2 = self._x * self._x + self._y * self._y
        self._r = sqrt(self._r2)
        self._vx += self._ax * self._timeStep
        self._vy += self._ay * self._timeStep
        self._ax = -MU * self._x / (self._r2 * self._r) # horiz acceleration
        self._ay = -MU * self._y / (self._r2 * self._r) # vertical acceleration

    def name(self): return(self._name)

    def mass(self): return self._mass

    def perihelion(self): return self._perihelion

    def maxSpeed(self): return self._maxSpeed

    def color(self): return(self._color)

    def position(self): return (self._x, self._y)

    def radius(self): return self._r

    def radiusSquared(self): return self._r2

    def velocity(self): return (self._vx, self._vy)

    def acceleration(self): return (self._ax, self._ay)

    def timeStep(self): return self._timeStep

#---------------------------------------------------------------------------------
# Simulated planets - global constants.

# The Sun is at (0,0).
# The starting position for a planet is (perihilionDistance, 0).
# The starting velocity is (0, sMax).
# The default time step in the simulation: 1000s.

# Inner planets - rocky:
MERCURY = SimulatedPlanet("Mercury",
                 MERCURY_DATA["mass"],
                 MERCURY_DATA["perihelion"],
                 MERCURY_DATA["max speed"],
                 "grey"
                )
VENUS   = SimulatedPlanet("Venus",
                 VENUS_DATA["mass"],
                 VENUS_DATA["perihelion"],
                 VENUS_DATA["max speed"],
                 "white"
                )
EARTH   = SimulatedPlanet("Earth",
                 EARTH_DATA["mass"],
                 EARTH_DATA["perihelion"],
                 EARTH_DATA["max speed"],
                 "turquoise"
                )
MARS    = SimulatedPlanet("Mars",
                 MARS_DATA["mass"],
                 MARS_DATA["perihelion"],
                 MARS_DATA["max speed"],
                 "red"
                )
# Outer planets - gas giants:
JUPITER = SimulatedPlanet("Jupiter",
                 JUPITER_DATA["mass"],
                 JUPITER_DATA["perihelion"],
                 JUPITER_DATA["max speed"],
                 "orange"
                )
SATURN  = SimulatedPlanet("Saturn",
                 SATURN_DATA["mass"],
                 SATURN_DATA["perihelion"],
                 SATURN_DATA["max speed"],
                 "gold"
                )
URANUS  = SimulatedPlanet("Uranus",
                 URANUS_DATA["mass"],
                 URANUS_DATA["perihelion"],
                 URANUS_DATA["max speed"],
                 "lightBlue1"
                )
NEPTUNE = SimulatedPlanet("Neptune",
                 NEPTUNE_DATA["mass"],
                 NEPTUNE_DATA["perihelion"],
                 NEPTUNE_DATA["max speed"],
                 "lightBlue2"
                )

ME = EARTH_DATA["mass"]
PE = EARTH_DATA["perihelion"]
S1_0 = sqrt(G * (MASS_SUN + ME) / PE) # the speed of a made-up planet
# with the same mass and perihelion as Earth, but with a circular orbit.

# Made-up planets, for computational experiments (default pencolor=green)
PLANET0_7 = SimulatedPlanet("Made-up planet 0.7", ME, PE, 0.7*S1_0) 
PLANET0_8 = SimulatedPlanet("Made-up planet 0.8", ME, PE, 0.8*S1_0)
PLANET0_9 = SimulatedPlanet("Made-up planet 0.9", ME, PE, 0.9*S1_0)
PLANET1_0 = SimulatedPlanet("Made-up planet 1.0", ME, PE,     S1_0)
PLANET1_1 = SimulatedPlanet("Made-up planet 1.1", ME, PE, 1.1*S1_0)
PLANET1_2 = SimulatedPlanet("Made-up planet 1.2", ME, PE, 1.2*S1_0)
PLANET1_3 = SimulatedPlanet("Made-up planet 1.3", ME, PE, 1.3*S1_0)

del ME
del PE
del S1_0

# Note.
# The inner planets and the made-up planets above are close to the Sun,
# and outer ones - are at vast distances.

#=================================================================================

