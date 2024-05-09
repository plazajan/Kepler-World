"""
Kepler's World

A simple simulation of movements of celestial bodies around the Sun
for the purpose of studying Kepler's laws 
and how they emerge from Newton's law of gravity.
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

The visualization shows orbits to scale,
but disregards some details not relevant to Kepler's laws.
Namely, in this simulation (unlike in reality):
* the orbits of all planets are in the same plane, and
* the foci of elliptical orbits of all planets are on the same straight line.
The period in which the turtle completes the orbit is proportional to the
actual orbital period of the planet, but the scaling factor is not the same
as for the orbit's size.

A slideshow is also available: Discovering the Mechanics of the Solar System.

These matrials can be used as a basis for a dicussion of concepts of
methodology of science, as outlined in README.md.
"""

# TO DO:
# experiment1() - planets with masses m, 2m, 3m and the same perihelion.
# experiment2() - plenets with the same perihelion and max speed v, 2v, 3v
# ...
# planets with the same orbtital period and different perihelions

#==================================================================================

from turtle import getscreen, Turtle
from math import sqrt, pi, sin, cos

#==================================================================================
# ASTRONOMICAL DATA FROM OBSERVATIONS AND CALCULATIONS 

# Astronomical data in kilograms, meters, seconds.
# From Wikipedia
# or JPL-NASA: https://ssd.jpl.nasa.gov
# Orbit viewer: https://ssd.jpl.nasa.gov/tools/orbit_viewer.html

#---------------------------------------------------------------------------------
# Physical constants

G = 6.67430e-11 # m^3 / (kg * s^2). Gravitational constant
KM = 1000 # m
AU = 149_597_870_700 #m, Astronomical Unit = about the avg Earth-Sun distance
SEC = 1 # s
DAY = 24*60*60 #s (actual day is 86400.002 s because of a slowing rotation)
YR = 365.25*DAY # Julian year, used in astronomy

#---------------------------------------------------------------------------------
# Sun

MS = 1.98847e30 # Mass in kg.
RS = 696_000*KM # Radius in m.
MU = G*MS      # Gravitational parameter with respect to the Sun

#---------------------------------------------------------------------------------

def bcsFromAPm(A,P,m):
    """A - aphelion in m, P - perihelion in m, m - mass in kg.
       Returns semi-minor axis, linear eccentricity and max and min speed.
    """
    a = (A + P) / 2 # semi-major axis in m.
    b = sqrt(A * P) # semi-minor axis in m.
    c = (A - P) / 2 # linear eccentricity = center-to-focus distance, in m.
    mu = G * (MS + m) # gravitational parameter (with respect to the Sun)
    sMax = sqrt(mu * (2/P - 1/a)) # speed at perihelion - max speed
    sMin = sqrt(mu * (2/A - 1/a)) # speed at aphelion - min speed
    # e = c / a; # eccentricity.
    # T = sqrt(4 * pi * pi * a * a * a / mu) # orbital period.
    return (b, c, sMax, sMin)
    
#---------------------------------------------------------------------------------

MERCURY_DATA = {}
MERCURY_DATA["mass"] = 3.3011e23
MERCURY_DATA["aphelion"]        = 69_820_000*KM 
MERCURY_DATA["perihelion"]      = 46_000_000*KM 
MERCURY_DATA["semi-major axis"] = 57_910_000*KM 
MERCURY_DATA["eccentricity"] = 0.205630      
MERCURY_DATA["avg speed"] = 47.36*KM/SEC  
MERCURY_DATA["sidereal orbital period"] = 87.9691*DAY
MERCURY_DATA["radius"] = 2439.7*KM

bcs = bcsFromAPm(MERCURY_DATA["aphelion"],
                 MERCURY_DATA["perihelion"],
                 MERCURY_DATA["mass"])

MERCURY_DATA["semi-minor axis"] = bcs[0]
MERCURY_DATA["linear eccentricity"] = bcs[1] # center-to-focus distance
MERCURY_DATA["max speed"] = bcs[2] # at the perihelion
MERCURY_DATA["min speed"] = bcs[3] # at the aphelion

#---------------------------------------------------------------------------------

VENUS_DATA = {}
VENUS_DATA["mass"] = 4.8675e24
VENUS_DATA["aphelion"]        = 108_940_000*KM
VENUS_DATA["perihelion"]      = 107_480_000*KM
VENUS_DATA["semi-major axis"] = 108_210_000*KM
VENUS_DATA["eccentricity"] = 0.006772   
VENUS_DATA["avg speed"] = 35.02*KM/SEC 
VENUS_DATA["sidereal orbital period"] = 224.701*DAY 
VENUS_DATA["radius"] = 6051.8*KM

bcs = bcsFromAPm(VENUS_DATA["aphelion"],
                 VENUS_DATA["perihelion"],
                 VENUS_DATA["mass"])

VENUS_DATA["semi-minor axis"] = bcs[0]
VENUS_DATA["linear eccentricity"] = bcs[1] # center-to-focus distance
VENUS_DATA["max speed"] = bcs[2] # at the perihelion
VENUS_DATA["min speed"] = bcs[3] # at the aphelion

#---------------------------------------------------------------------------------

EARTH_DATA = {}
EARTH_DATA["mass"] = 5.972168e24
EARTH_DATA["aphelion"]        = 152_097_597*KM
EARTH_DATA["perihelion"]      = 147_098_450*KM
EARTH_DATA["semi-major axis"] = 149_598_023*KM
EARTH_DATA["eccentricity"] = 0.0167086
EARTH_DATA["avg speed"] = 29.7827*KM/SEC
EARTH_DATA["sidereal orbital period"] = 365.2563630040*DAY 
EARTH_DATA["radius"] = 6371*KM

bcs = bcsFromAPm(EARTH_DATA["aphelion"],
                 EARTH_DATA["perihelion"],
                 EARTH_DATA["mass"])

EARTH_DATA["semi-minor axis"] = bcs[0]
EARTH_DATA["linear eccentricity"] = bcs[1] # center-to-focus distance
EARTH_DATA["max speed"] = bcs[2] # at the perihelion
EARTH_DATA["min speed"] = bcs[3] # at the aphelion

#---------------------------------------------------------------------------------

MARS_DATA = {}
MARS_DATA["mass"] = 6.4171e23
MARS_DATA["aphelion"]        = 249_261_000*KM
MARS_DATA["perihelion"]      = 206_650_000*KM
MARS_DATA["semi-major axis"] = 227_939_366*KM
MARS_DATA["eccentricity"] = 0.0934
MARS_DATA["avg speed"] = 24.07*KM/SEC
MARS_DATA["sidereal orbital period"] = 686.980*DAY
MARS_DATA["radius"] = 3389.5*KM

bcs = bcsFromAPm(MARS_DATA["aphelion"],
                 MARS_DATA["perihelion"],
                 MARS_DATA["mass"])

MARS_DATA["semi-minor axis"] = bcs[0]
MARS_DATA["linear eccentricity"] = bcs[1] # center-to-focus distance
MARS_DATA["max speed"] = bcs[2] # at the perihelion
MARS_DATA["min speed"] = bcs[3] # at the aphelion

#---------------------------------------------------------------------------------

# requires a bigger scaling constant in simulations!

JUPITER_DATA = {}
JUPITER_DATA["mass"] = 1.8982e27
JUPITER_DATA["aphelion"]        = 816_363_000*KM
JUPITER_DATA["perihelion"]      = 740_595_000*KM
JUPITER_DATA["semi-major axis"] = 778_479_000*KM
JUPITER_DATA["eccentricity"] = 0.0489
JUPITER_DATA["avg speed"] = 13.07*KM/SEC
JUPITER_DATA["sidereal orbital period"] = 4_332.59*DAY 
JUPITER_DATA["color"] = (176,127,53) # Mandalay -yellowish brown, "0xB07F35"
JUPITER_DATA["radius"] = 69_911*KM

bcs = bcsFromAPm(JUPITER_DATA["aphelion"],
                 JUPITER_DATA["perihelion"],
                 JUPITER_DATA["mass"])

JUPITER_DATA["semi-minor axis"] = bcs[0]
JUPITER_DATA["linear eccentricity"] = bcs[1] # center-to-focus distance
JUPITER_DATA["max speed"] = bcs[2] # at the perihelion
JUPITER_DATA["min speed"] = bcs[3] # at the aphelion

#---------------------------------------------------------------------------------

# requires a bigger scaling constant in simulations!

SATURN_DATA = {}
SATURN_DATA["mass"] = 5.6834e26
SATURN_DATA["aphelion"]        = 1_514_500_000*KM
SATURN_DATA["perihelion"]      = 1_352_550_000*KM
SATURN_DATA["semi-major axis"] = 1_433_530_000*KM
SATURN_DATA["eccentricity"] = 0.0565
SATURN_DATA["avg speed"] = 9.68*KM/SEC
SATURN_DATA["sidereal orbital period"] = 10_755.70*DAY 
SATURN_DATA["radius"] = 58_232*KM

bcs = bcsFromAPm(SATURN_DATA["aphelion"],
                 SATURN_DATA["perihelion"],
                 SATURN_DATA["mass"])

SATURN_DATA["semi-minor axis"] = bcs[0]
SATURN_DATA["linear eccentricity"] = bcs[1] # center-to-focus distance
SATURN_DATA["max speed"] = bcs[2] # at the perihelion
SATURN_DATA["min speed"] = bcs[3] # at the aphelion

#---------------------------------------------------------------------------------

# requires a bigger scaling constant in simulations!

URANUS_DATA = {}
URANUS_DATA["mass"] = 8.6810e25
URANUS_DATA["aphelion"]        = 3_006_390_000*KM
URANUS_DATA["perihelion"]      = 2_735_560_000*KM
URANUS_DATA["semi-major axis"] = 2_870_972_000*KM
URANUS_DATA["eccentricity"] = 0.04717
URANUS_DATA["avg speed"] = 6.80*KM/SEC
URANUS_DATA["sidereal orbital period"] = 30_688.5*DAY
URANUS_DATA["radius"] = 25_362*KM

bcs = bcsFromAPm(URANUS_DATA["aphelion"],
                 URANUS_DATA["perihelion"],
                 URANUS_DATA["mass"])

URANUS_DATA["semi-minor axis"] = bcs[0]
URANUS_DATA["linear eccentricity"] = bcs[1] # center-to-focus distance
URANUS_DATA["max speed"] = bcs[2] # at the perihelion
URANUS_DATA["min speed"] = bcs[3] # at the aphelion

#---------------------------------------------------------------------------------

# requires a bigger scaling constant in simulations!

NEPTUNE_DATA = {}
NEPTUNE_DATA["mass"] = 1.02409e26
NEPTUNE_DATA["aphelion"]        = 4_540_000_000*KM
NEPTUNE_DATA["perihelion"]      = 4_460_000_000*KM
NEPTUNE_DATA["semi-major axis"] = 4_500_000_000*KM
NEPTUNE_DATA["eccentricity"] = 0.008678
NEPTUNE_DATA["avg speed"] = 5.43*KM/SEC
NEPTUNE_DATA["sidereal orbital period"] = 60_195*DAY 
NEPTUNE_DATA["radius"] = 24_622*KM

bcs = bcsFromAPm(NEPTUNE_DATA["aphelion"],
                 NEPTUNE_DATA["perihelion"],
                 NEPTUNE_DATA["mass"])

NEPTUNE_DATA["semi-minor axis"] = bcs[0]
NEPTUNE_DATA["linear eccentricity"] = bcs[1] # center-to-focus distance
NEPTUNE_DATA["max speed"] = bcs[2] # at the perihelion
NEPTUNE_DATA["min speed"] = bcs[3] # at the aphelion

#---------------------------------------------------------------------------------

# requires a bigger scaling constant in simulations!

PLUTO_DATA = {}
PLUTO_DATA["mass"] = 1.303e22
PLUTO_DATA["aphelion"]        = 7_375_930_000*KM
PLUTO_DATA["perihelion"]      = 4_436_820_000*KM
PLUTO_DATA["semi-major axis"] = 5_906_380_000*KM
PLUTO_DATA["eccentricity"] = 0.2488
PLUTO_DATA["avg speed"] = 4.743*KM/SEC
PLUTO_DATA["sidereal orbital period"] = 90_560*DAY 
PLUTO_DATA["radius"] = 1_188.3*KM

bcs = bcsFromAPm(PLUTO_DATA["aphelion"],
                 PLUTO_DATA["perihelion"],
                 PLUTO_DATA["mass"])

PLUTO_DATA["semi-minor axis"] = bcs[0]
PLUTO_DATA["linear eccentricity"] = bcs[1] # center-to-focus distance
PLUTO_DATA["max speed"] = bcs[2] # at the perihelion
PLUTO_DATA["min speed"] = bcs[3] # at the aphelion

# While the planets' orbits are close to the ecliptic,
# Pluto's orbit has a significant inclination.
# This program shows orbits on a single plane. 
# Showing the planets and Pluto's orbits on the same plane would be
# greatly inaccurate and lead to misconceptions, for instance
# the visualization would wshow that Pluto's orbit intersects Uranus' orbit.
#---------------------------------------------------------------------------------

# requires a bigger scaling constant in simulations!

# Halley's Comet, unlike the planets, moves clockwise around the Sun
# when viewed form the north pole.

HALLEY_COMET_DATA = {}
HALLEY_COMET_DATA["mass"] = 2.2e14
HALLEY_COMET_DATA["aphelion"]        = 35.14*AU
HALLEY_COMET_DATA["perihelion"]      =  0.59278*AU
HALLEY_COMET_DATA["semi-major axis"] = 17.737*AU
HALLEY_COMET_DATA["eccentricity"] = 0.96658 
# HALLEY_COMET_DATA["avg speed"] = 
HALLEY_COMET_DATA["sidereal orbital period"] =  74.7*YR
HALLEY_COMET_DATA["radius"] = 5.5*KM

bcs = bcsFromAPm(HALLEY_COMET_DATA["aphelion"],
                 HALLEY_COMET_DATA["perihelion"],
                 HALLEY_COMET_DATA["mass"])

HALLEY_COMET_DATA["semi-minor axis"] = bcs[0]
HALLEY_COMET_DATA["linear eccentricity"] = bcs[1] # center-to-focus distance
HALLEY_COMET_DATA["max speed"] = bcs[2] # at the perihelion
HALLEY_COMET_DATA["min speed"] = bcs[3] # at the aphelion

# Halley's comet has a significant inclination with respect to the excliptic.
# and is also far from the plane of Pluto's orbit.

#=================================================================================
# PLANETS

class Planet(object):

    def __init__(self, name, mass, perihelionDistance, maxSpeed,
                 color="green", timeStep : int = 1000):
        """mass in kg, perihelionDistance in m, maxSpeed in m/s,
           simulation timeStep in s.
           Assumes that the Sun is at (0, 0) in a coordinate system.
           Puts the planet/self at coordinates (perihelionDistance, 0)
           giving it velocity (0, maxSpeed).
        """
        self._mass = mass
        self._x = perihelionDistance
        self._y = 0
        self._r2 = self._x*self._x + self._y*self._y # radius squared
        # we store radius squared to avoid loosing precision while recalculating it
        self._r = sqrt(self._r2) # radius = distance from Sun's center.
        self._vx = 0 # horizontal speed in m/s 
        self._vy = maxSpeed # vertical speed in m/s
        self._ax = -MU*self._x/(self._r2*self._r) # horizontal acceleration m/s^2
        self._ay = -MU*self._y/(self._r2*self._r) # vertical acceleration m/s^2
        self._timeStep = timeStep
        self._name = name
        self._color = color

    def move(self, timeStep="the default for this planet"):
        """Updates position, velocity, acceleration, radius
           to those after timeStep.
           timeStep is an optional parameter.
           If timeStep is is not provided, uses the default value
           with which the Planet object was created.
        """
        if timeStep == "the default for this planet":
            timeStep = self._timeStep
        self._x += self._vx * timeStep
        self._y += self._vy * timeStep
        self._r2 = self._x*self._x + self._y*self._y
        self._r = sqrt(self._r2)
        self._vx += self._ax * timeStep
        self._vy += self._ay * timeStep
        self._ax = -MU*self._x/(self._r2*self._r) # horizontal acceleration m/s^2
        self._ay = -MU*self._y/(self._r2*self._r) # vertical acceleration m/s^2

    def position(self, scaleDownFactor: int = 1):
        return (self._x/scaleDownFactor, self._y/scaleDownFactor)

    def velocity(self):
        return (self._vx, self._vy)

    def acceleration(self):
        return (self._ax, self._ay)

    def mass(self):
        return self._mass

    def radius(self):
        return self._r

    def radiusSquared(self):
        return self._r2

    def timeStep(self):
        return self._timeStep

    def name(self):
        return(self._name)

    def color(self):
        return(self._color)

    def setColor(self, color: str): 
        self._color = color

    def timeStep(self):
        return self._timeStep

    def setTimeStep(self, timeStep : int):
        self._timeStep = timeStep

#---------------------------------------------------------------------------------
# Simulated planets - global constants.

# The Sun is at (0,0).
# The starting position for a planet is (perihilionDistance, 0).
# The starting velocity is (0, sMax).
# The default time step in the simulation: 1000s.

# Inner planets - rocky:
MERCURY = Planet("Mercury",
                 MERCURY_DATA["mass"],
                 MERCURY_DATA["perihelion"],
                 MERCURY_DATA["max speed"],
                 "grey"
                )
VENUS   = Planet("Venus",
                 VENUS_DATA["mass"],
                 VENUS_DATA["perihelion"],
                 VENUS_DATA["max speed"],
                 "white"
                )
EARTH   = Planet("Earth",
                 EARTH_DATA["mass"],
                 EARTH_DATA["perihelion"],
                 EARTH_DATA["max speed"],
                 "turquoise"
                )
MARS    = Planet("Mars",
                 MARS_DATA["mass"],
                 MARS_DATA["perihelion"],
                 MARS_DATA["max speed"],
                 "red"
                )
# Outer planets - gas giants:
JUPITER = Planet("Jupiter",
                 JUPITER_DATA["mass"],
                 JUPITER_DATA["perihelion"],
                 JUPITER_DATA["max speed"],
                 "orange"
                )
SATURN  = Planet("Saturn",
                 SATURN_DATA["mass"],
                 SATURN_DATA["perihelion"],
                 SATURN_DATA["max speed"],
                 "gold"
                )
URANUS  = Planet("Uranus",
                 URANUS_DATA["mass"],
                 URANUS_DATA["perihelion"],
                 URANUS_DATA["max speed"],
                 "lightBlue1"
                )
NEPTUNE = Planet("Neptune",
                 NEPTUNE_DATA["mass"],
                 NEPTUNE_DATA["perihelion"],
                 NEPTUNE_DATA["max speed"],
                 "lightBlue2"
                )

ME = EARTH_DATA["mass"]
PE = EARTH_DATA["perihelion"]
SE = EARTH_DATA["max speed"]
S10 = sqrt(G*(MS+ME)/PE) # the speed of a made-up planet
# with the same mass and perihelion as Earth, but with a circular orbit.

# Made-up planets, for computational experiments (default pencolor=green)
PLANET07 = Planet("Planet 0.7", ME, PE, 0.7*S10) 
PLANET08 = Planet("Planet 0.8", ME, PE, 0.8*S10)
PLANET09 = Planet("Planet 0.9", ME, PE, 0.9*S10)
PLANET10 = Planet("Planet 1.0", ME, PE,     S10)
PLANET11 = Planet("Planet 1.1", ME, PE, 1.1*S10)
PLANET12 = Planet("Planet 1.2", ME, PE, 1.2*S10)
PLANET13 = Planet("Planet 1.3", ME, PE, 1.3*S10)

# Note.
# Inner planets and made up planets above are close to the Sun,
# and outer ones - are at vast distances.

#=================================================================================
# AUXILIARY FUNCTIONS

def sky(skyColor="black", showSun=True):
    """Create black canvas with the white Sun at (0,0).
       The Sun size is not to scale; it is shown much bigger."
    """
    screen = getscreen()
    screen.clear() # remove turtle image
    screen.screensize(1000,1000)
    screen.title("Kepler's world")
    screen.bgcolor(skyColor)
    if showSun:
        t = Turtle(visible=False)
        t.dot(10, "white") # SUN at (0,0) and is white.
        print("The Sun is not to scale; it is shown much bigger.")

#sky()

#---------------------------------------------------------------------------------
   
def drawEllipse(semiMajorAxis, semiMinorAxis, leftShift=0,
                color="white", focusColor="white"):
    """Draws an ellipse centered at (-leftShift,0) and shows the foci;
       the foci are on the x-axis.
       For an ellipse centered at (0,0) use leftShift=0.
       For an ellipse with the right focus at (0,0), use leftShift = c,
       where c is the linear eccentricity i.e. center-to-focus distance:
       c = sqrt(semiMajorAxis*semiMajorAxis - semiMinorAxis*semiMinorAxis)
       Note: make sure to create canvas before this function is called.
       Note: the parameters are in turtle canvas units, not in meters.
    """
    if semiMajorAxis < semiMinorAxis:
        raise ValueError(
            "semi-major axis cannot be smaller than semi-minor axis.")
    t = Turtle(visible=False)
    t.speed("fastest")
    t.pendown()
    t.pensize(1)
    t.pencolor(color)
    c = sqrt(semiMajorAxis*semiMajorAxis - semiMinorAxis*semiMinorAxis)
    # c is the linear eccentricity, i.e. center-to-focus distance.
    t.teleport(c-leftShift,0) # left focus
    t.dot(6, focusColor) 
    t.teleport(-c-leftShift,0) # right focus
    t.dot(6, focusColor) 
    t.teleport(semiMajorAxis-leftShift,0) # left vertex of the ellipse
    n = 200
    for i in range(n+1):
        angle = 2*pi*(i/n)
        x = semiMajorAxis*cos(angle)
        y = semiMinorAxis*sin(angle)
        t.goto(x-leftShift,y)

#sky(showSun=False) # run/uncomment this before running drawEllipse
#drawEllipse(200, 100, 0, "pink")
     
#---------------------------------------------------------------------------------

def simulate(planet: Planet, scaleDownFactor =  1_000_000_000):
    # Just draws the orbit, does not test Kepler's laws.
    # This function is not used by the top level functions in the program.
    # It is given here as a stepping stone to understand simulateAndTest below.
    """Precondition: planet position (x,y) must have x>0, y=0,
                     and velocity (vx,vy) must have vx=0, vy>0.
       So, the planet must be in the right vertex of its elliptical orbit.
       This function calculates the orbit resulting from the continuing local
       effect of the Newton's law of gravity, assuming the Sun is at (0,0).
       Every step in the simulation is done for the actual planet
       (such as Mars, with its actual mass, perihelion and actual max speed)
       and only later it is scaled down to be displayed in turtle graphics.
       Real distances in meters will be divided by scaleDownFactor
       before being given to the turtle.
       Note: make sure to create canvas before this function is called.
    """
    # Prepare turtle.
    t = Turtle(visible=False)
    t.speed("fastest")
    t.pendown()
    t.pensize(3)
    t.pencolor(planet.color())
    t.teleport(*planet.position(scaleDownFactor))

    # The planet starts from its perihelion.
    
    # Upper half of the orbit:
    done = False
    while True: 
        for j in range(100): # update turtle every 100 moves
            planet.move()
            if planet.position()[1]<0: # if y<0
                done = True
                break
        t.goto(*planet.position(scaleDownFactor))
        if done: break

    # The planet is now at its aphelion.

    # Lower half of the orbit:
    done = False
    while True: 
        for j in range(100): # update turtle every 100 moves
            planet.move()
            if planet.position()[1]>=0: # if y>=0
                done = True
                break
        t.goto(*planet.position(scaleDownFactor))
        if done: break

    # The planet is back at the perihelion.

#sky() # run/uncomment this before running simulate!
#simulate(PLANET12)

#-----------------------------------------------------------------------------

# Under construction
def simulateAndTest(planet: Planet, scaleDownFactor):
    """Precondition: planet position (x,y) must have x>0, y=0,
                     and velocity (vx,vy) must have vx=0, vy>0.
       So, the planet must be in the right vertex of its elliptical orbit.
       This function calculates the orbit resulting from the continuing local
       effect of the Newton's law of gravity, assuming the Sun is at (0,0).
       Every step in the simulation is done for the actual planet
       (such as Mars, with its actual mass, perihelion and actual max speed)
       and only later it is scaled down to be displayed in turtle graphics.
       Real distances in meters will be divided by scaleDownFactor
       before being given to the turtle.
       Tests Kepler's laws ...
       Returns the orbital period in seconds, ...
       Note: make sure to create canvas before this function is called.
    """
    # Prepare turtle t to draw an orbit
    t = Turtle(visible=False)
    t.speed("fastest")
    t.pendown()
    t.pensize(3)
    t.pencolor(planet.color())
    t.teleport(*planet.position(scaleDownFactor))

    # Prepare turtle t2 to draw a graph of area speed. 
    t2 = Turtle(visible=False)
    t2.speed("fastest")
    t2.pendown()
    t2.pensize(1)
    t2.pencolor(planet.color())
    

    # Concerning Kepler's 1st law
    x,y = planet.position() # perihelion,
                            # the right vertex of the elliptical orbit 
    maxX = x # perihelion distance from the Sun (Sun is at (0,0))
    minXsoFar = maxX            
    maxYsoFar = y
    # These will be used to find the semi-major axis and semi-minor axis.

    # Concerning Kepler's 2nd law
    vx,vy = planet.velocity()

    areaScaleDownFactor = scaleDownFactor*vy
    # area/areaScaleDownFactor ~ radius/scaleDownFactor
    
    area0 = x*vy-y*vx # The determinant of the matrix of column vectors r,v =
                      # = vector cross product  r x v.
                      # = area of the parallelogram spanned by vectors r,v =
                      # = twice the area of a triangle spanned by vectors r,v.
                      # = twice the area speed in (m^2)/s - area swept per sec.
                      # Notice that angular momentum is  r x mv.

    #t2.teleport(x/scaleDownFactor, area0/areaScaleDownFactor/2)

    minAreaSoFar = area0
    maxAreaSoFar = area0
    # If the difference between minArea and maxArea is small,
    # Kepler's 3rd law will be confirmed.

    # Concerning Kepler's 3rd law
    TsoFar = 0 # sidereal orbital period in seconds.

    # The planet starts from its perihelion.

    # Upper half of the orbit:
    done = False 
    while True: 
        for j in range(100): # update turtle every 100 moves
            planet.move()
            TsoFar += planet.timeStep() # update T
            x, y = planet.position() # update minX, maxY
            if x < minXsoFar: minXsoFar = x 
            if y > maxYsoFar: maxYsoFar = y
            vx, vy = planet.velocity() # update minArea, maxArea
            area = x*vy-y*vx
            if area < minAreaSoFar: minAreaSoFar = area 
            if area > maxAreaSoFar: maxAreaSoFar = area
            if y < 0:
                done = True
                break
        t.goto(*planet.position(scaleDownFactor))
        #t2.goto(x/scaleDownFactor, area/areaScaleDownFactor/2)
        if done: break

    # The planet is now at its aphelion.

    area1= area

    # Lower half of the orbit:
    done = False 
    while True: 
        for j in range(100): # update turtle every 100 moves
            planet.move()
            TsoFar += planet.timeStep() # update T
            x,y = planet.position()
            vx,vy = planet.velocity() # update minArea, maxArea
            area = x*vy-y*vx
            if area < minAreaSoFar: minAreaSoFar = area 
            if area > maxAreaSoFar: maxAreaSoFar = area
            if y >= 0:
                done = True
                break
        t.goto(*planet.position(scaleDownFactor))
        if done: break

        area2 = area

    # The planet is back at the perihelion.

    # Concerning Kepler's 1st law:
    minX = minXsoFar 
    maxY = maxYsoFar
    a = (maxX - minX)/2 # semi-major axis - from the simulation
    b = maxY # semi-minor axis - from the simulation
    # Notice that maxX is the distance from (0,0) to the perihelion.
    orbitCenter = maxX - a
    c = sqrt(a*a - b*b) # linear eccentricity = center to focus distance.
    t.teleport((orbitCenter + c)/scaleDownFactor,0) # draw right focus
    t.dot(6, planet.color())
    t.teleport((orbitCenter - c)/scaleDownFactor,0) # draw left focus
    t.dot(6, planet.color()) 

    # Concerning Kepler's 2nd law: 
    minArea = minAreaSoFar
    maxArea = maxAreaSoFar
    kepler2discrepancy = (maxArea - minArea) / area0
    kepler2discrepancyPercent = round(kepler2discrepancy*100, 2)
    print(kepler2discrepancyPercent, "%")
    avgArea = (minArea + maxArea)/2

    #print()
    #print("Area:")
    #print(minArea)
    #print(area0, 0) # starting area (perihelion) Why is minArea=area0 ?
    #print(area1, 1) # aphelion
    #print(area2, 2) # ending area (perihelion)
    #print(avgArea)
    #print(maxArea)
    #print("====")

    # Concerning Kepler's 3rd law:
    T = TsoFar # sidereal orbital period - from the simulation.
    # T^2 / a^2 = 4*pi^2 / G(M+m) - does the simulation support this?
    #print("Kepler3:")
    lhs = T*T/(a*a*a)
    #print(lhs)
    rhs = 4*pi*pi / (G*(MS+planet._mass))
    #print(rhs)    
    #print(abs((lhs-rhs)/rhs))
    #print(round(abs((lhs-rhs)/rhs)*100,2), "%")

    #print(area0/areaScaleDownFactor)
    #print(area1/areaScaleDownFactor)
    
    return a,b,kepler2discrepancyPercent,T

#sky() # run/uncomment this before running simulateAndTest!    
#simulateAndTest(PLANET12)
#simulate(PLANET11)
#simulate(PLANET10)
#simulate(PLANET09)
#simulate(PLANET08)
#simulate(PLANET07)

#=================================================================================
# FUNCTIONS CALLED BY main

# Under construction
def simmulationSummary(planetData: dict, planet: Planet, scaleDownFactor):
    print("\n")
    print(planet.name(), "orbital period in days:")
    T = planetData["sidereal orbital period"]
    print(round(T/DAY,2), "- actual")
    drawEllipse(planetData["semi-major axis"]/scaleDownFactor,
                planetData["semi-minor axis"]/scaleDownFactor,
                planetData["linear eccentricity"]/scaleDownFactor)
    TS = simulateAndTest(planet, scaleDownFactor)[3]
    print(round(abs(TS-T)*100/T, 2), "% error")

def planets(n: int = 4):
    """A computer simulation of orbits of n planets, resulting from the
       continuing local effect of the Newton's law of gravity.
       Tests if the simulated planets obey (global) Kepler's laws 1 and 3.
       Every step in the simulation is done for the actual planet
       (such as Mars, with its actual mass, perihelion and actual max speed)
       and only later it is scaled down to be displayed in turtle graphics.
       This function creates its own canvas-sky with the Sun,
       and draws an ellipse of the actual orbit with the foci,
       before starting the simulation.
       Besides displaying turtle graphics, it produces output.
    """
    print("""
The orbits displayed in turtle graphics are to scale.
and planets move counter-clockwise as seen from the north pole.
Other characteristics of the orbits are not modelled here -
all ellipses are shown in the same plane with the major axis on the x-axis
and perturbations due to gravitational influence of other planets is not shown.
         """)
    if n <= 5:
        scaleDownFactor =  1_000_000_000
    else:
        scaleDownFactor = 10_000_000_000
    # Real distances in meters will be divided by scaleDownFactor
    # before being given to the turtle.
    
    sky()

    if n>=9:
        print(
"""There are only 8 known planets in the Solar System.
Pluto is the biggest known dwarf planet
and the biggest known trans-Neptunian object.
It is bigger than all the asteroids in the belt between Mars and Jupiter,
but smaller than the Moon. This simulation does not include Pluto."""
             )
    
    print("\nTesting Kepler's 1st and 3rd laws")
    
    # Default time step = 1000 s
    if n>=1:
        simmulationSummary(MERCURY_DATA, MERCURY, scaleDownFactor)
    if n>=2:
        simmulationSummary(VENUS_DATA, VENUS, scaleDownFactor)
    if n>=3:
        simmulationSummary(EARTH_DATA, EARTH, scaleDownFactor)
    if n>=4:
        simmulationSummary(MARS_DATA, MARS, scaleDownFactor)
    if n==5:
        simmulationSummary(JUPITER_DATA, JUPITER, scaleDownFactor)
    if n>5:
        JUPITER.setTimeStep(10_000)
        simmulationSummary(JUPITER_DATA, JUPITER, scaleDownFactor)
    if n>=6:
        SATURN.setTimeStep(10_000)
        simmulationSummary(SATURN_DATA, SATURN, scaleDownFactor)
    if n>=7:
        URANUS.setTimeStep(10_000)
        simmulationSummary(URANUS_DATA, URANUS, scaleDownFactor)
    if n>=8:
        NEPTUNE.setTimeStep(10_000)
        simmulationSummary(NEPTUNE_DATA, NEPTUNE, scaleDownFactor)

#-------
#innerPlanets() # uncomment this to simulate the inner planets.

# Output:

# The orbits are to scale.
# and planets move counter-clockwise as seen from the north pole.
# The orientation of the major axes of orbits is not modeled here:
# all ellipses are shown with the major axis on the x-axis
# and the Sun in the right focus.
# The Sun is not to scale; shown much bigger.
# 
# Testing Kepler's 1st and 3rd laws
# 
# Mercury orbital period in days:
# 87.97 - actual
# 88.77 - from the simulation.
# 0.91 % error
# 
# Venus - orbital period in days:
# 224.7 - actual
# 225.46 - from the simulation.
# 0.34 % error
# 
# Earth - orbital period in days:
# 365.26 - actual
# 365.97 - from the simulation.
# 0.2 % error
# 
# Mars - orbital period in days:
# 686.98 - actual
# 687.85 - from the simulation.
# 0.13 % error
   
#---------------------------------------------------------------------------------

# Under construction
def testKepler(planet: Planet, scaleDownFactor = 1_000_000_000):
    """Precondition: planet position (x,y) must have x>0, y=0,
                    and velocity (vx,vy) must have vx=0, vy>0.
       So, the planet must be in the right vertex of its elliptical orbit.A computer simulation of the orbit of the planet, resulting from the
       continuing local effect of the Newton's law of gravity.
       Tests if the simulated planet obey (global) Kepler's laws 1 and 3.
       Every step in the simulation is done for the actual planet
       (such as Mars, with its actual mass, perihelion and actual max speed)
       and only later it is scaled down to be displayed in turtle graphics.
       Real distances in meters will be divided by scaleDownFactor
       before being given to the turtle.
       This function creates its own canvas-sky with the Sun,
       and draws the predicted ellipse of the orbit with the foci,
       before starting simulation.
       It produces output, besides displaying turtle graphics.
    """
    print("The orbit displayed in turtle graphics is to scale.")
    sky()
    print("\nTesting Kepler's 1st and 3rd laws")
    m = planet.mass()  # mass
    vmax = planet.velocity()[1] # maximal speed (at Perihelion)
    P = planet.position()[0]  # Perihelion distance (shortest distance from Sun)
    mu = G*(MS+m) # gravitational parameter
    a = P*mu / (2*mu-P*vmax*vmax) # semi-major axis:
    c = a - P # the linear eccentricity, i.e. center-to-focus distance
    A = a + c # Aphelion distance (biggest distance from the Sun)
    b = sqrt(A*P)  # semi-minor axis
    T = sqrt(4*pi*pi*a*a*a / mu) # orbital period in s
    print("Planet's orbital period in days:")
    print(round(T/DAY,2), "- predicted by the theory")
    drawEllipse(a/scaleDownFactor, b/scaleDownFactor, c/scaleDownFactor) # predicted orbit
    TS = simulateAndTest(planet, scaleDownFactor)[3] # simulation. TS - orbital period form simulation.
    print(round(TS/DAY,2), "- from the simulation")
    print( round(abs(TS-T)*100/T, 2), "% discrepancy")

#testKepler(PLANET12)

#=================================================================================

def main():
    print("Kepler's World")

    print("""
This program shows the orbit predicted by a formula,
then simulates the celestial body's movement and tests Kepler's laws.
In the case of actual celestial bodies
it also compares the simulated orbit to the data from astronomical tables."""
          )

    while True:
        print("""
1. Inner planets: Mercury, Venus, Earth and Mars.
   
2. Eight planets: Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune.
   
3. A made-up celestial body with the same mass and perihelion as Earth
   but with the maximal speed 20% bigger than Earth.

4. Exit.
          """)
        choice = input("Enter your choice (1-4): ")
        if choice == "":
            return
        if   choice[0] == "1":
            planets(4) # The 4 inner planets.
        elif choice[0] == "2":
            planets(8) # All 8 planets.
        elif choice[0] == "3":
            testKepler(PLANET12)
        else:
            print("Bye")
            return

if __name__ == "__main__": main()

#=================================================================================

#sky()
#simulateAndTest(PLANET12)
#innerPlanets()
#simulateAndTest(JUPITER)
#simulateAndTest(SATURN)
#simulateAndTest(URANUS)
#simulateAndTest(NEPTUNE)
