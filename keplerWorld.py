"""
A simple simulation of movements of planets around the Sun.
Uses a discrete approximation of the real-world continuous
application of Newton's law of gravity (F=GMm/rr).
The calculations are done for actual planets,
and results are later scaled down to be displayed in turtle graphics.
The simulation tests if the resulting orbits and planets' movements conform
to the Kepler laws, and to their enhancements, such as the Vis Viva Equation.
This provides an example of continuing local effect (of the law of gravity),
having a global effect (of elliptical orbits with details specified by Kepler.) 
"""

#==================================================================================
from turtle import getscreen, Turtle
from math import sqrt, pi, sin, cos

#==================================================================================
# ASTRONOMICAL DATA FROM OBSERVATIONS AND CALCUATIONS 

# Astronomical data in kilograms, meters, seconds, from Wikipedia.

#---------------------------------------------------------------------------------

G = 6.67430e-11 # m^3 / (kg * s^2). Gravitational constant
DAY = 24*60*60  # s
KM = 1000 # m

#---------------------------------------------------------------------------------

def bcsFromAPm(A,P,m):
    """A - aphelion in m, P - perihelion in m, m - mass in kg."""
    a = (A + P) / 2 # semi-major axis in m.
    b = sqrt(A * P) # semi-minor axis in m.
    c = (A - P) / 2 # Linear eccentricity = center-to-focus distance, in m.
    mu = G * (M_S + m) # Gravitational parameter.
    sMax = sqrt(mu * (2/P - 1/a)) # speed at perihelion (max speed)
    # e = c / a; # eccentricity.
    # T = sqrt(4 * pi * pi * a * a * a / mu) # Orbital period.
    return b, c, sMax
    
#print(bcsFromAPm(A_Me, P_Me, m_Me))
#print(bcsFromAPm(A_V,  P_V,  m_V ))
#print(bcsFromAPm(A_E,  P_E,  m_E ))
#print(bcsFromAPm(A_Ma, P_Ma, m_Ma))

#---------------------------------------------------------------------------------
# Sun

# white
M_S = 1.98847e30 # Mass in kg.
R_S = 696000*KM  # Radius in m.

#---------------------------------------------------------------------------------
# Mercury

# grey
A_Me = 69820000*KM # Aphelion in m.
P_Me = 46000000*KM # Perihlion in m.
a_Me = 57910000*KM # Semi-major axis in m.
e_Me = 0.205630    # Eccentricity.
T_Me = 87.9691*DAY # Sidereal orbital period in s.
s_Me = 47360       # Average orbital speed in m/s.
m_Me = 3.3011e23   # Mass in kg.

bcs = bcsFromAPm(A_Me, P_Me, m_Me)

b_Me    = bcs[0]   # semi-minor axis in m.
c_Me    = bcs[1]   # linear eccentricity = center-to-focus distance, in m.
sMax_Me = bcs[2]   # Max speed (at the prihelion) in m/s.

#---------------------------------------------------------------------------------
# Venus

# yellow/white
A_V = 108940000*KM # Aphelion in m.
P_V = 107480000*KM # Perihlion in m.
a_V = 108210000*KM # Semi-major axis in m.
e_V = 0.006772     # Eccentricity.
T_V = 224.701*DAY   # Sidereal orbital period in s.
s_V = 35020        # Average orbital speed in m/s.
m_V = 4.8675e24    # Mass in kg.

bcs = bcsFromAPm(A_V, P_V, m_V)

b_V    = bcs[0]    # semi-minor axis in m.
c_V    = bcs[1]    # linear eccentricity = center-to-focus distance, in m.
sMax_V = bcs[2]    # Max speed (at the prihelion) in m/s.

#---------------------------------------------------------------------------------
# Earth

# blue
A_E = 152097597*KM       # Aphelion in m.
P_E = 147098450*KM       # Perihlion in m.
a_E = 149598023*KM       # Semi-major axis in m.
e_E = 0.0167086          # Eccentricity
T_E = 365.2563630040*DAY # Sidereal orbital period in s. 
S_E = 29782.7            # Average orbital speed in m/s.
m_E = 5.972168e24        # Mass in kg.

bcs = bcsFromAPm(A_E, P_E, m_E)

b_E    = bcs[0]     # semi-minor axis in m.
c_E    = bcs[1]     # linear eccentricity = center-to-focus distance, in m.
sMax_E = bcs[2]     # Max speed (at the prihelion) in m/s.

#---------------------------------------------------------------------------------
# Mars

# red
A_Ma = 249261000*KM # Aphelion in m.
P_Ma = 206650000*KM # Perihlion in m.
a_Ma = 227939366*KM # Semi-major axis in m.
e_Ma = 0.0934       # Eccentricity.
T_Ma = 686.980*DAY  # Sidereal orbital period in s.
s_Ma = 24070        # Average orbital speed in m/s.
m_Ma = 6.4171e23    # Mass in kg.

bcs = bcsFromAPm(A_Ma, P_Ma, m_Ma)

b_Ma    = bcs[0]    # semi-minor axis in m.
c_Ma    = bcs[1]    # linear eccentricity = center-to-focus distance, in m.
sMax_Ma = bcs[2]    # Max speed (at the prihelion) in m/s.

#---------------------------------------------------------------------------------
# Jupiter

# requires a bigger scaling constant in simulations!

# yellow
A_J = 816363000*KM # Aphelion in m.
P_J = 740595000*KM # Perihlion in m.
a_J = 778479000*KM # Semi-major axis in m.
e_J = 0.0489       # Eccentricity.
T_J = 4332.59*DAY    # Sidereal orbital period in s.
s_J = 13070        # Average orbital speed in m/s.
m_J = 1.8982e27    # Mass in kg.

bcs = bcsFromAPm(A_J, P_J, m_J)

b_J    = bcs[0]    # semi-minor axis in m.
c_J    = bcs[1]    # linear eccentricity = center-to-focus distance, in m.
sMax_J = bcs[2]    # Max speed (at the prihelion) in m/s.

#---------------------------------------------------------------------------------
# Saturn

# requires a bigger scaling constant in simulations!

# yellow
A_S = 1514500000*KM # Aphelion in m.
P_S = 1352550000*KM # Perihlion in m.
a_S = 1433530000*KM # Semi-major axis in m.
e_S = 0.0565        # Eccentricity.
T_S = 10755.70*DAY  # Sidereal orbital period in s. 
s_S = 9680          # Average orbital speed in m/s.
m_S = 5.6834e26     # Mass in kg.

bcs = bcsFromAPm(A_S, P_S, m_S)

b_S    = bcs[0]    # semi-minor axis in m.
c_S    = bcs[1]    # linear eccentricity = center-to-focus distance, in m.
sMax_S = bcs[2]    # Max speed (at the prihelion) in m/s.

#---------------------------------------------------------------------------------
# Uranus

# requires a bigger scaling constant in simulations!

# white
A_U = 3006390000*KM # Aphelion in m.
P_U = 2735560000*KM # Perihlion in m.
a_U = 2870972000*KM # Semi-major axis in m.
e_U = 0.04717       # Eccentricity.
T_U = 30,688.5*DAY    # Sidereal orbital period in s.
s_U = 6800          # Average orbital speed in m/s.
m_U = 8.6810e25    # Mass in kg.

bcs = bcsFromAPm(A_U, P_U, m_U)

b_U    = bcs[0]    # semi-minor axis in m.
c_U    = bcs[1]    # linear eccentricity = center-to-focus distance, in m.
sMax_U = bcs[2]    # Max speed (at the prihelion) in m/s.

#---------------------------------------------------------------------------------
# Neptune

# requires a bigger scaling constant in simulations!

# blue
A_N = 4540000000*KM # Aphelion in m.
P_N = 4460000000*KM # Perihlion in m.
a_N = 4500000000*KM # Semi-major axis in m.
e_N = 0.008678      # Eccentricity.
T_N = 60195*DAY    # Sidereal orbital period in s.
s_N = 5.43          # Average orbital speed in m/s.
m_N = 1.02409e26    # Mass in kg.

bcs = bcsFromAPm(A_N, P_N, m_N)

b_N    = bcs[0]    # semi-minor axis in m.
c_N    = bcs[1]    # linear eccentricity = center-to-focus distance, in m.
sMax_N = bcs[2]    # Max speed (at the prihelion) in m/s.

#=================================================================================
# SIMULATION PARAMETERS

# A scaling factor suitable for showing orbits of the 4 inner planets.
# (outer planets require a factor 100 times bigger)
SCALING = 1000000000 # Real distances in meters will be divided by this
                     # before being given to the turtle.

TIME_STEP = 1000 # seconds.
        # Position, velocity, acceleration will be updated every TIME_STEP.

#=================================================================================
# PLANETS

class Planet(object):

    def __init__(self, mass, x, y, vx, vy, timeStep,
                 name=None, color="green", pensize=3):
        """mass in kg,
           x,y coodinates in m from the sun center,
           vx, vy - velocity in m/s,
           timeStep in s,
           size - size of the turtle pen while drawing the orbit.
        """
        self._mass = mass
        self._x = x
        self._y = y
        self._r2 = self._x*self._x + self._y*self._y # radius squared
        self._r = sqrt(self._r2) # radius = distance between sun and self.
        self._vx = vx # horizontal speed in m/s 
        self._vy = vy # vertical   speed in m/s
        self._ax = -G*M_S*x/(self._r2*self._r) # horizontal acceleration in m/s^2
        self._ay = -G*M_S*y/(self._r2*self._r) # vertical acceleration in m/s^2
        self._timeStep = timeStep
        self._name = name
        self._color = color
        self._pensize = pensize

    def move(self, timeStep=None):
        """Updates position, velocity, radius to those after timeStep.
           If timeStep is None, uses the default value
           with which the Planet object was created.
        """
        if timeStep is None: timeStep = self._timeStep
        self._x += self._vx * timeStep
        self._y += self._vy * timeStep
        self._r2 = self._x*self._x + self._y*self._y
        self._r = sqrt(self._r2)
        self._vx += self._ax * timeStep
        self._vy += self._ay * timeStep
        self._ax = -G*M_S*self._x/(self._r2*self._r) # horizontal acceleration in m / s^2
        self._ay = -G*M_S*self._y/(self._r2*self._r) # vertical acceleration in m / s^2

    def position(self):
        return (self._x, self._y)

    def turtlePosition(self):
        return (self._x/SCALING, self._y/SCALING)

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

    def pensize(self):
        return(self._pensize)

#---------------------------------------------------------------------------------
# Planets - global constants.

# The sun is at (0,0)
# The starting position for a planet is (Perihilion, 0)
# The starting velocity is (0, sMax), 
MERCURY = Planet(m_Me, P_Me, 0, 0, sMax_Me, TIME_STEP, "Mercury", "grey",      )
VENUS   = Planet(m_V,  P_V,  0, 0, sMax_V,  TIME_STEP, "Venus"  , "gold",      )
EARTH   = Planet(m_E,  P_E,  0, 0, sMax_E,  TIME_STEP, "Earth"  , "DeepSkyBlue")
MARS    = Planet(m_Ma, P_Ma, 0, 0, sMax_Ma, TIME_STEP, "Mars"   , "red",       )
JUPITER = Planet(m_J,  P_J,  0, 0, sMax_J,  TIME_STEP, "Jupiter", "yellow",    )
SATURN  = Planet(m_S,  P_S,  0, 0, sMax_S,  TIME_STEP, "Saturn" , "orange",    )
URANUS  = Planet(m_U,  P_U,  0, 0, sMax_U,  TIME_STEP, "Uranus" , "orange",    )
NEPTUN  = Planet(m_N,  P_N,  0, 0, sMax_N,  TIME_STEP, "Neptune", "blue",      )

S10 = sqrt(G*(M_S+m_E)/P_E)
# the speed of a made-up planet with the same mass and perihelion as Earth,
# but with a circular orbit.

# Made-up planets, for computational experiments
# with a default pencolor=green
PLANET07 = Planet(m_E, P_E, 0, 0, 0.7*S10, TIME_STEP, "Planet07") 
PLANET08 = Planet(m_E, P_E, 0, 0, 0.8*S10, TIME_STEP, "Planet08")
PLANET09 = Planet(m_E, P_E, 0, 0, 0.9*S10, TIME_STEP, "Planet09")
PLANET10 = Planet(m_E, P_E, 0, 0,     S10, TIME_STEP, "Planet10")
PLANET11 = Planet(m_E, P_E, 0, 0, 1.1*S10, TIME_STEP, "Planet11")
PLANET12 = Planet(m_E, P_E, 0, 0, 1.2*S10, TIME_STEP, "Planet12")
PLANET13 = Planet(m_E, P_E, 0, 0, 1.3*S10, TIME_STEP, "Planet13")

#=================================================================================
# AUXILIARY FUNCTIONS

def sky(skyColor="black", showSun=True):
    """Create black canvas with the white Sun at (0,0).
       The Sun is not to scale; it is shown much bigger."
    """
    # The Sun is always white. 
    screen = getscreen()
    screen.clear() # remove turtle image
    screen.screensize(5000,1000)
    screen.title("Kepler's world")
    screen.bgcolor(skyColor)
    t = Turtle(visible=False)
    if showSun:
        t.dot(10, "white") # SUN at (0,0)
        t.dot(1)
        print("The Sun is not to scale; it is shown much bigger.")
    del t

#sky()

#---------------------------------------------------------------------------------
   
def drawEllipse(semiMajorAxis, semiMinorAxis, leftShift=0, focusColor="white"):
    """Drwas an ellipse centered at (-leftShift,0) and the the foci;
       the foci are on the x-axis.
       For an ellipse centered at (0,0) use leftShift=0.
       For an ellipse with (0,0) in its right focus, use leftShift = c
       where c is the linear eccentricity i.e. center-to-focus distance.
       Note: make sure to create canvas before this function is called.
    """
    if semiMajorAxis < semiMinorAxis:
        print("You gave semi-major axis smaller than semi-minor axis.")
    t = Turtle(visible=False)
    t.speed("fastest")
    t.pendown()
    t.pensize(1)
    t.pencolor("white")
    c = sqrt(semiMajorAxis*semiMajorAxis - semiMinorAxis*semiMinorAxis)
    # c is the linear eccentricity, i.e. center-to-focus distance.
    t.teleport(c-leftShift,0) # left focus
    t.dot(6,"white") 
    t.teleport(-c-leftShift,0) # right focus
    t.dot(6,focusColor) 
    t.teleport(semiMajorAxis-leftShift,0) # left vertex of the ellipse
    for i in range(100+1):
        angle = 2*pi*(i/100)
        x = semiMajorAxis*cos(angle)
        y = semiMinorAxis*sin(angle)
        t.goto(x-leftShift,y)

#sky(showSun=False)
#drawEllipse(200,100)
        
#---------------------------------------------------------------------------------

def simulate0(planet: Planet):
    # version 0, just draw the orbit, do not verify laws, formulas.
    # This version is not used by the top level functions.
    # It is given here as a stepping stone to understand the version below.
    """Precondition: planet position (x,y) must have x>0, y=0,
                    and velocity (vx,vy) must have vx=0, vy>0.
       So, the planet must be in the right vertex of its elliptical orbit.
       This function calculats the orbit resulting from the continuing local
       effect of the Newton's law of gravity, assuming the Sun is at (0,0).
       Every step in the simulation is done for the actual planet
       (such as Mars, with its actual mass, perihelion and actual max speed)
       and only later it is scaled down to be displayed in turtle graphics.
       Note: make sure to create canvas before this function is called.
    """
    t = Turtle(visible=False)
    t.speed("fastest")
    t.pendown()
    t.pensize(planet.pensize())
    t.pencolor(planet.color())
    t.teleport(*planet.turtlePosition())
    rightVertex = planet.position() # the right vertex of the elliptical orbit,
                                    # assuming the the planet starts from the
                                    # x-axis with a vertical velocity.
    # Upper half of the orbit:
    done = False
    while True: 
        for j in range(100): # update turtle every 100 moves
            planet.move()
            if planet.position()[1]<0:
                done = True
                break
        t.goto(*planet.turtlePosition())
        if done: break

    # Lower half of the orbit:
    done = False
    while True: 
        for j in range(100): # update turtle every 100 moves
            planet.move()
            if planet.position()[1]>=0:
                done = True
                break
        t.goto(*planet.turtlePosition())
        if done: break

#sky()
#simulate0(PLANET13)

#-----------------------------------------------------------------------------

def simulate(planet: Planet):
    """Precondition: planet position (x,y) must have x>0, y=0,
                    and velocity (vx,vy) must have vx=0, vy>0.
       So, the planet must be in the right vertex of its elliptical orbit.
       This function calculats the orbit resulting from the continuing local
       effect of the Newton's law of gravity, assuming the Sun is at (0,0).
       Every step in the simulation is done for the actual planet
       (such as Mars, with its actual mass, perihelion and actual max speed)
       and only later it is scaled down to be displayed in turtle graphics.
       Returns the orbital period in seconds.
       Note: make sure to create canvas before this function is called.
    """
    t = Turtle(visible=False)
    t.speed("fastest")
    t.pendown()
    t.pensize(planet.pensize())
    t.pencolor(planet.color())
    t.teleport(*planet.turtlePosition())
    
    rightVertex = planet.position() # the right vertex of the elliptical orbit,
                                    # assuming the the planet starts from the
                                    # x-axis with a vertical velocity.
    maxX = rightVertex[0]
    minX = maxX           # will be updated 
    maxY = rightVertex[1] # will be updated
    
    x,y = planet.position()
    vx,vy = planet.velocity()

    area0 = x*vy-y*vx # The determinant of the matrix of column vectors r,v =
                        # area of the parallelogram spanned by vectors r, v =
                        # twice the area of a triangle spanned by vectors r, v.
    minArea = area0 # will be updated
    maxArea = area0 # will be updated
    # If the difference between minArea and maxArea is small,
    # Kepler's 3rd law will be confirmed.

    T = 0 # will be updated -- sidereal orbital period in seconds

    # Upper half of the orbit:
    done = False 
    while True: 
        for j in range(100): # update turtle every 100 moves
            planet.move()
            T += planet.timeStep() # update T
            x, y = planet.position() # update minX, maxY
            if x < minX: minX = x 
            if y > maxY: maxY = y
            vx, vy = planet.velocity() # update minArea, maxArea
            area = x*vy-y*vx
            if area < minArea: minArea = area 
            if area > maxArea: maxArea = area
            if planet.position()[1]<0: 
                done = True
                break
        t.goto(*planet.turtlePosition())
        if done: break

    # Lower half of the orbit:
    done = False 
    while True: 
        for j in range(100): # update turtle every 100 moves
            planet.move()
            T += planet.timeStep() # update T
            x,y = planet.position()
            vx,vy = planet.velocity() # update minArea, maxArea
            area = x*vy-y*vx
            if area < minArea: minArea = area 
            if area > maxArea: maxArea = area
            if planet.position()[1]>=0:
                done = True
                break
        t.goto(*planet.turtlePosition())
        if done: break
        
    a = (maxX - minX)/2 # semi-major axis
    b = maxY # semi-minor axis

    kepler3discrepancy = (maxArea - minArea) / area0
    kepler3discrepancyPercent = round(kepler3discrepancy*100, 2)

    print(kepler3discrepancyPercent, "%")
    
    return T #, a, b, kepler3discrepancy # orbital period, semi-major, semi-minor

#sky() # run/uncomment this before running simulate!    
#simulate(PLANET12)
#simulate(PLANET11)
#simulate(PLANET10)
#simulate(PLANET09)
#simulate(PLANET08)
#simulate(PLANET07)

#=================================================================================
# FUNCTIONS CALLED BY main

def innerPlanets():
    """A computer simulation of orbits of 4 inner planets, resulting from the
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
    print("The orbits displayed in turtle graphics are to scale.")
    print("and planets move counter-clockwise as seen from the north pole.")
    print("The orientation of the major axes of orbits is not modeled here:")
    print("all ellipses are shown with the major axis on the x-axis")
    print("and the sun in the right focus.")

    sky()
    
    print("\nTesting Kepler's 1st and 3rd laws")  

    print("\nMercury orbital period in days:")
    T = T_Me
    print(round(T/DAY,2), "- actual")
    drawEllipse(a_Me/SCALING, b_Me/SCALING, c_Me/SCALING, "gray")
    TS = simulate(MERCURY)
    print(round(abs(TS-T)*100/T, 2), "% error")

    print("\nVenus - orbital period in days:")
    T = T_V
    print(round(T/DAY,2), "- actual")
    drawEllipse(a_V/SCALING, b_V/SCALING, c_V/SCALING, "orange")
    TS = simulate(VENUS)
    print(round(abs(TS-T)*100/T, 2), "% error")

    print("\nEarth - orbital period in days:")
    T = T_E
    print(round(T/DAY,2), "- actual")
    drawEllipse(a_E/SCALING, b_E/SCALING, c_E/SCALING, "DeepSkyBlue")
    TS = simulate(EARTH)
    print(round(abs(TS-T)*100/T, 2), "% error")

    print("\nMars - orbital period in days:")
    T = T_Ma
    print(round(T/DAY,2), "- actual")
    drawEllipse(a_Ma/SCALING, b_Ma/SCALING, c_Ma/SCALING, "red")
    TS = simulate(MARS)
    print(round(abs(TS-T)*100/T, 2), "% error")

#innerPlanets() # uncomment this to simulate the inner planets.

#-------
#Output:

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

def testKepler(planet):
    """Precondition: planet position (x,y) must have x>0, y=0,
                    and velocity (vx,vy) must have vx=0, vy>0.
       So, the planet must be in the right vertex of its elliptical orbit.A computer simulation of the orbit of the planet, resulting from the
       continuing local effect of the Newton's law of gravity.
       Tests if the simulated planet obey (global) Kepler's laws 1 and 3.
       Every step in the simulation is done for the actual planet
       (such as Mars, with its actual mass, perihelion and actual max speed)
       and only later it is scaled down to be displayed in turtle graphics.
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
    P = planet.position()[0]  # Perihelion (shortest distance from the Sun)
    mu = G*(M_S+m) # gravitational parameter
    a = P*mu / (2*mu-P*vmax*vmax) # semi-major axis:
    c = a - P # the linear eccentricity, i.e. center-to-focus distance
    A = a + c # Aphelion (biggest distance from the Sun)
    b = sqrt(A*P)  # semi-minor axis
    T = sqrt(4*pi*pi*a*a*a / mu) # orbital period in s
    print("Planet's orbital period in days:")
    print(round(T/DAY,2), "- predicted by the theory")
    drawEllipse(a/SCALING, b/SCALING, c/SCALING, "green") # predicted orbit
    TS = simulate(planet) # simulation. TS - orbital period form simulation.
    print(round(TS/DAY,2), "- from the simulation")
    print( round(abs(TS-T)*100/T, 2), "% discrepancy")

testKepler(PLANET12)

#=================================================================================

def main():
    print("Welcome to Kepler's World.")

    while True:
        print(""" 
1. Inner planets: Mercury, Venus, Earth and Mars.
   For each planet, show the actual orbit, then simulate the movement.
   Print actual orbital period and that calculated in the simulation.
   
2. Test Kepler's laws: formulas vs. simulation.
   Use a made-up celestial body with the same mass and perihelion as Venus
   but with the maximal speed 1.2 times that of Venus.

3. Exit.
          """)
        choice = input("Enter your choice (1 or 2 or 3): ")
        if choice == "":
            return
        if choice[0] == "1":
            innerPlanets()
        elif choice[0] == "2":
            testKepler(PLANET12)
        else:
            print("Bye")
            return

if __name__ == "__main__":
    main()

#=================================================================================

