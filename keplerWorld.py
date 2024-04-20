"""
A simple simulation of movements of planets around the Sun.
It uses a discrete approximation of the real-world continuing
effect of Newton's law of gravity (F=GMm/rr) on the moving planet.

The calculations are done for actual planets,
and the results are later scaled down to be displayed in turtle graphics.

The simulation tests if the resulting orbits and planets' movements
obay the Kepler laws.

This provides an example of a local effects (of the law of gravity),
having a global effect (of elliptical orbits with details specified by Kepler.

The simulation shows planets' orbits to scale, 
but disregards some details not relevant to Kepler's laws.
Namely, in this simulation (unlike in reality):
* the orbits of all planets are in the same plane, and
* for each planet, the foci of its elliptical orbit are on the x-axis
  with the Sun in the right focus.
"""

#==================================================================================
from turtle import getscreen, Turtle
from math import sqrt, pi, sin, cos

#==================================================================================
# ASTRONOMICAL DATA FROM OBSERVATIONS AND CALCULATIONS 

# Astronomical data in kilograms, meters, seconds.
# From Wikipedia.

#---------------------------------------------------------------------------------

G = 6.67430e-11 # m^3 / (kg * s^2). Gravitational constant
DAY = 24*60*60  #s (actual day is 86400.002 s because of a slowing rotation)
KM = 1_000 # m
SEC = 1 # s 

#---------------------------------------------------------------------------------

def bcsFromAPm(A,P,m):
    """A - aphelion in m, P - perihelion in m, m - mass in kg."""
    a = (A + P) / 2 # semi-major axis in m.
    b = sqrt(A * P) # semi-minor axis in m.
    c = (A - P) / 2 # linear eccentricity = center-to-focus distance, in m.
    mu = G * (M_S + m) # gravitational parameter (with respect to the Sun)
    sMax = sqrt(mu * (2/P - 1/a)) # speed at perihelion - max speed
    sMin = sqrt(mu * (2/A - 1/a)) # speed at aphelion - min speed
    # e = c / a; # eccentricity.
    # T = sqrt(4 * pi * pi * a * a * a / mu) # orbital period.
    return b, c, sMax, sMin
    
#print(bcsFromAPm(A_Me, P_Me, m_Me))
#print(bcsFromAPm(A_V,  P_V,  m_V ))
#print(bcsFromAPm(A_E,  P_E,  m_E ))
#print(bcsFromAPm(A_Ma, P_Ma, m_Ma))

#---------------------------------------------------------------------------------
# Sun

# white
M_S = 1.98847e30 # Mass in kg.
R_S = 696_000*KM # Radius in m.

#---------------------------------------------------------------------------------
# Mercury

# grey
m_Me = 3.3011e23     # Mass in kg.
A_Me = 69_820_000*KM # Aphelion distance in m.
P_Me = 46_000_000*KM # Perihlion in m.
a_Me = 57_910_000*KM # Semi-major axis in m.
e_Me = 0.205630      # Eccentricity.
s_Me = 47.36*KM/SEC  # Average orbital speed in m/s.
T_Me = 87.9691*DAY   # Sidereal orbital period in s.

bcs = bcsFromAPm(A_Me, P_Me, m_Me)

b_Me    = bcs[0]     # Semi-minor axis in m.
c_Me    = bcs[1]     # Linear eccentricity = center-to-focus distance, in m.
sMax_Me = bcs[2]     # Max speed (at the prihelion) in m/s.
sMin_Me = bcs[3]     # Min speed (at the aphelion) in m/s.

#---------------------------------------------------------------------------------
# Venus

# yellow/white
m_V = 4.8675e24      # Mass in kg.
A_V = 108_940_000*KM # Aphelion distance in m.
P_V = 107_480_000*KM # Perihlion in m.
a_V = 108_210_000*KM # Semi-major axis in m.
e_V = 0.006772       # Eccentricity.
s_V = 35.02*KM/SEC   # Average orbital speed in m/s.
T_V = 224.701*DAY    # Sidereal orbital period in s.

bcs = bcsFromAPm(A_V, P_V, m_V)

b_V    = bcs[0]      # Semi-minor axis in m.
c_V    = bcs[1]      # Linear eccentricity = center-to-focus distance, in m.
sMax_V = bcs[2]      # Max speed (at the prihelion) in m/s.
sMin_V = bcs[3]      # Min speed (at the aphelion) in m/s.

#---------------------------------------------------------------------------------
# Earth

# blue
m_E = 5.972168e24        # Mass in kg.
A_E = 152_097_597*KM     # Aphelion distance in m.
P_E = 147_098_450*KM     # Perihlion in m.
a_E = 149_598_023*KM     # Semi-major axis in m.
e_E = 0.0167086          # Eccentricity
S_E = 29.7827*KM/SEC     # Average orbital speed in m/s.
T_E = 365.2563630040*DAY # Sidereal orbital period in s. 

bcs = bcsFromAPm(A_E, P_E, m_E)

b_E    = bcs[0]          # Semi-minor axis in m.
c_E    = bcs[1]          # Linear eccentricity = center-to-focus distance, in m.
sMax_E = bcs[2]          # Max speed (at the prihelion) in m/s.
sMin_E = bcs[3]          # Min speed (at the aphelion) in m/s.

#---------------------------------------------------------------------------------
# Mars

# red
m_Ma = 6.4171e23      # Mass in kg.
A_Ma = 249_261_000*KM # Aphelion distance in m.
P_Ma = 206_650_000*KM # Perihlion in m.
a_Ma = 227_939_366*KM # Semi-major axis in m.
e_Ma = 0.0934         # Eccentricity.
s_Ma = 24.07*KM/SEC   # Average orbital speed in m/s.
T_Ma = 686.980*DAY    # Sidereal orbital period in s.

bcs = bcsFromAPm(A_Ma, P_Ma, m_Ma)

b_Ma    = bcs[0]      # Semi-minor axis in m.
c_Ma    = bcs[1]      # Linear eccentricity = center-to-focus distance, in m.
sMax_Ma = bcs[2]      # Max speed (at the prihelion) in m/s.
sMin_Ma = bcs[3]      # Min speed (at the aphelion) in m/s.

#---------------------------------------------------------------------------------
# Jupiter

# requires a bigger scaling constant in simulations!

# yellow
m_J = 1.8982e27      # Mass in kg.
A_J = 816_363_000*KM # Aphelion distance in m.
P_J = 740_595_000*KM # Perihlion in m.
a_J = 778_479_000*KM # Semi-major axis in m.
e_J = 0.0489         # Eccentricity.
s_J = 13.07*KM/SEC   # Average orbital speed in m/s.
T_J = 4332.59*DAY    # Sidereal orbital period in s.

bcs = bcsFromAPm(A_J, P_J, m_J)

b_J    = bcs[0]    # Semi-minor axis in m.
c_J    = bcs[1]    # Linear eccentricity = center-to-focus distance, in m.
sMax_J = bcs[2]    # Max speed (at the prihelion) in m/s.
sMin_J = bcs[3]    # Min speed (at the aphelion) in m/s.

#---------------------------------------------------------------------------------
# Saturn

# requires a bigger scaling constant in simulations!

# yellow
m_S = 5.6834e26        # Mass in kg.
A_S = 1_514_500_000*KM # Aphelion distance in m.
P_S = 1_352_550_000*KM # Perihlion in m.
a_S = 1_433_530_000*KM # Semi-major axis in m.
e_S = 0.0565           # Eccentricity.
s_S = 9.68*KM/SEC      # Average orbital speed in m/s.
T_S = 10755.70*DAY     # Sidereal orbital period in s. 


bcs = bcsFromAPm(A_S, P_S, m_S)

b_S    = bcs[0]    # Semi-minor axis in m.
c_S    = bcs[1]    # Linear eccentricity = center-to-focus distance, in m.
sMax_S = bcs[2]    # Max speed (at the prihelion) in m/s.
sMin_S = bcs[3]    # Min speed (at the aphelion) in m/s.

#---------------------------------------------------------------------------------
# Uranus

# requires a bigger scaling constant in simulations!

# white
m_U = 8.6810e25        # Mass in kg.
A_U = 3_006_390_000*KM # Aphelion distance in m.
P_U = 2_735_560_000*KM # Perihlion in m.
a_U = 2_870_972_000*KM # Semi-major axis in m.
e_U = 0.04717          # Eccentricity.
s_U = 6.80*KM/SEC      # Average orbital speed in m/s.
T_U = 30,688.5*DAY     # Sidereal orbital period in s.

bcs = bcsFromAPm(A_U, P_U, m_U)

b_U    = bcs[0]        # Semi-minor axis in m.
c_U    = bcs[1]        # Linear eccentricity = center-to-focus distance, in m.
sMax_U = bcs[2]        # Max speed (at the prihelion) in m/s.
sMin_U = bcs[3]        # Min speed (at the aphelion) in m/s.

#---------------------------------------------------------------------------------
# Neptune

# requires a bigger scaling constant in simulations!

# blue
m_N = 1.02409e26       # Mass in kg.
A_N = 4_540_000_000*KM # Aphelion distance in m.
P_N = 4_460_000_000*KM # Perihlion in m.
a_N = 4_500_000_000*KM # Semi-major axis in m.
e_N = 0.008678         # Eccentricity.
s_N = 5.43*KM/SEC      # Average orbital speed in m/s.
T_N = 60195*DAY        # Sidereal orbital period in s.

bcs = bcsFromAPm(A_N, P_N, m_N)

b_N    = bcs[0]    # Semi-minor axis in m.
c_N    = bcs[1]    # Linear eccentricity = center-to-focus distance, in m.
sMax_N = bcs[2]    # Max speed (at the prihelion) in m/s.
sMin_N = bcs[3]    # Min speed (at the aphelion) in m/s.

#=================================================================================
# SIMULATION PARAMETERS

# A scaling factor suitable for showing orbits of the 4 inner planets
# and made-up planets 0.7 - 1.3
# Outer planets require a factor 10 times bigger.
SCALING = 1000000_000 # Real distances in meters will be divided by this
                     # before being given to the turtle.

# Position, velocity, acceleration will be updated every TIME_STEP.
TIME_STEP = 1_000 # seconds
# Suitable for inner plenets. Outer planets can use a step 10 times bigger.
# the bigger the time step, the faster the turtle moves.
        
#=================================================================================
# PLANETS

class Planet(object):

    def __init__(self, name, mass, perihelionDistance, maxSpeed, timeStep,
                 color="green"):
        """mass in kg,
           perihelionDistance in m.
           maxSpeed in m/s,
           simulation timeStep in s.
           Assumes that the Sun is at (0,0).
           Puts the planet/self at (perihelionDistance, 0)
           giving it velocity (0, maxSpeed).
        """
        self._mass = mass
        self._x = perihelionDistance
        self._y = 0
        self._r2 = self._x*self._x + self._y*self._y # radius squared
                   # we store radius squared 
        self._r = sqrt(self._r2) # radius = distance from Sun's center.
        self._vx = 0 # horizontal speed in m/SEC 
        self._vy = maxSpeed # vertical speed in m/s
        self._ax = -G*M_S*self._x/(self._r2*self._r) # horizontal acceleration m/s^2
        self._ay = -G*M_S*self._y/(self._r2*self._r) # vertical acceleration m/s^2
        self._timeStep = timeStep
        self._name = name
        self._color = color

    def move(self, timeStep=None):
        """Updates position, velocity, acceleration, radius
           to those after timeStep.
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
        self._ax = -G*M_S*self._x/(self._r2*self._r) # horizontal acceleration m/s^2
        self._ay = -G*M_S*self._y/(self._r2*self._r) # vertical acceleration m/s^2

    def position(self, scaleDownFactor=1):
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

#---------------------------------------------------------------------------------
# Planets - global constants.

# The Sun is at (0,0)
# The starting position for a planet is (perihilionDistance, 0)
# The starting velocity is (0, sMax),

# Inner planets - rocky:
MERCURY = Planet("Mercury", m_Me, P_Me, sMax_Me, TIME_STEP, "grey"       )
VENUS   = Planet("Venus"  , m_V,  P_V,  sMax_V,  TIME_STEP, "gold"       )
EARTH   = Planet("Earth"  , m_E,  P_E,  sMax_E,  TIME_STEP, "DeepSkyBlue")
MARS    = Planet("Mars"   , m_Ma, P_Ma, sMax_Ma, TIME_STEP, "red"        )
# Outer planets - gas giants:
JUPITER = Planet("Jupiter", m_J,  P_J,  sMax_J,  TIME_STEP, "yellow"     )
SATURN  = Planet("Saturn" , m_S,  P_S,  sMax_S,  TIME_STEP, "orange"     )
URANUS  = Planet("Uranus" , m_U,  P_U,  sMax_U,  TIME_STEP, "orange"     )
NEPTUNE = Planet("Neptune", m_N,  P_N,  sMax_N,  TIME_STEP, "blue"       )

S10 = sqrt(G*(M_S+m_E)/P_E)
# the speed of a made-up planet with the same mass and perihelion as Earth,
# but with a circular orbit.

# Made-up planets, for computational experiments (default pencolor=green)
PLANET07 = Planet("Planet 0.7", m_E, P_E, 0.7*S10, TIME_STEP) 
PLANET08 = Planet("Planet 0.8", m_E, P_E, 0.8*S10, TIME_STEP)
PLANET09 = Planet("Planet 0.9", m_E, P_E, 0.9*S10, TIME_STEP)
PLANET10 = Planet("Planet 1.0", m_E, P_E,     S10, TIME_STEP)
PLANET11 = Planet("Planet 1.1", m_E, P_E, 1.1*S10, TIME_STEP)
PLANET12 = Planet("Planet 1.2", m_E, P_E, 1.2*S10, TIME_STEP)
PLANET13 = Planet("Planet 1.3", m_E, P_E, 1.3*S10, TIME_STEP)

# Note.
# Inner planets and made up planets above are close to the Sun
# and outer ones are at vast distances.
# It is not practical to show them all on the same canvas.
# Either show the inner planets together with made up planets
# or the outer planets.

#=================================================================================
# AUXILIARY FUNCTIONS

def sky(skyColor="black", showSun=True):
    """Create black canvas with the white Sun at (0,0).
       The Sun size is not to scale; it is shown much bigger."
    """
    screen = getscreen()
    screen.clear() # remove turtle image
    screen.screensize(5000,1000)
    screen.title("Kepler's world")
    screen.bgcolor(skyColor)
    t = Turtle(visible=False)
    if showSun:
        t.dot(10, "white") # SUN at (0,0) and is white.
        print("The Sun is not to scale; it is shown much bigger.")
    del t

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
       Note: the parameters are in turtle canvas units, (not in meters).
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
    for i in range(100+1):
        angle = 2*pi*(i/100)
        x = semiMajorAxis*cos(angle)
        y = semiMinorAxis*sin(angle)
        t.goto(x-leftShift,y)

#sky(showSun=False) # run/uncomment this before running drawEllipse
#drawEllipse(200,100, 0, "pink")
     
#---------------------------------------------------------------------------------

def simulate(planet: Planet):
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
       Note: make sure to create canvas before this function is called.
    """
    # Prepare turtle.
    t = Turtle(visible=False)
    t.speed("fastest")
    t.pendown()
    t.pensize(3)
    t.pencolor(planet.color())
    t.teleport(*planet.position(scaleDownFactor=SCALING))

    # The planet starts from its perihelion.
    
    # Upper half of the orbit:
    done = False
    while True: 
        for j in range(100): # update turtle every 100 moves
            planet.move()
            if planet.position()[1]<0: # if y<0
                done = True
                break
        t.goto(*planet.position(scaleDownFactor=SCALING))
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
        t.goto(*planet.position(scaleDownFactor=SCALING))
        if done: break

    # The planet is back at the perihelion.

#sky() # run/uncomment this before running simulate!
#simulate(PLANET12)

#-----------------------------------------------------------------------------

def simulateAndTest(planet: Planet):
    """Precondition: planet position (x,y) must have x>0, y=0,
                     and velocity (vx,vy) must have vx=0, vy>0.
       So, the planet must be in the right vertex of its elliptical orbit.
       This function calculates the orbit resulting from the continuing local
       effect of the Newton's law of gravity, assuming the Sun is at (0,0).
       Every step in the simulation is done for the actual planet
       (such as Mars, with its actual mass, perihelion and actual max speed)
       and only later it is scaled down to be displayed in turtle graphics.
       Tests Kepler's laws ...
       Returns the orbital period in seconds, ...
       Note: make sure to create canvas before this function is called.
    """
    # Prepare turtle
    t = Turtle(visible=False)
    t.speed("fastest")
    t.pendown()
    t.pensize(3)
    t.pencolor(planet.color())
    t.teleport(*planet.position(scaleDownFactor=SCALING))

    # Concerning Kepler's 1st law
    x,y = planet.position() # perihelion,
                            # the right vertex of the elliptical orbit 
    maxX = x # perihelion distance from the Sun (Sun is at (0,0))
    minXsoFar = maxX            
    maxYsoFar = y
    # These will be used to find the semi-major axis and semi-minor axis.

    # Concerning Kepler's 2nd law
    vx,vy = planet.velocity()
    area0 = x*vy-y*vx # The determinant of the matrix of column vectors r,v =
                      # = area of the parallelogram spanned by vectors r,v =
                      # = twice the area of a triangle spanned by vectors r,v.
                      # = vector cross product  r x v.
                      # Notice that angular momentum is  r x mv.
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
        t.goto(*planet.position(scaleDownFactor=SCALING))
        if done: break

    # The planet is now at its aphelion.

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
        t.goto(*planet.position(scaleDownFactor=SCALING))
        if done: break

    # The planet is back at the perihelion.

    # Concerning Kepler's 1st law:
    minX = minXsoFar 
    maxY = maxYsoFar
    a = (maxX - minX)/2 # semi-major axis - from the simulation
    b = maxY # semi-minor axis - from the simulation
    # Notice that maxX is the distance from (0,0) to the perihelion.
    orbitCenter = maxX - a
    c = sqrt(a*a - b*b) # linear eccentricity = center to focus distance.
    t.teleport((orbitCenter + c)/SCALING,0) # draw right focus
    t.dot(6, planet.color())
    t.teleport((orbitCenter - c)/SCALING,0) # draw left focus
    t.dot(6, planet.color()) 

    # Concerning Kepler's 2nd law: 
    minArea = minAreaSoFar
    maxArea = maxAreaSoFar
    kepler2discrepancy = (maxAreaSoFar - minArea) / area0
    kepler2discrepancyPercent = round(kepler2discrepancy*100, 2)
    print(kepler2discrepancyPercent, "%")    

    # Concerning Kepler's 3rd law:
    T = TsoFar # sidereal orbital period - from the simulation.
    # T^2 / a^2 = 4*pi^2 / G(M+m) - does the simulation support this?
    #print("Kepler3:")
    lhs = T*T/(a*a*a)
    #print(lhs)
    rhs = 4*pi*pi / (G*(M_S+planet._mass))
    #print(rhs)    
    #print(abs((lhs-rhs)/rhs))
    #print(round(abs((lhs-rhs)/rhs)*100,2), "%")
    
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
    drawEllipse(a_Me/SCALING, b_Me/SCALING, c_Me/SCALING)
    TS = simulateAndTest(MERCURY)[3]
    print(round(abs(TS-T)*100/T, 2), "% error")

    print("\nVenus - orbital period in days:")
    T = T_V
    print(round(T/DAY,2), "- actual")
    drawEllipse(a_V/SCALING, b_V/SCALING, c_V/SCALING)
    TS = simulateAndTest(VENUS)[3]
    print(round(abs(TS-T)*100/T, 2), "% error")

    print("\nEarth - orbital period in days:")
    T = T_E
    print(round(T/DAY,2), "- actual")
    drawEllipse(a_E/SCALING, b_E/SCALING, c_E/SCALING)
    TS = simulateAndTest(EARTH)[3]
    print(round(abs(TS-T)*100/T, 2), "% error")

    print("\nMars - orbital period in days:")
    T = T_Ma
    print(round(T/DAY,2), "- actual")
    drawEllipse(a_Ma/SCALING, b_Ma/SCALING, c_Ma/SCALING)
    TS = simulateAndTest(MARS)[3]
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
    P = planet.position()[0]  # Perihelion distance (shortest distance from Sun)
    mu = G*(M_S+m) # gravitational parameter
    a = P*mu / (2*mu-P*vmax*vmax) # semi-major axis:
    c = a - P # the linear eccentricity, i.e. center-to-focus distance
    A = a + c # Aphelion distance (biggest distance from the Sun)
    b = sqrt(A*P)  # semi-minor axis
    T = sqrt(4*pi*pi*a*a*a / mu) # orbital period in s
    print("Planet's orbital period in days:")
    print(round(T/DAY,2), "- predicted by the theory")
    drawEllipse(a/SCALING, b/SCALING, c/SCALING) # predicted orbit
    TS = simulateAndTest(planet)[3] # simulation. TS - orbital period form simulation.
    print(round(TS/DAY,2), "- from the simulation")
    print( round(abs(TS-T)*100/T, 2), "% discrepancy")

#testKepler(PLANET12)

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

if __name__ == "__main__": main()

#=================================================================================


#innerPlanets()
#simulateAndTest(JUPITER)
#simulateAndTest(SATURN)
#simulateAndTest(URANUS)
#simulateAndTest(NEPTUNE)
