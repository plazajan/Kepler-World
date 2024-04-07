
import turtle
import math

#----------------------------------------------------------------------------------

# Astronomical data in kg,m,sec, from Wikipedia.

#----------------------------------------------------------------------------------

G = 6.67430e-11 # m^3 / (kg * s^2). Gravitational constant
DAY = 24*60*60  # s

#----------------------------------------------------------------------------------
# Sun

# white
M_S = 1.98847e30 # Mass in kg.
# R_S = 696000000  # Radius in m.

#----------------------------------------------------------------------------------
# Mercury

# grey
A_Me = 69820000000 # Aphelion in m.
P_Me = 46000000000 # Perihlion in m.
a_Me = 57910000000 # Semi-major axis in m.
e_Me = 0.205630    # Eccentricity.
T_Me = 7600521.6   # Sidereal orbital period in s. (87.9691 days)
s_Me = 47360       # Average orbital speed in m/s.
m_Me = 3.3011e23   # Mass in kg.

# Calculated from those above:
b_Me = 56672038961 # Semi-minor axis in m.
c_Me = 11910000000 # Linear eccentricity = center-to-focus distance, in m.
sMax_Me = 58979    # Max speed (at the prihelion) in m/s.

#----------------------------------------------------------------------------------
# Venus

# yellow/white
A_V = 108940000000 # Aphelion in m.
P_V = 107480000000 # Perihlion in m.
a_V = 108210000000 # Semi-major axis in m.
e_V = 0.006772     # Eccentricity.
T_V = 19414166.4   # Sidereal orbital period in s. (224.701 days)
s_V = 35020        # Average orbital speed in m/s.
m_V = 4.8675e24    # Mass in kg.

# Calculated from those above:
b_V = 108207537630 # Semi-minor axis in m.
c_V = 730000000    # Linear eccentricity = center-to-focus distance, in m.
sMax_V = 35258     # Max speed (at the prihelion) in m/s.

#----------------------------------------------------------------------------------
# Earth

# blue
A_E = 152097597000     # Aphelion in m.
P_E = 147098450000     # Perihlion in m.
a_E = 149598023000     # Semi-major axis in m.
e_E = 0.0167086        # Eccentricity
T_E = 31558149.7635456 # Sidereal orbital period in s. (365.2563630040 days)
S_E = 29782.7          # Average orbital speed in m/s.
m_E = 5.972168e24      # Mass in kg.

# Calculated from those above:
b_E = 149577139855 # Semi-minor axis in m.
c_E = 2499573500   # Linear eccentricity = center-to-focus distance, in m.
sMax_E = 30287     # Max speed (at the prihelion) in m/s.

#----------------------------------------------------------------------------------
# Mars

# red
A_Ma = 249261000000 # Aphelion in m.
P_Ma = 206650000000 # Perihlion in m.
a_Ma = 227939366000 # Semi-major axis in m.
e_Ma = 0.0934       # Eccentricity.
T_Ma = 59355072     # Sidereal orbital period in s. (686.980 dyas)
s_Ma = 24070        # Average orbital speed in m/s.
m_Ma = 6.4171e23    # Mass in kg.

# Calculated from those above:
b_Ma = 226957673697 # Semi-minor axis in m.
c_Ma = 21305500000  # Linear eccentricity = center-to-focus distance, in m.
sMax_Ma = 26500     # Max speed (at the prihelion) in m/s.

#----------------------------------------------------------------------------------

def calculateFromAPm(A,P,m):
    """A - aphelion in m, P - perihelion in m, m - mass in kg."""
    a= (A+P)/2 # semi-major axis in m.
    b=math.sqrt(A*P) # semi-minor axis in m.
    c=(A-P)/2 # Linear eccentricity = center-to-focus distance, in m.
    e=c/a; # eccentricity.
    mu=G*(M_S+m) # Gravitational parameter.
    sMax=math.sqrt(mu*(2/P-1/a)) # speed at perihelion (max speed)
    T=math.sqrt(4*a*a*a*math.pi*math.pi/mu) # Orbital period.
#    print("a", a)
#    print("e", e)
#    print("T", T)
#    print()
    print("b", b)
    print("c", c)
    print("sMax", sMax)
    
#calculateFromAPm(A_Me, P_Me, m_Me)
#calculateFromAPm(A_V,  P_V,  m_V )
#calculateFromAPm(A_E,  P_E,  m_E )
#calculateFromAPm(A_Ma, P_Ma, m_Ma)

#----------------------------------------------------------------------------------
# Simulation parameters

SCALING = 1000000000 # Real distances in meters will be divided by this
                     # before being given to the turtle.

TIME_STEP = 1000 # seconds.
        # Position, velocity, acceleration will be updated every TIME_STEP.

#----------------------------------------------------------------------------------

class Planet(object):

    def __init__(self, mass, x, y, vx, vy, timeStep,
                 name=None, color="white", pensize=1):
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
        self._r = math.sqrt(self._r2) # radius = distance between sun and self.
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
        self._r = math.sqrt(self._r2)
        self._vx -= self._ax * timeStep
        self._vy -= self._ay * timeStep
        self._ax = G*M_S*self._x/(self._r2*self._r) # horizontal acceleration in m / s^2
        self._ay = G*M_S*self._y/(self._r2*self._r) # vertical acceleration in m / s^2

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

    def timeStep(self):
        return self._timeStep

    def name(self):
        return(self._name)

    def color(self):
        return(self._color)

    def pensize(self):
        return(self._pensize)



#----------------------------------------------------------------------------------
# Sample planets - global constants.

# The sun is at (0,0)
# The starting position for a planet is (Perihilion, 0)
# The starting velocity is (0, sMax), 
MERCURY = Planet(m_Me, P_Me, 0, 0, sMax_Me, TIME_STEP, "Mercury", "grey",        3)
VENUS   = Planet(m_V,  P_V,  0, 0, sMax_V,  TIME_STEP, "Venus"  , "gold",        3)
EARTH   = Planet(m_E,  P_E,  0, 0, sMax_E,  TIME_STEP, "Earth"  , "DeepSkyBlue", 3)
MARS    = Planet(m_Ma, P_Ma, 0, 0, sMax_Ma, TIME_STEP, "Mars"   , "red",         3)

PLANET07 = Planet(m_V, P_V, 0, 0, 0.7*sMax_V, TIME_STEP, "Planet07", "green", 3) 
PLANET08 = Planet(m_V, P_V, 0, 0, 0.8*sMax_V, TIME_STEP, "Planet08", "green", 3)
PLANET09 = Planet(m_V, P_V, 0, 0, 0.9*sMax_V, TIME_STEP, "Planet09", "green", 3)
PLANET10 = VENUS
PLANET11 = Planet(m_V, P_V, 0, 0, 1.1*sMax_V, TIME_STEP, "Planet11", "green", 3)
PLANET12 = Planet(m_V, P_V, 0, 0, 1.2*sMax_V, TIME_STEP, "Planet12", "green", 3)
PLANET13 = Planet(m_V, P_V, 0, 0, 1.3*sMax_V, TIME_STEP, "Planet13", "green", 3)

#----------------------------------------------------------------------------------

def sky(skyColor="black", showSun=True):
    """Create black canvas with the white Sun at (0,0).
       The Sun is not to scale; it is shown much bigger."
    """
    # The Sun is always white. 
    screen = turtle.getscreen()
    screen.clear() # remove turtle image
    screen.screensize(5000,1000)
    screen.title("Kepler's world")
    screen.bgcolor(skyColor)
    t = turtle.Turtle(visible=False)
    if showSun:
        t.dot(10, "white") # SUN at (0,0)
        t.dot(1)
        print("The Sun is not to scale; shown much bigger.")
    del t

#sky()

#----------------------------------------------------------------------------------
   
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
    t = turtle.Turtle(visible=False)
    t.speed("fastest")
    t.pendown()
    t.pensize(1)
    t.pencolor("white")
    c = math.sqrt(semiMajorAxis*semiMajorAxis - semiMinorAxis*semiMinorAxis)
    # c is the linear eccentricity, i.e. center-to-focus distance.
    t.teleport(c-leftShift,0) # left focus
    t.dot(6,"white") 
    t.teleport(-c-leftShift,0) # right focus
    t.dot(6,focusColor) 
    t.teleport(semiMajorAxis-leftShift,0) # left verttex of the ellipse
    for i in range(100+1):
        angle = 2*math.pi*(i/100)
        x = semiMajorAxis*math.cos(angle)
        y = semiMinorAxis*math.sin(angle)
        t.goto(x-leftShift,y)

#sky(showSun=False)
#drawEllipse(200,100)
        
#----------------------------------------------------------------------------------

def drawOrbit(planet: Planet):
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
    t = turtle.Turtle(visible=False)
    t.speed("fastest")
    t.pendown()
    t.pensize(planet.pensize())
    t.pencolor(planet.color())
    t.teleport(*planet.turtlePosition())
    rightVertex = planet.position() # the right vertex of the elliptical orbit,
                                    # assuming the the planet starts from the
                                    # x-axis with a vertical velocity.
    maxX = minX = rightVertex[0]
    maxY = rightVertex[1]
    T = 0 # sidereal orbital period in seconds
    while planet.position()[1]>=0:
        for j in range(10):
            planet.move()
            T += planet.timeStep()
            currentX, currentY = planet.position()
            if currentX < minX: minX = currentX 
            if currentY > maxY: maxY = currentY
        t.goto(*planet.turtlePosition())
    while t.pos()[1] <= 0:
        for j in range(10):
            planet.move()
            T += planet.timeStep()
        t.goto(*planet.turtlePosition())
    #print(round(T/DAY,2), "- from the simulation, in days.")
    return T # orbital period in seconds

#sky() # run/uncomment this before running drawOrbit! 
    
#drawOrbit(MERCURY)
#drawOrbit(PLANET12)
#drawOrbit(PLANET11)
#drawOrbit(PLANET10)
#drawOrbit(PLANET09)
#drawOrbit(PLANET08)
#drawOrbit(PLANET07)

#----------------------------------------------------------------------------------

def innerPlanets():
    """A computer simulation of orbits of 4 inner planets, resulting from the
       continuing local effect of the Newton's law of gravity.
       Tests if the simulated planets obey (global) Kepler's laws 1 and 3.
       Every step in the simulation is done for the actual planet
       (such as Mars, with its actual mass, perihelion and actual max speed)
       and only later it is scaled down to be displayed in turtle graphics.
       This function creates its own canvas-sky with the Sun,
       and draws an ellipse of the actual orbit with the foci,
       before starting simulation.
       It produces output, besides displaying turtle graphics.
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
    TS = drawOrbit(MERCURY)
    print(round(abs(TS-T)*100/T, 2), "% error")

    print("\nVenus - orbital period in days:")
    T = T_V
    print(round(T/DAY,2), "- actual")
    drawEllipse(a_V/SCALING, b_V/SCALING, c_V/SCALING, "orange")
    TS = drawOrbit(VENUS)
    print(round(abs(TS-T)*100/T, 2), "% error")

    print("\nEarth - orbital period in days:")
    T = T_E
    print(round(T/DAY,2), "- actual")
    drawEllipse(a_E/SCALING, b_E/SCALING, c_E/SCALING, "DeepSkyBlue")
    TS = drawOrbit(EARTH)
    print(round(abs(TS-T)*100/T, 2), "% error")

    print("\nMars - orbital period in days:")
    T = T_Ma
    print(round(T/DAY,2), "- actual")
    drawEllipse(a_Ma/SCALING, b_Ma/SCALING, c_Ma/SCALING, "red")
    TS = drawOrbit(MARS)
    print(round(abs(TS-T)*100/T, 2), "% error")

#innerPlanets() # uncomment this to simulate the inner planets.

#-------
#Output:

# The orbits are to scale.
# and planets move counter-clockwise as seen from the north pole.
# The orientation of the major axes of orbits is not modeled here:
# all ellipses are shown with the major axis on the x-axis
# and the sun in the right focus.
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
   
#----------------------------------------------------------------------------------

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
    b = math.sqrt(A*P)  # semi-minor axis
    T = math.sqrt(4*math.pi*math.pi*a*a*a / mu) # orbital period in s
    print("Planet's orbital period in days:")
    print(round(T/DAY,2), "- predicted by the theory")
    drawEllipse(a/SCALING, b/SCALING, c/SCALING, "green") # predicted orbit
    TS = drawOrbit(planet) # simulation. TS - orbital period form simulation.
    print(round(TS/DAY,2), "- from the simulation")
    print( round(abs(TS-T)*100/T, 2), "% discrepancy")

testKepler(PLANET12)

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

