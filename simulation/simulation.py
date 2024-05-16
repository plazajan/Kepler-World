
#==================================================================================

from planets import *

#=================================================================================
# AUXILIARY FUNCTIONS

def sky(skyColor="black", showSun=True):
    """Create black canvas with the white Sun at (0,0).
       The Sun size is not to scale; it is shown much bigger."
    """
    screen = getscreen()
    screen.clear() # remove turtle image
    screen.screensize(2000,900)
    screen.title("Kepler's world")
    screen.bgcolor(skyColor)
    if showSun:
        turtle = Turtle(visible=False)
        turtle.dot(10, "white") # SUN is at (0,0) and is white.
        print("The Sun is not to scale; it is shown much bigger.")

#sky()

#---------------------------------------------------------------------------------
   
def drawEllipse(semiMajorAxis, semiMinorAxis, leftShift=0,
                color="pink"):
    """Draws an ellipse centered at (-leftShift,0) and shows the foci;
       the foci are on the x-axis.
       For an ellipse centered at (0,0) use leftShift=0.
       For an ellipse with the right focus at (0,0), use leftShift = c,
       where c is the linear eccentricity i.e. center-to-focus distance:
       c = sqrt(semiMajorAxis*semiMajorAxis - semiMinorAxis*semiMinorAxis)
       Note: make sure to create canvas before this function is called.
       Note: the parameteRADIUS_SUN are in turtle canvas units, not in meteRADIUS_SUN.
    """
    if semiMajorAxis < semiMinorAxis:
        raise ValueError(
            "semi-major axis cannot be smaller than semi-minor axis.")
    turtle = Turtle(visible=False)
    turtle.speed("fastest")
    turtle.pendown()
    turtle.pensize(1)
    turtle.pencolor(color)
    c = sqrt(semiMajorAxis*semiMajorAxis - semiMinorAxis*semiMinorAxis)
    # c is the linear eccentricity, i.e. center-to-focus distance.
    turtle.teleport(c-leftShift,0) # left focus
    turtle.dot(6, color) 
    turtle.teleport(-c-leftShift,0) # right focus
    turtle.dot(6, color) 
    turtle.teleport(semiMajorAxis-leftShift,0) # left vertex of the ellipse
    n = 500
    for i in range(n+1):
        angle = 2*pi*(i/n)
        x = semiMajorAxis*cos(angle)
        y = semiMinorAxis*sin(angle)
        turtle.goto(x-leftShift,y)

#sky(showSun=False) # run/uncomment this before running drawEllipse
#drawEllipse(200, 100, 0, "pink")


#---------------------------------------------------------------------------------

def predictByTheory(planet: SimulatedPlanet, displayParameters: dict):

    orbitScaleDownFactor = displayParameters["orbitScaleDownFactor"]
    graphStartX = displayParameters["graphStartX"]
    sweepSpeedScaleDownFactor = displayParameters["sweepSpeedScaleDownFactor"]
    timeScaleDownFactor = displayParameters["timeScaleDownFactor"]
    
##    t = planetData["orbital period"]
##    a = planetData["semi-major axis"]
##    b = planetData["semi-minor axis"]
##    c = planetData["linear eccentricity"]
##    vmax = planetData["max speed"]
##    perihelion = planetData["perihelion"]
    
    mu = G*(MASS_SUN+planet.mass()) # gravitational parameter
    perihelion = planet.perihelion()
    vmax = planet.maxSpeed()
    a = perihelion*mu / (2*mu-perihelion*vmax*vmax) # semi-major axis:
    c = a - perihelion # the linear eccentricity, i.e. center-to-focus distance
    aphelion = a + c # predicted aphelion distance
    b = sqrt(aphelion*perihelion)  # predicted semi-minor axis
    t = 2*pi*sqrt(a*a*a/mu) # predicted orbital period in s
    sweepSpeed0 = perihelion * vmax / 2

    drawEllipse(a/orbitScaleDownFactor, b/orbitScaleDownFactor,
                c/orbitScaleDownFactor)

    turtle2 = Turtle(visible=False)
    turtle2.speed("fastest")
    turtle2.pendown()
    turtle2.pensize(1)
    turtle2.pencolor("pink")

    turtle2.teleport(graphStartX,
                     sweepSpeed0 / sweepSpeedScaleDownFactor)
    turtle2.speed(3)

    # draw a pink hirizontal line with a blunt end 
    turtle2.goto(graphStartX+t/timeScaleDownFactor,
                     sweepSpeed0 /sweepSpeedScaleDownFactor)
    turtle2.left(90)
    turtle2.forward(3)
    turtle2.back(6)

#sky(showSun=False) # run/uncomment this before running predictByTheory

#---------------------------------------------------------------------------------

def simulate(planet: SimulatedPlanet):
    # Just draws the orbit, does not test Kepler's laws.
    # This function is not used by the top level functions in the program.
    # It is given here as a stepping stone to undeRADIUS_SUNtand simulateAndTest below.
    """Precondition: planet position (x,y) must have x>0, y=0,
                     and velocity (vx,vy) must have vx=0, vy>0.
       So, the planet must be in the right vertex of its elliptical orbit.
       This function calculates the orbit resulting from the continuing local
       effect of the Newton's law of gravity, assuming the Sun is at (0,0).
       Every step in the simulation is done for the actual planet
       (such as Mars, with its actual mass, perihelion and actual max speed)
       and only later it is scaled down to be displayed in turtle graphics.
       Real distances in meteRADIUS_SUN will be divided by orbitScaleDownFactor
       before being given to the turtle.
       Note: make sure to create canvas before this function is called.
    """
    # Prepare turtle.
    turtle = Turtle(visible=False)
    turtle.speed("fastest")
    turtle.pendown()
    turtle.pensize(1)
    turtle.pencolor(planet.color())

    orbitScaleDownFactor =  1_000_000_000

    (x,y) = planet.position() # the right vertex of the elliptical orbit
    turtle.teleport(x/orbitScaleDownFactor, y/orbitScaleDownFactor)

    # The planet starts from its perihelion.
    
    # Upper half of the orbit:
    done = False
    while True: 
        for j in range(100): # update turtle every 100 moves
            planet.move()
            (x,y) = planet.position()
            if y <= 0:
                done = True
                break
        turtle.goto(x/orbitScaleDownFactor, y/orbitScaleDownFactor)
        if done: break

    # The planet is now at its aphelion (or a little past)

    # Lower half of the orbit:
    done = False
    while True: 
        for j in range(100): # update turtle every 100 moves
            planet.move()
            (x,y) = planet.position()
            if y >= 0:
                done = True
                break
        turtle.goto(x/orbitScaleDownFactor, y/orbitScaleDownFactor)
        if done: break

    # The planet is back at the perihelion (or a little past)

#sky() # run/uncomment this before running simulate!
#simulate(PLANET1_2)

#-----------------------------------------------------------------------------

# Under construction
def simulateAndTest(planet: SimulatedPlanet, displayParameters):
    # The sweepSpeed is half the area spanned by r and v;
    # it will be graphed for every t from 0 to T
    # where T is the orbital period.
    # If the line in the graph is horizontal, it will confirm Kepler's 2nd law.
    
    """Precondition: planet position (x,y) must have x>0, y=0,
                     and velocity (vx,vy) must have vx=0, vy>0.
       So, the planet must be in the right vertex of its elliptical orbit.
       This function calculates the orbit resulting from the continuing local
       effect of the Newton's law of gravity, assuming the Sun is at (0,0).
       Every step in the simulation is done for the actual planet
       (such as Mars, with its actual mass, perihelion and actual max speed)
       and only later it is scaled down to be displayed in turtle graphics.
       Real distances in meteRADIUS_SUN will be divided by orbitScaleDownFactor
       before being given to the turtle.
       Tests Kepler's laws ...
       Returns the orbital period in seconds, ...
       Note: make sure to create canvas before this function is called.
    """
    # Prepare turtle to draw an orbit
    turtle = Turtle(visible=False)
    turtle.speed("fastest")
    turtle.pendown()
    turtle.pensize(3) # change to 1 to see the error
    turtle.pencolor(planet.color())
    
    # Prepare turtle2 to draw a graph of sweepSpeed. 
    turtle2 = Turtle(visible=False)
    turtle2.speed("fastest")
    turtle2.pendown()
    turtle2.pensize(1)
    turtle2.pencolor(planet.color())
    
    orbitScaleDownFactor = displayParameters["orbitScaleDownFactor"]
    graphStartX = displayParameters["graphStartX"]
    sweepSpeedScaleDownFactor = displayParameters["sweepSpeedScaleDownFactor"]
    timeScaleDownFactor = displayParameters["timeScaleDownFactor"]

    # JUPITER.setTimeStep(10_000)

    # Concerning Kepler's 1st law
    (x,y) = planet.position() # the right vertex of the elliptical orbit
    turtle.teleport(x/orbitScaleDownFactor, y/orbitScaleDownFactor)
    perihelion = x # shortest distance from Sun (Sun is at (0,0))
    maxX = x  
    minXsoFar = x # will be used to find the semi-major and semi-minor axes.
    maxYsoFar = y # will be used to find the semi-major and semi-minor axes.

    # Concerning Kepler's 2nd law
    (vx,vy) = planet.velocity()
    vmax = vy # maximum speed on the orbit.

    #sweepSpeedScaleDownFactor = orbitScaleDownFactor*vy
    # sweepSpeed/sweepSpeedScaleDownFactor ~ radius/orbitScaleDownFactor

    sweepSpeed0 = (x*vy-y*vx)/2
    # x*vy-y*vx is the determinant of matrix of column vectoRADIUS_SUN r,v =
    # = vector cross product  r x v
    # = the area of the parallelogram spanned by vectoRADIUS_SUN r,v =
    # = twice the area of a triangle spanned by vectoRADIUS_SUN r,v.
    # = twice the sweep speed in (m^2)/s (area swept per second).
    # Notice that the angular momentum is  r x mv.

    minSweepSpeedSoFar = sweepSpeed0
    maxSweepSpeedSoFar = sweepSpeed0
    # If the difference between minSweepSpeed and maxSweepSpeed is small,
    # Kepler's 3rd law will be confirmed.

    # Concerning Kepler's 3rd law
    tSoFar = 0 # simulated orbital period, in seconds.

    turtle2.teleport(graphStartX,
                     sweepSpeed0/sweepSpeedScaleDownFactor)
    
    # The planet starts from its perihelion.

    # Upper half of the orbit:
    done = False
    planet.setTimeStep(10)
    while True: 
        for j in range(10000): # update turtle every 100 moves
            planet.move()
            tSoFar += planet.timeStep() # update T
            (x,y) = planet.position() # update minX, maxY
            if x < minXsoFar: minXsoFar = x 
            if y > maxYsoFar: maxYsoFar = y
            (vx,vy) = planet.velocity() # update minSweepSpeed, maxSweepSpeed
            sweepSpeed = (x*vy-y*vx)/2
            if sweepSpeed < minSweepSpeedSoFar: minSweepSpeedSoFar = sweepSpeed 
            if sweepSpeed > maxSweepSpeedSoFar: maxSweepSpeedSoFar = sweepSpeed
            if y <= 0:
                done = True
                break
        turtle.goto(x/orbitScaleDownFactor, y/orbitScaleDownFactor)
        turtle2.goto(graphStartX+tSoFar/timeScaleDownFactor,
                     sweepSpeed/sweepSpeedScaleDownFactor)
        #t2.goto(x/orbitScaleDownFactor, sweepSpeed/sweepSpeedScaleDownFactor)
        if done: break

    # The planet is now at its aphelion (or a little past)

    sweepSpeed1 = sweepSpeed

    # Lower half of the orbit:
    done = False 
    while True: 
        for j in range(10000): # update turtle every 100 moves
            planet.move()
            tSoFar += planet.timeStep() # update T
            (x,y) = planet.position()
            (vx,vy) = planet.velocity() # update minSweepSpeed, maxSweepSpeed
            sweepSpeed = (x*vy-y*vx)/2
            if sweepSpeed < minSweepSpeedSoFar: minSweepSpeedSoFar = sweepSpeed 
            if sweepSpeed > maxSweepSpeedSoFar: maxSweepSpeedSoFar = sweepSpeed
            if y >= 0:
                done = True
                break
        turtle.goto(x/orbitScaleDownFactor, y/orbitScaleDownFactor)
        turtle2.goto(graphStartX+tSoFar/timeScaleDownFactor,
                     sweepSpeed/sweepSpeedScaleDownFactor)
        if done: break

    turtle2.left(90)
    turtle2.forward(3)
    turtle2.back(6)

    sweepSpeed2 = sweepSpeed

    # The planet is back at the perihelion (or a little past)

    # Concerning Kepler's 1st law:
    minX = minXsoFar 
    maxY = maxYsoFar
    a = (maxX - minX)/2 # semi-major axis - from the simulation
    b = maxY # semi-minor axis - from the simulation
    # Notice that maxX is the distance from (0,0) to the perihelion.
    orbitCenter = maxX - a
    c = sqrt(a*a - b*b) # linear eccentricity = center to focus distance.
    turtle.teleport((orbitCenter + c)/orbitScaleDownFactor,0) # draw right focus
    turtle.dot(6, planet.color())
    turtle.teleport((orbitCenter - c)/orbitScaleDownFactor,0) # draw left focus
    turtle.dot(6, planet.color()) 

    # Concerning Kepler's 2nd law: 
    minSweepSpeed = minSweepSpeedSoFar
    maxSweepSpeed = maxSweepSpeedSoFar
    kepler2discrepancy = (maxSweepSpeed - minSweepSpeed) / sweepSpeed0
    kepler2discrepancyPercent = round(kepler2discrepancy*100, 2)
    print(kepler2discrepancyPercent, "%")
    avgSweepSpeed = (minSweepSpeed + maxSweepSpeed)/2

    #print()
    #print("SweepSpeed:")
    #print(minSweepSpeed)
    #print(sweepSpeed0, 0) # starting sweepSpeed (perihelion)
                           # Why is minSweepSpeed=sweepSpeed0 ?
    #print(sweepSpeed1, 1) # aphelion
    #print(sweepSpeed2, 2) # ending sweepSpeed (perihelion)
    #print(avgSweepSpeed)
    #print(maxSweepSpeed)
    #print("====")

    # Concerning Kepler's 3rd law:
    t = tSoFar # sidereal orbital period - from the simulation.
    # T^2 / a^2 = 4*pi^2 / G(M+m) - does the simulation support this?
    #print("Kepler3:")
    lhs = t * t / (a * a * a)
    #print(lhs)
    rhs = 4 * pi * pi / (G*(MASS_SUN + planet._mass))
    #print(rhs)    
    #print(abs((lhs-rhs)/rhs))
    #print(round(abs((lhs-rhs)/rhs)*100,2), "%")

    #print(sweepSpeed0/sweepSpeedScaleDownFactor)
    #print(sweepSpeed1/sweepSpeedScaleDownFactor)
    
    return a, b, kepler2discrepancyPercent, t

#sky() # run/uncomment this before running simulateAndTest!    
#simulateAndTest(PLANET1_3, 1_000_000_000)
#simulateAndTest(PLANET1_2, 1_000_000_000)
#simulateAndTest(PLANET1_1, 1_000_000_000)
#simulateAndTest(PLANET1_0, 1_000_000_000)
#simulateAndTest(PLANET0_9, 1_000_000_000)
#simulateAndTest(PLANET0_8, 1_000_000_000)
#simulateAndTest(PLANET0_7, 1_000_000_000)

#=================================================================================
