
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
    screen.screensize(2100,950)
    screen.title("Kepler's world")
    screen.bgcolor(skyColor)
    if showSun:
        turtle = Turtle(visible=False)
        turtle.dot(10, "white") # The Sun is at (0,0) and is white.

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
       Note: the parameters are in turtle canvas units, not in meters.
       Note: the resulting ellipse is not perfect for two reasons:
             1. floating point calculation has limited precision;
             2. the computer screen has limited resolution.
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
    """Extracts perihelion and maxSpeed from the simulated planet,
       predicts and draws in pink its orbit; predicts the orbital period,
       predicts that sweepSpeed is constant,
       draws a graph of such a constant function depending on time,
       for time from 0 to the orbital predicted period.
       Returns predicted parameters of the ellipse.
    """
    # Display parameters
    orbitScaleDownFactor = displayParameters["orbitScaleDownFactor"]
    graphStartX = displayParameters["graphStartX"]
    sweepSpeedScaleDownFactor = displayParameters["sweepSpeedScaleDownFactor"]
    timeScaleDownFactor = displayParameters["timeScaleDownFactor"]

    # Calculate orbit parameters
    perihelion = planet.perihelion()
    vmax = planet.maxSpeed()
    mu = G * (MASS_SUN + planet.mass()) # gravitational parameter
    a = perihelion * mu / (2 * mu-perihelion * vmax * vmax) # semi-major axis
    c = a - perihelion # the linear eccentricity, i.e. center-to-focus distance
    aphelion = a + c # predicted aphelion distance
    b = sqrt(aphelion * perihelion)  # predicted semi-minor axis
    t = 2 * pi * sqrt(a * a * a / mu) # predicted sidereal orbital period in s
    sweepSpeedAtPerihelion = perihelion * vmax / 2
    minSpeed = sqrt(mu * (2/aphelion - 1/a))

    predictedOrbit = {}
    predictedOrbit["perihelion"] = perihelion
    predictedOrbit["max speed"] = vmax
    predictedOrbit["semi-major axis"] = a
    predictedOrbit["semi-minor axis"] = b
    predictedOrbit["aphelion"] = aphelion
    predictedOrbit["orbital period"] = t # sidereal
    predictedOrbit["min speed"] = minSpeed
    predictedOrbit["speed at aphelion"] = minSpeed

    # Draw scaled down predicted orbit
    drawEllipse(a/orbitScaleDownFactor, b/orbitScaleDownFactor,
                c/orbitScaleDownFactor)

    # Draw a horizontal segment with length proportional to the orbital period
    # at the hight propostional to sweepSpeedAtPerihelion.
    turtle = Turtle(visible=False)
    turtle.speed("fastest")
    turtle.pendown()
    turtle.pensize(1)
    turtle.speed(3)
    turtle.pencolor("pink")
    # Put turtle at the starting position
    turtle.teleport(graphStartX,
                    sweepSpeedAtPerihelion / sweepSpeedScaleDownFactor)
    # Draw left end:
    turtle.left(90)
    turtle.forward(3)
    turtle.back(6)
    turtle.forward(3)
    turtle.right(90)
    # Draw horizontal line:
    turtle.goto(graphStartX + t / timeScaleDownFactor,
                sweepSpeedAtPerihelion / sweepSpeedScaleDownFactor)
    # Draw right end:
    turtle.left(90)
    turtle.forward(3)
    turtle.back(6)

    return predictedOrbit

#sky(showSun=False) # run/uncomment this before running predictByTheory

#---------------------------------------------------------------------------------

def simulate(planet: SimulatedPlanet):
    # Just draws the orbit, does not test Kepler's laws.
    # This function is not used by the top level functions in the program.
    # It is given here as a stepping stone to understand simulateAndTest below.
    """Precondition: planet position (x,y) must have x>0, y=0,
                     and velocity (vx,vy) must have vx=0, vy>0.
       So, the planet must be in the apsis (perihelion or aphelion)
       and moving up / counter-clockwise.
       This function calculates the orbit resulting from the continuing local
       effect of the Newton's law of gravity, assuming the Sun is at (0,0).
       Every step in the simulation is done for the actual planet
       (such as Mars, with its actual mass, perihelion and actual max speed)
       and only later it is scaled down to be displayed in turtle graphics.
       Real distances in meters will be divided by orbitScaleDownFactor
       before being given to the turtle.
       Note: make sure to create canvas before this function is called.
    """
    # Prepare turtle.
    turtle = Turtle(visible=False)
    turtle.speed("fastest")
    turtle.pendown()
    turtle.pensize(1)
    turtle.pencolor(planet.color())

    orbitScaleDownFactor =  1_000_000_000 # suitable for inner planets

    planet.reset()
    (x,y) = planet.position() # the right vertex of the elliptical orbit
    turtle.teleport(x/orbitScaleDownFactor, y/orbitScaleDownFactor)

    # The planet starts from an apsis (perihelion or aphelion)
    
    # UPPER PART OF THE ORBIT:
    done = False
    while True: 
        for j in range(10_000): # update turtle every 10_000 moves
            planet.move()
            (x,y) = planet.position()
            if y <= 0:
                done = True
                break
        turtle.goto(x/orbitScaleDownFactor, y/orbitScaleDownFactor)
        if done: break

    # The planet is now at apsis2 (aphelion or perihelion) or a little past.

    # LOWER PART OF THE ORBIT:
    done = False
    while True: 
        for j in range(10_000): # update turtle every 10_000 moves
            planet.move()
            (x,y) = planet.position()
            if y >= 0:
                done = True
                break
        turtle.goto(x/orbitScaleDownFactor, y/orbitScaleDownFactor)
        if done: break

    # The planet is back at the initial apsis or a little past.

#sky() # run/uncomment this before running simulate!
#simulate(PLANET1_2)

#-----------------------------------------------------------------------------

def simulateAndTest(planet: SimulatedPlanet, displayParameters):
    """Precondition: planet position (x,y) must have x>0, y=0,
                     and velocity (vx,vy) must have vx=0, vy>0.
       The planet must be at the perihelion and move up & counter-clockwise.
       This function calculates the orbit resulting from the continuing local
       effect of the Newton's law of gravity, assuming the Sun is at (0,0).
       Every step in the simulation is done for the actual planet
       (such as Mars, with its actual mass, perihelion and actual max speed)
       and only later it is scaled down to be displayed in turtle graphics.
       Real distances in meters will be divided by orbitScaleDownFactor
       before being given to the turtle.
       The function tests Kepler's laws and provides output in turtle graphics
       and returns a dictionary containing simulation results.
       Note: make sure to create canvas before this function is called.
    """
    # Prepare turtle for drawing the orbit and turtle2 for sweep speed graph.
    turtle = Turtle(visible=False)
    turtle.speed("fastest")
    turtle.pendown()
    turtle.pensize(3) # You can change it to 1 to see the error
    turtle.pencolor(planet.color())
    turtle2 = Turtle(visible=False)
    turtle2.speed("fastest")
    turtle2.pendown()
    turtle2.pensize(1)
    turtle2.pencolor(planet.color())

    # Get display parameters
    orbitScaleDownFactor = displayParameters["orbitScaleDownFactor"]
    graphStartX = displayParameters["graphStartX"]
    sweepSpeedScaleDownFactor = displayParameters["sweepSpeedScaleDownFactor"]
    timeScaleDownFactor = displayParameters["timeScaleDownFactor"]

    planet.reset()

    # Concerning Kepler's 1st law
    (x,y) = planet.position() # the right vertex of the elliptical orbit
    if y != 0: raise ValueError("Planet's starting position must have y==0.")
    if x <= 0: raise ValueError("Planet's starting position must have x>0.")
    turtle.teleport(x/orbitScaleDownFactor, y/orbitScaleDownFactor)
    maxX = perihelion = x # shortest distance from the Sun, which is at (0,0)
    minXsoFar = x # will be used to find the semi-major axis.
    maxYsoFar = y # will be used to find the semi-minor axis.

    # Concerning Kepler's 2nd law
    (vx,vy) = planet.velocity()
    if vx != 0: raise ValueError("Planet's starting velocity must have vx==0.")
    if vy <= 0: raise ValueError("Planet's starting velocity must have vy>0.")
    # With the initial vy>0 the planet moves counter-clockwise araind the Sun;
    # otherwise it would move clockwise and the sweep speed would be negative.
    vAtPerihelion = vy

    sweepSpeedAtPerihelion = (x*vy-y*vx)/2
    # x*vy-y*vx is the determinant of a matrix of column vectors r,v =
    # = vector cross product:  r x v
    # = the area of the parallelogram spanned by vectors r,v =
    # = twice the area of a triangle spanned by vectors r,v.
    # = twice the sweep speed in (m^2)/s (area swept per second).
    # Notice that the angular momentum is  r x mv.

    minSweepSpeedSoFar = sweepSpeedAtPerihelion
    maxSweepSpeedSoFar = sweepSpeedAtPerihelion
    # If the difference between min and max sweep speed remains small,
    # Kepler's 3rd law will be confirmed.

    # Concerning Kepler's 3rd law
    time = 0 # simulated orbital period, in seconds.

    # Draw left end of the sweep speed graph:
    turtle2.teleport(graphStartX,
                     sweepSpeedAtPerihelion/sweepSpeedScaleDownFactor)
    turtle2.left(90)
    turtle2.forward(3)
    turtle2.back(6)
    turtle2.forward(3)
    turtle2.right(90)
    
    # The planet starts from its perihelion.

    # UPPER PART OF THE ORBIT:
    done = False
    while True: 
        for j in range(10_000): # update turtle every 10_000 moves
            planet.move()
            time += planet.timeStep() # update T
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
        turtle2.goto(graphStartX+time/timeScaleDownFactor,
                     sweepSpeed/sweepSpeedScaleDownFactor)
        if done: break

    # The planet is now at its aphelion (or a little past)

    vAtAphelion = vy

    # LOWER PART OF THE ORBIT:
    done = False 
    while True: 
        for j in range(10_000): # update turtle every 10_000 moves
            planet.move()
            time += planet.timeStep() # update T
            (x,y) = planet.position()
            (vx,vy) = planet.velocity() # update minSweepSpeed, maxSweepSpeed
            sweepSpeed = (x*vy-y*vx)/2
            if sweepSpeed < minSweepSpeedSoFar: minSweepSpeedSoFar = sweepSpeed 
            if sweepSpeed > maxSweepSpeedSoFar: maxSweepSpeedSoFar = sweepSpeed
            if y >= 0:
                done = True
                break
        turtle.goto(x/orbitScaleDownFactor, y/orbitScaleDownFactor)
        turtle2.goto(graphStartX+time/timeScaleDownFactor,
                     sweepSpeed/sweepSpeedScaleDownFactor)
        if done: break

    # Draw the right end of the sweep speed graph:
    turtle2.left(90)
    turtle2.forward(3)
    turtle2.back(6)

    # SIMULATION FINISHED

    # The planet is back at the perihelion (or a little past)

    # Concerning Kepler's 1st law:
    minX = minXsoFar
    maxY = maxYsoFar
    aphelion = abs(minX) # maximal distance from the Sun
    a = (maxX - minX) / 2 # semi-major axis - from the simulation
    b = maxY # semi-minor axis - from the simulation
    orbitCenter = maxX - a # this is the x-coordiante, the y-coordinate = 0
    c = sqrt(a*a - b*b) # linear eccentricity = center to focus distance.
    turtle.teleport((orbitCenter + c)/orbitScaleDownFactor,0) # draw right focus
    turtle.dot(6, planet.color())
    turtle.teleport((orbitCenter - c)/orbitScaleDownFactor,0) # draw left focus
    turtle.dot(6, planet.color())

    # Concerning Kepler's 2nd law: 
    minSweepSpeed = minSweepSpeedSoFar
    maxSweepSpeed = maxSweepSpeedSoFar
    sweepSpeedRelativeError = \
        (maxSweepSpeed - minSweepSpeed) / sweepSpeedAtPerihelion

    # Concerning Kepler's 3rd law:
    t = time # sidereal orbital period - from the simulation.
    
    #print(round(abs((lhs-rhs)/rhs)*100,2), "%")

    simulationResults = {}
    simulationResults["semi-major axis"] = a
    simulationResults["semi-minor axis"] = b
    simulationResults["orbital period"] = t
    simulationResults["sweep speed relative error"] = \
                             sweepSpeedRelativeError
    simulationResults["aphelion"] = aphelion
    simulationResults["speed at aphelion"] = abs(vAtAphelion)

    return simulationResults

#sky() # run/uncomment this before running simulateAndTest!    
#simulateAndTest(PLANET1_3, 1_000_000_000)
#simulateAndTest(PLANET1_2, 1_000_000_000)
#simulateAndTest(PLANET1_1, 1_000_000_000)
#simulateAndTest(PLANET1_0, 1_000_000_000)
#simulateAndTest(PLANET0_9, 1_000_000_000)
#simulateAndTest(PLANET0_8, 1_000_000_000)
#simulateAndTest(PLANET0_7, 1_000_000_000)

#=================================================================================
