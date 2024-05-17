
from simulation import *

#=================================================================================
# FUNCTIONS CALLED BY main

def simmulationSummary(planetData: dict, planet: SimulatedPlanet,
                       displayParameters):

    print("\n")
    print(planet.name(), "orbital period in days:")

    orbitScaleDownFactor = displayParameters["orbitScaleDownFactor"]
    graphStartX = displayParameters["graphStartX"]
    sweepSpeedScaleDownFactor = displayParameters["sweepSpeedScaleDownFactor"]
    timeScaleDownFactor = displayParameters["timeScaleDownFactor"]
    
    t = planetData["orbital period"]
    a = planetData["semi-major axis"]
    b = planetData["semi-minor axis"]
    c = planetData["linear eccentricity"]
    vmax = planetData["max speed"]
    perihelion = planetData["perihelion"]

    sweepSpeed0 = perihelion * vmax / 2
    
    print(round(t/DAY,2), "- actual")
    
    predictByTheory(planet, displayParameters)

    #mu = G*(MASS_SUN+planet.mass()) # gravitational parameter
    #a = perihelion*mu / (2*mu-perihelion*vmax*vmax) # semi-major axis:
    #c = a - perihelion # the linear eccentricity, i.e. center-to-focus distance
    #aphelion = a + c # predicted aphelion distance
    #b = sqrt(aphelion*perihelion)  # predicted se1mi-minor axis
    #t = 2*pi*sqrt(a*a*a/mu) # predicted orbital period in s

    turtle2 = Turtle(visible=False)
    turtle2.speed("fastest")
    turtle2.pendown()
    turtle2.pensize(1)
    turtle2.pencolor("pink")

    turtle2.teleport(graphStartX,
                     sweepSpeed0 / sweepSpeedScaleDownFactor)
    turtle2.speed(5)

    # draw a pink hirizontal line with a blunt end 
    turtle2.goto(graphStartX+t/timeScaleDownFactor,
                     sweepSpeed0 /sweepSpeedScaleDownFactor)
    turtle2.left(90)
    turtle2.forward(3)
    turtle2.back(6)
    
    ts = simulateAndTest(planet, displayParameters)[3]
    print(round(abs(ts-t)*100/t, 2), "% error")


def planets(n: int, displayParameters):
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
    
    orbitScaleDownFactor = displayParameters["orbitScaleDownFactor"]
    graphStartX = displayParameters["graphStartX"]
    sweepSpeedScaleDownFactor = displayParameters["sweepSpeedScaleDownFactor"]
    timeScaleDownFactor = displayParameters["timeScaleDownFactor"]
    
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
        simmulationSummary(MERCURY_DATA, MERCURY, displayParameters)
    if n>=2:
        simmulationSummary(VENUS_DATA, VENUS, displayParameters)
    if n>=3:
        simmulationSummary(EARTH_DATA, EARTH, displayParameters)
    if n>=4:
        simmulationSummary(MARS_DATA, MARS, displayParameters)
    if n>=5:
        simmulationSummary(JUPITER_DATA, JUPITER, displayParameters)
    if n>=6:
        simmulationSummary(SATURN_DATA, SATURN, displayParameters)
    if n>=7:
        simmulationSummary(URANUS_DATA, URANUS, displayParameters)
    if n>=8:
        simmulationSummary(NEPTUNE_DATA, NEPTUNE, displayParameters)

#-------
#inneRADIUS_SUNimulatedPlanets() # uncomment this to simulate the inner planets.

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
def testKepler(planet: SimulatedPlanet, displayParameters):
    """Precondition: planet position (x,y) must have x>0, y=0,
                    and velocity (vx,vy) must have vx=0, vy>0.
       So, the planet must be in the right vertex of its elliptical orbit.A computer simulation of the orbit of the planet, resulting from the
       continuing local effect of the Newton's law of gravity.
       Tests if the simulated planet obey (global) Kepler's laws 1 and 3.
       Every step in the simulation is done for the actual planet
       (such as Mars, with its actual mass, perihelion and actual max speed)
       and only later it is scaled down to be displayed in turtle graphics.
       Real distances in meteRADIUS_SUN will be divided by orbitScaleDownFactor
       before being given to the turtle.
       This function creates its own canvas-sky with the Sun,
       and draws the predicted ellipse of the orbit with the foci,
       before starting simulation.
       It produces output, besides displaying turtle graphics.
    """
    
    orbitScaleDownFactor = displayParameters["orbitScaleDownFactor"]
    graphStartX = displayParameters["graphStartX"]
    sweepSpeedScaleDownFactor = displayParameters["sweepSpeedScaleDownFactor"]
    timeScaleDownFactor = displayParameters["timeScaleDownFactor"]
    
    print("The orbit displayed in turtle graphics is to scale.")
    sky()
    print("\nTesting Kepler's 1st and 3rd laws")
    m = planet.mass()  # mass
    vmax = planet.velocity()[1] # maximal speed (at Perihelion)
    perihelion = planet.position()[0]  # shortest distance from Sun
    mu = G*(MASS_SUN+m) # gravitational parameter
    a = perihelion*mu / (2*mu-perihelion*vmax*vmax) # semi-major axis:
    c = a - perihelion # the linear eccentricity, i.e. center-to-focus distance
    aphelion = a + c # Aphelion distance (biggest distance from the Sun)
    b = sqrt(aphelion*perihelion)  # semi-minor axis
    t = sqrt(4*pi*pi*a*a*a / mu) # orbital period in s
    print("Planet's orbital period in days:")
    print(round(t/DAY,2), "- predicted by the theory")
    predictByTheory(planet, displayParameters)
    #drawEllipse(a/orbitScaleDownFactor, b/orbitScaleDownFactor, c/orbitScaleDownFactor) # predicted orbit
    ts = simulateAndTest(planet, displayParameters)[3] # simulation. TS - orbital period form simulation.
    print(round(ts/DAY,2), "- from the simulation")
    print( round(abs(ts-t)*100/t, 2), "% discrepancy")

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

    displayParameters = {
        "orbitScaleDownFactor" : 1_000_000_000,
        "graphStartX" : 200,
        "sweepSpeedScaleDownFactor" : 10_000_000_000_000,
        "timeScaleDownFactor" : 100_000
    }

    while True:
        print("""
1. Inner planets: Mercury, Venus, Earth and Mars.
   
2. A made-up celestial body with the same mass and perihelion as Earth
   but with the maximal speed 20% bigger than Earth.

3. Exit.
          """)
#2. Eight planets: Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune.

        choice = input("Enter your choice (1-3): ")
        if choice == "":
            return
        if   choice[0] == '1':
            planets(4, displayParameters) # The 4 inner planets.
        #elif choice[0] == '2':
        #    planets(8, displayParameters) # All 8 planets.
        elif choice[0] == '2':
            testKepler(PLANET1_2, displayParameters)
        else:
            print("Bye")
            return

if __name__ == "__main__": main()

#=================================================================================

#sky()
#simulateAndTest(PLANET1_2)
#simulateAndTest(JUPITER)
#simulateAndTest(SATURN)
#simulateAndTest(URANUS)
#simulateAndTest(NEPTUNE)
