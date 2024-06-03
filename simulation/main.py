
from simulation import *

#=================================================================================
# FUNCTIONS CALLED BY main

#---------------------------------------------------------------------------------

def simmulationSummary(planetTableData: dict, planet: SimulatedPlanet,
                       displayParameters: dict):

    orbitScaleDownFactor = displayParameters["orbitScaleDownFactor"]
    graphStartX = displayParameters["graphStartX"]
    sweepSpeedScaleDownFactor = displayParameters["sweepSpeedScaleDownFactor"]
    timeScaleDownFactor = displayParameters["timeScaleDownFactor"]

    print()
    print("-" * 79)
    print(planet.name())
    
    # Data from astronomical tables (observed and calculated)
    a0 = planetTableData["semi-major axis"]
    b0 = planetTableData["semi-minor axis"]
    t0 = planetTableData["orbital period"]
    sweepSpeed0 = planetTableData["perihelion"]*planetTableData["max speed"]/2
    
    # Predict orbit using theory
    predictedOrbit = predictByTheory(planet, displayParameters)
    a1 = predictedOrbit["semi-major axis"]
    b1 = predictedOrbit["semi-minor axis"]
    t1 = predictedOrbit["orbital period"] # sidereal
    v1 = predictedOrbit["speed at aphelion"]
    
    # Run simulation
    simulationResults = simulateAndTest(planet, displayParameters)
    a2 = simulationResults["semi-major axis"]
    b2 = simulationResults["semi-minor axis"]
    t2 = simulationResults["orbital period"]
    sweepSpeedRelativeError2 = simulationResults["sweep speed relative error"]
    v2 = simulationResults["speed at aphelion"]

    print()
    print("Semi-major axis, in km:")
    print(round(a0/KM), "- from astronomical tables,")
    print(round(a1/KM), "- predicted using theory,",
          round(abs(a1-a0)*100/a0, 3), "% relative error,")
    print(round(a2/KM), "- from simulation,",
          round(abs(a2-a0)*100/a0, 3), "% relative error.")
    print()
    print("Semi-minor axis, in km:")
    print(round(b0/KM), "- from astronomical tables,")
    print(round(b1/KM), "- predicted using theory,",
          round(abs(b1-b0)*100/b0, 3), "% relative error,")
    print(round(b2/KM), "- from simulation,",
          round(abs(b2-b0)*100/b0, 3), "% relative error.")
    print()
    print("Orbital period, in days:")
    print(round(t0/DAY), "- from astronomical tables,")
    print(round(t1/DAY), "- predicted using theory,",
          round(abs(t1-t0)*100/t0, 3), "% relative error,")
    print(round(t2/DAY), "- from simulation,",
          round(abs(t2-t0)*100/t0, 3), "% relative error.")
    print()
    print("Speed at the aphelion, in km/s:")
    print(round(v1/(KM/SEC)), "- predicted using theory,")
    print(round(v2/(KM/SEC)), "- from simulation,",
          round(abs(v2-v1)*100/v1, 3), "% relative error.")
    print()
    print("Area sweep speed, in (km^2)/s:")
    print(round(sweepSpeed0/(KM*KM/SEC)), "- a constant value by the theory,")
    print(round(sweepSpeedRelativeError2*100, 3),
          "% relative error in the simulation.")

#---------------------------------------------------------------------------------

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
    print(
"""The orbits displayed in turtle graphics are to scale,
and planets move counter-clockwise as seen from the north pole.
The Sun is not to scale; it is shown much bigger.
Other characteristics of the orbits are not modeled here -
all ellipses are shown in the same plane with the major axis on the x-axis, and
perturbations due to gravitational influence of other planets are not shown."""
         )
    
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

#planets(4) # uncomment this to simulate the inner planets.
   
#---------------------------------------------------------------------------------

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
    
    sky()

    print()
    print("-" * 79)
    print(planet.name())
        
    # Predict orbit using theory
    predictedOrbit = predictByTheory(planet, displayParameters)
    a1 = predictedOrbit["semi-major axis"]
    b1 = predictedOrbit["semi-minor axis"]
    t1 = predictedOrbit["orbital period"] # sidereal
    v1 = predictedOrbit["speed at aphelion"]
    sweepSpeed1 = predictedOrbit["perihelion"]*predictedOrbit["max speed"]/2
    
    # Run simulation
    simulationResults = simulateAndTest(planet, displayParameters)
    a2 = simulationResults["semi-major axis"]
    b2 = simulationResults["semi-minor axis"]
    t2 = simulationResults["orbital period"]
    sweepSpeedRelativeError2 = simulationResults["sweep speed relative error"]
    v2 = simulationResults["speed at aphelion"]

    print()
    print("Semi-major axis, in km:")
    print(round(a1/KM), "- predicted using theory,")
    print(round(a2/KM), "- from simulation,",
          round(abs(a2-a1)*100/a1, 3), "% relative error.")
    print()
    print("Semi-minor axis, in km:")
    print(round(b1/KM), "- predicted using theory,")
    print(round(b2/KM), "- from simulation,",
          round(abs(b2-b1)*100/b1, 3), "% relative error.")
    print()
    print("Orbital period, in days:")
    print(round(t1/DAY), "- predicted using theory,")
    print(round(t2/DAY), "- from simulation,",
          round(abs(t2-t1)*100/t1, 3), "% relative error.")
    print()
    print("Speed at the aphelion, in km/s:")
    print(round(v1/(KM/SEC)), "- predicted using theory,")
    print(round(v2/(KM/SEC)), "- from simulation,",
          round(abs(v2-v1)*100/v1, 3), "% relative error.")
    print()
    print("Area sweep speed, in (km^2)/s:")
    print(round(sweepSpeed1/(KM*KM/SEC)), "- a constant value by the theory,")
    print(round(sweepSpeedRelativeError2*100, 3),
          "% relative error in the simulation.")

#testKepler(PLANET12)

#=================================================================================

def main():
    print("Kepler's World")

    print("""
This program simulates planet movement and tests Kepler's laws.

The program draws in pink the elliptical orbit with the ellipses' foci
predicted by Kepler's laws,
then it simulates the planet movement based on Newton's laws
and draws the scaled down path of the planet,
after which it calculates and draws the foci.
If the two graphs coincide, this confirms Kepler's 1st law:
Every planet's orbit is an ellipse with the Sun in one focus.

The 2nd law says that the planet sweeps equal areas in equal amounts of time.
This is equivalent to having a constant area sweep speed.
The 3rd law says that the square of the orbital period is proportional to
the cube of the semi-major axis of the elliptical orbit
(and the constant of proportionality is known).
The program tests these laws as well.
First it draws in pink a graph showing the prediction from Kepler's laws 
of how the sweep speed depends on time,
for times from 0 to the predicted orbital period
-- this graph is (a segment of) a straight horizontal line.
Then, during the simulation, it graphs scaled down sweep speed
from the beginning to the end of a complete orbital period.
If the two graphs coincide, this confirms Kepler's 2nd and 3rd laws.

The program also prints information concerning the error of the simulation;
notice that the scaled down graphs concerning Kepler's 2nd and 3rd laws
coincide perfectly only because of limited resolution of computer screens,
while in reality there exists an error of the simulation.

In the case of actual celestial bodies
it also compares the simulated orbit to the data from astronomical tables
(which result form obeservations and theory, but are not purely observational.)
"""
)

    displayParameters = {
        "orbitScaleDownFactor" : 1_000_000_000,
        "graphStartX" : 200,
        "sweepSpeedScaleDownFactor" : 10_000_000_000_000,
        "timeScaleDownFactor" : 100_000
    }

    while True:
        print("-" * 79)
        print("""
1. Inner planets: Mercury, Venus, Earth and Mars.

2. A made-up celestial body with the same mass and perihelion as Earth
   but with the maximal speed 20% bigger than Earth.

3. Exit.
          """)

# add this menu option:
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

#---------------------------------------------------------------------------------

if __name__ == "__main__": main()

#=================================================================================
