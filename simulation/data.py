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

A slideshow is also available: Discovering the Mechanics of the Solar System.

These materials can be used as a basis for a discussion of concepts of
methodology of science, as outlined in README.md.
"""

# TO DO:
# experiment1() - planets with masses m, 2m, 3m and the same perihelion.
# experiment2() - plenets with the same perihelion and max speed v, 2v, 3v
# planets with the same orbtital period and different perihelions
# ...

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

MASS_SUN = 1.98847e30   # Mass in kg.
RADIUS_SUN = 696_000*KM # Radius in m - NOT USED
MU = G*MASS_SUN         # Gravitational parameter with respect to the Sun only.

#---------------------------------------------------------------------------------

def bcsFromAPm(aphelion, perihelion, m):
    """aphelion in m, perihelion in m, m - mass in kg.
       Returns semi-minor axis, linear eccentricity and max and min speed.
    """
    a = (aphelion + perihelion) / 2 # semi-major axis in m.
    b = sqrt(aphelion * perihelion) # semi-minor axis in m.
    c = (aphelion - perihelion) / 2 # linear eccentricity = center-to-focus distance, in m.
    mu = G * (MASS_SUN + m) # gravitational parameter (with respect to the Sun)
    sMax = sqrt(mu * (2/perihelion - 1/a)) # speed at perihelion - max speed
    sMin = sqrt(mu * (2/aphelion - 1/a)) # speed at aphelion - min speed
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
MERCURY_DATA["orbital period"] = 87.9691*DAY # sidereal
MERCURY_DATA["radius"] = 2439.7*KM

BCS = bcsFromAPm(MERCURY_DATA["aphelion"],
                 MERCURY_DATA["perihelion"],
                 MERCURY_DATA["mass"])

MERCURY_DATA["semi-minor axis"] = BCS[0]
MERCURY_DATA["linear eccentricity"] = BCS[1] # center-to-focus distance
MERCURY_DATA["max speed"] = BCS[2] # at the perihelion
MERCURY_DATA["min speed"] = BCS[3] # at the aphelion

#---------------------------------------------------------------------------------

VENUS_DATA = {}
VENUS_DATA["mass"] = 4.8675e24
VENUS_DATA["aphelion"]        = 108_940_000*KM
VENUS_DATA["perihelion"]      = 107_480_000*KM
VENUS_DATA["semi-major axis"] = 108_210_000*KM
VENUS_DATA["eccentricity"] = 0.006772   
VENUS_DATA["avg speed"] = 35.02*KM/SEC 
VENUS_DATA["orbital period"] = 224.701*DAY # sidereal
VENUS_DATA["radius"] = 6051.8*KM

BCS = bcsFromAPm(VENUS_DATA["aphelion"],
                 VENUS_DATA["perihelion"],
                 VENUS_DATA["mass"])

VENUS_DATA["semi-minor axis"] = BCS[0]
VENUS_DATA["linear eccentricity"] = BCS[1] # center-to-focus distance
VENUS_DATA["max speed"] = BCS[2] # at the perihelion
VENUS_DATA["min speed"] = BCS[3] # at the aphelion

#---------------------------------------------------------------------------------

EARTH_DATA = {}
EARTH_DATA["mass"] = 5.972168e24
EARTH_DATA["aphelion"]        = 152_097_597*KM
EARTH_DATA["perihelion"]      = 147_098_450*KM
EARTH_DATA["semi-major axis"] = 149_598_023*KM
EARTH_DATA["eccentricity"] = 0.0167086
EARTH_DATA["avg speed"] = 29.7827*KM/SEC
EARTH_DATA["orbital period"] = 365.2563630040*DAY # sidereal
EARTH_DATA["radius"] = 6371*KM

BCS = bcsFromAPm(EARTH_DATA["aphelion"],
                 EARTH_DATA["perihelion"],
                 EARTH_DATA["mass"])

EARTH_DATA["semi-minor axis"] = BCS[0]
EARTH_DATA["linear eccentricity"] = BCS[1] # center-to-focus distance
EARTH_DATA["max speed"] = BCS[2] # at the perihelion
EARTH_DATA["min speed"] = BCS[3] # at the aphelion

#---------------------------------------------------------------------------------

MARS_DATA = {}
MARS_DATA["mass"] = 6.4171e23
MARS_DATA["aphelion"]        = 249_261_000*KM
MARS_DATA["perihelion"]      = 206_650_000*KM
MARS_DATA["semi-major axis"] = 227_939_366*KM
MARS_DATA["eccentricity"] = 0.0934
MARS_DATA["avg speed"] = 24.07*KM/SEC
MARS_DATA["orbital period"] = 686.980*DAY # sidereal
MARS_DATA["radius"] = 3389.5*KM

BCS = bcsFromAPm(MARS_DATA["aphelion"],
                 MARS_DATA["perihelion"],
                 MARS_DATA["mass"])

MARS_DATA["semi-minor axis"] = BCS[0]
MARS_DATA["linear eccentricity"] = BCS[1] # center-to-focus distance
MARS_DATA["max speed"] = BCS[2] # at the perihelion
MARS_DATA["min speed"] = BCS[3] # at the aphelion

#---------------------------------------------------------------------------------

JUPITER_DATA = {}
JUPITER_DATA["mass"] = 1.8982e27
JUPITER_DATA["aphelion"]        = 816_363_000*KM
JUPITER_DATA["perihelion"]      = 740_595_000*KM
JUPITER_DATA["semi-major axis"] = 778_479_000*KM
JUPITER_DATA["eccentricity"] = 0.0489
JUPITER_DATA["avg speed"] = 13.07*KM/SEC
JUPITER_DATA["orbital period"] = 4_332.59*DAY # sidereal
JUPITER_DATA["radius"] = 69_911*KM

BCS = bcsFromAPm(JUPITER_DATA["aphelion"],
                 JUPITER_DATA["perihelion"],
                 JUPITER_DATA["mass"])

JUPITER_DATA["semi-minor axis"] = BCS[0]
JUPITER_DATA["linear eccentricity"] = BCS[1] # center-to-focus distance
JUPITER_DATA["max speed"] = BCS[2] # at the perihelion
JUPITER_DATA["min speed"] = BCS[3] # at the aphelion

#---------------------------------------------------------------------------------

SATURN_DATA = {}
SATURN_DATA["mass"] = 5.6834e26
SATURN_DATA["aphelion"]        = 1_514_500_000*KM
SATURN_DATA["perihelion"]      = 1_352_550_000*KM
SATURN_DATA["semi-major axis"] = 1_433_530_000*KM
SATURN_DATA["eccentricity"] = 0.0565
SATURN_DATA["avg speed"] = 9.68*KM/SEC
SATURN_DATA["orbital period"] = 10_755.70*DAY # sidereal 
SATURN_DATA["radius"] = 58_232*KM

BCS = bcsFromAPm(SATURN_DATA["aphelion"],
                 SATURN_DATA["perihelion"],
                 SATURN_DATA["mass"])

SATURN_DATA["semi-minor axis"] = BCS[0]
SATURN_DATA["linear eccentricity"] = BCS[1] # center-to-focus distance
SATURN_DATA["max speed"] = BCS[2] # at the perihelion
SATURN_DATA["min speed"] = BCS[3] # at the aphelion

#---------------------------------------------------------------------------------

URANUS_DATA = {}
URANUS_DATA["mass"] = 8.6810e25
URANUS_DATA["aphelion"]        = 3_006_390_000*KM
URANUS_DATA["perihelion"]      = 2_735_560_000*KM
URANUS_DATA["semi-major axis"] = 2_870_972_000*KM
URANUS_DATA["eccentricity"] = 0.04717
URANUS_DATA["avg speed"] = 6.80*KM/SEC
URANUS_DATA["orbital period"] = 30_688.5*DAY # sidereal 
URANUS_DATA["radius"] = 25_362*KM

BCS = bcsFromAPm(URANUS_DATA["aphelion"],
                 URANUS_DATA["perihelion"],
                 URANUS_DATA["mass"])

URANUS_DATA["semi-minor axis"] = BCS[0]
URANUS_DATA["linear eccentricity"] = BCS[1] # center-to-focus distance
URANUS_DATA["max speed"] = BCS[2] # at the perihelion
URANUS_DATA["min speed"] = BCS[3] # at the aphelion

#---------------------------------------------------------------------------------

NEPTUNE_DATA = {}
NEPTUNE_DATA["mass"] = 1.02409e26
NEPTUNE_DATA["aphelion"]        = 4_540_000_000*KM
NEPTUNE_DATA["perihelion"]      = 4_460_000_000*KM
NEPTUNE_DATA["semi-major axis"] = 4_500_000_000*KM
NEPTUNE_DATA["eccentricity"] = 0.008678
NEPTUNE_DATA["avg speed"] = 5.43*KM/SEC
NEPTUNE_DATA["orbital period"] = 60_195*DAY # sidereal 
NEPTUNE_DATA["radius"] = 24_622*KM

BCS = bcsFromAPm(NEPTUNE_DATA["aphelion"],
                 NEPTUNE_DATA["perihelion"],
                 NEPTUNE_DATA["mass"])

NEPTUNE_DATA["semi-minor axis"] = BCS[0]
NEPTUNE_DATA["linear eccentricity"] = BCS[1] # center-to-focus distance
NEPTUNE_DATA["max speed"] = BCS[2] # at the perihelion
NEPTUNE_DATA["min speed"] = BCS[3] # at the aphelion

#---------------------------------------------------------------------------------

# Compare sweepSpeed speeds of the closest and farthest planet
#print(NEPTUNE_DATA["perihelion"]*NEPTUNE_DATA["max speed"]/
#      MERCURY_DATA["perihelion"]/MERCURY_DATA["max speed"])
# 9.01

# Compare sweepSpeed speeds of the Mercury and Mars
#print(MARS_DATA["perihelion"]*MARS_DATA["max speed"]/
#      MERCURY_DATA["perihelion"]/MERCURY_DATA["max speed"])
# 2.02

# The quantity from Kepler's 3rd law:
#P = MARS_DATA["orbital period"]
#A = MARS_DATA["semi-major axis"]
#print(P*P*P/A/A)
# 4.02

#---------------------------------------------------------------------------------
# Pluto's data is not used in the current veRADIUS_SUNion of the program.

# While the planets' orbits are close to the ecliptic,
# Pluto's orbit has a significant inclination.
# This program shows the orbits on the same plane. 
# Showing the planets and Plutothis way would be
# greatly inaccurate and lead to misconceptions, for instance
# the visualization would show that Pluto's orbit inteRADIUS_SUNects Uranus' orbit.

PLUTO_DATA = {}
PLUTO_DATA["mass"] = 1.303e22
PLUTO_DATA["aphelion"]        = 7_375_930_000*KM
PLUTO_DATA["perihelion"]      = 4_436_820_000*KM
PLUTO_DATA["semi-major axis"] = 5_906_380_000*KM
PLUTO_DATA["eccentricity"] = 0.2488
PLUTO_DATA["avg speed"] = 4.743*KM/SEC
PLUTO_DATA["orbital period"] = 90_560*DAY # sidereal
PLUTO_DATA["radius"] = 1_188.3*KM

BCS = bcsFromAPm(PLUTO_DATA["aphelion"],
                 PLUTO_DATA["perihelion"],
                 PLUTO_DATA["mass"])

PLUTO_DATA["semi-minor axis"] = BCS[0]
PLUTO_DATA["linear eccentricity"] = BCS[1] # center-to-focus distance
PLUTO_DATA["max speed"] = BCS[2] # at the perihelion
PLUTO_DATA["min speed"] = BCS[3] # at the aphelion

#---------------------------------------------------------------------------------
# Halley comet's data is not used in the current veRADIUS_SUNion of the program.

# Halley's comet has a significant inclination with respect to the excliptic.
# and significant inclination to the plane of Pluto's orbit.

# Halley's Comet, unlike the planets, moves clockwise around the Sun
# when viewed form the north pole.

HALLEY_COMET_DATA = {}
HALLEY_COMET_DATA["mass"] = 2.2e14
HALLEY_COMET_DATA["aphelion"]        = 35.14*AU
HALLEY_COMET_DATA["perihelion"]      =  0.59278*AU
HALLEY_COMET_DATA["semi-major axis"] = 17.737*AU
HALLEY_COMET_DATA["eccentricity"] = 0.96658 
# HALLEY_COMET_DATA["avg speed"] = 
HALLEY_COMET_DATA["orbital period"] =  74.7*YR # sidereal
HALLEY_COMET_DATA["radius"] = 5.5*KM

BCS = bcsFromAPm(HALLEY_COMET_DATA["aphelion"],
                 HALLEY_COMET_DATA["perihelion"],
                 HALLEY_COMET_DATA["mass"])

HALLEY_COMET_DATA["semi-minor axis"] = BCS[0]
HALLEY_COMET_DATA["linear eccentricity"] = BCS[1] # center-to-focus distance
HALLEY_COMET_DATA["max speed"] = BCS[2] # at the perihelion
HALLEY_COMET_DATA["min speed"] = BCS[3] # at the aphelion

#---------------------------------------------------------------------------------

del BCS

#=================================================================================
