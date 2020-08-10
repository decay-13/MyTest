import ephem
from pyorbital import tlefile
import sys
import datetime
import math
import time
from datetime import timedelta
import matplotlib.pyplot as plt
import numpy as np
import Stepper
from matplotlib.pyplot import rc, grid, figure, plot, rcParams, savefig
from math import radians
import RPi.GPIO as GPIO

#from spacePi import Azimuth_Stepper_Motor

GPIO.setmode(GPIO.BCM)
Azimuth_Stepper_Motor = Stepper.Motor([21, 20, 16, 26])
Elevation_Stepper_Motor = Stepper.Motor([4,17,27,22])
Azimuth_Stepper_Motor.rpm = 5
Elevation_Stepper_Motor.rpm = 5
Azimuth_Stepper_Motor.deactivate()
Elevation_Stepper_Motor.deactivate()

def convert_radians_to_degrees(radians) :
    return radians * 180 / math.pi
def dms_to_dd(d, m, s):
    dd = d + float(m)/60 + float(s)/3600
    return dd

def convert(tude):
   # multiplier = 1 if tude[-1] in ['N', 'E'] else -1
    return sum(float(x) / 60 ** n for n, x in enumerate(tude[:-1].split('-')))

def get_planet_colour(planet):
    if (planet.name == 'Moon'):
        return 'blue'
    if (planet.name == 'Mars'):
        return 'red'
    if (planet.name == 'Saturn'):
        return 'cyan'
    if (planet.name == 'Jupiter'):
        return 'purple'
    if (planet.name == 'Neuptun'):
        return 'black'
    if (planet.name == 'Venus'):
        return 'yellow'
    if (planet.name == 'Mercury'):
        return 'silver'

    return 'white'


def predict_planet_orbit(obs, planet):
    planet_path = list()
    min_alt = 10. * math.pi / 180
    obs.date = datetime.datetime.utcnow()
    planet.compute(obs);
    #if planet.alt > min_alt:
    #    print planet.name, " is visable: Alt %s Azimuth %s Time %s" % (planet.alt, planet.az, ephem.localtime(obs.date))
    next_setting = obs.next_setting(planet)
    while obs.date < next_setting:
            obs.date += ephem.minute
            planet.compute(obs)
            print(planet.alt, planet.az)
            temp_tuple = (planet.alt), (planet.az)
            planet_path.append(temp_tuple)

    return planet_path



u = ephem.Uranus()
u.compute('1781/3/13')
print('%s %s %s' % (u.ra, u.dec, u.mag))
#5:35:45.28 23:32:54.1 5.6
print(ephem.constellation(u))

u.compute()
obs = ephem.Observer()
obs.lat = '-34.73'
obs.long = '138.66'
print('%s %s %s' % (u.ra, u.dec, u.mag))

d = ephem.Date('1984/12/21 15:00')
d = ephem.Date(datetime.datetime.now())
print (ephem.Date)
ephem.localtime(d)
print(ephem.localtime(d).ctime())

#using utcnow as ephem uses utc, convert using localtime
obs.date = datetime.datetime.utcnow()
print (obs.date)
sun, moon = ephem.Sun(), ephem.Moon()
venus, mercury, saturn, jupiter = ephem.Venus(), ephem.Mercury(), ephem.Saturn(), ephem.Jupiter()
print ("Sun sets at %s " % (ephem.localtime(obs.next_setting(sun))))
print ("Sunrise at %s " % (ephem.localtime(obs.next_rising(sun))))

print ("venus next rising %s " % ephem.localtime(obs.next_rising(venus)))
print ("mercury next rising %s " % ephem.localtime(obs.next_rising(mercury)))
print ("saturn next rising %s " % ephem.localtime(obs.next_rising(saturn)))
print ("jupiter next rising %s " % ephem.localtime(obs.next_rising(jupiter)))

#lets track the iss
print("tracking iss")
track_point = 0
track_point_limit = 3600
tlefile.TLE_URLS = (
    'http://celestrak.com/NORAD/elements/stations.txt',
)
iss_tle = tlefile.read('ISS (ZARYA)')
iss = ephem.readtle("ISS (ZARYA)", iss_tle.line1, iss_tle.line2)
Azimuth_Stepper_Motor.set_rpm(20);
Elevation_Stepper_Motor.set_rpm(40);
#insert test, small cog full rev, large cog full rev, 0, 90, 180, 270,360
Azimuth_Stepper_Motor.move_steps(Azimuth_Stepper_Motor.steps_per_full_rev)
Elevation_Stepper_Motor.move_steps(Elevation_Stepper_Motor.steps_per_full_rev)
Azimuth_Stepper_Motor.move_steps(Azimuth_Stepper_Motor.full_large_cog_rev)
Elevation_Stepper_Motor.move_steps(Elevation_Stepper_Motor.full_large_cog_rev)
Azimuth_Stepper_Motor.move_steps(0-Azimuth_Stepper_Motor.full_large_cog_rev)
Elevation_Stepper_Motor.move_steps(0-Elevation_Stepper_Motor.full_large_cog_rev)

while track_point < track_point_limit :
    obs.date = datetime.datetime.utcnow()

    iss.compute(obs)
    Azimuth_Stepper_Motor.move_to(convert_radians_to_degrees(float(iss.az)))
    Elevation_Stepper_Motor.move_to(convert_radians_to_degrees(float(iss.alt)))
    time.sleep(1)
    track_point += 1
    print("ISS: Track Point: %s  az: %s ele: %s " % (track_point,convert_radians_to_degrees(float(iss.az)),
                                                            convert_radians_to_degrees(float(iss.alt))))
#lets track the sun
track_point = 0
track_point_limit = 3600
while track_point < track_point_limit :
    obs.date = datetime.datetime.utcnow()
    sun.compute(obs)
    Azimuth_Stepper_Motor.move_to(convert_radians_to_degrees(float(sun.az)))
    Elevation_Stepper_Motor.move_to(convert_radians_to_degrees(float(sun.alt)))
    time.sleep(1)
    track_point += 1
    print("Sun Track Point: %s  az: %s ele: %s " % (track_point,convert_radians_to_degrees(float(sun.az)),
                                                convert_radians_to_degrees(float(sun.alt))))
#range = 24
for p in range(24):
    sun.compute(obs)
    moon.compute(obs)
    venus.compute(obs)
    mercury.compute(obs)
    jupiter.compute(obs)
    saturn.compute(obs)
    print ("Time %s " % ephem.localtime(obs.date))
    print ("Time %s " % obs.date)
    print("Sun Alt: %s Sun Azimuth %s" % (sun.alt, sun.az))
    print("Moon Alt %s Moon Azimuth %s" % (moon.alt, moon.az))
    print("Venus Alt %s Venus Azimute %s" % (venus.alt, venus.az))
    print("Mercury Alt %s Mercury Azimute %s" % (mercury.alt, mercury.az))
    print("Jupiter Alt %s Jupiter Azimute %s" % (jupiter.alt, jupiter.az))
    print("Saturn Alt %s Saturn Azimute %s" % (saturn.alt, saturn.az))
    obs.date = datetime.datetime.now()  + datetime.timedelta(hours=p)

planets = [
    ephem.Moon(),
    ephem.Mercury(),
    ephem.Venus(),
    ephem.Mars(),
    ephem.Jupiter(),
    ephem.Saturn()
    ]

min_alt = 10. * math.pi / 180
obs.date = datetime.datetime.utcnow()
for time in range(12):
    for planet in planets:
        planet.compute(obs)
        if planet.alt > min_alt:
            print planet.name, " is visable: Alt %s Azimuth %s Time %s" % (planet.alt, planet.az, ephem.localtime(obs.date))


obs.date = datetime.datetime.utcnow()
moon.compute(obs)
print (moon.alt, moon.az, ephem.localtime(obs.date))
obs.date = datetime.datetime.utcnow()
moon.compute(obs)
moon_setting = obs.next_setting(moon)
moon_location = list()
print ("Moon set at %s " % ephem.localtime(moon_setting))
while obs.date < moon_setting:
        obs.date += ephem.minute
        moon.compute(obs)
        print(moon.alt, moon.az)
        temp_tuple = (moon.alt), (moon.az)
        moon_location.append(temp_tuple)
        plt.plot([(moon.alt)],[(moon.az)])
t = np.arange(0.,90.,0.1)

jupiter_location = list()
obs.date = datetime.datetime.utcnow()
jupiter.compute(obs)
#calling it a night here......was goint to calculate jupitar but it'll be easier with a for loop current_time > planet.rising...
#while obs.date
#plt.axis([0,90,1,360])
#plt.show()
#plt.clf();
#generate_satellite_plot(12,14)

rc('grid', color='#316931', linewidth=1, linestyle='-')
rc('xtick', labelsize=15)
rc('ytick', labelsize=15)

# force square figure and square axes looks better for polar, IMO
width, height = rcParams['figure.figsize']
size = min(width, height)
# make a square figure
fig = figure(figsize=(size, size))

ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True, facecolor='#d5de9c')
ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)
legend = ax.legend(loc='upper center', shadow=True)

print(len(moon_location))
for planet in planets:
    planet_location = predict_planet_orbit(obs, planet)
    planet_color = get_planet_colour(planet)

    for (E, Az) in planet_location:
        ra,dec = obs.radec_of(Az,E)
        if (np.degrees(E) > min_alt):
            ax.scatter(Az,90 - np.degrees(E),color=planet_color)


ax.set_yticks(range(0, 90 + 10, 10))  # Define the yticks
yLabel = ['90', '', '', '60', '', '', '30', '', '', '']
ax.set_yticklabels(yLabel)
grid(True)

plt.show(block=True)
fig.show(block=True)
fig.hold()




GPIO.cleanup()