"""
Script for locating the Earth's shadow for a given epoch
"""

import argparse as ap
import numpy as np
from datetime import datetime
import ephem
from skyfield.api import load
from astropy import units as u
from astropy.coordinates import SkyCoord, Longitude, Latitude
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

RA_CORRECTION = Longitude(12., u.hourangle)
RAD_TO_DEG = 180. / np.pi

SITE_LATITUDE = 28.7603135
SITE_LONGITUDE = -17.8796168
SITE_ELEVATION = 2387

R_SUN = 695508000.
R_EARTH = 6371000.
R_GEO = 36000000.
CONE_ANGLE = 0.018867

def argParse():
    """
    Argument parser settings
    """
    
    parser = ap.ArgumentParser()
    
    parser.add_argument('epoch',
                        help='format "yyyy/MM/dd hh:mm:ss" or "now"',
                        type=str)
    
    parser.add_argument('ra',
                        help='Right ascension of target [deg]',
                        type=float)
    
    parser.add_argument('dec',
                        help='Declination og target [deg]',
                        type=float)
    
    return parser.parse_args()

def getSunPosition(epoch):
    """
    Compute position of Sun relative to observer
    """
    
    obs = ephem.Observer()
    
    obs.lat = SITE_LATITUDE
    obs.lon = SITE_LONGITUDE
    obs.elevation = SITE_ELEVATION
    obs.pressure = 0
    
    obs.date = epoch
    
    sun = ephem.Sun(obs)
    print('Sun is at ', SkyCoord(ra=sun.ra*u.rad, dec=sun.dec*u.rad, unit=u.deg))
    
    return SkyCoord(ra=sun.ra*u.rad, dec=sun.dec*u.rad, unit=u.deg)

def getAntiSolarPoint(coord):
    """
    Compute antisolar point given solar coordinates
    """
    
    ra_asp = Longitude(coord.ra - RA_CORRECTION, unit=u.deg)
    dec_asp = Latitude(-coord.dec, unit=u.deg)
    
    return SkyCoord(ra=ra_asp, dec=dec_asp)

def earthSunDistance(epoch):
    """
    Compute the Earth-Sun distance at the desired epoch
    """
    
    ts = load.timescale()
    epoch = datetime.strptime(epoch, '%Y/%m/%d %H:%M:%S')
    t = ts.tt(epoch.year,epoch.month,epoch.day,epoch.hour,epoch.minute,epoch.second)
    
    planets = load('de421.bsp')
    
    ephemeris = planets['earth'].at(t).observe(planets['sun'])
    
    _, _, distance = ephemeris.radec()
    
    return distance.m

def shadowRadius(altitude):
    """
    Compute the radius of the inner shadow
    """
    
    base_diam = 2.*R_EARTH
    wing_diam = 2.*R_GEO*np.tan(0.5*CONE_ANGLE)
    
    r_shadow = (base_diam - wing_diam) / 2.
    
    ang_radius = 2.*np.arctan(r_shadow / (2.*R_GEO))
    
    return ang_radius * RAD_TO_DEG

def penumbraRadius(d_sun):
    """
    Compute the radius of the outer (penumbral) shadow
    """
    
    theta = np.arctan((R_EARTH + R_SUN) / d_sun)
    
    d_sun_x = R_SUN * np.tan(np.pi / 2 - theta)
    
    r_penumbra = (R_EARTH + R_SUN) / d_sun * (d_sun - d_sun_x + R_EARTH + R_GEO)
    
    ang_radius = 2.*np.arctan(r_penumbra / (2.*R_GEO))
    
    return ang_radius * RAD_TO_DEG

def main(epoch, ra, dec, altitude=R_GEO):
    """
    Plot the Earth's inner and penumbral shadows
    """
    
    sun_coord = getSunPosition(epoch)
    asp = getAntiSolarPoint(sun_coord)
    
    d_s = earthSunDistance(epoch)
    r_s = shadowRadius(altitude)
    r_p = penumbraRadius(d_s)
    
    # plot the observer's target position
    fig, ax = plt.subplots()
    plt.plot(sun_coord.ra.deg, sun_coord.dec.deg,'rx')
    plt.plot(ra, dec, 'kx')
    
    plt.xlim(0, 360)
    plt.ylim(-90,90)
    
    c_i = Circle(xy=(asp.ra.deg,asp.dec.deg), radius=r_s, 
               facecolor="black", edgecolor='black', alpha=0.8)
    c_i.set_facecolor('black')
    c_i.set_edgecolor('black')
    ax.add_artist(c_i)
    
    c_p = Circle(xy=(asp.ra.deg,asp.dec.deg), radius=r_p, 
               facecolor="black", edgecolor='black', alpha=0.6)
    c_p.set_facecolor('black')
    c_p.set_edgecolor('black')
    ax.add_artist(c_p)
    
    plt.grid(True)
    
    plt.title('Epoch: ' + epoch)
    
    plt.xlabel('ra [deg]')
    plt.ylabel('dec [deg]')
    
    plt.show()
    plt.close(fig)

if __name__ == "__main__":
    
    args = argParse()
    
    main(args.epoch, args.ra, args.dec)
