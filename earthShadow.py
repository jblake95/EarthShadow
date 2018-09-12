"""
Script for locating the Earth's shadow for a given epoch
"""

import argparse as ap
import numpy as np
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import (
    Angle,
    Latitude,
    Longitude,
    SkyCoord,
    EarthLocation,
    solar_system_ephemeris, 
    get_body
    )
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

RA_CORRECTION = Longitude(12., u.hourangle)
CONE_ANGLE = 0.018867
RAD_TO_DEG = 180. / np.pi
RIGHT_ANGLE = np.pi / 2

SITE_LATITUDE = 28.7603135
SITE_LONGITUDE = -17.8796168
SITE_ELEVATION = 2387

SITE_LOCATION = EarthLocation(lat=SITE_LATITUDE*u.deg,
                              lon=SITE_LONGITUDE*u.deg,
                              height = SITE_ELEVATION*u.m)

R_SUN = 695508000.
R_EARTH = 6371000.
R_GEO = 36000000.

def argParse():
    """
    Argument parser settings
    """
    
    parser = ap.ArgumentParser()
    
    parser.add_argument('utc',
                        help='format "yyyy-MM-ddThh:mm:ss" or "now"',
                        type=str)
    
    parser.add_argument('ra_h',
                        help='Right ascension of target [hr]',
                        type=float)
    
    parser.add_argument('ra_m',
                        help='Right ascension of target [min]',
                        type=float)
    
    parser.add_argument('ra_s',
                        help='Right ascension of target [sec]',
                        type=float)
    
    parser.add_argument('dec_d',
                        help='Declination of target [deg]',
                        type=float)
                        
    parser.add_argument('dec_m',
                        help='Declination of target [deg]',
                        type=float)
    
    parser.add_argument('dec_s',
                        help='Declination of target [hr]',
                        type=float)
    
    return parser.parse_args()

def parseInput(args):
    """
    Read input epoch as an astropy Time object and input angles as 
    astropy Angle objects
    """
    
    epoch = Time(args.utc, format='isot', scale='utc')
    
    ra = Longitude((args.ra_h, args.ra_m, args.ra_s), u.hourangle)
    dec = Latitude((args.dec_d, args.dec_m, args.dec_s), u.deg)
    
    return epoch, ra, dec

def getSunEphemeris(epoch, loc=SITE_LOCATION):
    """
    Compute position of Sun relative to observer [epoch in utc]
    """
    
    sun_ephem = get_body('sun', epoch, loc)
    
    print('Sun is at ', (sun_ephem.ra, sun_ephem.dec))
    
    return sun_ephem

def getAntiSolarPoint(ephemeris):
    """
    Compute antisolar point given solar coordinates
    """
    
    ra_asp = Longitude(ephemeris.ra - RA_CORRECTION, unit=u.hourangle)
    dec_asp = Latitude(-ephemeris.dec, unit=u.deg)
    
    return SkyCoord(ra=ra_asp, dec=dec_asp)

def umbraRadius(altitude):
    """
    Compute the radius of the inner (umbral) shadow
    """
    
    base_diam = 2. * R_EARTH
    wing_diam = 2. * (altitude + R_EARTH) * np.tan(0.5 * CONE_ANGLE)
    
    r_umbra = (base_diam - wing_diam) / 4.
    
    ang_radius = 2. * np.arctan(r_umbra / (2. * (altitude - SITE_ELEVATION)))
    
    return ang_radius * RAD_TO_DEG

def penumbraRadius(d_sun):
    """
    Compute the radius of the outer (penumbral) shadow
    """
    
    theta = np.arctan((R_EARTH + R_SUN) / d_sun)
    d_sun_x = R_SUN * np.tan(RIGHT_ANGLE - theta)
    
    r_penumbra = (R_EARTH + R_SUN) / d_sun * (d_sun - d_sun_x + R_EARTH + R_GEO)
    
    ang_radius = 2. * np.arctan(r_penumbra / (2. * (R_GEO - SITE_ELEVATION)))
    
    return ang_radius * RAD_TO_DEG

def sexagesimal(angle_str):
    """
    Convert an angle string in degrees to a sexagesimal longitude
    """
    
    angle = Angle(angle_str, u.deg)
    sexagesimal_str = angle.to_string(u.hourangle)
    
    return sexagesimal_str

def main(args, altitude=R_GEO):
    """
    Plot Earth's umbral and penumbral shadows against target position
    """
    
    epoch, target_ra, target_dec = parseInput(args)
    
    sun_ephem = getSunEphemeris(epoch)
    asp = getAntiSolarPoint(sun_ephem)
    
    d_s = sun_ephem.distance
    r_u = umbraRadius(altitude)
    r_p = penumbraRadius(d_s.value) 
    
    # plot the observer's target position
    fig, ax = plt.subplots()
    
    plt.plot(sun_ephem.ra.deg, sun_ephem.dec.deg,'rx')
    plt.plot(target_ra, target_dec, 'gx')
    """
    # present right ascension axis in units of hour angle
    fig.canvas.draw()
    
    new_xlabels = []
    old_xlabels = [item.get_text() for item in ax.get_xticklabels()]
    for l, label in enumerate(old_xlabels):
        print(l, label)
        new_xlabels.append(sexagesimal(old_xlabels[l]))
    
    ax.set_xticklabels(new_xlabels) 
    """
    # add circles to represent umbral and penumbral shadows
    c_u = Circle(xy=(asp.ra.deg,asp.dec.deg), radius=r_u, 
               facecolor="black", edgecolor='black', alpha=0.8)
    c_u.set_facecolor('black')
    c_u.set_edgecolor('black')
    ax.add_artist(c_u)
    
    c_p = Circle(xy=(asp.ra.deg,asp.dec.deg), radius=r_p, 
               facecolor="black", edgecolor='black', alpha=0.6)
    c_p.set_facecolor('black')
    c_p.set_edgecolor('black')
    ax.add_artist(c_p)
    
    # final touches    
    plt.title('Epoch: ' + str(epoch))
    
    plt.xlabel('ra [deg]')
    plt.ylabel('dec [deg]')
    
    plt.xlim(0, 360)
    plt.ylim(-90,90)
    
    plt.grid(True)
    
    plt.show()
    plt.close(fig)

if __name__ == "__main__":
    
    args = argParse()
    
    main(args)
