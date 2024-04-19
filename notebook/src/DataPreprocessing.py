### all the data pre-processing code for the FGKM project
### contributions from ZI, Brandon Sally, Sid Chaini, Bruno Dominguez, and Karlo Mrakovcic 

import numpy as np
import astropy.units as u
import astropy.coordinates as ac
_ = ac.galactocentric_frame_defaults.set('v4.0') #set the default Astropy Galactocentric frame parameters to the values adopted in Astropy v4.0
from astropy.coordinates import (CartesianRepresentation, CartesianDifferential, Galactic)
import astropy.units as u

import healpy as hp
import healpy_utils as hpu 


def preprocessData(df):

    # solar motion parameters 
    Rsol, v_LSR, vX_sun, vY_sun, vZ_sun = getSolarMotion()

    ### first trivial transformations
    df['dist'] = df['r_med_photogeo']
    df['Dkpc'] = df['dist']/1000.0
    df['DM'] = 5*np.log10(df['r_med_photogeo']/10.0)
    # SDSS colors
    df['ug'] = df['psfmag_u'] - df['psfmag_g']
    df['gr'] = df['psfmag_g'] - df['psfmag_r']
    df['ri'] = df['psfmag_r'] - df['psfmag_i']
    df['iz'] = df['psfmag_i'] - df['psfmag_z']
    df['gi'] = df['gr'] + df['ri']

    
    ### For transforming heliocentric frame to galactocentric frame
    ra = np.array(df['ra'])
    dec = np.array(df['dec'])
    dist = np.array(df['Dkpc']) 
    pmra = np.array(df['pmra'])
    pmdec = np.array(df['pmdec'])
    rv = np.array(df['radial_velocity'])
    c_ecu = ac.ICRS(ra=ra*u.degree, dec=dec*u.degree, distance=dist*u.kpc, pm_ra_cosdec=pmra*u.mas/u.yr, pm_dec=pmdec*u.mas/u.yr, radial_velocity=rv*u.km/u.s)
    c_gal = c_ecu.transform_to(ac.Galactic())

    ### Cartesian (right-handed) galactocentric coordinates (in kpc): 
    # X: galactic center to Earth
    # Y: perpendicular to X and Z
    # Z: Earth to north galactic pole
    df['X'] = Rsol - dist*np.cos(np.radians(df['l']))*np.cos(np.radians(df['b']))
    df['Y'] = -dist*np.sin(np.radians(df['l']))*np.cos(np.radians(df['b']))
    df['Z'] = dist*np.sin(np.radians(df['b']))
    df['R'] = np.sqrt(df['X']**2 + df['Y']**2)
    df['phi'] = np.arctan(df['Y']/df['X']) # N.B. in radians
    df['pmL'] = c_gal.pm_l_cosb.value
    df['pmB'] = c_gal.pm_b.value
   
    ### change velocity frame from equatorial to galactic
    A_v = 4.74  # mu/(mas yr^-1) D/kpc km/s
    v_l = A_v*dist*c_gal.pm_l_cosb.value
    v_b = A_v*dist*c_gal.pm_b.value
    v_rad = c_gal.radial_velocity.value
    df['vL'] = v_l 
    df['vB'] = v_b 
    df['vradGal'] = c_gal.radial_velocity.value  

    # Cartesian galactocentric velocities
    df['vX_obs'] = -v_rad*np.cos(np.radians(df['l']))*np.cos(np.radians(df['b'])) \
                   + v_b * np.cos(np.radians(df['l'])) * np.sin(np.radians(df['b'])) \
                   + v_l*np.sin(np.radians(df['l']))
    df['vY_obs'] = -v_rad*np.sin(np.radians(df['l']))*np.cos(np.radians(df['b'])) \
                   + v_b * np.sin(np.radians(df['l'])) * np.sin(np.radians(df['b'])) \
                   - v_l*np.cos(np.radians(df['l']))
    df['vZ_obs'] = v_rad * np.sin(np.radians(df['b'])) + v_b * np.cos(np.radians(df['b']))

    # correct galactocentric velocities for solar motion
    df['v_X'] = df['vX_obs'] + vX_sun
    df['v_Y'] = df['vY_obs'] - v_LSR + vY_sun
    df['v_Z'] = df['vZ_obs'] + vZ_sun

    # velocity in cylindrical coordinates
    df['v_R'] = df['v_X']*df['X']/df['R'] + df['v_Y']*df['Y']/df['R']
    df['v_phi'] = -df['v_X']*df['Y']/df['R'] + df['v_Y']*df['X']/df['R']
    
    # N.B.
    # Obs: retrograde rotation is indicated by vφ > 0.
    # stars with vR > 0 move away from the Galactic center, and
    # stars with vZ > 0 move toward the North Galactic Pole.

    print('done!')
    return 



def preprocessDataGiants(df):

    # solar motion parameters 
    Rsol, v_LSR, vX_sun, vY_sun, vZ_sun = getSolarMotion()

    ### first trivial transformations
    df['piSNR'] = df['parallax']/df['parallax_error']
    df['Dkpc'] = 1.0/df['parallax']
    df['dist'] = 1000*df['Dkpc']
    df['DM'] = 5*np.log10(df['dist']/10.0)

    # Gaia photometry etc 
    df['BpRp'] = df['phot_bp_mean_mag'] - df['phot_rp_mean_mag']
    df['G'] = df['phot_g_mean_mag']
    df['FeH'] = df['mh_xgboost']
    
    ### For transforming heliocentric frame to galactocentric frame
    ra = np.array(df['ra'])
    dec = np.array(df['dec'])
    dist = np.array(df['Dkpc']) 
    pmra = np.array(df['pmra'])
    pmdec = np.array(df['pmdec'])
    rv = np.array(df['radial_velocity'])
    c_ecu = ac.ICRS(ra=ra*u.degree, dec=dec*u.degree, distance=dist*u.kpc, pm_ra_cosdec=pmra*u.mas/u.yr, pm_dec=pmdec*u.mas/u.yr, radial_velocity=rv*u.km/u.s)
    c_gal = c_ecu.transform_to(ac.Galactic())

    ### Cartesian (right-handed) galactocentric coordinates (in kpc): 
    # X: galactic center to Earth
    # Y: perpendicular to X and Z
    # Z: Earth to north galactic pole
    df['l'] = c_gal.l.value
    df['b'] = c_gal.b.value
    df['X'] = Rsol - dist*np.cos(np.radians(df['l']))*np.cos(np.radians(df['b']))
    df['Y'] = -dist*np.sin(np.radians(df['l']))*np.cos(np.radians(df['b']))
    df['Z'] = dist*np.sin(np.radians(df['b']))
    df['R'] = np.sqrt(df['X']**2 + df['Y']**2)
    df['phi'] = np.arctan(df['Y']/df['X']) # N.B. in radians
    df['pmL'] = c_gal.pm_l_cosb.value
    df['pmB'] = c_gal.pm_b.value
   
    ### change velocity frame from equatorial to galactic
    A_v = 4.74  # mu/(mas yr^-1) D/kpc km/s
    v_l = A_v*dist*c_gal.pm_l_cosb.value
    v_b = A_v*dist*c_gal.pm_b.value
    v_rad = c_gal.radial_velocity.value
    df['vL'] = v_l 
    df['vB'] = v_b 
    df['vradGal'] = c_gal.radial_velocity.value  

    # Cartesian galactocentric velocities
    df['vX_obs'] = -v_rad*np.cos(np.radians(df['l']))*np.cos(np.radians(df['b'])) \
                   + v_b * np.cos(np.radians(df['l'])) * np.sin(np.radians(df['b'])) \
                   + v_l*np.sin(np.radians(df['l']))
    df['vY_obs'] = -v_rad*np.sin(np.radians(df['l']))*np.cos(np.radians(df['b'])) \
                   + v_b * np.sin(np.radians(df['l'])) * np.sin(np.radians(df['b'])) \
                   - v_l*np.cos(np.radians(df['l']))
    df['vZ_obs'] = v_rad * np.sin(np.radians(df['b'])) + v_b * np.cos(np.radians(df['b']))

    # correct galactocentric velocities for solar motion
    df['v_X'] = df['vX_obs'] + vX_sun
    df['v_Y'] = df['vY_obs'] - v_LSR + vY_sun
    df['v_Z'] = df['vZ_obs'] + vZ_sun

    # velocity in cylindrical coordinates
    df['v_R'] = df['v_X']*df['X']/df['R'] + df['v_Y']*df['Y']/df['R']
    df['v_phi'] = -df['v_X']*df['Y']/df['R'] + df['v_Y']*df['X']/df['R']

    # velocity in spherical coordinates (Bond+2010, eq. 25)
    df['Rgc'] = np.sqrt(df['R']**2+df['Z']**2)
    df['v_r'] = df['v_R']*df['R']/df['Rgc'] + df['v_Z']*df['Z']/df['Rgc']
    df['v_th'] = df['v_R']*df['Z']/df['Rgc'] - df['v_Z']*df['R']/df['Rgc']

    
    # N.B.
    # Obs: retrograde rotation is indicated by vφ > 0.
    # stars with vR > 0 move away from the Galactic center, and
    # stars with vZ > 0 move toward the North Galactic Pole.

    print('done with giants!')
    return 


### Solar motion from Dehnen & Binney (1998) 
def getSolarMotion():
    ### numerical assumptions about Solar motion
    # standard vLSR
    vLSR = 220.0   
    # peculiar solar motion 
    vXo = -10.0 
    vYo =  -5.3    
    vZo =   7.2  
    # distance from the galactic center, from Eisenhauer+2003 (7.94 +- 0.42 kpc)  
    Rsun = 8.0  # kpc 
    return Rsun, vLSR, vXo, vYo, vZo


#### healpy projection tools

# Functions for creating the healbins for the map pixels, but now using rms as the applied function.
def sigG(x):
    '''
    Calculates the interquartile range and then adjusts it with the 0.741 parameter. 
    This creates a more robust version of the standard deviation.

    Arguments:
        x(array) = data we want to calculate the interquartile range with
    '''
    return 0.741*(np.percentile(x,75)-np.percentile(x,25))


def getHealpyBins(df, nside=64):
    
    STDfunc=sigG  # robust standard deviation 
    
    # Creating the bins to be used for mapping.
    pmLongM = hpu.healbin(np.array(df['l']), np.array(df['b']), np.array(df['pmL']), nside=nside)
    pmLatM = hpu.healbin(np.array(df['l']), np.array(df['b']), np.array(df['pmB']), nside=nside)
    radvelM = hpu.healbin(np.array(df['l']), np.array(df['b']), np.array(df['radial_velocity']), nside=nside)

    pmLongS = hpu.healbin(np.array(df['l']), np.array(df['b']), np.array(df['pmL']), nside=nside, reduce_func=STDfunc)
    pmLatS = hpu.healbin(np.array(df['l']), np.array(df['b']), np.array(df['pmB']), nside=nside, reduce_func=STDfunc)
    radvelS = hpu.healbin(np.array(df['l']), np.array(df['b']), np.array(df['radial_velocity']), nside=nside, reduce_func=STDfunc)

    return hp.ma(pmLongM), hp.ma(pmLongS), hp.ma(pmLatM), hp.ma(pmLatS), hp.ma(radvelM), hp.ma(radvelS)


# similar to getHealpyBins, except that SDSS-based model values from Bond et al. (2010) are
# subtracted from observed proper motions and radial velocity for means, and further normalized
# by expected dispersions to compute standard deviation  
def getHealpyBinsModels(df, nside=64):

    Z = df['Z'] # Z in kpc
    
    ### Bond+(2010) models and equations
    # eq. 17: mean vPhi for disk stars
    vPhiModel = -205 + 19.2*np.abs(Z)**1.25

    # eq. 18: vPhi dispersion for disk stars
    vPhiModelDisp = 30 + 3.0*Z**2.0

    # eq. 21: vR dispersion for disk stars
    vRModelDisp = 40 + 5.0*np.abs(Z)**1.5
    # sigmaM = 18 + 4.0*(Zg/1000)**1.5

    # eq. 24: vZ dispersion for disk stars
    vZModelDisp = 25 + 4.0*np.abs(Z)**1.5
    # vZModelDisp = 18 + 4.0*np.abs(Z)**1.5

    ### now get model-based proper motion components and radVel 
    pmLmodel, pmBmodel, radVelModel = vModel2obs(df['X'], df['Y'], df['Z'], 0*vPhiModel, vPhiModel, 0*vPhiModel)
    pmLdispP, pmBdispP, radVelDispP = vModel2obs(df['X'], df['Y'], df['Z'], vRModelDisp, vPhiModel+vPhiModelDisp, vZModelDisp)
    pmLdisp = np.std(pmLdispP - pmLmodel)
    pmBdisp = np.std(pmBdispP - pmBmodel)
    radVelDisp = np.std(radVelDispP - radVelModel)
    
    ## healpy binning... 
    STDfunc=sigG  # robust standard deviation 

    # Creating the bins to be used for mapping.
    pmLongM = hpu.healbin(np.array(df['l']), np.array(df['b']), np.array(df['pmL']-pmLmodel), nside=nside)
    pmLatM = hpu.healbin(np.array(df['l']), np.array(df['b']), np.array(df['pmB']-pmBmodel), nside=nside)
    radvelM = hpu.healbin(np.array(df['l']), np.array(df['b']), np.array(df['radial_velocity']-radVelModel), nside=nside)

    pmLongS = hpu.healbin(np.array(df['l']), np.array(df['b']), np.array((df['pmL']-pmLmodel)/pmLdisp), nside=nside, reduce_func=STDfunc)
    pmLatS = hpu.healbin(np.array(df['l']), np.array(df['b']), np.array((df['pmB']-pmBmodel)/pmBdisp), nside=nside, reduce_func=STDfunc)
    radvelS = hpu.healbin(np.array(df['l']), np.array(df['b']), np.array((df['radial_velocity']-radVelModel)/radVelDisp), nside=nside, reduce_func=STDfunc)

    return hp.ma(pmLongM), hp.ma(pmLongS), hp.ma(pmLatM), hp.ma(pmLatS), hp.ma(radvelM), hp.ma(radvelS)
   


# translated old SM code from the Bond et al. work...
# given (X, Y, Z) and (vR, vPhi, vZ), return (pmLong, pmLat, radVel) 
def vModel2obs(X, Y, Z, vR, vPhi, vZ):

    ### numerical assumptions about Solar motion
    Rsun, vLSR, vXo, vYo, vZo = getSolarMotion()

    ## first transform cylindrical velocities to cartesian velocities
    R = np.sqrt(X**2 + Y**2)
    vX = (vR*X - vPhi*Y)/R
    vY = (vPhi*X + vR*Y)/R

    ## add correction for LSR and solar peculiar velocity
    vXhelio = vX - vXo
    vYhelio = vY - vYo + vLSR
    vZhelio = vZ - vZo

    ## and now translate v(XYZ)helio to proper motions and radial velocity
    # Vxyz2obs vXhelio vYhelio vZhelio X Y Z vradial vL vB
    vL, vB, radVel = Vxyz2obs(X, Y, Z, vXhelio, vYhelio, vZhelio, Rsun=Rsun) 

    # distance in kpc
    Dkpc = np.sqrt((X-Rsun)**2 + Y**2 + Z**2)  
    # proper motions
    pmLong = vL / 4.74 / Dkpc
    pmLat  = vB / 4.74 / Dkpc

    return pmLong, pmLat, radVel 

 
# given Cartesian heliocentric coordinates and (traditional UVW) velocities, return tangential
# and radial velocity components 
def Vxyz2obs(X, Y, Z, vX, vY, vZ, Rsun=8.0):

    Xloc = X - Rsun
    zerovec = X*0.
    onevec = X*0.+1.

    runitx, runity, runitz = makeunitvec2(Xloc, Y, Z)
    radVel = dotproduct2(vX, vY, vZ, runitx, runity, runitz) 
    tempx, tempy, tempz = crossproduct2(runitx, runity, runitz, vX, vY, vZ)
    pmx, pmy, pmz = crossproduct2(runitx, runity, runitz, tempx, tempy, tempz)
    lx, ly, lz = crossproduct2(runitx, runity, runitz, zerovec, zerovec, onevec)
    lux, luy, luz = makeunitvec2(lx, ly, lz)
    vLong = dotproduct2(lux, luy, luz, pmx, pmy, pmz)
    bx, by, bz = crossproduct2(runitx, runity, runitz, lux, luy, luz)
    bux, buy, buz = makeunitvec2(bx, by, bz)
    vLat = dotproduct2(bux, buy, buz, pmx, pmy, pmz)
    
    return vLong, vLat, radVel


    
###########################################################
### vector tools from Nick Bond (translated from SM)
##	makeunitvec2
##	dotproduct2 
## 	crossproduct2 
def vectormag2(v1, v2, v3):                                                                    
    # returns the magnitude of a list of 3D vectors
    return np.sqrt(v1**2 + v2**2 + v3**2)

def dotproduct2(v1, v2, v3, w1, w2, w3):
	# returns the dot product of two 3D vectors
    return v1*w1+v2*w2+v3*w3

def crossproduct2(v1, v2, v3, w1, w2, w3):
    # returns the cross product of two 3D vectors
    p1 = v2*w3-v3*w2
    p2 = v3*w1-v1*w3
    p3 = v1*w2-v2*w1
    return p1, p2, p3

def makeunitvec2(v1, v2, v3): 
	# returns a unit 3D vector in the direction of the inputted vector
    magvec = vectormag2(v1, v2, v3)
    return v1/magvec, v2/magvec, v3/magvec





