### Utility functions for GNSS data:

import georinex as gr
import numpy as np
from scipy.optimize import root
from scipy.spatial.transform import Rotation as R

## Constants: 
mu = 3986004.418E+08         # mu = G* M_E (Units - m^3/s^2)
omega_e_rps = 7.2921159e-05  # Units - radians per sec
WEEKSEC = 604800
C = 3E+08

def keplerEccentricEquation(x, e, M_k):

    return x - e*np.sin(x) - M_k 

def getValidIndex(toe_array):

    for idx in range(0, len(toe_array)):

        if not np.isnan(toe_array[idx]):
            return idx

    return -1

def getSatelliteCoordsFromRinex(t, svid, rinex_filename):

    dat = gr.load(rinex_filename)

    sv_g01 = dat.sel(sv=f'G{svid:02d}')

    validIdx = getValidIndex(sv_g01['Toe'].values)

    print(' Valid index is: ', validIdx)

    if validIdx < 0:
        return np.zeros((3, 1))

    # Step 2: Determine the time offset from the Issue of Data Ephemeris 
    Toe = sv_g01['Toe'].values[validIdx]

    t_k = t*1E-09 - Toe

    print('T_k = ', t_k)

    # Step 3: Determine the mean anomaly
    sqrtA = sv_g01['sqrtA'].values[validIdx]
    delta_n = sv_g01['DeltaN'].values[validIdx]
    M0 = sv_g01['M0'].values[validIdx]

    Mk = M0 + (np.sqrt(mu)/ np.power(sqrtA, 3) + delta_n) * t_k

    # print(' Mean anomaly is: ', Mk)

    # Step 4: Determine the true anomaly:
    e = sv_g01['Eccentricity'].values[validIdx]

    res = root(keplerEccentricEquation, Mk, args=(e, Mk))

    Ek = res.x[0]

    # print(' Eccentric anomaly is', Ek)

    numerator = np.sqrt(1 - np.power(e,2)) * np.sin(Ek)
    denominator = np.cos(Ek) - e
    vk = np.arctan(numerator/denominator)  # True anomaly

    # print(' True anomaly is ', vk)

    # Step 5: Compute the radial distance of the satellite from the focus of the ellipse (Take into account harmonic corrections)
    Crc = sv_g01['Crc'].values[validIdx]
    Crs = sv_g01['Crs'].values[validIdx]
    omega = sv_g01['omega'].values[validIdx]

    angle_r = 2* (omega + vk)

    rk = sqrtA**2 * (1 -e* np.cos(Ek)) + Crc * np.cos(angle_r) + Crs * np.sin(angle_r)

    # print(' Radial distance is: ', rk)

    # Step 6: Compute the argument of latitude (uk) using the argument of perigee, true anomaly and harmonic corrections
    Cuc = sv_g01['Cuc'].values[validIdx]
    Cus = sv_g01['Cus'].values[validIdx]
    uk = omega + vk + Cuc * np.cos(angle_r) + Cus * np.sin(angle_r)

    # print(' Argument of latitude is: ', uk)

    # Step 7: Compute the inclination of the satellite orbital plane (ik) at the current time. 
    idot = sv_g01['IDOT'].values[validIdx]
    i0   = sv_g01['Io'].values[validIdx]

    Cic = sv_g01['Cic'].values[validIdx]
    Cis = sv_g01['Cis'].values[validIdx]

    ik = i0 + idot * t_k + Cic * np.cos(angle_r) + Cis * np.sin(angle_r)

    # print(' Inclination of orbit plane : ', ik)

    # Step 8: Compute the longitude of ascending node (lk) with respect to Greenwich
    Omega0 = sv_g01['Omega0'].values[validIdx]
    OmegaDot = sv_g01['OmegaDot'].values[validIdx]

    lk = Omega0 + (OmegaDot - omega_e_rps) * t_k - omega_e_rps * Toe

    # print(' Longitude of ascending node is: ', lk)

    # Step 9: Compute rotation matrices

    Ru = R.from_euler('Z', -uk).as_matrix()

    Ri = R.from_euler('X', -ik).as_matrix()

    Rl = R.from_euler('Z', -lk).as_matrix()

    # print(Ru)
    # print(Ri)
    # print(Rl)

    # Step 10: Compute ECEF coordinates:
    r_pf = np.zeros((3,1))
    r_pf[0] = rk

    r_ecef = Rl @ Ri @ Ru @ r_pf

    print(' ECEF coordinates are: ')
    print(r_ecef)

    # print('Radial distance is: ', np.linalg.norm(r_ecef))

    return r_ecef

def computePseudoRange(raw_gps_measurements, svid, utcMillis):

    timeNanos = raw_gps_measurements[utcMillis][svid]['TimeNanos']
    biasNanos = raw_gps_measurements[utcMillis][svid]['BiasNanos']
    fullBiasNanos = raw_gps_measurements[utcMillis][svid]['FullBiasNanos']

    receivedSvTimeNanos = raw_gps_measurements[utcMillis][svid]['ReceivedSvTimeNanos']
    timeOffsetNanos = raw_gps_measurements[utcMillis][svid]['TimeOffsetNanos']

    t_RX_GPS = timeNanos - (fullBiasNanos + biasNanos)

    t_RX_GPS_mod = t_RX_GPS % (WEEKSEC *1E+09)

    t_TX_GPS = receivedSvTimeNanos + timeOffsetNanos

    delta_T = t_RX_GPS_mod - t_TX_GPS

    # print(' Delta _T = ', delta_T * 1E-06, ' millseconds')

    rho = C * delta_T * 1E-09

    return rho
