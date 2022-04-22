import georinex as gr
import numpy as np
from scipy.optimize import root
from scipy.spatial.transform import Rotation as R

def keplerEccentricEquation(x, e, M_k):

    return x - e*np.sin(x) - M_k 

## Constants: 
mu = 3986004.418E+08         # mu = G* M_E (Units - m^3/s^2)
omega_e_rps = 7.2921159e-05  # Units - radians per sec

dat = gr.load('data/BREW00USA_R_20220311200_01H_GN.rnx')

sv_g01 = dat.sel(sv='G01')
print(sv_g01['sqrtA'].values)
print(sv_g01['sqrtA'].coords)

print(sv_g01['M0'].values)
print(sv_g01['M0'].coords)

print(sv_g01)

# Step 2: Determine the time offset from the Issue of Data Ephemeris 
Toe = sv_g01['Toe'].values[1]
t_curr = Toe + 45

t_k = t_curr - Toe

print('T_k = ', t_k)

# Step 3: Determine the mean anomaly
sqrtA = sv_g01['sqrtA'].values[1]
delta_n = sv_g01['DeltaN'].values[1]
M0 = sv_g01['M0'].values[1]

Mk = M0 + (np.sqrt(mu)/ np.power(sqrtA, 3) + delta_n) * t_k

print(' Mean anomaly is: ', Mk)

# Step 4: Determine the true anomaly:
e = sv_g01['Eccentricity'].values[1]

res = root(keplerEccentricEquation, Mk, args=(e, Mk))

Ek = res.x[0]

print(' Eccentric anomaly is', Ek)

numerator = np.sqrt(1 - np.power(e,2)) * np.sin(Ek)
denominator = np.cos(Ek) - e
vk = np.arctan(numerator/denominator)  # True anomaly

print(' True anomaly is ', vk)

# Step 5: Compute the radial distance of the satellite from the focus of the ellipse (Take into account harmonic corrections)
Crc = sv_g01['Crc'].values[1]
Crs = sv_g01['Crs'].values[1]
omega = sv_g01['omega'].values[1]

angle_r = 2* (omega + vk)

rk = sqrtA**2 * (1 -e* np.cos(Ek)) + Crc * np.cos(angle_r) + Crs * np.sin(angle_r)

print(' Radial distance is: ', rk)

# Step 6: Compute the argument of latitude (uk) using the argument of perigee, true anomaly and harmonic corrections
Cuc = sv_g01['Cuc'].values[1]
Cus = sv_g01['Cus'].values[1]
uk = omega + vk + Cuc * np.cos(angle_r) + Cus * np.sin(angle_r)

print(' Argument of latitude is: ', uk)

# Step 7: Compute the inclination of the satellite orbital plane (ik) at the current time. 
idot = sv_g01['IDOT'].values[1]
i0   = sv_g01['Io'].values[1]

Cic = sv_g01['Cic'].values[1]
Cis = sv_g01['Cis'].values[1]

ik = i0 + idot * t_k + Cic * np.cos(angle_r) + Cis * np.sin(angle_r)

print(' Inclination of orbit plane : ', ik)

# Step 8: Compute the longitude of ascending node (lk) with respect to Greenwich
Omega0 = sv_g01['Omega0'].values[1]
OmegaDot = sv_g01['OmegaDot'].values[1]

lk = Omega0 + (OmegaDot - omega_e_rps) * t_k - omega_e_rps * Toe

print(' Longitude of ascending node is: ', lk)

# Step 9: Compute rotation matrices

Ru = R.from_euler('Z', -uk).as_matrix()

Ri = R.from_euler('X', -ik).as_matrix()

Rl = R.from_euler('Z', -lk).as_matrix()

print(Ru)
print(Ri)
print(Rl)

# Step 10: Compute ECEF coordinates:
r_pf = np.zeros((3,1))
r_pf[0] = rk

r_ecef = Rl @ Ri @ Ru @ r_pf

print(' ECEF coordinates are: ', r_ecef)

print('Radial distance is: ', np.linalg.norm(r_ecef))