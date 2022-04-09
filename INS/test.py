import numpy as np

# Constants: 

# 1. Radius of Earth (meters)
Re_m = 6,371,000

# 2. Rotation rate of earth (rad/sec) 
omega_e_rps = 7.2921159e-05


def vec2skew3d(vec):

    res = np.zeros((3,3))

    res[0, 0] = 0.0;      res[0,1] = -vec[2];   res[0,2] = vec[1]

    res[1, 0] = vec[2];   res[1,1] = 0.0;       res[1,2] = -vec[0]

    res[2, 0] = -vec[1];  res[2,1] = vec[0];    res[2,2] = vec[1]
 
    return res

def deg2Rad(angle_deg):

    return angle_deg * (np.pi / 180.0)

def calcTransportRateVec(v_n, lat_rad, long_rad, height_m):

    lat_rate_rps  = v_n[1] / (Re_m + height_m)
    long_rate_rps = v_n[0] / (Re_m + height_m)

    omega_e2n_n = np.zeros((3,1))

    omega_e2n_n[0] = -lat_rate_rps
    omega_e2n_n[1] =  long_rate_rps * np.cos(lat_rad)
    omega_e2n_n[2] =  0.0

    return omega_e2n_n


def latLongToDCM(lat_rad, long_rad):

    C_e2n = np.zeros((3,3))

    C_e2n[0, 0] = np.cos(long_rad);                         C_e2n[0, 1] = 0.0;                  C_e2n[0, 2] = -np.sin(long_rad) 

    C_e2n[1, 0] = -np.sin(long_rad) * np.sin(lat_rad);      C_e2n[1, 1] = np.cos(lat_rad);      C_e2n[1, 2] = -np.cos(long_rad) * np.sin(lat_rad) 

    C_e2n[2, 0] = np.sin(long_rad) * np.cos(lat_rad);       C_e2n[2, 1] = np.sin(lat_rad);      C_e2n[2, 2] = np.cos(long_rad) * np.cos(lat_rad)

    return C_e2n

def DCMToLatLong(C_e2n):

    lat  = np.arcsin(C_e2n[2, 1])
    long = np.arctan( C_e2n[2,0] / ( C_e2n[0,0] * C_e2n[1,1] - C_e2n[0,1] * C_e2n[1,0]))

    return lat, long

## Setup time vector and simulated measurements: 
fs = 100
t_0 = 0.0
t_f = 10.0
n_steps = (t_f - t_0) * fs
t_vec = np.linspace(t_0, t_f, n_steps)


## Given a simulated set of accel and gyroscope measurements:
accel_meas = np.zeros((3, n_steps))
gyro_meas  = np.zeros((3, n_steps))

accel_meas[0,:] = 1.5 * np.sin(2 * np.pi * 0.05 * t_vec)
accel_meas[2,:] = 0.005 * np.ones_like(t_vec)

gyro_meas[2,:] = 0.0001 * np.ones_like(t_vec)

## Forward propagate the: 
# 1. Navigation to Earth Frame DCM
# 2. Velocity vector in navigation frame
# 3. Body to navigation frame DCM

# 1. Initialize the required state:
lat_init  = deg2Rad(37.230000) 
long_init = deg2Rad(-80.417778)

v_n = np.zeros((3,1))
C_b2n = np.eye(3)

C_e2n = latLongToDCM(lat_init, long_init)
