import numpy as np

# Constants: 

# 1. Radius of Earth (meters)
Re_m = 6371000

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

def radToDeg(angle_rad):

    return angle_rad * (180 / np.pi)

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

def calcEarthRotRateInNavFrame(omega_e_rps, lat_rad):

    res = np.zeros((3,1))

    res[0] = omega_e_rps * np.cos(lat_rad)
    res[1] = 0.0
    res[2] = - omega_e_rps * np.sin(lat_rad)

    return res

def caclGravityInNavFrame(lat_rad, omega_e_rps):

    g_n = np.zeros((3,))

    g_n[2] = -9.81

    return g_n

## Setup time vector and simulated measurements: 
fs = 100
dt = 1.0/fs
t_0 = 0.0
t_f = 10.0
n_steps = int ((t_f - t_0) * fs)
t_vec = np.linspace(t_0, t_f, n_steps)


## Given a simulated set of accel and gyroscope measurements:
accel_meas = np.zeros((3, n_steps))
gyro_meas  = np.zeros((3, n_steps))

accel_meas[0,:] = 1.5 * np.sin(2 * np.pi * 0.05 * t_vec)
accel_meas[2,:] = 0.005 * np.ones_like(t_vec) - 9.81

gyro_meas[2,:] = 0.0001 * np.ones_like(t_vec)

## Forward propagate the: 
# 1. Velocity vector in navigation frame (v_n)
# 2. Earth Frame to navigation frame DCM (C_e2n)
# 3. Body frame to navigation frame DCM (C_b2n)
# 4. Altitude (height)

# 1. Initialize the required state  (Blacksburg coordinates)
lat_init_rad  = deg2Rad(37.230000) 
long_init_rad = deg2Rad(-80.417778)
height_init_m = 500.0

# Reserve memory and intialize the states at time t = 0
v_n = np.zeros((3,n_steps))
C_b2n = np.zeros((3,3, n_steps))
C_e2n = np.zeros((3, 3, n_steps))
height_m = np.zeros((1, n_steps))

C_b2n[:, :, 0] = np.eye(3)
C_e2n[:, :, 0] = latLongToDCM(lat_init_rad, long_init_rad)
height_m[:,0] = height_init_m 

## Going to use simple Euler integration for now:
## TODO: Put into a for loop. For now just writing down all the equations:
for idx in range(0, n_steps - 1):

    # Extract Lat, long and height from the state: 
    (lat_rad, long_rad) = DCMToLatLong(C_e2n[:, :, idx])

    # Compute helper variables: 
    omega_en_n = calcTransportRateVec(v_n[:,idx], lat_rad, long_rad, height_m[:, idx])
    Omega_ie_n = calcEarthRotRateInNavFrame(omega_e_rps, lat_rad)
    Omega_en_n = vec2skew3d(omega_en_n)

    Omega_ib_b = vec2skew3d(gyro_meas[:, idx])

    Omega_n2b_b = Omega_ib_b - np.transpose(C_b2n[:, :, idx]) @ (Omega_en_n + Omega_ie_n) @ C_b2n[:, :, idx]

    g_n = caclGravityInNavFrame(lat_rad, omega_e_rps)

    # Equation #1:  
    t1 = np.matmul(C_b2n[:, :, idx] , accel_meas[:, idx])
    t2 = np.matmul( (Omega_en_n + 2 * Omega_ie_n) , v_n[:, idx]) 
    v_n_dot = np.matmul(C_b2n[:, :, idx] , accel_meas[:, idx]) - np.matmul( (Omega_en_n + 2 * Omega_ie_n) , v_n[:, idx]) - g_n

    # print(t1.shape)
    # print(t2.shape)
    # print(g_n.shape)

    # Equation #2: 
    C_e2n_dot = -1.0 * Omega_en_n @ C_e2n[:, :, idx]

    # Equation #3:
    C_b2n_dot = C_b2n[:, :, idx] @ Omega_n2b_b

    # Equation #4:  
    h_dot = v_n[2, idx]

    # Forward propagate the variables by simple euler integration:
    # 1. Velocity in navigation frame
    v_n[:, idx+1] = v_n[:, idx] + v_n_dot * dt

    # 2. Earth to Nav frame DCM: 
    C_e2n[:, :, idx + 1] = C_e2n[:, :, idx] + C_e2n_dot * dt

    # 3. Body to Nav frame DCM: 
    C_b2n[:, :, idx + 1] = C_b2n[:, :, idx] + C_b2n_dot * dt

    # 4. Altitude: 
    height_m[:, idx + 1] = height_m[:, idx] + h_dot * dt



print(" [Info] Finished integrations")

(lat_rad, lon_rad) = DCMToLatLong(C_e2n[:, :, -1])

print(" [Info] Final latitude = ", radToDeg(lat_rad), " Final longitude = " , radToDeg(lon_rad))




