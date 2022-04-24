import numpy as np

## Coordinate frame conversions: 

def deg2Rad(angle_deg):

    return angle_deg * (np.pi / 180.0)

def radToDeg(angle_rad):

    return angle_rad * (180 / np.pi)

# 1. LLH to ECEF:
# Params obtained from: https://en.wikipedia.org/wiki/World_Geodetic_System#A_new_World_Geodetic_System:_WGS_84 

def LLH2ECEF(lat, long, height):

    """
    Params: 
    1. R_N - Radius of curvature in the prime vertical (m)
    2. a - Semi-major axis of Earth (a.k.a Equatorial radius )according to WGS-84 model
    3. e - Eccentricity of the earth according to WGS-84 model 
    """ 

    a = 6378137  # meters
    e =  0.006694379990141316
    e_sq = e**2
    
    sin_lat = np.sin(lat)
    R_N = a/ np.sqrt(1 - e_sq * sin_lat**2)

    res = np.zeros((3,1))

    res[0] = (R_N + height) * np.cos(lat) * np.cos(long)
    res[1] = (R_N + height) * np.cos(lat) * np.sin(long)
    res[2] = ( (1 - e_sq) *R_N + height) * np.sin(lat)

    return res


# 2. ECEF to LLH: 

def ECEF2LLH(pos_ecef):

    x = pos_ecef[0]
    y = pos_ecef[1]
    z = pos_ecef[2]

    a = 6378137  # meters
    e =  0.006694379990141316
    e_sq = e**2

    long = np.arctan(y/x)

    p = np.sqrt(x**2 + y**2)

    lat_init = np.arctan( (z/p) /(1 - e_sq))
    h = 0
    lat_prev = 10000000000
    lat = lat_init

    itr_count = 1

    print('Outisde loop. Lat diff = ' , abs(lat- lat_prev) )
    
    while ( abs(lat - lat_prev) > 0.0000001):

        print('Itr number ', itr_count)

        lat_prev = lat

        sin_lat = np.sin(lat)
        R_N = a / np.sqrt(1 - e_sq * (sin_lat**2))

        h = p/np.cos(lat) - R_N

        denom = 1 - (R_N/(R_N + h)) * e_sq

        lat = np.arctan((z/p) / denom)

        itr_count = itr_count + 1

    return (lat, long, h)



## Test:

chennai_lat = deg2Rad(13.067439)
chennai_long = deg2Rad(80.237617)
chennai_elev = 6.7   # meters


ecef_chennai = LLH2ECEF(chennai_lat, chennai_long, chennai_elev)

print(' ECEF coords are: ', ecef_chennai)

(res_lat, res_long, res_elev) = ECEF2LLH(ecef_chennai)

print(' Latitude is: ', radToDeg(res_lat))
print(' Longitude is: ', radToDeg(res_long))
print(' Elevation is: ', res_elev)


