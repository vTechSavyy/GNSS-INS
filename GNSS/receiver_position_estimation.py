import numpy as np
import pprint
from gnss_utils import getSatelliteCoordsFromRinex, computePseudoRange
from utils import LLH2ECEF

def deg2Rad(angle_deg):

    return angle_deg * (np.pi / 180.0)

def radToDeg(angle_rad):

    return angle_rad * (180 / np.pi)


# Step 1: Read in the sensor log raw measurements file: 
def readSensorLog(filename): 

    raw_measurements = {}


    with open(filename) as f: 

        lines = f.readlines()

        for l in lines: 

            if l[0:3] == 'Raw':

                entries = l.split(',')

                utcMillisStr = entries[1]

                # print(entries)
                
                # FullBiasNanos,BiasNanos,BiasUncertaintyNanos,DriftNanosPerSecond,DriftUncertaintyNanosPerSecond,HardwareClockDiscontinuityCount,Svid,TimeOffsetNanos,State,ReceivedSvTimeNanos,ReceivedSvTimeUncertaintyNanos,Cn0DbHz,PseudorangeRateMetersPerSecond,PseudorangeRateUncertaintyMetersPerSecond,AccumulatedDeltaRangeState,AccumulatedDeltaRangeMeters,AccumulatedDeltaRangeUncertaintyMeters,CarrierFrequencyHz,CarrierCycles,CarrierPhase,CarrierPhaseUncertainty,MultipathIndicator,SnrInDb,ConstellationType,AgcDb,BasebandCn0DbHz,FullInterSignalBiasNanos,FullInterSignalBiasUncertaintyNanos,SatelliteInterSignalBiasNanos,SatelliteInterSignalBiasUncertaintyNanos,CodeType,ChipsetElapsedRealtimeNanos

                if entries[1] in raw_measurements:
                    # hey there
                    a = 1
                else: 
                    raw_measurements[utcMillisStr] = {}

                SvidStr = entries[11]
                raw_measurements[utcMillisStr][SvidStr] = {}



                raw_measurements[utcMillisStr][SvidStr]['TimeNanos'] = float(entries[2])
                raw_measurements[utcMillisStr][SvidStr]['LeapSecond'] = float(entries[3])
                # raw_measurements[utcMillisStr]['TimeUncertaintyNanos'] = float(entries[4])
                raw_measurements[utcMillisStr][SvidStr]['FullBiasNanos'] = float(entries[5])
                raw_measurements[utcMillisStr][SvidStr]['BiasNanos'] = float(entries[6])
                raw_measurements[utcMillisStr][SvidStr]['BiasUncertaintyNanos'] = float(entries[7])
                raw_measurements[utcMillisStr][SvidStr]['DriftNanosPerSecond'] = float(entries[8])
                raw_measurements[utcMillisStr][SvidStr]['DriftUncertaintyNanosPerSecond'] = float(entries[9])

                raw_measurements[utcMillisStr][SvidStr]['Svid'] = int(entries[11])
                raw_measurements[utcMillisStr][SvidStr]['TimeOffsetNanos'] = float(entries[12])

                raw_measurements[utcMillisStr][SvidStr]['ReceivedSvTimeNanos'] = float(entries[14])
                raw_measurements[utcMillisStr][SvidStr]['ReceivedSvTimeUncertaintyNanos'] = float(entries[15])

                raw_measurements[utcMillisStr][SvidStr]['ConstellationType'] = int(entries[28])


    return raw_measurements


gps_raw_meas = readSensorLog('data/gnss_log_2022_04_18_06_09_51.txt')

print('Raw measurements are: ')
pprint.pprint(gps_raw_meas)
print(' -      --------------------- ')

# Step 2: Read in the nav file using geo-rinex:
r_ecef_03 = getSatelliteCoordsFromRinex(88810344022050, 3, 'data/IISC00IND_R_20221080000_01H_GN.rnx')
# r_ecef_04 = getSatelliteCoordsFromRinex(88810345405825, 4, 'data/IISC00IND_R_20221080000_01H_GN.rnx')
r_ecef_04 = getSatelliteCoordsFromRinex(88810345405825, 1, 'data/IISC00IND_R_20221080000_01H_GN.rnx')
r_ecef_07 = getSatelliteCoordsFromRinex(88810333994443, 7, 'data/IISC00IND_R_20221080000_01H_GN.rnx')
# r_ecef_08 = getSatelliteCoordsFromRinex(88810343306344, 8, 'data/IISC00IND_R_20221080000_01H_GN.rnx')
r_ecef_08 = getSatelliteCoordsFromRinex(88810343306344, 9, 'data/IISC00IND_R_20221080000_01H_GN.rnx')

# Step 3: Compute the pseudo range for each GPS satellite:  
rho_03 = computePseudoRange(gps_raw_meas, str(3), str(1650242392415))
rho_04 = computePseudoRange(gps_raw_meas, str(4), str(1650242392415))
rho_07 = computePseudoRange(gps_raw_meas, str(7), str(1650242392415))
rho_08 = computePseudoRange(gps_raw_meas, str(8), str(1650242392415))

Rho_meas = np.zeros((4,1))
Rho_meas[0] = rho_03
Rho_meas[1] = rho_04
Rho_meas[2] = rho_07
Rho_meas[3] = rho_08

# print(' Rho SV 03 = ', rho_03)
# print(' Rho SV 04 = ', rho_04)
# print(' Rho SV 07 = ', rho_07)
# print(' Rho SV 08 = ', rho_08)

# Step 4: Set an initial guess:
chennai_lat = deg2Rad(13.067439)
chennai_long = deg2Rad(80.237617)
chennai_elev = 6.7   # meters

r_ecef_0 = LLH2ECEF(chennai_lat, chennai_long, chennai_elev)

# Step 5: Solve the linearized system of equations:

norm_03 = np.linalg.norm(r_ecef_03 - r_ecef_0)
norm_04 = np.linalg.norm(r_ecef_04 - r_ecef_0)
norm_07 = np.linalg.norm(r_ecef_07 - r_ecef_0)
norm_08 = np.linalg.norm(r_ecef_08 - r_ecef_0)

print(' Rho ref SV 03 = ', norm_03)
print(' Rho ref SV 04 = ', norm_04)
print(' Rho ref SV 07 = ', norm_07)
print(' Rho ref SV 08 = ', norm_08)

Rho_ref = np.zeros((4,1))
Rho_ref[0] = norm_03
Rho_ref[1] = norm_04
Rho_ref[2] = norm_07
Rho_ref[3] = norm_08

G = np.zeros((4,4))

G[0, 0:3] = (-1.0 /norm_03) * np.transpose(r_ecef_03 - r_ecef_0)
G[1, 0:3] = (-1.0 /norm_04) * np.transpose(r_ecef_04 - r_ecef_0)
G[2, 0:3] = (-1.0 /norm_07) * np.transpose(r_ecef_07 - r_ecef_0)
G[3, 0:3] = (-1.0 /norm_08) * np.transpose(r_ecef_08 - r_ecef_0)

G[:, 3] = 3E+08 * np.ones((1,4))

Rho_diff = Rho_meas - Rho_ref

delta_x = np.linalg.inv(G) @ Rho_diff

print(' The corretion in ECEF frame is: ')
print(delta_x)


# import datetime

# time_in_millis = 1650242392415
# dt = datetime.datetime.fromtimestamp(time_in_millis / 1000.0, tz=datetime.timezone.utc)

# print (dt)

# # Rough computation: 
# WEEKSEC = 604800
# C = 3E+08
# timeNanos = 283854068000000  # nanoseconds
# biasNanos = -0.3886222839355469  # sub-nanoseconds bias
# fullBiasNanos = -1333993756347467500

# t_RX_GPS = timeNanos - (fullBiasNanos + biasNanos)

# print(' t_RX_GPS =  ', t_RX_GPS)

# t_RX_GPS_mod = t_RX_GPS % (WEEKSEC *1E+09)

# print(' t_RX_GPS_mod =  ', t_RX_GPS_mod)

# receivedSvTimeNanos = 88810344022050
# timeOffsetNanos = 0.0

# t_TX_GPS = receivedSvTimeNanos + timeOffsetNanos

# print(' t_TX_GPS =  ', t_TX_GPS)

# delta_T = t_RX_GPS_mod - t_TX_GPS

# print(' Delta _T = ', delta_T * 1E-06, ' millseconds')

# rho = C * delta_T * 1E-09

# print(' Pseudo range is: ', rho)


# a = 6378137
# f = 1/298.257223563

# b = a * (1- f)

# e_sq = (a**2 - b**2)/ a**2

# print('E_sq = ', e_sq)