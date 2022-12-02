import numpy as np

time_nanos = 636005573000000
time_offset_nanos = 0.0

bias_nanos =  -0.9884319305419922
full_bias_nanos = -1329943102840388004

num_nanos_per_week = 6.048e+14
week_number_nanos = np.floor(-full_bias_nanos/num_nanos_per_week) * num_nanos_per_week 

t_RX_GNSS = time_nanos + time_offset_nanos - (bias_nanos + full_bias_nanos)

t_RX = np.mod(t_RX_GNSS, num_nanos_per_week)

print ("Reception time is: ", t_RX)

svid = 2

t_TX = 19108334548071

TX_RX_diff = t_RX - t_TX

print (" Time difference in nanoseconds is: ", TX_RX_diff)

c = 299792458
pseudorange_m = (TX_RX_diff /1e09) * c 

print (" Pseudo range is: ", pseudorange_m)