## Script to read GPS nav messages in RINEX format and compute the satellite position and velocity in ECEF frame. 

import re

## The RINEX nav messages usually speicfy the data for the Keplerian orbit params

def splitEntriesByNegative(entries):

    split_entries = []

    for e in entries: 

        negative_matches = re.finditer(r'(\d)\-(\d)', e)

        if negative_matches is None:
            print('Itr is None')
            split_entries.append(e)
        else: 
            start_idx = 0
            nm = {}
            for nm in negative_matches:

                print(e[start_idx:nm.start()+1])

                split_entries.append(float(e[start_idx:nm.start()+1]))

                start_idx = nm.start() + 1

            # print('Itr is not none. Start idx is: ', start_idx)
            # print(e)
            # print(e[start_idx:])
            split_entries.append(e[start_idx:])
        
    return split_entries

nav_data = {}

# Step 1: Read in the NAV data file
with open('data/BREW00USA_R_20220311200_01H_GN.rnx') as f:
    lines = f.readlines()

sat_name = ''
g = 0

for l in lines: 

    entries = l.split(' ') 
    print(entries)
    if (l[0] == 'G'):
        g = 1
        sat_name = str(entries[0])

    if g == 1:
        entries = splitEntriesByNegative(entries)
        nav_data[sat_name] = {}
        nav_data[sat_name]['sv_clock_bias'] = float(entries[7])
        nav_data[sat_name]['sv_clock_drift'] = float(entries[8])
        g = g + 1
    elif g == 2:
        entries = splitEntriesByNegative(entries)
        nav_data[sat_name]['iode'] = float(entries[5])
        nav_data[sat_name]['M0'] = float(entries[8])
        g = g + 1
    elif g ==3:
        entries = splitEntriesByNegative(entries)
        nav_data[sat_name]['e'] = float(entries[5])
        nav_data[sat_name]['sqrt_A'] = float(entries[7])
        g = g + 1

print('---')
print(nav_data)
print('---')

# entries = ['4.433109425008E-04-9.663381206337E-12-8.663381206337E-12', '5.433109425008E-04', '7.433109425008E-04-9.433109425008E-04']

# arr = splitEntriesByNegative(entries)
# print('---')
# print(arr)

# first = string[0:split_str.start() + 1]
# second = string[split_str.end() - 2:]

# print(split_str.start())
# print(split_str.end())

# print (first)
# print (second)
# Step 2: Determine the time offset from the Issue of Data Ephemeris 

# Step 3: Determine the mean anomaly

# Step 4: Determine the true anomaly

# Step 5: Determine position in perifocal frame

# Step 6: Transform position to ECEF frame
