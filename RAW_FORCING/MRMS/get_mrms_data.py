#!/usr/bin/env python3

import os
import sys
import wget
import datetime
from tqdm import tqdm

output_dir = 'archived_grib2'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

base_url = 'https://mtarchive.geol.iastate.edu'
st = datetime.datetime(2018, 4, 1, 0, 0)
et = datetime.datetime(2018, 9, 30, 23, 0)
ts_min = 60

total_steps = (et-st) / datetime.timedelta(minutes=ts_min)

pbar = tqdm(total=total_steps)
errors = []
while st < et:
    url = f'{base_url}/{st.strftime("%Y/%m/%d")}/' + \
           'mrms/ncep/GaugeCorr_QPE_01H'

    next_day = st + datetime.timedelta(days=1)

    # loop over hours
    while st.day != next_day.day:

        fname = 'GaugeCorr_QPE_01H_00.00_' + \
                f'{st.strftime("%Y%m%d-%H%M%S")}.grib2.gz'

        # download if file doesn't already exist
        fpath = os.path.join(output_dir, fname)
        if not os.path.exists(fpath):
            try:
                wget.download(f'{url}/{fname}', fpath, bar=None)

            except Exception as e:
                errors.append(f'URL: {url}/{fname}')

        # increment by timestep
        pbar.update()
        st += datetime.timedelta(minutes=ts_min)

if len(errors) > 0:
    print('Failed to located the following files:')
    for e in errors:
        print(e)

