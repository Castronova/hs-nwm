#!/usr/bin/env python3

import os
import sys
import subprocess
import datetime
from tqdm import tqdm
from getpass import getpass


uname = input("Enter GES Username: ").strip()
passw = getpass()


output_dir = "archived_grib2"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

with open("urls.txt", "r") as f:
    urls = f.readlines()


pbar = tqdm(total=len(urls))
errors = []
for url in urls:

    # build the url
    url = url.strip()

    fname = os.path.basename(url)
    # download if file doesn't already exist
    fpath = os.path.join(output_dir, fname)
    if not os.path.exists(fpath):
        cmd = (
            "wget --quiet "
            "--auth-no-challenge on "
            "--content-disposition "
            f"--user {uname} "
            f"--password {passw} "
            f"-P {output_dir} "
            f"{url}"
        )

        proc = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
        output, error = proc.communicate()

        if error is not None:
            errors.append(f"URL: {url}")

    # increment by timestep
    pbar.update()

if len(errors) > 0:
    print("Failed to located the following files:")
    for e in errors:
        print(e)
