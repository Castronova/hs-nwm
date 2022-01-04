#!/usr/bin/env python3

import os
import gzip
import glob
import shutil

archive_path = "archived_grib2"
output_path = "input_files"
for p in [archive_path, output_path]:
    if not os.path.exists(p):
        os.makedirs(p)


for mrms_file in glob.glob(f"{archive_path}/*gz"):
    print(f"processing: {mrms_file}")
    outpath = os.path.join(
        output_path, ".".join(os.path.basename(mrms_file).split(".")[:-1])
    )
    with gzip.open(mrms_file, "r") as f_in:
        with open(outpath, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
#    shutil.move(mrms_file, os.path.join(archive_path, mrms_file))
