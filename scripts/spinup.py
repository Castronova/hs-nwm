#!/usr/bin/env python3

import os
import tqdm
import glob
import time
import docker
import shutil
import f90nml
import datetime
import pickle
import plot
import argparse
import numpy
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.dates as mdates
import logging
import xarray

from pandas.plotting import register_matplotlib_converters

register_matplotlib_converters()

log = logging.getLogger(__name__)

pwd = os.path.abspath(os.path.dirname(__file__))
domain = os.path.join(pwd, "DOMAIN")
output = os.path.join(pwd, "OUTPUT")
forcing = os.path.join(pwd, "FORCING")
restarts = os.path.join(pwd, "RESTARTS")
configs = os.path.join(pwd, "configs")
hydro = os.path.join(configs, "hydro.namelist")
hrldas = os.path.join(configs, "namelist.hrldas")
target = "/home/docker/RUN"
pkl = "spinup/streamflow.pkl"


def prepare():
    """
    Preparation steps for simulation run
    """

    # clean the output
    if os.path.exists(output):
        # rename using last modified date
        mod_time = datetime.datetime.strptime(
            time.ctime(max(os.path.getmtime(root) for root, _, _ in os.walk(output))),
            "%c",
        ).strftime("%Y%m%d_%H_%M")

        print("Found existing output directory")
        print(f"Renaming to {output}.{mod_time}")
        shutil.move(output, f"{output}.{mod_time}")
    os.mkdir(output)

    # make a directory for the spinup namelists and restarts
    if not os.path.exists("spinup"):
        log.info("creating directory for spinup files")
        os.mkdir("spinup")
    else:
        log.info("using existing spinup directory")

    # create output file
    if not os.path.exists("spinup/streamflow.pkl"):
        with open(pkl, "wb") as f:
            pickle.dump({}, f)


def run(n, image, outlet_id):

    for i in range(n):

        # execute preparation steps
        prepare()

        # read the lsm namelist
        lsm_nml = f90nml.read(hrldas)
        lsm = lsm_nml["noahlsm_offline"]
        lsm_nml_path = os.path.join(pwd, "spinup", "namelist.hrldas")
        hyd_nml = f90nml.read(hydro)
        hyd_nml_path = os.path.join(pwd, "spinup", "hydro.namelist")

        # get the latest restarts
        lsm_restarts = glob.glob(f"{restarts}/*RESTART*")
        hyd_restarts = glob.glob(f"{restarts}/*HYDRO_RST*")

        if len(lsm_restarts) == 0:
            # first run: comment out restart lines
            log.info("could not find existing LSM restart file, skipping")
            lsm_nml["noahlsm_offline"].pop("RESTART_FILENAME_REQUESTED")
        else:
            lsm_latest = max(lsm_restarts, key=os.path.getctime)
            log.info(f"Using LSM restart: {os.path.basename(lsm_latest)}")
            lsm_nml["noahlsm_offline"][
                "RESTART_FILENAME_REQUESTED"
            ] = f"./RESTARTS/{os.path.basename(lsm_latest)}"

        if len(hyd_restarts) == 0:
            # first run: comment out restart lines
            log.info("could not find existing HYDRO restart file, skipping")
            hyd_nml["hydro_nlist"].pop("RESTART_FILE")
        else:
            hyd_latest = max(hyd_restarts, key=os.path.getctime)
            log.info(f"Using HYDRO restart: {os.path.basename(hyd_latest)}")
            hyd_nml["hydro_nlist"][
                "RESTART_FILE"
            ] = f"./RESTARTS/{os.path.basename(hyd_latest)}"
            hyd_nml["hydro_nlist"]["GW_RESTART"] = 1

        # !! DEBUGGING !!
        # lsm_nml["noahlsm_offline"]["khour"] = 24

        # write the new namelist files
        log.info(f"writing namelists to: {lsm_nml_path}, {hyd_nml_path}")
        if os.path.exists(lsm_nml_path):
            os.remove(lsm_nml_path)
        lsm_nml.write(lsm_nml_path)
        if os.path.exists(hyd_nml_path):
            os.remove(hyd_nml_path)
        hyd_nml.write(hyd_nml_path)

        # determine the start and end time of the simulation
        start_dt = datetime.datetime(
            lsm["start_year"],
            lsm["start_month"],
            lsm["start_day"],
            lsm["start_hour"],
            lsm["start_min"],
        )
        curr_dt = start_dt
        if "kday" in lsm:
            end_dt = start_dt + datetime.timedelta(days=lsm["kday"])
        else:
            end_dt = start_dt + datetime.timedelta(hours=lsm["khour"])
        simulation_total_sec = (end_dt - start_dt).total_seconds()

        # make a directory for restarts
        if not os.path.exists("RESTARTS"):
            os.mkdir("RESTARTS")

        # docker args
        cmd = '-c "./wrf_hydro.exe && mv *_DOMAIN* OUTPUT"'  # >/dev/null 2>&1'
        mnt = {
            domain: {"bind": f"{target}/DOMAIN", "mode": "rw"},
            output: {"bind": f"{target}/OUTPUT", "mode": "rw"},
            forcing: {"bind": f"{target}/FORCING", "mode": "rw"},
            restarts: {"bind": f"{target}/RESTARTS", "mode": "rw"},
            hyd_nml_path: {"bind": f"{target}/hydro.namelist", "mode": "rw"},
            lsm_nml_path: {"bind": f"{target}/namelist.hrldas", "mode": "rw"},
        }

        # spawn the simulation container
        client = docker.from_env()
        container = client.containers.run(
            image, command=cmd, volumes=mnt, remove=False, detach=True
        )

        log.info("Beginning Simulation")
        log.info(f"Start: {start_dt}")
        log.info(f"End: {end_dt}")

        # initialize the progress bar
        pbar = tqdm.tqdm(
            total=simulation_total_sec,
            bar_format=("{l_bar}{bar}|[{elapsed}" "<{remaining}{postfix}"),
        )
        pbar.set_description(f"Simulation {i+1}:")

        # wait for container to start
        time.sleep(5)

        # begin simulation loop
        while len(client.containers.list(filters={"id": container.id})) > 0:

            # get the last LDAS file that was created
            out_files = glob.glob(f"{output}/*LDAS*")
            latest = os.path.basename(max(out_files, key=os.path.getctime))

            # update the progress bar
            old_dt = curr_dt
            curr_dt = datetime.datetime.strptime(latest.split(".")[0], "%Y%m%d%H%M")
            diff = (curr_dt - old_dt).total_seconds()
            pbar.set_postfix(sim_time=curr_dt)
            pbar.update(diff)
            time.sleep(1)
        pbar.close()

        # save simulation data (i.e. outlet predictions and Restarts)
        # to the spinup directory
        save_to_spinup_dir(outlet_id)


def save_to_spinup_dir(outlet_id=None):
    dt_str = datetime.datetime.now().strftime("%Y.%m.%d.%H.%M.%S")

    # move the last restart files
    lsm_restarts = glob.glob(f"{output}/RESTART*")
    hyd_restarts = glob.glob(f"{output}/HYDRO_RST*")

    lsm_latest = max(lsm_restarts, key=os.path.getctime)
    hyd_latest = max(hyd_restarts, key=os.path.getctime)
    shutil.copyfile(lsm_latest, f"{restarts}/{dt_str}.{os.path.basename(lsm_latest)}")
    shutil.copyfile(hyd_latest, f"{restarts}/{dt_str}.{os.path.basename(hyd_latest)}")

    # read data for reach
    # save to pickle
    # plot when completely finished spinup run

    # save outlet streamflow
    if outlet_id is not None:
        try:
            p = plot.Plot("DOMAIN", "OUTPUT")
            x, y = p.get_reach_data(outlet_id)

            with open(pkl, "rb") as f:
                p = pickle.load(f)
            p[dt_str] = dict(x=x, y=y)
            with open(pkl, "wb") as f:
                pickle.dump(p, f)
        except Exception:
            print(f'ERROR saving results for reachID {outlet_id}') 


def make_plot():
    if not os.path.exists(pkl):
        log.critical(f"ERROR could not find: {pkl}")
        return None

    # open the pickle file
    with open(pkl, "rb") as f:
        dat = pickle.load(f)

    fig, ax = plt.subplots()
    colors = matplotlib.cm.jet(numpy.linspace(0, 1, len(dat.keys())))
    i = 0

    for k in sorted(dat):
        x = dat[k]["x"]
        y = dat[k]["y"]
        ax.plot_date(x, y, label=i, linestyle="-", marker="", color=colors[i])
        i += 1
    ax.set(xlabel="date")
    ax.grid()
    ax.legend(loc="best")

    # format the ticks
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m/%d/%Y"))
    fig.autofmt_xdate()

    plt.savefig("spinup/all_streamflow.png")


def list_gauges(domain_dir):

    # load routelink data, i.e. stream definitions
    routelink = os.path.join(domain_dir, "Route_Link.nc")
    route_df = xarray.open_dataset(routelink).to_dataframe()
    route_df.gages = route_df.gages.str.decode("utf-8").str.strip()

    filtered = route_df.loc[route_df["gages"] != ""]
    cols = ["order", "link", "gages", "lat", "lon"]
    data = filtered[cols].sort_values(by=["order"])
    data = data.rename(
        index=str,
        columns={
            "order": "StreamOrder",
            "link": "ReachID",
            "gages": "USGS_GageID",
            "lat": "lat-midpoint",
            "lon": "lon-midpoint",
        },
    )
    data.reset_index(inplace=True)
    return data


def select_outlet_reach():
    res = ""
    outlet_id = None
    while res != "q":
        data = list_gauges(domain_dir=args.domain_dir)
        print("\nFound the following reaches with USGS gauges:")
        print(data)
        res = input(
            'Enter the RowID of the reach you would like to save data for, or enter "q" to exit this menu: '
        )

        # parse response
        res = res.strip().lower()
        if res == "q":
            break

        if res.isdigit():
            idx = int(res)
        else:
            print(f"{res} is an invalid input.")
            continue

        if idx in data.index:
            outlet_id = data.iloc[idx]["ReachID"]
            break
        else:
            print(f"{res} is an invalid index.")
            continue

    return outlet_id


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # WRF-Hydro Directory Arguments
    parser.add_argument(
        "--domain-dir", default="DOMAIN", help="path to WRF-Hydro DOMAIN data"
    )

    # WRF-Hydro Simulation Arguments
    parser.add_argument("-n", type=int, default=0, help="number of simulations to run")
    parser.add_argument("-p", action="store_true", help="generate plots")
    parser.add_argument(
        "-i", default="cuahsi/wrfhydro-nwm:5.0.3", help="docker image to use"
    )
    parser.add_argument("-q", action="store_true", help="run in quiet mode")
    parser.add_argument(
        "-s",
        "--save-at-outlet",
        nargs=1,
        help="save output data and results for outlet id.",
    )

    # options: save all output (default delete)
    # plot specific run ids

    args = parser.parse_args()

    # set logging to INFO by default. Change to CRITICAL
    # if running in quiet mode
    log.setLevel(logging.INFO)
    if args.q:
        log.setLevel(logging.CRITICAL)
    else:
        log.addHandler(logging.StreamHandler())

    if args.save_at_outlet is None:
        # save was not provided.
        res = input(
            "Save-at-Outlet option was not provided, would you like to select reach to save [y/N]?"
        )
        if res.lower() == "y":
            outlet_id = select_outlet_reach()
            args.save_at_outlet = outlet_id

    if args.n:
        run(args.n, args.i, args.save_at_outlet)
    if args.p:
        make_plot()
#    if args.save_at_outlet is not None:
#        save_to_spinup_dir(args.save_at_outlet)
