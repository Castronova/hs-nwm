#!/usr/bin/env python3

import cdo
import xarray
from glob import glob
from pathlib import Path
from pynhd import NLDI, WaterData
import pandas
import numpy
import folium
import jenkspy
import matplotlib.pyplot as plt
import requests
import pytz
from datetime import datetime


class Simulation:
    def __init__(self, simulation_label, domain_dir, output_dir):
        self.simulation_label = simulation_label
        self.domain_dir = domain_dir
        self.output_dir = output_dir
        self.routelink_path = Path(self.domain_dir) / "Route_Link.nc"

        self.consolidated_chrout = f"{self.simulation_label}_chrtout_merged.nc"
        self.consolidated_gwout = f"{self.simulation_label}_gwout_merged.nc"

        self._chrtout = None
        self._gwout = None
        self._routelink = None
        self._reach_colors = None

    @property
    def routelink(self):
        if self._routelink is None:
            # load routelink nc as a dataframe
            rtlink = xarray.open_dataset(self.routelink_path).to_dataframe()
            rtlink.gages = rtlink.gages.str.decode("utf-8").str.strip()
            self.routelink = rtlink
        return self._routelink

    @routelink.setter
    def routelink(self, value):
        self._routelink = value

    @property
    def chrtout(self):
        if self._chrtout is None:
            if not Path(self.consolidated_chrout).exists():
                print('Building consolidated CHRTOUT...', end='', flush=True)
                chrtout_files = glob(f"{self.output_dir}/*CHRTOUT*")
                cdo.Cdo().cat(input=chrtout_files, output=self.consolidated_chrout)
                print('done')
            else: 
                print('Loading cached CHRTOUT (consolidated)')

            self.chrtout = xarray.open_dataset(self.consolidated_chrout)
        return self._chrtout

    @chrtout.setter
    def chrtout(self, value):
        self._chrtout = value

    @property
    def gwout(self):
        if self._gwout is None:
            if not Path(self.consolidated_gwout).exists():
                print('Building consolidated GWOUT...', end='', flush=True)
                gwout_files = glob(f"{self.output_dir}/*GWOUT*")
                cdo.Cdo().cat(input=gwout_files, output=self.consolidated_gwout)
                print('done')
            else:
                print('Loading cached GWOUT (consolidated)')

            self.gwout = xarray.open_dataset(self.consolidated_gwout)
        return self._gwout

    @gwout.setter
    def gwout(self, value):
        self._gwout = value

    @property
    def usgs_gages(self):
        """
        Lists usgs gauges that exist in the model domain

        returns:
            pandas dataframe containing USGS gauges and their corresponding linkIDs
        """

        filtered = self.routelink.loc[self.routelink["gages"] != ""]
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

    @property
    def reaches(self):
        """
        Lists all reaches within the domain

        returns:
            pandas dataframe containing simulation reaches and their upstream and downstream counterparts.
        """

        cols = ["order", "link", "gages", "lat", "lon", "to"]
        data = self.routelink[cols].sort_values(by=["order"])
        data = data.rename(
            index=str,
            columns={
                "order": "stream_order",
                "link": "comid",
                "to": "downstream_comid",
                "gages": "usgs_gageid",
                "lat": "lat-midpoint",
                "lon": "lon-midpoint",
            },
        )
        data.reset_index(inplace=True)
        return data

    def __prepare_reach_colors(self,
                               start_time="1900-01-01",
                               end_time="2030-01-01"):

        # slice the data to the specified time range
        sliced_reach_ds = self.chrtout.sel(time=slice(start_time, end_time))

        # combine total streamflow and feature id so we can use it to stylize the reaches
        self._reach_colors = sliced_reach_ds.streamflow.sum(axis=0).to_dataframe()
        breaks = jenkspy.jenks_breaks(self._reach_colors['streamflow'].values,
                                      nb_class=5)

        self._reach_colors['color'] = pandas.cut(self._reach_colors['streamflow'],
                                                 bins=breaks,
                                                 labels=['lightblue',
                                                         'cadetblue',
                                                         'blue',
                                                         'darkblue',
                                                         'purple'],
                                                 include_lowest=True)
        self._reach_colors.color = self._reach_colors.color.astype('str')

        # replace the color for reaches with zero flow
        self._reach_colors.loc[self._reach_colors.streamflow == 0.0,
                               'color'] = 'red'

    def __style_reach(self, feature):
        # define classes to style the data
        fid = int(feature['properties']['comid'])
        return {'color': self._reach_colors.loc[fid].color}

    def plot_domain(self,
                    station_id,
                    start_time="1900-01-01",
                    end_time="2030-01-01"):
        """
        args:
            station_id: USGS station id as a string, e.g. "05454300"
            start_time: YYYY-MM-DD string used to clip temporally
            end_time:   YYYY-MM-DD string used to clip temporally
        returns:
            m: Folium Map Object
        """

        # get NWM reaches and convert comid to int64 to prevent
        # type mismatch during future merges
        nldi = NLDI()
        basin = nldi.get_basins(station_id)
        mr = WaterData("nhdflowline_network")
        nhd = mr.bybox(basin.geometry[0].bounds)
        reaches = self.reaches
        nhd.comid = nhd.comid.astype("int64")

        # generate the reach colors
        self.__prepare_reach_colors(start_time=start_time, end_time=end_time)

        # merge NWM reaches with NHD+ reaches
        nhd = pandas.merge(nhd, reaches, how="inner", left_on="comid", right_on="comid")
        nhd = nhd.join(self._reach_colors, on="comid")

        # create empty lists in upstream_comids column
        nhd["upstream_comids"] = numpy.empty((len(nhd), 0)).tolist()

        # loop through each row and append its comid to its downstream id
        for idx, row in nhd.iterrows():
            # get the comid for the current reach
            comid = row.comid

            # get the downstream comid for the current reach
            downstream_comid = row.downstream_comid

            # append the current comid to the upstream_comids column of the downstream comid
            # todo: there's got to be a better way of doing this!!!
            if (
                len(nhd.loc[nhd.comid == downstream_comid, "upstream_comids"].values)
                > 0
            ):
                nhd.loc[nhd.comid == downstream_comid, "upstream_comids"].values[
                    0
                ].append(comid)

        # set upstream comids to -9999 if none exist
        nhd.loc[nhd["upstream_comids"].str.len() == 0, "upstream_comids"] = "-9999"

        # make comids strings so they appear without commas in popup
        nhd.comid = nhd.comid.astype(str)
        nhd.downstream_comid = nhd.downstream_comid.astype(str)

        usgs_gages = self.usgs_gages
        m = folium.Map(tiles="Stamen terrain", zoom_start=11)

        # add watershed
        # watershed = geopandas.GeoDataFrame.from_file('../spatial/watershed.shp')
        watershed_json = basin.to_crs(epsg="4326").to_json()
        w = folium.features.GeoJson(data=watershed_json)
        m.add_child(w)

        nhd_json = nhd.to_crs(epsg="4326").to_json()
        w = folium.features.GeoJson(
            data=nhd_json,
            style_function=self.__style_reach,
            popup=folium.features.GeoJsonPopup(
                fields=[
                    "comid",
                    "gnis_id",
                    "gnis_name",
                    "stream_order",
                    "downstream_comid",
                    "upstream_comids",
                    "streamflow",
                ]
            ),
            highlight_function=lambda feature: {"fillcolor": "green", "color": "green"},
        )
        m.add_child(w)

        for idx in range(0, len(usgs_gages)):
            dat = usgs_gages.iloc[idx]
            html = f"""
                   <html>
                    <b>NWM ID: </b>{dat.ReachID}<br>
                    <b>USGS ID: </b>{dat.USGS_GageID}<br>
                    <b>Stream Order: </b>{dat.StreamOrder}<br>
                   </html>
                   """
            popup = folium.Popup(folium.Html(html, script=True), max_width=2650)
            folium.Marker(
                location=[float(dat["lat-midpoint"]), float(dat["lon-midpoint"])],
                icon=folium.Icon(color="blue"),
                popup=popup,
            ).add_to(m)

        # Set the map extent (bounds) to the extent of the sites
        m.fit_bounds(m.get_bounds())

        # display map details
        print(f"\nWatershed outlet station ID: {station_id}")
        print(f"Reaches are symbolized based on the total flow within the time range:")
        print(f"  Start Date: {start_time}")
        print(f"  End Date: {end_time}\n")

        # convert comid back to int64 for later operations
        nhd.comid = nhd.comid.astype("int64")

        return m

    def plot_chrtout(self,
                     wrfhydro_link_ids=[],
                     usgs_station_ids=[],
                     start_time="1900-01-01",
                     end_time="2030-01-01",
                     tz='US/Eastern'):
        """
        args:
            wrfhydro_link_ids: list of nwm link ids to plot
            usgs_station_ids: list of USGS station ids to plot, e.g. "05454300"
            start_time: YYYY-MM-DD string used to clip temporally
            end_time:   YYYY-MM-DD string used to clip temporally
        returns:
            m: matplotlib object
        """
        
        plt.ioff()
        
        # create plot
        fig = plt.figure(figsize=(8, 5))
        ax = plt.gca()
        
        tz = pytz.timezone(tz)
        st = datetime.strptime(start_time, '%Y-%m-%d').replace(tzinfo=tz)
        et = datetime.strptime(end_time, '%Y-%m-%d').replace(tzinfo=tz)

        if len(usgs_station_ids) == 0:
            print('No USGS gauges provided, skipping')
        else:
            # plot observed data
            url = (
                   'https://nwis.waterservices.usgs.gov/nwis/iv/?'
                   f'format=json&sites={",".join(usgs_station_ids)}'
                   f'&startDT={st.strftime("%Y-%m-%dT%H:%M%z")}'
                   f'&endDT={et.strftime("%Y-%m-%dT%H:%M%z")}'
                   '&parameterCd=00060&siteStatus=all'
                  )
            res = requests.get(url)
            if res.status_code != 200:
                print(f'Error collecting USGS data; url={url}')
                return None
            usgs_dat = res.json()

            usgs_series = {}
            for i in range(0, len(usgs_dat['value']['timeSeries'])):
                ts = usgs_dat['value']['timeSeries'][i]
                label = (f"{ts['sourceInfo']['siteCode'][0]['value']}/"
                         f"{ts['variable']['variableCode'][0]['value']} - "
                         f"{ts['variable']['variableCode'][0]['vocabulary']}"
                         )
                x = []
                y = []
                for entry in ts['values'][0]['value']:
                    x.append(datetime.strptime(entry['dateTime'],
                             '%Y-%m-%dT%H:%M:%S.%f%z'))
                    # convert from cfs to cms
                    y.append(float(entry['value']) / 35.314666212661)
                usgs_series[label] = dict(x=x, y=y)
        
            # plot USGS data
            for label in usgs_series.keys():
                ax.plot_date(usgs_series[label]['x'],
                             usgs_series[label]['y'],
                             label=label, 
                             linestyle='-',
                             marker='')

        if len(wrfhydro_link_ids) == 0:
            print('No WRF-Hydro link IDs provided, skipping')
        else:
            # plot WRF-Hydro data
            for link in wrfhydro_link_ids:
                try:
                    timerange = slice(start_time, end_time)
                    (self.chrtout.sel(feature_id=int(link),
                                      time=timerange).streamflow).plot(ax=ax, label=f'{link} - WRF-Hydro')
                except Exception:
                    print(f'Failed to load WRF-Hydro reach {link}')

        # finish configuring the plot
        ax.set_ylabel('Streamflow (cms)')
        ax.set_xlabel('Date')
        _ = ax.legend()

        fig.autofmt_xdate()
        plt.tight_layout()
        
        return (fig, ax, plt)

