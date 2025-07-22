# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 14:56:23 2020

@author: raphael.dutrieux

"""


import numpy as np
import datetime
import random
import geopandas as gp
import xarray as xr
import click
import matplotlib
from rasterio import Affine, transform
import dask.array as da

from fordead.import_data import TileInfo, import_stackedmaskedVI, import_stress_data, import_stacked_anomalies, import_coeff_model, import_binary_raster, import_first_detection_date_index, import_dieback_data, import_soil_data
from fordead.model_vegetation_index import compute_HarmonicTerms

from fordead.results_visualisation import select_and_plot_time_series
import warnings
warnings.filterwarnings("ignore", message='The specified chunks separate the stored chunks along dimension "x" starting at index 100. This could degrade performance. Instead, consider rechunking after loading.')
warnings.filterwarnings("ignore", message="Converting non-nanosecond precision datetime values to nanosecond precision. This behavior can eventually be relaxed in xarray, as it is an artifact from pandas which is now beginning to support non-nanosecond precision values. This warning is caused by passing non-nanosecond np.datetime64 or np.timedelta64 values to the DataArray or Variable constructor; it can be silenced by converting the values to nanosecond precision ahead of time.")

@click.group()
def graph_series():
    """Main entrypoint."""
    
@graph_series.command(name='graph_series')
@click.option("-o", "--data_directory", type = str, help = "Path of the directory containing results from the region of interest")
@click.option("-x","--x", type = float, default = None, help = "x coordinate in the Sentinel-2 data CRS of the pixel of interest. Not used if shape_path parameter is used.")
@click.option("-y","--y", type = float, default = None, help = "y coordinate in the Sentinel-2 data CRS of the pixel of interest. Not used if shape_path parameter is used.")
@click.option("--shape_path", type = str, help = "Path to shapefile containing points whose data will be plotted. If None, indexes or coordinates for x and y can be given")
@click.option("--name_column", type = str, default = "id", help = "Name of the column containing the name of the point, used to name the exported image. Not used if pixel is selected from indexes or coordinates")
@click.option("--ymin", type = float, default = 0, help = "ymin limit of graph")
@click.option("--ymax", type = float, default = 2, help = "ymax limit of graph")
@click.option("--chunks", type = int, default = None, help = "Chunk length to import data as dask arrays and save RAM, advised if computed area in data_directory is large")
def cli_vi_series_visualisation(data_directory, x = None, y = None, shape_path = None, name_column = "id", ymin = 0, ymax = 2, chunks = None):
    """
    From previously computed results, graphs the results for specific pixels showing the vegetation index for each dates, the model and the detection.
    By specifying 'shape_path' and 'name_column' parameters, it can be used with a shapefile containing points with a column containing a unique ID used to name the exported image.
    By specifying 'x' and 'y' parameters, it can be used with coordinates in the Sentinel-2 data CRS.
    If neither shape_path or x and y parameters are specified, the user will be prompted to give coordinates in the system of projection of the tile. Graphs can also be plotted for random pixels inside the forest mask.
    The user can also choose to specify pixels by their indices from the top left hand corner of the computed area (If only a small region of interest was computed (for example by using extent_shape_path parameter in the step 01_compute_masked_vegetationindex (https://fordead.gitlab.io/fordead_package/docs/user_guides/english/01_compute_masked_vegetationindex/)), then create a timelapse on this whole region of interest (https://fordead.gitlab.io/fordead_package/docs/user_guides/Results_visualization/#creer-des-timelapses), then these indices correspond to the indices in the timelapse) 
    The graphs are exported in the data_directory/TimeSeries directory as png files.
    See details here : https://fordead.gitlab.io/fordead_package/docs/user_guides/Results_visualization/
    \f
    Parameters
    ----------
    data_directory
    x
    y
    shape_path
    name_column
    ymin
    ymax
    chunks

    Returns
    -------

    """
    vi_series_visualisation(data_directory, x, y, shape_path, name_column, ymin, ymax, chunks)



def vi_series_visualisation(data_directory, x= None, y = None, shape_path = None, name_column = "id", ymin = 0, ymax = 2, chunks = None):
    
    """
    From previously computed results, graphs the results for specific pixels showing the vegetation index for each dates, the model and the detection.
    By specifying 'shape_path' and 'name_column' parameters, it can be used with a shapefile containing points with a column containing a unique ID used to name the exported image.
    By specifying 'x' and 'y' parameters, it can be used with coordinates in the Sentinel-2 data CRS.
    If neither shape_path or x and y parameters are specified, the user will be prompted to give coordinates in the system of projection of the tile. Graphs can also be plotted for random pixels inside the forest mask.
    The user can also choose to specify pixels by their indices from the top left hand corner of the computed area (If only a small region of interest was computed (for example by using extent_shape_path parameter in the step 01_compute_masked_vegetationindex (https://fordead.gitlab.io/fordead_package/docs/user_guides/english/01_compute_masked_vegetationindex/)), then create a timelapse on this whole region of interest (https://fordead.gitlab.io/fordead_package/docs/user_guides/Results_visualization/#creer-des-timelapses), then these indices correspond to the indices in the timelapse) 
    The graphs are exported in the data_directory/TimeSeries directory as png files.
    See details here : https://fordead.gitlab.io/fordead_package/docs/user_guides/Results_visualization/
    
    Parameters
    ----------
    data_directory : str
        Path of the directory containing results from the region of interest
    x : int
        x coordinate in the Sentinel-2 data CRS of the pixel of interest. Not used if shape_path parameter is used.
    y : int
        y coordinate in the Sentinel-2 data CRS of the pixel of interest. Not used if shape_path parameter is used.
    shape_path: str
        Path to shapefile containing points whose data will be plotted. Not used if pixel is selected from x and y parameters
    name_column: str
        Name of the column containing the name of the point, used to name the exported image. Not used if pixel is selected from x and y parameters
    ymin: float
        ymin limit of graph
    ymax: float
        ymax limit of graph
    chunks: int
        Chunk length to import data as dask arrays and save RAM, advised if computed area in data_directory is large

    """
    tile = TileInfo(data_directory)
    tile = tile.import_info()

    # IMPORTING ALL DATA
    stack_vi, stack_masks = import_stackedmaskedVI(tile,chunks = chunks)
    stack_vi["DateNumber"] = ("Time", np.array([(datetime.datetime.strptime(date, '%Y-%m-%d')-datetime.datetime.strptime('2015-01-01', '%Y-%m-%d')).days for date in np.array(stack_vi["Time"])]))
    coeff_model = import_coeff_model(tile.paths["coeff_model"],chunks = chunks)
    if tile.parameters["stress_index_mode"] is not None:
        stress_data = import_stress_data(tile.paths,chunks = chunks)
    else:
        stress_data = None

    first_detection_date_index = import_first_detection_date_index(tile.paths["first_detection_date_index"],chunks = chunks)
    if tile.parameters["soil_detection"]:
        soil_data = import_soil_data(tile.paths,chunks = chunks)
    else:
        soil_data=None
    dieback_data = import_dieback_data(tile.paths,chunks = chunks)
    forest_mask = import_binary_raster(tile.paths["forest_mask"],chunks = chunks)
    tile.getdict_datepaths("Anomalies",tile.paths["AnomaliesDir"])
    anomalies = import_stacked_anomalies(tile.paths["Anomalies"],chunks = chunks)
    
    tile.add_dirpath("series", tile.data_directory / "TimeSeries")
    tile.save_info()
    xx = np.array(range(int(stack_vi.DateNumber.min()), int(stack_vi.DateNumber.max())))
    xxDate = [np.datetime64(datetime.datetime.strptime("2015-01-01", '%Y-%m-%d').date()+ datetime.timedelta(days=int(day))) for day in xx]
    # xxDate = [date for date in xxDate if tile.parameters["ignored_period"] is None or (date.astype(str)[5:] > min(tile.parameters["ignored_period"]) and date.astype(str)[5:] < max(tile.parameters["ignored_period"]))]
    # xx = np.array([(date - np.datetime64(datetime.datetime.strptime("2015-01-01", '%Y-%m-%d').date())).astype(int) for date in xxDate])
    
    harmonic_terms = np.array([compute_HarmonicTerms(DateAsNumber) for DateAsNumber in xx])
    harmonic_terms = xr.DataArray(harmonic_terms, coords={"Time" : xxDate, "coeff" : [1,2,3,4,5]},dims=["Time","coeff"])
    
    if tile.parameters["correct_vi"]:
        stack_vi = stack_vi + tile.correction_vi
  
    stack_vi.coords["Time"] = tile.dates.astype("datetime64[D]")
    stack_masks.coords["Time"] = tile.dates.astype("datetime64[D]")

    if shape_path is not None:
        matplotlib.use('Agg')
        shape = gp.read_file(shape_path)
        shape = shape.to_crs(crs = tile.raster_meta["crs"])
        
        for point_index in range(len(shape)):
            id_point = shape.iloc[point_index][name_column]
            geometry_point = shape.iloc[point_index]["geometry"]
            print(id_point)
            row, col = transform.rowcol(tile.raster_meta["transform"],geometry_point.x,geometry_point.y)

            select_and_plot_time_series(col,row, forest_mask, harmonic_terms, coeff_model, first_detection_date_index, soil_data, dieback_data, stack_masks, stack_vi, anomalies, stress_data, tile, ymin, ymax, name_file = str(id_point))
    elif (x is not None) and (y is not None):
        row, col = transform.rowcol(tile.raster_meta["transform"],x,y)
        select_and_plot_time_series(col, row, forest_mask, harmonic_terms, coeff_model, first_detection_date_index, soil_data, dieback_data, stack_masks, stack_vi, anomalies, stress_data, tile, ymin, ymax)

    else:
        #Initialiser x,y
        # matplotlib.use('TkAgg')
        PixelsToChoose = np.where(forest_mask)
        
        mode = input("Type 'c' to input coordinates as coordinates in the system of projection of the tile\nType 'i' to input coordinates by positional indexing (x = col, y = row) \n[c/i]?")
        if mode not in ["c","i"]: raise Exception("Index or coordinate mode incorrect. Type 'c' for coordinate mode, 'i' for index mode")
      
        x=0
        while x!=-1:
            
            x=input("x ? ")
            if x=="":
                #PIXEL ALEATOIRE DANS LE MASQUE FORET
                PixelID=random.randint(0,PixelsToChoose[0].shape[0])
                x=PixelsToChoose[1][PixelID]
                y=PixelsToChoose[0][PixelID]

                select_and_plot_time_series(x,y, forest_mask, harmonic_terms, coeff_model, first_detection_date_index, soil_data, dieback_data, stack_masks, stack_vi, anomalies, stress_data, tile, ymin, ymax)


            elif x=="-1":
                #ARRET SI x = -1
                break
            else:
                x=int(x)
                y=int(input("y ? "))
                
                if mode == "c": y, x = transform.rowcol(tile.raster_meta["transform"],x,y)

                select_and_plot_time_series(x,y, forest_mask, harmonic_terms, coeff_model, first_detection_date_index, soil_data, dieback_data, stack_masks, stack_vi, anomalies, stress_data, tile, ymin, ymax)

    
if __name__ == '__main__':
    # print(dictArgs)
    cli_vi_series_visualisation()