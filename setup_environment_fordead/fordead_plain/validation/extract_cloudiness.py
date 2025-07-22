# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 16:57:21 2023

@author: rdutrieux
"""

# import time
import click
# import geopandas as gp
from pathlib import Path
import pandas as pd
import numpy as np
from fordead.import_data import TileInfo, get_band_paths, get_cloudiness

@click.command(name='extract_cloudiness')
@click.option("--sentinel_dir", type = str,default = None, help = "Path of the directory containing Sentinel-2 data.", show_default=True)
@click.option("--export_path", type = str,default = None, help = "Path to write csv file with extracted cloudiness", show_default=True)
@click.option("-t","--tile_selection", multiple=True, type = str, default = None, help = "List of tiles from which to extract reflectance (ex : -t T31UFQ -t T31UGQ). If None, all tiles are extracted.", show_default=True)
@click.option("--sentinel_source", type = str,default = "THEIA", help = "Source of data, can be 'THEIA' et 'Scihub' et 'PEPS'", show_default=True)
def cli_extract_cloudiness(sentinel_dir, export_path, tile_selection, sentinel_source):
    """
    
    \f

    """
    
    extract_cloudiness(**locals())


def extract_cloudiness(sentinel_dir, export_path, tile_selection = None, sentinel_source = "THEIA"):
    """
    
    For each acquisition, extracts percentage of pixels in the mask provided by the Sentinel-2 data provider.
    For THEIA, all pixels different to 0 in the CLM mask are considered cloudy
    For Scihub and PEPS, all pixels different to 4 or 5 in the SCL mask are considered cloudy
    The results are exported in a csv file, with the columns "area_name", "Date" and "cloudiness", containing the name of the Sentinel-2 tile, the date of acquisition, and the percentage of cloudy pixels.
 
    Parameters
    ----------
    sentinel_dir : str
        Path of the directory containing Sentinel-2 data.
    export_path : str
        Path to write csv file with extracted reflectance.
    tile_selection : list
        List of tiles from which to extract reflectance (ex : ["T31UFQ", "T31UGQ"]). If None, all tiles are extracted.
    sentinel_source : str
        Source of data, can be 'THEIA' et 'Scihub' et 'PEPS'
    """
    
    sentinel_dir = Path(sentinel_dir) ; export_path = Path(export_path)
    
    # if export_path.exists():
    #     extracted_cloudiness = pd.read_csv(export_path)
    
    list_dir = [x for x in sentinel_dir.iterdir() if x.is_dir()]
    cloudiness_list = []
    for directory in list_dir :
        if tile_selection is None or directory.stem in tile_selection:
        # for directory in list_dir :
            tile = TileInfo(directory)
            tile.getdict_datepaths("Sentinel",directory) #adds a dictionnary to tile.paths with key "Sentinel" and with value another dictionnary where keys are ordered and formatted dates and values are the paths to the directories containing the different bands
            tile.paths["Sentinel"] = get_band_paths(tile.paths["Sentinel"]) #Replaces the paths to the directories for each date with a dictionnary where keys are the bands, and values are their paths
            if len(tile.paths["Sentinel"]) >= 1:
                area_name = directory.stem
                # dict_example_raster[directory.stem] = list(list(tile.paths["Sentinel"].values())[0].values())[0]
                cloudiness = get_cloudiness(sentinel_dir / area_name / "cloudiness", tile.paths["Sentinel"], sentinel_source) #Returns dictionnary with cloud percentage for each date, except if lim_perc_cloud is set as 1, in which case cloud percentage is -1 for every date so source mask is not used and every date is used 
        
                area_cloudiness = pd.DataFrame.from_dict({"area_name" : len(cloudiness)*[area_name], "Date" : cloudiness.keys(), "cloudiness" : cloudiness.values()})
                cloudiness_list += [area_cloudiness]
    
    total_cloudiness = pd.concat(cloudiness_list, ignore_index=True)
    total_cloudiness.to_csv(export_path, mode='w', index=False, header=True)
    
    print("Cloudiness extracted from " + ', '.join(list(np.unique(total_cloudiness.area_name))))
    
if __name__ == '__main__':

        extract_cloudiness(
            sentinel_dir = "D:/fordead/fordead_data/sentinel_data/validation_tutorial/sentinel_data", 
            # sentinel_dir = "D:/fordead/Data/Sentinel",
            export_path = "D:/fordead/fordead_data/output/cloudiness_tuto.csv",
            tile_selection = ["T31UGP"])
        