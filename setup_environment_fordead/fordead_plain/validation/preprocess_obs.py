# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 10:52:12 2022

@author: rdutrieux
"""
import time
import click
import geopandas as gp
import pandas as pd
from pathlib import Path
from fordead.reflectance_extraction import attribute_id_to_obs, buffer_obs

@click.command(name='preprocess_obs')
@click.option("--obs_path", type = str,default = None, help = "Path of vector file containing observation points or polygons to preprocess", show_default=True)
@click.option("--export_path", type = str,default = None, help = "Path used to export the resulting preprocessed observation points or polygons", show_default=True)
@click.option("--buffer", type = int,default = None, help = "Length in meters of the buffer used to dilate (positive integer) or erode (negative integer) the observations. If None, no buffer is applied. Some observations may disappear completely if a negative buffer is applied", show_default=True)
@click.option("--name_column", type = str,default = "id", help = "Name of the column used to identify observations. If the column doesn't already exists, it is added as an integer between 1 and the number of observations", show_default=True)
def cli_preprocess_obs(obs_path, export_path, buffer, name_column):
    """
    Used as a preprocessing function for a vector file containing observation points or polygons. Can add an ID column if one does not already exist, and can also apply a buffer to erode or dilate observations.
    \f

    """

    preprocess_obs(**locals())
    
    
def preprocess_obs(obs_path, export_path, buffer = None, name_column = "id"):
    """
    Used as a preprocessing function for a vector file containing observation points or polygons. Can add an ID column if one does not already exist, and can also apply a buffer to erode or dilate observations.

    Parameters
    ----------
    obs_path : str
        Path of vector file containing observation points or polygons to preprocess.
    export_path : str
        Path used to export the resulting preprocessed observation points or polygons.
    buffer : int, optional
        Length in meters of the buffer used to dilate (positive integer) or erode (negative integer) the observations. If None, no buffer is applied. Some observations may disappear completely if a negative buffer is applied. The default is None.
    name_column : str, optional
        Name of the column used to identify observations. If the column doesn't already exists, it is added as an integer between 1 and the number of observations. The default is "id".

    """
    
    obs = gp.read_file(obs_path)
    obs = attribute_id_to_obs(obs, name_column)

    if buffer is not None:
        obs = buffer_obs(obs, buffer, name_column)
        
    if Path(export_path).exists():
        raise Exception(str(export_path) + " already exists.")
    obs.to_file(export_path)
    

if __name__ == '__main__':

    preprocess_obs(obs_path = "D:/fordead/Data/Test_programme/export_reflectance/vectors/Export_57_2022.shp",
                   export_path = "D:/fordead/Data/Test_programme/export_reflectance/reflectance_scolytes/preprocessed_Export_57_2022.shp", 
                   buffer = None, name_column = "id")