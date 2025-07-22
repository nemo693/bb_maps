# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 10:52:12 2022

@author: rdutrieux
"""

import click
import geopandas as gp
import pandas as pd
from pathlib import Path
from fordead.reflectance_extraction import get_grid_points, process_points, get_polygons_from_sentinel_dirs
from fordead.stac.stac_module import getItemCollection, get_bbox, get_harmonized_planetary_collection, get_polygons_from_sentinel_planetComp

@click.command(name='obs_to_s2_grid')
@click.option("--obs_path", type = str,default = None, help = "Path to a vector file containing observation points or polygons, must have an ID column corresponding to name_column parameter", show_default=True)
@click.option("--sentinel_source", type = str,default = None, help = "Can be either 'Planetary', in which case the Sentinel-2 grid is infered from Microsoft Planetary Computer stac catalogs, or the path of the directory containing Sentinel-2 data.", show_default=True)
@click.option("--export_path", type = str,default = None, help = "Path used to write resulting vector file, with added 'epsg','area_name' and 'id_pixel' columns", show_default=True)
@click.option("--name_column", type = str,default = "id", help = "Name of the ID column", show_default=True)
@click.option("-t","--tile_selection", type = list, default = None, help = "A list of names of Sentinel-2 directories. (ex : -t T31UFQ -t T31UGQ). If None, all tiles are used.", show_default=True)
@click.option("--overwrite",  is_flag=True, help = "Overwrites file at obs_path", show_default=True)
def cli_obs_to_s2_grid(obs_path, sentinel_source, export_path, name_column, tile_selection, overwrite):
    """
    Attributes intersecting Sentinel-2 tiles to observation points or polygons, adding their epsg and name. If polygons are used, they are converted to grid points located at the centroid of Sentinel-2 pixels.
    If points or polygons intersect several Sentinel-2 tiles, they are duplicated for each of them.
    If some intersect no Sentinel-2 tiles, they are removed and their IDs are printed.

    \f

    """
    
    obs_to_s2_grid(**locals())

def obs_to_s2_grid(obs_path, sentinel_source, export_path, name_column = "id", tile_selection = None, overwrite = False):
    """
    Attributes intersecting Sentinel-2 tiles to observation points or polygons, adding their epsg and name. If polygons are used, they are converted to grid points located at the centroid of Sentinel-2 pixels.
    If points or polygons intersect several Sentinel-2 tiles, they are duplicated for each of them.
    If some intersect no Sentinel-2 tiles, they are removed and their IDs are printed.

    Parameters
    ----------
    obs_path : str
        Path to a vector file containing observation points or polygons, must have an ID column corresponding to name_column parameter.
    sentinel_source : str
        Can be either 'Planetary', in which case the Sentinel-2 grid is infered from Microsoft Planetary Computer stac catalogs, or the path of the directory containing Sentinel-2 data.
    export_path : str
        Path used to write resulting vector file, with added "epsg","area_name" and "id_pixel" columns.
    name_column : str, optional
        Name of the ID column. The default is "id".
    tile_selection : list
        A list of names of Sentinel-2 directories. If this parameter is used, extraction is limited to those directories.
    overwrite : bool
        If True, allows overwriting of file at obs_path

    """
    
    # sentinel_dir = Path(sentinel_dir) 
    export_path = Path(export_path)
    # if sentinel_source == "Planetary":
    #     collection = get_harmonized_planetary_collection(sentinel_source, start_date, end_date, get_bbox(tile_obs), lim_perc_cloud, tile)
    #     get_polygons_from_sentinel_planetComp(item_collection)
        
    obs = gp.read_file(obs_path)

    if export_path.exists():
        print(str(export_path) + " already exists")
        if overwrite:
            print("It will be overwritten\n")
        else:
            raise Exception("Set 'overwrite' parameter as True to overwrite")
    
    geom_type = obs.geom_type.drop(columns = "geometry")
    points = obs[(geom_type == 'Point') | (geom_type == 'MultiPoint')]
    polygons = obs[(geom_type == 'Polygon') | (geom_type == 'MultiPolygon')]
    
    if sentinel_source == "Planetary":
        collection = get_harmonized_planetary_collection("2015-01-01", "2024-01-01", get_bbox(obs), 0.01)
        sen_polygons = get_polygons_from_sentinel_planetComp(collection, tile_selection)
    else:
        sen_polygons = get_polygons_from_sentinel_dirs(sentinel_source, tile_selection)
    
    # sen_polygons.to_file("D:/fordead/05_SUBPROJECTS/03_stac/03_RESULTS/sen_polygons.shp")
    points_from_points = process_points(points, sen_polygons, name_column) if len(points) != 0 else None
    points_from_poly = get_grid_points(polygons, sen_polygons, name_column) if len(polygons) != 0 else None

    total_points = pd.concat([points_from_poly,points_from_points])
    total_points.to_file(export_path)

    
if __name__ == '__main__':
        
    #Planetary
    obs_to_s2_grid(
        obs_path = "D:/fordead/fordead_data/vector/observations_tuto.shp",
        sentinel_source = "Planetary", 
        export_path = "D:/fordead/05_SUBPROJECTS/03_stac/03_RESULTS/pp_observations_tuto_plan.shp",
        name_column = "id",
        # tile_selection = ["T32ULU"],
        overwrite = True)
    
    #Theia locally
    # obs_to_s2_grid(
    #     obs_path = "D:/fordead/fordead_data/vector/observations_tuto.shp",
    #     sentinel_source = "D:/fordead/fordead_data/sentinel_data/validation_tutorial/sentinel_data", 
    #     export_path = "D:/fordead/05_SUBPROJECTS/03_stac/03_RESULTS/pp_observations_tuto_theia.shp",
    #     name_column = "id",
    #     overwrite = True)