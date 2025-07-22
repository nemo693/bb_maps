import pandas as pd
from shapely.geometry import Polygon, Point
import rasterio
import geopandas as gp
from pathlib import Path
import numpy as np
from rasterio.crs import CRS

from fordead.import_data import TileInfo, get_band_paths, get_raster_metadata
from dask.diagnostics import ProgressBar

# import rasterio.sample
# =============================================================================
# PREPROCESS POLYGONS
# =============================================================================

def attribute_id_to_obs(obs, name_column):
    """
    Adds an ID column if it doesn't already exists. If column named after name_column parameter does not exist in the geodataframe, adds one with integers from 1 to the number of observations.

    Parameters
    ----------
    obs : geopandas GeoDataFrame
        Observation points or polygons
    name_column : str
        Name of the ID column.

    Returns
    -------
    obs : geopandas GeoDataFrame
        Observation points or polygons with added column named using parameter name_column if it doesn't already exist.

    """
    
    if name_column not in obs.columns:
        print("Creating " + name_column + " column")
        obs.insert(1, name_column, range(1,len(obs)+1))
    return obs

def buffer_obs(obs, buffer, name_column):
    """
    Applies a buffer to dilate or erode observations. Names of discarded observations too small compared to a eroding buffer are printed.

    Parameters
    ----------
    obs : geopandas GeoDataFrame
        Observation points or polygons with a column named name_column used to identify observations.
    buffer : int
        Length in meters of the buffer used to dilate (positive integer) or erode (negative integer) the observations. Some observations may disappear completely if a negative buffer is applied.
    name_column : str
        Name of the column used to identify observations.

    Returns
    -------
    obs : geopandas GeoDataFrame
        Observation polygons with the buffer applied.

    """
    
    obs['geometry']=obs.geometry.buffer(buffer)
    empty_obs = obs[obs.is_empty]
    if len(empty_obs) != 0:
        print(str(len(empty_obs)) + " observations were too small for the chosen buffer. \nIds are the following : \n" + ', '.join(map(str, empty_obs[name_column])) + " \n")
    obs=obs[~(obs.is_empty)]
    return obs


# =============================================================================
#   MAKE GRID POINTS
# =============================================================================

def get_grid_points(obs_polygons, sen_polygons, name_column):
    """
    Generates points in a grid corresponding to the centroids of Sentinel-2 pixels inside the polygons.
    

    Parameters
    ----------
    obs_polygons : geopandas GeoDataFrame
        Observation polygons with a column named name_column used to identify observations.
    sentinel_dir : str
        Path of directory containing Sentinel-2 data.
    name_column : str
        Name of the column used to identify observations.
    list_tiles : list
        A list of names of Sentinel-2 directories. If this parameter is used, extraction is  limited to those directories.

    Returns
    -------
    grid_points : geopandas GeoDataFrame
        Points corresponding to the centroids of the Sentinel-2 pixels of each Sentinel-2 tile intersecting with the polygons.

    """
    
    sen_obs_intersection = get_sen_obs_intersection(obs_polygons, sen_polygons, name_column)
    
    grid_points = polygons_to_grid_points(sen_obs_intersection, name_column)
    
    return grid_points
    

def get_sen_obs_intersection(obs_polygons, sen_polygons, name_column):
    """
    Observations polygons are intesected with Sentinel-2 tiles extent vector. Intersections where the observation polygon did not entirely fit in the sentinel-2 tile are removed.
    Observation polygons which intersect no Sentinel-2 tiles are removed and their IDs are printed.


    Parameters
    ----------
    obs_polygons : geopandas GeoDataFrame
        Polygons of observations.
    sen_polygons : geopandas GeoDataFrame
        Polygons of Sentinel-2 tiles extent with 'area_name' and 'epsg' columns corresponding to the name of the tile, and the projection system respectively.
    name_column : str
        Name of the ID column.

    Returns
    -------
    geopandas GeoDataFrame
        Intersection of obs_polygons and sen_polygons, with incomplete intersections removed.

    """
    # obs_polygons = obs_polygons.to_crs("2154")
    # obs_polygons = obs_polygons.to_crs("2154")
    # obs_polygons.crs
    sen_polygons = sen_polygons.to_crs(obs_polygons.crs)
    # obs_polygons = obs_polygons.to_crs(sen_polygons.crs)
    obs_area_tot = obs_polygons[[name_column]]
    # obs_area_tot = obs_polygons[name_column]
    obs_area_tot.insert(1, "area_tot", obs_polygons.area)

    obs_polygons = obs_polygons.merge(obs_area_tot, on= name_column, how='left')
    obs_intersection = gp.overlay(obs_polygons, sen_polygons)
    
    obs_intersection.insert(1, "area_intersect", obs_intersection.area)
    obs_intersection = obs_intersection[obs_intersection["area_tot"].round(3) == obs_intersection['area_intersect'].round(3)]
    
    unvalid_obs = obs_polygons[~obs_polygons[name_column].isin(obs_intersection[name_column])]
    if len(unvalid_obs) != 0:
        print("Observations not fitting entirely on one tile found.\nThey are removed from the dataset. \nIds are the following : \n" + ', '.join(map(str, unvalid_obs[name_column])) + " \n")
    
    # return obs_intersection[[name_column,"area_name","epsg"]]
    return obs_intersection.drop(columns = ["area_intersect","area_tot"])
    
    # return sen_obs_intersection

def get_polygons_from_sentinel_dirs(sentinel_dir, list_tiles):
    """
    Creates a Sentinel-2 tiles extent vector from existing Sentinel-2 data

    Parameters
    ----------
    sentinel_dir : str
        Path of directory containing Sentinel-2 data.
    list_tiles : list
        A list of names of Sentinel-2 directories. If this parameter is used, extraction is  limited to those directories.

    Returns
    -------
    concat_areas : geopandas GeoDataFrame
        Vector containing the extent of Sentinel-2 tiles contained in sentinel_dir directory, with 'epsg' and 'area_name' columns corresponding to the projection system and the name of the tile derived from the name of the directory containing its data.

    """
    
    sentinel_dir = Path(sentinel_dir)
    
    list_dir = [x for x in sentinel_dir.iterdir() if x.is_dir()]
    dict_example_raster = {}
    for directory in list_dir :
        if list_tiles is None or directory.stem in list_tiles:
        # for directory in list_dir :
            tile = TileInfo(directory)
            tile.getdict_datepaths("Sentinel",directory) #adds a dictionnary to tile.paths with key "Sentinel" and with value another dictionnary where keys are ordered and formatted dates and values are the paths to the directories containing the different bands
            tile.paths["Sentinel"] = get_band_paths(tile.paths["Sentinel"]) #Replaces the paths to the directories for each date with a dictionnary where keys are the bands, and values are their paths
            if len(tile.paths["Sentinel"]) >= 1:
                dict_example_raster[directory.stem] = list(list(tile.paths["Sentinel"].values())[0].values())[0]
    if len(dict_example_raster) == 0 :
        raise Exception("No Sentinel-2 data to extract")
        
    for area_index,area in enumerate(dict_example_raster):
        raster_metadata = get_raster_metadata(dict_example_raster[area])
        
        lon_point_list = [raster_metadata["extent"][0],raster_metadata["extent"][2],raster_metadata["extent"][2],raster_metadata["extent"][0]]
        lat_point_list = [raster_metadata["extent"][1],raster_metadata["extent"][1],raster_metadata["extent"][3],raster_metadata["extent"][3]]
        # lon_point_list = [raster_metadata["transform"][2], raster_metadata["transform"][2]+raster_metadata["sizes"]["x"]*10, raster_metadata["transform"][2]+raster_metadata["sizes"]["x"]*10, raster_metadata["transform"][2], raster_metadata["transform"][2]]
        # lat_point_list = [raster_metadata["transform"][5], raster_metadata["transform"][5], raster_metadata["transform"][5]-10*raster_metadata["sizes"]["y"], raster_metadata["transform"][5]-10*raster_metadata["sizes"]["y"],raster_metadata["transform"][5]]
        polygon_geom = Polygon(zip(lon_point_list, lat_point_list))
        polygon = gp.GeoDataFrame(index=[0], crs=raster_metadata["crs"], geometry=[polygon_geom])
        polygon.insert(1, "area_name", area)
        polygon.insert(2, "epsg", raster_metadata["crs"].to_epsg())
        
        if area_index == 0:
            concat_areas = polygon
        else:
            if concat_areas.crs != polygon.crs:
                polygon = polygon.to_crs(concat_areas.crs)
            concat_areas = pd.concat([concat_areas,polygon])
    # concat_areas.to_file("E:/fordead/Data/Vecteurs/aa_concat_areas.shp")
    return concat_areas





# =============================================================================


def polygons_to_grid_points(polygons, name_column):
    """
    Converts polygons to points corresponding to the centroids of the Sentinel-2 pixels of each Sentinel-2 tile intersecting with the polygons.
    Prints polygons with no pixels centroids inside of them.  
    
    Parameters
    ----------
    polygons : geopandas GeoDataFrame
        Observation polygons with 'area_name' and 'epsg' columns, corresponding to the name of a Sentinel-2 tile and its CRS respectively.
    name_column : str
        Name of the ID column.

    Returns
    -------
    grid_points : geopandas GeoDataFrame
        Points corresponding to the centroids of the Sentinel-2 pixels of each Sentinel-2 tile intersecting with the polygons.

    """
    
    # polygons = obs_polygons.merge(polygon_area_crs, on= name_column, how='inner') 
    
    grid_list = []
    for epsg in np.unique(polygons.epsg):
        print("epsg : " + str(epsg))
        
        polygons_epsg = polygons[polygons.epsg == epsg]
        polygons_epsg = polygons_epsg.to_crs(epsg = epsg)
        
        for area_name in np.unique(polygons.area_name):
            polygons_area = polygons_epsg[polygons_epsg.area_name == area_name]
        
            for polygon_index in range(polygons_area.shape[0]):
                polygon = polygons_area.iloc[polygon_index:(polygon_index+1)]
                
                bounds = get_bounds(polygon)
                
                obs_grid = gp.GeoDataFrame(
                    geometry=[
                        Point(x, y)
                        for x in np.arange(bounds[0], bounds[2], 10)
                        for y in np.arange(bounds[1], bounds[3], 10)
                    ],
                    crs=CRS.from_epsg(epsg),
                ).to_crs(CRS.from_epsg(epsg))
                
                obs_grid = gp.clip(obs_grid, polygon)
                
                obs_grid.insert(0,name_column,polygon[name_column].iloc[0])
                obs_grid.insert(1,"id_pixel",range(1,len(obs_grid)+1))
                obs_grid.insert(0,"area_name",area_name)
                obs_grid.insert(0,"epsg",epsg)
                obs_grid = obs_grid.to_crs(polygons.crs)
                
                if len(obs_grid) == 0:
                    print("Polygon " + str(polygon.iloc[0][name_column]) + " contains no pixel centroid")
                
                grid_list += [obs_grid]
                
        
    try:
        grid_points = gp.GeoDataFrame( pd.concat(grid_list, ignore_index=True), crs=grid_list[0].crs)
    except ValueError:
        raise Exception("It is probable that none of the observations intersect with the Sentinel-2 data available")
        
    return grid_points

def get_bounds(obs):
    """
    Get bounds of around of a polygons so it matches the limits of Sentinel-2 pixels
    
    Parameters
    ----------
    obs : geopandas GeoDataFrame
        A polygon

    Returns
    -------
    bounds : 1D array
        Bounds around the polygon

    """
    
    bounds = obs["geometry"].total_bounds 
    bounds[[0,1]] = bounds[[0,1]] - bounds[[0,1]] % 10 - 5
    bounds[[2,3]] = bounds[[2,3]] - bounds[[2,3]] % 10 + 15
    bounds = bounds.astype(int)
    return bounds


# =============================================================================
#   GET REFLECTANCE AT POINTS
# =============================================================================


def get_reflectance_at_points(grid_points,sentinel_dir, extracted_reflectance, name_column, bands_to_extract, tile_selection):
    reflectance_list = []
    
    if tile_selection is not None:
        grid_points = grid_points[grid_points["area_name"].isin(tile_selection)]
        if len(grid_points) == 0:
            raise Exception("No observations in selected tiles")
    
    for epsg in np.unique(grid_points.epsg):
        print("epsg : " + str(epsg))
        points_epsg = grid_points[grid_points.epsg == epsg]
        points_epsg = points_epsg.to_crs(epsg = epsg)
        
        for area_name in np.unique(points_epsg.area_name):
            print("area_name : " + area_name)
            points_area = points_epsg[points_epsg.area_name == area_name]
            
            if extracted_reflectance is not None:
                extracted_reflectance_area = extracted_reflectance[extracted_reflectance["area_name"] == area_name]
            else:
                extracted_reflectance_area = None
            
            raster_values = (points_area, sentinel_dir / area_name, extracted_reflectance_area, name_column, bands_to_extract)
            
            if raster_values is not None:
                reflectance_list += [raster_values]
    # reflectance = gp.GeoDataFrame(pd.concat(reflectance_list, ignore_index=True), crs=reflectance_list[0].crs)
    if len(reflectance_list) != 0:
        reflectance = pd.concat(reflectance_list, ignore_index=True)
        if "Mask" in bands_to_extract:
            reflectance["Mask"] = reflectance["Mask"].astype(np.uint8)
    else:
        reflectance = None
        
    return reflectance


def extract_points(x, df, **kwargs):
    """_summary_

    Parameters
    ----------
    x : xarray.DataArray or xarray.Dataset
    df : pandas.DataFrame
        Coordinates of the points

    Returns
    -------
    xarray.DataArray or xarray.Dataset
        The points values.

    Examples
    --------
    >>> import xarray as xr
    >>> import pandas as pd
    >>> import dask
    >>> da = xr.DataArray(
    ... # np.random.random((100,200)),
    ... dask.array.random.random((100,200), chunks=10),
    ... coords = [('x', range(100)), ('y', range(200))]
    ... )
    >>> df = pd.DataFrame(
    ...    dict(
    ...        x=np.random.permutation(range(100))[:100]+np.random.random(100),
    ...        y=np.random.permutation(range(100))[:100]+np.random.random(100),
    ...        z=range(100),
    ...    )
    ... )
    >>> extract_points(da, df, method="nearest", tolerance=.5)

    """
    # x = da
    xk = x.dims
    coords_cols = [c for c in df.keys() if c in xk]
    coords = df[coords_cols]
    points = x.sel(coords.to_xarray(), **kwargs)
    return points

def extract_raster_values(points, tile_coll, bands_to_extract, extracted_reflectance, export_path=None):
    """
    Sample raster values for each XY points

    Parameters
    ----------
    points : geodataframe
        Observation points
    tile_coll : pystac.ItemCollection
        Item_collection of a unique MGRS Tile
    bands_to_extract : list of strings
        List of bands to extract
    extracted_reflectance: pandas.DataFrame
        Table of already sampled raster values, not to extract again.
        Expected columns are "area_name", "Date" and an ID columns to merge with points dataframe.
    export_path : str, optional
        Path to export the result
        
    Returns
    -------
    pandas DataFrame or None
        The result is written in file `export_path`.
    """
    if extracted_reflectance is not None:
        er_keys = extracted_reflectance.keys()
        if "area_name" not in er_keys:
            raise ValueError("area_name not in extracted_reflectance")
        if "Date" not in er_keys:
            raise ValueError("Date not in extracted_reflectance")


    arr = tile_coll.to_xarray(xy_coords="center", assets=bands_to_extract).rename("value")
    if not points.crs.equals(arr.rio.crs):
        points = points.to_crs(arr.rio.crs)
    points.index.rename("id_point", inplace=True)
    coords = points.get_coordinates().join(points.drop(columns='geometry'))

    # To respect original implementation: extract only the points not already extracted
    # not sure it is faster than extracting directly the all the points...
    if extracted_reflectance is None or extracted_reflectance.empty:
        p = extract_points(arr, coords, method="nearest", tolerance=arr.rio.resolution()[0]/2)
    else:
        p_temp = extract_points(arr, coords, method="nearest", tolerance=arr.rio.resolution()[0]/2)
        # copy coords for each dates
        coords2 = coords.reset_index(drop=False).merge(pd.Series(arr.time.values, name="time"), how="cross")
        coords2["Date"] = coords2.time.dt.date.astype(str)
        # keep only the points and dates not in extracted_reflectance
        to_extract = coords2.merge(extracted_reflectance, on=extracted_reflectance.keys(), how="outer", indicator=True)
        to_extract = to_extract.query("_merge == 'left_only'").drop("_merge", axis=1)
        if to_extract.empty:
            # nothing to extract
            return
        # subset only the points to extract
        p = extract_points(p_temp, to_extract)
    with ProgressBar():
        p = p.to_dataframe()
    p.reset_index(drop=False, inplace=True)

    # it may have a problem if two images at the same time,
    # in that case `id` should be added to the index arg,
    # and a solution should be found to choose between these values
    index = ["id_point", "time"]
    p = p.pivot(columns="band", values="value", index=index)

    # drop rows with NA values
    if p.isna().any().any():
        from warnings import warn
        warn("Found NA values, dropping them...")
        p.dropna(inplace=True)

    # change type to int
    extractions = p.astype("int").reset_index(drop=False)
    # join extractions with points
    points = points.drop(columns='geometry').reset_index(drop=False)
    # reformat result
    res = points.merge(extractions, on=["id_point"])
    res.time = res.time.dt.strftime("%Y-%m-%d")
    res = res.rename(columns={"time": "Date"}).drop(columns="id_point")

    if export_path is None:
        return res
    res.to_csv(export_path, mode='a', index=False, header=not(export_path.exists()))
    # epsg,area_name,id,id_pixel,Date,B2,B3,B4,B5,B6,B7,B8,B8A,B11,B12,Mask
    # 32631,T31UGP,0,1,2019-02-17,118,172,129,267,797,986,1224,1265,427,203,0


# ==========================================================================================================================================================

# =============================================================================
#   process_points
# =============================================================================

def process_points(points, sen_polygons, name_column):
    """
    Intersects observation points with Sentinel-2 tiles extent vector.
    Adds an 'id_pixel' column filled with 0 so the resulting vector can be used in export_reflectance function.
    Points outside of available Sentinel-2 tiles are detected and their IDs are printed.
    
    Parameters
    ----------
    points : geopandas GeoDataFrame
        Observation points used for intersection
    sen_polygons : geopandas GeoDataFrame
        Polygons of Sentinel-2 tiles extent with 'area_name' and 'epsg' columns corresponding to the name of the tile, and the projection system respectively.
    name_column : str
        Name of the ID column in points

    Returns
    -------
    obs_intersection : geopandas GeoDataFrame
        Intersection of points and sen_polygons, with added 'id_pixel' columns

    """
    
    sen_polygons = sen_polygons.to_crs(points.crs)
    obs_intersection = gp.overlay(points, sen_polygons)
    # obs_intersection["id_pixel"] = 0
    obs_intersection.insert(1,"id_pixel",1) #Insert column
    
    outside_points = points[~points[name_column].isin(obs_intersection[name_column])]
    
    if len(outside_points) != 0:
        print(str(len(outside_points)) + " point observations were outside available Sentinel-2 tiles.\nThey are removed from the dataset. \nIds are the following : \n" + ', '.join(map(str, outside_points[name_column])) + " \n")

    if len(obs_intersection) == 0 :
        raise Exception("No available Sentinel-2 tiles at point locations")
    
    return obs_intersection[["epsg","area_name",name_column,"id_pixel","geometry"]]
    
    # return sen_intersection_points
    

# def get_sen_intersection_points(points, sen_polygons, name_column):

    

def get_already_extracted(export_path, name_column, bands_to_extract):
    """
    
    Returns already extracted acquisition dates for each observation and each tile.

    Parameters
    ----------
    export_path : str
        Path to reflectance.csv
    name_column : str
        Name of the ID column in obs
    bands_to_extract : list
        List of bands to extract

    Returns
    -------
    pandas DataFrame
        Already extracted acquisition dates for each observation and each tile
        Columns are "area_name", name_column, "Date".

    """
    
    if not export_path.exists():
        return
    reflectance = pd.read_csv(export_path)
    if name_column not in reflectance.columns:
        # raise Exception("name_column '"+ name_column + "' not in reflectance.csv found in " + str(export_dir))
        raise Exception("name_column '"+ name_column + "' not in " + str(export_path))
    if not all([b in reflectance.keys() for b in bands_to_extract]):
        raise ValueError("Some bands to extract are not in reflectance.csv found in " + str(export_path))

    extracted_reflectance = reflectance[["area_name", name_column,"Date"]].drop_duplicates()
    
    # obs_extracted = np.unique(extracted_reflectance[name_column])
    # obs_to_extract = np.unique(obs[name_column])
    
    # diff_nb_obs = len(obs_extracted) != len(obs_to_extract)
    # missing_obs = not(all(item in obs_to_extract for item in obs_extracted))
    # added_obs = not(all(item in obs_extracted for item in obs_to_extract))
    
    # if diff_nb_obs or missing_obs or added_obs:
    #     print("Changes to "+ str(obs_path) + " have been detected since last extracting reflectances to " + str(export_path))
    #     if diff_nb_obs:
    #         print("Different number of observations detected")
    #     if missing_obs:
    #         print("Some observations are missing, already extracted reflectance will stay untouched")
    #     if added_obs:
    #         print("Some observations were added, their reflectance will be extracted")
    #     answer = input("Do you want to continue ? (yes/no)")
    #     if answer != "yes":
    #         raise Exception("Extraction of reflectance stopped")

    return extracted_reflectance    
    

# def _vrt(points,sentinel_dir):
#     """Must have the same crs"""

# from osgeo import gdal

#     tile = TileInfo(sentinel_dir)
#     tile.getdict_datepaths("Sentinel",sentinel_dir) #adds a dictionnary to tile.paths with key "Sentinel" and with value another dictionnary where keys are ordered and formatted dates and values are the paths to the directories containing the different bands
#     tile.paths["Sentinel"] = get_band_paths(tile.paths["Sentinel"]) #Replaces the paths to the directories for each date with a dictionnary where keys are the bands, and values are their paths
    
#     coord_list = [(x,y) for x,y in zip(points['geometry'].x , points['geometry'].y)]
#     dates = list(tile.paths["Sentinel"])
#     bands = list(tile.paths["Sentinel"][list(tile.paths["Sentinel"])[0]])

#     dict_data = {"id_obs" : [id_obs for id_obs in points.id_obs for date in dates],
#                  "id_pixel" : [id_pixel for id_pixel in points.id_pixel for date in dates],
#                  "Date" : [date for id_obs in points.id_obs for date in dates]}
#     for band in bands:
#         print(band)
#         dict_data[band] = []
#         path_list = [str(tile.paths["Sentinel"][date][band]) for date in dates]
#         gdal.BuildVRT(str(sentinel_dir / "vrt.vrt"), path_list, separate=True)
#         with rasterio.open(str(sentinel_dir / "vrt.vrt")) as raster:
#             reflect_band = raster.sample(coord_list)
#             # dict_data[band] = list(reflect_band)
#             # list(reflect_band)
#             for x in reflect_band:
#                 dict_data[band] += list(x)

#     reflectance = pd.DataFrame(data=dict_data)
#     return reflectance
    

# def _chunks(points,sentinel_dir):
#     """Must have the same crs"""
#     tile = TileInfo(sentinel_dir)
#     tile.getdict_datepaths("Sentinel",sentinel_dir) #adds a dictionnary to tile.paths with key "Sentinel" and with value another dictionnary where keys are ordered and formatted dates and values are the paths to the directories containing the different bands
#     tile.paths["Sentinel"] = get_band_paths(tile.paths["Sentinel"]) #Replaces the paths to the directories for each date with a dictionnary where keys are the bands, and values are their paths
    
#     coord_list = [(x,y) for x,y in zip(points['geometry'].x , points['geometry'].y)]

    # dates = list(tile.paths["Sentinel"])
#     bands = list(tile.paths["Sentinel"][list(tile.paths["Sentinel"])[0]])
    
    
#     list_bands=[]
#     for band in bands:
#         print(band)
#         list_bands += [xr.open_mfdataset(path_list,concat_dim = "Time", combine = "nested", chunks = 1280, parallel=True).assign_coords(Time=dates).band_data]
#         # stack_bands=stack_bands.assign_coords(band=bands)
#     stack_bands=xr.concat(list_bands,dim="band").assign_coords(band=bands)
#     # stack_bands = stack_bands.chunk(1280)  
    
#     list_point = [stack_bands.sel(x=point[0], y=point[1], method="nearest") for point in coord_list]
#     data=xr.concat(list_point,dim="id")
#     data.to_dataframe()

