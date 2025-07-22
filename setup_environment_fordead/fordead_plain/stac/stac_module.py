# author(s) : Kenji Ose / Raphael Dutrieux

import geopandas as gp
import pandas as pd
import numpy as np
# from shapely.geometry import Polygon, Point
from pyproj import Proj, transform
import pystac_client
import pystac
import planetary_computer
from pathlib import Path
from fordead.stac.theia_collection import build_theia_collection, ItemCollection
from datetime import datetime
import fordead.stac.stac_module as st
from shapely.geometry import Polygon
from pyproj import Transformer

def get_vectorBbox(vector_path):
    """
    get layer bounding box

    - vector_path: <string> vectorfile path

    return: <list> of bbox coordinates in lat/lon 
            [lon_min, lat_min, lon_max, lat_max]
    """
    obs = gp.read_file(vector_path)
    obs = obs.to_crs("EPSG:4326")
    bbox = list(obs.total_bounds)
    return bbox

def get_bbox(vector):
    """
    get layer bounding box

    - vector: <geopanda dataframe> vector

    return: <list> of bbox coordinates in lat/lon 
            [lon_min, lat_min, lon_max, lat_max]
    """
    vector = vector.to_crs("EPSG:4326")
    bbox = list(vector.total_bounds)
    return bbox


def getItemCollection(startdate, enddate, bbox, cloud_nb = 100):
    """
    get item collection from Microsoft Planetary Computer

    - startdate: <string> start date yyyy-mm-dd
    - enddate: <string> end date yyyy-mm-dd
    - bbox: <list> bounding box [lon_min, lat_min, lon_max, lat_max]
    - cloud_nb: <int> cloud percentage [default: 0]

    return pystac itemcollection object
    """
    catalog = pystac_client.Client.open(
    "https://planetarycomputer.microsoft.com/api/stac/v1",
    modifier=planetary_computer.sign_inplace,
    )
    time_range = f"{startdate}/{enddate}" #"2021-07-01/2021-12-01"
    bbox = bbox
    cloud_nb = cloud_nb

    # recherche avec filtre de type 'query'
    search = catalog.search(
        collections=["sentinel-2-l2a"],
        bbox=bbox,
        datetime=time_range,
        query={"eo:cloud_cover": {"lt": float(cloud_nb)}},
        sortby="datetime",
        # seems to solve the duplicates issue if less than 1000 items
        # see https://github.com/microsoft/PlanetaryComputer/issues/163
        limit=1000
    )

    # itemsColl = search.pages()
    itemsColl = search.item_collection()
    return ItemCollection(itemsColl)

# Conversion de tous les items découverts par la recherche, 
# on crée un objet ItemcCollection pour avoir la liste des items. 
# from sys import sys.exit
def get_s2_epsgXY(item_collection):
    """
    Using an item collection, returns a list with all tilenames, the epsg of the corresponding sentinel-2 data, as well as the list of longitude and latitude coordinates of the tile in epsg 4326.
    """
    
        
    corner_names = []
    x_coordinates = []
    y_coordinates = []
    s2_mgrs_tiles = []
    epsg_list = []
    for item in item_collection:
        # print(item.properties["eo:cloud_cover"])
        # Get properties from the item
        s2_mgrs_tile = "T" + item.properties.get('s2:mgrs_tile')
        epsg = item.properties["proj:epsg"]
        # print(epsg)
        # Extract the coordinates
        coordinates = item.geometry['coordinates'][0]
        
        for i, (x, y) in enumerate(coordinates):
            if len(coordinates) == 5:
                corner_name = f'Corner{i+1}'  # Create corner names (Corner1, Corner2, etc.)
                corner_names.append(corner_name)
                x_coordinates.append(x)
                y_coordinates.append(y)
                s2_mgrs_tiles.append(s2_mgrs_tile)
                epsg_list.append(epsg)
    
    # Create a Pandas DataFrame
    data = {
        'tilename': s2_mgrs_tiles,
        'epsg' : epsg_list,
        'Corner': corner_names,
        'x': x_coordinates,
        'y': y_coordinates
    }
    df = pd.DataFrame(data)
    
    # Calculate the min and max for each corner of each tile
    grouped = df.groupby(['tilename', 'Corner',"epsg"])
    min_max = grouped.agg({'x': ['min', 'max'], 'y': ['min', 'max']})
    min_max = grouped.agg({'x': ['min', 'max'], 'y': ['min', 'max']})

    min_max.reset_index(inplace=True)
    min_max.columns = ['tilename', 'Corner', "epsg",'xmin', 'xmax', 'ymin', 'ymax']
    
    out_list = []
    for tile in np.unique(min_max.tilename):
        # print(tile)
        # for corner in ["Corner1","Corner2","Corner3","Corner4","Corner5"]:
        corner1 = min_max[(min_max["tilename"] == tile) & (min_max["Corner"] == "Corner1")]
        corner2 = min_max[(min_max["tilename"] == tile) & (min_max["Corner"] == "Corner2")]
        corner3 = min_max[(min_max["tilename"] == tile) & (min_max["Corner"] == "Corner3")]
        corner4 = min_max[(min_max["tilename"] == tile) & (min_max["Corner"] == "Corner4")]
        # corner5 = min_max[(min_max["tilename"] == tile) & (min_max["Corner"] == "Corner5")]
        
        lon_point_list = [corner1.xmin.iloc[0],corner2.xmin.iloc[0],corner3.xmax.iloc[0],corner4.xmax.iloc[0]]
        lat_point_list = [corner1.ymax.iloc[0],corner2.ymin.iloc[0],corner3.ymin.iloc[0],corner4.ymax.iloc[0]]
        tileinfo = [tile, corner1.epsg.iloc[0], lon_point_list, lat_point_list]
        out_list.append(tileinfo)
    return out_list


def tile_filter(item_collection, tile_name):
    """
    Filter S2 item_collection according to Tile ID

    - item_collection: <pystac object> input S2 item_collection
    - tile_name: <string> S2 tile ID 

    return: <pystac object> filtered item_collection
    """
    filt_item = []
    for item in item_collection:
        if item.properties["s2:mgrs_tile"] == tile_name[1:]:
            filt_item.append(item)
    new_coll = pystac.item_collection.ItemCollection(filt_item)
    return new_coll



def get_items_tiles(item_collection):
    """
    Link S2 tile with corresponding EPSG

    - item_collection: <pystac object> S2 item collection

    return: <list> of list with following elements [tile name, epsg code]
    """
    out_list = []
    for item in item_collection:
        tile = item.properties['s2:mgrs_tile']
        epsg = item.properties["proj:epsg"]
        info = [tile, epsg]
        if info not in out_list:
            out_list.append(info)
    return out_list



    
S2_THEIA_BANDS = [f"B{i+1}" for i in range(12)]+["B8A"]
S2_SEN2COR_BANDS = [f"B{i+1:02}" for i in range(12)]+["B8A"]
def get_harmonized_planetary_collection(start_date, end_date, obs_bbox, lim_perc_cloud, tile = None):
    """
    Get the planetary item collection with band names
    and offset harmonized with THEIA format.

    Parameters
    ----------
    start_date : str
        Start of the period.
    end_date : str
        End of the period.
    obs_bbox : list
        bounding box [lon_min, lat_min, lon_max, lat_max]
    lim_perc_cloud : float
        Max cloud cover between 0 et 1.
    tile : str, optional
        The name of a single tile used to filter the collection. If None, all tiles are kept in the collection. The default is None.

    Returns
    -------
    fordead.stac.stac_module.ItemCollection
        The item collection.
    
    Notes
    -----
    The offset is harmonized with `harmonize_sen2cor_offet`.
    """
    
    corresp_keys = {'B01' : 'B1', 'B02' : 'B2', 'B03' : "B3", 'B04' : "B4", 'B05' : "B5", 
                    'B06' : "B6", 'B07' : "B7", 'B08' : "B8", 'B09' : "B9", 'SCL' : "Mask"}
    
    collection = getItemCollection(start_date, end_date, obs_bbox, int(lim_perc_cloud*100))
    # item_list = []
    for item in collection:
        # item.properties["cloud_cover"] = item.properties['eo:cloud_cover']
        item.properties["tilename"] = "T" + item.properties['s2:mgrs_tile']
        #Changing B01 to B1, up to B9, and SCL to Mask in item asset names

        for key_change in corresp_keys:
            item.assets[corresp_keys[key_change]] = item.assets.pop(key_change)
            
    # for item in collection:
    #     print(item.properties["cloud_cover"])
    if tile is not None:
        print("Extracting Sentinel-2 data from "+ tile + " collected from Microsoft Planetary Computer")
        collection = collection.filter(filter=f"tilename = '{tile}'")
    
    harmonize_sen2cor_offet(collection, bands=S2_THEIA_BANDS, inplace=True)
    collection.drop_duplicates(inplace=True)

    return collection

def harmonize_sen2cor_offet(collection, bands=set(S2_THEIA_BANDS + S2_SEN2COR_BANDS), inplace=False):
    if not inplace:
        collection = collection.copy()
    for item in collection:
        for asset in bands:
            if asset in item.assets:
                if item.properties["datetime"] >= "2022-01-25":
                    item.assets[asset].extra_fields["raster:bands"] = [dict(offset=-1000)]
                else:
                    item.assets[asset].extra_fields["raster:bands"] = [dict(offset=0)]
    if inplace:
        return collection

def get_harmonized_theia_collection(sentinel_source, tile_cloudiness, start_date, end_date, lim_perc_cloud, tile):
    # lim_cloud_cover = int(lim_perc_cloud*100)
    corresp_keys = {'CLM' : "Mask"}
    
    sentinel_source = Path(sentinel_source)
    sen_dir_path =  sentinel_source / tile
    
    if sentinel_source.is_dir():
        if sen_dir_path.is_dir():
            print("Extracting Sentinel-2 data from " + tile + " stored locally")

            collection = build_theia_collection(sen_dir_path)
            
            
            for item in collection:
                
                #Adding cloud_cover to theia_collection 
                if tile_cloudiness is not None:
                    try:
                        date_datetime = datetime.strptime(item.properties["datetime"], "%Y-%m-%dT%H:%M:%S.%fZ")
                    except ValueError:
                        date_datetime = datetime.strptime(item.properties["datetime"], "%Y-%m-%dT%H:%M:%SZ")
                    try:
                        acqu_cloudiness = tile_cloudiness[tile_cloudiness.Date == date_datetime.strftime("%Y-%m-%d")]["cloudiness"].values[0]
                    except Exception as e:
                        # traceback_str = traceback.format_exc()
                        print(f"Error: {e}\nCloudiness does not seem to have been extracted for this tile")
                    item.properties["cloud_cover"] = acqu_cloudiness
                
                #Changing 'CLM' to 'Mask' 
                for key_change in corresp_keys:
                    item.assets[corresp_keys[key_change]] = item.assets.pop(key_change) #Might be changed in the future

            collection = collection.filter(datetime=f"{start_date}/{end_date}")
            
            if tile_cloudiness is not None : collection = collection.filter(filter = f'cloud_cover < {lim_perc_cloud}') 
            
            
            if len(collection) == 0:
                print("Collection is empty")
            
        else:
            raise Exception("No directory found at " + str(sen_dir_path))
    else:
        raise Exception("Unrecognized sentinel_source")
    
    return collection
    
def get_polygons_from_sentinel_planetComp(item_collection, tile_selection = None, outputfile = None):
    """
    get image footprint with specified tile crs
    
    - item_collection: pystac item_collection object

    return: geopandas dataframe
    """    
    
    # for item in item_collection:
    #     print(item.properties)
    
    tilesinfo = st.get_s2_epsgXY(item_collection)
    # tilesinfo = st.coord2epsg(tile_extents_4326, in_epsg = "4326", out_epsg = str(epsg))   
    
    if tile_selection is not None:
        tilesinfo = [tile for tile in tilesinfo if tile[0] in tile_selection]
        
    # i=0
    for i in range(len(tilesinfo)):
        # print(i)
        # tile.xmin
        
        # lon_point_list = [tile[2][0][0],tile[2][1][0],tile[2][1][0],tile[2][0][0]]
        # lon_point_list = [tile.xmin,tile.xmax,tile.xmax,tile.xmin]
        # lat_point_list = [tile.ymin,tile.ymin,tile.ymax,tile.ymax]
        #polygon_geom = gp.points_from_xy(lon_point_list, lat_point_list)
        polygon_geom = Polygon(zip(tilesinfo[i][2], tilesinfo[i][3]))
        polygon = gp.GeoDataFrame(index=[0], crs='epsg:4326', geometry=[polygon_geom])
        polygon.insert(1, "area_name", tilesinfo[i][0])
        polygon.insert(2, "epsg", f'{tilesinfo[i][1]}')
        
        if i == 0:
            concat_areas = polygon
        else:
            # if concat_areas.crs != polygon.crs:
            #         polygon = polygon.to_crs(concat_areas.crs)
            concat_areas = pd.concat([concat_areas,polygon])
        i +=1
        
    
    if outputfile != None:
        concat_areas.to_file(outputfile) # a sécuriser
    
    return concat_areas


if __name__ == '__main__':
    obs_path = "D:/PROJETS/PRJ_FORDEAD/TEST_STAC/data/observations_tuto.shp"
    obs_bbox = get_vectorBbox(obs_path)
    print(obs_bbox)

    coll = getItemCollection("2021-07-01", "2021-12-01", obs_bbox, 20)
    print(coll.items)
    
