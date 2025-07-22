# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 10:30:53 2020

@author: Raphael Dutrieux


Crée un timelapse à partir des résultats calculés


"""

#%% =============================================================================
#   LIBRAIRIES
# =============================================================================

# from glob import glob
import geopandas as gp
import click
from plotly.offline import plot
from zipfile import ZipFile, ZIP_DEFLATED

# %%=============================================================================
#  IMPORT LIBRAIRIES PERSO
# =============================================================================

from fordead.results_visualisation import CreateTimelapse, intersects_area
from fordead.import_data import TileInfo

@click.group()
def timelapse():
    """Main entrypoint."""
    
@timelapse.command(name='timelapse')
@click.option("-o", "--data_directory", type = str, help = "Path of the directory containing results from the region of interest")
@click.option("--shape_path", type = str, help = "Path of the shapefile of the area, or points, to convert to timelapse. Not used if timelapse made from x and y coordinates")
@click.option("--name_column", type = str, default = "id", help = "Name of the column containing the name of the export. Not used if timelapse made from x and y coordinates", show_default=True)
@click.option("-x", "--x", type = int, help = "Coordinate x in the crs of the Sentinel-2 tile. Not used if timelapse is made using a shapefile")
@click.option("-y", "--y", type = int, help = "Coordinate y in the crs of the Sentinel-2 tile. Not used if timelapse is made using a shapefile")
@click.option("--buffer", type = int, default = 100, help = "Buffer around polygons or points for the extent of the timelapse", show_default=True)
@click.option("--vector_display_path", type = str, help = "Path of the shapefile with ground observations")
@click.option("-d","--hover_column_list", multiple=True, type=str, default=None, help="List of columns to display when hovering mouse over vectors of vector_display_path", show_default=True)
@click.option("--max_date", type = str, default = None, help = "Last date used in the timelapse")
@click.option("--show_confidence_class",  is_flag=True, help = "If True, detected dieback is shown with the confidence class at the last date used, as vectorized in the step [05_export_results](https://fordead.gitlab.io/fordead_package/docs/user_guides/english/05_export_results/)", show_default=True)
@click.option("--zip_results",  is_flag=True, help = "If True, puts timelapses in a zip file", show_default=True)
def cli_create_timelapse(data_directory, shape_path = None, name_column = "id", x = None, y = None, buffer = 100, vector_display_path = None, hover_column_list = None, max_date = None, show_confidence_class = False, zip_results = False):
    """
    Create timelapse allowing navigation through Sentinel-2 dates with detection results superimposed.
    By specifying 'shape_path' and 'name_column' parameters, it can be used with a shapefile containing one or multiple polygons or points with a column containing a unique ID used to name the export. 
    By specifying 'x' and 'y' parameters, it can be used by specifying coordinates in the system of projection of the tile. 
    The timelapse is exported in the data_directory/Timelapses directory as an html file.
    See details https://fordead.gitlab.io/fordead_package/docs/user_guides/Results_visualization/
    \f
    Parameters
    ----------
    data_directory
    shape_path
    name_column
    x
    y
    buffer
    vector_display_path
    hover_column_list
    max_date
    show_confidence_class
    zip_results


    """
    create_timelapse(data_directory, shape_path, name_column, x, y, buffer, vector_display_path, hover_column_list, max_date, show_confidence_class, zip_results)

#%% =============================================================================
#   MAIN CODE
# =============================================================================
       


def create_timelapse(data_directory, shape_path = None, name_column = "id",  x = None, y = None, buffer = 100, vector_display_path = None, hover_column_list = None, max_date = None, show_confidence_class = False, zip_results = False):
    """
    Create timelapse allowing navigation through Sentinel-2 dates with detection results superimposed.
    By specifying 'shape_path' and 'name_column' parameters, it can be used with a shapefile containing one or multiple polygons with a column containing a unique ID used to name the export. 
    By specifying 'x' and 'y' parameters, it can be used by specifying coordinates in the system of projection of the tile. 
    The timelapse is exported in the data_directory/Timelapses directory as an html file.
    See details https://fordead.gitlab.io/fordead_package/docs/user_guides/Results_visualization/

    Parameters
    ----------
    data_directory : str
        Path of the directory containing results from the region of interest.
    shape_path : str, optional
        Path of the shapefile of the area or points to convert to timelapse. Not used if timelapse made from x and y coordinates. The default is None.
    name_column : str, optional
        Name of the column containing the name of the export. Not used if timelapse made from x and y coordinates. The default is "id".
    x : int, optional
        Coordinate x in the crs of the Sentinel-2 tile. Not used if timelapse is made using a shapefile. The default is None.
    y : int, optional
        Coordinate y in the crs of the Sentinel-2 tile. Not used if timelapse is made using a shapefile. The default is None.
    buffer : int, optional
        Buffer around polygons or points for the extent of the timelapse. The default is 100.
    vector_display_path : str, optional
        Optionnal, Path of a vector to display in the timelapse, ground observations, landmarks... The default is None.
    hover_column_list : list, optional
        Optionnal, list of columns to display when hovering mouse over vectors of vector_display_path.
    max_date: str
        Last date used in the timelapse, if None all dates available are used.
    show_confidence_class: bool
        If True, detected dieback is shown with the confidence class at the last date used, as vectorized in the step [05_export_results](https://fordead.gitlab.io/fordead_package/docs/user_guides/english/05_export_results/)
    zip_results: bool
        If True, transfers the timelapse to a zip file "Timelapses.zip" in data_directory/Timelapses directory
    """
    #%%
    
    tile = TileInfo(data_directory)
    tile = tile.import_info()
    # tile.add_parameters({"shape_path" : shape_path})
    # if tile.parameters["Overwrite"] : tile.delete_dirs("timelapse") #Deleting previous detection results if they exist
    tile.add_dirpath("timelapse", tile.data_directory / "Timelapses")
    tile.save_info()
    #%% Importing/creating ROI
    
    if (x is not None) and (y is not None):
        print("Timelapse created from coordinates")
        shape_interest = gp.GeoDataFrame({"id" : [str(x)+"_"+str(y)]},geometry = gp.points_from_xy([x], [y]),crs =tile.raster_meta["crs"] )
        name_column = "id"
    elif shape_path is not None:
        print("Timelapse(s) created from " + shape_path)
        shape_interest=gp.read_file(shape_path)
        shape_interest=shape_interest.to_crs(crs = tile.raster_meta["crs"])
    else:
        raise Exception("No shape_path or coordinates")
    #%% Creating zip file
    
    if zip_results:
        if not (tile.paths["timelapse"] / "Timelapses.zip").exists():
            zipObj = ZipFile(tile.paths["timelapse"] / "Timelapses.zip", 'w', compression = ZIP_DEFLATED, compresslevel = 6)
        else:
            zipObj = ZipFile(tile.paths["timelapse"] / "Timelapses.zip", 'a', compression = ZIP_DEFLATED, compresslevel = 6)
    #%%Creating timelapses
 
    for shapeindex in range(shape_interest.shape[0]):
        shape=shape_interest.iloc[shapeindex:(shapeindex+1)]
        try:
            NameFile = str(shape[name_column].iloc[0])
        except KeyError:
            raise Exception("No column "+name_column+" in " + str(shape_path))
        
        # if not((tile.paths["timelapse"] / (NameFile + ".html")).exists()): 

        if intersects_area(shape.geometry.buffer(buffer), tile.raster_meta["extent"], tile.raster_meta["crs"]):
            print("Creating timelapse | Id : " + NameFile)
            fig = CreateTimelapse(shape.geometry.buffer(buffer),tile, vector_display_path, hover_column_list, max_date, show_confidence_class)
            plot(fig,filename=str(tile.paths["timelapse"] / (NameFile + ".html")),auto_open=False)
            if zip_results: 
                zipObj.write(str(tile.paths["timelapse"] / (NameFile + ".html")),NameFile + ".html") #Adding to zipfile
                (tile.paths["timelapse"] / (NameFile + ".html")).unlink() #Removing html file
        else:
            print("Chosen area for the timelapse does not intersect extent of computed area")
    if zip_results: zipObj.close()

if __name__ == '__main__':
    cli_create_timelapse()

