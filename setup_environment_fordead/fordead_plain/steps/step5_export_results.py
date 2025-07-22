# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 11:32:57 2020

@author: Raphaël Dutrieux
"""

import click
from fordead.import_data import import_dieback_data, TileInfo, import_binary_raster, import_soil_data, import_stress_data, import_stress_index
from fordead.writing_data import vectorizing_confidence_class, get_bins, convert_dateindex_to_datenumber, get_periodic_results_as_shapefile, get_state_at_date, union_confidence_class, write_tif
import numpy as np
import geopandas as gp
import pandas as pd
import warnings
warnings.filterwarnings("ignore", message="`unary_union` returned None due to all-None GeoSeries. In future, `unary_union` will return 'GEOMETRYCOLLECTION EMPTY' instead.")

@click.command(name='export_results')
@click.option("-o", "--data_directory",  type=str, help="Path of the output directory", show_default=True)
@click.option("--start_date",  type=str, default='2015-06-23',
                    help="Start date for exporting results", show_default=True)
@click.option("--end_date",  type=str, default="2030-01-02",
                    help="End date for exporting results", show_default=True)
@click.option("--frequency",  type=str, default='M',
                    help="Frequency used to aggregate results, if value is 'sentinel', then periods correspond to the period between sentinel dates used in the detection, or it can be the frequency as used in pandas.date_range. e.g. 'M' (monthly), '3M' (three months), '15D' (fifteen days)", show_default=True)
@click.option("--multiple_files",  is_flag=True,
                    help="If True, one shapefile is exported for each period containing the areas in dieback at the end of the period. Else, a single shapefile is exported containing diebackd areas associated with the period of dieback", show_default=True)
@click.option("-t", "--conf_threshold_list", type = float, multiple=True, default = None, help = "List of thresholds used as bins to discretize the confidence index into several classes", show_default=True)
@click.option("-c", "--conf_classes_list", type = str, multiple=True, default = None, help = "List of classes names, if conf_threshold_list has n values, conf_classes_list must have n+1 values", show_default=True)
def cli_export_results(
    data_directory,
    start_date = '2015-06-23',
    end_date = "2030-01-02",
    frequency = 'M',
    multiple_files = False,
    conf_threshold_list = None,
    conf_classes_list = None,
    ):
    """
    Export results to a vectorized shapefile format.
    \f

    """
    export_results(data_directory, start_date, end_date, frequency, multiple_files, conf_threshold_list, conf_classes_list)



def export_results(
    data_directory,
    start_date = '2015-06-23',
    end_date = "2030-01-02",
    frequency = 'M',
    multiple_files = False,
    conf_threshold_list = None,
    conf_classes_list = None,
    ):
    """
    Writes results in the chosen period, form and using chosen frequency.
    See details here : https://fordead.gitlab.io/fordead_package/docs/user_guides/english/05_export_results/
    \f

    Parameters
    ----------
    data_directory : str
        Path of the output directory
    start_date : str
        Start date for exporting results (format : 'YYYY-MM-DD')
    end_date : str
        End date for exporting results (format : 'YYYY-MM-DD')
    frequency : str
        Frequency used to aggregate results, if value is 'sentinel', then periods correspond to the period between sentinel dates used in the detection, or it can be the frequency as used in pandas.date_range. e.g. 'M' (monthly), '3M' (three months), '15D' (fifteen days)
    multiple_files : bool
        If True, one shapefile is exported for each period containing the areas in dieback at the end of the period. Else, a single shapefile is exported containing diebackd areas associated with the period of dieback
    conf_threshold_list : list
        List of thresholds used as bins to discretize the confidence index into several classes
    conf_classes_list : list
        List of classes names, if conf_threshold_list has n values, conf_classes_list must have n+1 values
    
    Returns
    -------

    """    
    tile = TileInfo(data_directory)
    tile = tile.import_info()
    dieback_data = import_dieback_data(tile.paths, chunks= None)
    tile.add_parameters({"start_date" : start_date,"end_date" : end_date, "frequency" : frequency, "multiple_files" : multiple_files, "conf_threshold_list": conf_threshold_list, "conf_classes_list" : conf_classes_list})
    if tile.parameters["Overwrite"] : 
        tile.delete_dirs("periodic_results_dieback","result_files") #Deleting previous detection results if they exist
        tile.delete_attributes("last_date_export")
        
    tile.add_path("confidence_index", tile.data_directory / "Results" / "confidence_index.tif")
    # tile.add_path("stress_periods", tile.data_directory / "Results" / "stress_periods.shp")
    tile.add_path("stress_periods", tile.data_directory / "Results" / "stress_periods.shp")

    exporting = (tile.dates[-1] != tile.last_date_export) if hasattr(tile, "last_date_export") else True
    # exporting = True
    if exporting:
        print("Exporting results")
        bins_as_date, bins_as_datenumber = get_bins(start_date,end_date,frequency,tile.dates)
        first_date_number = convert_dateindex_to_datenumber(dieback_data.first_date,dieback_data.state, tile.dates)
    
        if tile.parameters["soil_detection"]:
            soil_data = import_soil_data(tile.paths, chunks= 1280)
            first_date_number_soil = convert_dateindex_to_datenumber(soil_data.first_date, soil_data.state, tile.dates)
            

        forest_mask = import_binary_raster(tile.paths["forest_mask"], chunks= 1280)
        sufficient_coverage_mask = import_binary_raster(tile.paths["sufficient_coverage_mask"], chunks= 1280)
        if tile.parameters["stress_index_mode"] is not None:
            too_many_stress_periods_mask = import_binary_raster(tile.paths["too_many_stress_periods_mask"], chunks= 1280)
            relevant_area = forest_mask & sufficient_coverage_mask & too_many_stress_periods_mask
        else:
            relevant_area = forest_mask & sufficient_coverage_mask
      
        #EXPORTING STRESS RESULTS
        if tile.parameters["stress_index_mode"] is not None:
            if conf_threshold_list is None or conf_classes_list is None or len(conf_threshold_list) == 0 or len(conf_classes_list) == 0:
                print("Parameters conf_threshold_list and conf_classes_list are not provided, stress results can't be exported")
            else:
                stress_index = import_stress_index(tile.paths["stress_index"])
                stress_data = import_stress_data(tile.paths)
                stress_list = []
                for period in range(tile.parameters["max_nb_stress_periods"]):
                    stress_period = stress_index.isel(period = period)
                    conf_threshold_list = np.array(conf_threshold_list).astype(float)
                    stress_class = vectorizing_confidence_class(stress_period, stress_data.nb_dates.isel(period = period), (relevant_area & (period < stress_data["nb_periods"])).compute(), conf_threshold_list, np.array(conf_classes_list), tile.raster_meta["attrs"])
                    if stress_class.size != 0:
                        stress_start_date = stress_data["date"].isel(change = period*2)
                        stress_start_date_number = convert_dateindex_to_datenumber(stress_start_date, period < stress_data["nb_periods"], tile.dates)
                        vector_stress_start = get_periodic_results_as_shapefile(stress_start_date_number, bins_as_date, bins_as_datenumber, relevant_area, forest_mask.attrs)
                        # stress_class.to_file(tile.paths["stress_periods"] / ("stress_class" + str(period+1) + ".shp"))
                        # vector_stress_start.to_file(tile.paths["stress_periods"] / ("stress_period" + str(period+1) + ".shp"))
                        stress_list += [gp.overlay(vector_stress_start, stress_class, how='intersection',keep_geom_type = True)]
                # for period_date in pd.unique(vector_stress_start["period"]):
                #     for class_stress in pd.unique(stress_class["class"]):
                #         stress_list += [gp.overlay(vector_stress_start[vector_stress_start["period"]==period_date], stress_class[stress_class["class"] == class_stress], how='intersection',keep_geom_type = True)]
                
                if len(stress_list) != 0:
                    stress_total = gp.GeoDataFrame( pd.concat(stress_list, ignore_index=True), crs=stress_list[0].crs)
                    stress_total.to_file(tile.paths["stress_periods"])
                else:
                    print("No stress periods detected")

        if multiple_files:
            tile.add_dirpath("result_files", tile.data_directory / "Results")
            for date_bin_index, date_bin in enumerate(bins_as_date):
                state_code = first_date_number <= bins_as_datenumber[date_bin_index]
                if tile.parameters["soil_detection"]:
                    state_code = state_code + 2*(first_date_number_soil <= bins_as_datenumber[date_bin_index])
                period_end_results = get_state_at_date(state_code,relevant_area,dieback_data.state.attrs)
                if not(period_end_results.empty):
                    period_end_results.to_file(tile.paths["result_files"] / (date_bin.strftime('%Y-%m-%d')+".shp"))
                    
        else:
            #EXPORTING DIEBACK 
            tile.add_path("periodic_results_dieback", tile.data_directory / "Results" / "periodic_results_dieback.shp")
            periodic_results = get_periodic_results_as_shapefile(first_date_number, bins_as_date, bins_as_datenumber, relevant_area, dieback_data.state.attrs)
          
            if conf_threshold_list is not None and conf_classes_list is not None and len(conf_threshold_list) != 0 and len(conf_classes_list) != 0:
                if tile.parameters["stress_index_mode"] is None:
                    print("Confidence index was not saved, parameters conf_threshold_list and conf_classes_list are ignored. Change stress_index_mode parameter in step 3 to compute confidence index")
                else:
                    stress_data = import_stress_data(tile.paths)
                    stress_index = import_stress_index(tile.paths["stress_index"])
                    confidence_area = relevant_area & dieback_data["state"] & ~soil_data["state"] if tile.parameters["soil_detection"] else relevant_area & dieback_data["state"]
               
                    confidence_index = stress_index.sel(period = (stress_data["nb_periods"]+1).where(stress_data["nb_periods"]<=tile.parameters["max_nb_stress_periods"],tile.parameters["max_nb_stress_periods"])) #The selection probably makes no sense for pixels with nb_periods higher that max_nb_stress_periods, but it doesn't matter since they are excluded from result exports, but it removes bugs of inexistant period values.
                    nb_dates = stress_data["nb_dates"].sel(period = (stress_data["nb_periods"]+1).where(stress_data["nb_periods"]<=tile.parameters["max_nb_stress_periods"],tile.parameters["max_nb_stress_periods"]))  #The selection probably makes no sense for pixels with nb_periods higher that max_nb_stress_periods, but it doesn't matter since they are excluded from result exports, but it removes bugs of inexistant period values.
                   
                    write_tif(confidence_index.where(confidence_area,0), forest_mask.attrs,nodata = 0, path = tile.paths["confidence_index"])
                   
                    confidence_class = vectorizing_confidence_class(confidence_index, nb_dates, confidence_area.compute(), conf_threshold_list, np.array(conf_classes_list), tile.raster_meta["attrs"])

                    periodic_results = union_confidence_class(periodic_results, confidence_class)
                    
            if not(periodic_results.empty):
                periodic_results.to_file(tile.paths["periodic_results_dieback"],index = None)
            del periodic_results
            
            #EXPORTING SOIL DETECTION
            if tile.parameters["soil_detection"]:
                tile.add_path("periodic_results_soil", tile.data_directory / "Results" / "periodic_results_soil.shp")
                periodic_results = get_periodic_results_as_shapefile(first_date_number_soil, bins_as_date, bins_as_datenumber, relevant_area, soil_data.state.attrs)
                if not(periodic_results.empty):
                    periodic_results.to_file(tile.paths["periodic_results_soil"])
                del periodic_results
        tile.last_date_export = tile.dates[-1]
        tile.save_info()
        dieback_data.close()
    else:
        print("Results already exported")

if __name__ == '__main__':
    cli_export_results()
    