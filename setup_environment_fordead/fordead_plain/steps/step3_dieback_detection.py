# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 09:25:23 2020

@author: Raphaël Dutrieux
"""


import click
from fordead.import_data import import_coeff_model, import_dieback_data, import_stress_data, initialize_dieback_data, initialize_stress_data, import_masked_vi, import_first_detection_date_index, TileInfo, import_binary_raster
from fordead.writing_data import write_tif
from fordead.dieback_detection import detection_anomalies, detection_dieback, save_stress
from fordead.model_vegetation_index import prediction_vegetation_index, correct_vi_date

@click.command(name='dieback_detection')
@click.option("-o", "--data_directory",  type=str, help="Path of the output directory")
@click.option("-s", "--threshold_anomaly",  type=float, default=0.16,
                    help="Minimum threshold for anomaly detection", show_default=True)
@click.option("--max_nb_stress_periods",  type=int, default=5,
                    help="Maximum number of stress periods", show_default=True)
@click.option("--stress_index_mode",  type=str, default=None,
                    help="Chosen stress index, if 'mean', the index is the mean of the difference between the vegetation index and the predicted vegetation index for all unmasked dates after the first anomaly subsequently confirmed. If 'weighted_mean', the index is a weighted mean, where for each date used, the weight corresponds to the number of the date (1, 2, 3, etc...) from the first anomaly. If None, the stress periods are not detected, and no informations on stress periods are saved.", show_default=True)
@click.option("--vi",  type=str, default=None,
                    help="Chosen vegetation index, only useful if step1 was skipped", show_default=True)
@click.option("--path_dict_vi",  type=str, default=None,
                    help="Path of text file to add vegetation index formula, only useful if step1 was skipped", show_default=True)
def cli_dieback_detection(
    data_directory,
    threshold_anomaly=0.16,
    max_nb_stress_periods = 5,
    stress_index_mode = None,
    vi = None,
    path_dict_vi = None
    
    ):
    """
    Detects anomalies by comparing the vegetation index and its prediction from the model. 
    Detects pixels suffering from dieback when there are 3 successive anomalies. If pixels detected as suffering from dieback have 3 successive dates without anomalies, they are considered healthy again.
    If stress_index_mode parameter is given, Those periods between detection and return to normal are saved, with the date of first anomaly, date of return to normal, number of dates, and an associated stress index.
    Anomalies and dieback data are written in the data_directory
    See details here : https://fordead.gitlab.io/fordead_package/docs/user_guides/english/03_dieback_detection/
    \f
    """
    dieback_detection(data_directory, threshold_anomaly, max_nb_stress_periods, stress_index_mode, vi, path_dict_vi)


def dieback_detection(
    data_directory,
    threshold_anomaly=0.16,
    max_nb_stress_periods = 5,
    stress_index_mode = None,
    vi = None,
    path_dict_vi = None
    ):
    """
    Detects anomalies by comparing the vegetation index and its prediction from the model. 
    Detects pixels suffering from dieback when there are 3 successive anomalies. If pixels detected as suffering from dieback have 3 successive dates without anomalies, they are considered healthy again.
    If stress_index_mode parameter is given, information on those periods between detection and return to normal are saved, with the date of first anomaly, date of return to normal, number of dates, an associated stress index, and the total number of those periods.
    Anomalies and dieback data are written in the data_directory
    See details here : https://fordead.gitlab.io/fordead_package/docs/user_guides/english/03_dieback_detection/
    
    \f
    Parameters
    ----------
    data_directory : str
        Path of the output directory
    threshold_anomaly : float
        Minimum threshold for anomaly detection
    max_nb_stress_periods : int
        Maximum number of stress periods, if this number is reached, the pixel is masked in the too_many_stress_periods, thus removed from future exports. Only used if stress_index_mode is not None.
    stress_index_mode : str
        Chosen stress index, if 'mean', the index is the mean of the difference between the vegetation index and the predicted vegetation index for all unmasked dates after the first anomaly subsequently confirmed.
        If 'weighted_mean', the index is a weighted mean, where for each date used, the weight corresponds to the number of the date (1, 2, 3, etc...) from the first anomaly.
        If None, the stress periods are not detected, and no informations on stress periods are saved.
    vi : str
        Chosen vegetation index, only useful if step1 was skipped
    path_dict_vi : str
        Path of text file to add vegetation index formula, only useful if step1 was skipped

    Returns
    -------

    """
    tile = TileInfo(data_directory)
    tile = tile.import_info()
    tile.add_parameters({"threshold_anomaly" : threshold_anomaly, "max_nb_stress_periods" : max_nb_stress_periods, "stress_index_mode" : stress_index_mode})
    if tile.parameters["Overwrite"] : 
        tile.delete_dirs("AnomaliesDir","state_dieback","periodic_results_dieback","result_files","timelapse","series","nb_periods_stress") #Deleting previous detection results if they exist
        tile.delete_files("too_many_stress_periods_mask")
        tile.delete_attributes("last_computed_anomaly","last_date_export")

    if vi==None : vi = tile.parameters["vi"]
    if path_dict_vi==None : path_dict_vi = tile.parameters["path_dict_vi"] if "path_dict_vi" in tile.parameters else None
    
    tile.add_dirpath("AnomaliesDir", tile.data_directory / "DataAnomalies") #Choose anomalies directory
    
    #added
    tile.add_dirpath("AnomaliesValuesDir", tile.data_directory / "ValuesAnomalies") #Choose anomalies directory
    
    tile.getdict_datepaths("Anomalies",tile.paths["AnomaliesDir"]) # Get paths and dates to previously calculated anomalies
    tile.search_new_dates() #Get list of all used dates
    
    tile.add_path("too_many_stress_periods_mask", tile.data_directory / "TimelessMasks" / "too_many_stress_periods_mask.tif")
    
    tile.add_path("state_dieback", tile.data_directory / "DataDieback" / "state_dieback.tif")
    tile.add_path("first_date_dieback", tile.data_directory / "DataDieback" / "first_date_dieback.tif")
    tile.add_path("first_date_unconfirmed_dieback", tile.data_directory / "DataDieback" / "first_date_unconfirmed_dieback.tif")
    tile.add_path("count_dieback", tile.data_directory / "DataDieback" / "count_dieback.tif")
    
    tile.add_path("dates_stress", tile.data_directory / "DataStress" / "dates_stress.tif")
    tile.add_path("nb_periods_stress", tile.data_directory / "DataStress" / "nb_periods_stress.tif")
    tile.add_path("cum_diff_stress", tile.data_directory / "DataStress" / "cum_diff_stress.tif")
    tile.add_path("nb_dates_stress", tile.data_directory / "DataStress" / "nb_dates_stress.tif")
    tile.add_path("stress_index", tile.data_directory / "DataStress" / "stress_index.tif")

    #Verify if there are new SENTINEL dates
    #new_dates = tile.dates[tile.dates > tile.last_computed_anomaly] if hasattr(tile, "last_computed_anomaly") else tile.dates[tile.dates >= tile.parameters["min_last_date_training"]]
    
       ###### modded, comment added
    new_dates = tile.dates[tile.dates > tile.last_computed_anomaly] if hasattr(tile, "last_computed_anomaly") else tile.dates[tile.dates >= tile.parameters["min_last_date_training"]] # get new dates or all dates after the start of detection
    
    #############
    
    if  len(new_dates) == 0:
        print("Dieback detection : no new dates")
    else:
        print("Dieback detection : " + str(len(new_dates))+ " new dates")
        
        #IMPORTING DATA
        first_detection_date_index = import_first_detection_date_index(tile.paths["first_detection_date_index"])
        coeff_model = import_coeff_model(tile.paths["coeff_model"])
        
        if tile.paths["state_dieback"].exists():
            dieback_data = import_dieback_data(tile.paths)
        else:
            dieback_data = initialize_dieback_data(first_detection_date_index.shape,first_detection_date_index.coords)
        if tile.paths["nb_dates_stress"].exists(): 
            stress_data = import_stress_data(tile.paths)
        else:
            stress_data = initialize_stress_data(first_detection_date_index.shape,first_detection_date_index.coords, max_nb_stress_periods)
   
        if tile.parameters["correct_vi"]:
            forest_mask = import_binary_raster(tile.paths["forest_mask"])
            
        #dieback DETECTION
        for date_index, date in enumerate(tile.dates):
            if date in new_dates:
                vegetation_index, mask = import_masked_vi(tile.paths,date)
                if tile.parameters["correct_vi"]:
                    vegetation_index, tile.correction_vi = correct_vi_date(vegetation_index, mask,forest_mask, tile.large_scale_model, date, tile.correction_vi)

                mask = mask | (date_index < first_detection_date_index) #Masking pixels where date was used for training
                
                predicted_vi=prediction_vegetation_index(coeff_model,[date])
                
                anomalies, diff_vi = detection_anomalies(vegetation_index, mask, predicted_vi, threshold_anomaly, 
                                                vi = vi, path_dict_vi = path_dict_vi)
                                
                dieback_data, changing_pixels = detection_dieback(dieback_data, anomalies, mask, date_index)
                
                if stress_index_mode is not None: stress_data = save_stress(stress_data, dieback_data, changing_pixels, diff_vi, mask, stress_index_mode) 

                write_tif(anomalies, first_detection_date_index.attrs, tile.paths["AnomaliesDir"] / str("Anomalies_" + date + ".tif"),nodata=0)
                print('\r', date, " | ", len(tile.dates)-date_index-1, " remaining         ", sep='', end='', flush=True) if date_index != (len(tile.dates) -1) else print('\r', "                                              ", sep='', end='\r', flush=True) 
                del vegetation_index, mask, predicted_vi, anomalies, changing_pixels, diff_vi

        tile.last_computed_anomaly = new_dates[-1]
  
        write_tif(dieback_data["state"], first_detection_date_index.attrs,tile.paths["state_dieback"],nodata=0)
        write_tif(dieback_data["first_date"], first_detection_date_index.attrs,tile.paths["first_date_dieback"],nodata=0)
        write_tif(dieback_data["first_date_unconfirmed"], first_detection_date_index.attrs,tile.paths["first_date_unconfirmed_dieback"],nodata=0)
        write_tif(dieback_data["count"], first_detection_date_index.attrs,tile.paths["count_dieback"],nodata=0)
        del dieback_data
        
        if stress_index_mode is not None:
            # valid_model = import_binary_raster(tile.paths["sufficient_coverage_mask"])
            # valid_model = valid_model.where(stress_data["nb_periods"]<=max_nb_stress_periods,False)
            too_many_stress_periods_mask = stress_data["nb_periods"]<=max_nb_stress_periods
            write_tif(too_many_stress_periods_mask, first_detection_date_index.attrs,tile.paths["too_many_stress_periods_mask"],nodata=0) 

            if stress_index_mode == "mean":
                stress_index = stress_data["cum_diff"]/stress_data["nb_dates"]
            elif stress_index_mode == "weighted_mean":
                stress_index = stress_data["cum_diff"]/(stress_data["nb_dates"]*(stress_data["nb_dates"]+1)/2)
            else:
                raise Exception("Unrecognized stress_index_mode")
                   
            write_tif(stress_index, first_detection_date_index.attrs,tile.paths["stress_index"],nodata=0)
            del stress_index
            #Writing dieback data to rasters
            write_tif(stress_data["date"], first_detection_date_index.attrs,tile.paths["dates_stress"],nodata=0)
            write_tif(stress_data["nb_periods"], first_detection_date_index.attrs,tile.paths["nb_periods_stress"],nodata=0)
            write_tif(stress_data["cum_diff"], first_detection_date_index.attrs,tile.paths["cum_diff_stress"],nodata=0)
            write_tif(stress_data["nb_dates"], first_detection_date_index.attrs,tile.paths["nb_dates_stress"],nodata=0)
            del stress_data

        # print("Détection du déperissement")
    tile.save_info()



if __name__ == '__main__':
    # print(dictArgs)
    # start_time = time.time()
    cli_dieback_detection()
    # print("Temps d execution : %s secondes ---" % (time.time() - start_time))
