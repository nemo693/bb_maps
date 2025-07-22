# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 13:01:07 2023

@author: Raphael Dutrieux
"""
import pandas as pd
import time
from fordead.validation_process import add_diff_vi_to_vi, add_status_to_vi, fill_periods, detect_state_changes, detection_anomalies_dataframe, compute_stress_index, prediction_vegetation_index_dataframe
from pathlib import Path
import click

@click.command(name='calval_dieback_detection')
@click.option("--masked_vi_path", type = str, help = "Path of the csv containing the vegetation index for each pixel of each observation, for each valid Sentinel-2 acquisition.")
@click.option("--pixel_info_path", type = str, help = "Path used to write the csv containing pixel info such as the validity of the model and its coefficients.")
@click.option("--periods_path", type = str, help = "Path of the csv containing pixel periods")
@click.option("--name_column", type = str,default = "id", help = "Name of the ID column", show_default=True)
@click.option("--update_masked_vi",  is_flag=True, help = "If True, updates the csv at masked_vi_path with the columns 'period_id', 'state', 'predicted_vi', 'diff_vi' and 'anomaly'", show_default=True)
@click.option("-s", "--threshold_anomaly",  type=float, default=0.16, help="Minimum threshold for anomaly detection", show_default=True)
@click.option("--stress_index_mode", type=str, default = None, help="Chosen stress index, if 'mean', the index is the mean of the difference between the vegetation index and the predicted vegetation index for all unmasked dates after the first anomaly subsequently confirmed. If 'weighted_mean', the index is a weighted mean, where for each date used, the weight corresponds to the number of the date (1, 2, 3, etc...) from the first anomaly. If None, no stress period index is computed.", show_default=True)
@click.option("--vi", type = str,default = "CRSWIR", help = "Chosen vegetation index", show_default=True)
@click.option("--path_dict_vi", type = str,default = None, help = "Path to a text file used to add potential vegetation indices. If not filled in, only the indices provided in the package can be used (CRSWIR, NDVI, NDWI). The file [ex_dict_vi.txt](https://gitlab.com/fordead/fordead_package/-/blob/master/docs/examples/ex_dict_vi.txt) gives an example for how to format this file. One must fill the index's name, formula, and '+' or '-' according to whether the index increases or decreases when anomalies occur.", show_default=True)
def cli_dieback_detection_from_dataframe(masked_vi_path, pixel_info_path, periods_path, name_column, update_masked_vi, threshold_anomaly, vi, path_dict_vi):
    """
    Updates of the csv file containing periods, so for each pixel, the whole time series is covered with the first and last unmasked Sentinel-2 acquisition date of each period, and its associated state.
    The 'state' column can now hold the following values :
        - Training : Period used in training the harmonic model
        - Healthy : Period detected as healthy, with no stress, dieback or bare ground detected
        - Stress : Period beginning with 3 successive anomalies, ending with the last anomaly before three successive non-anomalies of the beginning of a 'Healthy' period.
        - Dieback : Period beginning with 3 successive anomalies, ending with the last available acquisition, or the beggining of a Bare ground period.
        - Invalid : The pixel is invalid, there were not enough valid acquisitions to compute a harmonic model
    A new column 'anomaly_intensity' is also added, it is a weighted mean of the difference between the calculated vegetation indices and their predicted value for the period. The weight is the number of the date within that period (1+2+3+...+ nb_dates). It is only calculated for 'Healthy', 'Stress' and 'Dieback' periods

    if 'update_masked_vi' is True, this function also updates the csv at 'masked_vi_path' with the following columns:
        - period_id : id of the period the acquisition is associated with
        - state : Status of of the associated period, can be 'Training', 'Healthy', 'Stress', 'Dieback' or 'Invalid'.
        - predicted_vi : The prediction of the vegetation index using the harmonic model
        - diff_vi : Difference between the vegetation and its prediction, in the expected direction of anomalies for the vegetation index
        - anomaly : True if 'diff_vi' exceeds 'threshold_anomaly', else False
    
    See additional information [here](https://fordead.gitlab.io/fordead_package/docs/user_guides/english/validation_tools/07_dieback_detection_from_dataframe/)
    """
    dieback_detection_from_dataframe(**locals())
    

def dieback_detection_from_dataframe(masked_vi_path, pixel_info_path, periods_path, name_column = "id", update_masked_vi = False,
                                     threshold_anomaly = 0.16, stress_index_mode = None, vi = "CRSWIR", path_dict_vi = None):
    """
    
    Updates of the csv file containing periods, so for each pixel, the whole time series is covered with the first and last unmasked Sentinel-2 acquisition date of each period, and its associated state.
    The 'state' column can now hold the following values :
        - Training : Period used in training the harmonic model
        - Healthy : Period detected as healthy, with no stress, dieback or bare ground detected
        - Stress : Period beginning with 3 successive anomalies, ending with the last anomaly before three successive non-anomalies of the beginning of a 'Healthy' period.
        - Dieback : Period beginning with 3 successive anomalies, ending with the last available acquisition, or the beggining of a Bare ground period.
        - Invalid : The pixel is invalid, there were not enough valid acquisitions to compute a harmonic model
    A new column 'anomaly_intensity' is also added, it is a weighted mean of the difference between the calculated vegetation indices and their predicted value for the period. The weight is the number of the date within that period (1+2+3+...+ nb_dates). It is only calculated for 'Healthy', 'Stress' and 'Dieback' periods

    if 'update_masked_vi' is True, this function also updates the csv at 'masked_vi_path' with the following columns:
        - period_id : id of the period the acquisition is associated with
        - state : Status of of the associated period, can be 'Training', 'Healthy', 'Stress', 'Dieback' or 'Invalid'.
        - predicted_vi : The prediction of the vegetation index using the harmonic model
        - diff_vi : Difference between the vegetation and its prediction, in the expected direction of anomalies for the vegetation index
        - anomaly : True if 'diff_vi' exceeds 'threshold_anomaly', else False
    
    See additional information [here](https://fordead.gitlab.io/fordead_package/docs/user_guides/english/validation_tools/07_dieback_detection_from_dataframe/)

    Parameters
    ----------
    masked_vi_path : str
        Path of the csv containing the vegetation index for each pixel of each observation, for each valid Sentinel-2 acquisition. Must have the following columns : epsg, area_name, id, id_pixel, Date and vi.
    pixel_info_path : str
        Path of the the csv containing pixel info.
    periods_path : str
        Path of the csv containing pixel periods.
    name_column : str
        Name of the ID column. The default is 'id'.
    update_masked_vi : bool, optional
        If True, updates the csv at masked_vi_path with the columns 'period_id', 'state', 'predicted_vi', 'diff_vi' and 'anomaly'. The default is False.
    threshold_anomaly : float, optional
        Threshold at which the difference between the actual and predicted vegetation index is considered as an anomaly. The default is 0.16.
    stress_index_mode : str, optional
        Chosen stress index, if 'mean', the index is the mean of the difference between the vegetation index and the predicted vegetation index for all unmasked dates after the first anomaly subsequently confirmed.
        If 'weighted_mean', the index is a weighted mean, where for each date used, the weight corresponds to the number of the date (1, 2, 3, etc...) from the first anomaly.
        If None, no stress period index is computed.
    vi : str, optional
        Chosen vegetation index. The default is "CRSWIR".
    path_dict_vi : str, optional
        Path to a text file used to add potential vegetation indices. If not filled in, only the indices provided in the package can be used (CRSWIR, NDVI, NDWI). The file [ex_dict_vi.txt](https://gitlab.com/fordead/fordead_package/-/blob/master/docs/examples/ex_dict_vi.txt) gives an example for how to format this file. One must fill the index's name, formula, and "+" or "-" according to whether the index increases or decreases when anomalies occur. The default is None.


    """
    
    
    masked_vi = pd.read_csv(masked_vi_path)
    pixel_info = pd.read_csv(pixel_info_path)
    periods = pd.read_csv(periods_path)
    
    if "anomaly_intensity" in periods.columns:
        raise Exception(str(periods_path) + " already exists.")
        
    pixel_info["coeff"] = pixel_info[["coeff1","coeff2","coeff3", "coeff4","coeff5"]].values.tolist()
    
    # masked_vi = remove_bare_ground(masked_vi,periods)
    
    
    if "bare_ground" in masked_vi.columns:
        # period_bare_ground = masked_vi[masked_vi["bare_ground"]]
        # period_bare_ground = period_bare_ground.groupby(by = ["area_name", name_column,"id_pixel"]).first().reset_index()
        # period_bare_ground["state"] = "Bare ground"
        # period_bare_ground = period_bare_ground.rename(columns={"Date": "start_date"}).drop(columns = ["vi", "bare_ground", "epsg"])
        detection_dates = masked_vi[~masked_vi["bare_ground"]]
    else:
        detection_dates = masked_vi
    
    # test = masked_vi[(masked_vi[name_column] == 5) & (masked_vi["id_pixel"] == 41)]
    
    detection_dates = detection_dates.merge(pixel_info, on=["epsg", "area_name",name_column,"id_pixel"], how='left')
   
    # training_periods = get_training_period(masked_vi, name_column)
    
    detection_dates = detection_dates[detection_dates["Date"] > detection_dates["last_training_date"]].reset_index()

    detection_dates["predicted_vi"] = prediction_vegetation_index_dataframe(detection_dates, pixel_info, name_column)
    
    detection_dates["anomaly"], detection_dates["diff_vi"] = detection_anomalies_dataframe(detection_dates, threshold_anomaly, vi = vi, path_dict_vi = path_dict_vi)

    dated_changes = detect_state_changes(detection_dates, name_column)
 
    periods = fill_periods(periods,dated_changes, masked_vi, name_column)

    if stress_index_mode is not None: periods = compute_stress_index(detection_dates, periods, name_column, stress_index_mode)
    
    if update_masked_vi:
        masked_vi = add_diff_vi_to_vi(masked_vi, pixel_info, threshold_anomaly, vi, path_dict_vi, name_column)
        masked_vi = add_status_to_vi(masked_vi, periods, name_column, stress_index_mode)
        masked_vi.to_csv(masked_vi_path, mode='w', index=False,header=True)
        
        
    # periods = compute_stress_index(detection_dates, periods, name_column)
    # dieback_data = dieback_detection(detection_dates, dated_changes, name_column)


    periods.to_csv(periods_path, mode='w', index=False,header=True)
    # dieback_data.to_csv(Path(export_directory) / 'dieback_data.csv', mode='w', index=False,header=True)

    
if __name__ == '__main__':
    start_time_debut = time.time()
    #Exemple tuto
    output_dir = Path("D:/fordead/fordead_data/output/calval_tuto")

    dieback_detection_from_dataframe(masked_vi_path = output_dir / "mask_vi_tuto.csv",
                                      pixel_info_path = output_dir / "pixel_info_tuto.csv",
                                        periods_path = output_dir / "periods_tuto.csv",
                                        name_column = "id",
                                        update_masked_vi = True)
    # output_dir = Path("D:/fordead/fordead_data/output/calval_tuto")

    # dieback_detection_from_dataframe(masked_vi_path = output_dir / "masked_vi_scolytes.csv",
    #                                    pixel_info_path = output_dir / "pixel_info_tuto.csv",
    #                                    periods_path = output_dir / "periods_scolytes.csv",
    #                                    update_masked_vi = True,
    #                                    name_column = "Id")
    print("Temps de calcul : %s secondes ---" % (time.time() - start_time_debut))
