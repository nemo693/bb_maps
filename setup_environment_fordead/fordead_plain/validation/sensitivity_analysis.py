# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 18:55:28 2023

@author: Raphael Dutrieux
"""

import time
import click
import shutil
from pathlib import Path
from itertools import product
import pandas as pd

from fordead.validation_process import get_args, filter_args, combine_validation_results, get_test_id, get_default_args
from fordead.validation.mask_vi_from_dataframe import mask_vi_from_dataframe
from fordead.validation.train_model_from_dataframe import train_model_from_dataframe
from fordead.validation.dieback_detection_from_dataframe import dieback_detection_from_dataframe

def get_args_dataframe(comb_dict):
    args_dict = get_default_args(mask_vi_from_dataframe)
    d2 = get_default_args(train_model_from_dataframe)
    d3 = get_default_args(dieback_detection_from_dataframe)
    args_dict.update(d2)
    args_dict.update(d3)
    args_dict.update(comb_dict)
    for key in args_dict:
        args_dict[key] = str(args_dict[key])
    args_dataframe = pd.DataFrame(data=args_dict, index=[0])
    
    return args_dataframe

@click.command(name='sensitivity_analysis')
@click.option("--testing_directory", type = str,default = None, help = "Directory where the results will be exported", show_default=True)
@click.option("--reflectance_path", type = str,default = None, help = "Path of the csv file with extracted reflectance.", show_default=True)
@click.option("--cloudiness_path", type = str,default = None, help = "Path of a csv with the columns 'area_name','Date' and 'cloudiness', can be calculated by the [extract_cloudiness function](https://fordead.gitlab.io/fordead_package/docs/Tutorials/Validation/03_extract_cloudiness/). Not used if not given.", show_default=True)
@click.option("--args_to_test", type = str,default = None, help = "Path to a text file where each line begins parameter name, then each value is separated with a space. All combinations will be tested. See an example [here](https://gitlab.com/fordead/fordead_package/-/blob/master/docs/examples/ex_dict_args.txt)", show_default=True)
@click.option("--update_masked_vi",  is_flag=True, help = "If True, updates the csv at masked_vi_path with the columns 'period_id', 'state', 'predicted_vi', 'diff_vi' and 'anomaly'", show_default=True)
@click.option("--name_column", type = str,default = "id", help = "Name of the ID column", show_default=True)
def cli_sensitivity_analysis(testing_directory, reflectance_path, cloudiness_path, args_to_test, update_masked_vi, name_column):
    """
    Allows the testing of many parameter combinations, running three detection steps mask_vi_from_dataframe, train_model_from_dataframe and dieback_detection_from_dataframe using default parameters as well as user defined parameter combinations.
    A 'test_info.csv' is written in 'testing_directory', where each test_id is associated with the value of all parameters used in the iteration.
    Each test results can be found in a newly created directory whose name is the test_id.
    
    See additional information [here](https://fordead.gitlab.io/fordead_package/docs/user_guides/english/validation_tools/08_sensitivity_analysis/)

    """
    
    sensitivity_analysis(**locals())



def sensitivity_analysis(testing_directory, reflectance_path, args_to_test, cloudiness_path = None, update_masked_vi = False, name_column = "id"):
    """
    Allows the testing of many parameter combinations, running three detection steps mask_vi_from_dataframe, train_model_from_dataframe and dieback_detection_from_dataframe using default parameters as well as user defined parameter combinations.
    A 'test_info.csv' is written in 'testing_directory', where each test_id is associated with the value of all parameters used in the iteration.
    Each test results can be found in a newly created directory whose name is the test_id.

    See additional information [here](https://fordead.gitlab.io/fordead_package/docs/user_guides/english/validation_tools/08_sensitivity_analysis/)

    Parameters
    ----------
    testing_directory : str
        Directory where the results will be exported.
    reflectance_path : str
        Path of the csv file with extracted reflectance.
    args_to_test : dict or str
        Either a dict where each key is any argument to functions mask_vi_from_dataframe, train_model_from_dataframe or dieback_detection_from_dataframe, and the values are list of values to test. All combinations will be tested.
        If a str is given, it is interpreted as the path to a text file where each line begins parameter name, then each value is separated with a space. See an example [here](https://gitlab.com/fordead/fordead_package/-/blob/master/docs/examples/ex_dict_args.txt)
    cloudiness_path : str, optional
        Path of a csv with the columns 'area_name','Date' and 'cloudiness' used to filter dates, can be calculated by the [extract_cloudiness function](https://fordead.gitlab.io/fordead_package/docs/Tutorials/Validation/03_extract_cloudiness/). Not used if None.
    update_masked_vi : bool, optional
        If True, updates the csv containing the vegetation for each acquisition with the columns 'period_id', 'state', 'predicted_vi', 'diff_vi' and 'anomaly'. Doesn't apply if 'overwrite = True'. The default is False.
    name_column : str, optional
        Name of the ID column. The default is "id".

    """
    
    
    testing_directory = Path(testing_directory)
    test_info_path = testing_directory / "test_info.csv"
        
    args_dict = get_args(args_to_test)
    
    for combs in product(*args_dict.values()):
        print(combs)
        args_dataframe = get_args_dataframe(dict(zip(args_dict.keys(), combs)))
        test_id = get_test_id(test_info_path, args_dataframe)
    
        if test_id is not None:
        
            # if overwrite:
            #     dir_path = testing_directory / "temp_dir"
            # else:
            #     # dir_path = testing_directory / '__'.join([list(args_dict.keys())[i] + "_" + str(combs[i]) for i in range(len(combs))])
            #     dir_path = testing_directory / str(test_id)
                
            dir_path = testing_directory / str(test_id)
            dir_path.mkdir(parents=True, exist_ok=True) 
    
            masked_vi_path = dir_path / "masked_vi.csv"
            pixel_info_path = dir_path / "pixel_info.csv"
            periods_path = dir_path / "periods.csv"
            
            mask_vi_from_dataframe(reflectance_path = reflectance_path,
                                    masked_vi_path = masked_vi_path,
                                    periods_path = periods_path,
                                    name_column = name_column,
                                    cloudiness_path = cloudiness_path,
                                    **filter_args(mask_vi_from_dataframe, args_dict, combs))
         
            
            train_model_from_dataframe(masked_vi_path = masked_vi_path,
                                        pixel_info_path = pixel_info_path,
                                        periods_path = periods_path,
                                        name_column = name_column,
                                        **filter_args(train_model_from_dataframe, args_dict, combs)
                                        )
            
            dieback_detection_from_dataframe(masked_vi_path = masked_vi_path, 
                                             pixel_info_path = pixel_info_path, 
                                             periods_path = periods_path, 
                                             name_column = name_column,
                                             update_masked_vi = update_masked_vi,
                                             **filter_args(dieback_detection_from_dataframe, args_dict, combs))
            

            # combine_validation_results(csv_path_list = [pixel_info_path, 
            #                                             periods_path,
            #                                             masked_vi_path],
            #                             merged_csv_path_list = [testing_directory / "merged_pixel_info.csv", 
            #                                             testing_directory / "merged_periods.csv",
            #                                             testing_directory / "merged_masked_vi.csv"],
            #                             test_info_path = test_info_path,
            #                             args_dataframe = args_dataframe, test_id = test_id)
            
            args_dataframe["test_id"] = test_id
            args_dataframe.to_csv(test_info_path, mode='a', index=False,header=not(test_info_path.exists()))
    
    # if (testing_directory / "temp_dir").exists():
    #     shutil.rmtree(testing_directory / "temp_dir")

    
if __name__ == '__main__':
    start_time_debut = time.time()
    # args_to_test = {"stress_index_mode" : ["mean","weighted_mean"], "threshold_anomaly" : [0.15,0.16], "list_bands" :  ["B2","B3","B4", "B8A", "B11","B12"]}
    args_to_test = {"threshold_anomaly" : [0.16], 
                    "list_bands" :  [["B2","B3","B4", "B8A", "B11","B12"]], 
                    "vi" : ["CRSWIR"]}

    # args_to_test = "D:/fordead/Data/Validation/results_from_raster/Feuillu/args_dict.txt"
    
    sensitivity_analysis(testing_directory = "D:/fordead/fordead_data/output",
                       reflectance_path = "D:/fordead/fordead_data/output/reflectance_tuto.csv",
                       cloudiness_path = "D:/fordead/fordead_data/output/cloudiness_tuto.csv",
                       # tiles = ["T31UGP"],
                       args_to_test = args_to_test,
                       overwrite = True)
    
    print("Temps de calcul : %s secondes ---" % (time.time() - start_time_debut))