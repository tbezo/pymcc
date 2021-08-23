"""
array.py

@author: tbezo
"""
import pandas as pd
import numpy as np

  
# Class that holds 729 Array measurement data in a Dataframe
class SEVEN29:
    
    def __init__(self, mcc_list: list):  
        self.mcclist = mcc_list
        self.dataframe = self.merge_profiles()
  
    def merge_profiles(self) -> pd.DataFrame:
        """
        merges all single xyProfiles into one 2D Dataframe with aditional
        interpolation between the profiles. The final resolution is 0.1 mm.
        
        Returns
        -------
        array : pd.DataFrame()
            returns an interpolated and upsampled dataframe.

        """
        y = 2
        array = self.mcclist[0].dataframe.drop('reference', 1)
        del self.mcclist[0]
        for i in self.mcclist:
            for j in range(1, 100):
                array.insert(y, "int" + str(y-1), np.nan)
                y = y + 1                
            array.insert(y, "meas_values" + str(y-1), i.dataframe.meas_values)
            y = y + 1
        
        array.interpolate(method='linear', axis=1, limit=None, inplace=True)
        
        return array


# Class that holds Starcheck measurement data in a list of XyProfiles
class STARCHECK:
    
    def __init__(self, mcc_list: list):  
        self.mcclist = mcc_list
        self.dataframe = self.merge_profiles()  
  
    # def merge_profiles(self) -> pd.DataFrame:
    #     """
    #     merges all single xyProfiles into one 2D Dataframe with aditional
    #     interpolation between the profiles. The final resolution is 0.1 mm.
        
    #     Returns
    #     -------
    #     array : pd.DataFrame()
    #         returns an interpolated and upsampled dataframe.

    #     """
    #     y = 2
    #     array = self.mcclist[0].dataframe.drop('reference', 1)
    #     del self.mcclist[0]
    #     for i in self.mcclist:
    #         for j in range(1, 100):
    #             array.insert(y, "int" + str(y-1), np.nan)
    #             y = y + 1                
    #         array.insert(y, "meas_values" + str(y-1), i.dataframe.meas_values)
    #         y = y + 1
        
    #     #array.interpolate(method='linear', axis=1, limit=None, inplace=True)
        
    #     return array