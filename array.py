"""
array.py

@author: tbezo
"""
import pandas as pd

  
# Class that holds 729 Array measurement data in a Dataframe
class ARRAY:
    
    def __init__(self, mcc_list: list):  
        self.mcclist = mcc_list
        self.array = self.merge_profiles()
        return
  
    def merge_profiles(self) -> pd.DataFrame:
        """
        
        Returns
        -------
        array : TYPE
            DESCRIPTION.

        """
        y = 2
        array = self.mcclist[0].dataframe.drop('reference', 1)
        for i in self.mcclist:
            array.insert(y, "meas_values" + str(y), i.dataframe.meas_values)
            y = y + 1
            
        return array
    
    def upsample_rows(self) -> pd.DataFrame:
        y = 1
        return
