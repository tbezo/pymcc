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
        Merges all xyProfiles into one 2D Dataframe with aditional
        interpolation between the columns. The final resolution is 0.1 mm 
        (which might too much)
        
        Returns
        -------
        array : pd.DataFrame()
            returns an interpolated and upsampled dataframe.

        """
        # starting point for inserting columns
        y = 1
        
        # remove reference values, they are not needed
        array = self.mcclist[0].dataframe.drop('reference', 1)
        
        # remove first profile element from list since it is already contained 
        # in DataFrame array
        del self.mcclist[0]
             
        for i in self.mcclist:
            # add NaN columns for interpolation
            for j in range(1, 100):
                array.insert(y+1, "int" + str(y), np.nan)
                y = y + 1                
            array.insert(y+1, "meas_values" + str(y), i.dataframe.meas_values)
            #print(i.dataframe.meas_values[1300])
            y = y + 1
        
        # interpolate between columns
        array.interpolate(method='linear', axis=1, limit=None, inplace=True)
        # use position as index, the final transposed df will have the
        # position at both index and column names
        array.set_index('position', inplace=True)
        # use index also for column names before transposing
        array.set_axis(array.index, axis='columns', inplace=True)

        # return transposed array                        
        return array.transpose()
 
    
    def downsample(self, ds_val: int = 100) -> pd.DataFrame:
        """
        function that return a subsample of the dataframe. The default is to
        only return the original data values.

        Parameters
        ----------
        ds_val : int, optional
            Downsample value. Only every ds_val row and column gets returned.
            The default is 100.

        Returns
        -------
        pd.DataFrame
            by ds_val downsampled dataframe.

        """
        return self.dataframe.iloc[::ds_val, ::ds_val]


# Class that holds Starcheck measurement data in a list of XyProfiles
class STARCHECK:
    
    def __init__(self, mcc_list: list):  
        self.mcclist = mcc_list
        
  
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