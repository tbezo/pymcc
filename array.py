"""
array.py

@author: tbezo
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

  
# Class that holds Octavius 729 Array measurement data in a Dataframe
class OCT729:
    
    def __init__(self, mcc_list: list):  
        self.mcclist = mcc_list
        self.dataframe = self.merge_profiles()
        self.dpmm = 0.1
  
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


    def plot(self):
        """
        Plots the DataFrame content with plt.imshow() and a few x/y-ticks

        Returns
        -------
        None.

        """
        plt.imshow(self.dataframe)
        plt.yticks(np.arange(0.5, len(self.dataframe.index), 400),
                   self.dataframe.index[::400])
        plt.xticks(np.arange(0.5, len(self.dataframe.columns), 400),
                   self.dataframe.columns[::400])

# Class that holds Octavius 1500 Array measurement data in a Dataframe
class OCT1500(OCT729):
    
    def __init__(self, mcc_list: list):
        super().__init__(mcc_list)

  
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
        y = 0
        
        # remove reference values, they are not needed
        array = self.mcclist[0].dataframe.drop('reference', 1)
        # use position as index, the final transposed df will have the
        # position at both index and column names
        array.set_index('position', inplace=True)
        
        # remove first profile element from list since it is already contained 
        # in DataFrame array
        del self.mcclist[0]
             
        for i in self.mcclist:
            i.dataframe.set_index('position', inplace=True)
            # add NaN columns for interpolation
            for j in range(1, 50):
                array.insert(y+1, "int" + str(y), np.nan)
                y = y + 1                
            array.insert(y+1, "meas_values" + str(y), i.dataframe.meas_values)
            #print(i.dataframe.meas_values[1300])
            y = y + 1
        
        # interpolate between columns
        array.interpolate(method='linear', axis=1, limit=None, inplace=True)
        # use index also for column names before transposing
        array.set_axis(array.index, axis='columns', inplace=True)

        # return transposed array                        
        return array.transpose()


# Class that holds Starcheck measurement data in a XyProfiles
class STARCHECK:
    
    def __init__(self, mcc_list: list):  
        self.mcclist = mcc_list
        self.center_crossplane = self.mcclist[3]
        self.center_inplane = self.mcclist[10]
        # target left gun right
        self.diagonal_tlgr = self.mcclist[14]
        # target right gun left
        self.diagonal_trgl = self.mcclist[15]
  
    def analyze_center(self) -> dict:
        """
        calculates the results for the center profiles (crossplane/inplane) 
        using the calc_results function from the corresponding profile.

        Returns
        -------
        dict
            dict with results from profile analysis.

        """
        mcc_dict = {
            "CROSSPLANE_PROFILE": self.center_crossplane.calc_results(),
            "INPLANE_PROFILE": self.center_inplane.calc_results(),
            }
        
        return mcc_dict
    
    def analyze_diagonal(self) -> dict:
        """
        calculates the results for the diagonal profiles (tlgr/trgl) 
        using the calc_results function from the corresponding profile.

        Returns
        -------
        dict
            dict with results from profile analysis.

        """
        mcc_dict = {
            "TLGR_PROFILE": self.diagonal_tlgr.calc_results(),
            "TRGL_PROFILE": self.diagonal_trgl.calc_results(),
            }
        
        return mcc_dict     
        