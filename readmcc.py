# -*- coding: utf-8 -*-
"""
readmcc.py
Created on Fri Aug 01 13:25:39 2021

@author: bez0t
"""


import re
import numpy as np
import pandas as pd
from pymcc.wtscans import XyProfile, PDD


def read_file(filepath: str) -> list:
    """Read a mcc file and create an object from each data part.
    To access the filename in QATrack+ use FILE.name
    """
    # default values for our water tank scans
    nominal_fs = 100.0
    scan_depth = 100.0
    ssd = 900.0
    offset = 0.0
    filter = "FF"
    modality = "X"

    # list of mcc data objects
    data_obj = []

    with open(filepath) as meas_file:  # elegant way to parse file

        # empty list for data lines "Position\t\tMeasurementValue\t\tReference"
        lines = []
        # mark for data values, they get only edded if copy = true
        copy = False
        data_type = ""

        for line in meas_file:  # QATrack+ needs "in FILE"
            # strip() removes whitespace characters bevor and after text
            line = line.strip()
            # detect the profile type
            if line.split('=')[0] == "SCAN_CURVETYPE":
                data_type = line.split('=')[1]

            # elektrons or photons?
            if line.split('=')[0] == "MODALITY":
                modality = line.split('=')[1]

            # get field offset information
            if line.split('=')[0] == "COLL_OFFSET_INPLANE":
                offset = float(line.split('=')[1])
            if line.split('=')[0] == "COLL_OFFSET_CROSSPLANE":
                offset = float(line.split('=')[1])

            # get nominal field size
            if line.split('=')[0] == "FIELD_INPLANE":
                nominal_fs = float(line.split('=')[1])

            if line.split('=')[0] == "FILTER":
                filter = line.split('=')[1]

            # get geometry information
            if line.split("=")[0] == "ISOCENTER":
                isocenter = float(line.split("=")[1])
            if line.split("=")[0] == "SSD":
                ssd = float(line.split("=")[1])
            if line.split("=")[0] == "SCAN_DEPTH":
                scan_depth = float(line.split("=")[1])

            # find the data block
            if line == "BEGIN_DATA":
                copy = True
            elif line == "END_DATA":
                copy = False

                if data_type == "INPLANE_PROFILE":
                    data = conv_data(lines)
                    data_obj.append(XyProfile(modality, data_type, offset, 
                                             nominal_fs, filter, isocenter, 
                                             ssd, scan_depth, data))
                    lines = [] # empty line buffer
                if data_type == "CROSSPLANE_PROFILE":
                    data = conv_data(lines)
                    data_obj.append(XyProfile(modality, data_type, offset, 
                                             nominal_fs, filter, isocenter, 
                                             ssd, scan_depth, data))
                    lines = [] # empty line buffer
                if data_type == "PDD":
                    data = conv_data(lines)
                    data_obj.append(PDD(modality, data_type, offset, nominal_fs,
                        filter, isocenter, ssd, scan_depth, data))
                    lines = [] # empty line buffer

            elif copy:
                lines.append(line)

    return data_obj


def conv_data(lines: list) -> pd.DataFrame:
    """Function that converts text lines to numeric data with regular
    expressions and returns a DataFrame
    """

    # \s for unicode (str) patterns, matches whitespace characters
    # \d for unicode (str) patterns, matches decimal digit
    ptw_pattern = re.compile(r'(?P<position>\S+)\s{2}(?P<meas_values>\S+)\s{,2}#*(?P<reference>\S*)')
    data = []
    for line in lines:
        match = ptw_pattern.search(line)
        if match is not None:
            # Create Dict from match and append to data
            data_line = match.groupdict()
            data.append(data_line)

    data = pd.DataFrame(data)  # create DataFrame from list of Dicts
    data = data.apply(pd.to_numeric)  # Convert Strings to np.floats

    return linearize_data(data)


def linearize_data(dataframe: pd.DataFrame) -> pd.DataFrame:
    """Linearize and upsample the measurement data"""

    # ndarray with 10th mm resolution between first an last measurement
    # position, beginning at the first index (could be smaller than 0!)
    start = dataframe.position[dataframe.first_valid_index()] * 10
    # last index
    stop = dataframe.position[dataframe.last_valid_index()] * 10

    interp_arr = np.arange(start, stop)

    # create dict with single entry to construct dataframe with column names
    dfinter = {"position": interp_arr / 10}
    # convert dict to 1-D dataframe
    dataint = pd.DataFrame(dfinter)

    # connect data frames
    dataint = pd.merge_ordered(dataint, dataframe, on="position",
                                how="outer")
    dataint.meas_values = dataint.meas_values.interpolate()
    #dataint.meas_values = dataint.meas_values.interpolate(method='cubic')

    return dataint

# if __name__ == "__main__":
#     import sys
#     #print(sys.argv[1])
#     mymcc = read_file(sys.argv[1])
#     for i in mymcc:
#         if i.curve_type == "PDD":
#             print(i.calc_pdd())
#         elif i.curve_type == "INPLANE_PROFILE":
#             print(i.calc_profile())
#         elif i.curve_type == "CROSSPLANE_PROFILE":
#             print(i.calc_profile())