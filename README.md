# pymcc
Module that reads mephisto mcc files from watertank scans or Array files. pymcc relies on Pandas and uses a Pandas DataFrame to store the measurement values inside the class objects. 

At the moment the tests implemented are the ones used at the University Medical Center Mainz (mainly the default PTW Data Analyze Varian profile ones). pymcc is able to handle Photon (FF and FFF) as well as electron water tank measurements and Starcheck files (at the moment a list of XyProfiles). Octavius 729 measurement files can be read and the measurement is put into a interpolated Pandas DataFrame (No further analysis possible).

pymcc was created to analyze mcc files from within QATrack+ through an upload test that creates a dict (of dicts), followed by composite tests that grab the individual results from that dict.

Only symmetric fields are tested, asymmetric fields should work.

## Examples for QATrack+

Using pymcc.wtscans profile and pdd curves can be analyzed. 

### Water Tank Profile
*File Upload Test: read_mcc_6x_10x10*
```Python
import pymcc

# read mcc file and return list of measurement objects (PDD and/or Profiles)
mymcc = pymcc.readmcc.read_file(FILE.name)

mcc_dict = {}
for i in mymcc:
    mcc_dict[i.curve_type] = i.calc_results()

# provide object dict for composite tests
read_mcc_6x_10x10 = mcc_dict
```

*Access single result from dict in composite test: flat_10x10_x_py (makro name)*
```Python
flat_10x10_x_py = read_mcc_6x_10x10["CROSSPLANE_PROFILE"]["Flatness"]
```

A complete mcc_dict can looks like this (here electrons, 20 x 20 cm²):

{'PDD': {'Type': 'PDD', 'R80': 19.986482881280565, 'R50': {'R50 (DIN)': 24.940544317729188, 'R50': 23.788279990789775}, 'Rp': 29.794591153821518}, 
'INPLANE_PROFILE': {'Type': 'INPLANE_PROFILE', 'CaxDev': 0.3168300690775254, 'FWHM': {'fwhm (nominal)': 204.02442621763075, 'fwhm': 201.42603042514637}, 'Flatness': 1.9405358627975489, 'Symmetry': 0.9219008965810869}, 
'CROSSPLANE_PROFILE': {'Type': 'CROSSPLANE_PROFILE', 'CaxDev': 0.49851535225756294, 'FWHM': {'fwhm (nominal)': 203.39409918216086, 'fwhm': 200.8037310515953}, 'Flatness': 1.6065593462425587, 'Symmetry': -0.2758778417071142}}


### Starcheck
(without using the starcheck class)
*File Upload Test: read_mcc_0_6x*
```Python
import pymcc

# read mcc file and return list of measurement objects (PDD and/or Profiles)
mymcc = pymcc.readmcc.read_file(FILE.name)

mcc_dict = {
    "CROSSPLANE_PROFILE": mymcc[3].calc_results(),
    "INPLANE_PROFILE": mymcc[10].calc_results(),
}

# provide object dict for composite tests
read_mcc_0_6x = mcc_dict
```
*Access single result from dict in composite test: fs_0_y_py (macro name)*
```Python
fs_0_y_py = read_mcc_0_6x["INPLANE_PROFILE"]["FWHM"]["fwhm"]
```

The complete dict for a Starcheck mcc file contains 16 profiles of variable length. #3 and # 10 are crossplane and inplane profiles through the central axis. #14 and #15 are the diagonal profiles.
