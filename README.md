# pymcc
Module that reads mephisto mcc files from watertank scans or array files. pymcc relies on Pandas and uses a Pandas DataFrame to store the measurement values inside the class objects. 

At the moment the tests implemented are the ones used at the University Medical Center Mainz (mainly the default PTW Data Analyze Varian profile ones). pymcc is able to handle Photon (FF and FFF) as well as electron water tank measurements, Starcheck (files with options to analyze the four main profiles) as well as Octavius 729 measurement files (measurement is put into a interpolated Pandas DataFrame - no further analysis possible).

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

For a single photon pdd the dict contains the following:
{'PDD': {'Type': 'PDD', 'Q Index': 0.6671976734598629, 'Surface Dose': 47.71301033726637, 'D100': 67.11913960530438, 'D200': 38.52406807977445}}

### Starcheck
*File Upload Test: read_mcc_0_6x*
```Python
import pymcc

# read Starcheck mcc file
mystar = pymcc.readmcc.read_file(FILE.name)
# analyze center inplane and crossplane profiles
mcc_dict = mystar.analyze_center()

# provide object dict for composite tests
read_mcc_0_6x = mcc_dict

```
*Access single result from dict in composite test: fs_0_y_py (macro name)*
```Python
fs_0_y_py = read_mcc_0_6x["INPLANE_PROFILE"]["FWHM"]["fwhm"]
```

The two diagonals are labeled TLGR_PROFILE and TRGL_PROFILE. The dict from mystar.analyze_diagonal() looks like this:

{'TLGR_PROFILE': {'Type': 'INPLANE_PROFILE', 'CaxDev': -1.7500384951353993, 'FWHM': {'fwhm (nominal)': 276.49675566173937, 'fwhm': 276.49675566173937}, 'Flatness': 2.118534052326748, 'Symmetry': -0.7116183813665132}, 'TRGL_PROFILE': {'Type': 'INPLANE_PROFILE', 'CaxDev': 0.12798456264559377, 'FWHM': {'fwhm (nominal)': 275.69602807162545, 'fwhm': 275.69602807162545}, 'Flatness': 2.20840822848013, 'Symmetry': -0.8215006262255542}}
