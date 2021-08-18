# -*- coding: utf-8 -*-
"""
wtscans.py
Created on Fri Aug 13 13:20:39 2021

@author: bez0t
"""

import numpy as np
import pandas as pd
import pylinac
from scipy.stats import linregress


class PDD:
    """
    Calculate the field parameters asked for in DIN6847-5
    """
    def __init__(self, mod: str, curve_type: str, offset: float, 
                 nominal_fs: float, filter: str, isocenter: float, 
                 ssd: float, scan_depth: float, datafr: pd.DataFrame()) -> None:

        self.curve_type = curve_type
        self.dataframe = datafr
        self.isocenter = isocenter
        self.scan_depth = scan_depth
        self.offset = offset
        self.ssd = ssd
        self.nominal_fs = nominal_fs
        self.filter = filter
        self.modality = mod


    @staticmethod
    def interp_value(point_a: np.float64, point_b: np.float64,
                     interp_x: np.float64) -> np.float64:
        """Calculate an interpolated value between two points."""

        slope = (point_b[1]-point_a[1]) / (point_b[0]-point_a[0])

        intersect = point_a[1] - slope*point_a[0]

        # y value at x
        interp_y = (interp_x-intersect) / slope

        return interp_y
    

    def depth_max(self, interp=True) -> np.float64:
        """function that returns the position of the maximum of the pdd
        """
        if interp:
                wl = int(self.dataframe.position.size / 100)
                rol = self.dataframe.rolling(window=wl, min_periods=None, center=True)
                max_depth = self.dataframe.position[rol.meas_values.mean().idxmax()]
        else:
                idx_max = self.dataframe['meas_values'].idxmax()
                max_depth = self.dataframe.position.iat[idx_max]

        return max_depth


    def depth_x(self, x: np.float) -> np.float64:
        """function that returns the position of x% of the max of the pdd with
        linear interpolation between data points (in case of a steep gradient)
        """

        max_dose = self.dataframe['meas_values'].max()
        x_val = max_dose * x/100
        # searchsorted is a binary search and expects a sorted series
        # hence we need invert the curve.
        idx_x_depth = self.dataframe.iloc[::-1].meas_values.searchsorted(x_val)

        b_1 = (self.dataframe.iloc[::-1].position.iat[idx_x_depth-1],
               self.dataframe.iloc[::-1].meas_values.iat[idx_x_depth-1])
        b_2 = (self.dataframe.iloc[::-1].position.iat[idx_x_depth],
               self.dataframe.iloc[::-1].meas_values.iat[idx_x_depth])

        x_depth = self.interp_value(b_1, b_2, x_val)

        # also the returned index is from the inverted curve...
        #x_depth = self.dataframe.iloc[::-1].position.iat[idx_x_depth]

        return x_depth


    def dose_x(self, depth: np.float) -> np.float64:
        """ returns the percent dose value at x mm depth
        """

        idx = self.dataframe.position.searchsorted(depth)
        max_dose = self.dataframe['meas_values'].max()
        dx = self.dataframe.meas_values.iat[idx] / max_dose

        return dx

    def dose_100(self) -> np.float64:
        """ returns the percent dose value at 100 mm depth
        """

        idx = self.dataframe.position.searchsorted(100.0)
        max_dose = self.dataframe['meas_values'].max()
        d100 = self.dataframe.meas_values.iat[idx] / max_dose

        return d100


    def dose_200(self) -> np.float64:
        """ returns the percent dose value at 200 mm depth
        """

        idx = self.dataframe.position.searchsorted(200.0)
        max_dose = self.dataframe['meas_values'].max()
        d200 = self.dataframe.meas_values.iat[idx] / max_dose

        return d200


    def calc_R50_din(self) -> np.float64:
        """function that returns the R50 value according to german DIN 6800-2
        """
        # half maximum
        D50ion = self.dataframe['meas_values'].max() / 2

        # searchsorted is a binary search and expects a sorted series
        # hence we need invert the curve.
        idx_D50ion = self.dataframe.iloc[::-1].meas_values.searchsorted(D50ion)

        b_1 = (self.dataframe.iloc[::-1].position.iat[idx_D50ion-1],
               self.dataframe.iloc[::-1].meas_values.iat[idx_D50ion-1])
        b_2 = (self.dataframe.iloc[::-1].position.iat[idx_D50ion],
               self.dataframe.iloc[::-1].meas_values.iat[idx_D50ion])

        R50ion = self.interp_value(b_1, b_2, D50ion)

        R50 = (0.00171 * (R50ion*R50ion) + 1.00805 * R50ion - 0.00689)

        return R50


    def calc_surface_dose(self) -> np.float64:
        """function that returns the relative surface dose at 0.5 mm depth
        """
        
        # only works if values are already interpolated with stepsize 0.1 mm
        surface_dose = self.dataframe.meas_values.iat[5]
        max_dose = self.dataframe['meas_values'].max()

        rel_surface_dose = surface_dose / max_dose * 100

        return rel_surface_dose


    def calc_q_index(self) -> np.float64:
        """calculate Q_index following IAEA TRS 398
        """

        d100 = self.dose_100()
        d200 = self.dose_200()

        qi = (1.2661 * d200 / d100) - 0.0595

        return qi


    def calc_E0_mean(self) -> np.float64:
        """ Estimate of the mean Electron energy at the phantom surface
        """

        E0_mean = 2.33 * self.calc_R50_din() / 10

        return E0_mean


    def calc_Rp(self) -> np.float64:
        """Calculate practical Range of an electron beam"""

        max_dose = self.dataframe['meas_values'].max()
        
        # get 40/60% data points (x,y) for slope calculation
        a1 = (self.depth_x(60.0), 0.6 * max_dose)
        a2 = (self.depth_x(40.0), 0.4 * max_dose)

        # calculate slope:
        slope_50 = (a1[1]-a2[1])/(a1[0]-a2[0])
        
        # intercept
        a50 = (self.depth_x(50.0), 0.5 * max_dose)
        inter_50 = a50[1]- (a50[0] * slope_50)

        # define data range for linear regression, start with an estimated
        # practical Range (formula from PTW data analyze handbook)
        E0_mean = self.calc_E0_mean()
        Rp_est = (0.11 + 0.505 * E0_mean - 3E-4 * E0_mean**2)*100
        # index equals steps of 0.1 mm
                
        # add saftey distance behind Rp (slope should be shallow here)
        lin_start = int(Rp_est + 100)
        
        print(lin_start)
        print(self.dataframe.position.size - 2)
        
        # if there are enough data points do a linear regression
        if lin_start < ( self.dataframe.position.size - 2 ):
            #self.dataframe[lin_start:].meas_values.plot()
        
            # do linear regression on Dataframe Series
            lr = linregress(self.dataframe.position[lin_start:], 
                        self.dataframe.meas_values[lin_start:])
            inter_bs = lr.intercept
            slope_bs = lr.slope
        # else, use the last point as an approximation for the intercept
        # with 0 slope
        else:
            inter_bs = self.dataframe['meas_values'].iloc[-1]
            slope_bs = 0
                       
        # Rp is the depth of the point of intersection 
        Rp = (inter_bs - inter_50) / (slope_50 - slope_bs)
        
        return Rp
    
    
    def calc_results(self) -> dict:
        """Calculate all relevant values for inplane and crossplane data and
        return two dicts, one for inplane and one for crossplane.
        """
        if self.modality == "EL":
            results = {
                "Type": self.curve_type,
                "R80": self.depth_x(80.0),
                "R50": {"R50 (DIN)": self.calc_R50_din(), 
                        "R50": self.depth_x(50.0)},
                "Rp": self.calc_Rp(),
                }
        elif self.modality == "X":
            results = {
                "Type": self.curve_type,
                "Q Index": self.calc_q_index(),
                "Surface Dose": self.calc_surface_dose(),
                "D100": self.dose_100(),
                "D200": self.dose_200()
                }
        else:
            results = {}
        
        
        return results    
    
    

class XyProfile:
    """
    Calculate the field parameters asked for in DIN6847-5
    """
    def __init__(self, mod: str, curve_type: str, offset: float, 
                 nominal_fs: float, filter: str, isocenter: float, 
                 ssd: float, scan_depth: float, datafr: pd.DataFrame()) -> None:

        self.dataframe = datafr
        self.offset = offset
        self.isocenter = isocenter
        self.scan_depth = scan_depth
        self.ssd = ssd
        self.nominal_fs = nominal_fs
        self.filter = filter
        self.curve_type = curve_type
        self.modality = mod



    @staticmethod
    def interp_value(point_a: np.float64, point_b: np.float64,
                     interp_x: np.float64) -> np.float64:
        """Calculate an interpolated value between two points."""

        slope = (point_b[1]-point_a[1]) / (point_b[0]-point_a[0])

        intersect = point_a[1] - slope*point_a[0]

        # y value at x
        interp_y = (interp_x-intersect) / slope

        return interp_y

    def calc_fwhm(self, max_type: str = 'cax') -> dict:
        """Calculate the FWHM from data in the DataFrame and return the
        Fieldsize nominal and at isocenter distance.
        """

        half_max = self.calc_halfmax(max_type=max_type)

        # find the index where half max should be inserted before to keep the
        # series sorted. Hence value at _pos is always higher.
        left_pos = self.dataframe.meas_values.searchsorted(half_max)
        a_1 = (self.dataframe.position.iat[left_pos-1],
               self.dataframe.meas_values.iat[left_pos-1])
        a_2 = (self.dataframe.position.iat[left_pos],
               self.dataframe.meas_values.iat[left_pos])

        right_pos = self.dataframe.iloc[::-1].meas_values.searchsorted(half_max)
        b_1 = (self.dataframe.iloc[::-1].position.iat[right_pos-1],
               self.dataframe.iloc[::-1].meas_values.iat[right_pos-1])
        b_2 = (self.dataframe.iloc[::-1].position.iat[right_pos],
               self.dataframe.iloc[::-1].meas_values.iat[right_pos])

        fwhm = self.interp_value(b_1, b_2, half_max) - self.interp_value(a_1, a_2, half_max)

        # correct for measurement depth so that Fieldsize is returned at iso
        iso_corr = self.isocenter/(self.ssd+self.scan_depth)

        fwhm_results = {'fwhm (nominal)': fwhm, 'fwhm': fwhm*iso_corr}

        return fwhm_results


    def calc_dinflat(self) -> np.float64:
        """Calculate the flatness according to Varian"""

        # extract np.array from Dataframe and create Pylinac profile
        messwerte = self.dataframe.meas_values.values
        #profil = pylinac.core.profile.SingleProfile(messwerte) #alte pylinac Ver.
        profil = pylinac.core.profile.SingleProfile(messwerte,
                    None,pylinac.core.profile.Interpolation.NONE,False,0.1,10,
                    pylinac.core.profile.Normalization.NONE)

        profil_max = profil.field_calculation(0.8, 'max')
        #print(\"Profil 80% Max: \", profil_max)
        profil_min = profil.field_calculation(0.8, 'min')
        #print(\"Profil 80% Min: \", profil_min)
        center = profil.geometric_center()
        #print(\"Messwert Zentralstrahl: \", profil.values[int(y_center['index (exact)'])])


        #Wert im Zentrahlstrahl
        din100 = profil.values[int(center['index (exact)'])]

        #Flatness nach PTW Data Analyse Varian Protocol
        din_flat = (profil_max /din100 - profil_min / din100) / 2 * 100

        #Geometric Center: profil.geometric_center()
        #Dmax: \", \"%.2f\" % (profil_max / din100 * 100), \"%\")
        #Dmin: \", \"%.2f\" % (profil_min / din100 * 100), \"%\")

        #Wert im Zentrahlstrahl (über dataframe, sollte mit pylianc ident. sein)
        #ptw100 = self.dataframe.loc[(self.dataframe.position == 0.0), 'meas_values'].iat[0]

        return din_flat


    def calc_fff_unflat(self) -> np.float64:
        """Calculate the FFF unflatness according to Fogliata"""

        field_width = self.calc_fwhm()
        idx_center = self.dataframe.position.searchsorted(0.0)

        # Unflattness
        idx_80_left = self.dataframe.position.searchsorted(0.0 - 0.8 * field_width['fwhm (nominal)']/2)
        idx_80_right = self.dataframe.position.searchsorted(0.8 * field_width['fwhm (nominal)']/2)

        unflattness = 1.0
        for x in [idx_80_left, idx_80_right]:
            tmp = (self.dataframe.meas_values.iat[idx_center] /
                                self.dataframe.meas_values.iat[x])
            if tmp > unflattness: unflattness = tmp

        return unflattness


    def calc_fff_slopes_peak(self) -> np.float64:
        """Calculate the left and right slopes of the fff profiles"""

        # fff field with needed to find 1/3 and 2/3 points on slopes
        field_width = self.calc_fwhm()

        # renormalize slope to CAX (%)
        idx_center = self.dataframe.position.searchsorted(0.0)
        renorm_percent = (self.calc_fff_renorm10x()*100 /
                            self.dataframe.meas_values.iat[idx_center])

        # find point positions (index) on left side
        idx_a1 = self.dataframe.position.searchsorted(0.0 - field_width['fwhm (nominal)']/3)
        idx_a2 = self.dataframe.position.searchsorted(0.0 - field_width['fwhm (nominal)']/6)

        # find point positions (index) on right side
        idx_b1 = self.dataframe.position.searchsorted(field_width['fwhm (nominal)']/6)
        idx_b2 = self.dataframe.position.searchsorted(field_width['fwhm (nominal)']/3)

        #print("CAX-Index: ", idx_center)

        # get data points for slope calculation (no interpolation)
        a1 = (self.dataframe.position.iat[idx_a1],
                self.dataframe.meas_values.iat[idx_a1]*renorm_percent)

        a2 = (self.dataframe.position.iat[idx_a2],
                self.dataframe.meas_values.iat[idx_a2]*renorm_percent)

        b1 = (self.dataframe.position.iat[idx_b1],
                self.dataframe.meas_values.iat[idx_b1]*renorm_percent)

        b2 = (self.dataframe.position.iat[idx_b2],
                self.dataframe.meas_values.iat[idx_b2]*renorm_percent)


        # calculate slopes:
        slope_left = (a1[1]-a2[1])/(a1[0]-a2[0])
        slope_right = (b1[1]-b2[1])/(b1[0]-b2[0])

        # Intercepts, Peak Position
        i_left = a1[1]- (a1[0] * slope_left)
        i_right = b2[1]- (b2[0] * slope_right)

        peak_pos = (i_left - i_right) / (slope_right - slope_left)

        return slope_left, slope_right, peak_pos


    def calc_caxdev_pylinac(self) -> np.float64:
        """Calculate the distance between the CAX and the center of the
        field that has been calculated. (maybe not perfect for large
        FFF fields)
        """

        messwerte = self.dataframe.meas_values.values
        profil = pylinac.core.profile.SingleProfile(messwerte,
                    None,pylinac.core.profile.Interpolation.NONE,False,0.1,10,
                    pylinac.core.profile.Normalization.NONE)
        # index der Position 0 ermitteln:
        #cax = self.dataframe.loc[self.dataframe.position == 0.0]
        #cax = cax.index[0]

        # CAX Abweichung berechnen aus der Mitte des Profils - Zentrahlstrahl
        cax_dev = (profil.beam_center()['index (exact)'] -
                                profil.geometric_center()['index (exact)']) / 10

        return cax_dev

    def calc_caxdev(self, max_type: str = 'cax') -> np.float64:
        """Calculate the distance between the CAX and the center of the
        field that has been calculated.
        """

        half_max = self.calc_halfmax(max_type=max_type)

        # find the index where half max should be inserted before to keep the
        # series sorted. Hence value at _pos is always higher.
        left_pos = self.dataframe.meas_values.searchsorted(half_max)
        a_1 = (self.dataframe.position.iat[left_pos-1],
               self.dataframe.meas_values.iat[left_pos-1])
        a_2 = (self.dataframe.position.iat[left_pos],
               self.dataframe.meas_values.iat[left_pos])

        right_pos = self.dataframe.iloc[::-1].meas_values.searchsorted(half_max)
        b_1 = (self.dataframe.iloc[::-1].position.iat[right_pos-1],
               self.dataframe.iloc[::-1].meas_values.iat[right_pos-1])
        b_2 = (self.dataframe.iloc[::-1].position.iat[right_pos],
               self.dataframe.iloc[::-1].meas_values.iat[right_pos])
        
        # for debugging
        # print(a_1, b_1)
        # print(self.interp_value(b_1, b_2, half_max))
        # print(self.interp_value(a_1, a_2, half_max))        
        
        cax_dev = (self.interp_value(b_1, b_2, half_max) +
                    self.interp_value(a_1, a_2, half_max)) / 2
        
        return cax_dev


    def calc_sym(self) -> np.float64:
        """Calculates the symmetry of the field plane with the point difference
        method. (same for FF and FFF? should the result be rescaled?).
        """
        messwerte = self.dataframe.meas_values.values
        profil = pylinac.core.profile.SingleProfile(messwerte,
                    None,pylinac.core.profile.Interpolation.NONE,False,0.1,10,
                    pylinac.core.profile.Normalization.NONE)
        symmetry = pylinac.field_analysis.symmetry_point_difference(profil, 0.8)
        return symmetry


    def calc_halfmax(self, max_type: str = 'cax') -> np.float64:
        """Return half of the max value or half the cax value depending on
        the given max_type. Default is to return half the CAX value. For FFF
        fields the value gets renormalized to return the correct value
        """

        if max_type == 'cax':
            max_val = self.dataframe.loc[(self.dataframe.position == self.offset),
                                    'meas_values'].iat[0]
        else:
            max_val = self.dataframe['meas_values'].max()

        half_max = (max_val / 2)

        if self.filter == "FFF":
            return half_max/self.calc_fff_renorm10x()

        return half_max


    def calc_fff_renorm10x(self) -> np.float64:
        """Return the field size correction factor for FFF beams
        Formula is taken from Folgiata et. al.
        (https://doi.org/10.1118/1.4754799) - valid for Truebeam 10 FFF
        """

        renorm_factor = (89.08 + 2.4826 * self.nominal_fs / 10
                         + 0.1152 * self.scan_depth / 10) \
                         / (1 - 0.0078 * self.nominal_fs / 10
                            + 0.0011 * self.scan_depth / 10)

        return renorm_factor/100

    def calc_fff_renorm6x(self) -> np.float64:
        """Return the field size correction factor for FFF beams
        Formula is taken from Folgiata et. al.
        (https://doi.org/10.1118/1.4754799) - valid for Truebeam 6 FFF
        """

        renorm_factor = (95.60 + 0.6595 * self.nominal_fs / 10
                         + 0.1255 * self.scan_depth / 10) \
                         / (1 - 0.0099 * self.nominal_fs / 10
                            + 0.0013 * self.scan_depth / 10)

        return renorm_factor/100

    def calc_results(self) -> dict:
        """Calculate all relevant values for inplane or crossplane data and
        return dict
        """
        if self.filter == "FF":
            results = {
                "Type": self.curve_type,
                "CaxDev": self.calc_caxdev(),
                "FWHM": self.calc_fwhm(),
                "Flatness": self.calc_dinflat(),
                "Symmetry": self.calc_sym()
                }
        elif self.filter == "FFF":
            results = {
                "Type": self.curve_type,
                "CaxDev": self.calc_caxdev(),
                "FWHM": self.calc_fwhm(),
                "Symmetry": self.calc_sym(),
                "Peak": self.calc_fff_slopes_peak()[2]
                }
        else:
            results ={}

        return results    