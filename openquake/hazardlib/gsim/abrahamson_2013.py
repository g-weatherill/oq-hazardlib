# The Hazard Library
# Copyright (C) 2013 GEM Foundation
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
Module exports :class:`AbrahamsonEtAl2013`.
"""

from __future__ import division

import numpy as np

from openquake.hazardlib.gsim.base import GMPE, CoeffsTable
from openquake.hazardlib import const
from openquake.hazardlib.imt import PGA, SA


class AbrahamsonEtAl2013Inter(GMPE):
    """
    Implements the Subduction GMPE developed by Norman Abrahamson, Nicholas
    Gregor and Kofi Addo, otherwise known as the "BC Hydro" Model, published
    as "BC Hydro Ground Motion Prediction Equations For Subduction Earthquakes
    (2013, Earthquake Spectra, in press).
    This implements only the interface GMPE for forearc sites
    """

    #: Supported tectonic region type is subduction interface
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.SUBDUCTION_INTERFACE

    #: Supported intensity measure types are spectral acceleration,
    #: and peak ground acceleration, see table 1, page 1715
    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([
        PGA,
        SA
    ])

    #: Supported intensity measure component is the random horizontal
    #component :
    #attr:`~openquake.hazardlib.const.IMC.GEOMETRIC_MEAN`, see
    #paragraph 'Functional : Form', page 1706
    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.GEOMETRIC_MEAN

    #: Supported standard deviation types are inter-event, intra-event
    #: and total, see table 1, page 1715
    DEFINED_FOR_STANDARD_DEVIATION_TYPES = set([
        const.StdDev.TOTAL,
        const.StdDev.INTER_EVENT,
        const.StdDev.INTRA_EVENT,
        const.StdDev.SINGLE_STATION
    ])

    #: Required site parameters is Vs30, used to distinguish between NEHRP
    #: soil classes, see paragraph 'Functional Form', page 1706
    #: For the Abrahamson et al (2013) GMPE a new term is introduced to 
    #: determine whether a site is on the forearc with respect to the 
    #: subduction interface, or on the backarc. This boolean is a vector 
    #: containing True for a forearc site or False for a backarc site.

    REQUIRES_SITES_PARAMETERS = set(('vs30', 'forearc'))

    #: Required rupture parameters are magnitude and focal depth, see equation
    #: 1, page 1706
    REQUIRES_RUPTURE_PARAMETERS = set(('mag'))

    #: Required distance measure is closest distance to rupture, for 
    #: interface events
    REQUIRES_DISTANCES = set(('rrup'))
    

    def get_mean_and_stddevs(self, sites, rup, dists, imt, stddev_types):
        """
        See :meth:`superclass method
        <.base.GroundShakingIntensityModel.get_mean_and_stddevs>`
        for spec of input and result values.
        """
        # extract dictionaries of coefficients specific to required
        # intensity measure type and for PGA
        C = self.COEFFS[imt]
        C = self._get_delta_c1_coeffs(C)
        C_PGA = self.COEFFS[PGA()]

        # compute median pga on rock (vs30=1100), needed for site response
        # term calculation
        pga1000 = np.exp(self._compute_imt1000(C_PGA, sites, rup, dists))
        mean = (self._compute_base_term(C) +
                self._compute_magnitude_term(C, rup.mag) +
                self._compute_distance_term(C, rup.mag, dists) + 
                self._compute_focal_depth_term(C, rup) + 
                self._compute_forearc_backarc_term(C, sites, dists) +
                self._compute_site_response_term(C, sites, pga1000))

        stddevs = self._getstddevs(C, stddev_types, len(sites.vs30))
        return mean, stddevs

        
    def _get_delta_c1_coeffs(self, C):
        """
        Returns the choice of coefficients to use for Delta C1
        The original manuscript provides three sets of coefficients for the 
        large magnitude scaling term DeltaC_1. These are intended to represent 
        the epistemic uncertainty in the large magnitude scaling with 
        a "lower", "central" and "upper" set of coefficients. See Table 8.
        As the DeltaC_1 coefficients were only given for five spectral periods
        and PGA in the original manuscript the values for intermediate periods
        were derived by linear interpolation. For periods > 3.0 s the scaling
        values for Sa (3.0) were assumed.
        """
        C['delta_c1'] = C['dc1_cent']
        return C

    def _compute_base_term(self, C):
        """
        Computes the base model term described in equation (3)
        """
        return C['theta1'] + (C['theta4'] * C['delta_c1'])


    def _compute_focal_depth_term(self, C, rup):
        """
        Computes the hypocentral depth scaling term - as indicated by 
        equation (5)
        """
        # For interface events F_EVENT = 0.. so no depth scaling is returned
        return 0.


    def _compute_forearc_backarc_term(self, C, sites, dists):
        """
        Computes the forearc/backarc scaling term given by equation (6).
        """
        f_faba = np.zeros_like(dists, rrup)
        # Term only applies to backarc sites (F_FABA = 0. for forearc)
        idx = np.logical_not(sites.forearc)

        max_dist = dists.rrup[idx]
        max_dist[max_dist < 100.0] = 100.0
        f_faba[idx] = C['theta15'] + C['theta16'] * np.log(max_dist / 40.0)
        return f_faba


    def _compute_magnitude_term(self, C, mag):
        """
        Computes the magnitude scaling term given by equation (4). 
        """

        if mag > C['dc1_cent']
            f_mag = self.CONSTS['theta4'] * (mag - (7.8 + C['delta_c1'])) + \
                C['theta13'] * ((10. - mag) ** 2.)
        else:
            f_mag = self.CONSTS['theta5'] * (mag - (7.8 + C['delta_c1'])) + \
                C['theta13'] * ((10. - mag) ** 2.)
        return f_mag


    def _compute_distance_term(self, C, mag, dists):
        """
        Computes the distance scaling term, as contained within equation (3)
        """
        return (C['theta2'] + self.CONSTS['theta3'] * (mag - 7.8)) * \
            np.log(dists.rrup + self.CONSTS['c4'] * np.exp((mag - 6.) * 
            self.CONSTS['theta9'])) + (C['theta6'] * dists.rrup)


    def _compute_site_response_term(self, C, sites, pga1000):
        """
        Compute and return site response model term, the Abrahamson et al 
        (2013) GMPE adopts the same site response scaling model of 
        Walling et al (2008) as implemented in the Abrahamson & Silva (2008)
        GMPE. The functional form is retained here.
        """
        site_resp_term = np.zeros_like(sites.vs30)

        vs_star = sites.vs30.copy()
        vs_star[vs_star > 1000.0] = 1000.
        
        vlin, c, n = C['vlin'], self.CONSTS['c'], self.CONSTS['n']
        theta12, b = C['theta12'], C['b']

        idx = sites.vs30 < vlin
        arg = vs_star[idx] / vlin
        site_resp_term[idx] = (theta12 * np.log(arg) -
                               b * np.log(pga1000[idx] + c) +
                               b * np.log(pga1000[idx] + c * (arg ** n)))

        idx = sites.vs30 >= vlin
        site_resp_term[idx] = (theta12 + b * n) * np.log(vs_star[idx] / vlin)
        return site_resp_term


    def _compute_imt_1000(self, C, sites, rup, dists):
        """
        Compute and return mean imt value for rock conditions
        (vs30 = 1100 m/s)
        """
        mean = (self._compute_base_term(C) +
                self._compute_magnitude_term(C, rup.mag) +
                self._compute_distance_term(C, rup.mag, dists) + 
                self._compute_focal_depth_term(C, rup) + 
                self._compute_forearc_backarc_term(C, sites, dists))
        site_response = (C['theta12'] + C['b'] * C['n']) * 
            np.log(1000. / C['vlin'])
        return mean + site_response
        

    def _get_stddevs(self, C, stddev_types, num_sites):
        """
        Return standard deviations as defined in table 5, p. 200.
        """
        stddevs = []
        for stddev_type in stddev_types:
            assert stddev_type in self.DEFINED_FOR_STANDARD_DEVIATION_TYPES
            if stddev_type == const.StdDev.TOTAL:
                stddevs.append(C['sigma'] + np.zeros(num_sites))
            elif stddev_type == const.StdDev.INTRA_EVENT:
                stddevs.append(C['tau'] + np.zeros(num_sites))
            elif stddev_type == const.StdDev.INTER_EVENT:
                stddevs.append(C['phi'] + np.zeros(num_sites))
            elif stddev_type == const.StdDev.SINGLE_STATION:
                stddevs.append(C['sigma_ss'] + np.zeros(num_sites))
        return stddevs


    # Period-dependent coefficients (Table 5)
    COEFFS_SINTER=CoeffsTable(sa_damping5, table="""\
imt          vlin        b   theta1    theta2    theta6   theta7    theta8  theta10  theta11   theta12   theta13   theta14  theta15   theta16      phi     tau   sigma  sigma_ss  dc1_low  dc1_cent  dc1_high
pga      865.1000  -1.1860   4.2203   -1.3500   -0.0012   1.0988   -1.4200   3.1200   0.0130    0.9800   -0.0135   -0.4000   0.9969   -1.0000   0.6000  0.4300  0.7400    0.6000   0.0000    0.2000    0.4000
0.0200   865.1000  -1.1860   4.2203   -1.3500   -0.0012   1.0988   -1.4200   3.1200   0.0130    0.9800   -0.0135   -0.4000   0.9969   -1.0000   0.6000  0.4300  0.7400    0.6000   0.0000    0.2000    0.4000
0.0500  1053.5000  -1.3460   4.5371   -1.4000   -0.0012   1.2536   -1.6500   3.3700   0.0130    1.2880   -0.0138   -0.4000   1.1030   -1.1800   0.6000  0.4300  0.7400    0.6000   0.0000    0.2000    0.4000
0.0750  1085.7000  -1.4710   5.0733   -1.4500   -0.0012   1.4175   -1.8000   3.3700   0.0130    1.4830   -0.0142   -0.4000   1.2732   -1.3600   0.6000  0.4300  0.7400    0.6000   0.0000    0.2000    0.4000
0.1000  1032.5000  -1.6240   5.2892   -1.4500   -0.0012   1.3997   -1.8000   3.3300   0.0130    1.6130   -0.0145   -0.4000   1.3042   -1.3600   0.6000  0.4300  0.7400    0.6000   0.0000    0.2000    0.4000
0.1500   877.6000  -1.9310   5.4563   -1.4500   -0.0014   1.3582   -1.6900   3.2500   0.0130    1.8820   -0.0153   -0.4000   1.2600   -1.3000   0.6000  0.4300  0.7400    0.6000   0.0000    0.2000    0.4000
0.2000   748.2000  -2.1880   5.2684   -1.4000   -0.0018   1.1648   -1.4900   3.0300   0.0129    2.0760   -0.0162   -0.3500   1.2230   -1.2500   0.6000  0.4300  0.7400    0.6000   0.0000    0.2000    0.4000
0.2500   654.3000  -2.3810   5.0594   -1.3500   -0.0023   0.9940   -1.3000   2.8000   0.0129    2.2480   -0.0172   -0.3100   1.1600   -1.1700   0.6000  0.4300  0.7400    0.6000   0.0000    0.2000    0.4000
0.3000   587.1000  -2.5180   4.7945   -1.2800   -0.0027   0.8821   -1.1800   2.5900   0.0128    2.3480   -0.0183   -0.2800   1.0500   -1.0600   0.6000  0.4300  0.7400    0.6000   0.0000    0.2000    0.4000
0.4000   503.0000  -2.6570   4.4644   -1.1800   -0.0035   0.7046   -0.9800   2.2000   0.0127    2.4270   -0.0206   -0.2300   0.8000   -0.7800   0.6000  0.4300  0.7400    0.6000  -0.0500    0.1500    0.3500
0.5000   456.6000  -2.6690   4.0181   -1.0800   -0.0044   0.5799   -0.8200   1.9200   0.0125    2.3990   -0.0231   -0.1900   0.6620   -0.6200   0.6000  0.4300  0.7400    0.6000  -0.1000    0.1000    0.3000
0.6000   430.3000  -2.5990   3.6055   -0.9900   -0.0050   0.5021   -0.7000   1.7000   0.0124    2.2730   -0.0256   -0.1600   0.5800   -0.5000   0.6000  0.4300  0.7400    0.6000  -0.1200    0.0800    0.2800
0.7500   410.5000  -2.4010   3.2174   -0.9100   -0.0058   0.3687   -0.5400   1.4200   0.0120    1.9930   -0.0296   -0.1200   0.4800   -0.3400   0.6000  0.4300  0.7400    0.6000  -0.1500    0.0500    0.2500
1.0000   400.0000  -1.9550   2.7981   -0.8500   -0.0062   0.1746   -0.3400   1.1000   0.0114    1.4700   -0.0363   -0.0700   0.3300   -0.1400   0.6000  0.4300  0.7400    0.6000  -0.2000    0.0000    0.2000
1.5000   400.0000  -1.0250   2.0123   -0.7700   -0.0064  -0.0820   -0.0500   0.7000   0.0100    0.4080   -0.0493    0.0000   0.3100    0.0000   0.6000  0.4300  0.7400    0.6000  -0.2500   -0.0500    0.1500
2.0000   400.0000  -0.2990   1.4128   -0.7100   -0.0064  -0.2821    0.1200   0.7000   0.0085   -0.4010   -0.0610    0.0000   0.3000    0.0000   0.6000  0.4300  0.7400    0.6000  -0.3000   -0.1000    0.1000
2.5000   400.0000   0.0000   0.9976   -0.6700   -0.0064  -0.4108    0.2500   0.7000   0.0069   -0.7230   -0.0711    0.0000   0.3000    0.0000   0.6000  0.4300  0.7400    0.6000  -0.3500   -0.1500    0.0500
3.0000   400.0000   0.0000   0.6443   -0.6400   -0.0064  -0.4466    0.3000   0.7000   0.0054   -0.6730   -0.0798    0.0000   0.3000    0.0000   0.6000  0.4300  0.7400    0.6000  -0.4000   -0.2000    0.0000
4.0000   400.0000   0.0000   0.0657   -0.5800   -0.0064  -0.4344    0.3000   0.7000   0.0027   -0.6270   -0.0935    0.0000   0.3000    0.0000   0.6000  0.4300  0.7400    0.6000  -0.4000   -0.2000    0.0000
5.0000   400.0000   0.0000  -0.4624   -0.5400   -0.0064  -0.4368    0.3000   0.7000   0.0005   -0.5960   -0.0980    0.0000   0.3000    0.0000   0.6000  0.4300  0.7400    0.6000  -0.4000   -0.2000    0.0000
6.0000   400.0000   0.0000  -0.9809   -0.5000   -0.0064  -0.4586    0.3000   0.7000  -0.0013   -0.5660   -0.0980    0.0000   0.3000    0.0000   0.6000  0.4300  0.7400    0.6000  -0.4000   -0.2000    0.0000
7.5000   400.0000   0.0000  -1.6017   -0.4600   -0.0064  -0.4433    0.3000   0.7000  -0.0033   -0.5280   -0.0980    0.0000   0.3000    0.0000   0.6000  0.4300  0.7400    0.6000  -0.4000   -0.2000    0.0000
10.0000  400.0000   0.0000  -2.2937   -0.4000   -0.0064  -0.4828    0.3000   0.7000  -0.0060   -0.5040   -0.0980    0.0000   0.3000    0.0000   0.6000  0.4300  0.7400    0.6000  -0.4000   -0.2000    0.0000""")

    CONSTS = {
        # Period-Independent Coefficients (Table 4)
        'n': 1.18,
        'c': 1.88,
        'theta3': 0.1,
        'theta4': 0.9,
        'theta5': 0.0,
        'theta9': 0.4,
        'c4': 10.0
        }


class AbrahamsonEtAl2013InterHigh(AbrahamsonEtAl2013Inter):
    """
    Defines the Abrahamson et al. (2013) scaling relation  assuming the upper
    values of the magnitude scaling for large slab earthquakes, as defined in
    table 8
    """

    def _get_delta_c1_coeffs(self, C):
        """
        Returns the choice of coefficients to use for Delta C1
        """
        C['delta_c1'] = C['dc1_high']
        return C


class AbrahamsonEtAl2013InterLow(AbrahamsonEtAl2013Inter):
    '''
    Defines the Abrahamson et al. (2013) scaling relation  assuming the upper
    values of the magnitude scaling for large slab earthquakes, as defined in
    table 8
    '''
    def _get_delta_c1_coeffs(self, C):
        """
        Returns the choice of coefficients to use for Delta C1
        """
        C['delta_c1'] = C['dc1_low']
        return C


class AbrahamsonEtAl2013SSlab(AbrahamsonEtAl2013Inter):
    """
    Implements the Subduction GMPE developed by Norman Abrahamson, Nicholas
    Gregor and Kofi Addo, otherwise known as the "BC Hydro" Model, published
    as "BC Hydro Ground Motion Prediction Equations For Subduction Earthquakes
    (2013, Earthquake Spectra, in press).
    This implements only the inslab GMPE. For inslab events the source is
    considered to be a point source located at the hypocentre. Therefore 
    the hypocentral distance metric is used in place of the rupture distance,
    and the hypocentral depth is used to scale the ground motion by depth
    """

    #: Required distance measure is hypocentral for in-slab events
    REQUIRES_DISTANCES = set(('rhypo'))
    #: In-slab events require constraint of hypocentral depth
    REQUIRES_RUPTURE_PARAMETERS = set(('mag', 'hypo_depth'))


    def _compute_base_term(self, C):
        """
        Computes the base model term described in equation (3)
        """
        return C['theta1'] + (C['theta4'] * C['delta_c1']) + C['theta10']
                                                             

    def _compute_focal_depth_term(self, C, rup):
        """
        Computes the hypocentral depth scaling term - as indicated by 
        equation (5)
        """
        return C['theta11'] * (rup.hypo_depth - 60.)

    
    def _compute_distance_term(self, C, mag, dists):
        """
        Computes the distance scaling term, as contained within equation (3)
        """
        return (C['theta2'] + C['theta14'] + self.CONSTS['theta3'] * 
            (mag - 7.8)) * np.log(dists.rhypo + self.CONSTS['c4'] * 
            np.exp((mag - 6.) * self.CONSTS['theta9'])) + \
            (C['theta6'] * dists.rrup)


    def _compute_forearc_backarc_term(self, C, sites, dists):
        """
        Computes the forearc/backarc scaling term given by equation (6).
        """
        f_faba = np.zeros_like(dists.rhypo)
        # Term only applies to backarc sites (F_FABA = 0. for forearc)
        idx = np.logical_not(sites.forearc)

        max_dist = dists.rhypo[idx]
        max_dist[max_dist < 85.0] = 85.0
        f_faba[idx] = C['theta7'] + C['theta8'] * np.log(max_dist / 40.0)
        return f_faba


class AbrahamsonEtAl2013SSlabHigh(AbrahamsonEtAl2013SSlab):
    """
    Defines the Abrahamson et al. (2013) scaling relation  assuming the upper
    values of the magnitude scaling for large slab earthquakes, as defined in
    table 8
    """
    def _get_delta_c1_coeffs(self, C):
        """
        Returns the choice of coefficients to use for Delta C1
        """
        C['delta_c1'] = C['dc1_high']
        return C


class AbrahamsonEtAl2013SSlabLow(AbrahamsonEtAl2013SSlab):
    """
    Defines the Abrahamson et al. (2013) scaling relation  assuming the upper
    values of the magnitude scaling for large slab earthquakes, as defined in
    table 8
    """
    def _get_delta_c1_coeffs(self, C):
        """
        Returns the choice of coefficients to use for Delta C1
        """
        C['delta_c1'] = C['dc1_low']
        return C

