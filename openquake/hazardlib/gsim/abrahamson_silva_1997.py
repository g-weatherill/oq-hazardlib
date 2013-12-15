# Copyright (C) 2012 GEM Foundation
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
Module exports :class:`AbrahamsonSilva2008`.
"""
from __future__ import division

import numpy as np

from openquake.hazardlib.gsim.base import GMPE, CoeffsTable
from openquake.hazardlib import const
from openquake.hazardlib.imt import PGA, SA


class AbrahamsonSilva1997(GMPE):
    """
    Implements GMPE developed by Norman Abrahamson and Walter Silva and
    published as "Summary of the Abrahamson & Silva NGA Ground-Motion
    Relations" (2008, Earthquakes Spectra, Volume 24, Number 1, pages 67-97).
    This class implements only the equations for mainshock/foreshocks/swarms
    type events, that is the aftershock term (4th term in equation 1, page 74)
    is set to zero. The constant displacement model (page 80) is also not
    implemented (that is equation 1, page 74 is used for all periods and no
    correction is applied for periods greater than the constant displacement
    period). This class implements also the corrections (for standard
    deviation and hanging wall term calculation) as described in:
    http://peer.berkeley.edu/products/abrahamson-silva_nga_report_files/
    AS08_NGA_errata.pdf
    """
    #: Supported tectonic region type is active shallow crust, see paragraph
    #: 'Data Set Selection', see page 68.
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.ACTIVE_SHALLOW_CRUST

    #: Supported intensity measure types are spectral acceleration, peak
    #: ground velocity and peak ground acceleration, see tables 5a and 5b
    #: pages 84, 85, respectively.
    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([
        PGA,
        SA
    ])

    #: Supported intensity measure component is average horizontal
    #: :attr:`~openquake.hazardlib.const.IMC.GMRotI50`
    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.AVERAGE_HORIZONTAL

    #: Supported standard deviation types are inter-event, intra-event
    #: and total, see paragraph "Equations for standard deviations", page 81.
    DEFINED_FOR_STANDARD_DEVIATION_TYPES = set([
        const.StdDev.TOTAL
    ])

    #: Required site parameters are Vs30
    #: Abrahamson & Silva (1997) define only a "rock" and a "soil" class
    #: based on Geomatrix site classification with Geomatrix Classes A and B
    #: classed as 'rock', and C and D as 'soil'. However, distinction
    #: between B and C is not based on Vs30, but on a qualitative description.
    #: For the p
    REQUIRES_SITES_PARAMETERS = set(('vs30'))

    #: Required rupture parameters are magnitude and rake
    #: The GMPE defines three rupture types: Reverse, Reverse/Oblique and
    #: otherwise.
    REQUIRES_RUPTURE_PARAMETERS = set(('mag', 'rake'))

    #: Required distance measures are Rrup and Rx.
    REQUIRES_DISTANCES = set(('rrup', 'rx'))

    def get_mean_and_stddevs(self, sites, rup, dists, imt, stddev_types):
        """
        See :meth:`superclass method
        <.base.GroundShakingIntensityModel.get_mean_and_stddevs>`
        for spec of input and result values.
        """
        # extract dictionaries of coefficients specific to required
        # intensity measure type and for PGA
        C = self.COEFFS[imt]
        C_PGA = self.COEFFS[PGA()]

        # compute median pga on rock (vs30=1100), needed for site response
        # term calculation
        pga_rock = self._get_pga_on_rock(C_PGA, rup, dists)

        mean = (self._get_magnitude_scaling_term(C, rup.mag) +
                self._get_distance_scaling_term(C, rup.mag, dists.rrup) +
                self._get_style_of_faulting_term(C, rup) +
                self._get_hanging_wall_term(C, rup.mag, dists) +
                self._get_site_response_term(C, pga_rock, sites.vs30))
        npts = len(dists.rrup)
        stddevs = self._get_stddevs(C, rup.mag, stddev_types, npts)

        return mean, stddevs

    def _get_pga_on_rock(self, C_PGA, rup, dists):
        """
        """
        return np.exp(
            self._get_magnitude_scaling_term(C_PGA, rup.mag) +
            self._get_distance_scaling_term(C_PGA, rup.mag, dists.rrup) +
            self._get_style_of_faulting_term(C_PGA, rup) +
            self._get_hanging_wall_term(C_PGA, rup.mag, dists))

    def _get_magnitude_scaling_term(self, C, mag):
        """
        Returns the magnitude scaling term described in Equation 4 (page 105)
        """
        if mag <= self.CONSTS['c1']:
            return C['a1'] + self.CONSTS['a2'] * (mag - self.CONSTS['c1']) +\
                C['a12'] * ((8.5 - mag) ** self.CONSTS['n'])
        else:
            return C['a1'] + self.CONSTS['a4'] * (mag - self.CONSTS['c1']) +\
                C['a12'] * ((8.5 - mag) ** self.CONSTS['n'])

    def _get_distance_scaling_term(self, C, mag, rrup):
        """
        Returns the distance scaling term, described in Equation 4 (page 105)
        """
        return (C['a3'] + self.CONSTS['a13'] * (mag - self.CONSTS['c1'])) *\
            np.log(np.sqrt(rrup ** 2. + C['c4'] ** 2.))

    def _get_style_of_faulting_term(self, C, rup):
        """
        Returns the style of faulting term described in Equation 6 (page 106)
        """
        if rup.mag >= self.CONSTS['c1']:
            f3m = C['a6']
        elif rup.mag <= 5.8:
            f3m = C['a5']
        else:
            f3m = C['a5'] + ((C['a6'] - C['a5']) / (self.CONSTS['c1'] - 5.8))
        return f3m * self._get_style_of_faulting_dummy_variable(rup.rake)

    def _get_style_of_faulting_dummy_variable(self, rake):
        """
        """
        if np.fabs(rake - 90.) < 30.0:
            return 1.0
        elif (np.fabs(rake - 45.0) < 15.0) or (np.fabs(rake - 135) < 15.0):
            return 0.5
        else:
            return 0.

    def _get_hanging_wall_term(self, C, mag, dists):
        """
        Returns the hanging wall term
        """
        if dists.rx <= 0.:
            # Is not on the hanging wall side
            return 0.

        # Get hanging wall magnitude term - Equation 8
        if mag <= 5.5:
            fhwm = 0.
        elif mag >= 6.5:
            fhwm = 1.
        else:
            fhwm = mag - 5.5
        return fhwm * self._get_hanging_wall_distance_term(C, dists.rrup)

    def _get_hanging_wall_distance_term(self, C, rrup):
        """
        Returns the hanging wall distance term - Equation 9
        """
        fhwr = np.zeros_like(rrup, dtype=float)
        idx = np.logical_and(rrup > 4.0, rrup <= 8.0)
        fhwr[idx] = C['a9'] * ((rrup[idx] - 4.0) / 4.0)
        idx = np.logical_and(rrup > 8.0, rrup <= 18.0)
        fhwr[idx] = C['a9']
        idx = np.logical_and(rrup > 18.0, rrup <= 24.0)
        fhwr[idx] = C['a9'] * (1.0 - ((rrup - 18.0) / 7.0))
        return fhwr

    def _get_site_response_term(self, C, pga_rock, vs30):
        """
        Returns the site response term described in Equation 10 (page 106)
        """
        if vs30 < 400:
            return C['a10'] + C['a11'] * np.log(pga_rock + self.CONSTS['c5'])
        else:
            return 0.

    def _get_stddevs(self, C, mag, stddev_types, npts):
        """
        Returns only the total sigma term as given in Equation 13
        """
        stddevs = []
        for stddev_type in stddev_types:
            assert stddev_type in self.DEFINED_FOR_STANDARD_DEVIATION_TYPES
            if stddev_type == const.StdDev.TOTAL:
                vals = np.ones(npts, dtype=float)
                if mag <= 5.0:
                    sig_total = C['b5'] * vals
                elif mag >= 7.0:
                    sig_total = (C['b5'] - 2.0 * C['b6']) * vals
                else:
                    sig_total = (C['b5'] - C['b6'] * (mag - 5.)) * vals
                stddevs.append(sig_total)
        return stddevs

    #: Coefficient tables from Table 3 and Table 4 (pages 108 and 109
    #: repectively)
    COEFFS = CoeffsTable(sa_damping=5, table="""\
    IMT     c4       a1      a3       a5       a6      a9      a10      a11       a12      b5      b6
    pga     5.600    1.640   -1.145   0.610    0.260   0.370   -0.417   -0.230    0.0000   0.700   0.135
    0.010   5.600    1.640   -1.145   0.610    0.260   0.370   -0.417   -0.230    0.0000   0.700   0.135
    0.020   5.600    1.640   -1.145   0.610    0.260   0.370   -0.417   -0.230    0.0000   0.700   0.135
    0.030   5.600    1.690   -1.145   0.610    0.260   0.370   -0.470   -0.230    0.0143   0.700   0.135
    0.040   5.600    1.780   -1.145   0.610    0.260   0.370   -0.555   -0.251    0.0245   0.705   0.135
    0.050   5.600    1.870   -1.145   0.610    0.260   0.370   -0.620   -0.267    0.0280   0.713   0.135
    0.060   5.600    1.940   -1.145   0.610    0.260   0.370   -0.665   -0.280    0.0300   0.720   0.135
    0.075   5.580    2.037   -1.145   0.610    0.260   0.370   -0.628   -0.280    0.0300   0.728   0.135
    0.090   5.540    2.100   -1.145   0.610    0.260   0.370   -0.609   -0.280    0.0300   0.735   0.135
    0.100   5.500    2.160   -1.145   0.610    0.260   0.370   -0.598   -0.280    0.0280   0.739   0.135
    0.120   5.390    2.272   -1.145   0.610    0.260   0.370   -0.591   -0.280    0.0180   0.746   0.135
    0.150   5.270    2.407   -1.145   0.610    0.260   0.370   -0.577   -0.280    0.0050   0.754   0.135
    0.170   5.200    2.430   -1.135   0.610    0.260   0.370   -0.522   -0.265   -0.0040   0.759   0.135
    0.200   5.100    2.406   -1.115   0.610    0.260   0.370   -0.445   -0.245   -0.0138   0.765   0.135
    0.240   4.970    2.293   -1.079   0.610    0.232   0.370   -0.350   -0.223   -0.0238   0.772   0.135
    0.300   4.800    2.114   -1.035   0.610    0.198   0.370   -0.219   -0.195   -0.0360   0.780   0.135
    0.360   4.620    1.955   -1.005   0.610    0.170   0.370   -0.123   -0.173   -0.0460   0.787   0.135
    0.400   4.520    1.860   -0.988   0.610    0.154   0.370   -0.065   -0.160   -0.0518   0.791   0.135
    0.460   4.380    1.717   -0.965   0.592    0.132   0.370    0.020   -0.136   -0.0594   0.796   0.132
    0.500   4.300    1.615   -0.952   0.581    0.119   0.370    0.085   -0.121   -0.0635   0.799   0.130
    0.600   4.120    1.428   -0.922   0.557    0.091   0.370    0.194   -0.089   -0.0740   0.806   0.127
    0.750   3.900    1.160   -0.885   0.528    0.057   0.331    0.320   -0.050   -0.0862   0.814   0.123
    0.850   3.810    1.020   -0.865   0.512    0.038   0.309    0.370   -0.028   -0.0927   0.819   0.121
    1.000   3.700    0.828   -0.838   0.490    0.013   0.281    0.423    0.000   -0.1020   0.825   0.118
    1.500   3.550    0.260   -0.772   0.438   -0.049   0.210    0.600    0.040   -0.1200   0.840   0.110
    2.000   3.500   -0.150   -0.725   0.400   -0.094   0.160    0.610    0.040   -0.1400   0.851   0.105
    3.000   3.500   -0.690   -0.725   0.400   -0.156   0.089    0.630    0.040   -0.1726   0.866   0.097
    4.000   3.500   -1.130   -0.725   0.400   -0.200   0.039    0.640    0.040   -0.1956   0.877   0.092
    5.000   3.500   -1.460   -0.725   0.400   -0.200   0.000    0.664    0.040   -0.2150   0.885   0.087
    """)

    CONSTS = {
        # coefficients in table 4, page 84
        'c1': 6.75,
        'a2': 0.512,
        'a4': -0.144,
        'a13': 0.17,
        'n': 2.0,
        'c5': 0.03
    }
