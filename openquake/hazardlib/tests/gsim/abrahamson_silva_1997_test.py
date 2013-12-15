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

from openquake.hazardlib.gsim.abrahamson_silva_1997 import \
    AbrahamsonSilva1997
from openquake.hazardlib.tests.gsim.utils import BaseGSIMTestCase

class AbrahamsonSilva1997TestCase(BaseGSIMTestCase):
    GSIM_CLASS = AbrahamsonSilva1997
    
    # Test data constructed using the Matlab implementation by Jack Baker
    # available at http://www.stanford.edu/~bakerjw/GMPEs.html
    def test_mean(self):
        self.check('AS1997/AS97_MEAN.csv',
                   max_discrep_percentage=0.5)

    def test_std_total(self):
        self.check('AS1997/AS97_TOTAL.csv',
                   max_discrep_percentage=0.5)
