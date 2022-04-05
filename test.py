import os
import sys

# dir_path = os.path.dirname(os.path.abspath(__file__))
# root_dir_path = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..'))
# sys.path.append(root_dir_path)
from importing import *

import PARSEC_functions
import airfoil_handling
# import SolidWorks_macro
import optimization_NSGA2
import Beta_dist
import prop_simulate
import xfoil_interface

x = PARSEC_functions.PFoil(random = True)
xx, yy, zz = xfoil_interface.get_curve_com_default(100000, 0, 10, 1)
x.plot()
plt.plot(xx, yy)
plt.plot(xx, zz)
plt.show()