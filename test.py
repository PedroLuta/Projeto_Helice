# import os
# import sys

# dir_path = os.path.dirname(os.path.abspath(__file__))
# root_dir_path = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..'))
# sys.path.append(root_dir_path)
from importing import *

# import PARSEC_functions
# import airfoil_handling
# # import SolidWorks_macro
# import optimization_NSGA2
# import Beta_dist
import prop_simulate
# import xfoil_interface

# x = PARSEC_functions.PFoil(random = True)
# xx, yy, zz = xfoil_interface.get_curve_com_default(100000, 0, 10, 1)
# x.plot()
# plt.plot(xx, yy)
# plt.plot(xx, zz)
# plt.show()

# airfoil = "airfoils\\s1223.txt"

# foil = PARSEC_functions.PFoil(selig_file = airfoil)
# x, y = PARSEC_functions.read_foil(airfoil, PARSEC_functions.count_header_lines(airfoil))
# x_lo, y_lo, x_u, y_u = PARSEC_functions.split_upper_lower(x, y)

# foil.plot_diff(x_lo, y_lo, x_u, y_u)

Betas = [60, 50, 40, 30, 20]
chord_dist = [0.25, 0.3, 0.275, 0.175, 0.05]
r_vector = [0.08, 0.1, 0.15, 0.2, 0.25]
airfoil = "airfoils\\sunnysky.dat"

Vax = 10
rpm = 3000
radps = rpm*2*pi/60
Blades = 6
R = 0.3

dT_vector_q, dQ_vector_q, _, _, _, _, _ = prop_simulate.qprop_fixed_pitch(Vax, radps, Blades, R, r_vector, Betas, chord_dist, airfoil = airfoil)
dT_vector_m, dQ_vector_m = prop_simulate.momentum_Ftip_fixed_pitch(Vax, radps, Blades, R, r_vector, Betas, chord_dist, airfoil = airfoil)

plt.plot(r_vector, dT_vector_q, label = 'qprop')
plt.plot(r_vector, dT_vector_m, label = 'momentum')
plt.legend()
plt.show()