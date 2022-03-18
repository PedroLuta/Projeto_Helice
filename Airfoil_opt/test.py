from importing import *
import simulate
from pyqprop import *
import Chord_dist

rpm = 2500
radps = rpm*2*pi/60
Radius = 0.3
p = 0.2
Blades = 2
V = 16
alphas = [2, 8, 0.25]
sections = 21

r_vector1, chord1, Beta1 = Chord_dist.blade_design_vortex_v1_2(V, radps, Blades, p, Radius, sections = sections, Cl_ref = 1.2, Cd_ref = 0.014, a_ref = 4, Prescribed_power = 520, Prescribed_thrust = 5, init_disp = 0)
r_vector2, chord2, Beta2 = Chord_dist.blade_design_vortex_v1_1(V, radps, Blades, p, Radius, airfoil = "airfoils\\sunnysky.dat", sections = sections, alphas = alphas, Prescribed_power = 520, Prescribed_thrust = 5, init_disp = 1)
#fixed_pitch_qprop_v2(16, radps, 2, [0.1, 0.2, 0.3], [45, 30, 15], [0.06, 0.04, 0], airfoil = 'airfoils\\sunnysky.dat')

R_vec1 = [r/Radius for r in r_vector1]
R_vec2 = [r/Radius for r in r_vector2]

plt.plot(R_vec1, chord1, 'r', label = "fixed Cl and Cd")
plt.plot(R_vec2, chord2, 'b', label = "iterative Cl and Cd")
plt.legend()
plt.show()

# absolute_limits = nup.array([[0, 5], [0, 3]])

# sampling = LHS(xlimits = absolute_limits)

# num = 50
# x = sampling(num)

# plt.plot(x[:, 0], x[:, 1], "o")
# plt.xlabel("x")
# plt.ylabel("y")
# plt.show()

