from importing import *
import airfoil_handling
import Chord_dist
import Beta_dist
import Auxiliary
import simulate
import xfoil_interface
import Cube_compatibility

truncate = Auxiliary.truncate_float

Diameter_inches = 12
Radius_meters = Diameter_inches*0.0254/2
R = Radius_meters 
#print(f'maximum chord recommended: {Radius_meters*0.2}')
V = 12
Blades = 5
rpm = 12000
radps = rpm*2*pi/60


airfoil = 'airfoils\\NACA6407_modBF.txt'
anchor_position = 0.3 #Naca 4-digit standard max thickness point = 0.3
extra_thickness_at_root = 3
offset_at_root = 0.5

rho = 1.09
kvisc = 1.69*(10**-5)

min_pitch = 0
max_pitch = (R*2*39.37)
pitches_to_test = 201

p1 = 0.1
p3 = 0.35
#sections1 = 8
sections2 = 21
#r_vector1 = nup.linspace(p1*R, p3*R, num = sections1, endpoint = False)
#r_vector2 = nup.linspace(p3*R, R, num = sections2)
r_vector2 = nup.linspace(p1*R, R, num = sections2)

#chord_vector, Beta_vector_from_chord_script = Chord_dist.blade_design_vortex_v1_1_1(V, rpm*2*pi/60, Blades, R, r_vector2, airfoil = airfoil, Prescribed_power = 580, init_disp = 1, rho = rho, alphas = [2, 8, 1])
#chord_vector = [0.04125642828317191, 0.04102379549791816, 0.038872677163378666, 0.0363426138260957, 0.03391140449781531, 0.03155694056252, 0.02936883477375551, 0.027386481892879972, 0.02557174580752842, 0.023920888026726247, 0.022409328483390587, 0.021020572401307187, 0.01973756813161086, 0.01851953186403174, 0.021909436767059424, 0.020233042645877854, 0.0184778640564797, 0.016446542776239966, 0.011186922289564678, 0.002, 0.002]
#chord_vector = [0.03602142101403248, 0.03379152332925716, 0.03169101962201139, 0.029747626059909434, 0.027953475594706734, 0.0263001102785749, 0.02477865816439229, 0.0233676799847905, 0.022067144265344698, 0.02085684174244585, 0.019684232037484525, 0.023466722295084754, 0.021911546136561273, 0.02029523081921109, 0.018529898043166324, 0.016466483150171903, 0.011170244691841679, 0.002, 0.002]
# coefs_chord = nup.polyfit(r_vector2, chord_vector, 5)
# chord_poly = [coefs_chord[0]*x**5 + coefs_chord[1]*x**4 + coefs_chord[2]*x**3 + coefs_chord[3]*x**2 + coefs_chord[4]*x + coefs_chord[5] for x in r_vector2]
# chord_poly[-1] = 0.005

Re_vec = []
Bladesv = [6, 1]
rpmv = [2000, 12000]
for i in range(2):
    rpm, Blades = rpmv[i], Bladesv[i]
    x = 0.75*R
    chord = 0.06
    Vrrr = rpmv[i]*2*pi*x/60
    ww = ((Vrrr**2) + (V**2))**0.5
    Reynolds = chord*ww/kvisc
    Re_vec.append(Reynolds)
print(Re_vec)

# plt.plot(r_vector2, chord_poly)
# plt.ylim(0, 0.1)
# plt.show()

#Beta_vector_momentum = Beta_dist.momentum_Ftip_highest_efficiency(V, rpm*2*pi/60, Blades, R, r_vector, chord_vector, airfoil = airfoil, alphas = [2, 8, 1])

#Beta_vector_qprop = Beta_dist.qprop_highest_efficiency(V, rpm*2*pi/60, Blades, R, r_vector2, chord_poly, airfoil = airfoil, alphas = [2, 8, 1])
# Beta_vector_simple_pitch = Beta_dist.simple_pitch_highest_efficiency(V, rpm*2*pi/60, r_vector2, Blades, R, chord_poly, min_pitch, max_pitch, pitches_to_test, airfoil = airfoil)

# Beta_vector_simple_pitch[-1] = 0
# Beta_vector_simple_pitch[-2] = 0

# chord3 = chord_poly[0]
# chord1 = (2*pi*p1*R/Blades)/1.3
# der_chord3 = 0#((chord_poly[1] - chord_poly[0])/(r_vector2[1] - r_vector2[0]))*0.5
# der_chord1 = 0
# c_chord = Cube_compatibility.compatibility_function(p1, p3, chord1, chord3, der_chord1, der_chord3)
# aa, bb, cc, dd = c_chord[0], c_chord[1], c_chord[2], c_chord[3] 
# chord_vector_root = [aa*(x/R)**3 + bb*(x/R)**2 + cc*(x/R) + dd for x in r_vector1]
# chord_vector_root.extend(chord_poly)
# chord_vector = chord_vector_root.copy()

# Pitch3 = Beta_vector_simple_pitch[0]
# Pitch1 = 0
# der_Pitch3 = 0#((Beta_vector_simple_pitch[1] - Beta_vector_simple_pitch[0])/(r_vector2[1] - r_vector2[0]))*0.5
# der_Pitch1 = 0
# c_Pitch = Cube_compatibility.compatibility_function(p1, p3, Pitch1, Pitch3, der_Pitch1, der_Pitch3)
# aa, bb, cc, dd = c_Pitch[0], c_Pitch[1], c_Pitch[2], c_Pitch[3] 
# Beta_vector_root = [aa*(x/R)**3 + bb*(x/R)**2 + cc*(x/R) + dd for x in r_vector1]
# Beta_vector_root.extend(Beta_vector_simple_pitch)
# Beta_vector_simple = Beta_vector_root.copy()

# offset3 = anchor_position
# offset1 = 0.5
# der_off3 = 0
# der_off1 = 0
# c_offset = Cube_compatibility.compatibility_function(p1, p3, offset1, offset3, der_off1, der_off3)
# aa, bb, cc, dd = c_offset[0], c_offset[1], c_offset[2], c_offset[3]
# offset_vector_root = [aa*(x/R)**3 + bb*(x/R)**2 + cc*(x/R) + dd for x in r_vector1]

# thickness3 = 1
# thickness1 = extra_thickness_at_root
# der_tck3 = 0
# der_tck1 = 0
# c_thick = Cube_compatibility.compatibility_function(p1, p3, thickness1, thickness3, der_tck1, der_tck3)
# aa, bb, cc, dd = c_thick[0], c_thick[1], c_thick[2], c_thick[3]
# thickness_vector_root = [aa*(x/R)**3 + bb*(x/R)**2 + cc*(x/R) + dd for x in r_vector1]

# Pitch3 = Beta_vector_qprop[0]
# Pitch1 = 0
# der_Pitch3 = 0#((Beta_vector_qprop[1] - Beta_vector_qprop[0])/(r_vector2[1] - r_vector2[0]))*0.5
# der_Pitch1 = 0
# c_Pitch = Cube_compatibility.compatibility_function(p1, p3, Pitch1, Pitch3, der_Pitch1, der_Pitch3)
# aa, bb, cc, dd = c_Pitch[0], c_Pitch[1], c_Pitch[2], c_Pitch[3] 
# Beta_vector_root = [aa*(x/R)**3 + bb*(x/R)**2 + cc*(x/R) + dd for x in r_vector1]
# Beta_vector_root.extend(Beta_vector_qprop)
# Beta_vector_qprop = Beta_vector_root.copy()

# r_vector = nup.append(r_vector1, r_vector2)

# x_orig, y_orig = airfoil_handling.read_foil(airfoil, header_lines = airfoil_handling.count_header_lines(airfoil))


# Beta_vector = Beta_vector_simple

# for i in range(len(r_vector)):
#     output_file = f"outputs\\simple_{i + 1}.txt"
#     x, y = x_orig.copy(), y_orig.copy()
#     if i < len(r_vector1):
#         thickness = thickness_vector_root[i]
#         offset = offset_vector_root[i]
#         x, y = airfoil_handling.apply_thickness_on_camber_easy(x, y, thickness)
#     else:
#         offset = anchor_position
#     x = airfoil_handling.apply_offset_x(x, -offset)
#     x, y = airfoil_handling.apply_chord(x, y, chord_vector[i])
#     x, y = airfoil_handling.apply_rotation_deg(x, y, Beta_vector[i])
#     with open(output_file, 'w') as out:
#         for ii in range(len(x)):
#             out.write(f"{truncate(r_vector[i], 5)} {truncate(x[ii], 5)} {truncate(y[ii], 5)}\n")

# Beta_vector = Beta_vector_qprop

# for i in range(len(r_vector)):
#     output_file = f"outputs\\qprop_{i + 1}.txt"
#     x, y = x_orig.copy(), y_orig.copy()
#     if i < len(r_vector1):
#         thickness = 5*(1 - tck)*(r_vector[i]/R) + ((3*tck) - 1)/2
#         x, y = airfoil_handling.apply_thickness_on_camber_easy(x, y, thickness)
#     x = airfoil_handling.apply_offset_x(x, -anchor_position)
#     x, y = airfoil_handling.apply_chord(x, y, chord_vector[i])
#     x, y = airfoil_handling.apply_rotation_deg(x, y, Beta_vector[i])
#     with open(output_file, 'w') as out:
#         for ii in range(len(x)):
#             out.write(f"{truncate(r_vector[i], 5)} {truncate(x[ii], 5)} {truncate(y[ii], 5)}\n")
