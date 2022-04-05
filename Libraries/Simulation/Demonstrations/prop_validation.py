from importing import *
import Beta_dist

vi = 150
Pmax = 560
rho = 1.225
dvisc = 1.8/100000
kvisc = dvisc/rho

Blades = 2
p = 0.2
R = 3.048/2
D = R*2
airfoil = 'airfoils\\naca16-009.txt'#NACA 10-(3)(062)-045 propeller
#airfoil = 'airfoils\\SD8000-089-88.txt'
sections = 9

rpm = 2300
rps = rpm/60
radps = 2*pi*rpm/60

#Beta_vector = [61.30475642056321, 53.698149009210034, 59.90488514583296, 55.17941843793663, 51.57797055905809, 48.29090872404107, 45.879909897651444, 42.832671953954474, 40.63210106205273, 38.52547131333741, 25.32804121807452, 34.1138776411056, 32.597182165668045, 31.218485002500472, 29.960345510442345, 28.81298051644318, 27.767412340577742, 26.814608083810256, 25.946695026156938, 25.15810588121797, 24.442261023764274, 23.794891683150713, 23.21219439020902, 22.69393887616421, 22.242126679640656, 21.871075957844546, 21.611828352903267, 22.258068493412765, 22.96372655158519, 3.172526986476278]
#Beta_vector = [77, 72.75, 68.5, 65, 61.5, 58.75, 56, 53.5, 51, 49, 47, 45, 43, 42, 41, 39, 37] #NACA 10-(3)(062)-045 propeller
Beta_vector = [77, 68.5, 61.5, 56, 51, 47, 43, 41, 37] #NACA 10-(3)(062)-045 propeller
#Beta_vector = [68, 65, 61, 56, 51, 48, 45, 43, 42, 41] 

#Betas = [15.873628, 29.830716, 40.20649549, 40.20649549, 35.16022374, 31.12177578, \
#    27.84693006, 25.15422488, 22.9103745, 21.01722617, 19.4018361, 18.00934099, 16.79793758, \
#        15.73535684, 14.79637473, 13.96103831, 13.21338756, 12.54052339, 11.93191978]

#dT_vector1, dQ_vector1, Beta_vector, r_vector = simulate.momentum_Ftip(vi, radps, Blades, p, R, airfoil = airfoil, sections = sections, rho = rho, dvisc = dvisc, alphas = [0, 10, 0.5])
#dT_vector, dQ_vector, Beta_vector, r_vector = simulate.vortex(vi, radps, Blades, p, R, airfoil = airfoil, sections = sections, rho = rho, dvisc = dvisc, alphas = [0, 10, 0.5])

Ct_vec = []
Cp_vec = []
eff_vec = []
J_vec = [1.4, 1.6, 1.8, 2.0, 2.2, 2.4]
for J in J_vec:
    #J = vi/(rps*D)
    rps = vi/(J*D)
    radps = rps*2*pi
    #vi = J*rps*D
    dT_vector, dQ_vector, r_vector, Re_dist, WA_dist, Cl_dist, Cd_dist = Beta_dist.fixed_pitch_qprop(vi, radps, Blades, p, R, Beta_vector, airfoil = airfoil, sections = sections, rho = rho, dvisc = dvisc)

# Thrust1 = Auxiliary.area_under_curve(r_vector, dT_vector1)
# Torque1 = Auxiliary.area_under_curve(r_vector, dQ_vector1)
    Thrust = Auxiliary.area_under_curve(r_vector, dT_vector)
    Torque = Auxiliary.area_under_curve(r_vector, dQ_vector)
    Power = Torque*radps
    Cp = Power/(rho*(rps**3)*(D**5))
    Ct = Thrust/(rho*(rps**2)*(D**4))
    eff = (Ct/Cp)*J
    Ct_vec.append(Ct)
    Cp_vec.append(Cp)
    eff_vec.append(eff)

oficial_Ct = [0.145, 0.133, 0.11, 0.08, 0.067, 0.035]
oficial_Cp = [0.31, 0.28, 0.24, 0.18, 0.135, 0.075]
oficial_eff = [0.665, 0.765, 0.845, 0.885, 0.92, 0.875]
# print(J_vec)
# print(Ct_vec)
# print(oficial_Ct)
# plt.plot(J_vec, Ct_vec, label = "analysis")
# plt.plot(J_vec, oficial_Ct, label = "experimental")
# plt.xlabel("Advance Ratio")
# plt.ylabel("Ct")
# plt.legend()
# plt.show()

plt.plot(J_vec, Cp_vec, label = "analysis")
plt.plot(J_vec, oficial_Cp, label = "experimental")
plt.xlabel("Advance Ratio")
plt.ylabel("Cp")
plt.legend()
plt.show()

plt.plot(J_vec, eff_vec, label = "analysis")
plt.plot(J_vec, oficial_eff, label = "experimental")
plt.xlabel("Advance Ratio")
plt.ylabel("eff")
plt.legend()
plt.show()


#print(r_vector)
# print([R*simulate.chord_distribution(x/R) for x in r_vector])
# print(Beta_vector)
# print(Thrust1)
# print(Torque1)
#print(Thrust)
#print(Torque)
# print(Re_dist)
# print(WA_dist)
# print(Cl_dist)
# print(Cd_dist)
# print(dT_vector1)
# print(dQ_vector1)


#print(dT_vector)
#print(dQ_vector)
#print(J)
