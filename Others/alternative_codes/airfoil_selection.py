#default imports
import numpy as nup
#from scipy import linalg
import math
import os
import subprocess as sp
#import time
from numba import njit, jit

pi = nup.pi
euler = nup.e

#functions import
from xfoil_interface import *
from induction import *
from PARSEC_functions import *


def cut_half(list):
    i = 0
    accumulator = 0
    while i < len(list):
        #print(list[i][3])
        if math.isnan(list[i][3]):
            i += 1
            continue
        accumulator += list[i][3]
        i += 1
    mean_eff = accumulator/len(list)

    i = 0
    list2 = []
    while i < len(list):
        if list[i][3] > mean_eff:
            list2.append(list[i])
        i += 1
    return list2



#operating conditions
rpm = 2300
rho = 1.0856 #[kg/(m^3)]
kvisc = 1.8/100000 #[kg/(m*s)]
vi = 14 #[m/s]
Pmax = 700 #[N*m/s]
radps = (rpm*2*pi)/60 #[1/s]
Available_Torque = Pmax/radps






#propeller characterization
Blades = 4
sections = 5 #more = more precision
p = 0.15 #radius (percentage) at which there's no more geometric constraints regarding the hub
R = 0.30225
D = 2*R
J = vi/(radps*D)






#initialization





def chord_distribution(x, factor = 1):
    if x > 0.4:
        return factor*(-1.8598247908*(x**4) + 4.7243249915*(x**3) - 4.4099344990*(x**2) + 1.5997758465*x + 0.0076143297)
    else:
        return factor*(-775.7226922959*(x**6) + 1170.2815273553*(x**5) - 660.9381028488*(x**4) + 161.9663843811*(x**3) - 13.9599756954*(x**2) + 0.1950071591*x + 0.1005241407)




#main starts here
induction_axial = []
induction_radial = []
a1 = 0
a2 = 10
astep = 1


list_at = os.listdir(r"C:\Users\PEDRO\Desktop\IC_DE_HÉLICE\git_control\git_prop\airfoils\query_airfoiltools")

best_airfoil = ''
Torque_remember = 1000
Thrust_remember = 0
airfoil_list1 = []
Thrust_dist_list = []
Torque_dist_list = []
eff_vec = []
count_vec = []
count = 0

for airfoil in list_at:
    Thrust = 0
    Torque = 0
    Beta_vector = []
    radius_vector = []
    Thrust_vector = []
    Torque_vector = []
    first = True
    #foil = PFoil(selig_file = f"C:\\Users\\PEDRO\\Desktop\\IC_DE_HÉLICE\\git_control\\git_prop\\airfoils\\query_airfoiltools\\{airfoil}")
    #check = write_from_coeffs(foil.aup, foil.alo, coord_file = f"C:\\Users\\PEDRO\\Desktop\\IC_DE_HÉLICE\\git_control\\git_prop\\airfoils\\airfoil.txt")
    #if airfoil == "ag03-il.dat":
    #    break

    for rr in nup.linspace(p*R, R, num = sections):
        dT = 0
        dQ = 0

        chord = R*chord_distribution(rr/R)
        Vr = radps*rr
        V = ((Vr**2)+(vi**2))**0.5
        Re = ((V*chord)/kvisc) 
        alpha, Cl, Cd, _ = communicate_range_flexible(Re, a1, a2, astep, afile = f'query_airfoiltools\\{airfoil}')
        try:
            dT, dQ, phi, aai0, aai = jitted_induction_Ftip(radps, rr, Cl, Cd, Blades, rho, R, chord, vi)
        except:
            dT, dQ, phi, aai0, aai = 0, 0, 0, 0, 0 
        #induction_axial.append(aai*vi)
        #induction_radial.append(radps*rr*aai0)
        if first or (dT == 0 and dQ == 0) or math.isnan(dT) or math.isnan(dQ):
            dQ_remember = dQ
            dT_remember = dT
            r_remember = rr
            first = False
        else:
            Beta_vector.append(alpha + math.degrees(phi))
            radius_vector.append(rr)
            Torque_vector.append(dQ)
            Thrust_vector.append(dT)
            Thrust += ((dT + dT_remember)/2)*(rr - r_remember)
            Torque += ((dQ + dQ_remember)/2)*(rr - r_remember)
            dQ_remember = dQ
            dT_remember = dT
            r_remember = rr
    print(f"{airfoil} calculated")
    if Torque == 0:
        pass
    elif True:
        airfoil_list1.append([airfoil, Thrust, Torque, Thrust/Torque])
        Thrust_dist_list.append(Thrust_vector)
        Torque_dist_list.append(Torque_vector)
        kt = Thrust/(rho*(radps**2)*(D**4))
        kq = Torque/(rho*(radps**2)*(D**5))
        eff = (1/(2*pi))*(kt/kq)*J
        count_vec.append(count)
        count += 1
        eff_vec.append(eff)

#plt.plot(eff_vec)
plt.scatter(count_vec, eff_vec, s = 5)
plt.show()

#list = cut_half(airfoil_list1)
#with open("out.txt", 'w') as o:
#    i = 0
#    while i < len(airfoil_list1):
#        o.write(f"{airfoil_list1[i]}\n")
#        o.write(f"{Thrust_dist_list[i]}\n")
#        o.write(f"{Torque_dist_list[i]}\n\n")
#        i += 1
    
#list = cut_half(list)
#list = cut_half(list)
#list = cut_half(list)
#print(list)

