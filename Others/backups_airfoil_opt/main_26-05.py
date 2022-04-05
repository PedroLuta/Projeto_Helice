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
from PARSEC_airfoilgen import *



#operating conditions
rpm = 2300
rho = 1.0856 #[kg/(m^3)]
kvisc = 1.8/100000 #[kg/(m*s)]
vi = 0.01 #[m/s]
Pmax = 700 #[N*m/s]
radps = (rpm*2*pi)/60 #[1/s]
Available_Torque = Pmax/radps






#propeller characterization
Blades = 4
sections = 21 #more = more precision
p = 0.15 #radius at which there's no more geometric constraints regarding the hub
R = 0.30225






#initialization
Thrust = 0
Torque = 0
Beta_vector = []
radius_vector = []
Thrust_vector = []
Torque_vector = []
first = True




def chord_distribution(x):
    if x > 0.4:
        return 0.8*(-1.8598247908*(x**4) + 4.7243249915*(x**3) - 4.4099344990*(x**2) + 1.5997758465*x + 0.0076143297)
    else:
        return 0.8*(-775.7226922959*(x**6) + 1170.2815273553*(x**5) - 660.9381028488*(x**4) + 161.9663843811*(x**3) - 13.9599756954*(x**2) + 0.1950071591*x + 0.1005241407)




#main starts here
induction_axial = []
induction_radial = []
a1 = 0
a2 = 10
astep = 1



while True:
    #Parsec_Vector = [ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo]
    Pv = [-6.4434, 9.4527, 2.3314E-004, 0.003, 8.1099E-003, 0.4308, 6.2920E-002, -0.4254, 8.4989E-003, 0.3440, -5.8901E-002, 0.7057]
    run_check = airfoilgen(Pv[0], Pv[1], Pv[2], Pv[3], Pv[4], Pv[5], Pv[6], Pv[7], Pv[8], Pv[9], Pv[10], Pv[11])

    if not run_check:
        print("An airfoil was not generated")
        break

    for rr in nup.linspace(0, R, num = sections):
        dT = 0
        dQ = 0

        if rr <= (p*R):
            Beta_vector.append(0)
            radius_vector.append(rr)
            Torque_vector.append(0)
            Thrust_vector.append(0)

        else:
            chord = R*chord_distribution(rr/R)

            Vr = radps*rr
            V = ((Vr**2)+(vi**2))**0.5
            Re = ((V*chord)/kvisc) 
            #print(Re) 
            #print(int(round(Re)))

            Cl, Cd, alpha = get_properties(Re, a1, a2, astep)

            dT, dQ, phi, aai0, aai = jitted_induction(radps, rr, Cl, Cd, Blades, rho, R, chord, vi)

            induction_axial.append(aai*vi)
            induction_radial.append(radps*rr*aai0)
            Beta_vector.append(alpha + math.degrees(phi))
            radius_vector.append(rr)
            Torque_vector.append(dQ)
            Thrust_vector.append(dT)

            if first:
                dQ_remember = dQ
                dT_remember = dT
                r_remember = rr
                first = False
            else:
                Thrust += ((dT + dT_remember)/2)*(rr - r_remember)
                Torque += ((dQ + dQ_remember)/2)*(rr - r_remember)
                dQ_remember = dQ
                dT_remember = dT
                r_remember = rr

    if Torque == 0:
        print('hell no')
    else:
        #print(induction_axial)
        #print(induction_radial)
        #print(Available_Torque)
        print(Torque)
        print(Thrust)  
        #print(Beta_vector)
        #print(chord_vector)
        #print(radius_vector)
        #print(Thrust_vector)
        #print(Torque_vector)
    
    print("loop")
    break