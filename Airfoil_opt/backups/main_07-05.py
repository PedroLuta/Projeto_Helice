#default imports
import numpy as nup
#from scipy import linalg
import math
import os
import subprocess as sp
#import time
from numba import njit, jit

#functions import
from xfoil_interface import *
from induction import *

xfoilpath = r'C:\Users\PEDRO\Desktop\python xfoil scripts\xfoil.exe'
xfoil_outfile = r'C:\Users\PEDRO\Desktop\python xfoil scripts\xfoil_output.txt'
#coeff_file = r'C:\Users\PEDRO\Desktop\python xfoil scripts\Coeff.txt' 
#coord_file = r'C:\Users\PEDRO\Desktop\python xfoil scripts\Coor_XY.txt' 

pi = nup.pi
euler = nup.e






#operating conditions
rpm = 2300
rho = 1.0856 #[kg/(m^3)]
kvisc = 1.8/100000 #[kg/(m*s)]
vi = 15 #[m/s]
Pmax = 700 #[N*m/s]
radps = (rpm*2*pi)/60 #[1/s]
Available_Torque = Pmax/radps







#xfoil processing constants
itr = 150
a1 = 2
a2 = 10
astep = 1






#propeller characterization
Blades = 4
sections = 20 #more = more precision
p = 0.15 #radius at which there's no more geometric constraints regarding the hub
R = 0.30225






#initialization
Thrust = 0
Torque = 0
chord_vector = []
Beta_vector = []
radius_vector = []
Thrust_vector = []
Torque_vector = []
dQ_remember = 0
dT_remember = 0
first = True






#subprocess stuff for quality of life
stout = 0
startupinfo = sp.STARTUPINFO()
startupinfo.dwFlags |= sp.STARTF_USESHOWWINDOW





#functions import
#from xfoil_interface import *


def chord_distribution(x):
    if x > 0.4:
        return 0.8*(-1.8598247908*(x**4) + 4.7243249915*(x**3) - 4.4099344990*(x**2) + 1.5997758465*x + 0.0076143297)
    else:
        return 0.8*(-775.7226922959*(x**6) + 1170.2815273553*(x**5) - 660.9381028488*(x**4) + 161.9663843811*(x**3) - 13.9599756954*(x**2) + 0.1950071591*x + 0.1005241407)


















#main starts here

induction_axial = []
induction_radial = []

for rr in nup.linspace(0, R, num = sections):
    dT = 0
    dQ = 0

    if rr <= (p*R):
        Beta_vector.append(0)
        chord_vector.append(0)
        radius_vector.append(rr)

    else:
        chord = R*chord_distribution(rr/R)

        Vr = radps*rr
        V = ((Vr**2)+(vi**2))**0.5
        Re = ((V*chord)/kvisc)  #*(10**-4)

        Cl, Cd, alpha = get_properties(Re)

        #exp_func1 = (-Blades/2)*((R-rr)/rr)
        #dT1 = (rho*Blades*chord)/2
        #dQ1 = (rho*Blades*chord*rr)/2           
        #denominator1 = 4*pi*rr*rho
        #denominator2 = 4*pi*(rr**3)*rho

        #dT, dQ, phi = correct_for_induction(radps, rr, exp_func1, dT1, dQ1, Cl, Cd, denominator1, denominator2) 
        dT, dQ, phi, aai0, aai = jitted_induction(radps, rr, Cl, Cd, Blades, rho, R, chord, vi)
        
        induction_axial.append(aai*vi)
        induction_radial.append(radps*rr*aai0)

        chord_vector.append(chord)
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