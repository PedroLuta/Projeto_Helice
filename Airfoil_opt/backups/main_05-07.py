import numpy as nup
import math

pi = nup.pi
euler = nup.e

#functions import
from xfoil_interface import *
from induction import *
from PARSEC_functions import *

def chord_distribution(x, factor = 1):
    if x > 0.4:
        return factor*(-1.8598247908*(x**4) + 4.7243249915*(x**3) - 4.4099344990*(x**2) + 1.5997758465*x + 0.0076143297)
    else:
        return factor*(-775.7226922959*(x**6) + 1170.2815273553*(x**5) - 660.9381028488*(x**4) + 161.9663843811*(x**3) - 13.9599756954*(x**2) + 0.1950071591*x + 0.1005241407)

def area_under_curve(x, curve):
    i = 0
    area = 0
    while i < len(x):
        if i == 0:
            val_remember = curve[i]
            x_remember = x[i]
            i += 1
        area += (x[i] - x_remember)*(curve[i] + val_remember)/2
        val_remember = curve[i]
        x_remember = x[i]
        i += 1
    return area

#operating conditions
rpm = 2300
rho = 1.0856 #[kg/(m^3)]
dvisc = 1.68/100000 #[kg/(m*s)]
kvisc = dvisc/rho
vi = 14 #[m/s]
radps = (rpm*2*pi)/60 #[1/s]

#propeller characterization
Blades = 4
sections = 5 #more = more precision
p = 0.15 #radius (percentage) at which there's no more geometric constraints regarding the hub
R = 0.30225
D = 2*R
J = vi/(radps*D)

#main starts here
#induction_axial = []
#induction_radial = []
a1 = 0
a2 = 10
astep = 1
#Torque_remember = 1000
#Thrust_remember = 0
Thrust = 0
Torque = 0
Beta_vector = []
radius_vector = []
Thrust_vector = []
Torque_vector = []
first = True

for rr in nup.linspace(p*R, R, num = sections):
    dT = 0
    dQ = 0
    chord = R*chord_distribution(rr/R)
    Vr = radps*rr
    V = ((Vr**2)+(vi**2))**0.5
    Re = ((V*chord)/kvisc) 
    alpha, Cl, Cd, _ = communicate_range_flexible(Re, a1, a2, astep, afile = f'airfoils\\clarky.txt')
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
if Torque == 0:
    pass
else:
    kt = Thrust/(rho*(radps**2)*(D**4))
    kq = Torque/(rho*(radps**2)*(D**5))
    eff = (1/(2*pi))*(kt/kq)*J

