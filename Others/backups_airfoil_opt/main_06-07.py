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
    area = 0
    i = 1
    while i < len(x):
        area += (x[i] - x[i - 1])*(curve[i] + curve[i - 1])/2
        i += 1
    return area

def calculate_T_Q(vi, rpm, Blades, p, R, airfoil = 'airfoils\\airfoil.txt', sections = 21, rho = 1.225, dvisc = 1.8/100000, alphas = [0, 10, 1]):
    a1, a2, astep = alphas[0], alphas[1], alphas[2]
    kvisc = dvisc/rho
    Beta_vector = []
    r_vector = []
    dT_vector = []
    dQ_vector = []
    radps = (rpm*2*pi)/60
    for rr in nup.linspace(p*R, R, num = sections):
        chord = R*chord_distribution(rr/R)
        Vr = radps*rr
        V = ((Vr**2)+(vi**2))**0.5
        Re = ((V*chord)/kvisc) 
        alpha, Cl, Cd, Cm = communicate_range_flexible(Re, a1, a2, astep, afile = airfoil)
        try:
            dT, dQ, phi, aai0, aai = jitted_induction_Ftip(radps, rr, Cl, Cd, Blades, rho, R, chord, vi)
        except:
            dT, dQ, phi, aai0, aai = 0, 0, 0, 0, 0 
        Beta_vector.append(alpha + math.degrees(phi))
        r_vector.append(rr)
        dQ_vector.append(dQ)
        dT_vector.append(dT)

    Thrust = area_under_curve(r_vector, dT_vector)
    Torque = area_under_curve(r_vector, dQ_vector)

    return Thrust, Torque

##operating conditions
#rpm = 2300
#rho = 1.0856 #[kg/(m^3)]
#dvisc = 1.68/100000 #[kg/(m*s)]
#kvisc = dvisc/rho
#vi = 14 #[m/s]
#radps = (rpm*2*pi)/60 #[1/s]
#
##propeller characterization
#Blades = 4
#sections = 5 #more = more precision
#p = 0.15 #radius (percentage) at which there's no more geometric constraints regarding the hub
#R = 0.30225
#
##main starts here
#a1, a2, astep = 0, 10, 1
#Beta_vector = []
#r_vector = []
#dT_vector = []
#dQ_vector = []
#
#for rr in nup.linspace(p*R, R, num = sections):
#    chord = R*chord_distribution(rr/R)
#    Vr = radps*rr
#    V = ((Vr**2)+(vi**2))**0.5
#    Re = ((V*chord)/kvisc) 
#    alpha, Cl, Cd, _ = communicate_range_flexible(Re, a1, a2, astep, afile = f'airfoils\\clarky.txt')
#    try:
#        dT, dQ, phi, aai0, aai = jitted_induction_Ftip(radps, rr, Cl, Cd, Blades, rho, R, chord, vi)
#    except:
#        dT, dQ, phi, aai0, aai = 0, 0, 0, 0, 0 
#    Beta_vector.append(alpha + math.degrees(phi))
#    r_vector.append(rr)
#    dQ_vector.append(dQ)
#    dT_vector.append(dT)
#
#Thrust = area_under_curve(r_vector, dT_vector)
#Torque = area_under_curve(r_vector, dQ_vector)
Thrust, Torque = calculate_T_Q(14, 2300, 4, 0.15, 0.30225, airfoil = 'airfoils\\clarky.txt', sections = 5, rho = 1.0856, dvisc = 1.68/100000)
print(Thrust)
print(Torque)