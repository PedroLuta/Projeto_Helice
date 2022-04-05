from numpy.lib.function_base import _i0_1
from importing import *
import xfoil_interface
import induction

def chord_distribution(x, factor = 1):
    if x > 0.4:
        return factor*(-1.8598247908*(x**4) + 4.7243249915*(x**3) - 4.4099344990*(x**2) + 1.5997758465*x + 0.0076143297)
    else:
        return factor*(-775.7226922959*(x**6) + 1170.2815273553*(x**5) - 660.9381028488*(x**4) + 161.9663843811*(x**3) - 13.9599756954*(x**2) + 0.1950071591*x + 0.1005241407)

def vortex(vi, radps, Blades, p, R, airfoil = 'airfoils\\airfoil.txt', sections = 21, rho = 1.225, dvisc = 1.8/100000, alphas = [2, 8, 1], speed_sound = 340):
    a1, a2, astep = alphas[0], alphas[1], alphas[2]
    kvisc = dvisc/rho
    r_vector = nup.linspace(p*R, R, num = sections)
    Beta_vector = []
    dT_vector = []
    dQ_vector = []
    for rr in r_vector:
        chord = R*chord_distribution(rr/R)
        Vr = radps*rr
        V = ((Vr**2)+(vi**2))**0.5
        M = V/speed_sound
        Re = ((V*chord)/kvisc) 
        alpha_c, Cl_c, Cd_c = xfoil_interface.get_curve_com_default(Re, a1, a2, astep, afile = airfoil, M = M)
        alpha, Cl, Cd = xfoil_interface.calculate_most_eff_alpha(alpha_c, Cl_c, Cd_c)
        WA, WT = induction.induction_qprop_adapted(radps, rr, Blades, Cl, R, chord, vi)
        W = (WA**2 + WT**2)**0.5
        phi = math.atan(WA/WT)
        dT, dQ = induction.dT_dQ_from_Cl_Cd(Cl, Cd, phi, W, rho, Blades, chord, rr)
        Beta_vector.append(alpha + math.degrees(phi))
        dQ_vector.append(dQ)
        dT_vector.append(dT)

    return dT_vector, dQ_vector, Beta_vector, r_vector

def fixed_pitch_qprop(vi, radps, Blades, p, R, Beta_dist, airfoil = 'airfoils\\airfoil.txt', sections = 21, rho = 1.225, dvisc = 1.8/100000):
    kvisc = dvisc/rho
    r_vector = nup.linspace(p*R, R, num = sections)
    dT_vector = []
    dQ_vector = []
    Re_vector = []
    WA_vector = []
    WT_vector = []
    Cl_vector = []
    Cd_vector = []
    if len(Beta_dist) != len(r_vector):
        print("Beta vector and radius vector mismatch, please adjust the sections number")
    for i in range(len(r_vector)):
        rr = r_vector[i]
        Beta = Beta_dist[i]
        chord = 2*R*chord_distribution(rr/R)
        Vr = radps*rr
        V = ((Vr**2)+(vi**2))**0.5
        Re = ((V*chord)/kvisc) 
        alpha_c, Cl_c, Cd_c = xfoil_interface.get_curve_com_default(Re, -5, 10, 1, afile = airfoil)
        WA, WT, Cl, Cd = induction.induction_qprop_list(radps, rr, Blades, alpha_c, Cl_c, Cd_c, Beta, R, chord, vi)
        W = (WA**2 + WT**2)**0.5
        phi = math.atan(WA/WT)
        dT, dQ = induction.dT_dQ_from_Cl_Cd(Cl, Cd, phi, W, rho, Blades, chord, rr)
        dQ_vector.append(dQ)
        dT_vector.append(dT)
        Re_vector.append(Re)
        WA_vector.append(WA)
        WT_vector.append(WT)
        Cl_vector.append(Cl)
        Cd_vector.append(Cd)
    return dT_vector, dQ_vector, r_vector, Re_vector, WA_vector, Cl_vector, Cd_vector




