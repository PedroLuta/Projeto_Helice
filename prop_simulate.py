import xfoil_interface
import math
import numpy as np
import pandas as pd
pi = np.pi
euler = np.euler

def qprop_fixed_pitch(vi, radps, Blades, R, r_vector, Beta_dist, chord_dist, airfoil = 'airfoils\\airfoil.txt', rho = 1.225, dvisc = 1.8/100000, alphas = [-5, 10, 1]):
    #CHORD DISTRIBUTION TAKEN AS INPUT
    a1, a2, astep = alphas[0], alphas[1], alphas[2]
    kvisc = dvisc/rho

    dT_vector = []
    dQ_vector = []
    Re_vector = []
    WA_vector = []
    WT_vector = []
    Cl_vector = []
    Cd_vector = []

    for i in range(len(r_vector)):
        rr = r_vector[i]
        Beta = Beta_dist[i]
        chord = chord_dist[i]

        Vr = radps*rr
        V = ((Vr**2)+(vi**2))**0.5
        Re = ((V*chord)/kvisc) 

        alpha_c, Cl_c, Cd_c = xfoil_interface.get_curve_com_default(Re, a1, a2, astep, afile = airfoil)
        WA, WT, Cl, Cd = induction_qprop_fixed_pitch(radps, rr, Blades, alpha_c, Cl_c, Cd_c, Beta, R, chord, vi)
        W = (WA**2 + WT**2)**0.5
        phi = math.atan(WA/WT)
        dT = (rho*Blades*chord)*(W**2)*(Cl*math.cos(phi) - Cd*math.sin(phi))/2
        dQ = (rho*Blades*chord*rr)*(W**2)*(Cl*math.sin(phi) + Cd*math.cos(phi))/2


        dQ_vector.append(dQ)
        dT_vector.append(dT)
        Re_vector.append(Re)
        WA_vector.append(WA)
        WT_vector.append(WT)
        Cl_vector.append(Cl)
        Cd_vector.append(Cd)
    return dT_vector, dQ_vector, r_vector, Re_vector, WA_vector, Cl_vector, Cd_vector

def momentum_Ftip_fixed_pitch(vi, radps, Blades, R, r_vector, Beta_dist, chord_dist, airfoil = 'airfoils\\airfoil.txt', rho = 1.225, dvisc = 1.8/100000, alphas = [-5, 10, 1]):
    a1, a2, astep = alphas[0], alphas[1], alphas[2]
    kvisc = dvisc/rho

    dT_vector = []
    dQ_vector = []

    for i in range(len(r_vector)):
        rr = r_vector[i]
        chord = chord_dist[i]
        Beta = Beta_dist[i]

        Vr = radps*rr
        V = ((Vr**2)+(vi**2))**0.5
        Re = ((V*chord)/kvisc)

        alpha_c, Cl_c, Cd_c = xfoil_interface.get_curve_com_default(Re, a1, a2, astep, afile = airfoil)

        dT, dQ = induction_momentum_Ftip_fixed_pitch(radps, rr, Cl_c, Cd_c, alpha_c, Beta, Blades, rho, R, chord, vi) 
        dQ_vector.append(dQ)
        dT_vector.append(dT)

    return dT_vector, dQ_vector

def BEMT_PrePolars(vi, radps, Blades, R, r_vector, Beta_dist, chord_dist, CLPolar = 'TestCLPolar.dat', CDPolar = 'TestCDPolar.dat', rho = 1.225, dvisc = 1.8/100000):
    kvisc = dvisc/rho

    dT_vector = []
    dQ_vector = []

    for i in range(len(r_vector)):
        rr = r_vector[i]
        chord = chord_dist[i]
        Beta = Beta_dist[i]

        Vr = radps*rr
        V = ((Vr**2)+(vi**2))**0.5
        Re = ((V*chord)/kvisc)

        # alpha_c, Cl_c, Cd_c = xfoil_interface.get_curve_com_default(Re, a1, a2, astep, afile = airfoil)

        dT, dQ = induction_momentum_Ftip_fixed_pitch(radps, rr, Cl_c, Cd_c, alpha_c, Beta, Blades, rho, R, chord, vi) 
        dQ_vector.append(dQ)
        dT_vector.append(dT)

    return dT_vector, dQ_vector










#Induction methods
def induction_qprop_fixed_pitch(OMG, rr, BLDS, a_list, CL_list, CD_list, Beta, RAD, CHORD, VEL):
    EPS = 1E-06
    UA     = VEL   
    UT     = OMG*rr 

    WZ = (UA**2 + UT**2)**0.5

    PSImid = 0
    PSIup = math.radians(90)
    PSIlo = math.radians(-90)

    first = True

    while True:
        if first:
            RESup, _, _, _, _ = calculate_residual_fixed_pitch(UA, UT, WZ, Beta, PSIup, CHORD, a_list, CL_list, CD_list, rr, BLDS, RAD)
            RESlo, _, _, _, _ = calculate_residual_fixed_pitch(UA, UT, WZ, Beta, PSIlo, CHORD, a_list, CL_list, CD_list, rr, BLDS, RAD)
            first = False
        else:
            PSImid = (PSIup + PSIlo)/2

        RESmid, WAmid, WTmid, CL, CD = calculate_residual_fixed_pitch(UA, UT, WZ, Beta, PSImid, CHORD, a_list, CL_list, CD_list, rr, BLDS, RAD)

        if(abs(PSIup - PSIlo) < EPS):
            return WAmid, WTmid, CL, CD

        if RESup*RESmid < 0:
            RESlo = RESmid
            PSIlo = PSImid
        elif RESlo*RESmid < 0:
            RESup = RESmid
            PSIup = PSImid
        else:
            #print(f"Induction failed, section at radial position {(rr/RAD)*100}% will be assumed as simple flow")
            CL, CD = find_alpha_interval_return_CL_CD(a_list, CL_list, CD_list, math.degrees(math.atan(UA/UT)) - Beta)
            return UA, UT, CL, CD

def induction_momentum_Ftip_fixed_pitch(radps, rr, Cl_c, Cd_c, alpha_c, Beta, Blades, rho, R, chord, Vax_before):
    check = 0
    ai = 0.1
    ai0 = 0.01

    while True:
        Vr = radps*rr*(1-ai0)
        Vax = Vax_before*(1+ai)
        V = ((Vr**2)+(Vax**2))**0.5

        phi = math.atan(Vax/Vr)
        alpha = Beta - math.degrees(phi)
        Cl, Cd = find_alpha_interval_return_CL_CD(alpha_c, Cl_c, Cd_c, alpha)

        exp_func = (-Blades/2)*((R-rr)/rr)*(V/Vax)
        F_tip = (2/pi)*math.acos(euler**exp_func)

        dT = (rho*(V**2)/2)*(Cl*math.cos(phi)-Cd*math.sin(phi))*(Blades*chord)
        dQ = (rho*(V**2)/2)*(Cl*math.sin(phi)+Cd*math.cos(phi))*(Blades*chord)*rr

        ai_new = dT/(4*pi*rr*rho*Vax*Vax_before*F_tip) #
        ai0_new = dQ/(4*pi*rr*rho*Vax*(radps*rr)*F_tip*rr) #

        ai_middle = (ai_new + ai)/2
        ai0_middle = (ai0_new + ai0)/2

        if ((abs(ai_middle - ai) < 1/100000) and (abs(ai0_middle - ai0) < 1/100000)) or check > 500:
            return dT, dQ

        ai = ai_middle
        ai0 = ai0_middle
        check += 1





#Residuals
def calculate_residual_fixed_pitch(UA, UT, WZ, Beta, PSI, CHORD, a_list, CL_list, CD_list, rr, BLDS, RAD):
    COSP = math.cos(PSI)
    SINP = math.sin(PSI)
    WA     = 0.5*UA     + 0.5*WZ    *SINP
    WT     = 0.5*UT     + 0.5*WZ    *COSP
    PHI = math.degrees(math.atan(WA/WT))
    alpha = Beta - PHI
    CL, CD = find_alpha_interval_return_CL_CD(a_list, CL_list, CD_list, alpha)
    if (WA <= 0.0):
        F     = 1.0
        ADW     = 0
    else:
        TSR = WT/WA * RAD/rr
        FARG     = 0.5*BLDS*(1.0-rr/RAD)*TSR
        FARG = min(FARG, 20.0 )   
        FEXP = euler**(-FARG) 
        if FEXP > 1 or FEXP < -1:
            F = 1
        else: 
            F = (2.0/pi) * math.acos(FEXP)
        ADW     =  1.0    /TSR
    VT     = UT     - WT
    QBI = 4.0/BLDS
    PIR = ((pi*rr)**2 + (QBI*RAD*ADW)**2)**0.5
    GAM     = QBI* F*VT                *PIR
    W = (WA**2 + WT**2)**0.5
    RES     = GAM     - 0.5*CHORD* CL*W
    return RES, WA, WT, CL, CD



#Auxiliary
def find_alpha_interval_return_CL_CD(a_list, cl_list, CD_list, alpha):
    #Exception cases
    if len(a_list) == 0:
        #print("Exception case 1: alpha list empty")
        return 0, 1
    if alpha < min(a_list):
        #print("Exception case 2: queried alpha is smaller than the minimum alpha -> returning values for lowest alpha")
        return cl_list[0], CD_list[0]
    if alpha > max(a_list):
        #print("Exception case 3: queried alpha is bigger than the maximum alpha -> returning values for highest alpha")
        return cl_list[-1], CD_list[-1]


    a_list_deducted = [x - alpha for x in a_list]
    for i in range(len(a_list) - 1):
        if a_list_deducted[i]*a_list_deducted[i + 1] < 0: 
            return linear_interpolate(a_list[i], a_list[i + 1], cl_list[i], cl_list[i + 1], alpha), linear_interpolate(a_list[i], a_list[i + 1], CD_list[i], CD_list[i + 1], alpha)
    #print("alpha interval not found")
    return 0, 1

def linear_interpolate(x0, x1, y0, y1, x):
    if x1 - x0 == 0:
        return (y0 + y1)/2
    return y0 + ((x - x0)*(y1 - y0)/(x1 - x0))


def calculate_most_eff_alpha(a_list, cl_list, cd_list):
    a, cl, cd = 0, 0, 1
    clcd_remember = 0
    for i in range(len(a_list)):
        if cd_list[i] == 0:
            continue
        clcd_try = cl_list[i]/cd_list[i]
        if clcd_try > clcd_remember:
            a, cl, cd = a_list[i], cl_list[i], cd_list[i]
            clcd_remember = clcd_try
    return a, cl, cd