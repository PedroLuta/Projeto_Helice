import numpy as nup
import math
import xfoil_interface

euler = nup.e
pi = nup.pi

def chord_distribution(x, factor = 1):
    if x > 0.4:
        return factor*(-1.8598247908*(x**4) + 4.7243249915*(x**3) - 4.4099344990*(x**2) + 1.5997758465*x + 0.0076143297)
    else:
        return factor*(-775.7226922959*(x**6) + 1170.2815273553*(x**5) - 660.9381028488*(x**4) + 161.9663843811*(x**3) - 13.9599756954*(x**2) + 0.1950071591*x + 0.1005241407)

def fixed_pitch_qprop_v1(vi, radps, Blades, p, R, Beta_dist, airfoil = 'airfoils\\airfoil.txt', sections = 21, rho = 1.225, dvisc = 1.8/100000):
    #CHORD DISTRIBUTION DEFINED BY A FUNCTION
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
        chord = R*chord_distribution(rr/R)

        Vr = radps*rr
        V = ((Vr**2)+(vi**2))**0.5
        Re = ((V*chord)/kvisc) 

        alpha_c, Cl_c, Cd_c = xfoil_interface.get_curve_com_default(Re, -5, 10, 1, afile = airfoil)
        WA, WT, Cl, Cd = induction_qprop_list(radps, rr, Blades, alpha_c, Cl_c, Cd_c, Beta, R, chord, vi)
        W = (WA**2 + WT**2)**0.5
        phi = math.atan(WA/WT)
        dT = (rho*Blades*chord)*(W**2)*(Cl*math.cos(phi)-Cd*math.sin(phi))/2
        dQ = (rho*Blades*chord*rr)*(W**2)*(Cl*math.sin(phi)+Cd*math.cos(phi))/2


        dQ_vector.append(dQ)
        dT_vector.append(dT)
        Re_vector.append(Re)
        WA_vector.append(WA)
        WT_vector.append(WT)
        Cl_vector.append(Cl)
        Cd_vector.append(Cd)
    return dT_vector, dQ_vector, r_vector, Re_vector, WA_vector, Cl_vector, Cd_vector

def fixed_pitch_qprop_v2(vi, radps, Blades, R, r_vector, Beta_dist, chord_dist, airfoil = 'airfoils\\airfoil.txt', rho = 1.225, dvisc = 1.8/100000):
    #CHORD DISTRIBUTION TAKEN AS INPUT
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

        alpha_c, Cl_c, Cd_c = xfoil_interface.get_curve_com_default(Re, -5, 10, 0.5, afile = airfoil)
        WA, WT, Cl, Cd = induction_qprop_list(radps, rr, Blades, alpha_c, Cl_c, Cd_c, Beta, R, chord, vi)
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














def induction_qprop_list(OMG, rr, BLDS, a_list, CL_list, CD_list, Beta, RAD, CHORD, VEL):
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
            RESup, _, _, _, _ = calculate_residual_curve(UA, UT, WZ, Beta, PSIup, CHORD, a_list, CL_list, CD_list, rr, BLDS, RAD)
            RESlo, _, _, _, _ = calculate_residual_curve(UA, UT, WZ, Beta, PSIlo, CHORD, a_list, CL_list, CD_list, rr, BLDS, RAD)
            first = False
        else:
            PSImid = (PSIup + PSIlo)/2

        RESmid, WAmid, WTmid, CL, CD = calculate_residual_curve(UA, UT, WZ, Beta, PSImid, CHORD, a_list, CL_list, CD_list, rr, BLDS, RAD)

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

def calculate_residual_curve(UA, UT, WZ, Beta, PSI, CHORD, a_list, CL_list, CD_list, rr, BLDS, RAD):
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
