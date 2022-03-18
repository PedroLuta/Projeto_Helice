from importing import *
import Auxiliary
#import xfoil_interface

def dT_dQ_from_Cl_Cd(Cl, Cd, phi, V, rho, Blades, chord, rr):
    dT = (rho*Blades*chord)*(V**2)*(Cl*math.cos(phi)-Cd*math.sin(phi))/2
    dQ = (rho*Blades*chord*rr)*(V**2)*(Cl*math.sin(phi)+Cd*math.cos(phi))/2
    return dT, dQ

def correct_for_induction_Ftip(radps, rr, Cl, Cd, Blades, rho, R, chord, vi):
    #pi = nup.pi
    #euler = nup.e
    check = 0
    ai = 0.1
    ai0 = 0.01

    exp_func1 = (-Blades/2)*((R-rr)/rr)
    dT1 = (rho*Blades*chord)/2
    dQ1 = (rho*Blades*chord*rr)/2           
    denominator1 = 4*pi*rr*rho
    denominator2 = 4*pi*(rr**3)*rho

    while True:
        Vr = radps*rr*(1-ai0)
        Vax = vi*(1+ai)
        V = ((Vr**2)+(Vax**2))**0.5

        phi = math.atan(Vax/Vr)
        exp_func = exp_func1*(V/Vax)

        F_tip = (2/pi)*math.acos(euler**exp_func)

        dT = dT1*(V**2)*(Cl*math.cos(phi)-Cd*math.sin(phi))
        dQ = dQ1*(V**2)*(Cl*math.sin(phi)+Cd*math.cos(phi))

        ai_new = dT/(denominator1*(V**2)*(1+ai)*F_tip)
        ai0_new = dQ/(denominator2*V*(1+ai)*radps*F_tip)

        ai_middle = (ai_new + ai)/2
        ai0_middle = (ai0_new + ai0)/2

        if ((abs(ai_middle - ai) < 1/100000) and (abs(ai0_middle - ai0) < 1/100000)) or check > 500:
            return dT, dQ, phi, Vr, Vax

        ai = ai_middle
        ai0 = ai0_middle
        check += 1

def correct_for_induction_Ftip_Beta(radps, rr, a_list, Cl_curve, Cd_curve, Beta, Blades, rho, R, chord, vi):
    #pi = nup.pi
    #euler = nup.e
    check = 0
    ai = 0.1
    ai0 = 0.01

    exp_func1 = (-Blades/2)*((R-rr)/rr)
    dT1 = (rho*Blades*chord)/2
    dQ1 = (rho*Blades*chord*rr)/2           
    denominator1 = 4*pi*rr*rho
    denominator2 = 4*pi*(rr**3)*rho

    while True:
        Vr = radps*rr*(1-ai0)
        Vax = vi*(1+ai)
        V = ((Vr**2)+(Vax**2))**0.5

        phi = math.degrees(math.atan(Vax/Vr))
        alpha = Beta - phi
        Cl, Cd = find_alpha_interval_return_CL_CD(a_list, Cl_curve, Cd_curve, alpha)
        exp_func = exp_func1*(V/Vax)

        F_tip = (2/pi)*math.acos(euler**exp_func)

        dT = dT1*(V**2)*(Cl*math.cos(phi)-Cd*math.sin(phi))
        dQ = dQ1*(V**2)*(Cl*math.sin(phi)+Cd*math.cos(phi))

        ai_new = dT/(denominator1*(V**2)*(1+ai)*F_tip)
        ai0_new = dQ/(denominator2*V*(1+ai)*radps*F_tip)

        ai_middle = (ai_new + ai)/2
        ai0_middle = (ai0_new + ai0)/2

        if ((abs(ai_middle - ai) < 1/100000) and (abs(ai0_middle - ai0) < 1/100000)) or check > 500:
            return dT, dQ, phi, Vr, Vax

        ai = ai_middle
        ai0 = ai0_middle
        check += 1

def correct_for_induction_optimized(radps, rr, Cl, Cd, Blades, rho, chord, vi):
    pi = nup.pi
    check = 0
    ai = 0.1
    ai0 = 0.01

    dT1 = (rho*Blades*chord)/2
    dQ1 = (rho*Blades*chord*rr)/2           
    denominator1 = 4*pi*rr*rho
    denominator2 = 4*pi*(rr**3)*rho

    while True:
        Vr = radps*rr*(1-ai0)
        Vax = vi*(1+ai)
        V = ((Vr**2)+(Vax**2))**0.5
        phi = math.atan(Vax/Vr)

        dT = dT1*(V**2)*(Cl*math.cos(phi)-Cd*math.sin(phi))
        dQ = dQ1*(V**2)*(Cl*math.sin(phi)+Cd*math.cos(phi))

        ai_new = dT/(denominator1*(V**2)*(1+ai))
        ai0_new = dQ/(denominator2*V*(1+ai)*radps)

        ai_middle = (ai_new + ai)/2
        ai0_middle = (ai0_new + ai0)/2

        if ((abs(ai_middle - ai) < 1/100000) and (abs(ai0_middle - ai0) < 1/100000)) or check > 500:
            return dT, dQ, phi, ai0, ai

        ai = ai_middle
        ai0 = ai0_middle
        check += 1

def correct_for_induction_simple(radps, rr, Cl, Cd, Blades, rho, chord, vi):
    pi = nup.pi
    check = 0
    ai = 0.1
    ai0 = 0.01

    while True:
        Vr = radps*rr*(1-ai0)
        Vax = vi*(1+ai)
        V = ((Vr**2)+(Vax**2))**0.5
        phi = math.atan(Vax/Vr)

        dT = (rho*Blades*chord)*(V**2)*(Cl*math.cos(phi)-Cd*math.sin(phi))/2
        dQ = (rho*Blades*chord*rr)*(V**2)*(Cl*math.sin(phi)+Cd*math.cos(phi))/2

        ai_new = dT/(4*pi*rr*rho*(V**2)*(1+ai))
        ai0_new = dQ/(4*pi*(rr**3)*rho*V*(1+ai)*radps)

        ai_middle = (ai_new + ai)/2
        ai0_middle = (ai0_new + ai0)/2

        if ((abs(ai_middle - ai) < 1/100000) and (abs(ai0_middle - ai0) < 1/100000)) or check > 500:
            return dT, dQ, phi, ai0, ai

        ai = ai_middle
        ai0 = ai0_middle
        check += 1
        
def calculate_residual(UA, UT, WZ, PSI, CHORD, CL, rr, BLDS, RAD):
    COSP = math.cos(PSI)
    SINP = math.sin(PSI)
    WA     = 0.5*UA     + 0.5*WZ    *SINP
    WT     = 0.5*UT     + 0.5*WZ    *COSP
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
    WSQ = WA**2 + WT**2
    W = (WSQ)**0.5
    RES     = GAM     - 0.5*CHORD* CL*W
    return RES, WA, WT

def induction_qprop_adapted(OMG, rr, BLDS, CL, RAD, CHORD, VEL):
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
            RESup, _, _ = calculate_residual(UA, UT, WZ, PSIup, CHORD, CL, rr, BLDS, RAD)
            RESlo, _, _ = calculate_residual(UA, UT, WZ, PSIlo, CHORD, CL, rr, BLDS, RAD)
            first = False
        else:
            PSImid = (PSIup + PSIlo)/2
        RESmid, WAmid, WTmid = calculate_residual(UA, UT, WZ, PSImid, CHORD, CL, rr, BLDS, RAD)

        if(abs(PSIup - PSIlo) < EPS):
            return WAmid, WTmid

        if RESup*RESmid < 0:
            RESlo = RESmid
            PSIlo = PSImid
        elif RESlo*RESmid < 0:
            RESup = RESmid
            PSIup = PSImid
        else:
            #print(f"Induction failed, section at radial position {(rr/RAD)*100}% will be assumed as simple flow")
            return UA, UT

def find_alpha_interval_return_CL_CD(a_list, cl_list, CD_list, alpha):
    if len(a_list) == 0:
        return 0, 1000
    if alpha < min(a_list):
        return cl_list[0], CD_list[0]
    if alpha > max(a_list):
        return cl_list[-1], CD_list[-1]
    length = len(a_list)
    for i in range(length - 1):
        if a_list[i]*a_list[i + 1] <= 0:
            return Auxiliary.linear_interpolate(a_list[i], a_list[i + 1], cl_list[i], cl_list[i + 1], alpha), Auxiliary.linear_interpolate(a_list[i], a_list[i + 1], CD_list[i], CD_list[i + 1], alpha)
    return 0, 1000

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
    WSQ = WA**2 + WT**2
    W = (WSQ)**0.5
    RES     = GAM     - 0.5*CHORD* CL*W
    return RES, WA, WT, CL, CD

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


#s = time.time()
jitted_induction_Ftip = njit()(correct_for_induction_Ftip)
jitted_induction_simple = njit()(correct_for_induction_simple)
jitted_induction_optimized = njit()(correct_for_induction_optimized)
#e = time.time()
#
#start1 = time.time()
#_, _, _, _, _ = jitted_induction_simple(1, 1, 1, 1, 1, 1, 1, 1)
#end1 = time.time()
#
#start2 = time.time()
#_, _, _, _, _ = jitted_induction_optimized(1, 1, 1, 1, 1, 1, 1, 1)
#end2 = time.time()
#
#print(f"Normal: {end1 - start1} \nOtimizado: {end2 - start2} \nCompilação: {e - s}")