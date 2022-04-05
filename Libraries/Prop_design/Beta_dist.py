from importing import *
import xfoil_interface
import prop_simulate
#import Auxiliary

def qprop_highest_efficiency(vi, radps, Blades, R, r_vector, chord_dist, airfoil = 'airfoils\\airfoil.txt', rho = 1.225, dvisc = 1.8/100000, alphas = [2, 8, 1]):
    a1, a2, astep = alphas[0], alphas[1], alphas[2]
    kvisc = dvisc/rho

    Beta_vector = []

    for i in range(len(r_vector)):
        rr = r_vector[i]
        chord = chord_dist[i]

        Vr = radps*rr
        V = ((Vr**2)+(vi**2))**0.5
        Re = ((V*chord)/kvisc) 

        alpha_c, Cl_c, Cd_c = xfoil_interface.get_curve_com_default(Re, a1, a2, astep, afile = airfoil)
        alpha, Cl, Cd = calculate_most_eff_alpha(alpha_c, Cl_c, Cd_c)

        WA, WT = induction_qprop_original_adapted(radps, rr, Blades, Cl, R, chord, vi)
        phi = math.atan(WA/WT)
        Beta_vector.append(alpha + math.degrees(phi))

    return Beta_vector

def momentum_Ftip_highest_efficiency(vi, radps, Blades, R, r_vector, chord_dist, airfoil = 'airfoils\\airfoil.txt', rho = 1.225, dvisc = 1.8/100000, alphas = [2, 8, 1]):
    a1, a2, astep = alphas[0], alphas[1], alphas[2]
    kvisc = dvisc/rho

    Beta_vector = []

    for i in range(len(r_vector)):
        rr = r_vector[i]
        chord = chord_dist[i]

        Vr = radps*rr
        V = ((Vr**2)+(vi**2))**0.5
        Re = ((V*chord)/kvisc)

        alpha_c, Cl_c, Cd_c = xfoil_interface.get_curve_com_default(Re, a1, a2, astep, afile = airfoil)
        alpha, Cl, Cd = calculate_most_eff_alpha(alpha_c, Cl_c, Cd_c)
        try:
            _, _, phi, _, _ = jitted_induction_momentum_Ftip(radps, rr, Cl, Cd, Blades, rho, R, chord, vi)
        except:
            _, _, phi, _, _ = 0, 0, 0, 0, 0 
        Beta_vector.append(alpha + math.degrees(phi))

    return Beta_vector

def simple_pitch_inches1(r_vector_inches, Pitch_inches): #r_vector and pitch in inches
    r_vector = [rr*0.0254 for rr in r_vector_inches]
    Pitch = Pitch_inches*0.0254
    return simple_pitch(r_vector, Pitch)

def simple_pitch_inches2(r_vector, Pitch_inches): #only pitch in inches
    Pitch = Pitch_inches*0.0254
    return simple_pitch(r_vector, Pitch)

def simple_pitch(r_vector, Pitch): #All values in meters
    Beta_dist = []
    for rr in r_vector:
        Beta_dist.append(math.degrees(math.atan(Pitch/(2*pi*rr))))
    return Beta_dist

def simple_pitch_highest_efficiency(vi, radps, r_vector, Blades, R, chord_dist, min_pitch, max_pitch, quantity_to_test,  airfoil = 'airfoils\\airfoil.txt', rho = 1.225, dvisc = 1.8/100000):
    kvisc = dvisc/rho

    airfoil_data = []
    for i in range(len(r_vector)):
        rr = r_vector[i]
        chord = chord_dist[i]

        Vr = radps*rr
        V = ((Vr**2)+(vi**2))**0.5
        Re = ((V*chord)/kvisc) 

        alpha_c, Cl_c, Cd_c = xfoil_interface.get_curve_com_default(Re, -5, 10, 0.5, afile = airfoil)

        airfoil_data.append([alpha_c, Cl_c, Cd_c])
    #print("Airfoil Data acquired")

    best_eff = 0
    best_pitch = []
    for pitch in nup.linspace(min_pitch, max_pitch, num = quantity_to_test):
        Beta_dist = simple_pitch_inches2(r_vector, pitch)
        #print(f"Trying pitch {pitch}")
        dT_vector = []
        dQ_vector = []
        for i in range(len(r_vector)):
            rr = r_vector[i]
            Beta = Beta_dist[i]
            chord = chord_dist[i]

            alpha_c, Cl_c, Cd_c = airfoil_data[i]

            WA, WT, Cl, Cd = induction_qprop_fixed_pitch(radps, rr, Blades, alpha_c, Cl_c, Cd_c, Beta, R, chord, vi)
            W = (WA**2 + WT**2)**0.5
            phi = math.atan(WA/WT)
            dT = (rho*Blades*chord)*(W**2)*(Cl*math.cos(phi) - Cd*math.sin(phi))/2
            dQ = (rho*Blades*chord*rr)*(W**2)*(Cl*math.sin(phi) + Cd*math.cos(phi))/2

            dQ_vector.append(dQ)
            dT_vector.append(dT)
        Thrust = Auxiliary.area_under_curve(r_vector, dT_vector)
        Torque = Auxiliary.area_under_curve(r_vector, dQ_vector)
        n = radps/(2*pi)
        Ct = Thrust/(rho*(n**2)*((2*R)**4))
        Cp = Torque*radps/(rho*(n**3)*((2*R)**5))
        J = vi/(n*2*R)
        eff = Ct*J/Cp
        if eff > best_eff:
            #print(f"-------------Best pitch changed to {pitch}---------------")
            best_eff = eff
            best_pitch = Beta_dist
    return best_pitch
        






#Induction methods
def induction_qprop_original_adapted(OMG, rr, BLDS, CL, RAD, CHORD, VEL):
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
            RESup, _, _ = calculate_residual_original(UA, UT, WZ, PSIup, CHORD, CL, rr, BLDS, RAD)
            RESlo, _, _ = calculate_residual_original(UA, UT, WZ, PSIlo, CHORD, CL, rr, BLDS, RAD)
            first = False
        else:
            PSImid = (PSIup + PSIlo)/2
        RESmid, WAmid, WTmid = calculate_residual_original(UA, UT, WZ, PSImid, CHORD, CL, rr, BLDS, RAD)

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

def induction_momentum_Ftip(radps, rr, Cl, Cd, Blades, rho, R, chord, vi):
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

jitted_induction_momentum_Ftip = njit()(induction_momentum_Ftip)



#Residuals
def calculate_residual_original(UA, UT, WZ, PSI, CHORD, CL, rr, BLDS, RAD):
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