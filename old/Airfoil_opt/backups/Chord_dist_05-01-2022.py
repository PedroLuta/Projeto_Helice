from importing import *
import xfoil_interface
import induction


def determine_section_Cl(W_times_chord_times_Cl, kvisc, alphas, airfoil, M, Cl_ini = 1):
    Cl_new = 0
    tolerance = 0.01
    tries = 0
    a1, a2, astep = alphas[0], alphas[1], alphas[2]
    Cl_ref = Cl_ini
    while True:
        W_times_chord = W_times_chord_times_Cl/Cl_ref
        Re = ((W_times_chord)/kvisc)
        if Re == 0:
            return 0, 0, 1
        #print(Re)
        alpha_c, Cl_c, Cd_c = xfoil_interface.get_curve_com_default(Re, a1, a2, astep, afile = airfoil, M = M)
        alpha, Cl_new, Cd = xfoil_interface.calculate_most_eff_alpha(alpha_c, Cl_c, Cd_c)
        if (abs(Cl_ref - Cl_new) < tolerance):
            break
        if tries >= 20:
            return 0, 0, 1
        tries += 1
        Cl_ref = (Cl_new + Cl_ref*3)/4
    #print("Done Cl")
    return alpha, Cl_new, Cd



def blade_design_vortex_v1_1(vi, radps, Blades, p, R, airfoil = 'airfoils\\airfoil.txt', sections = 21, rho = 1.225, dvisc = 1.8/100000, alphas = [2, 8, 1], speed_sound = 340, Prescribed_power = 520, Prescribed_thrust = 50, init_disp = 0):
    #From Design of Optimum Propellers, Charles N. Adkins*, Falls Church, Virginia 22042 and Robert H. Liebeckt, Douglas Aircraft Company, Long Beach, California 90846
    
    #With an iterative process to find section Cl

    kvisc = dvisc/rho
    lmb = vi/(radps*R)
    Pc = 2*Prescribed_power/(rho*(vi**3)*pi*(R**2)) 
    Tc = 2*Prescribed_thrust/(rho*(vi**2)*pi*(R**2))
    r_vector = nup.linspace(p*R, R, num = sections)

    disp = init_disp
    tolerance_percentage = 0.001
    disp_new = disp
    disp_new2 = disp
    while True:
        disp = disp_new #(disp_new + disp*2)/3
        #disp = disp_new2 #(disp_new2 + disp*2)/3
        print(f"Try disp = {disp}")
        I1_vec = []
        I2_vec = []
        J1_vec = []
        J2_vec = []
        Beta_vector = []
        chord_vector = []


        # Cl_before = Cl_ref
        # a_before = a_ref
        # Cd_before = Cd_ref
        phi_tip = math.atan(lmb*(1 + (disp/2)))
        for rr in r_vector:
            print("----NEW SECTION----")
            Csi = rr/R
            f = (Blades/2)*(1 - rr/R)/math.sin(phi_tip)
            F = (2/pi)*math.acos(euler**(-f))
            phi = math.atan(math.tan(phi_tip)/Csi)
            speed_ratio = radps*rr/vi
            G = F*math.cos(phi)*math.sin(phi)*speed_ratio
            placeholder_V = ((vi**2) + ((radps*rr)**2))**0.5
            M = placeholder_V/speed_sound
            W_times_chord_times_Cl = 4*pi*lmb*G*vi*R*disp/Blades
            alpha, Cl, Cd = determine_section_Cl(W_times_chord_times_Cl, kvisc, alphas, airfoil, M) 
            # if Cl == 0:
            #     Cl = Cl_before
            #     Cd = Cd_before
            #     alpha = a_before
            # else:
            #     Cl_before = Cl
            #     Cd_before = Cd
            #     a_before = alpha
            if Cl == 0:
                aa = 0
                aa_ = 0
                W = vi*(1 + aa)/(math.sin(phi))
                chord = 0
                Beta = 0
                dI1 = 0
                dI2 = 0
                dJ1 = 0
                dJ2 = 0
            else:
                eps = Cd/Cl
                aa = (disp/2)*(math.cos(phi)**2)*(1 - eps*math.tan(phi)) 
                aa_ = (disp*lmb/2)*math.cos(phi)*math.sin(phi)*(1 + eps/math.tan(phi)) 
                W = vi*(1 + aa)/(math.sin(phi))
                chord = W_times_chord_times_Cl/(Cl*W)
                Beta = math.radians(alpha) + phi
                dI1 = 4*Csi*G*(1 - (eps*math.tan(phi)))
                dI2 = lmb*(dI1/(2*Csi))*(1 + (eps/math.tan(phi)))*math.sin(phi)*math.cos(phi)
                dJ1 = 4*Csi*G*(1 + (eps/math.tan(phi)))
                dJ2 = (dJ1/2)*(1 - (eps*math.tan(phi)))*(math.cos(phi)**2)

            Beta_vector.append(math.degrees(Beta))
            chord_vector.append(chord)
            I1_vec.append(dI1)
            I2_vec.append(dI2)
            J1_vec.append(dJ1)
            J2_vec.append(dJ2)
        I1 = Auxiliary.area_under_curve(r_vector, I1_vec)
        I2 = Auxiliary.area_under_curve(r_vector, I2_vec)
        J1 = Auxiliary.area_under_curve(r_vector, J1_vec)
        J2 = Auxiliary.area_under_curve(r_vector, J2_vec)

        disp_new = -(J1/(2*J2)) + (((J1/(2*J2))**2) + (Pc/J2))**0.5
        disp_new2 = (I1/(2*I2)) - ((((I1/(2*I2))**2) - (Tc/I2))**0.5)
        Tc_after = (I1*disp) - (I2*(disp**2))
        Pc_after = (J1*disp) + (J2*(disp**2))

        if abs((disp/disp_new) - 1) < tolerance_percentage:
            print(f"Final disp: {disp}")
            print(f"Thrust gotten: {Tc_after*(rho*(vi**2)*pi*(R**2))/2}")
            #print(f"Power needed: {Pc_after*(rho*(vi**3)*pi*(R**2))/2}")
            break
        disp = (disp_new + disp*2)/3

    return r_vector, chord_vector, Beta_vector

def blade_design_vortex_v1_2(vi, radps, Blades, p, R, sections = 21, rho = 1.225, dvisc = 1.8/100000, speed_sound = 340, Cl_ref = 1, a_ref = 5, Cd_ref = 0.01, Prescribed_power = 520, Prescribed_thrust = 50, init_disp = 0):
    #From Design of Optimum Propellers, Charles N. Adkins*, Falls Church, Virginia 22042 and Robert H. Liebeckt, Douglas Aircraft Company, Long Beach, California 90846
    
    #With fixed alpha, Cl and Cd defined by user 

    kvisc = dvisc/rho
    lmb = vi/(radps*R)
    Pc = 2*Prescribed_power/(rho*(vi**3)*pi*(R**2)) 
    Tc = 2*Prescribed_thrust/(rho*(vi**2)*pi*(R**2))
    r_vector = nup.linspace(p*R, R, num = sections)

    disp = init_disp
    tolerance_percentage = 0.00000001
    disp_new = disp
    disp_new2 = disp
    while True:
        disp = disp_new #(disp_new + disp*2)/3
        #disp = disp_new2 #(disp_new2 + disp*2)/3
        print(f"Try disp = {disp}")
        I1_vec = []
        I2_vec = []
        J1_vec = []
        J2_vec = []
        Beta_vector = []
        chord_vector = []

        phi_tip = math.atan(lmb*(1 + (disp/2)))
        for rr in r_vector:
            print("----NEW SECTION----")
            Csi = rr/R
            f = (Blades/2)*(1 - Csi)/math.sin(phi_tip)
            F = (2/pi)*math.acos(euler**(-f))
            phi = math.atan(math.tan(phi_tip)/Csi)
            speed_ratio = radps*rr/vi
            G = F*math.cos(phi)*math.sin(phi)*speed_ratio #CHECK PAPER
            W_times_chord = 4*pi*lmb*G*vi*R*disp/(Cl_ref*Blades)
            Re = ((W_times_chord)/kvisc)
            print(f"Reynolds of section: {Re}")
            eps = Cd_ref/Cl_ref
            aa = (disp/2)*(math.cos(phi)**2)*(1 - eps*math.tan(phi)) 
            aa_ = (disp*lmb/2)*math.cos(phi)*math.sin(phi)*(1 + eps/math.tan(phi)) 
            print(f"phi:{math.degrees(phi)}\ncalculated phi:{math.degrees(math.atan(vi*(1 + aa)/(radps*rr*(1 - aa_))))}")
            W = vi*(1 + aa)/(math.sin(phi))
            # print(W)
            # print(aa)
            chord = W_times_chord/W
            Beta = math.radians(a_ref) + phi
            dI1 = 4*Csi*G*(1 - (eps*math.tan(phi)))
            dI2 = lmb*(dI1/(2*rr/R))*(1 + (eps/math.tan(phi)))*math.sin(phi)*math.cos(phi)
            dJ1 = 4*Csi*G*(1 + (eps/math.tan(phi)))
            dJ2 = (dJ1/2)*(1 - (eps*math.tan(phi)))*(math.cos(phi)**2)

            Beta_vector.append(math.degrees(Beta))
            chord_vector.append(chord)
            I1_vec.append(dI1)
            I2_vec.append(dI2)
            J1_vec.append(dJ1)
            J2_vec.append(dJ2)
        I1 = Auxiliary.area_under_curve(r_vector, I1_vec)
        I2 = Auxiliary.area_under_curve(r_vector, I2_vec)
        J1 = Auxiliary.area_under_curve(r_vector, J1_vec)
        J2 = Auxiliary.area_under_curve(r_vector, J2_vec)

        disp_new = -(J1/(2*J2)) + ((((J1/(2*J2))**2) + (Pc/J2))**0.5)
        disp_new2 = (I1/(2*I2)) - ((((I1/(2*I2))**2) - (Tc/I2))**0.5)
        Tc_after = (I1*disp) - (I2*(disp**2))
        Pc_after = (J1*disp) + (J2*(disp**2))

        if abs((disp/disp_new) - 1) < tolerance_percentage:
            print(f"Final disp: {disp}")
            print(f"Thrust gotten: {Tc_after*(rho*(vi**2)*pi*(R**2))/2}")
            #print(f"Power needed: {Pc_after*(rho*(vi**3)*pi*(R**2))/2}")
            break

    return r_vector, chord_vector, Beta_vector

# def blade_design_vortex_v1_3(vi, radps, Blades, p, R, airfoil = 'airfoils\\airfoil.txt', sections = 21, rho = 1.225, dvisc = 1.8/100000, alphas = [2, 8, 1], speed_sound = 340, Prescribed_power = 520, Prescribed_thrust = 50, chord_at_75 = 0.06, init_disp = 0):
#     #From Design of Optimum Propellers, Charles N. Adkins*, Falls Church, Virginia 22042 and Robert H. Liebeckt, Douglas Aircraft Company, Long Beach, California 90846
    
#     #With only one XFoil query, at 75% of radius

#     kvisc = dvisc/rho
#     Vr_at_75 = radps*0.75*R
#     V_at_75 = ((Vr_at_75**2)+(vi**2))**0.5
#     M_at_75 = V_at_75/speed_sound
#     Re_at_75 = ((V_at_75*chord_at_75)/kvisc) 
#     alpha_c, Cl_c, Cd_c = xfoil_interface.get_curve_com_default(Re_at_75, alphas[0], alphas[1], alphas[2], afile = airfoil, M = M_at_75)
#     a_ref, Cl_ref, Cd_ref = xfoil_interface.calculate_most_eff_alpha(alpha_c, Cl_c, Cd_c)
#     lmb = vi/(radps*R)

#     Pc = 2*Prescribed_power/(rho*(vi**3)*pi*(R**2)) 
#     Tc = 2*Prescribed_thrust/(rho*(vi**2)*pi*(R**2))
#     r_vector = nup.linspace(p*R, R, num = sections)

#     disp = init_disp
#     tolerance_percentage = 0.00000001
#     disp_new = disp
#     disp_new2 = disp
#     while True:
#         disp = disp_new #(disp_new + disp*2)/3
#         #disp = disp_new2 #(disp_new2 + disp*2)/3
#         print(f"Try disp = {disp}")
#         I1_vec = []
#         I2_vec = []
#         J1_vec = []
#         J2_vec = []
#         Beta_vector = []
#         chord_vector = []

#         phi_tip = math.atan(lmb*(1 + (disp/2)))
#         for rr in r_vector:
#             print("----NEW SECTION----")
#             Csi = rr/R
#             f = (Blades/2)*(1 - Csi)/math.sin(phi_tip)
#             F = (2/pi)*math.acos(euler**(-f))
#             phi = math.atan(math.tan(phi_tip)/Csi)
#             speed_ratio = radps*rr/vi
#             G = F*math.cos(phi)*math.sin(phi)*speed_ratio #CHECK PAPER
#             W_times_chord = 4*pi*lmb*G*vi*R*disp/(Cl_ref*Blades)
#             Re = ((W_times_chord)/kvisc)
#             print(f"Reynolds of section: {Re}")
#             eps = Cd_ref/Cl_ref
#             aa = (disp/2)*(math.cos(phi)**2)*(1 - eps*math.tan(phi)) 
#             aa_ = (disp*lmb/2)*math.cos(phi)*math.sin(phi)*(1 + eps/math.tan(phi)) 
#             W = vi*(1 + aa)/(math.sin(phi))
#             chord = W_times_chord/W
#             Beta = math.radians(a_ref) + phi
#             dI1 = 4*Csi*G*(1 - (eps*math.tan(phi)))
#             dI2 = lmb*(dI1/(2*Csi))*(1 + (eps/math.tan(phi)))*math.sin(phi)*math.cos(phi)
#             dJ1 = 4*Csi*G*(1 + (eps/math.tan(phi)))
#             dJ2 = (dJ1/2)*(1 - (eps*math.tan(phi)))*(math.cos(phi)**2)

#             Beta_vector.append(math.degrees(Beta))
#             chord_vector.append(chord)
#             I1_vec.append(dI1)
#             I2_vec.append(dI2)
#             J1_vec.append(dJ1)
#             J2_vec.append(dJ2)
#         I1 = Auxiliary.area_under_curve(r_vector, I1_vec)
#         I2 = Auxiliary.area_under_curve(r_vector, I2_vec)
#         J1 = Auxiliary.area_under_curve(r_vector, J1_vec)
#         J2 = Auxiliary.area_under_curve(r_vector, J2_vec)

#         disp_new = -(J1/(2*J2)) + ((((J1/(2*J2))**2) + (Pc/J2))**0.5)
#         disp_new2 = (I1/(2*I2)) - ((((I1/(2*I2))**2) - (Tc/I2))**0.5)
#         Tc_after = (I1*disp) - (I2*(disp**2))
#         Pc_after = (J1*disp) + (J2*(disp**2))

#         if abs((disp/disp_new) - 1) < tolerance_percentage:
#             print(f"Final disp: {disp}")
#             print(f"Thrust gotten: {Tc_after*(rho*(vi**2)*pi*(R**2))/2}")
#             #print(f"Power needed: {Pc_after*(rho*(vi**3)*pi*(R**2))/2}")
#             break

#     return r_vector, chord_vector, Beta_vector







# def blade_design_vortex_v2(vi, radps, Blades, p, R, airfoil = 'airfoils\\airfoil.txt', sections = 21, rho = 1.225, dvisc = 1.8/100000, alphas = [2, 8, 1], speed_sound = 340, Cl_ref = 1, Prescribed_power = 520):
#     #Adapted from Design of Optimum Propellers, Charles N. Adkins*, Falls Church, Virginia 22042 and Robert H. Liebeckt, Douglas Aircraft Company, Long Beach, California 90846
#     kvisc = dvisc/rho
#     a1, a2, astep = alphas[0], alphas[1], alphas[2]
#     chord_at_75 = 0.06
#     Vr_at_75 = radps*0.75*R
#     V_at_75 = ((Vr_at_75**2)+(vi**2))**0.5
#     M_at_75 = V_at_75/speed_sound
#     Re_at_75 = ((V_at_75*chord_at_75)/kvisc) 
#     alpha_c, Cl_c, Cd_c = xfoil_interface.get_curve_com_default(Re_at_75, a1, a2, astep, afile = airfoil, M = M_at_75)
#     alpha_ref, Cl_ref, Cd_ref = xfoil_interface.calculate_most_eff_alpha(alpha_c, Cl_c, Cd_c)
#     WA, WT = induction.induction_qprop_adapted(radps, 0.75*R, Blades, Cl_ref, R, chord_at_75, vi)
#     disp_v1 = (WA - vi)/vi
#     print(disp_v1)
#     disp_v2 = WA - vi
#     lmb = vi/(radps*R)


#     r_vector = nup.linspace(p*R, R, num = sections)

#     Beta_vector = []
#     chord_vector = []
#     for rr in r_vector:
#         print("----NEW SECTION----")
#         Vr = radps*rr
#         phi_t = math.atan(lmb*(1 + (disp_v1/2)))
#         f = (Blades/2)*(1 - rr/R)/math.sin(phi_t)
#         if f > 1 or f < -1:
#             F = 1
#         else: 
#             F = (2/pi)*math.acos(euler**(-f))
#         phi = math.atan(math.tan(phi_t)/(rr/R))
#         speed_ratio = Vr/vi
#         G = F*math.cos(phi)*math.sin(phi)*speed_ratio
#         placeholder_V = ((vi**2) + (Vr**2))**0.5
#         M = placeholder_V/speed_sound
#         alpha, Cl, Cd = determine_section_Cl(Cl_ref, disp_v1, Blades, lmb, G, kvisc, R, vi, a1, a2, astep, airfoil, M)
#         eps = Cl/Cd
#         aa = (disp_v1/2)*(math.cos(phi)**2)*(1 - eps*math.tan(phi)) #(disp/2)*(math.cos(phi*(1 - eps*math.tan(phi)))**2)
#         print(aa)
#         W = vi*(1 + aa)/(math.sin(phi))
#         print(lmb)
#         print(G)
#         print(disp_v1)
#         print(Cl)
#         print(W)
#         chord = (4*pi*lmb*G*vi*R*disp_v1/(Cl*Blades)) / W


#         Beta = math.radians(alpha) + phi
#         aa_ = (disp_v1*lmb/2)*math.cos(phi)*math.sin(phi)*(1 + eps/math.tan(phi)) #(disp*lmb/2)*math.cos(phi)*math.sin(phi*(1 + eps/math.tan(phi)))
#         Vrr = Vr*(1 - aa_)
#         Vax = vi*(1 + aa)
#         Beta_vector.append(math.degrees(Beta))
#         chord_vector.append(chord)
#     WA, WT = W*math.sin(phi), W*math.cos(phi)

#     return chord_vector, Beta_vector






# def blade_design_vortex_v3(vi, radps, Blades, p, R, airfoil = 'airfoils\\airfoil.txt', sections = 21, rho = 1.225, dvisc = 1.8/100000, alphas = [2, 8, 1], speed_sound = 340, Cl_ref = 1, Prescribed_power = 520):
#     #Adapted from Design of Optimum Propellers, Charles N. Adkins*, Falls Church, Virginia 22042 and Robert H. Liebeckt, Douglas Aircraft Company, Long Beach, California 90846
#     kvisc = dvisc/rho
#     a1, a2, astep = alphas[0], alphas[1], alphas[2]
#     chord_at_75 = 0.06
#     Vr_at_75 = radps*0.75*R
#     V_at_75 = ((Vr_at_75**2)+(vi**2))**0.5
#     M_at_75 = V_at_75/speed_sound
#     Re_at_75 = ((V_at_75*chord_at_75)/kvisc) 
#     alpha_c, Cl_c, Cd_c = xfoil_interface.get_curve_com_default(Re_at_75, a1, a2, astep, afile = airfoil, M = M_at_75)
#     alpha_ref, Cl_ref, Cd_ref = xfoil_interface.calculate_most_eff_alpha(alpha_c, Cl_c, Cd_c)
#     WA, WT = induction.induction_qprop_adapted(radps, 0.75*R, Blades, Cl_ref, R, chord_at_75, vi)
#     disp_v2 = WA
#     lmb = vi/(radps*R)
#     tolerance = 0.001
#     print(disp_v2)


#     r_vector = nup.linspace(p*R, R, num = sections)

#     Beta_vector = []
#     chord_vector = []
#     wa_before = 0
#     for rr in r_vector:
#         print("----NEW SECTION----")
#         delta = 0.01
#         chord = 0.0001
#         while True:
#             Vr = radps*rr
#             V = ((Vr**2)+(vi**2))**0.5
#             M = V/speed_sound
#             Re = ((V*chord)/kvisc) 
#             alpha_c, Cl_c, Cd_c = xfoil_interface.get_curve_com_default(Re, a1, a2, astep, afile = airfoil, M = M)
#             alpha, Cl, Cd = xfoil_interface.calculate_most_eff_alpha(alpha_c, Cl_c, Cd_c)
#             WA, WT = induction.induction_qprop_adapted(radps, rr, Blades, Cl, R, chord, vi)
#             print(f"chord: {chord}")
#             print(WA - disp_v2)
#             if abs(WA - disp_v2) < tolerance:
#                 break
#             if wa_before*WA < 0:
#                 print("Delta diminished")
#                 chord -= delta
#                 delta /= 10
#             chord += delta
#             wa_before = WA

#         phi = math.atan(WA/WT)
#         Beta = math.radians(alpha) + phi
#         Beta_vector.append(math.degrees(Beta))
#         chord_vector.append(chord)

#     return chord_vector, Beta_vector

