import random
import numpy as np
import matplotlib.pyplot as plt
import prop_design

# r_R = np.linspace(0, 1, 101, endpoint = True)
# # ChordVec_m_R = [-1]
# # i = 0
# # while i < 50:
# #     a = random.uniform(-1, 5)
# #     b = random.uniform(-5, 5)
# #     c = random.uniform(-5, 5)



# #     ChordSquaredVec_m2_R2 = [(a*(r**2)) + (b*r) + c for r in r_R]
# #     if any(chord < 0 for chord in ChordSquaredVec_m2_R2):
# #         continue
# #     ChordVec_m_R = [ChordSquared**0.5 for ChordSquared in ChordSquaredVec_m2_R2]
# #     if any(chord < 0.05 for chord in ChordVec_m_R) or any(chord > 0.5 for chord in ChordVec_m_R):
# #         continue
# #     plt.plot(r_R, ChordVec_m_R)
# #     i += 1

# # print(a)
# # print(b)
# # print(c)
# # plt.ylim([0,1])
# # plt.show()
# Ctip_adim = 0.1
# Croot_adim = 0.2
# Cmax_adim = 0.3
# rmax_adim = 0.7

# # b = (((Croot_adim**2)*((rmax_adim**2) - 1)) - ((Ctip_adim**2)*(rmax_adim**2)) + (Cmax_adim**2))/(rmax_adim*(1 - rmax_adim))
# # a = (Ctip_adim**2) - (Croot_adim**2) - b
# # c = Croot_adim**2

# # plt.plot(r_R, [((a*(r**2)) + (b*r) + c)**0.5 for r in r_R])
# # # plt.ylim([0,1])
# # print([((a*(r**2)) + (b*r) + c)**0.5 for r in r_R])
# i = 0
# while i < 50:
#     Ctip_adim = random.uniform(0, 0.3)
#     Croot_adim = random.uniform(0, 0.3)
#     Cmax_adim = random.uniform(0, 0.5)
#     rmax_adim = random.uniform(0, 1)

#     a = -(2*Cmax_adim**2*rmax_adim - Cmax_adim**2 + Croot_adim**2*rmax_adim**2 - 2*Croot_adim**2*rmax_adim + Croot_adim**2 - Ctip_adim**2*rmax_adim**2)/(rmax_adim**2*(rmax_adim - 1)**2)
#     b = -(- 3*Cmax_adim**2*rmax_adim**2 + Cmax_adim**2 - 2*Croot_adim**2*rmax_adim**3 + 3*Croot_adim**2*rmax_adim**2 - Croot_adim**2 + 2*Ctip_adim**2*rmax_adim**3)/(rmax_adim**2*(rmax_adim - 1)**2)
#     c = -(3*Cmax_adim**2*rmax_adim - 2*Cmax_adim**2 + Croot_adim**2*rmax_adim**3 - 3*Croot_adim**2*rmax_adim + 2*Croot_adim**2 - Ctip_adim**2*rmax_adim**3)/(rmax_adim*(rmax_adim - 1)**2)
#     d = Croot_adim**2

#     DerivativeZero = c/(2*(d**0.5))
#     if DerivativeZero < 0:
#         continue
#     DerivativeOne = (3*a + 2*b + c)/(2*((a + b + c + d)**0.5))
#     if DerivativeOne > 0:
#         continue

#     ChordVec_m_R = [((a*(r**3)) + (b*(r**2)) + (c*r) + d)**0.5 for r in r_R]
#     plt.plot(r_R, ChordVec_m_R)
#     i += 1

# plt.show()

def MotorTorque(Omega_rad_s, Voltage_V):
    #inputs: Omega, Voltagem, Parametros do motor
    #outputs: Torque, Torque_OMG, Torque_VOLT, Current, Current_OMG, Current_VOLT
    MotorKV_rad_sV = MotorKV_rpm_V*2*np.pi/60

    MotorVoltage_V = Omega_rad_s/MotorKV_rad_sV
    dMotorVoltageByOmega_Vs_rad = 1/MotorKV_rad_sV

    MotorCurrent_A = (Voltage_V - MotorVoltage_V)/MotorResistance_ohm
    dMotorCurrentByOmega_As_rad = - dMotorVoltageByOmega_Vs_rad/MotorResistance_ohm
    dMotorCurrentByVoltage_A_V = 1/MotorResistance_ohm

    MotorTorque_Nm = (MotorCurrent_A - MotorNoLoadCurrent_A)/MotorKV_rad_sV
    dMotorTorqueByOmega_sNm_rad = dMotorCurrentByOmega_As_rad/MotorKV_rad_sV
    dMotorTorqueByVoltage_Nm_V = dMotorCurrentByVoltage_A_V/MotorKV_rad_sV
    
    return MotorTorque_Nm, dMotorTorqueByOmega_sNm_rad, dMotorTorqueByVoltage_Nm_V, MotorCurrent_A, dMotorCurrentByOmega_As_rad, dMotorCurrentByVoltage_A_V

def MotorVoltage(Omega_rad_s, PropellerTorque_Nm):
    #inputs: Omega, Torque, Parametros do motor
    #outputs: Volt, Volt_OMG, Volt_TORQUE, Amp, Amp_OMG, Amp_TORQUE
    MotorKV_rad_sV = MotorKV_rpm_V*2*np.pi/60

    MotorCurrent_A = (PropellerTorque_Nm*MotorKV_rad_sV) + MotorNoLoadCurrent_A
    Voltage_V = (MotorCurrent_A*MotorResistance_ohm) + (Omega_rad_s/MotorKV_rad_sV)

    for _ in range(50):
        MotorTorque_Nm, dMotorTorqueByOmega_sNm_rad, dMotorTorqueByVoltage_Nm_V, MotorCurrent_A, dMotorCurrentByOmega_As_rad, dMotorCurrentByVoltage_A_V = MotorTorque(Omega_rad_s, Voltage_V)
        Residue = MotorTorque_Nm - PropellerTorque_Nm
        dResidueByVoltage_1_V = dMotorTorqueByVoltage_Nm_V
        dVoltage_V = -Residue/dResidueByVoltage_1_V
        if dVoltage_V < (10**(-8))*max(1, abs(Voltage_V)):
            break
        Voltage_V += dVoltage_V
        
    dResidueByOmega_s_rad = dMotorTorqueByOmega_sNm_rad
    dResidueByTorque_1_Nm = -1

    dVoltageByOmega_sV_rad = -dResidueByOmega_s_rad/dResidueByVoltage_1_V
    dVoltageByTorque_V_Nm = -dResidueByTorque_1_Nm/dResidueByVoltage_1_V

    dCurrentByOmega_sA_rad = (dMotorCurrentByVoltage_A_V*dVoltageByOmega_sV_rad) + dMotorCurrentByOmega_As_rad
    dCurrentByTorque_A_Nm = dMotorCurrentByVoltage_A_V*dVoltageByTorque_V_Nm

    return Voltage_V, dVoltageByOmega_sV_rad, dVoltageByTorque_V_Nm, MotorCurrent_A, dCurrentByOmega_sA_rad, dCurrentByTorque_A_Nm


VelocityVector_m_s = [0, 3, 6, 9, 12, 15, 18]
MotorKV_rpm_V = 450
MotorNoLoadCurrent_A = 1.4
MaxCurrent_A = 54
MotorResistance_ohm = 0.032
BatteryVoltage_V = 22.2 #DEIXAR VOLTAGEM LIVRE -> SE PASSAR DE 22.2 DESCARTA
Power_W = 700

CurrentOutput_A = Power_W/BatteryVoltage_V

Diameter_m = 3
Croot_adim = 0.2
Cmax_adim = 0.3
rmax_adim = 0.5
Ctip_adim = 0.2
Pitch_m = 2
Colective_deg = 1

#c_R = (a*r**3 + b*r**2 + c*r + d)**0.5
a = -(2*Cmax_adim**2*rmax_adim - Cmax_adim**2 + Croot_adim**2*rmax_adim**2 - 2*Croot_adim**2*rmax_adim + Croot_adim**2 - Ctip_adim**2*rmax_adim**2)/(rmax_adim**2*(rmax_adim - 1)**2)
b = -(- 3*Cmax_adim**2*rmax_adim**2 + Cmax_adim**2 - 2*Croot_adim**2*rmax_adim**3 + 3*Croot_adim**2*rmax_adim**2 - Croot_adim**2 + 2*Ctip_adim**2*rmax_adim**3)/(rmax_adim**2*(rmax_adim - 1)**2)
c = -(3*Cmax_adim**2*rmax_adim - 2*Cmax_adim**2 + Croot_adim**2*rmax_adim**3 - 3*Croot_adim**2*rmax_adim + 2*Croot_adim**2 - Ctip_adim**2*rmax_adim**3)/(rmax_adim*(rmax_adim - 1)**2)
d = Croot_adim**2

r_R = np.linspace(0, 1, 11, endpoint = False)
ChordDistribution_m = [((a*r**3 + b*r**2 + c*r + d)**0.5)*Diameter_m/2 for r in r_R]
TwistDistribution_deg = prop_design.simple_pitch(r_R*Diameter_m/2, Pitch_m)
TwistDistributionCollective_deg = [twist + Colective_deg for twist in TwistDistribution_deg]
for AxialVelocity_m_s in [0, 3, 6, 9, 12, 15]:
    Omega_rad_s = 1
    PropellerTorque_Nm = 0
    MotorVoltage_V, _, _, _, _, _ =  MotorVoltage(Omega_rad_s, PropellerTorque_Nm)
    for _ in range(50):
        #calcular torque requerido pela hélice no OMG atual
        #MotorTorque(Omega_rad_s, Voltage)
        #Resíduos
        #Se voltagem for maior que a da bateria, fazer algo
        #Iterar
        pass

#CurrentOutput_A = (BatteryVoltage_V - (RPM/KV_rpm_V))/MotorResistance_ohm

    