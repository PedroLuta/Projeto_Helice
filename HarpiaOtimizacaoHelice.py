import prop_simulate
import prop_design
import optimization_NSGA2
import numpy as np

VelocityVector_m_s = [0, 3, 6, 9, 12, 15, 18]
MotorKV_rpm_V = 450
MotorNoLoadCurrent_A = 1.4
MaxCurrent_A = 54
MotorInternalResistance_ohm = 0.032
BatteryVoltage_V = 22.2

#Genes Contínuos -> Diametro, Croot, Cmax, rmax, Ctip, Passo, Coletivo
#Genes Discretos -> Numero de pas

def Evaluation(Chromossome):
    Diameter_m = Chromossome[0]
    Croot_adim = Chromossome[1]
    Cmax_adim = Chromossome[2]
    rmax_adim = Chromossome[3]
    Ctip_adim = Chromossome[4]
    Pitch_m = Chromossome[5]
    Colective_deg = Chromossome[6]

    #c_R = (a*r**3 + b*r**2 + c*r + d)**0.5
    a = -(2*Cmax_adim**2*rmax_adim - Cmax_adim**2 + Croot_adim**2*rmax_adim**2 - 2*Croot_adim**2*rmax_adim + Croot_adim**2 - Ctip_adim**2*rmax_adim**2)/(rmax_adim**2*(rmax_adim - 1)**2)
    b = -(- 3*Cmax_adim**2*rmax_adim**2 + Cmax_adim**2 - 2*Croot_adim**2*rmax_adim**3 + 3*Croot_adim**2*rmax_adim**2 - Croot_adim**2 + 2*Ctip_adim**2*rmax_adim**3)/(rmax_adim**2*(rmax_adim - 1)**2)
    c = -(3*Cmax_adim**2*rmax_adim - 2*Cmax_adim**2 + Croot_adim**2*rmax_adim**3 - 3*Croot_adim**2*rmax_adim + 2*Croot_adim**2 - Ctip_adim**2*rmax_adim**3)/(rmax_adim*(rmax_adim - 1)**2)
    d = Croot_adim**2

    r_R = np.linspace(0, 1, 11, endpoint = False)
    ChordDistribution_m = [((a*r**3 + b*r**2 + c*r + d)**0.5)*Diameter_m/2 for r in r_R]

    TwistDistribution_deg = prop_design.simple_pitch(r_R*Diameter_m/2, Pitch_m)
    TwistDistributionCollective_deg = [twist + Colective_deg for twist in TwistDistribution_deg]

    #CurrentOutput_A = (BatteryVoltage_V - (RPM/KV_rpm_V))/MotorResistance_ohm
    #PowerOutput_W = BatteryVoltage_V*CurrentOutput_A

def Validation(Chromossome):
    Diameter_m = Chromossome[0]
    Croot_adim = Chromossome[1]
    Cmax_adim = Chromossome[2]
    rmax_adim = Chromossome[3]
    Ctip_adim = Chromossome[4]
    Pitch_m = Chromossome[5]
    Colective_deg = Chromossome[6]

    #c_R = (a*r**3 + b*r**2 + c*r + d)**0.5
    a = -(2*Cmax_adim**2*rmax_adim - Cmax_adim**2 + Croot_adim**2*rmax_adim**2 - 2*Croot_adim**2*rmax_adim + Croot_adim**2 - Ctip_adim**2*rmax_adim**2)/(rmax_adim**2*(rmax_adim - 1)**2)
    b = -(- 3*Cmax_adim**2*rmax_adim**2 + Cmax_adim**2 - 2*Croot_adim**2*rmax_adim**3 + 3*Croot_adim**2*rmax_adim**2 - Croot_adim**2 + 2*Ctip_adim**2*rmax_adim**3)/(rmax_adim**2*(rmax_adim - 1)**2)
    c = -(3*Cmax_adim**2*rmax_adim - 2*Cmax_adim**2 + Croot_adim**2*rmax_adim**3 - 3*Croot_adim**2*rmax_adim + 2*Croot_adim**2 - Ctip_adim**2*rmax_adim**3)/(rmax_adim*(rmax_adim - 1)**2)
    d = Croot_adim**2

    DerivativeZero = c/(2*(d**0.5))
    if DerivativeZero < 0:
        return False
    DerivativeOne = (3*a + 2*b + c)/(2*((a + b + c + d)**0.5))
    if DerivativeOne > 0:
        return False
