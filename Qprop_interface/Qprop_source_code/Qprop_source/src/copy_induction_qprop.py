import math
import numpy as nup

#INPUTS
#CHORD   local blade chord
#BETA    local blade angle  (radians)
#rr       local radius
#BLDS    number of blades
#RAD     tip radius  (for Prandtl's factor)
#VEL     forward flight speed
#OMG     rotational speed  (radians/time)
#VSO     speed of sound
#MCRIT

PI = nup.pi
euler = nup.e
EPS = 1.0E-6 
MSQMAX = 0.9

rr = 0.1
VEL = 10
rpm = 2300
OMG = rpm*2*PI/60
BLDS = 4
RAD = 0.3
BETA = math.radians(30)
VSO = 340
CHORD = 0.05

CL0 = 0.56 
DCLDA = 5.58 
CLMIN = -0.3 
CLMAX = 1.37 
MCRIT = 0.7


U0A = 0.
U0T = 0.
UA     = VEL   + U0A
UT     = OMG*rr - U0T

WZ = (UA**2 + UT**2)**0.5

PSI2 = BETA + CL0/DCLDA
PSI = math.atan(UA/UT)
PSI = max(PSI, PSI2)



while True:
    COSP = math.cos(PSI)
    SINP = math.sin(PSI)

    WA     = 0.5*UA     + 0.5*WZ    *SINP
    WA_PSI =              0.5*WZ    *COSP
    WT     = 0.5*UT     + 0.5*WZ    *COSP
    WT_PSI =            - 0.5*WZ    *SINP

    if (WA == 0.0):
        F     = 1.0
        F_PSI = 0
        ADW     = 0
        ADW_PSI = 0

    else:
        TSR = WT/WA * RAD/rr
        TSR_PSI = (WT_PSI * RAD/rr - TSR*WA_PSI)/WA
        FARG     = 0.5*BLDS*(1.0-rr/RAD)*TSR
        FARG_TSR = 0.5*BLDS*(1.0-rr/RAD)
        FARG = min( FARG , 20.0 )   
        FEXP = euler**(-FARG)
        FEXP_TSR = -FEXP*FARG_TSR   
        F = (2.0/PI) * math.acos(FEXP)
        F_TSR = -(2.0/PI) / ((1.0 - FEXP**2)**0.5) * FEXP_TSR  
        F_PSI = F_TSR*TSR_PSI
        ADW     =  1.0    /TSR
        ADW_PSI = -TSR_PSI/TSR**2

    VA     = WA     - UA
    VA_PSI = WA_PSI
    VT     = UT     - WT
    VT_PSI =        - WT_PSI

    QBI = 4.0/BLDS
    PIR = ((PI*rr)**2 + (QBI*RAD*ADW)**2)**0.5
    PIR_ADW = (QBI*RAD)**2*ADW/PIR
    PIR_PSI = PIR_ADW*ADW_PSI

    GAM     = QBI* F*VT                *PIR
    GAM_PSI = QBI*(F*VT_PSI + F_PSI*VT)*PIR + QBI*F*VT*PIR_PSI

    WSQ = WA**2 + WT**2
    W = (WSQ)**0.5
    W_PSI = (WA*WA_PSI + WT*WT_PSI)/W

    A = BETA - math.atan(WA/WT)
    A_PSI = (-WT*WA_PSI + WA*WT_PSI)/WSQ

    MSQ = WSQ / VSO**2
    MSQ_PSI = 2.0*W*W_PSI / VSO**2

    if (MSQ > MSQMAX):
        MSQ = MSQMAX
        MSQ_PSI = 0.

    PG = 1.0 / (1.0 - MSQ)**0.5
    PG_MSQ = 0.5*PG / (1.0 - MSQ)
    PG_PSI = PG_MSQ*MSQ_PSI

    CL     = (DCLDA*A + CL0)*PG
    CL_PSI =  DCLDA*A_PSI   *PG + (DCLDA*A + CL0)*PG_PSI  

    STALL = False
    if (CL > CLMAX):
        STALL = True
        ACL0 = CL0/DCLDA
        CL   =  CLMAX*math.cos(A-ACL0)
        CL_A = -CLMAX*math.sin(A-ACL0)
        CL_PSI = CL_A*A_PSI

    elif (CL < CLMIN):
        STALL = True
        ACL0 = CL0/DCLDA
        CL   =  CLMIN*math.cos(A-ACL0)
        CL_A = -CLMIN*math.sin(A-ACL0)
        CL_PSI = CL_A*A_PSI

    RES     = GAM     - 0.5*CHORD* CL*W
    RES_PSI = GAM_PSI - 0.5*CHORD*(CL*W_PSI + CL_PSI*W)

    DPSI = -RES/RES_PSI
    DPSI = max(-0.1, min(0.1, DPSI))

    if(abs(DPSI) < EPS):
        LCONV = True
        break

    PSI = PSI + DPSI

print(f"Res, a, CL = {RES} {A} {CL}")