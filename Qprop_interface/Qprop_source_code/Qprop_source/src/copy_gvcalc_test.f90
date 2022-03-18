PROGRAM GVCALC
    
    REAL, parameter :: PI = 3.14159265 
    REAL, parameter :: EPS = 1.0E-6 
    REAL, parameter :: MSQMAX = 0.9    


    !IMPLICIT REAL (A-H,M,O-Z)
    LOGICAL :: STALL, LCONV
    INTEGER :: BLDS
    REAL :: CHORD,BETA,R,RAD,VEL,RPM,OMG,VSO
    REAL :: CL0,DCLDA,CLMIN,CLMAX,MCRIT
    REAL :: GAM
    REAL :: VA
    REAL :: VT
    REAL :: CL
    REAL :: WSQ, MSQ, MSQ_PSI

    BLDS = 4

    CHORD = 0.05
    BETA = 0.523599
    R = 0.1
    RAD = 0.3
    VEL = 10.0
    RPM = 2300.0
    OMG = RPM*2*PI/60
    VSO = 340.0

                          
    CL0 = 0.56 
    DCLDA = 5.58 
    CLmin = -0.3 
    CLmax = 1.37 
    MCRIT = 0.7

    U0A = 0.
    U0T = 0.
    UA     = VEL   + U0A
    UT     = OMG*R - U0T
    WZ = SQRT(UA**2 + UT**2)
    PSI1 = ATAN2(UA,UT)
    PSI2 = BETA + CL0/DCLDA
    PSI = MAX( PSI1 , PSI2 )


    LCONV = .FALSE.
    DO ITER = 1, 20
       COSP = COS(PSI)
       SINP = SIN(PSI)

       WA     = 0.5*UA     + 0.5*WZ    *SINP
       WA_PSI =              0.5*WZ    *COSP

       WT     = 0.5*UT     + 0.5*WZ    *COSP
       WT_PSI =            - 0.5*WZ    *SINP

       WSQ = WA**2 + WT**2
       W = SQRT(WSQ)
       W_PSI = (WA*WA_PSI + WT*WT_PSI)/W

       VA     = WA     - UA
       VA_PSI = WA_PSI

       VT     = UT     - WT
       VT_PSI =        - WT_PSI

       A = BETA - ATAN2(WA,WT)
       A_PSI = (-WT*WA_PSI + WA*WT_PSI)/WSQ

       MSQ = WSQ / VSO**2
       MSQ_PSI = 2.0*W*W_PSI / VSO**2

       IF(MSQ .GT. MSQMAX) THEN
        MSQ = MSQMAX
        MSQ_PSI = 0.
       ENDIF

       PG = 1.0 / SQRT(1.0 - MSQ)
       PG_MSQ = 0.5*PG / (1.0 - MSQ)
       PG_PSI = PG_MSQ*MSQ_PSI

       CL     = (DCLDA*A + CL0)*PG
       CL_PSI =  DCLDA*A_PSI   *PG + (DCLDA*A + CL0)*PG_PSI 
         
       STALL = .FALSE.
       IF    (CL.GT.CLMAX) THEN
        STALL = .TRUE.
        ACL0 = CL0/DCLDA
        CL   =  CLMAX*COS(A-ACL0)
        CL_A = -CLMAX*SIN(A-ACL0)
        CL_PSI = CL_A*A_PSI
       ELSEIF(CL.LT.CLMIN) THEN
        STALL = .TRUE.
        ACL0 = CL0/DCLDA
        CL   =  CLMIN*COS(A-ACL0)
        CL_A = -CLMIN*SIN(A-ACL0)
        CL_PSI = CL_A*A_PSI
       ENDIF

       IF(WA.EQ.0.0) THEN
        F     = 1.0
        F_PSI = 0.

        ADW     = 0.
        ADW_PSI = 0.
       ELSE
        TSR = WT/WA * RAD/R
        TSR_PSI = (WT_PSI * RAD/R - TSR*WA_PSI)/WA

        FARG     = 0.5*BLDS*(1.0-R/RAD)*TSR
        FARG_TSR = 0.5*BLDS*(1.0-R/RAD)
        FARG = MIN( FARG , 20.0 )

        FEXP = EXP(-FARG)
        FEXP_TSR = -FEXP*FARG_TSR

        F = (2.0/PI) * ACOS(FEXP)
        F_TSR = -(2.0/PI) / SQRT(1.0 - FEXP**2) * FEXP_TSR

        F_PSI = F_TSR*TSR_PSI

        ADW     =  1.0    /TSR
        ADW_PSI = -TSR_PSI/TSR**2
       ENDIF

       QBI = 4.0/BLDS
       PIR = SQRT((PI*R)**2 + (QBI*RAD*ADW)**2)
       PIR_ADW = (QBI*RAD)**2*ADW/PIR
       PIR_PSI = PIR_ADW*ADW_PSI

       GAM     = QBI* F*VT                *PIR
       GAM_PSI = QBI*(F*VT_PSI + F_PSI*VT)*PIR + QBI*F*VT*PIR_PSI

       RES     = GAM     - 0.5*CHORD* CL*W
       RES_PSI = GAM_PSI - 0.5*CHORD*(CL*W_PSI + CL_PSI*W)

       DPSI = -RES/RES_PSI
       DPSI = MAX( -0.1 , MIN ( 0.1 , DPSI ) )


       IF(ABS(DPSI) .LT. EPS) THEN
        LCONV = .TRUE.
        GO TO 50
       ENDIF

       PSI = PSI + DPSI
    ENDDO

    WRITE(*,*) 'GVCALC: Not converged.  Res a CL =', RES, A, CL

50  CONTINUE

    print *, "RES, A, CL", RES, A, CL

END program

