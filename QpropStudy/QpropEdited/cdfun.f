SUBROUTINE CDFUN(CL,   RE,   MA, CLCD0,CD0,CD2,REREF,REEXP, MCRIT, CD,CD_CL,CD_RE,CD_MA )
    REAL MA, MCRIT
    DATA CDMF / 10.0 /
    DATA IEXP / 3 /


    FAC = (RE/REREF)**REEXP
    FAC_RE = REEXP*FAC/RE

    CLB = CL - CLCD0
    CD    = (CD0 + CD2*CLB**2 )*FAC
    CD_CL = (      CD2*CLB*2.0)*FAC
    CD_RE = (CD0 + CD2*CLB**2 )*FAC_RE
    CD_MA = 0.

    IF(MA .GT. MCRIT) THEN
        CD    = CD    + CDMF*(MA-MCRIT)** IEXP
        CD_MA = CD_MA + CDMF*(MA-MCRIT)**(IEXP-1) * FLOAT(IEXP)
    ENDIF

    RETURN
END ! CDFUN



SUBROUTINE DCDFUN(BE,WA,WT, CLCD0,CL0,DCLDA,DCD,DCD_BE,DCD_WA,DCD_WT )
    A    = BE - ATAN2(WA,WT)
    A_BE = 1.0
    A_WA = -WT/(WA**2 + WT**2)
    A_WT =  WA/(WA**2 + WT**2)

    ACD0 = (CLCD0-CL0)/DCLDA
    DCD   = 2.0*SIN(A-ACD0)**2
    DCD_A = 4.0*SIN(A-ACD0)*COS(A-ACD0)

    DCD_BE = DCD_A*A_BE
    DCD_WA = DCD_A*A_WA
    DCD_WT = DCD_A*A_WT

    RETURN
END ! DCDFUN
