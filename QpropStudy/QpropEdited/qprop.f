PROGRAM QPROP
 !===============VARIABLES=============================================================================================================================
	IMPLICIT REAL (A-H,M,O-Z)

    PARAMETER (IRDIM=81)
    REAL WORK(IRDIM)
    REAL RB(IRDIM), CB(IRDIM), BB(IRDIM)
    REAL CL0B(IRDIM), DCLDAB(IRDIM), CLMINB(IRDIM), CLMAXB(IRDIM)
    REAL CD0B(IRDIM), CD2UB(IRDIM), CD2LB(IRDIM), CLCD0B(IRDIM)
    REAL REREFB(IRDIM), REEXPB(IRDIM), MCRITB(IRDIM)

    PARAMETER (IDIM=25)
    REAL R(IDIM), C(IDIM), B(IDIM), DR(IDIM)
    REAL CL0(IDIM), DCLDA(IDIM), CLMIN(IDIM), CLMAX(IDIM)
    REAL CD0(IDIM), CD2U(IDIM), CD2L(IDIM), CLCD0(IDIM)
    REAL REREF(IDIM), REEXP(IDIM), MCRIT(IDIM)
    REAL VA(IDIM), VT(IDIM), CL(IDIM), CD(IDIM)
    LOGICAL STALL(IDIM)

    REAL TP_C(IDIM), TP_B(IDIM), QP_C(IDIM), QP_B(IDIM)

    PARAMETER (NMPDIM=10)
    REAL PARMOT(NMPDIM)
    CHARACTER*32 PMLAB(NMPDIM)

    CHARACTER*1 CHARF, ANS
    CHARACTER*80 PNAME, MNAME
    CHARACTER*80 ARGP1, ARGP2, ARGP3, ARGP4, ARGP5, ARGP6, ARGP7, ARGP8, ARGP9, ARGP10
    CHARACTER*80 FILNAM
    CHARACTER*128 LINE

    LOGICAL LRDUMP
    LOGICAL LRPMSET, LVOLTSET, LTHRUSET, LTORQSET, LAMPSSET, LPELESET
    LOGICAL ERROR

    INCLUDE 'QDEF.INC'

    REAL RVAL(15)
    INTEGER IVAL(15)

    DATA PI / 3.14159265 /
    DATA EPS / 1.0E-8 /
    DATA VERSION / 1.22 /

    MCRIT0 = 0.70
 
 !===============SETUP ARGUMENTS=============================================================================================================================

    CALL GETARG0(1,ARGP1)
    CALL GETARG0(2,ARGP2)
    CALL GETARG0(3,ARGP3)
    CALL GETARG0(4,ARGP4)
    CALL GETARG0(5,ARGP5)
    CALL GETARG0(6,ARGP6)
    CALL GETARG0(7,ARGP7)
    CALL GETARG0(8,ARGP8)
    CALL GETARG0(9,ARGP9)
    CALL GETARG0(10,ARGP10)

    IF(ARGP1.EQ.' ') THEN
        WRITE(*,1005)
        1005  FORMAT(
        	& /' QPROP usage:'
        	&//' % qprop propfile motorfile Vel Rpm ',
        	&             '[ Volt dBeta Thrust Torque Amps Pele ]',
        	&  '   (single-point)'
        	&//' % qprop propfile motorfile Vel1,Vel2,dVel Rpm ["]           ',
        	&  '   (multi-point 1-parameter sweep over Vel, Rpm set)'
        	&//' % qprop propfile motorfile Vel1,Vel2,dVel 0 Volt ["]        ',
        	&  '   (multi-point 1-parameter sweep over Vel, Volt set)'
        	&//' % qprop propfile motorfile Vel1,Vel2,dVel Rpm1,Rpm2,dRpm ["]',
        	&  '   (multi-point 2-parameter sweep over Vel and Rpm)'
        	&//' % qprop propfile motorfile runfile                          ',
        	&  '   (multi-point, via file specification)'
        	&       )
        WRITE(*,*)
        WRITE(*,*) 'Run with default inputs?  Y'
        READ(*,1000) ANS
        IF(INDEX('Nn',ANS) .NE. 0) STOP
        WRITE(*,*) 
    ENDIF

    RHO = RHO1   ! density
    RMU = RMU1   ! viscosity
    VSO = VSO1   ! speed of sound

    CALL QCGET(RHO,RMU,VSO)

    PNAME = 'Graupner CAM 6x3 folder'
    BLDS = 2.0

    NR = 7

    DO IR = 1, NR
        CL0B(IR) = 0.5
        DCLDAB(IR) = 5.8
        CLMINB(IR) = -0.4
        CLMAXB(IR) = 1.2
	ENDDO

    DO IR = 1, NR
        CD0B(IR) = 0.028
        CD2UB(IR) = 0.050
        CD2LB(IR) = 0.050
        CLCD0B(IR) = 0.5
        REREFB(IR) = 70000.0
        REEXPB(IR) = -0.7
        MCRITB(IR) = MCRIT0
    ENDDO

    RFAC = 0.0254
    RADD = 0.
    RB(1) = 0.75  
    RB(2) = 1.00  
    RB(3) = 1.50  
    RB(4) = 2.00  
    RB(5) = 2.50  
    RB(6) = 2.875 
    RB(7) = 3.00  

    CFAC = 0.0254
    CADD = 0.
    CB(1) = 0.66 
    CB(2) = 0.69 
    CB(3) = 0.63 
    CB(4) = 0.55 
    CB(5) = 0.44 
    CB(6) = 0.30 
    CB(7) = 0.19 

    BFAC = 1.0
    BADD = 0.
    BB(1) = 27.5
    BB(2) = 22.0
    BB(3) = 15.2
    BB(4) = 10.2
    BB(5) =  6.5
    BB(6) =  4.6
    BB(7) =  4.2

    RAD = RB(NR)

    MNAME = "Speed-400 3321 (6V) direct drive"
    IMOTYPE = 1
    PARMOT(1) = 0.31    ! Rmotor  (Ohms)
    PARMOT(2) = 0.77    ! Io      (Amps)
    PARMOT(3) = 2760.0  ! Kv      (rpm/Volt)
    PMLAB(1) = 'R  (Ohm)'
    PMLAB(2) = 'Io (Amp)'
    PMLAB(3) = 'Kv (rpm/Volt)'
    NMPAR = 3

    VEL1 =  0.0
    VEL2 = 10.0
    NVEL = 6

    RPM1 = 10000.0
    RPM2 = 16000.0
    NRPM = 7

    VOLT1 = 6.0
    VOLT2 = 9.0
    NVOLT = 4

    DBET1 = -2.0
    DBET2 =  2.0
    NDBET = 5

    THRU1 = 0.
    THRU2 = 0.
    NTHRU = 1

    AMPS1 = 0.
    AMPS2 = 0.
    NAMPS = 1

    PELE1 = 0.
    PELE2 = 0.
    NPELE = 1

    LRDUMP = .FALSE.

    1000 FORMAT(A)

    FILNAM = ARGP1
    IF(FILNAM.EQ.' ') GO TO 18

    LU = 1
    OPEN(LU,FILE=FILNAM,STATUS='OLD',ERR=18)

    ILINE = 0


    CALL FREAD(LU,LINE,ILINE,IERR,PNAME)
    IF(IERR.EQ.+1) GO TO 900
    IF(IERR.EQ.-1) GO TO 950


    NVAL = 2
    CALL RREAD(LU,LINE,ILINE,IERR,NVAL,RVAL)
    IF(IERR.EQ.+1) GO TO 900
    IF(IERR.EQ.-1) GO TO 950
    IF(NVAL.LT. 1) GO TO 980
    BLDS = RVAL(1)
    IF(NVAL.GE.2) THEN
        RAD = RVAL(2)
    ELSE
        RAD = 0.
    ENDIF

    NVAL = 2
    CALL RREAD(LU,LINE,ILINE,IERR,NVAL,RVAL)
    IF(IERR.EQ.+1) GO TO 900
    IF(IERR.EQ.-1) GO TO 950
    IF(NVAL.LT. 2) GO TO 980
    DO IR = 1, IRDIM
        CL0B(IR)   = RVAL(1)
        DCLDAB(IR) = RVAL(2)
    ENDDO

    NVAL = 2
    CALL RREAD(LU,LINE,ILINE,IERR,NVAL,RVAL)
    IF(IERR.EQ.+1) GO TO 900
    IF(IERR.EQ.-1) GO TO 950
    IF(NVAL.LT. 2) GO TO 980
    DO IR = 1, IRDIM
        CLMINB(IR) = RVAL(1)
        CLMAXB(IR) = RVAL(2)
    ENDDO

    NVAL = 4
    CALL RREAD(LU,LINE,ILINE,IERR,NVAL,RVAL)
    IF(IERR.EQ.+1) GO TO 900
    IF(IERR.EQ.-1) GO TO 950
    IF(NVAL.LT. 3) GO TO 980
    DO IR = 1, IRDIM
        CD0B(IR)   = RVAL(1)
        CD2UB(IR)  = RVAL(2)
        CD2LB(IR)  = RVAL(3)
        CLCD0B(IR) = RVAL(4)
    ENDDO

    NVAL = 2
    CALL RREAD(LU,LINE,ILINE,IERR,NVAL,RVAL)
    IF(IERR.EQ.+1) GO TO 900
    IF(IERR.EQ.-1) GO TO 950
    IF(NVAL.LT. 2) GO TO 980
    DO IR = 1, IRDIM
      	REREFB(IR) = RVAL(1)
      	REEXPB(IR) = RVAL(2)
    ENDDO


    NVAL = 3
    CALL RREAD(LU,LINE,ILINE,IERR,NVAL,RVAL)
    IF(IERR.EQ.+1) GO TO 900
    IF(IERR.EQ.-1) GO TO 950
    IF(NVAL.LT. 3) GO TO 980
    RFAC = RVAL(1)
    CFAC = RVAL(2)
    BFAC = RVAL(3)

    NVAL = 3
    CALL RREAD(LU,LINE,ILINE,IERR,NVAL,RVAL)
    IF(IERR.EQ.+1) GO TO 900
    IF(IERR.EQ.-1) GO TO 950
    IF(NVAL.LT. 3) GO TO 980
    RADD = RVAL(1)
    CADD = RVAL(2)
    BADD = RVAL(3)

    KR = 0

    14   CONTINUE

    NVAL = 13
    CALL RREAD(LU,LINE,ILINE,IERR,NVAL,RVAL)
    IF(IERR.EQ.+1) GO TO 900
    IF(IERR.EQ.-1) GO TO 16
    IF(NVAL.LT. 3) GO TO 980

    KR = KR + 1
    IR = MIN( KR , IRDIM )
    RB(IR) = RVAL(1)
    CB(IR) = RVAL(2)
    BB(IR) = RVAL(3)

    MCRITB(IR) = MCRIT0

    IF(NVAL.GE. 4) CL0B(IR)   = RVAL( 4)
    IF(NVAL.GE. 5) DCLDAB(IR) = RVAL( 5)
    IF(NVAL.GE. 6) CLMINB(IR) = RVAL( 6)
    IF(NVAL.GE. 7) CLMAXB(IR) = RVAL( 7)
    IF(NVAL.GE. 8) CD0B(IR)   = RVAL( 8)
    IF(NVAL.GE. 9) CD2UB(IR)  = RVAL( 9)
    IF(NVAL.GE.10) CD2LB(IR)  = RVAL(10)
    IF(NVAL.GE.11) CLCD0B(IR) = RVAL(11)
    IF(NVAL.GE.12) REREFB(IR) = RVAL(12)
    IF(NVAL.GE.13) REEXPB(IR) = RVAL(13)

	GO TO 14

 	16   CONTINUE
    CLOSE(LU)

    IF(KR.GT.0) THEN
     	NR = IR
     	IF(KR.GT.NR) THEN
      		WRITE(*,*) 'Array overflow.  Increase IRDIM to', KR
      		STOP
     	ENDIF
    ENDIF
    GO TO 19

 	18   CONTINUE
    WRITE(*,*)
    WRITE(*,*) 'Prop file not found:  ', FILNAM(1:48)
    WRITE(*,*) 'Default prop used  :  ', PNAME


 	19   CONTINUE

    IF(NR.LE.1) THEN
    	WRITE(*,*)
       	WRITE(*,*) '*** Must define at least two radial stations'
       	STOP
    ENDIF


    DO IR = 1, NR
    	RB(IR) =  RB(IR)*RFAC + RADD
        CB(IR) =  CB(IR)*CFAC + CADD
        BB(IR) = (BB(IR)*BFAC + BADD)* PI / 180.0
    ENDDO

    IF(RAD .EQ. 0.0) THEN
    	RAD = RB(NR)
	ENDIF

    DO IR = 1, NR-1
    	IF(CB(IR) .LE. 0.0) STOP 'Chords must be positive'
        IF(RB(IR) .LT. 0.0) STOP 'Radii must be nonnegative'
        IF(RB(IR) .GE. RB(IR+1)) STOP 'Radii must increase monotonically'
    ENDDO

    IF(RAD .LT. RB(NR)) THEN
       	WRITE(*,1050) RAD, RB(NR)
 		1050  FORMAT(/' Given on line 2:  R =', G12.4,
     		&        /' Last r station :  r =', G12.4,
     		&       //' Must have  R > r' / )
       	STOP
    ENDIF

    FILNAM = ARGP2
    IF(FILNAM.EQ.' ') GO TO 28

    LU = 2
    OPEN(LU,FILE=FILNAM,STATUS='OLD',ERR=28)


    DO IMPAR = 1, NMPDIM
      	PARMOT(IMPAR) = 0.0
      	PMLAB(IMPAR) = ' '
    ENDDO

    ILINE = 0


    CALL FREAD(LU,LINE,ILINE,IERR,MNAME)
    IF(IERR.EQ.+1) GO TO 900
    IF(IERR.EQ.-1) GO TO 950


    NVAL = 1
    CALL IREAD(LU,LINE,ILINE,IERR,NVAL,IVAL)
    IF(IERR.EQ.+1) GO TO 900
    IF(IERR.EQ.-1) GO TO 950
    IF(NVAL.LT. 1) GO TO 980
    IMOTYPE = IVAL(1)


    DO IMPAR = 1, NMPDIM+1
      	NVAL = 1
      	CALL RREAD(LU,LINE,ILINE,IERR,NVAL,RVAL)
      	IF(IERR.EQ.+1) GO TO 900
      	IF(IERR.EQ.-1) GO TO 25
      	IF(NVAL.LT. 1) GO TO 980
      	IF(IMPAR.EQ.NMPDIM+1) THEN
       		WRITE(*,*) '* Motor parameter array overflow. Increase NMPDIM'
       		STOP
      	ENDIF
      	PARMOT(IMPAR) = RVAL(1)
      	KEX = INDEX(LINE,'!')
      	IF(KEX.GE.1) THEN
       		PMLAB(IMPAR) = LINE(KEX+1:80)
      	ELSE
       		PMLAB(IMPAR) = ' '
      	ENDIF
    ENDDO

 	25   CONTINUE
    NMPAR = IMPAR-1

    CLOSE(LU)
    GO TO 29

 	28   CONTINUE
	WRITE(*,*)
	WRITE(*,*) 'Motor file not found:  ', FILNAM(1:48)
	WRITE(*,*) 'Default motor used  :  ', MNAME

 	29   CONTINUE

    FILNAM = ARGP3
    IF(FILNAM.EQ.' ') GO TO 80

    LU = 4
    OPEN(LU,FILE=FILNAM,STATUS='OLD',ERR=31)

    ILINE = 0

    NVAL = 3
    CALL RREAD(LU,LINE,ILINE,IERR,NVAL,RVAL)
    IF(IERR.EQ.+1) GO TO 900
    IF(IERR.EQ.-1) GO TO 950
    IF(NVAL.LT. 3) GO TO 980
    VEL1 = RVAL(1)
    VEL2 = RVAL(2)
    NVEL = INT( RVAL(3) + 0.01 )

    NVAL = 3
    CALL RREAD(LU,LINE,ILINE,IERR,NVAL,RVAL)
    IF(IERR.EQ.+1) GO TO 900
    IF(IERR.EQ.-1) GO TO 950
    IF(NVAL.LT. 3) GO TO 980
    RPM1 = RVAL(1)
    RPM2 = RVAL(2)
    NRPM = INT( RVAL(3) + 0.01 )

    NVAL = 3
    CALL RREAD(LU,LINE,ILINE,IERR,NVAL,RVAL)
    IF(IERR.EQ.+1) GO TO 900
    IF(IERR.EQ.-1) GO TO 950
    IF(NVAL.LT. 3) GO TO 980
    VOLT1 = RVAL(1)
    VOLT2 = RVAL(2)
    NVOLT = INT( RVAL(3) + 0.01 )

    NVAL = 3
    CALL RREAD(LU,LINE,ILINE,IERR,NVAL,RVAL)
    IF(IERR.EQ.+1) GO TO 900
    IF(IERR.EQ.-1 .OR. NVAL.LT. 3) THEN
     	DBET1 = 0.0
     	DBET2 = 0.0
     	NDBET = 0
    ELSE
     	DBET1 = RVAL(1)
     	DBET2 = RVAL(2)
     	NDBET = INT( RVAL(3) + 0.01 )
    ENDIF
    NDBET = MAX( 1 , NDBET )

    NVAL = 3
    CALL RREAD(LU,LINE,ILINE,IERR,NVAL,RVAL)
    IF(IERR.EQ.+1) GO TO 900
    IF(IERR.EQ.-1 .OR. NVAL.LT. 3) THEN
     	THRU1 = 0.0
     	THRU2 = 0.0
     	NTHRU = 0
    ELSE
     	THRU1 = RVAL(1)
     	THRU2 = RVAL(2)
     	NTHRU = INT( RVAL(3) + 0.01 )
    ENDIF
    NTHRU = MAX( 1 , NTHRU )

    CLOSE(LU)
    GO TO 82

 	31   CONTINUE

    CALL PPARSE(ARGP3,VEL1,VEL2,NVEL,IERR)
    IF(IERR.EQ.+1) GO TO 80
    IF(IERR.EQ.-1) GO TO 80


    RPM1 = 0.
    RPM2 = 0.
    NRPM = 0

    VOLT1 = 0.
    VOLT2 = 0.
    NVOLT = 0

    DBET1 = 0.
    DBET2 = 0.
    NDBET = 1

    THRU1 = 0.
    THRU2 = 0.
    NTHRU = 1

    TORQ1 = 0.
    TORQ2 = 0.
    NTORQ = 1

    AMPS1 = 0.
    AMPS2 = 0.
    NAMPS = 1

    PELE1 = 0.
    PELE2 = 0.
    NPELE = 1


    CALL PPARSE(ARGP4,RPM1,RPM2,NRPM,IERR)
    IF(IERR.EQ.+1 .OR. IERR.EQ.-1) THEN
    	RPM = 0.0
    	NRPM = 0
    ENDIF


    CALL PPARSE(ARGP5,VOLT1,VOLT2,NVOLT,IERR)
    IF(IERR.EQ.+1 .OR. IERR.EQ.-1) THEN
    	VOLT = 0.0
       	NVOLT = 0
	ENDIF


    CALL PPARSE(ARGP6,DBET1,DBET2,NDBET,IERR)
    IF(IERR.EQ.+1 .OR. IERR.EQ.-1) THEN
       	DBET = 0.0
       	NDBET = 1
    ENDIF


	CALL PPARSE(ARGP7,THRU1,THRU2,NTHRU,IERR)
    IF(IERR.EQ.+1 .OR. IERR.EQ.-1) THEN
       	THRU = 0.0
       	NTHRU = 1
    ENDIF


    CALL PPARSE(ARGP8,TORQ1,TORQ2,NTORQ,IERR)
    IF(IERR.EQ.+1 .OR. IERR.EQ.-1) THEN
    	TORQ = 0.0
       	NTHRU = 1
    ENDIF


	CALL PPARSE(ARGP9,AMPS1,AMPS2,NAMPS,IERR)
    IF(IERR.EQ.+1 .OR. IERR.EQ.-1) THEN
       	AMPS = 0.0
       	NAMPS = 1
    ENDIF


    CALL PPARSE(ARGP10,PELE1,PELE2,NPELE,IERR)
    IF(IERR.EQ.+1 .OR. IERR.EQ.-1) THEN
       	PELE = 0.0
       	NPELE = 1
    ENDIF


    LRDUMP = NVEL .LE.1 .AND. NRPM .LE.1 .AND. NVOLT.LE.1 .AND. NDBET.LE.1 .AND. NTHRU.LE.1 .AND. NTORQ.LE.1 .AND. NAMPS.LE.1 .AND. NPELE.LE.1

    GO TO 82


 	80   CONTINUE
    WRITE(*,*)
    WRITE(*,*) 'Run parameter file not found: ', FILNAM(1:48)
    WRITE(*,*) 'Default velocities, voltages, pitch used'

 	82   CONTINUE

    IF(NRPM .EQ.0 .AND. NVOLT.EQ.0 .AND. NTHRU.EQ.0 .AND.NAMPS.EQ.0 .AND.NPELE.EQ.0) THEN
       	WRITE(*,*)  'Must specify either Rpm or Volts or Thrust or Amps or Pele'
       	STOP
    ENDIF

    LRPMSET  = NRPM  .GT. 0 .AND. (RPM1  .NE. 0.0 .OR. RPM2  .NE. 0.0)
    LVOLTSET = NVOLT .GT. 0 .AND. (VOLT1 .NE. 0.0 .OR. VOLT2 .NE. 0.0)
    LTHRUSET = NTHRU .GT. 0 .AND. (THRU1 .NE. 0.0 .OR. THRU2 .NE. 0.0)
    LTORQSET = NTORQ .GT. 0 .AND. (TORQ1 .NE. 0.0 .OR. TORQ2 .NE. 0.0)
    LAMPSSET = NAMPS .GT. 0 .AND. (AMPS1 .NE. 0.0 .OR. AMPS2 .NE. 0.0)
    LPELESET = NPELE .GT. 0 .AND. (PELE1 .NE. 0.0 .OR. PELE2 .NE. 0.0)

    R0 = RB(1)
    R1 = RB(NR)

    N = IDIM
    DO I = 1, N
    	FRAC = (FLOAT(I)-0.5)/FLOAT(N)
    	R(I) = R0*(1.0-FRAC) + R1*FRAC
    	DR(I) = (R1-R0)/FLOAT(N)
    ENDDO

    CALL SPLINE(CB,WORK,RB,NR)
    DO I = 1, N
      	C(I) = SEVAL(R(I),CB,WORK,RB,NR)
    ENDDO

    CALL SPLINE(BB,WORK,RB,NR)
    DO I = 1, N
      	B(I) = SEVAL(R(I),BB,WORK,RB,NR)
    ENDDO

    CALL SPLINE(CL0B,WORK,RB,NR)
    DO I = 1, N
      	CL0(I) = SEVAL(R(I),CL0B,WORK,RB,NR)
    ENDDO

    CALL SPLINE(DCLDAB,WORK,RB,NR)
    DO I = 1, N
      	DCLDA(I) = SEVAL(R(I),DCLDAB,WORK,RB,NR)
    ENDDO

    CALL SPLINE(CLMINB,WORK,RB,NR)
    DO I = 1, N
      	CLMIN(I) = SEVAL(R(I),CLMINB,WORK,RB,NR)
    ENDDO

    CALL SPLINE(CLMAXB,WORK,RB,NR)
    DO I = 1, N
      	CLMAX(I) = SEVAL(R(I),CLMAXB,WORK,RB,NR)
    ENDDO

    CALL SPLINE(CD0B,WORK,RB,NR)
    DO I = 1, N
      	CD0(I) = SEVAL(R(I),CD0B,WORK,RB,NR)
    ENDDO

    CALL SPLINE(CD2UB,WORK,RB,NR)
    DO I = 1, N 
      	CD2U(I) = SEVAL(R(I),CD2UB,WORK,RB,NR)
    ENDDO

    CALL SPLINE(CD2LB,WORK,RB,NR)
    DO I = 1, N
      	CD2L(I) = SEVAL(R(I),CD2LB,WORK,RB,NR)
    ENDDO

    CALL SPLINE(CLCD0B,WORK,RB,NR)
    DO I = 1, N
      	CLCD0(I) = SEVAL(R(I),CLCD0B,WORK,RB,NR)
    ENDDO

    CALL SPLINE(REREFB,WORK,RB,NR)
    DO I = 1, N
      	REREF(I) = SEVAL(R(I),REREFB,WORK,RB,NR)
    ENDDO

    CALL SPLINE(REEXPB,WORK,RB,NR)
    DO I = 1, N
      	REEXP(I) = SEVAL(R(I),REEXPB,WORK,RB,NR)
    ENDDO

    CALL SPLINE(MCRITB,WORK,RB,NR)
    DO I = 1, N
      	MCRIT(I) = SEVAL(R(I),MCRITB,WORK,RB,NR)
    ENDDO


    ERROR = .FALSE.
    DO I = 1, N
      	IF(C(I) .LE. 0.0) THEN
       		WRITE(*,*) 'Negative chord at i =', I
       		ERROR = .TRUE.
      	ENDIF
      	IF(REREF(I) .LE. 0.0) THEN
       		WRITE(*,*) 'Negative Re_ref at i =', I
       		ERROR = .TRUE.
      	ENDIF
      	IF(MCRIT(I) .LE. 0.0) THEN
       		WRITE(*,*) 'Negative Mcrit at i =', I
       		ERROR = .TRUE.
      	ENDIF
    ENDDO

    IF(ERROR) THEN
        WRITE(*,*)
       	WRITE(*,1100) ' i   radius   chord     beta    Re_ref'
        DO I = 1, N
          	IRE = INT( REREF(I) )
          	WRITE(*,1070) I, R(I), C(I), B(I)*180.0/PI, IRE
 			1070     FORMAT(1X,I3, F9.4, F9.4, F9.3, I9)
        ENDDO
        WRITE(*,*)
        STOP
    ENDIF

    LU = 6
    WRITE(*,*)

 	1105 FORMAT('# QPROP Version', F5.2)
 	1100 FORMAT('# ', A,A,A,A)
 	1110 FORMAT('#  ', G12.5, 1X, A)
 	1120 FORMAT('#   rho =', G12.5,' kg/m^3'
     &      /'#   mu  =', G12.5,' kg/m-s'
     &      /'#   a   =', G12.5,' m/s   ' )

    WRITE(LU,1105) VERSION
    WRITE(LU,1100)
    WRITE(LU,1100) PNAME
    WRITE(LU,1100)
    WRITE(LU,1100) MNAME
    DO IMPAR=1, NMPAR
      	WRITE(LU,1110) PARMOT(IMPAR), PMLAB(IMPAR)
    ENDDO
    WRITE(LU,1100)
    WRITE(LU,1120) RHO, RMU, VSO
    WRITE(LU,1100)
    WRITE(LU,1100)' 1         2        3          4          5        ',' 6            7         8       9        10        11    ','    12          13        14        15      16          17   ','        18      19'
    WRITE(LU,1100)
    WRITE(LU,1100)' V(m/s)    rpm      Dbeta      T(N)       Q(N-m)   ',' Pshaft(W)    Volts     Amps    effmot   effprop   adv   ','    CT          CP        DV(m/s)   eff     Pelec       Pprop','        cl_avg  cd_avg'

    IF(LRDUMP) THEN
       	CHARF = '#'
    ELSE
       	CHARF = ' '
    ENDIF

    NVELM = MAX( 1 , NVEL-1 )
    DVEL = (VEL2-VEL1)/FLOAT(NVELM)

    NDBETM = MAX( 1 , NDBET-1 )
    DDBET = (DBET2-DBET1)/FLOAT(NDBETM)

    IF(LRPMSET) THEN
    	NRPMM = MAX( 1 , NRPM-1 )
    	PAR1 = RPM1
    	PAR2 = RPM2
    	DPAR = (RPM2-RPM1)/FLOAT(NRPMM)
    	NPAR = NRPM
    ELSEIF(LVOLTSET) THEN
    	NVOLTM = MAX( 1 , NVOLT-1 )
    	PAR1 = VOLT1
    	PAR2 = VOLT2
    	DPAR = (VOLT2-VOLT1)/FLOAT(NVOLTM)
    	NPAR = NVOLT
    ELSEIF(LTHRUSET) THEN
    	NTHRUM = MAX( 1 , NTHRU-1 )
    	PAR1 = THRU1
    	PAR2 = THRU2
    	DPAR = (THRU2-THRU1)/FLOAT(NTHRUM)
    	NPAR = NTHRU
    ELSEIF(LTORQSET) THEN
    	NTORQM = MAX( 1 , NTORQ-1 )
    	PAR1 = TORQ1
    	PAR2 = TORQ2
    	DPAR = (TORQ2-TORQ1)/FLOAT(NTORQM)
    	NPAR = NTORQ
    ELSEIF(LAMPSSET) THEN
    	NAMPSM = MAX( 1 , NAMPS-1 )
    	PAR1 = AMPS1
    	PAR2 = AMPS2
    	DPAR = (AMPS2-AMPS1)/FLOAT(NAMPSM)
    	NPAR = NAMPS
    ELSEIF(LPELESET) THEN
    	NPELEM = MAX( 1 , NPELE-1 )
    	PAR1 = PELE1
    	PAR2 = PELE2
    	DPAR = (PELE2-PELE1)/FLOAT(NPELEM)
    	NPAR = NPELE
    ELSE
    	WRITE(*,*) 'Additional parameter not specified'
    	WRITE(*,*) ' RPM   :',  lrpmset
    	WRITE(*,*) ' Volt  :',  lvoltset
    	WRITE(*,*) ' Thrust:',  lthruset
    	WRITE(*,*) ' Torque:',  ltorqset
    	WRITE(*,*) ' Amps  :',  lampsset
    	WRITE(*,*) ' Pelec :',  lpeleset
    	STOP
    ENDIF

    LRDUMP = NVEL .LE.1 .AND. NRPM .LE.1 .AND. NVOLT.LE.1 .AND. NDBET.LE.1 .AND. NTHRU.LE.1 .AND. NTORQ.LE.1 .AND. NAMPS.LE.1 .AND. NPELE.LE.1

 !===============CALCULATE=============================================================================================================================
    DO IDBET = 1, NDBET
        DBET = DBET1 + DDBET*FLOAT(IDBET-1)
        DBE = DBET * PI/180.0
        DO IPAR = 1, NPAR
          	PAR = PAR1 + DPAR*FLOAT(IPAR-1)
          	IF    (LRPMSET ) THEN
           		RPM = PAR
          	ELSEIF(LVOLTSET) THEN
           		VOLT = PAR
          	ELSEIF(LTHRUSET) THEN
           		THRU = PAR
          	ELSEIF(LTORQSET) THEN
           		TORQ = PAR
          	ELSEIF(LAMPSSET ) THEN
           		AMPS = PAR
          	ELSEIF(LPELESET) THEN
           		PELE = PAR
          	ENDIF

          	DO IVEL = 1, NVEL
            	VEL = VEL1 + DVEL*FLOAT(IVEL-1)
            	IF    (LRPMSET) THEN
             		OMG = RPM * PI/30.0
            	ELSEIF(LVOLTSET) THEN
             		I = MAX( 1 , (8*N)/10 )
             		RT = R(I)
             		BT = B(I) - CL0(I)/DCLDA(I) + DBE
             		BT = MAX( 0.02 , MIN( 0.45*PI , BT ) )
             		IF(VEL.EQ.0.0) THEN
              			OMG = 1.0
             		ELSE
              			OMG = VEL/(RT*TAN(BT))
             		ENDIF
            	ELSE
             		I = MAX( 1 , (8*N)/10 )
             		RT = R(I)
             		BT = B(I) - CL0(I)/DCLDA(I) + DBE
             		BT = MAX( 0.02 , MIN( 0.45*PI , BT ) )
             		IF(VEL.EQ.0.0) THEN
              			OMG = 1.0
             		ELSE
              			OMG = VEL/(RT*TAN(BT))
             		ENDIF
             		QP = 0.
             		CALL VOLTM(OMG,QP, IMOTYPE, PARMOT,NMPAR, VM,VM_OMG,VM_QP, AM,AM_OMG,AM_QP )
             		VOLT = VM
            	ENDIF

            	DO ITER = 1, 25
              		CALL TQCALC(N,C,B,R,DR,VA,VT,CL,CD,STALL,BLDS,RAD,VEL,OMG,DBE,RHO,RMU,VSO,CL0,DCLDA,CLMIN,CLMAX,MCRIT,CD0,CD2U,CD2L,CLCD0,REREF,REEXP,TP, TP_VEL, TP_OMG, TP_DBE, TP_C, TP_B,QP, QP_VEL, QP_OMG, QP_DBE, QP_C, QP_B )
              		IF(LRPMSET) THEN
               			RES = 0.
               			RES_OMG = 1.0

               			DOMG = -RES/RES_OMG
               			DVOLT = 0.

               			CALL VOLTM(OMG,QP, IMOTYPE, PARMOT,NMPAR,VM,VM_OMG,VM_QP,AM,AM_OMG,AM_QP )
               			VOLT = VM
               			AMPS = AM
              		ELSEIF(LVOLTSET) THEN
               			CALL MOTORQ(OMG,VOLT, IMOTYPE, PARMOT,NMPAR,QM,QM_OMG,QM_VOLT,AM,AM_OMG,AM_VOLT )
               			RES     = QP     - QM
               			RES_OMG = QP_OMG - QM_OMG

               			DOMG = -RES/RES_OMG
               			DVOLT = 0.
               			AMPS = AM
              		ELSEIF(LTHRUSET) THEN
               			CALL MOTORQ(OMG,VOLT, IMOTYPE, PARMOT,NMPAR,QM,QM_OMG,QM_VOLT,AM,AM_OMG,AM_VOLT )
               			RES1      = QP     - QM
               			RES1_OMG  = QP_OMG - QM_OMG
               			RES1_VOLT =        - QM_VOLT

               			RES2      = TP     - THRU
               			RES2_OMG  = TP_OMG
               			RES2_VOLT = 0.

               			A11 = RES1_OMG
               			A12 = RES1_VOLT
               			A21 = RES2_OMG
               			A22 = RES2_VOLT

               			DET   =   A11 *A22  - A12 *A21
               			DOMG  = -(RES1*A22  - A12 *RES2) / DET
               			DVOLT = -(A11 *RES2 - RES1*A21 ) / DET

               			AMPS = AM

              		ELSEIF(LTORQSET) THEN
               			CALL MOTORQ(OMG,VOLT, IMOTYPE, PARMOT,NMPAR, QM,QM_OMG,QM_VOLT,AM,AM_OMG,AM_VOLT )
            			RES1      = QP     - QM
            			RES1_OMG  = QP_OMG - QM_OMG
            			RES1_VOLT =        - QM_VOLT

            			RES2      = QP     - TORQ
            			RES2_OMG  = QP_OMG
            			RES2_VOLT = 0.

            			A11 = RES1_OMG
            			A12 = RES1_VOLT
            			A21 = RES2_OMG
            			A22 = RES2_VOLT

            			DET   =   A11 *A22  - A12 *A21
            			DOMG  = -(RES1*A22  - A12 *RES2) / DET
            			DVOLT = -(A11 *RES2 - RES1*A21 ) / DET

            			AMPS = AM

              		ELSEIF(LAMPSSET) THEN
               			CALL MOTORQ(OMG,VOLT, IMOTYPE, PARMOT,NMPAR,QM,QM_OMG,QM_VOLT,AM,AM_OMG,AM_VOLT )
               			RES1      = QP     - QM
               			RES1_OMG  = QP_OMG - QM_OMG
               			RES1_VOLT =        - QM_VOLT

               			RES2      = AM   - AMPS
               			RES2_OMG  = AM_OMG
               			RES2_VOLT = AM_VOLT

               			A11 = RES1_OMG
               			A12 = RES1_VOLT
               			A21 = RES2_OMG
               			A22 = RES2_VOLT

               			DET   =   A11 *A22  - A12 *A21
               			DOMG  = -(RES1*A22  - A12 *RES2) / DET
               			DVOLT = -(A11 *RES2 - RES1*A21 ) / DET

              		ELSEIF(LPELESET) THEN
               			CALL MOTORQ(OMG,VOLT, IMOTYPE, PARMOT,NMPAR,QM,QM_OMG,QM_VOLT,AM,AM_OMG,AM_VOLT )
            			RES1      = QP     - QM
            			RES1_OMG  = QP_OMG - QM_OMG
            			RES1_VOLT =        - QM_VOLT


            			RES2      = VOLT*AM   - PELE
            			RES2_OMG  = VOLT*AM_OMG
            			RES2_VOLT = VOLT*AM_VOLT + AM

            			A11 = RES1_OMG
            			A12 = RES1_VOLT
            			A21 = RES2_OMG
            			A22 = RES2_VOLT

            			DET   =   A11 *A22  - A12 *A21
            			DOMG  = -(RES1*A22  - A12 *RES2) / DET
            			DVOLT = -(A11 *RES2 - RES1*A21 ) / DET

            			AMPS = AM

              		ENDIF

              		RLX = 1.0
              		IF(RLX*DOMG  .GT.  1.0*OMG) RLX =  1.0*OMG/DOMG
              		IF(RLX*DOMG  .LT. -0.5*OMG) RLX = -0.5*OMG/DOMG

              		IF(RLX*DVOLT .GT.  2.0*VOLT) RLX =  2.0*VOLT/DVOLT
              		IF(RLX*DVOLT .LT. -0.5*VOLT) RLX = -0.5*VOLT/DVOLT

              		IF(ABS(DOMG) .LT. EPS*ABS(OMG)) GO TO 110

              		OMG  = OMG  + RLX*DOMG
              		VOLT = VOLT + RLX*DVOLT

 				ENDDO
 				110        CONTINUE

        		DTSUM = 0.
        		CLAVG = 0.
        		CDAVG = 0.

        		DO I = 1, N
        		  	WA = VEL + VA(I)
        		  	WT = OMG*R(I) - VT(I)
        		  	WSQ = WA**2 + WT**2
        		  	DTSUM = DTSUM + WSQ*C(I)*DR(I)
        		  	CLAVG = CLAVG + WSQ*C(I)*DR(I)*CL(I)
        		  	CDAVG = CDAVG + WSQ*C(I)*DR(I)*CD(I)
        		ENDDO
        		CLAVG = CLAVG / DTSUM
        		CDAVG = CDAVG / DTSUM

        		RPM = OMG*30.0/PI
        		PPROP = TP*VEL
        		POWER = QP*OMG

        		PINPUT = VOLT*AMPS

        		IF(POWER .NE. 0.0) THEN
         			EFFP = PPROP/POWER
        		ELSE
         			EFFP = 0.
        		ENDIF

        		IF(PINPUT .NE. 0.0) THEN
        		 	EFFM = POWER/PINPUT
        		ELSE
        		 	EFFM = 0.0
        		ENDIF

        		EFF = EFFM*EFFP

        		IF(ABS(OMG).GT.0.0) THEN
        		 	ADV = VEL/(OMG*RAD)
        		ELSE
        		 	ADV = 0.
        		ENDIF

        		IF(OMG .EQ. 0.0) THEN
        		 	WRI = 0.
        		ELSE
        		 	WRI = 1.0 / (OMG*RAD)
        		ENDIF

        		CT = TP * WRI**2 * 2.0 / (RHO * PI * RAD**2)
        		CP = QP * WRI**2 * 2.0 / (RHO * PI * RAD**3)
        		DV = SQRT(VEL**2 + TP * 2.0/(RHO*PI*RAD**2)) - VEL

        		WRITE(LU,2100) CHARF, VEL,   RPM,  DBET,   TP,     QP, POWER, VOLT,AMPS,  EFFM,  EFFP, ADV, CT, CP, DV, EFF, PINPUT, PPROP, CLAVG, CDAVG
 				2100       FORMAT(A,F8.3,  G12.4,  F7.3, G12.4, G12.4, G12.4,F8.3, F10.4,  F9.4, F9.4, F10.5, G12.4, G12.4, F9.4,F9.4, G12.4, G12.4, F9.4, G12.4)
      		ENDDO
      		IF(.NOT.LRDUMP .AND. NVEL.GT.1) WRITE(LU,1000)
    	ENDDO
    ENDDO

    IF(LRDUMP) THEN

    	WRITE(LU,1100)
    	WRITE(LU,1100)' radius   chord   beta      Cl       Cd       Re    Mach','     effi     effp    Wa(m/s)     Aswirl      adv_wake'

    	DO I = 1, N
    	 	WA = VEL + VA(I)
    	 	WT = OMG*R(I) - VT(I)
    	 	WSQ = WA**2 + WT**2
    	 	W = SQRT(WSQ)

    	 	AMA = W/VSO
    	 	IRE = INT( RHO*W*C(I)/RMU + 0.5 )

    	 	IF(WA.NE.0.0 .AND. WT.NE.0.0) THEN
    	  		ADW = (WA/WT) * (R(I)/RAD)
    	  		EFFI = (VEL/(OMG*R(I))) * (WT/WA)
    	  		EFFP = (CL(I) - CD(I)*WA/WT)/(CL(I) + CD(I)*WT/WA)
    	 	ELSE
    	  		ADW = 0.
    	  		EFFI = 0.
    	  		EFFP = 0.
    	 	ENDIF

    		EFFI = MAX( -99.0 , MIN( 99.0 , EFFI ) )
    		EFFP = MAX( -99.0 , MIN( 99.0 , EFFP ) )

    		RU = R(I)
    		CU = C(I)
    		BU = B(I) * 180.0/PI

    		ASWIRL = ATAN2( VT(I) , WA ) * 180.0/PI

    		WRITE(LU,3100)  RU,   CU,   BU, CL(I), CD(I),IRE,  AMA, EFFI, EFFP, WA, ASWIRL, ADW
 				3100    FORMAT(1X,F8.4, F8.4, F8.3, F9.4,  F9.5,I9, F7.3, F9.4, F9.4, G12.4, G12.4, G12.4)
    	ENDDO
    ENDIF

    STOP

 	900  CONTINUE
    WRITE(*,9000) FILNAM(1:64), ILINE, LINE
 	9000 FORMAT(/' Read error'/'   in file:  ', A/'   line',I3,':  ', A)
    STOP

 	950  CONTINUE
    WRITE(*,9500) FILNAM(1:64), ILINE
 	9500 FORMAT(/' Unexpected end-of-file reached'/'   in file:  ', A/'   line',I3)
    STOP

 	980  CONTINUE
    WRITE(*,9500) FILNAM(1:64), ILINE
 	9800 FORMAT(/' Fewer parameters than required'/'   in file:  ', A/'   line',I3)
    STOP

END ! QPROP
         
