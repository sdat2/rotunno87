c*****************************************************************************
c*****************************************************************************
c
c
c                    ROTUNNO-EMANUEL 1987 MODEL
c
c   Rotunno, R., and K.A. Emanuel, 1987:  An air-sea interaction theory for
c      tropical cyclones, Part II: Evolutionary study using axisymmetric
c      nonhydrostatic numerical model.  J. Atmos. Sci., 44, 542-561.
c
c   Modified where indicated by K. Emanuel; last modified July, 2020
c
c
c*****************************************************************************
c*****************************************************************************
c*****************************************************************************
c
      PROGRAM HURAXI    
C
      PARAMETER(M=400,N=80,MP1=M+1,MM1=M-1,NP1=N+1,NM1=N-1)     
      PARAMETER (NTM=5000)
C
      DIMENSION  QRAD(M,N), UB(N),
     $ U(MP1,NP1),V(MP1,NP1),W(MP1,NP1),P(MP1,NP1),
     $ T(MP1,NP1)  , QV(MP1,NP1) , QL(MP1,NP1),VTERM2D(M,N),
     $ DUM(M,N), RC2(N)  , ZC2(N), CPTDR(N), CPTDZ(N), QS(N),
     $ RHOTVB(N), RHOTVW(NP1), RHOT(N),RHOW(NP1),RTB(N),RQVB(N),
     $ DTSR(M),DTSZ(M,N),DTSRW(N),DTSZW(N),DTLR(M),DTLZ(N),DTSV(M),
     $ A(N), B(N), C(N), D(M,N), E(N), WS(M,NP1), PS(M,N) 
     $,UAV(MP1,N),VAV(M,N),WAV(M,NP1),TAV(M,N),PAV(M,N),QVAV(M,N)
     $,QLAV(M,N), XKMAV(M,N), SST(M), ESS(M)
       REAL UMN(NTM),VMX(NTM),VMXS(NTM),WMX(NTM),
     $ RMX(NTM),PMN(NTM), TTP(NTM)
       REAL VMAXZ(NTM,N),THEHOV(NTM,M)
       REAL VBHOV(NTM,M),UBHOV(NTM,M),VTHOV(NTM,M),UTHOV(NTM,M)
       REAL RPL(M),UPL(M),VPL(M),PPL(M),THEPL(M)
       REAL ZPL(N), DUMT(M)
       character*18 chtime
       character*2 anp
       CHARACTER*1 PDEP,EVAP
C
      COMMON/DIFF/      
     $ UA(MP1,N) , WA(M,NP1), VA(M,N)  , TA(M,N), QVA(M,N), QLA(M,N),
     $ U1(MP1,N) , V1(M,N)  , W1(M,NP1), T1(M,N), QV1(M,N), QL1(M,N),
     $ XKM(M,NP1), P1(M,N)  , R(MP1) ,RS(M),Z(NP1),ZS(N),TAD(M,N),
     $ RTBW(N),RQVBW(N),TAU(N),TB(N), QVB(N),PN(NP1), PD(NP1), 
     $ TSURF(M),QSURF(M),RDR,RDZ,RDR2,RDZ2,DTS,DTL,DZ,
     $ XVL2, XHL2, CD, CE, XLDCP, XKAPPA, A1, G, 
     $ CERS(M), CDRS(M), DISS, VCAP, IVMAXS
C 
        DATA  CSTAR/30./,EP/0.1/,NZTR/1/,QVBT/0./
C
C   Open Parameter Input File and Read Parameters
C
        OPEN(UNIT=11,FILE='hurr.in',STATUS='old')
C
        READ(11,*)
        READ(11,*)
        READ(11,*)PDS
        READ(11,*)FC
        READ(11,*)CD
        READ(11,*)CE
        READ(11,*)CD1
        READ(11,*)CDCAP
        READ(11,*)CKCAP
        READ(11,*)VCAP
        READ(11,*)PDEP
        READ(11,*)XHL
        READ(11,*)XVL
        READ(11,*)DISS
        READ(11,*)TAUR
        READ(11,*)RADMAX
        READ(11,*)ACT
        READ(11,*)VTERM
        READ(11,*)VTERMSNOW
        READ(11,*)EVAP
        READ(11,*)TMIN
        READ(11,*)
        READ(11,*)
        READ(11,*)
        READ(11,*)RMAX
        READ(11,*)RO
        READ(11,*)VMAX
        READ(11,*)TMID
        READ(11,*)RSST
        READ(11,*)
        READ(11,*)
        READ(11,*)
        READ(11,*)DT
        READ(11,*)NS
        READ(11,*)EPS
        READ(11,*)ALPHA
        READ(11,*)RB
        READ(11,*)ZB
        READ(11,*)NSPONGE
        READ(11,*)ETIME
        READ(11,*)TAVE
        READ(11,*)PLTIME
        READ(11,*)TIMMAX
        READ(11,*)TIMEPL
        READ(11,*)ROG
        READ(11,*)ZOG
c
        CLOSE(11)
c
c-----------------------------------------------------------------------
c
c   Read input sounding file
c
        DZTEMP=ZB*1000./FLOAT(N)
        CALL INTERPOLATE(N,PDS,DZTEMP,TB,QVB,DISS,TBS,VMAXTH)
c
c  Modify VMAX according to Emanuel and Rotunno, J. Atmos. Sci., 2011
c
        CR1=CE/CD
        IF(CR1.EQ.2.0)CR1=2.01 ! Avoid singularity (singular limit)
        VMAXTH=VMAXTH*SQRT(2.0*CR1)*(0.5*CR1)**(CR1/(4.-2.*CR1))
c                 
        TBSW=TBS
        TBS=(TBS+273.)*(1000./PDS)**0.286
        TBT=2.0*TB(N)-TB(N-1)
c
c     open output file    
c
        open(17,file='output/texout.txt',
     1   status='unknown',form='formatted')
c
c  Define and/or normalize parameters
c
        F=FC*1.0E-5
        CD=CD*0.001
        CD1=CD1*1.0E-5
        CE=CE*0.001
        CDCAP=CDCAP*0.001
        CKCAP=CKCAP*0.001
        RADT=1.0/(3600.*TAUR)
        RADMAX=2.*DT*RADMAX/(3600.*24.)
        RMAX=RMAX*1000.0
        RO=RO*1000.0
        RSST=RSST*1000.0
        MPL=(FLOAT(M)*ROG)/RB  
        NPL=(FLOAT(N)*ZOG)/ZB 
        RB=1000.*RB
        ZB=1000.*ZB
        ISTOP=ETIME*24.*3600./DT
        NZSP=N-NSPONGE+1
        IAVE=TAVE*3600.*24./DT
        IPLOT=PLTIME*3600.*24./DT
        ITMAX=TIMMAX*3600./DT
        IPRINT=ISTOP/2-1
        TIMEPL=TIMEPL*3600.*24.
        ISTART=1
        AVE=FLOAT(IAVE)
        IPL=ISTART-1+IPLOT
        PDEPA=1.0
        IF(PDEP.EQ.'N'.OR.PDEP.EQ.'n')PDEPA=0.0
        ACT=0.001*ACT ! Normalize autoconversion threshold (July, 2020)
c
c     Added  9/3/2004 by KE
c
c	iave=min(iave,ipl)
c
      DR  = RB / FLOAT(M)       
      DZ  = ZB / FLOAT(N)       
      DTS = DT / FLOAT(NS)      
      DTL = DT  
      DZ2 = DZ * DZ     
      DR2 = DR * DR     
      RDR = 1. / DR     
      RDZ = 1. / DZ     
      RDR2= 1. / DR2    
      RDZ2= 1. / DZ2    
      PI  = 4. * ATAN(1.)       
      C2  = 90000.      
      CP  = 1005.       
      CV  = 718.
      RD  = 287.
      XKAPPA= RD / CP   
      XLDCP = 2500.     
      G     = 9.81      
      A1    = 7.5 * LOG(10.)   
      PNS   = (PDS/1000.) ** XKAPPA     
      TAMB  = TBS*PNS
      ESSS   =  6.11 * EXP( A1* (TAMB-273.)/( TAMB-36.)  )
      QSP   = .622 * ESSS /( PDS -ESSS )
      QVBS  = QVB(1)
      TVBS  = TBS * ( 1. + .61 * QVBS ) 
c      CC1 = .61 * G * DZ /(2.* CP * TVBS )/(1. + .61 *QVBS   )
c      CC2 = 1. - G * DZ / (2. * CP * TVBS * PNS )
c      CC5 = - TBS / PNS
      CC3   =  G * DZ / CP      
c
c   Modified by K. Emanuel, 11/4/2010
c
c     XHL2  = .04 * DR2 
      XHL2  = XHL*XHL
c
      XVL2  = XVL * XVL 
      TIME  = DT * ( FLOAT(ISTART) - 1. )       
C       
C        DEFINE GRID ARRAYS     
C       
      DO 10 I = 1 , MP1 
      R(I) = ( FLOAT(I) - 1.0 ) * DR    
      IF( I .NE. MP1 )  RS(I) = ( FLOAT(I) - 0.5 ) * DR 
   10 CONTINUE  
      DO 20 J = 1 , NP1 
      Z(J) = ( FLOAT(J) - 1.0 ) * DZ    
      IF( J .NE. NP1)  ZS(J) = ( FLOAT(J) - 0.5 ) * DZ  
   20 CONTINUE  
      TSTOP  = DT * FLOAT(ISTOP)
      write(17,920)  TIME, TSTOP, ISTART, ISTOP, IPLOT, IPRINT,    
     $ DT, EPS, XVL, F, RB, ZB, M, N, CD, CE    
  920 FORMAT(1H ,       
     $  9X,'TSTART=',F10.3,/,10X,'TSTOP =',F10.3,/,     
     $ 10X,'ISTART=',I5   ,/,10X,'ISTOP =',I5   ,/,     
     $ 10X,'IPLOT =',I5   ,/,10X,'IPRINT=',I5   ,/,     
     $ 10X,'DT    =',F10.3,/,10X,'EPS   =',F10.3,/,     
     $ 10X,'XL    =',F10.3,/,10X,'F     =',F7.5 ,/,     
     $ 10X,'RB    =',F10.1,/,10X,'ZB    =',F10.3,/,     
     $ 10X,'M     =',I3   ,/,10X,'N     =',I3   ,/,     
     $ 10X,'CD    =',F10.5,/,10X,'CE    =',F10.5   )    
C       
C        BASE STATE (HYDROSTATIC EQUATION FOR PN)       
C       
      PN(1) = PNS -.5*CC3/(.5*(TB(1)+TBS)*(1.+.61*.5*(QVB(1)+QVBS)))    
      PD(1) = 1000. * PN(1) ** ( 1. / XKAPPA )  
      DO 40 J = 2 , N   
       PN(J) = PN(J-1) - CC3/(.5*(TB(J)+TB(J-1))*
     $                          (1. + .61*.5*(QVB(J)+QVB(J-1)) ) )      
         TBAR=0.5*(PN(J)*TB(J)+PN(J-1)*TB(J-1))
         IF(TBAR.LE.TMIN)THEN
          PN(J)=PN(J-1)*EXP(-CC3/TMIN)
          TB(J)=TMIN/PN(J)
         END IF
       PD(J) = 1000. * PN(J) ** ( 1. / XKAPPA )  
   40 CONTINUE  
C   mod by KE  8/17/04
      PN(NP1)=2.*PN(N)-PN(N-1)
        PD(NP1)=2.*PD(N)-PD(N-1)
C  end mod
        TBAR=0.5*(PN(N-1)*TB(N-1)+PN(N)*TB(N))
      PNT = PN(N) -.5*CC3/(.5*(TB(N)+TBT)*(1.+.61*.5*(QVB(N)+QVBT)))    
        IF(TBAR.LE.TMIN)THEN
          PNT=PN(N)*EXP(-.5*CC3/TMIN)
        END IF
      PDT = 1000. * PNT ** ( 1. / XKAPPA )      
      write(17,945) 
  945 FORMAT(1H ,'     Z          P          PI        PT       ',
     $'QVB         QVS') 
      write(17,940) ZB, PDT, PNT, TBT, QVBT, 0.    
      DO 50 J = N , 2, -1       
      ES =6.11*EXP(A1*(PN(J)*TB(J)- 273.)/(PN(J)*TB(J)-36.))
      QS(J) = .622 * ES /MAX(( PD(J)-ES ),ES) 
      write(17,940) ZS(J), PD(J), PN(J), TB(J), QVB(J), QS(J) 
  940 FORMAT(1H , 6(F10.4,1X) ) 
   50 CONTINUE  
      write(17,940) 0., PDS, PNS, TBS, QVBS, QSP   
C       
C       ARRAYS FOR SMALL TIME STEP AND BUOYANCY CALCULATION     
C       
      DO 60 J = 1 , N   
      RC2(J)    =  RDR*DTS*C2/(CP*TB(J)*(1.+.61*QVB(J)))
      CPTDR(J)  =  DTS*RDR * CP * TB(J)* (1. + .61 * QVB(J) )   
      RTB(J)    = .5  * G / TB(J)       
c      RQVB(J)   = .61 * G * .5 / ( 1. + .61 * QVB(J) )  
      RHOTVB(J) = (1000./RD) * PN(J) ** (CV /RD)
      RHOT(J) = 100.* RHOTVB(J)/( TB(J)*( 1.+.61*QVB(J) ) )  
      ZC2(J)    = RDZ*RC2(J)/RHOTVB(J)/RDR      
      IF( J.EQ.1 )  GO TO 60    
      PNW = .5 * ( PN(J) + PN(J-1) )    
      TVW = .5*(TB(J)+TB(J-1) )*(1.+.305*(QVB(J)+QVB(J-1)))  
      RHOTVW(J)=  (1000./RD) * PNW ** (CV /RD)  
      RHOW(J)  = 100. * RHOTVW(J) / TVW 
      RTBW(J) = G / ( TB(J) + TB(J-1) ) 
      RQVBW(J)= .305*G/( 1. + .305*(QVB(J)+QVB(J-1)) )  
      CPTDZ(J)=DTS*RDZ*CP*.5*(TB(J)+TB(J-1))*
     1   (1.+.305*(QVB(J)+QVB(J-1)))
   60 CONTINUE  
      RHOTVW(1  )= (1000./RD) * PNS ** (CV /RD) 
      RHOTVW(NP1)= (1000./RD) * PNT ** (CV /RD) 
      RHOW(1  ) = 100. * RHOTVW(1  ) / TVBS     
      RHOW(NP1)=  100. * RHOTVW(NP1) / TBT      
C       
C         SPONGE LAYER DAMPING COEFFICIENT      
C       
      DO 65 J = 1 , N
      ZSP = ( ZS(J) - ZS(NZSP) ) / ( ZS(N) - ZS(NZSP) )
      IF(ZSP.LT.0.) TAU(J) = 0. 
      IF(ZSP.GE.0..AND.ZSP.LE..5) TAU(J)=-.5*ALPHA*(1.-COS(ZSP*PI) )    
      IF(ZSP.GT..5) TAU(J)=-.5*ALPHA*( 1. + PI*(ZSP-.5) )       
   65 CONTINUE  
C       
C         ARRAYS TO INCREASE EFFICIENCY 
C       
      DTSF = .5 * DTS * F       
      DTSG = .5 * DTS * G   
      DTL2 = .5 * DTL   
      DO 70 I = 1 , M   
      DTLR(I)  = .5  * DTL * RDR /RS(I) 
      DTSV(I)  = .5  * DTS  / RS(I)     
      IF( I .NE. 1 )  DTSR(I)  = .25 * DTS * RDR / R(I) 
   70 CONTINUE  
      DO 75 J = 1 , N   
      DTLZ(J)  = .5  * DTL * RDZ / RHOT(J)      
      DTSZW(J) = .25 * DTS * RDZ / RHOW(J)      
      DTSRW(J) = .25 * DTS * RDR / RHOW(J)      
   75 CONTINUE  
      DO 80 J = 1 , N   
      DO 80 I = 2 , M   
      DTSZ(I,J)  = .25 * DTS * RDZ /( RHOT(J) * R(I) )  
   80 CONTINUE  
C       
C        ARRAYS FOR SEMI - IMPLICIT SMALL TIME STEP     
C       
      CC4=.25*(1.+EP)**2
      DO 85 J = 2 , N   
      A(J)=   CC4*CPTDZ(J)*RHOTVW(J+1)*    ZC2(J)       
      B(J)=1.+CC4*CPTDZ(J)*RHOTVW(J  )*(ZC2(J)+ZC2(J-1))
      C(J)=   CC4*CPTDZ(J)*RHOTVW(J-1)*   ZC2(J-1)      
   85 CONTINUE  
      E(1)=0.   
      DO 95 I = 1 , M
      D(I,1)=0.   
95    CONTINUE
      DO 90 J = 2 , N   
      E(J) = A(J)/(B(J)-C(J)*E(J-1))    
   90 CONTINUE  
C
C       INITIAL CONDITIONS      
C
C
C     TEMPERATURE AND VAPOR PRESSURE AT SEA SURFACE
C
      DO 5 I = 1 , M
      SST(I)=TAMB+ TMID/(1.+(RS(I)/RSST)**2)
      ESS(I)=  6.11 * EXP( A1* (SST(I)-273.)/(SST(I)-36.)  )
 5    CONTINUE
      zd=zs(nzsp)
      E1 = 3.   
      D2 = 2. * RMAX / (    RO + RMAX ) 
      DO 100 J = 1 , N  
      DO 100 I = 1 , M  
      D1 = 2. * RMAX / ( RS(I) + RMAX ) 
      IF( RS(I) .LT. RO )THEN
       VR = SQRT(VMAX**2 *(RS(I)/RMAX)**2 *  
     $ ( D1 ** E1 - D2 ** E1 )+.25*F*F*RS(I)**2) - .5 * F * RS(I)       
      ELSE
       VR = 0.       
      END IF
      V(I,J )   = VR * ZS(J)/ZS(NZTR)
      IF( ZS(J).GT.ZS(NZTR) ) V(I,J)=VR*(ZS(J)-ZD)/(ZS(NZTR)-ZD)
      IF( ZS(J).GT.ZS(NZSP) ) V(I,J) = 0.0
      V1(I,J)   = V(I,J)
      QL(I,J)   = 0.    
      QL1(I,J)  = 0.    
  100 CONTINUE  
      DO 105 J = 1 , N  
      P(M,J) = 0.       
      P1(M,J)=P(M,J)    
      DO 105 I = M , 2 , -1     
      P(I-1,J) = P(I,J) - ( DTS / CPTDR(J) ) * .5 *      
     $ ( V(I,J)*V(I,J)/RS(I) + V(I-1,J)*V(I-1,J)/RS(I-1)
     $   + F * ( V(I,J) + V(I-1,J) ))
      P1(I-1,J)= P(I-1,J)       
  105 CONTINUE  
      DO 106 I = 1 , M  
      DO 107 J = 2 , NM1
      T(I,J) = TB(J) + .25*( CPTDZ(J+1)*(P(I,J+1)-P(I,J  ))     
     $  +CPTDZ(J  )*(P(I,J  )-P(I,J-1)) )/(DTS*RTB(J))  
      T1(I,J)=T(I,J)    
  107 CONTINUE  
      T(I,1) = TB(1) + .5*CPTDZ(2)*(P(I,2)-P(I,1  ))/(DTS*RTB(1))       
      T(I,N) = TB(N) + .5*CPTDZ(N)*(P(I,N)-P(I,NM1))/(DTS*RTB(N))       
      T1(I,1)=T(I,1)    
      T1(I,N)=T(I,N)    
  106 CONTINUE  
C
C    GIVE 0'S TO THOSE ARRAY ELEMENTS THAT ARE OUTSIDE PHYSICAL BOUNDARIES (IE
C    TO ELEMENTS THAT WERE ADDED TO ARRAY DEFINITIONS TO AVOID INDEX OVERFLOW)
C
      DO 141 I=1,MP1
         V(I,NP1)=0.0
         P(I,NP1)=0.0
         T(I,NP1)=0.0
         QV(I,NP1)=0.0
         QL(I,NP1)=0.0
         U(I,NP1)=0.0
  141  CONTINUE
      W(MP1,NP1)=0.0
      DO 142 J=1,NP1
         V(MP1,J)=0.0
         P(MP1,J)=0.0
         T(MP1,J)=0.0
         QV(MP1,J)=0.0
         QL(MP1,J)=0.0
         W(MP1,J)=0.0
  142  CONTINUE
       U(MP1,NP1)=0.0
c
c
c    ADJUST Qv SO THAT THETAe= constant across vortex
c
c
      DO 108 J = 1 , N  
      DO 108 I = 1 , M  
c      QV(I,J)   = QVB(J)+PN(J)*(TB(J)-T(I,J))/XLDCP
      QV(I,J)   = QVB(J)
c       IF(R(I).LT.150000.0.AND.J.GT.1)THEN
c         PNST= PN(J) + P1(I,J) - P1(M,1)
c         PDST= 1000. * PNST ** (1./XKAPPA)
c         TEMP= PNST * T1(I,J)
c         ES  = 6.11 * EXP(A1 * ( TEMP  - 273. ) / ( TEMP  - 36. ) )
c         QV(I,J) = 0.98*.622 * ES /(PDST-ES)
c       END IF
      QV1(I,J)  = QV(I,J)
 108  CONTINUE  
      DO 110 J = 1 , N  
      DO 110 I = 1 , MP1
      U(I,J)  = 0.0     
      U1(I,J) = U(I,J)  
  110 CONTINUE  
      DO 120 J = 1 , NP1
      DO 120 I = 1 , M
      W(I,J)  = 0.0     
      W1(I,J) = W(I,J)  
  120 CONTINUE  
      DO 505 J = 1 , N  
      UA(1,J) = 0.      
  505 CONTINUE  
        TIPL=0.0
        NPT=0
C
C   Write some arrays
C
        DO I=1,MPL
         RPL(I)=R(I)*0.001
        END DO
        DO 821 I=1,NPL
         ZPL(I)=0.001*DZ*FLOAT(I-1)
  821   CONTINUE
        OPEN(UNIT=22,FILE='output/rgraph.out',STATUS='UNKNOWN')
        OPEN(UNIT=23,FILE='output/zgraph.out',STATUS='UNKNOWN')
        WRITE(22,1025)(RPL(I),I=1,MPL)
 1025   FORMAT(2X,800(F12.7,' '))
        WRITE(23,1030)(ZPL(I),I=1,NPL)
 1030   FORMAT(2X,800(F8.4,1X))
        CLOSE(22)
        CLOSE(23)
        NPLOT=0
C       
C        BEGIN TIME MARCH       
C       
      DO 500 ITIME = ISTART , ISTOP     
      TIME = TIME + DT  
      TIPL=TIPL + DT
C       
C       SURFACE TEMP AND WATER-VAPOR MIXING RATIO       
C       
      DO 515 I = 1 , M  
      PNSURF = PDEPA*(P1(I,1)-P1(M,1))/(1.-.5*G*DZ/(CP*SST(I)))
      PDST   = 1000. * (PNS +PNSURF) ** (1./XKAPPA)
      QSURF(I) = .622 * ESS(I) /(PDST -ESS(I))
      TSURF(I) = SST(I)/(PNS+PNSURF)
  515 CONTINUE  
      DO 6 I = 1 , M
        CDRS(I)=CD+CD1
     $    *SQRT(.25*(U1(I,1)+U1(I+1,1))**2+V1(I,1)**2) 
        CDRS(I)=MIN(CDRS(I),CDCAP)
        CERS(I)=CE+CD1
     $    *SQRT(.25*(U1(I,1)+U1(I+1,1))**2+V1(I,1)**2)
        CERS(I)=MIN(CERS(I),CKCAP)
C
    6 CONTINUE
C
      CALL DIFFUSE 
C
C     RADIATION (R < 1 DEG/DAY)                           
C     
      DO 150 J = 1 , N
      DO 150 I = 1 , M
      TDIF=T1(I,J)-TB(J)
        QRAD(I,J) = -DTL*TDIF*RADT
        QRAD(I,J)=MIN(QRAD(I,J),RADMAX)
        QRAD(I,J)=MAX(QRAD(I,J),-RADMAX)
  150 CONTINUE                             
C       
C        FORCING FOR U EQUATION 
C       
      DO 510 J = 1 , N  
      DO 510 I = 2 , M  
      IF(J.EQ.1)THEN
        UMINUS=U(I,J)
      ELSE
        UMINUS=U(I,J-1)
      END IF
      UA(I,J)    = UA(I,J)      
     1 - DTSR(I)   *(   
     2 ( R(I+1)*U(I+1,J)+R(I  )*U(I  ,J) )*( U(I+1,J)-U(I  ,J) )
     3+( R(I  )*U(I  ,J)+R(I-1)*U(I-1,J) )*( U(I  ,J)-U(I-1,J) ) )      
     4 - DTSZ(I,J) *(   
     5 RHOW(J+1)*(RS(I)*W(I,J+1)+RS(I-1)*W(I-1,J+1))*(U(I,J+1)-U(I,J))  
     6+RHOW(J  )*(RS(I)*W(I,J)+RS(I-1)*W(I-1,J))*(U(I,J)-UMINUS)) 
     6+ DTSV(I) * V(I,J)*V(I,J) +  DTSV(I-1) * V(I-1,J)*V(I-1,J)
     7  +  DTSF    * ( V(I,J) + V(I-1,J) )      
  510 CONTINUE  
C       
C       OUTER BOUNDARY  
C       
      DO 525 J = 1 , N  
      UA(MP1,J) =       
     $ -MAX( U(MP1,J) + CSTAR , 0. )*DTL*RDR*(U1(MP1,J)-U1(M,J))
     $+ DTL * (V(M,J)*V(M,J)/RS(M) +  F * V(M,J)   )    
  525 CONTINUE  
C       
C       DRAG AT J = 1   
C       
      UA(MP1,1) = UA(MP1,1)     
     $ - CD *DTL*RDZ*U1(MP1,1)*SQRT( U1(MP1,1)**2+V1(M,1)**2 ) 
C       
C       FORCING FOR W EQUATION  
C       
      DO 550 J = 2 , N  
      DO  I = 2 , MM1
      WA(I,J)    = WA(I,J)      
     1 -DTSRW(J)*(      
     2 (RHOT(J  )*U(I+1,J)+RHOT(J-1)*U(I+1,J-1))*(W(I+1,J)-W(I  ,J))    
     3+(RHOT(J  )*U(I  ,J)+RHOT(J-1)*U(I  ,J-1))*(W(I  ,J)-W(I-1,J)))  
     4 -DTSZW(J)*(      
     5 (RHOW(J+1)*W(I,J+1)+RHOW(J  )*W(I  ,J  ))*(W(I,J+1)-W(I,J  ))    
     6+(RHOW(J-1)*W(I,J-1)+RHOW(J  )*W(I  ,J  ))*(W(I,J  )-W(I,J-1)))  
     7   + DTS*RTBW(J)  * ( T(I,J) - TB(J)  + T(I,J-1)  - TB(J-1)  )    
     8   + DTS*RQVBW(J) * (QV(I,J) - QVB(J) + QV(I,J-1) - QVB(J-1) )    
     9   - DTSG * ( QL(I,J) +QL(I,J-1) )
      END DO
      I=1
      WA(I,J)    = WA(I,J)      
     1 -DTSRW(J)*(      
     2 (RHOT(J  )*U(I+1,J)+RHOT(J-1)*U(I+1,J-1))*(W(I+1,J)-W(I  ,J)))    
     4 -DTSZW(J)*(      
     5 (RHOW(J+1)*W(I,J+1)+RHOW(J  )*W(I  ,J  ))*(W(I,J+1)-W(I,J  ))    
     6+(RHOW(J-1)*W(I,J-1)+RHOW(J  )*W(I  ,J  ))*(W(I,J  )-W(I,J-1)))  
     7   + DTS*RTBW(J)  * ( T(I,J) - TB(J)  + T(I,J-1)  - TB(J-1)  )    
     8   + DTS*RQVBW(J) * (QV(I,J) - QVB(J) + QV(I,J-1) - QVB(J-1) )    
     9   - DTSG * ( QL(I,J) +QL(I,J-1) )
  550 CONTINUE
C       
C       OUTER BOUNDARY  
C       
      DO 520 J = 2 , N  
      UB(J) =      MAX( RHOT(J  ) * ( U(MP1,J) + U(M,J  ) )   
     $                   +RHOT(J-1) * ( U(M,J-1) + U(M,J-1) )  , 0.)  
  520 CONTINUE  
      I = M     
      DO 555 J = 2 , N  
      WA(I,J)    = WA(I,J)      
     1 -DTSRW(J) * UB(J) * ( W1(I  ,J) - W1(I-1,J) )    
     4 -DTSZW(J)*(      
     5 (RHOW(J+1)*W(I,J+1)+RHOW(J  )*W(I  ,J  ))*(W(I,J+1)-W(I,J  ))    
     6+(RHOW(J-1)*W(I,J-1)+RHOW(J  )*W(I  ,J  ))*(W(I,J  )-W(I,J-1)))  
     7   + DTS*RTBW(J)  * ( T(I,J) - TB(J)  + T(I,J-1)  - TB(J-1)  )    
     8   + DTS*RQVBW(J) * (QV(I,J) - QVB(J) + QV(I,J-1) - QVB(J-1) )    
     9   - DTSG * ( QL(I,J) +QL(I,J-1) )
  555 CONTINUE  
C       
C       FORCING FOR  V, T, QV, QL  EQUATIONS    
C       
      DO 540 J = 1 , N  
      DO  I = 2 , MM1
      IF(J.EQ.1)THEN
       VMINUS=V(I,J)
       TMINUS=T(I,J)
       TBMINUS=TB(J)
       QVMINUS=QV(I,J)
       QLMINUS=QL(I,J)
      ELSE
       VMINUS=V(I,J-1)
       TMINUS=T(I,J-1)
       TBMINUS=TB(J-1)
       QVMINUS=QV(I,J-1)
       QLMINUS=QL(I,J-1)
      END IF
      VA(I,J)    = VA(I,J)      
     1 -DTLR(I)* ( R(I+1) * U(I+1,J) * ( V(I+1,J  ) - V(I  ,J  ) )      
     2            +R(I  ) * U(I  ,J) * ( V(I  ,J  ) - V(I-1,J  ) ) )    
     3 -DTLZ(J)* ( RHOW(J+1)*W(I,J+1)* ( V(I  ,J+1) - V(I  ,J  ) )      
     4            +RHOW(J  )*W(I,J  )* ( V(I  ,J  ) - VMINUS ) )    
     5 -DTL2*(F+ V(I,J)/RS(I))*(R(I+1)*U(I+1,J)+R(I)*U(I,J))/RS(I)      
      TA(I,J)    = TA(I,J)      
     1 -DTLR(I)* ( R(I+1) * U(I+1,J) * ( T(I+1,J  ) - T(I  ,J  ) )      
     2            +R(I  ) * U(I  ,J) * ( T(I  ,J  ) - T(I-1,J  ) ) )    
c     3 -DTLZ(J)* ( RHOW(J+1)*W(I,J+1)* (
c     & T(I  ,J+1) -TB(J+1) - T(I  ,J  ) +TB(J) )      
c     4            +RHOW(J  )*W(I,J  )* ( 
c     & T(I  ,J  ) -TB(J)   - TMINUS +TBMINUS ) )    
     3 -DTLZ(J)* ( RHOW(J+1)*W(I,J+1)* ( T(I  ,J+1) - T(I  ,J  ) )      
     4            +RHOW(J  )*W(I,J  )* ( T(I  ,J  ) - TMINUS ) )    
     5 +QRAD(I,J)
      QVA(I,J)    = QVA(I,J)    
     1 -DTLR(I)* ( R(I+1) * U(I+1,J) * ( QV(I+1,J  ) - QV(I  ,J  ) )    
     2            +R(I  ) * U(I  ,J) * ( QV(I  ,J  ) - QV(I-1,J  ) ) )  
     3 -DTLZ(J)* ( RHOW(J+1)*W(I,J+1)* ( QV(I  ,J+1) - QV(I  ,J  ) )    
     4            +RHOW(J  )*W(I,J  )* ( QV(I  ,J  ) - QVMINUS ) )  
      QLA(I,J)    = QLA(I,J)    
     1 -DTLR(I)* ( R(I+1) * U(I+1,J) * ( QL(I+1,J  ) - QL(I  ,J  ) )    
     2            +R(I  ) * U(I  ,J) * ( QL(I  ,J  ) - QL(I-1,J  ) ) )  
     3 -DTLZ(J)* ( RHOW(J+1)*W(I,J+1)* ( QL(I  ,J+1) - QL(I  ,J  ) )    
     4            +RHOW(J  )*W(I,J  )* ( QL(I  ,J  ) - QLMINUS ) )  
      END DO
c
c  Inner boundary
c
      I=1
c
      IF(J.EQ.1)THEN
       VMINUS=V(I,J)
       TMINUS=T(I,J)
       TBMINUS=TB(J)
       QVMINUS=QV(I,J)
       QLMINUS=QL(I,J)
      ELSE
       VMINUS=V(I,J-1)
       TMINUS=T(I,J-1)
       TBMINUS=TB(J-1)
       QVMINUS=QV(I,J-1)
       QLMINUS=QL(I,J-1)
      END IF
      VA(I,J)    = VA(I,J)      
     1 -DTLR(I)* ( R(I+1) * U(I+1,J) * ( V(I+1,J  ) - V(I  ,J  ) ) ) 
     3 -DTLZ(J)* ( RHOW(J+1)*W(I,J+1)* ( V(I  ,J+1) - V(I  ,J  ) )      
     4            +RHOW(J  )*W(I,J  )* ( V(I  ,J  ) - VMINUS ) )    
     5 -DTL2*(F+ V(I,J)/RS(I))*(R(I+1)*U(I+1,J)+R(I)*U(I,J))/RS(I)      
      TA(I,J)    = TA(I,J)      
     1 -DTLR(I)* ( R(I+1) * U(I+1,J) * ( T(I+1,J  ) - T(I  ,J  ) ) )
c     3 -DTLZ(J)* ( RHOW(J+1)*W(I,J+1)* (
c     & T(I  ,J+1) -TB(J+1) - T(I  ,J  ) +TB(J) )      
c     4            +RHOW(J  )*W(I,J  )* ( 
c     & T(I  ,J  ) -TB(J)   - TMINUS +TBMINUS ) )    
     3 -DTLZ(J)* ( RHOW(J+1)*W(I,J+1)* ( T(I  ,J+1) - T(I  ,J  ) )      
     4            +RHOW(J  )*W(I,J  )* ( T(I  ,J  ) - TMINUS ) )    
     5 +QRAD(I,J)
      QVA(I,J)    = QVA(I,J)    
     1 -DTLR(I)* ( R(I+1) * U(I+1,J) * ( QV(I+1,J  ) - QV(I  ,J  ) ) ) 
     3 -DTLZ(J)* ( RHOW(J+1)*W(I,J+1)* ( QV(I  ,J+1) - QV(I  ,J  ) )    
     4            +RHOW(J  )*W(I,J  )* ( QV(I  ,J  ) - QVMINUS ) )  
      QLA(I,J)    = QLA(I,J)    
     1 -DTLR(I)* ( R(I+1) * U(I+1,J) * ( QL(I+1,J  ) - QL(I  ,J  ) ) )
     3 -DTLZ(J)* ( RHOW(J+1)*W(I,J+1)* ( QL(I  ,J+1) - QL(I  ,J  ) )    
     4            +RHOW(J  )*W(I,J  )* ( QL(I  ,J  ) - QLMINUS ) )  
c
  540 CONTINUE  
C       
C       OUTER BOUNDARY  
C       
      DO 530 J = 1 , N  
      UB(J) =  MAX( R(MP1) * U(MP1,J) + R(M) * U(M,J) , 0. )  
  530 CONTINUE  
      I = M     
      DO 545 J = 1 , N  
      IF(J.EQ.1)THEN
       VMINUS=V(I,J)
       TMINUS=T(I,J)
       TBMINUS=TB(J)
       QVMINUS=QV(I,J)
       QLMINUS=QL(I,J)
      ELSE
       VMINUS=V(I,J-1)
       TMINUS=T(I,J-1)
       TBMINUS=TB(J-1)
       QVMINUS=QV(I,J-1)
       QLMINUS=QL(I,J-1)
      END IF
      VA(I,J)    = VA(I,J)      
     1 -DTLR(I) * UB(J) * ( V1(I  ,J  ) - V1(I-1,J  ) ) 
     3 -DTLZ(J)* ( RHOW(J+1)*W(I,J+1)* ( V(I  ,J+1) - V(I  ,J  ) )      
     4            +RHOW(J  )*W(I,J  )* ( V(I  ,J  ) - VMINUS ) )    
     5 -DTL2*(F+ V(I,J)/RS(I))*(R(I+1)*U(I+1,J)+R(I)*U(I,J))/RS(I)      
      TA(I,J)    = TA(I,J)      
     1 -DTLR(I) * UB(J) * ( T1(I  ,J  ) - T1(I-1,J  ) ) 
c     3 -DTLZ(J)* ( RHOW(J+1)*W(I,J+1)* (
c     & T(I  ,J+1) -TB(J+1) - T(I  ,J  ) +TB(J) )      
c     4            +RHOW(J  )*W(I,J  )* ( 
c     & T(I  ,J  ) -TB(J)   - TMINUS +TBMINUS ) )    
     3 -DTLZ(J)* ( RHOW(J+1)*W(I,J+1)* ( T(I  ,J+1) - T(I  ,J  ) )      
     4            +RHOW(J  )*W(I,J  )* ( T(I  ,J  ) - TMINUS ) )    
     5 +QRAD(I,J)
      QVA(I,J)    = QVA(I,J)    
     1 -DTLR(I) * UB(J) * ( QV1(I  ,J  ) - QV1(I-1,J  ) )       
     3 -DTLZ(J)* ( RHOW(J+1)*W(I,J+1)* ( QV(I  ,J+1) - QV(I  ,J  ) )    
     4            +RHOW(J  )*W(I,J  )* ( QV(I  ,J  ) - QVMINUS ) )  
      QLA(I,J)    = QLA(I,J)    
     1 -DTLR(I) * UB(J) * ( QL1(I  ,J  ) - QL1(I-1,J  ) )       
     3 -DTLZ(J)* ( RHOW(J+1)*W(I,J+1)* ( QL(I  ,J+1) - QL(I  ,J  ) )    
     4            +RHOW(J  )*W(I,J  )* ( QL(I  ,J  ) - QLMINUS ) )  
  545 CONTINUE  
C
C         RAIN (AND SNOW) FALL
C
      DO 570 J = 1 , N
      DO 570 I = 1 , M
c
      VTERM2D(I,J)=VTERM
c
      IF(VTERM.LT.1.0E-6)THEN ! If VTERM = 0 calculate fall speed
c
       VTERM2D(I,J)=36.34*SQRT(1.1/RHOT(J))*(RHOT(J)*QL(I,J))**0.1364
       PNST= PN(J) + P1(I,J) - P1(M,1)
       PDST= 1000. * PNST ** (1./XKAPPA)
       TEMP= PNST * T1(I,J) -273.0
       IF(TEMP.LT.0.0)THEN  !  Weight toward set fall speed of snow ast T approaches -20C
         SWEIGHT=(TEMP+20.0)/20.0
         SWEIGHT=MAX(SWEIGHT,0.0)
         VTERM2D(I,J)=SWEIGHT*VTERM2D(I,J)+VTERMSNOW*(1.-SWEIGHT)
       END IF
       VTERM2D(I,J)=MIN(VTERM2D(I,J),12.0)  ! Put cap on fall speed of rain
c
      END IF
c
      DUM(I,J)=0.
      IF((QL(I,J) - ACT).GT.0.0) DUM(I,J)=VTERM2D(I,J)
  570 CONTINUE
      DO 575 J = 2 , NM1
      DO 575 I = 1 , M
      QLA(I,J) = QLA(I,J) + DTLZ(J) *
     $(DUM(I,J+1)*RHOT(J+1)*QL(I,J+1)-DUM(I,J-1)*RHOT(J-1)*QL(I,J-1)) 
  575 CONTINUE
c
      DO 580 I = 1 , M
      DUM(I,2)=0.
c
      VTERM2D(I,J)=VTERM
c
      IF(VTERM.LT.1.0E-6)THEN ! If VTERM = 0 calculate fall speed
c
      VTERM2D(I,2)=36.34*SQRT(1.1/RHOT(2))*(RHOT(2)*QL1(I,2))**0.1364
      VTERM2D(I,1)=36.34*SQRT(1.1/RHOT(1))*(RHOT(2)*QL1(I,1))**0.1364
c
      END IF
c
      DUM(I,2)=0.0
      IF((QL1(I,2) - ACT).GT.0.0) DUM(I,2)=VTERM2D(I,2)
      DUM(I,1)=0.
      IF((QL1(I,1) - ACT).GT.0.0) DUM(I,1)=VTERM2D(I,1)
c
      QLA(I,1) = QLA(I,1) + 2.*DTLZ(1)*
     $( DUM(I,2)*RHOT(2)*QL1(I,2)  - DUM(I,1)*RHOT(1)*QL1(I,1 ) )
c
  580 CONTINUE
C       
C        TIME SMOOTHER  
C       
c      cfnew= dts/dtl

      DO 160 J = 1 , N  
      DO 160 I = 1 , M  
      U(I+1,J)= U(I+1,J)+ EPS * ( U1(I+1,J)-2.*U(I+1,J) )       
      V(I,J)  = V(I,J)  + EPS * ( V1(I,J)  -2.*V(I,J)   )       
      W(I,J+1)= W(I,J+1)+ EPS * ( W1(I,J+1)-2.*W(I,J+1) )       
      P(I,J)  = P(I,J)  + EPS * ( P1(I,J)  -2.*P(I,J)   )       
      T(I,J)  = T(I,J)  + EPS * ( T1(I,J)  -2.*T(I,J)   )       
      QV(I,J) = QV(I,J) + EPS * ( QV1(I,J) -2.*QV(I,J)  )       
      QL(I,J) = QL(I,J) + EPS * ( QL1(I,J) -2.*QL(I,J)  )       
  160 CONTINUE  
C       
C         SMALL TIME STEP       
C       
      DO 170 NSMALL = 1 , NS    
      DO 171 J = 1 , N  
      DO 171 I = 2 , M  
c  Modified by K.E. 4/15/02
c  171 U1(I,J)=U1(I,J)-CPTDR(J)*(P1(I,J)-P1(I-1,J))+UA(I,J)      
        CFAC=DTS*RDR*CP*0.5*(T1(I,J)*(1.+0.61*QV(I,J))+T1(I-1,J)*
     1 (1.+0.61*QV(I-1,J)))
  171 U1(I,J)=U1(I,J)-CFAC*(P1(I,J)-P1(I-1,J))+UA(I,J)      
c      do 1761, j=1,N
c      do 1761, i=1,M
c        t1(i,j)= t1(i,j) +cfnew*( ta(i,j) 
c     & -DTLZ(J)* ( RHOW(J+1)*W1(I,J+1)* ( tb(j+1) - tb(j) )
c     &            +RHOW(J  )*W1(I,J  )* ( tb(j) - tb(j-1) ) ) )
c 1761 continue
      DO 172 J = 2 , N  
      DO 172 I = 1 , M  
c  Modified by K.E. 4/15/02
C  172 WS(I,J)=W1(I,J)-.5*(1.-EP)*CPTDZ(J)*(P1(I,J)-P1(I,J-1))+WA(I,J)   
        CFAC=DTS*RDZ*CP*0.5*(T1(I,J)+T1(I,J-1))*
     1   (1.+0.305*(QV(I,J)+QV(I,J-1)))
  172 WS(I,J)=W1(I,J)-.5*(1.-EP)*CFAC*(P1(I,J)-P1(I,J-1))+WA(I,J)   
c     7   + DTS*RTBW(J)  * ( T1(I,J) - TB(J)  + T1(I,J-1)  - TB(J-1))    
      DO 173 J = 1 , N  
      DO 173 I = 1 , M  
      PS(I,J)=P1(I,J)-RC2(J) *  
     1   ( R(I+1)*U1(I+1,J)  -R(I)*U1(I,J)  )/RS(I)     
     2               -.5*(1.-EP)*ZC2(J) *       
     3    ( RHOTVW(J+1)*W1(I,J+1)-RHOTVW(J)*W1(I,J))    
  173 CONTINUE  
      DO 175 J = 2 , N  
      DO 175 I = 1 , M  
c  Modified by K.E. 4/15/02
C      D(I,J) = (WS(I,J) - CPTDZ(J)*.5*(1.+EP)*(PS(I,J)-PS(I,J-1)) 
C     $                + C(J) * D(I,J-1) ) * E(J) /A(J)    
       CFAC=DTS*RDZ*CP*0.5*(T1(I,J)+T1(I,J-1))*
     1  (1.+0.305*(QV(I,J)+QV(I,J-1)))
      D(I,J) = (WS(I,J) - CFAC*.5*(1.+EP)*(PS(I,J)-PS(I,J-1)) 
     $                + C(J) * D(I,J-1) ) * E(J) /A(J)    
  175 CONTINUE  
      DO 174 J = N , 1 , -1     
       DO 174 I = 1 , M  
        W1(I,J) = E(J) * W1(I,J+1) + D(I,J) 
        P1(I,J) = PS(I,J) -.5*(1.+EP)*ZC2(J)*     
     $    ( RHOTVW(J+1)*W1(I,J+1)-RHOTVW(J)*W1(I,J))    
  174  CONTINUE  
  170 CONTINUE  
C       
C        OUTER BOUNDARY 
C       
      DO 195 J = 1 , N  
       U1(MP1,J) = U1(MP1,J) + UA(MP1,J) 
  195 CONTINUE  
C       
C        ADVANCE V , T , QV , QL
C       
      DO 490 J=1,N      
      DO 490 I=1,M      
      V1(I,J)  = V1(I,J)  + VA(I,J)     
      T1(I,J)  = T1(I,J)  + TA(I,J)     
      QV1(I,J) = QV1(I,J) + QVA(I,J)    
      QL1(I,J) = QL1(I,J) + QLA(I,J)    
      IF(QV1(I,J).LT.0.   ) QV1(I,J)=0. 
      IF(QL1(I,J).LT.0.   ) QL1(I,J)=0. 
C       
C        CONDENSATION / EVAPORATION     
C           
      PNST= PN(J) + P1(I,J) - P1(M,1)   
      PDST= 1000. * PNST ** (1./XKAPPA) 
      TEMP= PNST * T1(I,J)      
      ES  = 6.11 * EXP(A1 * ( TEMP  - 273. ) / ( TEMP  - 36. ) )
      QSS = .622 * ES /MAX((PDST-ES),ES)  
      IF( QV1(I,J) .LT. QSS .AND. QL1(I,J).LE. 1.E-8 ) GO TO 480
      R1 = 1./(1. + XLDCP * 237. * A1 * QSS /( TEMP - 36. )**2 )
      QVD = R1 * ( QV1(I,J) - QSS )     
c
      IF(QVD.LT.0.0)THEN
C
      IF( EVAP .EQ. 'N' .OR. EVAP .EQ. 'n')THEN ! No Evaporation
C
        QVD=0.
C
      ELSEIF(  EVAP .EQ. 'I' .OR. EVAP .EQ. 'i' )THEN ! Instantaneous evaporation (per RE(1987))
c
      IF ( (QL1(I,J)   + QVD) .LT. 0. ) THEN    
C
C         EVAPORATE     
C       
      T1 (I,J) = T1(I,J)  - XLDCP * QL1(I,J) /PNST      
      QV1(I,J) = QV1(I,J) + QL1(I,J)    
      QL1(I,J) = 0.     
c
      ENDIF
c
      ELSE
C
C         EVAPORATE     
C       
c  Added July 3 2020
c
      eamount=0.5*DTL*(1.-QV1(I,J)/QSS)*sqrt(max(QL1(I,J),0.0))/ ! Evaporation scheme (Emanuel, 1991)
     1  (2.0e3+1.0e4/(PDST*QSS))
      eamount=min(eamount,QL1(I,J))
c
      T1 (I,J) = T1(I,J)  - XLDCP *eamount  /PNST      
      QV1(I,J) = QV1(I,J) + eamount
      QL1(I,J) = QL1(I,J) - eamount
c
      END IF
c
      ELSE
c
C         CONDENSE      
C       
      T1(I,J)   = T1(I,J)  + XLDCP * QVD / PNST 
      QV1(I,J)  = QV1(I,J)         - QVD
      QL1(I,J)  = QL1(I,J)         + QVD
c
      END IF
c
  480 CONTINUE  
c
      IF(QV1(I,J).LT.0.   ) QV1(I,J)=0. 
      IF(QL1(I,J).LT.0.   ) QL1(I,J)=0. 
  490 CONTINUE  
C       
C        TIME FLIP      
C       
      DO 190 J=1,N      
      DO 190 I=1,M      
      UTEMP    = U(I+1,J) + EPS * U1(I+1,J)     
      VTEMP    = V(I,J)   + EPS * V1(I,J)       
      WTEMP    = W(I,J+1) + EPS * W1(I,J+1)     
      PTEMP    = P(I,J)   + EPS * P1(I,J)       
      TTEMP    = T(I,J)   + EPS * T1(I,J)       
      QVTEMP   = QV(I,J)  + EPS * QV1(I,J)      
      QLTEMP   = QL(I,J)  + EPS * QL1(I,J)      

      U(I+1,J) = U1(I+1,J)      
      V(I,J)   = V1(I,J)
      W(I,J+1) = W1(I,J+1)      
      P(I,J)   = P1(I,J)
      T(I,J)   = T1(I,J)
      QV(I,J)  = QV1(I,J)       
      QL(I,J)  = QL1(I,J)       

      U1(I+1,J)= UTEMP  
      V1(I,J)  = VTEMP  
      W1(I,J+1)= WTEMP  
      T1(I,J)  = TTEMP  
      P1(I,J)  = PTEMP  
      QV1(I,J) = QVTEMP 
      QL1(I,J) = QLTEMP 
  190 CONTINUE  
      IF( mod(itime,ipl-iave+1).eq.0) then
                NPLOT=NPLOT+1
        IF(NPLOT.NE.1)THEN
            nmplot=nplot-1
            WRITE(ANP,191)NMPLOT
  191 FORMAT(I2.2)        
      OPEN(UNIT=15,FILE='output/vcon'//anp//'.out',STATUS='UNKNOWN')
      OPEN(UNIT=16,FILE='output/ucon'//anp//'.out',STATUS='UNKNOWN')
      OPEN(UNIT=18,FILE='output/wcon'//anp//'.out',STATUS='UNKNOWN')
      OPEN(UNIT=19,FILE='output/tcon'//anp//'.out',STATUS='UNKNOWN')
      OPEN(UNIT=20,FILE='output/pcon'//anp//'.out',STATUS='UNKNOWN')
      OPEN(UNIT=21,FILE='output/tecon'//anp//'.out',STATUS='UNKNOWN')
      OPEN(UNIT=25,FILE='output/liqcon'//anp//'.out',STATUS='UNKNOWN')
      OPEN(UNIT=28,FILE='output/tescon'//anp//'.out',STATUS='UNKNOWN')
      OPEN(UNIT=29,FILE='output/qcon'//anp//'.out',STATUS='UNKNOWN')
      OPEN(UNIT=30,FILE='output/tfcon'//anp//'.out',STATUS='UNKNOWN')
        TT = TAV(1,1) * ( PN(1) + PAV(1,1) - PAV(1,1) )
        TEMAX=TAV(1,1) * EXP( XLDCP * QVAV(1,1) / TT )
        TEMAX=TEMAX+10.0
        DO J=1,NPL
         WRITE(15,1035)(VAV(I,J),I=1,MPL)
         DO I=1,MPL
          DUMT(I)=.5 * ( UAV(I+1,J) + UAV(I,J) )
         END DO
         WRITE(16,1035)(DUMT(I),I=1,MPL)
       DO I=1,MPL
          DUMT(I)=.5 * ( WAV(I,J+1) + WAV(I,J) )
         END DO
         WRITE(18,1035)(DUMT(I),I=1,MPL)
c       DO I=1,MPL
c          DUMT(I)=TAV(I,J)*(PN(J)+PAV(I,J)-PAV(M,1)) - TB(J) * PN(J)
c          DUMT(I)=TAV(I,J)*(PN(J)+PAV(I,J)-PAV(M,1)) 
c         END DO
c         WRITE(19,1035)(DUMT(I),I=1,MPL)
       DO I=1,MPL
          DUMT(I)=1000.*QLAV(I,J)
         END DO
         WRITE(25,1035)(DUMT(I),I=1,MPL)
c       DO I=1,MPL
c          DUMT(I)=1000.*PN(J)**(1./XKAPPA) * (PAV(I,J)-PAV(M,1))/XKAPPA
c         END DO
c         WRITE(20,1035)(DUMT(I),I=1,MPL)
       DO I=1,MPL
        TT = TAV(I,J) * ( PN(J) + PAV(I,J) - PAV(M,1) )
          DUMT(I)=TAV(I,J) * EXP( XLDCP * QVAV(I,J) / TT )
c          DUMT(I)=MIN(DUMT(I),TEMAX)
         END DO
       WRITE(21,1035)(DUMT(I),I=1,MPL)
       DO I=1,MPL
        TT = TAV(I,J) * ( PN(J) + PAV(I,J) - PAV(M,1) )
          TTC=TT-273.15
          ESTT=6.112*EXP(17.67*TTC/(243.5+TTC))
          PFULL=1000.*( PN(J) + PAV(I,J)- PAV(M,1))**(1./XKAPPA)
          QSTT=0.622*ESTT/MAX(1.0,(PFULL-ESTT))
          DUMT(I)=TAV(I,J) * EXP( XLDCP * QSTT / TT )
c          DUMT(I)=MIN(DUMT(I),TEMAX)
         END DO
         WRITE(28,1035)(DUMT(I),I=1,MPL)
         DO I=1,MPL
          TT = TAV(I,J) * ( PN(J) + PAV(I,J) - PAV(M,1) )
          DUMT(I)=TT
         END DO
         WRITE(30,1035)(DUMT(I),I=1,MPL)
        END DO
c
        DO I=1,M
         DO J=1,N
          DUMT(J)=TAV(I,J)*(PN(J)+PAV(I,J)-PAV(M,1))
         END DO
         WRITE(19,1035)(DUMT(J),J=1,N)
         DO J=1,N
          DUMT(J)=1000.*( PN(J) + PAV(I,J)- PAV(M,1))**(1./XKAPPA)
         END DO
         WRITE(20,1035)(DUMT(J),J=1,N)
         DO J=1,N
          DUMT(J)=1000.*QVAV(I,J)
         END DO
         WRITE(29,1035)(DUMT(J),J=1,N)
        END DO
c
        CLOSE(15)
        CLOSE(16)
        CLOSE(18)
        CLOSE(19)
        CLOSE(20)
        CLOSE(21)
        CLOSE(25)
        CLOSE(28)
        CLOSE(29)
        CLOSE(30)
        END IF
c
      DO 375 J = 1 , N
      DO 375 I = 1 , M
      UAV(I+1,J) = 0.0
      VAV(I,J)   = 0.0
      WAV(I,J+1) = 0.0
      TAV(I,J)   = 0.0
      PAV(I,J)   = 0.0
      QVAV(I,J)  = 0.0
      QLAV(I,J)  = 0.0
      XKMAV(I,J) = 0.0
 375  CONTINUE
      ipl=ipl+iplot
      end if
      if(itime.ge.ipl-iplot-iave+1.and.itime.le.ipl-iplot) then
      DO 370 I = 1 , M
      DO 370 J = 1 , N
      UAV(I+1,J) = UAV(I+1,J) + U1(I+1,J)/AVE
      VAV(I,J)   = VAV(I,J)   + V1(I,J)/AVE
      WAV(I,J+1) = WAV(I,J+1) + W1(I,J+1)/AVE
      TAV(I,J)   = TAV(I,J)   + T1(I,J)/AVE
      PAV(I,J)   = PAV(I,J)   + P1(I,J)/AVE
      QVAV(I,J)  = QVAV(I,J)  + QV1(I,J)/AVE
      QLAV(I,J)  = QLAV(I,J)  + QL1(I,J)/AVE
      XKMAV(I,J) = XKMAV(I,J) + XKM(I,J)/AVE
  370 CONTINUE
      end if
C       
C      COMPUTE WMAX , VMAX      
C       
      IWMAX=1   
      JWMAX=1   
      WMAX=0.   
      IVMAX=1   
      IVMAXS=1
      JVMAX=1   
      VMAX=0.   
      UMIN=0.0
       VMAXS=0.0
       DO I=1,M
         IF(V(I,1).GT.VMAXS)THEN
        VMAXS=V(I,1)
         IVMAXS=I
         END IF
       END DO
      DO 210 J = 2 , N  
      DO 210 I = 1 , M  
      IF(W(I,J).GT.WMAX) THEN   
      WMAX = W(I,J)     
      IWMAX = I 
      JWMAX = J 
      END IF    
      IF(V(I,J).GT.VMAX) THEN   
      VMAX = V(I,J)     
      IVMAX = I 
      JVMAX = J 
      END IF    
      UMIN=MIN(UMIN,U(I,2))
  210 CONTINUE  
      TINHRS = TIME / 3600.     
      IF( ITIME.EQ.ISTART ) GO TO 301
      IF(MOD(ITIME,ITMAX).EQ.0) 
     $ write(17,925) TINHRS,WMAX,IWMAX,JWMAX,VMAX,IVMAX,JVMAX      
  925 FORMAT(1H ,'TIME = ',F9.1,2X,'WMAX = ',E11.4,2X,'IW= ',I3,
     $2X,'JW= ',I2,2X,'VMAX = ',E11.4,2X,'IV= ',I3,2X,'JV= ',I3)
C       
C        PRINT FIELDS   
C       
      IF(MOD(ITIME,IPRINT).NE.0) GO TO 200      
  301 CONTINUE
      TINHRS = TIME / 3600.     
      WRITE(17, 915) TINHRS 
      WRITE(17, 900) 'U'    
      DO 300 J=N,1,-1   
  300 WRITE(17, 910) (U(I,J)  ,I=1,20,2)      
      WRITE(17, 900) 'W'    
      DO 305 J=N,2,-1   
  305 WRITE(17, 910) (W(I,J)  ,I=1,20,2)      
      WRITE(17, 900) 'V'    
      DO 306 J =N,1,-1  
  306 WRITE(17, 910) (V(I,J)  ,I=1,20,2)      
      WRITE(17, 900) 'THTA' 
      DO 310 J=N,1,-1   
  310 WRITE(17, 910) ( T(I,J) ,I=1,20,2 )       
      WRITE(17, 900) 'TSUR' 
      WRITE(17, 910) (TSURF(I),I=1,20,2)       
      WRITE(17, 900) 'TEMP' 
      DO 311 J=N,1,-1   
  311 WRITE(17, 910) ( T(I,J)*(PN(J)+P(I,J)-P(M,1))-TB(J)*PN(J) 
     &                ,I=1,20,2 )   
      WRITE(17, 900) 'QV'   
      DO 315 J=N,1,-1   
  315 WRITE(17, 910) (1000.*QV(I,J)  ,I=1,20,2)       
      WRITE(17, 900) 'QSUR' 
      WRITE(17, 910) (1000.*QSURF(I) ,I=1,20,2)       
      WRITE(17, 900) 'QL'   
      DO 316 J=N,1,-1   
  316 WRITE(17, 910) (1000.*QL(I,J)  ,I=1,20,2)       
      WRITE(17, 900) 'P'    
      DO 320 J=N,1,-1   
      FAC = PN(J)**1./XKAPPA  
  320 WRITE(17, 910) ( 1000.*FAC*(P(I,J)-P(M,1))/XKAPPA ,I=1,20,2)
      WRITE(17, 900) 'KM'   
      DO 325 J=N,1,-1   
  325 WRITE(17, 910) ( .01*XKM(I,J) ,I=1,20,2)
      WRITE(17, 900) 'TE'   
      DO 330 J =N,1,-1  
      TT = T(I,J) * ( PN(J) + P(I,J) - P(M,1) ) 
  330 WRITE(17, 910) (T(I,J) * (1. + XLDCP * QV(I,J) / TT ),I=1,20,2)
  200 CONTINUE  
C       
      TINHRS = TIME/3600.       
      write(chtime,812) tinhrs
  812 format('   TIME=',F10.4)
      IF(ITIME.GT.1) GO TO 499  
      DTL = 2. * DT     
      NS  = 2  * NS     
      DTL2 = .5 * DTL   
      DO 17 I = 1 , M   
   17 DTLR(I)  = .5  * DTL * RDR /RS(I) 
      DO 27 J = 1 , N   
   27 DTLZ(J)  = .5  * DTL * RDZ / RHOT(J)      
  499 CONTINUE
        IF(TIPL.GE.TIMEPL)THEN
         NPT=NPT+1
         TTP(NPT)=TIME/3600.0
         VMX(NPT)=VMAX
         VMXS(NPT)=VMAXS
         UMN(NPT)=UMIN
         WMX(NPT)=WMAX*100.
         PMN(NPT)= PDS+1000.*PN(1)**(1./XKAPPA) *
     1    (P(1,1)-P(M,1))/XKAPPA
         RMX(NPT)=0.001*DR*FLOAT(IVMAX-1)
          DO J=1,N
           VMAXZ(NPT,J)=0.0
           DO I=1,M
            VMAXZ(NPT,J)=MAX(VMAXZ(NPT,J),V(I,J))
           END DO
          END DO
           JOUT=15000.0/DZ
          DO I=1,MPL
           TT = T(I,1) * ( PN(1) + P(I,1) - P(M,1) )
           THEHOV(NPT,I)= T(I,1) * EXP( XLDCP * QV(I,1) / TT )
           VBHOV(NPT,I)=V(I,1)
           UBHOV(NPT,I)=U(I,1)
           VTHOV(NPT,I)=V(I,JOUT)
           UTHOV(NPT,I)=U(I,JOUT)
          END DO
         TIPL=0.0
        END IF
C
  500 CONTINUE  
C
        DO 820 I=1,MPL
         UPL(I)=UAV(I,1)
         VPL(I)=VAV(I,1)
         PPL(I)=PDS+1000.*PN(1)**(1./XKAPPA)*
     1    (PAV(I,1)-PAV(M,1))/XKAPPA
         TT = TAV(I,1) * ( PN(1) + PAV(I,1) - PAV(M,1) )
         THEPL(I)= TAV(I,1) * EXP( XLDCP * QVAV(I,1) / TT )
  820   CONTINUE
  900 FORMAT(/,1X,A4,/) 
  910 FORMAT(1H ,16F7.2)
  915 FORMAT(//,1H ,'TIME=',F10.3)      
         OPEN(UNIT=13,FILE='output/time.out',STATUS='UNKNOWN')
         OPEN(UNIT=14,FILE='output/radius.out',STATUS='UNKNOWN')
         OPEN(UNIT=15,FILE='output/vcon.out',STATUS='UNKNOWN')
         OPEN(UNIT=16,FILE='output/ucon.out',STATUS='UNKNOWN')
         OPEN(UNIT=18,FILE='output/wcon.out',STATUS='UNKNOWN')
         OPEN(UNIT=19,FILE='output/tcon.out',STATUS='UNKNOWN')
         OPEN(UNIT=20,FILE='output/pcon.out',STATUS='UNKNOWN')
         OPEN(UNIT=21,FILE='output/tecon.out',STATUS='UNKNOWN')
         OPEN(UNIT=24,FILE='output/vmaxz.out',STATUS='UNKNOWN')
         OPEN(UNIT=25,FILE='output/liqcon.out',STATUS='UNKNOWN')
         OPEN(UNIT=27,FILE='output/thehov.out',STATUS='UNKNOWN')
         OPEN(UNIT=28,FILE='output/tescon.out',STATUS='UNKNOWN')
         OPEN(UNIT=29,FILE='output/s.out',STATUS='UNKNOWN')
         OPEN(UNIT=30,FILE='output/tfcon.out',STATUS='UNKNOWN')
         OPEN(UNIT=31,FILE='output/xkcon.out',STATUS='UNKNOWN')
          OPEN(UNIT=32,FILE='output/vbhov.out',STATUS='UNKNOWN')
          OPEN(UNIT=33,FILE='output/vthov.out',STATUS='UNKNOWN')
          OPEN(UNIT=34,FILE='output/ubhov.out',STATUS='UNKNOWN')
          OPEN(UNIT=35,FILE='output/uthov.out',STATUS='UNKNOWN')
        DO 1010 I=1,NPT
         WRITE(13,1005)TTP(I),VMX(I),UMN(I),WMX(I),PMN(I),RMX(I),
     1    VMAXTH,VMXS(I)
          WRITE(24,2006)(VMAXZ(I,J),J=1,NPL)
 1005  FORMAT(2X,7(F12.6,' '),F12.6)
 1010        CONTINUE
 2006 FORMAT(50(1X,F12.6))
        DO 1020 I=1,MPL
         WRITE(14,1015)RPL(I),VPL(I),UPL(I),PPL(I),THEPL(I)
 1015         FORMAT(2X,4(F12.6,' '),F12.6)
 1020        CONTINUE
        TT = TAV(1,1) * ( PN(1) + PAV(1,1) - PAV(1,1) )
        TEMAX=TAV(1,1) * EXP( XLDCP * QVAV(1,1) / TT )
        TEMAX=TEMAX+10.0
c
        WRITE(29,1031)
 1031   FORMAT(9X,'Sea Surface Temperature (C)')
        WRITE(29,1032)
 1032   FORMAT(9X,'---------------------------')       
        WRITE(29,1033)TBSW
 1033   FORMAT(19X,F5.2)    
        WRITE(29,*)
        WRITE(29,1034)
 1034   FORMAT(4X,'p (hPa)         T (C)            q (g/Kg)')
        WRITE(29,1036)   
 1036   FORMAT(4X,'-----------------------------------------')      
        WRITE(29,*)   
C               
        DO 1040 J=1,NPL
         WRITE(15,1035)(VAV(I,J),I=1,MPL)
         DO I=1,MPL
          DUMT(I)=.5 * ( UAV(I+1,J) + UAV(I,J) )
         END DO
         WRITE(16,1035)(DUMT(I),I=1,MPL)
         DO I=1,MPL
          DUMT(I)=.5 * ( WAV(I,J+1) + WAV(I,J) )
         END DO
         WRITE(18,1035)(DUMT(I),I=1,MPL)
         DO I=1,MPL
          DUMT(I)=.5 * ( XKMAV(I,J+1) + XKMAV(I,J) )
         END DO
         WRITE(31,1035)(DUMT(I),I=1,MPL)
         DO I=1,MPL
          DUMT(I)=TAV(I,J)*(PN(J)+PAV(I,J)-PAV(M,1)) - TB(J) * PN(J)
         END DO
         WRITE(19,1035)(DUMT(I),I=1,MPL)
         DO I=1,MPL
          DUMT(I)=1000.*QLAV(I,J)
         END DO
         WRITE(25,1035)(DUMT(I),I=1,MPL)
         DO I=1,MPL
          DUMT(I)=1000.*PN(J)**(1./XKAPPA)*(PAV(I,J)-PAV(M,1))/XKAPPA
         END DO
         WRITE(20,1035)(DUMT(I),I=1,MPL)
         DO I=1,MPL
          TT = TAV(I,J) * ( PN(J) + PAV(I,J) - PAV(M,1) )
          DUMT(I)=TAV(I,J) * EXP( XLDCP * QVAV(I,J) / TT )
c          DUMT(I)=MIN(DUMT(I),TEMAX)
         END DO
         WRITE(21,1035)(DUMT(I),I=1,MPL)
         DO I=1,MPL
          TT = TAV(I,J) * ( PN(J) + PAV(I,J) - PAV(M,1) )
          TTC=TT-273.15
          ESTT=6.112*EXP(17.67*TTC/(243.5+TTC))
          PFULL=1000.*( PN(J) + PAV(I,J)- PAV(M,1))**(1./XKAPPA)
          QSTT=0.622*ESTT/MAX(1.0,(PFULL-ESTT))
          DUMT(I)=TAV(I,J) * EXP( XLDCP * QSTT / TT )
c          DUMT(I)=MIN(DUMT(I),TEMAX)
         END DO
         WRITE(28,1035)(DUMT(I),I=1,MPL)
         TT = TAV(M,J) * ( PN(J) + PAV(M,J) - PAV(M,1) )-273.15
         PFULL=1000.*( PN(J) + PAV(M,J)- PAV(M,1))**(1./XKAPPA)
          WRITE(29,1041) PFULL,TT,1000.*QVAV(M,J)
         DO I=1,MPL
          TT = TAV(I,J) * ( PN(J) + PAV(I,J) - PAV(M,1) )
          DUMT(I)=TT
         END DO
         WRITE(30,1035)(DUMT(I),I=1,MPL)
 1040    CONTINUE
        DO I=1,NPT
         WRITE(27,1035)(THEHOV(I,J),J=1,MPL)
         WRITE(32,1035)(vbhov(I,J),J=1,MPL)
         WRITE(33,1035)(vthov(I,J),J=1,MPL)
         WRITE(34,1035)(ubhov(I,J),J=1,MPL)
         WRITE(35,1035)(uthov(I,J),J=1,MPL)
        END DO
 1035   FORMAT(2X,F12.7,' ',800(F8.2,' '))
 1041   FORMAT(4X,F8.4,5X,F9.4,11X,F8.5)
C
      STOP      
      END       
C
      SUBROUTINE DIFFUSE   
C
      PARAMETER(M=400,N=80,MP1=M+1,MM1=M-1,NP1=N+1,NM1=N-1)     
      DIMENSION 
     $ TRR(M,N  )  , TTT(M,N)  , TZZ(M,NP1),    
     $ TRZ(MP1,NP1), TRT(MP1,N), TZT(M,NP1),    
     $ TR(MP1, N)  , TZ(M,NP1) ,
     $ QVR(MP1, N) , QVZ(M,NP1), QLR(MP1,N), QLZ(M,NP1),
     $ XKMH(M,NP1) , TE(M,N)   , UDR(M,N)  , WDZ(M,N)  , VDR(MP1,N), 
     $ UZWR(MP1,N) , VDZ(M,N)  , DEFH2(M,N), DEF2(M,N) , DTDZ(M,N)
      COMMON/DIFF/      
     $ UA(MP1,N),WA(M,NP1),VA(M,N),TA(M,N),QVA(M,N),QLA(M,N),
     $ U1(MP1,N),V1(M,N),W1(M,NP1),T1(M,N), QV1(M,N), QL1(M,N),
     $ XKM(M,NP1),P1(M,N),R(MP1),RS(M),Z(NP1),ZS(N),TAD(M,N),
     $ RTBW(N),RQVBW(N),TAU(N),TB(N),QVB(N),PN(NP1),PD(NP1),
     $ TSURF(M),QSURF(M),RDR,RDZ,RDR2,RDZ2,DTS,DTL,DZ,
     $ XVL2, XHL2, CD, CE, XLDCP, XKAPPA, A1, G, 
     $ CERS(M), CDRS(M), DISS, VCAP, IVMAXS
      DO 500 J = 1 , N  
      DO 500 I = 1 , M  
      TE(I,J)= T1(I,J) * (1.+ XLDCP * QV1(I,J) /
     $        ( T1(I,J) * ( PN(J) + P1(I,J) - P1(M,1) )  )  )   
  500 CONTINUE  
      DO 110 J = 1 , N  
      DO 110 I = 1 , M  
      UDR(I,J) = RDR * ( U1(I+1,J) - U1(I,J) )  
      WDZ(I,J) = RDZ * ( W1(I,J+1) - W1(I,J) )  
  110 CONTINUE  
      DO 120 J = 1 , N  
      DO 125 I = 2 , M  
      VDR(I,J) = R(I) * RDR * ( V1(I,J)/RS(I) - V1(I-1,J)/RS(I-1) )     
  125 CONTINUE  
      VDR(1  ,J) = 0.   
      VDR(MP1,J) = VDR(M,J)     
  120 CONTINUE  
      DO 130 J = 2 , N  
      DO 135 I = 2 , M  
      UZWR(I,J) = RDZ*(U1(I,J)-U1(I,J-1)) + RDR*(W1(I,J)-W1(I-1,J))     
  135 CONTINUE  
      UZWR(1  ,J) = 0.  
      UZWR(MP1,J) = UZWR(M,J)   
  130 CONTINUE  
      DO 140 J = 2 , N  
      DO 140 I = 1 , M  
      VDZ(I,J) = RDZ * ( V1(I,J) - V1(I,J-1) )  
  140 CONTINUE  
      DO 150 J = 2 , N  
      DO 150 I = 1 , M  
      DEFH2(I,J) = UDR(I,J)*UDR(I,J) + UDR(I,J-1)*UDR(I,J-1)    
     $ + ( .25/(RS(I)*RS(I)) )*(
     $  (U1(I+1,J  )+U1(I,J  )) *  (U1(I+1,J  )+U1(I,J  ))      
     $ +(U1(I+1,J-1)+U1(I,J-1)) *  (U1(I+1,J-1)+U1(I,J-1)) )    
     $ + .25*( VDR(I+1,J  )*VDR(I+1,J  ) + VDR(I,J  )*VDR(I,J  )
     $        +VDR(I+1,J-1)*VDR(I+1,J-1) + VDR(I,J-1)*VDR(I,J-1) )      
  150 CONTINUE  
      DO 160 J = 2 , N  
      DO 160 I = 1 , M  
      DEF2(I,J) = DEFH2(I,J) + 
     $    WDZ(I,J)*WDZ(I,J) + WDZ(I,J-1)*WDZ(I,J-1)     
     $ +.5*( UZWR(I+1,J)*UZWR(I+1,J) + UZWR(I,J)*UZWR(I,J) )    
     $ + VDZ(I,J) * VDZ(I,J)    
c
c   Line below added 9/29/2010 to change effective Richardson Number
c
c	DEF2(I,J)=2.*DEF2(I,J)
c
  160 CONTINUE  
      DO 170 J = 2 , N  
      DO 170 I = 1 , M  
      DTDZ(I,J) = RDZ*( 
     $ 2.*RTBW(J) * (T1(I,J) - T1(I,J-1)  )     
     $+2.*RQVBW(J)* (QV1(I,J)- QV1(I,J-1) ) )   
      IF( (QL1(I,J) * QL1(I,J-1)) .GE. 1.E-8 ) THEN
      TT = .5 * (  T1(I,J  ) * ( PN(J  )+P1(I,J  )-P1(M,1) )    
     $           + T1(I,J-1) * ( PN(J-1)+P1(I,J-1)-P1(M,1) )  ) 
      AA=(1. + 4362.   *(QV1(I,J)+QV1(I,J-1)) / TT      )       
     $  /(1. + 6738953.*(QV1(I,J)+QV1(I,J-1))/(TT*TT)   )       
     $  * (2. / (TB(J)+TB(J-1)) )       
      DTDZ(I,J) = RDZ* G * ( AA * ( TE(I,J) - TE(I,J-1) )       
     $ - ( QL1(I,J) + QV1(I,J) - QL1(I,J-1) - QV1(I,J-1) ) )    
      END IF    
  170 CONTINUE  
      DO 180 J = 2 , N  
      DO 180 I = 1 , M  
      XKM(I,J) = 0.     
c
c  Factor of 0.5 inserted here to mimic larger Richardson Number, 9/29/2010
c
      IF(DEF2(I,J).GT.DTDZ(I,J))
     $  XKM(I,J)= XVL2 * SQRT( DEF2(I,J) - DTDZ(I,J))   
      XKMH(I,J) = XHL2 * SQRT( DEFH2(I,J) )     
      IF(XKM(I,J).GE. .4 * DZ *DZ/DTL) XKM(I,J) = .4 *DZ*DZ/DTL 
      IF(XKMH(I,J).LT.XKM(I,J))   XKMH(I,J) = XKM(I,J)  
  180 CONTINUE  
      DO 190 I = 1 , M  
      XKM(I,  1) = XKM(I,2)     
      XKM(I,NP1) = XKM(I,N)     
      XKMH(I,  1) = XKMH(I,2)   
      XKMH(I,NP1) = XKMH(I,N)   
  190 CONTINUE  
C       
C        CALCULATE STRESS       
C       
C        TRR(M,N)       
C       
      DO 200 J = 1 , N  
      DO 200 I = 1 , M  
      TRR(I,J) =  ( XKMH(I,J+1) + XKMH(I,J) ) * UDR(I,J)
  200 CONTINUE  
C       
C        TTT(M,N)       
C       
      DO 210 J = 1 , N  
      TTT(1,J) = 0.     
      DO 210 I = 2 , M  
      TTT(I,J) =
     $.5*(XKMH(I,J+1)+XKMH(I,J)+XKMH(I-1,J+1)+XKMH(I-1,J))*U1(I,J)/R(I)
  210 CONTINUE  
C       
C        TRZ(M,NP1)     
C       
      DO 220 J=2,N      
      DO 220 I=2,M      
      TRZ(I,J)=.5* ( XKM(I-1,J) + XKM(I,J) ) *  UZWR(I,J) 
  220 CONTINUE  
      DO 230 I =2,M     
      TRZ(I,1  )=.5*(CDRS(I)+CDRS(I-1))
     $  *U1(I,1)*SQRT(U1(I,1)**2+.25*(V1(I,1)+V1(I-1,1))**2) 
      TRZ(I,NP1) = 0.   
  230 CONTINUE  
      DO 235 J = 1 , NP1
      TRZ(1  ,J) = 0.   
      TRZ(MP1,J) = TRZ(M,J) * R(M) / R(MP1)     
  235 CONTINUE  
C       
C        TRT(M,N)       
C       
      DO 245 J=1,N      
      DO 240 I=2,M      
      TRT(I,J) =
     $.25*( XKMH(I,J+1)+XKMH(I,J)+XKMH(I-1,J+1)+XKMH(I-1,J))* VDR(I,J)
  240 CONTINUE  
      TRT(1  ,J) = 0.   
      TRT(MP1,J) = TRT(M,J) * R(M)**2 / R(MP1)**2       
  245 CONTINUE  
C       
C        TZT(M,NP1)     
C       
      DO 250 I = 1 , M  
      DO 255 J = 2 , N  
      TZT(I,J) = XKM(I,J) * VDZ(I,J)   
  255 CONTINUE  
      TZT(I,1  )=CDRS(I)
     $   *V1(I,1)*SQRT(.25*(U1(I,1)+U1(I+1,1))**2+V1(I,1)**2)  
      TZT(I,NP1)= 0.    
  250 CONTINUE  
C       
C          TZZ(M,N)     
C       
      DO 260 J = 1 , N  
      DO 260 I = 1 , M  
      TZZ(I,J) = ( XKM(I,J+1) + XKM(I,J) )*WDZ(I,J) 
  260 CONTINUE  
C       
C       TEMPERATURE FLUX
C       
C          TR(MP1,N)    
C       
      DO 270 J=1,N      
      DO 275 I=2,M      
      TR(I,J) = .25*( XKMH(I,J+1)+XKMH(I,J)+XKMH(I-1,J+1)+XKMH(I-1,J))  
     $                          *RDR*(T1(I,J)-T1(I-1,J))
  275 CONTINUE  
      TR(1  ,J)=0.0     
      TR(MP1,J)= TR(M,J) * R(M) / R(MP1)
  270 CONTINUE  
C       
C         TZ(M,NP1)     
C       
      DO 280 I=1,M      
      DO 285 J=2,N      
      TZ(I,J) = XKM(I,J) * RDZ * (T1(I,J) - T1(I,J-1))  
c   Mod by K.E.  4/12/02
        TZ(I,J)=TZ(I,J)*PD(J)
  285 CONTINUE  
        VABS=SQRT(.25 * ( U1(I+1,1) + U1(I,1) ) **2 + V1(I,1)**2 )
c
        VABS=MIN(VABS, VCAP)
        TZ(I,1  ) =( T1(I,1)-TSURF(I))*CERS(I)*VABS
c
c
c   Experiment of 2/4/2016
c
c       IM=MAX(I-1,1)
c       VABSM=SQRT(.25 * ( U1(I,1) + U1(IM,1) ) **2 + V1(IM,1)**2 )
c       TZ(I,1  ) =( T1(IM,1)-TSURF(IM))*CERS(IM)*VABSM
c
c   Mod by K.E.  4/12/02
        TZ(I,1)=TZ(I,1)*PD(1)
      TZ(I,NP1) = 0.    
  280 CONTINUE  
C       
C       QV FLUX 
C       
C         QVR(MP1,N)    
C       
      DO 370 J=1,N      
      DO 375 I=2,M      
      QVR(I,J)= .25*(XKMH(I,J+1)+XKMH(I,J)+XKMH(I-1,J+1)+XKMH(I-1,J))   
     $                            *RDR*(QV1(I,J)-QV1(I-1,J))    
  375 CONTINUE  
      QVR(1  ,J)=0.0    
      QVR(MP1,J)= QVR(M,J) * R(M) / R(MP1)      
  370 CONTINUE  
C       
C         QVZ(M,NP1)    
C       
      DO 380 I=1,M      
      DO 385 J=2,N      
      QVZ(I,J) = XKM(I,J) * RDZ * (QV1(I,J) - QV1(I,J-1))       
c   Mod by K.E.  4/12/02
        QVZ(I,J)=QVZ(I,J)*PD(J)/PN(J)
  385 CONTINUE  
      VABS=SQRT(.25 * ( U1(I+1,1) + U1(I,1) ) **2 + V1(I,1)**2 )
c
      VABS=MIN(VABS, VCAP)
      QVZ(I,1  ) =( QV1(I,1)-QSURF(I))*CERS(I)*VABS      
c
c
c
c   Experiment of 2/4/2016
c
c       IM=MAX(I-1,1)
c       VABSM=SQRT(.25 * ( U1(IM,1) + U1(IM,1) ) **2 + V1(IM,1)**2 )
c       QVZ(I,1  ) =( QV1(IM,1)-QSURF(IM))*CERS(IM)*VABSM
c
c   Mod by K.E.  4/12/02
      QVZ(I,1)=QVZ(I,1)*PD(1)/PN(1)
      QVZ(I,NP1) = 0.   
  380 CONTINUE  
C       
C       QL FLUX 
C       
C        QLR(MP1,N)     
C       
      DO 570 J=1,N      
      DO 575 I=2,M      
      QLR(I,J)= .25*(XKMH(I,J+1)+XKMH(I,J)+XKMH(I-1,J+1)+XKMH(I-1,J))   
     $                          *RDR*(QL1(I,J)-QL1(I-1,J))      
  575 CONTINUE  
      QLR(1  ,J)=0.0    
      QLR(MP1,J)= QLR(M,J) * R(M) / R(MP1)      
  570 CONTINUE  
C       
C         QLZ(M,NP1)    
C       
      DO 580 I=1,M      
      DO 585 J=2,N      
      QLZ(I,J) = XKM(I,J) * RDZ * (QL1(I,J) - QL1(I,J-1))       
  585 CONTINUE  
      QLZ(I,1  ) = 0.   
      QLZ(I,NP1) = 0.   
  580 CONTINUE  
      DO 400 J = 1 , N  
      DO 400 I = 2 , M  
      UA(I,J)= -DTS*(   
     1 - RDR * ( RS(I) * TRR(I,J) - RS(I-1) * TRR(I-1,J) ) / R(I)       
     2 + TTT(I,J) / R(I)
     3 - RDZ * ( TRZ(I,J+1) - TRZ(I,J) )   )    
     4 + DTS * TAU(J) * U1(I,J) 
  400 CONTINUE  
      DO 410  J = 2 , N 
      DO 410  I = 1 , M 
      WA(I,J) = -DTS*(  
     1 - RDR * ( R(I+1)*TRZ(I+1,J) - R(I)*TRZ(I,J) ) / RS(I)    
     2 - RDZ * ( TZZ(I,J) - TZZ(I,J-1) )   )    
     3 + DTS * .5 * ( TAU(J) + TAU(J-1) ) * W1(I,J)     
  410 CONTINUE  
      DO 420  J = 1 , N 
      DO 420  I = 1 , M 
      VA(I,J) = -DTL * (
     1 - RDR * (R(I+1)*R(I+1)*TRT(I+1,J)-R(I)*R(I)*TRT(I,J))    
     2                                           /(RS(I)*RS(I)) 
     3 - RDZ * ( TZT(I,J+1) - TZT(I,J) )   )    
     4 +DTL * TAU(J) * V1(I,J)  
c   Mod by K.E.  4/12/02
      TA(I,J) = -DTL * (
     1 - RDR * ( R(I+1) * TR(I+1,J) - R(I) * TR(I,J) ) / RS(I) 
c     2 - RDZ * (  TZ(I,J+1) - TZ(I,J) )    )    
     2 - RDZ * (  TZ(I,J+1) - TZ(I,J) )/(0.5*(PD(J)+PD(J+1))))    
     3 + DTL * TAU(J) * ( T1(I,J) - TB(J) )     
c   Mod by K.E.  4/12/02
      QVA(I,J) = -DTL * (       
     1 - RDR * ( R(I+1) * QVR(I+1,J) - R(I) * QVR(I,J) ) / RS(I)
C     2 - RDZ * (  QVZ(I,J+1) - QVZ(I,J) )    )  
     2 - RDZ * (  QVZ(I,J+1) - QVZ(I,J) )/
     * (0.5*(PD(J)/PN(J)+PD(J+1)/PN(J+1)))   )  
     3 + DTL * TAU(J) *  QV1(I,J)       
      QLA(I,J) = -DTL * (       
     1 - RDR * ( R(I+1) * QLR(I+1,J) - R(I) * QLR(I,J) ) / RS(I)
     2 - RDZ * (  QLZ(I,J+1) - QLZ(I,J) )    )  
     3 + DTL * TAU(J) *  QL1(I,J)       
  420 CONTINUE  
c
c       Dissipative Heating
c
        DO 430 I=1,M-1
         DO 425 J=2,N-1
C         TA(I,J)=TA(I,J)-(0.25*(U1(I+1,J)+U1(I,J))*(UA(I+1,J)+UA(I,J))
C     1    +V1(I,J)*VA(I,J))/(PN(J)*1005.0)
         TAD(I,J)=DTL*(TZT(I,J+1)*VDZ(I,J+1)+XKM(I,J)*0.25*(UZWR
     1    (I,J)+UZWR(I+1,J))**2+XKMH(I,J+1)*(1./16.)*(VDR(I,J)+VDR(
     2   I+1,J)+VDR(I,J+1)+VDR(I+1,J+1))**2+XKMH(I,J+1)*0.25*(UDR(I,J)+
     3    UDR(I,J+1))**2)/(PN(J)*1005.0)
         TA(I,J)=TA(I,J)+DISS*TAD(I,J)
  425        CONTINUE        
         TAD(I,1)=DTL*(RDZ*CDRS(I)*(0.25*(U1(I,1)+U1(I+1,1))**2+
     1    V1(I,1)**2)**1.5+XKMH(I,2)*(1./16.)*(VDR(I,1)+VDR(
     2    I+1,1)+VDR(I,2)+VDR(I+1,2))**2+XKMH(I,2)*0.25*(UDR(I,1)+
     3    UDR(I,2))**2)/(PN(1)*1005.0)
         TA(I,1)=TA(I,1)+DISS*TAD(I,1)
  430	CONTINUE
      RETURN    
      END       
c-----------------------------------------------------------
C
        SUBROUTINE INTERPOLATE(NLEVELS,PDS,DZ,TZ,QZ,DISS,SST,VMAX)
        PARAMETER(NA=600)
        REAL TZ(NLEVELS),QZ(NLEVELS)
        REAL PZ(NLEVELS), TAZ(NLEVELS),QAZ(NLEVELS)
        REAL Z(NA),T(NA),Q(NA),P(NA)
        OPEN(UNIT=101,FILE='s.in',STATUS='OLD')
        RD=287.
        EPSI=1./0.622
        G=9.8
        READ(101,*)
        READ(101,*)
        READ(101,*)SST
        READ(101,21)
   21   FORMAT(3X,///)
c
        DO 40 I=1,2000
          READ(101,*,END=41)P(I),T(I),Q(I)
          NP=I
   40	CONTINUE
   41   CONTINUE
        CLOSE(101)
        DO I=1,NP
         T(I)=T(I)+273.15
         Q(I)=0.001*Q(I)
        END DO
        Z(1)=0.0
        DO 50 I=2,NP
         TV1=T(I-1)*(1.+EPSI*Q(I-1)-Q(I-1))
         TV2=T(I)*(1.+EPSI*Q(I)-Q(I))
         TVBAR=0.5*(TV1+TV2)
         PBAR=0.5*(P(I-1)+P(I))
         Z(I)=Z(I-1)+(RD*TVBAR/(G*PBAR))*(P(I-1)-P(I))
   50   CONTINUE
        DO 100 I=1,NLEVELS
         AZ=0.5*DZ+DZ*FLOAT(I-1)
         DO 70 J=1,NP
         K=J
         IF(Z(J).GE.AZ)GOTO 80 
   70   CONTINUE
   80   CONTINUE
         TM1=T(K-1)*(1000./P(K-1))**(0.286)
         TM=T(K)*(1000./P(K))**(0.286)
         TZ(I)=TM1+(TM-TM1)*(AZ-Z(K-1))/(Z(K)-Z(K-1))
         QZ(I)=Q(K-1)+(Q(K)-Q(K-1))*(AZ-Z(K-1))/(Z(K)-Z(K-1))
C
         QAZ(I)=1000.*QZ(I)
         TAZ(I)=T(K-1)+(T(K)-T(K-1))*(AZ-Z(K-1))/(Z(K)-Z(K-1))-273.15
         PZ(I)=P(K-1)+(P(K)-P(K-1))*(AZ-Z(K-1))/(Z(K)-Z(K-1))
C
  100   CONTINUE
C
        CALL PCMIN(SST,PDS,PZ,TAZ,QAZ,NLEVELS,NLEVELS-1,DISS,
     1        PMIN,VMAX,IFL) 
        RETURN
        END
C
       SUBROUTINE PCMIN(SST,PSL,P,T,R,NA,N,DISS,PMIN,VMAX,IFL)
C
C   ***   This subroutine calculates the maximum wind speed        ***
C   ***             and mimimum central pressure                   ***
C   ***    achievable in tropical cyclones, given a sounding       ***
C   ***             and a sea surface temperature.                 ***
C
C  INPUT:   SST: Sea surface temperature in C
C
C           PSL: Sea level pressure (mb)
C
C           P,T,R: One-dimensional arrays of dimension NA
C             containing pressure (mb), temperature (C),
C             and mixing ratio (g/kg). The arrays MUST be
C             arranged so that the lowest index corresponds
C             to the lowest model level, with increasing index
C             corresponding to decreasing pressure. The temperature
C             sounding should extend to at least the tropopause and 
C             preferably to the lower stratosphere, however the
C             mixing ratios are not important above the boundary
C             layer. Missing mixing ratios can be replaced by zeros.
C
C           NA: The dimension of P,T and R
C
C           N:  The actual number of points in the sounding
C                (N is less than or equal to NA)
C
C  OUTPUT:  PMIN is the minimum central pressure, in mb
C
C           VMAX is the maximum surface wind speed, in m/s
C                  (reduced to reflect surface drag)
C
C           IFL is a flag: A value of 1 means OK; a value of 0
C              indicates no convergence (hypercane); a value of 2
C              means that the CAPE routine failed.
C
C-----------------------------------------------------------------------------
       REAL T(NA), P(NA), R(NA)
C
C   ***   Adjustable constant: Ratio of C_k to C_D    ***
C
       CKCD=1.0
C
C   ***   Adjustable constant for buoyancy of displaced parcels:  ***
C   ***    0=Reversible ascent;  1=Pseudo-adiabatic ascent        ***
C
      SIG=0.0
C
C   ***  Exponent, b, in assumed profile of azimuthal velocity in eye,   ***
C   ***   V=V_m(r/r_m)^b. Used only in calculation of central pressure   ***
C
       b=2.0
C
C   *** Set level from which parcels lifted   ***
C
       NK=1
C
C   *** Factor to reduce gradient wind to 10 m wind
C
       VREDUC=1.0
C
C   ***   Normalize certain quantities   ***
C
       SSTK=SST+273.15
       ES0=6.112*EXP(17.67*SST/(243.5+SST))
       DO 40 I=1,N
        R(I)=R(I)*0.001
        T(I)=T(I)+273.15
   40       CONTINUE
C
C   ***   Default values   ***
C
      VMAX=0.0
       PMIN=PSL
       IFL=1
C
       NP=0
       PM=950.0
C
C   ***   Find environmental CAPE *** 
C
      TP=T(NK)
      RP=R(NK)
      PP=P(NK)
      CALL CAPE(TP,RP,PP,T,R,P,NA,N,SIG,CAPEA,TOA,IFLAG)
      IF(IFLAG.NE.1)IFL=2
C
C   ***   Begin iteration to find mimimum pressure   ***
C
  100 CONTINUE
C
C   ***  Find CAPE at radius of maximum winds   ***
C
      TP=T(NK)
      PP=PM
      RP=0.622*R(NK)*PSL/(PM*(0.622+R(NK))-R(NK)*PSL)
      CALL CAPE(TP,RP,PP,T,R,P,NA,N,SIG,CAPEM,TOM,IFLAG) 
      IF(IFLAG.NE.1)IFL=2
      RAT=SSTK/TOM
      RAT=RAT/(DISS+RAT*(1.-DISS))
C
C  ***  Find saturation CAPE at radius of maximum winds   ***
C
      TP=SSTK
      PP=PM
      RP=0.622*ES0/(PM-ES0)
      CALL CAPE(TP,RP,PP,T,R,P,NA,N,SIG,CAPEMS,TOMS,IFLAG)
      IF(IFLAG.NE.1)IFL=2
C
C  ***  Initial estimate of minimum pressure   ***
C
      RS0=RP
      TV1=T(1)*(1.+R(1)/0.622)/(1.+R(1))
       TVAV=0.5*(TV1+SSTK*(1.+RS0/0.622)/(1.+RS0))
       CAT=CAPEM-CAPEA+0.5*CKCD*RAT*(CAPEMS-CAPEM)
C       CAT=0.5*CKCD*RAT*(CAPEMS-CAPEM)
       PNEW=PSL*EXP(-CAT/(287.04*TVAV))
C
C   ***  Test for convergence   ***
C
       IF(ABS(PNEW-PM).GT.0.2)THEN
        PM=PNEW
        NP=NP+1
        IF(NP.GT.1000.OR.PM.LT.400.0)THEN
         PMIN=400.0
         IFL=0
         GOTO 900
        END IF
        GOTO 100
       ELSE
        CATFAC=0.5*(1.+1./b)
        CAT=CAPEM-CAPEA+CKCD*RAT*CATFAC*(CAPEMS-CAPEM)
C        CAT=CKCD*RAT*CATFAC*(CAPEMS-CAPEM)
        PMIN=PSL*EXP(-CAT/(287.04*TVAV))
       END IF
  900       CONTINUE
       FAC=MAX(0.0,(CAPEMS-CAPEM))
       VMAX=VREDUC*SQRT(CKCD*RAT*FAC)
C
C   ***  Renormalize sounding arrays   ***
C       
       DO 910 I=1,N
        R(I)=R(I)*1000.0
        T(I)=T(I)-273.15
  910       CONTINUE
C
       RETURN
       END
C        
      SUBROUTINE CAPE(TP,RP,PP,T,R,P,ND,N,SIG,CAPED,TO,IFLAG)
C
C     This subroutine calculates the CAPE of a parcel with pressure PP (mb), 
C       temperature TP (K) and mixing ratio RP (gm/gm) given a sounding
C       of temperature (T in K) and mixing ratio (R in gm/gm) as a function
C       of pressure (P in mb). ND is the dimension of the arrays T,R and P,
C       while N is the actual number of points in the sounding. CAPED is
C       the calculated value of CAPE and TO is the temperature at the
C       level of neutral buoyancy.  IFLAG is a flag
C       integer. If IFLAG = 1, routine is successful; if it is 0, routine did
C       not run owing to improper sounding (e.g.no water vapor at parcel level).
C       IFLAG=2 indicates that routine did not converge.                 
C
      REAL T(ND),R(ND),P(ND),TVRDIF(100)   
      REAL NA
C
C   ***   Default values   ***
C      
      CAPED=0.0
      TO=T(1)
      IFLAG=1
C
C   ***   Check that sounding is suitable    ***
C
      IF(RP.LT.1.0E-6.OR.TP.LT.200.0)THEN
       IFLAG=0
       RETURN
      END IF            
C
C   ***   Assign values of thermodynamic constants     ***
C
      CPD=1005.7
      CPV=1870.0
C      CL=4190.0
      CL=2500.0
      CPVMCL=CPV-CL
      RV=461.5
      RD=287.04
      EPS=RD/RV
      ALV0=2.501E6
C
C   ***  Define various parcel quantities, including reversible   ***
C   ***                       entropy, S.                         ***
C                           
      TPC=TP-273.15
      ESP=6.112*EXP(17.67*TPC/(243.5+TPC))
      EVP=RP*PP/(EPS+RP)
      RH=EVP/ESP
      ALV=ALV0-CPVMCL*TPC
      S=(CPD+RP*CL)*LOG(TP)-RD*LOG(PP-EVP)+
     1   ALV*RP/TP-RP*RV*LOG(RH)            
C
C   ***  Find lifted condensation pressure, PLCL   ***
C     
       CHI=TP/(1669.0-122.0*RH-TP)
       PLCL=PP*(RH**CHI)
C
C   ***  Begin updraft loop   ***
C
       NCMAX=0
        DO J=1,N
         TVRDIF(J)=0.0
        END DO
       DO 200 J=2,N
C
C    ***   Don't bother lifting parcel above 60 mb    ***
C
        IF(P(J).LT.59.0)GOTO 200
C
C    ***  Parcel quantities below lifted condensation level   ***
C        
        IF(P(J).GE.PLCL)THEN
         TG=TP*(P(J)/PP)**(RD/CPD)
         RG=RP
C
C   ***   Calculate buoyancy   ***
C  
         TLVR=TG*(1.+RG/EPS)/(1.+RG)
         TVRDIF(J)=TLVR-T(J)*(1.+R(J)/EPS)/(1.+R(J))
        ELSE
C
C   ***  Parcel quantities above lifted condensation level  ***
C        
         TG=T(J)          
         TJC=T(J)-273.15 
         ES=6.112*EXP(17.67*TJC/(243.5+TJC)) 
         RG=EPS*ES/(P(J)-ES)
C
C   ***  Iteratively calculate lifted parcel temperature and mixing   ***
C   ***                ratio for reversible ascent                    ***
C
         NC=0
  120         CONTINUE
         NC=NC+1
C
C   ***  Calculate estimates of the rates of change of the entropy    ***
C   ***           with temperature at constant pressure               ***
C  
         ALV=ALV0-CPVMCL*(TG-273.15)
         SL=(CPD+RP*CL+ALV*ALV*RG/(RV*TG*TG))/TG
         EM=RG*P(J)/(EPS+RG)
         SG=(CPD+RP*CL)*LOG(TG)-RD*LOG(P(J)-EM)+
     1      ALV*RG/TG
         IF(NC.LT.3)THEN
          AP=0.3
         ELSE
          AP=1.0
         END IF
         TGNEW=TG+AP*(S-SG)/SL  
C
C   ***   Test for convergence   ***
C
         IF(ABS(TGNEW-TG).GT.0.01)THEN
          TG=TGNEW
          TC=TG-273.15
          ENEW=6.112*EXP(17.67*TC/(243.5+TC))
C
C   ***   Bail out if things get out of hand   ***
C
          IF(NC.GT.500.OR.ENEW.GT.(P(J)-1.0))THEN
            IFLAG=2
            RETURN
          END IF
          RG=EPS*ENEW/(P(J)-ENEW)           
          GOTO 120
         END IF
         NCMAX=MAX(NC,NCMAX)
C
C   *** Calculate buoyancy   ***
C
          RMEAN=SIG*RG+(1.-SIG)*RP
         TLVR=TG*(1.+RG/EPS)/(1.+RMEAN)
         TVRDIF(J)=TLVR-T(J)*(1.+R(J)/EPS)/(1.+R(J))
        END IF
  200       CONTINUE
C
C  ***  Begin loop to find NA, PA, and CAPE from reversible ascent ***
C
       NA=0.0
       PA=0.0
C
C   ***  Find maximum level of positive buoyancy, INB    ***
C
       INB=1
       DO 550 J=N,1,-1
        IF(TVRDIF(J).GT.0.0)INB=MAX(INB,J)
  550       CONTINUE
       IF(INB.EQ.1)RETURN
C
C   ***  Find positive and negative areas and CAPE  ***
C
       IF(INB.GT.1)THEN
        DO 600 J=2,INB
         TVM=0.5*(TVRDIF(J)+TVRDIF(J-1))
         PMA=0.5*(P(J)+P(J-1))
         IF(TVM.LE.0.0)THEN
          NA=NA-RD*TVM*(P(J-1)-P(J))/PMA
         ELSE
          PA=PA+RD*TVM*(P(J-1)-P(J))/PMA
         END IF
  600        CONTINUE
C
C   ***   Find residual positive area above INB and TO  ***
C
       PAT=0.0
       TO=T(INB)
       IF(INB.LT.N)THEN
        PINB=(P(INB+1)*TVRDIF(INB)-P(INB)*TVRDIF(INB+1))/
     1   (TVRDIF(INB)-TVRDIF(INB+1))
        PAT=RD*TVRDIF(INB)*(P(INB)-PINB)/(P(INB)+PINB)
       TO=(T(INB)*(PINB-P(INB+1))+T(INB+1)*(P(INB)-PINB))/
     1    (P(INB)-P(INB+1))
       END IF
C
C   ***   Find CAPE  ***
C            
        CAPED=PA+PAT-NA
        CAPED=MAX(CAPED,0.0)
       END IF
C
       RETURN
       END

                  


