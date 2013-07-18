*** SUPPLEMENT TO eos12.f FOR THE CASE OF ARBITRARY MAGNETIC FIELD
*
*             EQUATION OF STATE FOR ELECTRON-ION PLASMAS
*
*                 USES SUBROUTINES FROM eos12.f
*               (!!! eos12.f needs to be linked !!!)
*
* A.Y.Potekhin & G.Chabrier, in preparation (2012)
* Please send comments/suggestions to Alexander Potekhin:
*                                                   palex@astro.ioffe.ru
* Last updated 02.10.12
**   L I S T   O F   S U B R O U T I N E S (functions go with "="):
* MAIN (normally commented-out) - example driving routine
* EOSMAG - EOS of a homogeneous plasma of electrons and one ion species
* CHEMAG8 - EOS of magnetized electron ideal gas (density as input)
* ELECTMAG - EOS of magnetized electron ideal gas (chem.pot. as input)
* MAGNION - EOS of nonrelativistic nondegenerate ion gas in magn.field.
* SPINION - contributions to the EOS of nonrelativistic nondegenerate
*           ion gas due to the spin degeneracy and magnetic moments.
* EOSFIM12 - non-ideal parts of thermodynamic functions in magn.field.
* EXCORM - exchange-correlation contribution for the electron gas
* FHARMAG - harmonic (including static-lattice and zero-point)
*          contributions to the free and internal energies, pressure,
*          entropy, heat capacity, derivatives of pressure over
*          logarithms of temperature and density for solid OCP.
* HLMAG - the same as FHARM12, but only for thermal contributions
* darcth=hyperbolic arc-tangent (inverse dtanh)
************************************************************************
*                           MAIN program:               Version 26.05.12
* This driving routine allows one to compile and run this code "as is".
* In practice, however, one usually needs to link subroutines from this
* file to another (external) code, therefore the MAIN program is
* normally commented-out.
C%C      implicit double precision (A-H), double precision (O-Z)
C%C      parameter (UN_B12=425.438,UN_T6=.3157746)
C%C      data LIQSOL/0/
C%C      write(*,'('' Input Z, A: ''$)')
C%C      read*,Zion,CMI
C%C      write(*,110) Zion,CMI
C%C      write(*,'('' Input B12, T6, and RHO: ''$)')
C%C      read*,B12,T6,RHO
C%C    1 continue
C%C      write(*,111) B12,T6
C%C      TEMP=T6/UN_T6 ! T [au]
C%C      GAMAG=B12*UN_B12
C%C    2 continue
C%C      call EOSMAG(Zion,CMI,RHO,TEMP,GAMAG,
C%C     *  DENS,GAMI,CHI,TPT,LIQSOL,PnkT,UNkT,SNk,CV,CHIR,CHIT)
C%C      Tnk=8.31447d13/CMI*RHO*T6 ! n_i kT [erg/cc]
C%C      P=PnkT*Tnk/1.d12 ! P [Mbar]
C%C*   --------------------   OUTPUT   --------------------------------   *
C%C* Here in the output we have:
C%C* RHO - mass density in g/cc
C%C* P - total pressure in Mbar (i.e. in 1.e12 dyn/cm^2)
C%C* PnkT=P/nkT, where n is the number density of ions, T temperature
C%C* CV - heat capacity at constant volume, divided by number of ions, /k
C%C* CHIT - logarithmic derivative of pressure \chi_T
C%C* CHIR - logarithmic derivative of pressure \chi_\rho
C%C* UNkT - internal energy divided by NkT, N being the number of ions
C%C* SNk - entropy divided by number of ions, /k
C%C* GAMI - ionic Coulomb coupling parameter
C%C* TPT=T_p/T, where T_p is the ion plasma temperature
C%C* CHI - electron chemical potential, divided by kT
C%C* LIQSOL = 0 in the liquid state, = 1 in the solid state
C%C      write(*,111) RHO,P,PnkT,CV,CHIT,CHIR,UNkT,SNk,GAMI,TPT,CHI,
C%C     *  LIQSOL
C%C      write(*,'('' Next RHO ( < 0 to new B and T): ''$)')
C%C      read*,RHO
C%C      if (RHO.gt.0.) goto 2
C%C      write(*,'('' Next B12,T6 ( < 0 to stop): ''$)')
C%C      read*,B12,T6
C%C      if (T6.le.0..or.B12.lt.0.) stop
C%C      write(*,'('' RHO: ''$)')
C%C      read*,RHO
C%C      goto 1
C%C  110 format(' Output of eosmag v.27.05.12'/
C%C     *  '   Zion      CMI'/1P,2E11.4/
C%C     *  ' rho [g/cc]  P [Mbar]  P/(n_i kT)  Cv/N_i',
C%C     *  '     chi_T      chi_r    U/(N_i kT)    S/N_i     Gamma_i',
C%C     *  '   T_p/T   chi_e   liq/sol'/'   B12      T6')
C%C  111 format(1P,11E11.3,I2)
C%C      end

      subroutine EOSMAG(Zion,CMI,RHO,TEMP,GAMAG,
     *  DENS,GAMI,CHI,TPT,LIQSOL,PnkT,UNkT,SNk,CV,CHIR,CHIT)
!f2py   intent(in) :: Zion,CMI,RHO,TEMP,GAMAG,LIQSOL
!f2py   intent(out) ::DENS,GAMI,CHI,TPT,LIQSOL,PnkT,UNkT,SNk,CV,CHIR,CHIT
*                                                       Version 06.07.12
*                                             slightly modified 13.09.12
* EOS of fully ionized electron-ion plasma with magnetic field
* Limitations: inapplicable in the regimes of (1) bound-state formation,
*      (2) quantum liquid, (3) presence of positrons, (4) mixtures
* Input: Zion and CMI - ion charge and mass numbers,
*        RHO - total mass density [g/cc]
*        TEMP - temperature [in a.u.=2Ryd=3.158e5 K]
*        GAMAG - magnetic field [in a.u. (=B/2.3505e9 G)]
*        LIQSOL - regulator of the solid vs liquid regime (see below)
* NB: instead of RHO, a true input is CHI, determined below in CHEMAG8.
*     Hence, disagreement between RHO and DENS is the CHEMAG8 error.
* NB': Photon-gas quantities (radiation pressure etc.) are NOT INCLUDED
*     in this subroutine (if needed, add them externally).
* Output: DENS - electron number density [in a.u.=6.7483346e24 cm^{-3}]
*         GAMI - ion-ion Coulomb coupling constant
*         CHI = mu_e/kT, where mu_e is the electron chem.potential
*         TPT - ionic quantum parameter (T_p/T)
*         LIQSOL - regulator of the solid vs liquid regime (see below)
*         SNk - total dimensionless entropy per 1 ion (assumes spin=1/2)
*         UNkT - internal energy per kT per ion
*         PnkT - pressure / n_i kT, where n_i is the ion number density
*         CV - heat capacity per ion, divided by the Boltzmann constant
*         CHIR - inverse compressibility -(d ln P / d ln V)_T ("\chi_r")
*         CHIT = (d ln P / d ln T)_V ("\chi_T")
* LIQSOL is the liquid/solid switch, which works as follows:
*     if LIQSOL=0 or 1 on the input, then:
*        if (either GAMI or 1/RS is below its critical value) then
*           liquid regime and LIQSOL=0 on the output
*        otherwise solid regime and LIQSOL=1 on the output
*     if LIQSOL=2 or 3 on the input, then:
*        if (LIQSOL=2) then liquid regime and LIQSOL=2 on the output
*        if (LIQSOL=3) then solid regime and LIQSOL=3 on the output
*     if LIQSOL=4 or 5 on the input, then
*         consider non-ideal free energy FC1:
*        if (FC1=min in liquid) then liquid regime, output LIQSOL=4
*        if (FC1=min in solid) then solid regime, output LIQSOL=5
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter(TINY=1.d-7,EPS=1.d-2)
      parameter (PI=3.141592653d0,AUM=1822.888d0, 
     *   CMN=0.6, 
     *   GAMIMELT=175., 
     *   GAML=1.d3, 
     *   GAMS=1.d0, 
     *   BOHR=137.036, 
     *   RSIMELT=140.) 
      parameter (PI2=PI**2,BOHR2=BOHR**2,BOHR3=BOHR2*BOHR)
      DENS=RHO/11.20587*Zion/CMI 
      B=GAMAG/BOHR**2 
      TEMR=TEMP/BOHR2 
      Z13=Zion**.33333333d0
      DENR=DENS/BOHR3 
      DENRI=DENR/Zion 
* Plasma parameters:
      RS=(.75d0/PI/DENS)**.33333333d0 
      RSI=RS*CMI*Zion**2*Z13*AUM 
      GAME=1.d0/RS/TEMP 
      GAMI=Zion**2*GAME/Z13 
      if (LIQSOL.eq.0.or.LIQSOL.eq.1) then
        if (GAMI.lt.GAMIMELT.or.RSI.lt.RSIMELT) then
           LIQSOL=0 
        else
           LIQSOL=1 
        endif
      endif
      TPT2=GAME**2/RS*3.d0/(AUM*CMI)*Zion
* In the case of a mixture, this estimate of the quantum ion parameter
* is as crude as the linear mixing model for Fharm employed below.
      TPT=dsqrt(TPT2) ! effective T_p/T - ion quantum parameter
* (1) ideal electron gas (including relativity, degeneracy, and m.f.):
      call CHEMAG8(B,TEMR,DENR,
     *  DENR1,CHI,FEid,PEid,UEid,SEid,CVE,CHITE,CHIRE,DDENR)
* NB: CHI can be used as true input instead of RHO or DENS. In that case
*     one should define density as "DENS=DENR1*BOHR3"
* (2) ideal ion gas (including relativity, degeneracy, and m.f.):
      GYRO=0.d0
      MULTI=1
      if (dabs(Zion-1.d0).lt.TINY) then ! assume protons
         GYRO=5.5857d0 ! = g_p
         MULTI=2
      endif
      if (LIQSOL.eq.0.or.LIQSOL.eq.2.or.LIQSOL.eq.4.or.LIQSOL.eq.5) then
         call MAGNION(B,TEMR,DENRI,Zion,CMI,GYRO,MULTI,
     *     FIid,PIid,UIid,SIid,CVI,CHITI,CHIRI)
        if (LIQSOL.eq.4.or.LIQSOL.eq.5) then ! remember
           FIidL=FIid
           PIidL=PIid
           UIidL=UIid
           SIidL=SIid
           CVIL=CVI
           CHITIL=CHITI
        endif
      endif
      if (LIQSOL.eq.1.or.LIQSOL.eq.3.or.LIQSOL.eq.4.or.LIQSOL.eq.5) then
         PIid=1.d0
         CHITI=1.d0
         CHIRI=1.d0
         ZETI=Zion/(CMI*AUM)*GAMAG/TEMP ! \hbar\omega_{ci}/kT
         call SPINION(GYRO,ZETI,MULTI,Fspin,Uspin,CVspin)
         FIid=1.5d0*dlog(TPT**2/GAMI)-1.323515d0+Fspin
* 1.3235=1+0.5*ln(6/pi); FIid = F_{id.ion gas}/(N_i kT)
* NB: for a mixture, this FIid must be modified (see FidION in EOSFI7)
         UIid=1.5d0+Uspin
         CVI=1.5d0+CVspin
         SIid=UIid-FIid
      endif
* --- corr.for zero-point Landau vibrations:
      UIZP=GAMAG/(2.d0*CMI*AUM*TEMP)
* (3) Coulomb+xc nonideal contributions:
      if (LIQSOL.eq.0.or.LIQSOL.eq.2) then ! liquid
         call EOSFIM12(0,CMI,Zion,RS,GAMI,TEMP,GAMAG,1,0,
     *     FC1,UC1,PC1,SC1,CV1,PDT1,PDR1,
     *     FC2,UC2,PC2,SC2,CV2,PDT2,PDR2)
      endif
      if (LIQSOL.eq.1.or.LIQSOL.eq.3) then ! solid
         call EOSFIM12(1,CMI,Zion,RS,GAMI,TEMP,GAMAG,1,0,
     *     FC1,UC1,PC1,SC1,CV1,PDT1,PDR1,
     *     FC2,UC2,PC2,SC2,CV2,PDT2,PDR2)
      endif
      if (LIQSOL.eq.4.or.LIQSOL.eq.5) then ! choose liquid or solid
         call EOSFIM12(0,CMI,Zion,RS,GAMI,TEMP,GAMAG,1,0,
     *     FC1L,UC1L,PC1L,SC1L,CV1L,PDT1L,PDR1L,
     *     FC2L,UC2L,PC2L,SC2L,CV2L,PDT2L,PDR2L)
         call EOSFIM12(1,CMI,Zion,RS,GAMI,TEMP,GAMAG,1,0,
     *     FC1,UC1,PC1,SC1,CV1,PDT1,PDR1,
     *     FC2,UC2,PC2,SC2,CV2,PDT2,PDR2)
        if (FC1L.lt.FC1) then ! switch to the liquid; recall it.
           LIQSOL=4
           FIid=FIidL
           PIid=PIidL
           UIid=UIidL
           SIid=SIidL
           CVI=CVIL
           CHITI=CHITIL
           FC1=FC1L
           UC1=UC1L
           PC1=PC1L
           SC1=SC1L
           CV1=CV1L
           PDT1=PDT1L
           PDR1=PDR1L
           FC2=FC2L
           UC2=UC2L
           PC2=PC2L
           SC2=SC2L
           CV2=CV2L
           PDT2=PDT2L
           PDR2=PDR2L
        else
           LIQSOL=5
        endif
      endif
* Calculate partial thermodynamic quantities and combine them together:
** Normalization factors:
      DTE=DENS*TEMP ! electron normalization [au]
      DENSI=DENS/Zion ! number density of all ions [au]
      DTI=DENSI*TEMP ! ion normalization [au]
** Electrons:
*** First-order TD functions:
      UINTE=UEid*DTE
      PRESSE=PEid*DTE ! P_e [a.u.]
      StotE=SEid*DENS ! electron entropy per unit volume [a.u.]
*** 2nd-order TD functions:
      CVtotE=CVE*DENS
      PDLTE=PRESSE*CHITE
      PDLRE=PRESSE*CHIRE
** Ions:
      if (LIQSOL.eq.0.or.LIQSOL.eq.2.or.LIQSOL.eq.4) then ! liquid
*** First-order TD functions:
         UINTI=DTI*(UIid+UC1) ! (i+ii+ie+ee)
         PRESSI=DTI*(PIid+PC1)
         StotI=DENSI*(SIid+SC1)
*** 2nd-order TD functions:
         CVtotI=DENSI*(CVI+CV1)
         PDLTI=DTI*(PIid*CHITI+PDT1)
         PDLRI=DTI*(PIid*CHIRI+PDR1)
      else ! solid
*** First-order TD functions:
         UINTI=DTI*UC2
         PRESSI=DTI*PC2
         StotI=DENSI*SC2
*** 2nd-order TD functions:
         CVtotI=DENSI*CV2
         PDLTI=DTI*PDT2
         PDLRI=DTI*PDR2
      endif
** Final sum-up:
      UINT=UINTE+UINTI ! int.en. per unit V (e+i+Coul.) [a.u.]
      PRESS=PRESSE+PRESSI ! pressure (e+i+Coul.) [a.u.]
      Stot=StotE+StotI
      CVtot=CVtotE+CVtotI
      PDLT=PDLTE+PDLTI
      PDLR=PDLRE+PDLRI
* Dimensionless ratios:
** First-order:
      PnkT=PRESS/DTI ! beta P / n_i
      UNkT=UINT/DTI ! beta U / N_i
      SNk=Stot/DENSI ! S / N_i k_B
** Second-order:
      CV=CVtot/DENSI ! C_V per ion
      CHIR=PDLR/PRESS ! d ln P / d ln\rho
      CHIT=PDLT/PRESS ! d ln P / d ln T
      return
      end

*   =========  CHEMICAL POTENTIAL and DENSITY DERIVATIVE  ==========   *
      subroutine CHEMAG8(B,TEMR,DENR,
     *  DENR1,CHI,FEid,PEid,UEid,SEid,CVE,CHITE,CHIRE,DDENR)
      implicit double precision (A-H), double precision (O-Z)
      save
*                                                       Version 04.03.09
*                                               Cosmetic change 02.09.11
*                                                EPS diminished 26.05.12
* Stems from CHEMPOT v.03.05.07
* Ideal electron-gas EOS
* Input: B - magnetic field [rel.un.] (=B/4.414e13 G)
*        TEMR - temperature [rel.un.] (=T/5.93e9 K)
*        DENR - electron number density n_e [rel.un.]
* Output: DENR1 - the density adjusted to the found CHI
*         CHI=\mu/T, where \mu is the chem.pot-l w/o the rest-energy
*               (= n_e/1.7366e31 cm^{-3})
*         FEid - free energy / N_e kT,
*         PEid - pressure (P_e) / n_e kT
*         UEid - internal energy / N_e kT
*         SEid - entropy, CVE - heat capacity [per 1 electron, in k_B],
*         CHITE=(d ln P_e / d ln T)_V, CHIRE=(d ln P_e / d ln n_e)_T
*         DDENR = (d n / d\mu)_T (rel.un.)
*  All quantities are by default in relativistic units
*  B - magnetic field (=h\omega_B/mc^2)
      parameter (BOHR=137.036,PI=3.141592653d0)
      parameter (PI2=PI**2,BOHR2=BOHR**2,BOHR3=BOHR2*BOHR)
      parameter (NLMAX1=50,PRECLN=9.,EPS=1.D-4,NI=50) 
      BT=B/TEMR
      if (BT.lt.10.) then 
         NLMAX=NLMAX1
      else
         NLMAX=NLMAX1/dsqrt(BT/10.)
         NLMAX=max0(10,NLMAX)
      endif
      CNLMAX=NLMAX
      PF0=(29.6088132d0*DENR)**.33333333d0 
      if (B.gt.0.) then
         PFB=PF0**3/(1.5*B)
         PF=dmin1(PF0,PFB) 
        if (PF.gt.1.d-7) then
           TF=dsqrt(1.d0+PF**2)-1.d0 
        else
           TF=0.5d0*PF**2
        endif
         DE=TF+PRECLN*TEMR 
         CNL0=DE*(1.d0+0.5d0*DE)/B 
      else
         CNL0=NLMAX+1
      endif
* Preliminary nonmagnetic calculation
      call CHEMFIT7(DENR,TEMR,CHI,CMU1,1,
     *  CMUDENR,CMUDT,CMUDTT)
      DDENR=1.d0/CMUDENR
      TEMP=TEMR*BOHR2
      call ELECT11(TEMP,CHI,
     *  DENSNM,FEidNM,PEidNM,UEidNM,SEidNM,CVENM,CHITENM,CHIRENM,DDNNM,
     *  DlnDT,DlnDHH,DlnDTT,DlnDHT)
      DENR1NM=DENSNM/BOHR3
      if (CNL0.gt.CNLMAX) then ! Nonmagnetic case
         DENR1=DENR1NM
         FEid=FEidNM
         PEid=PEidNM
         UEid=UEidNM
         SEid=SEidNM
         CVE=CVENM
         CHITE=CHITENM
         CHIRE=CHIRENM
         DDN=DDNNM
        return
      else ! Magnetic calculation:
         ZET=(dmin1(1.d0,PF0**2/B/1.5)*PF0)**2
         NL0=CNL0+.5
      endif
*   ----------  Calculation of chemical potential CMU   ------------   *'
      Kup=0.
      Klow=0.
      POWER=.666666667d0
      if (NL0.lt.1) POWER=2.d0
      DENR1=0.
      do 5 I=1,NI
      if (ZET.gt.3.d0*TEMR) then ! Degenerate plasma
         SHIFTMU=0.
      else ! Non-degenerate plasma
         SHIFTMU=1.5d0*TEMR*dlog(3.d0*TEMR/ZET)
         if (NL0.lt.2) SHIFTMU=SHIFTMU/3.d0
      endif
      if (SHIFTMU.lt..5d0) then
         CMU=sqrt(1.+ZET)*(1.-SHIFTMU)
      else
         CMU=sqrt(1.+ZET)-SHIFTMU
      endif
      CHI=(CMU-1.d0)/TEMR
      DENR1old=DENR1
      call ELECTMAG(B,TEMR,CHI,
     *  DENR1,FEid,PEid,UEid,SEid,CVE,CHITE,CHIRE,DDENR)
      if (DENR1.gt.DENR) then
        if (Kup.eq.0) Zup=ZET
         Kup=1
         Zup=dmin1(Zup,ZET)
      else
        if (Klow.eq.0) Zlow=ZET
         Klow=1
         Zlow=dmax1(Zlow,ZET)
      endif
      ZET1=ZET
      DENRMIN=dmin1(1.d-19,DENR/1.d2)
      DENR1=dmax1(DENRMIN,DENR1)
      ZET=ZET1*(DENR/DENR1)**POWER
      if (Kup.ne.0.and.(ZET.gt.Zup.or.(I.gt.NI/2.and.DENR1.lt.DENR)))
     &   ZET=(ZET1+Zup)/2.d0
      if (Klow.ne.0.and.(ZET.lt.Zlow.or.(I.gt.NI/2.and.DENR1.gt.DENR)))
     &   ZET=(ZET1+Zlow)/2.d0
      if (abs(DENR/DENR1-1.).lt.EPS) goto 6
      if (dabs(DENR1-DENR1old).lt.EPS*.1d0*DENR1.
     .   and.DENR1.gt.DENRMIN*(1.d0+EPS)) goto 55
    5 continue
   55 continue
    6 continue
* Interpolation in order to exclude jump at nonmagnetic/magnetic switch
      Y=(1.d0-CNL0/CNLMAX)**3
      X=1.d0-Y
      DENR1=X*DENR1NM+Y*DENR1
      FEid=X*FEidNM+Y*FEid
      PEid=X*PEidNM+Y*PEid
      UEid=X*UEidNM+Y*UEid
      SEid=X*SEidNM+Y*SEid
      CVE=X*CVENM+Y*CVE
      CHITE=X*CHITENM+Y*CHITE
      CHIRE=X*CHIRENM+Y*CHIRE
      return
      end

* ==============  IDEAL ELECTRONS  =================================== *
      subroutine ELECTMAG(B,TEMR,CHI,
     *  DENR,FEid,PEid,UEid,SEid,CVE,CHITE,CHIRE,DDENR)
*                                                       Version 04.03.09
* BLIN8 --> BLIN9                                               24.11.10
* Ideal magnetized electron-gas EOS.
* Input: B - magnetic field [rel.un.] (=B/4.414e13 G)
*        TEMR - temprature [rel.un.] (=T/5.93e9 K)
*        CHI=\mu/T, where \mu is the chem.pot-l w/o the rest-energy
* Output: DENR - electron number density n_e [rel.un.]
*               (= n_e/1.7366e31 cm^{-3})
*         FEid - free energy / N_e kT,
*         PEid - pressure (P_e) / n_e kT
*         UEid - internal energy / N_e kT
*         SEid - entropy, CVE - heat capacity [per 1 electron, in k_B],
*         CHITE=(d ln P_e / d ln T)_V, CHIRE=(d ln P_e / d ln n_e)_T
*         DDENR=(d n_e / d \mu)_T
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter (BOHR=137.036,PI=3.141592653d0)
      parameter (PI2=PI**2,BOHR2=BOHR**2,BOHR3=BOHR2*BOHR)
      parameter (NLMAX=200,PRECLN=9.,EPST=1.d-17,EPS14=3.d0,EPSR=.1d0)
      if (B.gt.0.) then
         DE=(dmax1(0.d0,CHI)+PRECLN)*TEMR ! (1+DE)=Upper E bound to int.
         NL0=DE*(1.+DE/2.)/B ! Max contributing Landau number to integr.
      else
         NL0=NLMAX+1
      endif
      if (NL0.gt.NLMAX) then ! nonmagnetic calc.
         TEMP=TEMR*BOHR2
         call ELECT11(TEMP,CHI,
     *     DENS,FEid,PEid,UEid,SEid,CVE,CHITE,CHIRE,DDN,
     *     DlnDT,DlnDHH,DlnDTT,DlnDHT)
         DENR=DENS/BOHR3
         DDENR=DENR*DDN/TEMR
        return
      endif
      if (CHI.lt.-40.d0) then ! nondegenerate
         CLE=dsqrt(2.d0*PI/TEMR) ! thermal length
         DENR=2.d0*dexp(CHI)/CLE**3
         call MAGNION(B,TEMR,DENR,1.d0,5.4858d-4,2.d0,2,
     *     FIid,PIid,UIid,SIid,CVI,CHITI,CHIRI)
         DDENR=DENR/TEMR
        return
      endif
* Magnetic calculation:
      SP=0.
      SN=0.
      SNDH=0.
      SNDT=0.
      SPDH=0.
      SPDT=0.
      SNDHH=0.
      SNDTT=0.
      SNDHT=0.
      SPDHH=0.
      SPDTT=0.
      SPDHT=0.
      do NL=0,NL0
         EN1=B*NL
         EN0=dsqrt(1.d0+2.d0*EN1)
        if (EN1.gt.1.d-7) EN1=EN0-1.d0
         CHIN=CHI-EN1/TEMR
         TN=TEMR/EN0
         SQREN=dsqrt(EN0)
         call BLIN9(TN,CHIN,
     *      W0,W0DX,W0DT,W0DXX,W0DTT,W0DXT,
     *      W1,W1DX,W1DT,W1DXX,W1DTT,W1DXT,
     *      W2,W2DX,W2DT,W2DXX,W2DTT,W2DXT,
     *      W0XXX,W0XTT,W0XXT)
         H14=EPS14
        if (CHIN.lt.14.d0.and.14.d0-CHIN.lt.H14) then ! regularization
           CHINM=14.d0-H14
           CHINP=14.d0
           XM=(CHIN-CHINM)/H14
           XP=1.d0-XM
           call BLIN9(TN,CHINM,
     *        W0M,W0DXM,W0DTM,W0DXXM,W0DTTM,W0DXTM,
     *        W1M,W1DXM,W1DTM,W1DXXM,W1DTTM,W1DXTM,
     *        W2M,W2DXM,W2DTM,W2DXXM,W2DTTM,W2DXTM,
     *        W0XXXM,W0XTTM,W0XXTM)
           call BLIN9(TN,CHINP,
     *        W0P,W0DXP,W0DTP,W0DXXP,W0DTTP,W0DXTP,
     *        W1P,W1DXP,W1DTP,W1DXXP,W1DTTP,W1DXTP,
     *        W2P,W2DXP,W2DTP,W2DXXP,W2DTTP,W2DXTP,
     *        W0XXXP,W0XTTP,W0XXTP)
           W0=XP*W0+XM*W0P
           W0DX=XP*W0DX+XM*W0DXP
           W0DT=XP*W0DT+XM*W0DTP
           W0DXX=XP*W0DXX+XM*W0DXXP
           W0DXT=XP*W0DXT+XM*W0DXTP
           W0DTT=XP*W0DTT+XM*W0DTTP
           W0XXX=XP*W0XXX+XM*W0XXXP
           W0XTT=XP*W0XTT+XM*W0XTTP
           W0XXT=XP*W0XXT+XM*W0XXTP
           W1=XP*W1+XM*W1P
           W1DX=XP*W1DX+XM*W1DXP
           W1DT=XP*W1DT+XM*W1DTP
           W1DXX=XP*W1DXX+XM*W1DXXP
           W1DXT=XP*W1DXT+XM*W1DXTP
           W1DTT=XP*W1DTT+XM*W1DTTP
           W2=XP*W2+XM*W2P
           W2DX=XP*W2DX+XM*W2DXP
           W2DT=XP*W2DT+XM*W2DTP
           W2DXX=XP*W2DXX+XM*W2DXXP
           W2DXT=XP*W2DXT+XM*W2DXTP
           W2DTT=XP*W2DTT+XM*W2DTTP
        endif
         SN=SN+W0DX*SQREN
         SP=SP+SQREN*W0
         SNDH=SNDH+W0DXX*SQREN
         SPDH=SPDH+SQREN*W0DX
         SNDHH=SNDHH+W0XXX*SQREN
         SPDHH=SPDHH+SQREN*W0DXX
         CHINDT=EN1/TEMR**2
         CHINDTT=-2.d0*EN1/TEMR**3
         SNDT=SNDT+(W0DX/(2.d0*TEMR)+W0DXX*CHINDT+W0DXT/EN0)*SQREN
         SPDT=SPDT+SQREN*(1.5d0*W0/TEMR+CHINDT*W0DX+W0DT/EN0)
         SNDTT=SNDTT+(-.25d0*W0DX/TEMR**2-CHINDT/TEMR*W0DXX+
     +     W0DXT/(TEMR*EN0)+CHINDT**2*W0XXX+2.d0*CHINDT/EN0*W0XXT+
     +     W0XTT/EN0**2)*SQREN
         SPDTT=SPDTT+SQREN*(.75d0/TEMR**2*W0+CHINDT/TEMR*W0DX+
     +     3.d0*W0DT/(EN0*TEMR)+CHINDT**2*W0DXX+2.d0*CHINDT/EN0*W0DXT+
     +     W0DTT/EN0**2)
         SNDHT=SNDHT+(.5d0*W0DXX/TEMR+CHINDT*W0XXX+W0XXT/EN0)*SQREN
         SPDHT=SPDHT+SQREN*(1.5d0*W0DX/TEMR+CHINDT*W0DXX+W0DXT/EN0)
         SCHT1=SPDT-SPDH*SNDT/SNDH
        if (NL.eq.0) then
           SN=SN/2.d0
           SP=SP/2.d0
           SNDH=SNDH/2.d0
           SPDH=SPDH/2.d0
           SNDT=SNDT/2.d0
           SPDT=SPDT/2.d0
           SNDHH=SNDHH/2.d0
           SPDHH=SPDHH/2.d0
           SNDTT=SNDTT/2.d0
           SPDTT=SPDTT/2.d0
           SNDHT=SNDHT/2.d0
           SPDHT=SPDHT/2.d0
           SCHT=SCHT1/2.d0
        else
           if(SCHT1.ge.SCHT.or.NL.lt.NL0-1) SCHT=SCHT1
        endif
      enddo ! next NL
    1 continue
      CN=dsqrt(2.d0*TEMR)/PI2*B
      DENR=CN*SN ! n_e [rel.un.]
      PRI=DENR*TEMR ! ideal-gas pressure = nkT = normalization factor
      CNP=CN*TEMR
      PR=CNP*SP ! P [rel.un.]
      PEid=PR/PRI ! P/nkT
      FR=CHI*PRI-PR ! F/V [rel.un.]
      FEid=FR/PRI ! F/NkT
* derivatives over chi at constant T:
      dndH=CN*SNDH ! (d n_e/d\chi)_T
      dPdH=CNP*SPDH ! (d P_e/d\chi)_T
      dndHH=CN*SNDHH ! (d^2 n_e/d\chi^2)_T
      dPdHH=CNP*SPDHH ! (d^2 P_e/d\chi^2)_T
      dFdH=(DENR+CHI*dndH)*TEMR-dPdH ! [d (F/V)/d\chi]_T
      dFdHH=(2.d0*dndH+CHI*dndHH)*TEMR-dPdHH ! [d^2 (F/V)/d\chi^2]_T
* derivatives over T at constant chi:
      dndT=CN*SNDT ! (d n_e/dT)_\chi
      dPdT=CNP*SPDT ! (dP/dT)_\chi
      dFdT=CHI*(DENR+TEMR*dndT)-dPdT
      dndTT=CN*SNDTT
      dPdTT=CNP*SPDTT
      dFdTT=chi*(2.d0*dndT+TEMR*dndTT)-dPdTT
* mixed derivatives
      dndHT=CN*SNDHT
      dPdHT=CNP*SPDHT
      dFdHT=DENR+CHI*dndH+TEMR*dndT+CHI*TEMR*dndHT-dPdHT
* normalized thermodynamic functions:
      G=dndT/dndH ! = - (d\chi/dT)_V
      GDT=(dndTT-dndT*dndHT/dndH)/dndH
      GDH=(dndHT-dndT*dndHH/dndH)/dndH
      S=dFdH*G-dFdT ! S=entropy per volume=-(dF/dT)_V
      SEid=S/DENR ! S/Nk
      UEid=FEid+SEid
      SDT=dFdHT*G+dFdH*GDT-dFdTT
      SDH=dFdHH*G+dFdH*GDH-dFdHT
      CV=TEMR*(SDT-SDH*G) ! C_V per volume = T(dS/dT)_V
      CVE=CV/DENR ! C_V/Nk
C      CHITE1=TEMR/PR*(dPdT-dPdH*G) ! d ln P / d ln T |_V
      CHITE=TEMR/PR*CNP*SCHT ! d ln P / d ln T |_V
      if (dabs(dPdT-dPdH*G).lt.EPST*dPdT) then ! accuracy loss; non-mag.
         TEMP=TEMR*BOHR2
         call ELECT11(TEMP,CHI,
     *     DENS,FEnm,PEnm,UEnm,SEnm,CVEnm,CHITE,CHIREnm,DDNnm,
     *     DlnDT,DlnDHH,DlnDTT,DlnDHT)
      endif
      CHIRE=DENR/PR*dPdH/dndH ! d ln P / d ln rho |_T
      DDENR=dndH/TEMR ! (d n / d\mu)_T (rel.un.)
      return
      end

* =============================   IONS   ============================= *
      subroutine MAGNION(B,TEMR,DENRI,Zion,CMI,GYRO,MULTI,
     *   FIid,PIid,UIid,SIid,CVI,CHITI,CHIRI)
*                                                       Version 02.10.12
* Previous version was 04.03.09. The difference: SPINION is separate now
* Thermodynamic functions of an ideal magnetized gas of atomic nuclei
*   (see the "EOS" written notes, page 129)
* Input: B - magnetic field [rel.un.] (=B/4.414e13 G)
*        TEMR - temperature [rel.un.] (=T/5.93e9 K)
*        DENRI - number density of ions [rel.un.]
*        Zion,CMI - atomic charge and mass numbers of the ions
*        GYRO - g-factor
*        MULTI - spin multiplicity (2S+1)
* Output: FIid - free energy per ion per kT,
*         PIid=1 - pressure /(n_i kT)
*         UIid - internal energy per ion per kT,
*         SIid - entropy /k, per ion,
*         CVI - specific heat /k, per ion,
*         CHITI=CHIRI=1 - log.derivatives of pressure.
* NB: zero-point term is not included.
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter (AUM=1822.888d0,PI=3.141592653d0,EPS=1.d-8)
      parameter (TwoPI=2.d0*PI)
      RMI=CMI*AUM ! nucleus-to-electron mass ratio
      CLI=dsqrt(TwoPI/(RMI*TEMR)) ! thermal length
      ZETI=Zion/RMI*B/TEMR ! \hbar\omega_{ci}/kT
      if (ZETI.lt.EPS) then ! non-magnetic case
         FIid=dlog(DENRI*CLI**3)-1.d0
         UIid=1.5d0         
         CVI=1.5d0
      else ! magnetic case
         EZ1=dexp(-ZETI)
         EZ=1.d0-EZ1
         FIid=dlog(TwoPI/B*CLI*DENRI*EZ/Zion)-1.d0+.5d0*ZETI
         ZETI1=ZETI*EZ1/EZ
         UIid=.5d0+ZETI1+.5d0*ZETI
         CVI=.5d0+ZETI1**2
      endif
      PIid=1.d0
      CHITI=1.d0
      CHIRI=1.d0
      call SPINION(GYRO,ZETI,MULTI,Fspin,Uspin,CVspin)
      FIid=FIid+Fspin
      UIid=UIid+Uspin
      CVI=CVI+CVspin
      SIid=UIid-FIid
      return
      end

      subroutine EOSFIM12(LIQSOL,CMI,Zion,RS,GAMI,TEMP,GAMAG,KWK,KDIR,
     *  FC1,UC1,PC1,SC1,CV1,PDT1,PDR1,
     *  FC2,UC2,PC2,SC2,CV2,PDT2,PDR2)
*                                                       Version 03.06.12
* Stems from EOSFIM9 v.17.01.11
* New (as of 2012): the magnetic effects on phonon gas are included,
*          spin multiplicity in solid is made consistent with the liquid
* Non-ideal parts of thermodynamic functions in the fully ionized plasma
* Input: LIQSOL=0/1(liquid/solid), 
*        Zion,CMI - ion charge and mass numbers,
*        RS=r_s (electronic density parameter),
*        GAMI=Gamma_i (ion coupling),
*        TEMP - temperature in a.u.
*        GAMAG=B/2.35e9 G (magn.field in a.u.),
*        KWK=1 switches quantum diffraction for ion liquid on,
*        KDIR=1 switches an alternative direction of the crystal on.
* Output: FC1 and UC1 - non-ideal "ii+ie+ee" contribution to the 
*         free and internal energies (per ion per kT),
*         PC1 - analogous contribution to pressure divided by (n_i kT),
*         CV1 - "ii+ie+ee" heat capacity per ion [units of k]
*         PDT1=(1/N_i kT)*(d P_C/d ln T)_V
*         PDR1=(1/N_i kT)*(d P_C/d ln\rho)_T
* FC2,UC2,PC2,SC2,CV2 - analogous to FC1,UC1,PC1,SC1,CV1, but including
* the part corresponding to the ideal ion gas WITHOUT magnetic field. 
* This is useful for preventing accuracy loss in some cases
*   (e.g., when SC2 << SC1).
* IDEAL-GAS MAGNETIC-FIELD CONTRIBUTION FOR THE IONS SHOULD BE ADDED
*   EXTERNALLY.
* FC2 does not take into account the ENTROPY OF MIXING S_{mix}: in a
*   mixture, S_{mix}/(N_i k) HAS TO BE ADDED EXTERNALLY.
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter(TINY=1.d-20) 
      parameter (AUM=1822.888d0) 
      GAME=GAMI/Zion**1.666667
      call EXCORM(RS,GAME,GAMAG,
     *  FXC,UXC,PXC,CVXC,SXC,PDTXC,PDRXC) 
* Calculate "ii" part:
      COTPT=dsqrt(3./AUM/CMI)/Zion**(7./6.) ! auxiliary coefficient
      TPT=GAMI/dsqrt(RS)*COTPT              ! T_p/T
      FidION=1.5*dlog(TPT**2/GAMI)-1.323515
* 1.3235=1+0.5*ln(6/pi); FidION = F_{id.ion gas}/(N_i kT), but without
* the term x_i ln x_i = -S_{mix}/(N_i k).
      FWK=0.
      UWK=0.
      PWK=0.
      CVWK=0.
      PDRWK=0.
      PDTWK=0.
      if (LIQSOL.eq.0) then                 ! liquid
         call FITION9(GAMI,FION,UION,PION,CVii,PDTii,PDRii)
         FItot=FION+FidION
         UItot=UION+1.5d0
         PItot=PION+1.d0
         CVItot=CVii+1.5d0
         SCItot=UItot-FItot
         PDTi=PDTii+1.d0
         PDRi=PDRii+1.d0
* Wigner-Kirkwood (quantum diffraction) terms:
        if (KWK.eq.1) then
           ZETI=Zion/(CMI*AUM)*GAMAG/TEMP ! \hbar\omega_{ci}/kT
           call QDMAG(TPT,ZETI,FWK,UWK,PWK,CVWK,PDTWK,PDRWK)
        endif
      else ! solid
         BP=Zion/(CMI*AUM)*GAMAG/TEMP/TPT ! \hbar\omega_{ci}/kT_p
         call FHARMAG(GAMI,TPT,BP,
     *     Fharm,Uharm,Pharm,CVharm,Sharm,PDTharm,PDRharm,EVAR)
        if (KDIR.eq.1) then ! less favorable orientation of the crystal
           Fharm=Fharm+EVAR
           Uharm=Uharm+EVAR
        endif
         call ANHARM8(GAMI,TPT,Fah,Uah,Pah,CVah,PDTah,PDRah) 
         FItot=Fharm+Fah ! spin term -.6931=-ln2 is excluded 28.04.12
         FION=FItot-FidION
         UItot=Uharm+Uah
         UION=UItot-1.5d0 ! minus 1.5=ideal-gas, in order to get "ii"
         PItot=Pharm+Pah
         PION=PItot-1.d0 ! minus 1=ideal-gas
         PDTi=PDTharm+PDTah
         PDRi=PDRharm+PDRah
         PDTii=PDTi-1.d0 ! minus 1=ideal-gas
         PDRii=PDRi-1.d0 ! minus 1=ideal-gas
         CVItot=CVharm+CVah
         SCItot=Sharm+Uah-Fah ! spin term .6931=ln2 is excluded 28.04.12
         CVii=CVItot-1.5d0 ! minus 1.5=ideal-gas
      endif
* Calculate "ie" part: NB: here the ie-scaling is abandoned!
      if (LIQSOL.eq.1) then
         call FSCRsol8(RS,GAMI,Zion,TPT,
     *     FSCR,USCR,PSCR,S_SCR,CVSCR,PDTSCR,PDRSCR)
      else
         call FSCRliq8(RS,GAME,Zion,
     *     FSCR,USCR,PSCR,CVSCR,PDTSCR,PDRSCR)
         S_SCR=USCR-FSCR
      endif
* Total excess quantities ("ii"+"ie"+"ee", per ion):
      FC0=FSCR+Zion*FXC+FWK
      UC0=USCR+Zion*UXC+UWK
      PC0=PSCR+Zion*PXC+PWK
      SC0=S_SCR+Zion*SXC+FWK
      CV0=CVSCR+Zion*CVXC+CVWK
      PDT0=PDTSCR+Zion*PDTXC+PDTWK
      PDR0=PDRSCR+Zion*PDRXC+PDRWK
      FC1=FION+FC0
      UC1=UION+UC0
      PC1=PION+PC0
      SC1=(UION-FION)+SC0
      CV1=CVii+CV0
      PDT1=PDTii+PDT0
      PDR1=PDRii+PDR0
* Total excess + ideal-ion quantities
      FC2=FItot+FC0
      UC2=UItot+UC0
      PC2=PItot+PC0
      SC2=SCItot+SC0
      CV2=CVItot+CV0
      PDT2=PDTi+PDT0
      PDR2=PDRi+PDR0
      return
      end

* ==============  ELECTRON EXCHANGE AND CORRELATION   ================ *
      subroutine EXCORM(RS,GAME,GAMAG,
     *  FXC,UXC,PXC,CVXC,SXC,PDTXC,PDRXC) 
*                                                       Version 06.05.08
*                          "d0" is added at some exact numerals 27.04.12
* Exchange-correlation contribution for the electron gas
* Stems from EXCOR7 v.09.06.07 and TANAKMAG v.27.04.98
* Input: RS - electron density parameter =electron-sphere radius in a.u.
*        GAME - electron Coulomb coupling parameter
*        GAMAG - atomic magnetic-field parameter
* Output: FXC - excess free energy of e-liquid per kT per one electron
*             according to Tanaka & Ichimaru 85-87 and Ichimaru 93 FXC,
*         UXC - free and internal energy contr.[per 1 electron, kT]
*         PXC - pressure contribution divided by (n_e kT)
*         CVXC - heat capacity divided by N_e k
*         SXC - entropy divided by N_e k
*         PDTXC,PDRXC = PXC + d ln PXC / d ln(T,\rho)
      implicit double precision(A-H),double precision(O-Z)
      save
      parameter(TINY=1.d-99)
      if (RS.le.0.) stop'wrong RS'
      if (GAME.le.0.) stop'wrong GAME'
      THETA1=.543*RS/GAME ! non-relativistic non-mag.degeneracy par.
*  ------   06.04.98:  magnetic scale factor on RS at high T   ------  *
      U=GAMAG/2.d0*GAME*RS ! hbar omega_c/(2T)...
      SCALE=1.d0
      SCADU=0.
      SCADUU=0.
      if (U.gt.1.d-4) then
         THU=dtanh(U)
         TU=THU/U
         S=dsqrt(1.d0-TU)
         CH2=2.d0
         TDU=-TU/U
         CDU=0.
         CDUU=0.
        if (U.lt.20.) then
           CHU2=dcosh(U)**2
           CH2=dcosh(2.d0*U)/CHU2
           VU=darcth(S)/S
           TDU=1.d0/(U*CHU2)+TDU
           TDUU=2.d0/U**2*TU-2./(U**2*CHU2)-2.d0*TU/CHU2
           CDU=2.d0*THU/CHU2
           CDUU=(6.d0/CHU2-4.d0)/CHU2
        else
           TDUU=2.d0/U**2*TU
           VU=dlog(4.d0*U-1.d0)/(2.d0-1.d0/U)
        endif
         VDU=.5d0*TDU/S**2*(VU-1.d0/TU)
         VDUU=VDU*TDUU/TDU+(VDU*TDU*1.5d0+(TDU/TU)**2*.5d0)/S**2
         SCALE=CH2*TU*VU ! =f_2 of Steinberg=f_1 of PCS99
         SCADU=CDU*TU*VU+CH2*TDU*VU+CH2*TU*VDU
         SCADUU=CDUU*TU*VU+CH2*TDUU*VU+CH2*TU*VDUU+
     +     2.d0*(CDU*TDU*VU+CDU*TU*VDU+CH2*TDU*VDU)
      endif
      UDG=GAMAG*THETA1*GAME/.543
      SCADH=SCADU*GAMAG/1.086*GAME**2
      SCADG=SCADU*UDG
      SCADHH=SCADUU*(GAMAG*GAME**2/1.086)**2
      SCADGG=SCADUU*UDG**2+SCADU*GAMAG*THETA1/.543
      SCADHG=GAMAG*GAME/.543*(SCADU+SCADUU*GAME/2.d0*UDG)
*  ------ Finished calculation of SCALE and its derivatives.
      GR4=(GAMAG*RS**2)**2+TINY
*   ----------------------   scaled theta and its derivatives:
      THETAM=.090063*GR4*RS/GAME ! magn.theta; .09=8/(9\pi^2)
      TMT=.16586*GR4 ! theta_m/theta_0
      TMTDH=4.d0*TMT/THETA1
      TMTDG=4.d0*TMT/GAME
      TMTDHH=12.d0*TMT/THETA1**2
      TMTDGG=12.d0*TMT/GAME**2
      TMTDHG=16.d0*TMT/(THETA1*GAME)
      SCUP=1.d0+TMT
      E=dexp(-1.d0/THETAM)
      EDH=5.d0*E/(THETA1*THETAM)
      EDG=4.d0*E/(THETAM*GAME)
      EDHH=EDH/THETA1*(5.d0/THETAM-6.d0)
      EDGG=EDG/GAME*(4.d0/THETAM-5.d0)
      EDHG=EDH/GAME*(4.d0/THETAM-4.d0)
      PM=SCALE*TMT*E
      PMDH=TMTDH*SCALE*E+TMT*E*SCADH+TMT*EDH*SCALE
      PMDG=TMTDG*SCALE*E+TMT*E*SCADG+TMT*EDG*SCALE
      PMDHH=TMTDHH*SCALE*E+TMT*E*SCADHH+TMT*EDHH*SCALE+
     +  2.d0*(TMTDH*SCADH*E+TMTDH*SCALE*EDH+TMT*SCADH*EDH)
      PMDGG=TMTDGG*SCALE*E+TMT*E*SCADGG+TMT*EDGG*SCALE+
     +  2.d0*(TMTDG*SCADG*E+TMTDG*SCALE*EDG+TMT*SCADG*EDG)
      PMDHG=TMTDHG*SCALE*E+TMT*E*SCADHG+TMT*EDHG*SCALE+
     +  E*(TMTDH*SCADG+TMTDG*SCADH)+SCALE*(TMTDH*EDG+TMTDG*EDH)+
     +  TMT*(SCADH*EDG+SCADG*EDH)
      SCDN=1.d0+PM
      THETA=THETA1*SCUP/SCDN ! SCALED degeneracy parameter
      THDH=(SCUP+THETA1*TMTDH-THETA*PMDH)/SCDN
      THDG=(THETA1*TMTDG-THETA*PMDG)/SCDN
      THDHH=(2.d0*(TMTDH+(THETA*PMDH-SCUP)/SCDN*PMDH)+
     +  THETA1*(TMTDHH-2.d0*TMTDH*PMDH/SCDN)-THETA*PMDHH)/SCDN
      THDGG=(THETA1*TMTDGG+
     +  2.d0*(THETA*PMDG-THETA1*TMTDG)*PMDG/SCDN-THETA*PMDGG)/SCDN
      THDHG=(TMTDG+THETA1*TMTDHG-THETA*PMDHG+
     +  ((2.d0*THETA*PMDH-SCUP)*PMDG-
     -  THETA1*(TMTDG*PMDH+TMTDH*PMDG))/SCDN)/SCDN
*   ----------------------------------------------------------------   *
      SQTH=dsqrt(THETA)
      THETA2=THETA**2
      THETA3=THETA2*THETA
      THETA4=THETA3*THETA
      if (THETA.gt..005) then
         CHT1=dcosh(1.d0/THETA)
         SHT1=dsinh(1.d0/THETA)
         CHT2=dcosh(1.d0/SQTH)
         SHT2=dsinh(1.d0/SQTH)
         T1=SHT1/CHT1 ! dtanh(1./THETA)
         T2=SHT2/CHT2 ! dtanh(1./dsqrt(THETA))
         T1DH=-1.d0/(THETA*CHT1)**2 ! d T1 / d\theta
         T1DHH=2.d0/(THETA*CHT1)**3*(CHT1-SHT1/THETA)
         T2DH=-.5d0*SQTH/(THETA*CHT2)**2
         T2DHH=(.75d0*SQTH*CHT2-.5d0*SHT2)/(THETA*CHT2)**3
      else
         T1=1.d0
         T2=1.d0
         T1DH=0.
         T2DH=0.
         T1DHH=0.
         T2DHH=0.
      endif
      A0=.75d0+3.04363*THETA2-.09227*THETA3+1.7035*THETA4
      A0DH=6.08726*THETA-.27681*THETA2+6.814*THETA3
      A0DHH=6.08726-.55362*THETA+20.442*THETA2
      A1=1.d0+8.31051*THETA2+5.1105*THETA4
      A1DH=16.62102*THETA+20.442*THETA3
      A1DHH=16.62102+61.326*THETA2
      A=.610887*A0/A1*T1 ! HF fit of Perrot and Dharma-wardana
      AH=A0DH/A0-A1DH/A1+T1DH/T1
      ADH=A*AH
      ADHH=ADH*AH+A*(A0DHH/A0-(A0DH/A0)**2-A1DHH/A1+(A1DH/A1)**2+
     +  T1DHH/T1-(T1DH/T1)**2)
      B0=.341308+12.070873d0*THETA2+1.148889d0*THETA4
      B0DH=24.141746d0*THETA+4.595556d0*THETA3
      B0DHH=24.141746d0+13.786668d0*THETA2
      B1=1.d0+10.495346d0*THETA2+1.326623*THETA4
      B1DH=20.990692d0*THETA+5.306492*THETA3
      B1DHH=20.990692d0+15.919476d0*THETA2
      B=SQTH*T2*B0/B1
      BH=.5d0/THETA+T2DH/T2+B0DH/B0-B1DH/B1
      BDH=B*BH
      BDHH=BDH*BH+B*(-.5d0/THETA2+T2DHH/T2-(T2DH/T2)**2+
     +  B0DHH/B0-(B0DH/B0)**2-B1DHH/B1+(B1DH/B1)**2)
      D0=.614925+16.996055d0*THETA2+1.489056*THETA4
      D0DH=33.99211d0*THETA+5.956224d0*THETA3
      D0DHH=33.99211d0+17.868672d0*THETA2
      D1=1.d0+10.10935*THETA2+1.22184*THETA4
      D1DH=20.2187*THETA+4.88736*THETA3
      D1DHH=20.2187+14.66208*THETA2
      D=SQTH*T2*D0/D1
      DH=.5d0/THETA+T2DH/T2+D0DH/D0-D1DH/D1
      DDH=D*DH
      DDHH=DDH*DH+D*(-.5d0/THETA2+T2DHH/T2-(T2DH/T2)**2+
     +  D0DHH/D0-(D0DH/D0)**2-D1DHH/D1+(D1DH/D1)**2)
          E0=.539409+2.522206*THETA2+.178484*THETA4
      E0DH=5.044412*THETA+.713936*THETA3
      E0DHH=5.044412+2.141808*THETA2
      E1=1.d0+2.555501*THETA2+.146319*THETA4
      E1DH=5.111002*THETA+.585276*THETA3
      E1DHH=5.111002+1.755828*THETA2
      E=THETA*T1*E0/E1
      EH=1.d0/THETA+T1DH/T1+E0DH/E0-E1DH/E1
      EDH=E*EH
      EDHH=EDH*EH+E*(T1DHH/T1-(T1DH/T1)**2+E0DHH/E0-(E0DH/E0)**2-
     -  E1DHH/E1+(E1DH/E1)**2-1.d0/THETA2)
      EXP1TH=dexp(-1.d0/THETA)
      C=(.872496+.025248*EXP1TH)*E
      CDH=.025248*EXP1TH/THETA2*E+C*EDH/E
      CDHH=.025248*EXP1TH/THETA2*(EDH+(1.d0-2.d0*THETA)/THETA2*E)+
     +  CDH*EDH/E+C*EDHH/E-C*(EDH/E)**2
      DISCR=dsqrt(4.d0*E-D**2)
      DIDH=.5d0/DISCR*(4.d0*EDH-2.d0*D*DDH)
      DIDHH=(-((2.d0*EDH-D*DDH)/DISCR)**2+2.d0*EDHH-DDH**2-D*DDHH)/DISCR
*   -----------------------------------------------
      S1=-C/E*GAME
      S1H=CDH/C-EDH/E
      S1DH=S1*S1H
      S1DHH=S1DH*S1H+S1*(CDHH/C-(CDH/C)**2-EDHH/E+(EDH/E)**2)
      S1DG=-C/E ! => S1DGG=0
      S1DHG=S1DG*(CDH/C-EDH/E)
      B2=B-C*D/E
      B2DH=BDH-(CDH*D+C*DDH)/E+C*D*EDH/E**2
      B2DHH=BDHH-(CDHH*D+2.d0*CDH*DDH+C*DDHH)/E+
     +  (2.d0*(CDH*D+C*DDH-C*D*EDH/E)*EDH+C*D*EDHH)/E**2
      SQGE=dsqrt(GAME)
      S2=-2.d0/E*B2*SQGE
      S2H=B2DH/B2-EDH/E
      S2DH=S2*S2H
      S2DHH=S2DH*S2H+S2*(B2DHH/B2-(B2DH/B2)**2-EDHH/E+(EDH/E)**2)
      S2DG=.5d0*S2/GAME
      S2DGG=-.5d0*S2DG/GAME
      S2DHG=.5d0*S2DH/GAME
      R3=E*GAME+D*SQGE+1.d0
      R3DH=EDH*GAME+DDH*SQGE
      R3DHH=EDHH*GAME+DDHH*SQGE
      R3DG=E+.5d0*D/SQGE
      R3DGG=-.25d0*D/(GAME*SQGE)
      R3DHG=EDH+.5d0*DDH/SQGE
      B3=A-C/E
      B3DH=ADH-CDH/E+C*EDH/E**2
      B3DHH=ADHH-CDHH/E+(2.d0*CDH*EDH+C*EDHH)/E**2-2.d0*C*EDH**2/E**3
      C3=(D/E*B2-B3)/E ! =D*B2/E**2-B3/E
      C3DH=(DDH*B2+D*B2DH+B3*EDH)/E**2-2.d0*D*B2*EDH/E**3-B3DH/E
      C3DHH=(-B3DHH+
     +  (DDHH*B2+2.d0*DDH*B2DH+D*B2DHH+B3DH*EDH+B3*EDHH+B3DH*EDH)/E-
     -  2.d0*((DDH*B2+D*B2DH+B3*EDH+DDH*B2+D*B2DH)*EDH+D*B2*EDHH)/E**2+
     +  6.d0*D*B2*EDH**2/E**3)/E
      S3=C3*dlog(R3)
      S3DH=S3*C3DH/C3+C3*R3DH/R3
      S3DHH=(S3DH*C3DH+S3*C3DHH)/C3-S3*(C3DH/C3)**2+
     +  (C3DH*R3DH+C3*R3DHH)/R3-C3*(R3DH/R3)**2
      S3DG=C3*R3DG/R3
      S3DGG=C3*(R3DGG/R3-(R3DG/R3)**2)
      S3DHG=(C3DH*R3DG+C3*R3DHG)/R3-C3*R3DG*R3DH/R3**2
      B4=2.d0-D**2/E
      B4DH=EDH*(D/E)**2-2.d0*D*DDH/E
      B4DHH=EDHH*(D/E)**2+2.d0*EDH*(D/E)**2*(DDH/D-EDH/E)-
     -  2.d0*(DDH**2+D*DDHH)/E+2.d0*D*DDH*EDH/E**2
      C4=2.d0*E*SQGE+D
      C4DH=2.d0*EDH*SQGE+DDH
      C4DHH=2.d0*EDHH*SQGE+DDHH
      C4DG=E/SQGE
      C4DGG=-.5d0*E/(GAME*SQGE)
      C4DHG=EDH/SQGE
      S4A=2.d0/E/DISCR
      S4AH=EDH/E+DIDH/DISCR
      S4ADH=-S4A*S4AH
      S4ADHH=-S4ADH*S4AH-
     -  S4A*(EDHH/E-(EDH/E)**2+DIDHH/DISCR-(DIDH/DISCR)**2)
      S4B=D*B3+B4*B2
      S4BDH=DDH*B3+D*B3DH+B4DH*B2+B4*B2DH
      S4BDHH=DDHH*B3+2.d0*DDH*B3DH+D*B3DHH+
     +  B4DHH*B2+2.d0*B4DH*B2DH+B4*B2DHH
      S4C=datan(C4/DISCR)-datan(D/DISCR)
      UP1=C4DH*DISCR-C4*DIDH
      DN1=DISCR**2+C4**2
      UP2=DDH*DISCR-D*DIDH
      DN2=DISCR**2+D**2
      S4CDH=UP1/DN1-UP2/DN2
      S4CDHH=(C4DHH*DISCR-C4*DIDHH)/DN1-
     -  UP1*2.d0*(DISCR*DIDH+C4*C4DH)/DN1**2-
     -  (DDHH*DISCR-D*DIDHH)/DN2+UP2*2.d0*(DISCR*DIDH+D*DDH)/DN2**2
      S4CDG=C4DG*DISCR/DN1
      S4CDGG=C4DGG*DISCR/DN1-2.d0*C4*DISCR*(C4DG/DN1)**2
      S4CDHG=(C4DHG*DISCR+C4DG*DIDH-
     -  C4DG*DISCR/DN1*2.d0*(DISCR*DIDH+C4*C4DH))/DN1
      S4=S4A*S4B*S4C
      S4DH=S4ADH*S4B*S4C+S4A*S4BDH*S4C+S4A*S4B*S4CDH
      S4DHH=S4ADHH*S4B*S4C+S4A*S4BDHH*S4C+S4A*S4B*S4CDHH+
     +  2.d0*(S4ADH*S4BDH*S4C+S4ADH*S4B*S4CDH+S4A*S4BDH*S4CDH)
      S4DG=S4A*S4B*S4CDG
      S4DGG=S4A*S4B*S4CDGG
      S4DHG=S4A*S4B*S4CDHG+S4CDG*(S4ADH*S4B+S4A*S4BDH)
      FXC=S1+S2+S3+S4
      FXCDH1=S1DH+S2DH+S3DH+S4DH
      FXCDG1=S1DG+S2DG+S3DG+S4DG
      FXCDHH1=S1DHH+S2DHH+S3DHH+S4DHH
      FXCDGG1=S2DGG+S3DGG+S4DGG
      FXCDHG1=S1DHG+S2DHG+S3DHG+S4DHG
* -------- Magnetic correction of derivatives
       FXCDH=FXCDH1*THDH
       FXCDG=FXCDG1+FXCDH1*THDG
       FXCDHH=FXCDHH1*THDH**2+FXCDH1*THDHH
       FXCDGG=FXCDGG1+FXCDH1*THDGG+(2.d0*FXCDHG1+FXCDHH1*THDG)*THDG
       FXCDHG=FXCDHG1*THDH+FXCDH1*THDHG+FXCDHH1*THDH*THDG
* --------
      PXC=(GAME*FXCDG-2.d0*THETA1*FXCDH)/3.d0
      UXC=GAME*FXCDG-THETA1*FXCDH
      SXC=(GAME*S2DG-S2+GAME*S3DG-S3+S4A*S4B*(GAME*S4CDG-S4C))-
     -  THETA1*FXCDH
      if (dabs(SXC).lt.1.d-9*dabs(THETA1*FXCDH)) SXC=0. ! accuracy loss
      CVXC=2.d0*THETA1*(GAME*FXCDHG-FXCDH)-
     -  THETA1**2*FXCDHH-GAME**2*FXCDGG
      if (dabs(CVXC).lt.1.d-9*dabs(GAME**2*FXCDGG)) CVXC=0. ! accuracy
      PDLH=THETA1*(GAME*FXCDHG-2.d0*FXCDH-2.d0*THETA1*FXCDHH)/3.d0
      PDLG=GAME*(FXCDG+GAME*FXCDGG-2.d0*THETA1*FXCDHG)/3.d0
      PDRXC=PXC+(PDLG-2.d0*PDLH)/3.d0
      PDTXC=GAME*(THETA1*FXCDHG-GAME*FXCDGG/3.d0)-
     -  THETA1*(FXCDH/.75d0+THETA1*FXCDHH/1.5d0)
      return
      end

      function darcth(X)
*                                                       Version 26.03.98
      implicit double precision (A-H), double precision (O-Z)
      if (dabs(X).ge.1.d0) stop'darcth: X out of range.'
      if (dabs(X).lt.1.d-9) then
         D=X
      else
         D=.5d0*dlog((1.d0+X)/(1.d0-X))
      endif
      darcth=D
      return
      end

* ==============   SUBROUTINES FOR THE SOLID STATE   ================= *
      subroutine HLMAG(TPT,BP,UMth,CMth,SMth,PMth,PDTth,PDRth,
     *  U1m,dLNu1,d2LNu1,DU1,UM,SM,CMU,CMS)
*                                                       Version 09.07.12
* Input: TPT=T_p/T
*        BP=\omega_c/T_p
* Output: UMth - thermal energy / NkT
*         CMth - specific heat/Nk
*         SMth - entropy/Nk
* Auxiliary output:
*         U1m=u1=<\omega_{phonon}/\omega_p> - normalized 1st moment
*             NB: zero-point cyclotron oscillations are excluded,
*             i.e., \omega_{phonon}=\omega_{phonon,total}-\omega_c/3.
*         dLNu1=d\ln(u1)/d\ln(BP)
*         d2LNu1=d^2\ln(u1)/[d\ln(BP)]^2
*         DU1 - variation of u1 with changing direction of the field
* Additional provisional output:
*         UM - a simpler fit to thermal energy / NkT
*         SM - a simpler fit to entropy / Nk
*         CMU - specific heat derived from UM
*         CMS - specific heat derived from SM, = d(SM)/dt,
*              where t=ln(TTP)=-ln(TPT)
* Useful internal variables:
*         dUMdt=d(UM)/dt, dUMdb=d(UM)/db, where b=ln(BP),
*         dUMdtt=d^2(UM)/dt^2, dUMdbb=d^2(UM)/db^2, dUMdtb=d^2(UM)/dtdb,
*         dSMdb=d(SM)/db, dCMSd=d(CMS)/dt=d^2(SM)/dt^2,
*         dCMSdt=d^2(SM)/dt^2, dSMdbb=d^2(SM)/db^2, dSMdtb=d^2(SM)/dtdb
      implicit double precision (A-H), double precision (O-Z)
      parameter(TINY=1.d-19)
      parameter(Zb=12.5d0,A1=24.d0,A3=15.d0,A4=119.d0,A6=0.4d0)
      parameter(A4S=9.8d0,A3S=7.9d0,A1S=42.d0)
      parameter(U1inf=0.1801,CU1=1.27)
      if (TPT.lt.0.) stop'HLMAG: TPT < 0'
      TTP=1.d0/(TPT+TINY)
      call HLfit12(TPT,F,U,CVth,Sth,U1,CWth,1)
      if (BP.lt.0.) stop'HLMAG: BP < 0'
      if (BP.eq.0.) then
         UMth=U
         CMth=CVth
         SMth=Sth
         UM=U
         CMU=CVth
         CMS=CVth
         SM=Sth
         PMth=U/2.d0
         PDTth=CVth/2.d0
         PDRth=.75d0*U-.25d0*CVth
         U1m=U1
         dLNu1=0.
         d2LNu1=0.
         DU1=0.
      else
         TA=TTP/BP
         TB=TTP*BP
         SqrtTB=sqrt(TB)
* U:
         TA3=3.d0*TA
         TQ3=(TA3)**2*dsqrt(TA3)
         Q=.5d0/(1.d0+TQ3)**A6
         UMH=U/(1.d0+Q)
         TV=(Zb*SqrtTB+A4*TB)*TB
         DV=1.d0+TV+A3*TB
         TW=A1*TTP**2
         DW=1.d0+TW
         UML=TV/DV/dsqrt(DW)
         UM=UML+UMH
* S:
         TVS=Zb*(SqrtTB/.6d0+A4S*TB)*TB
         DVS=1.d0+A3S*TB
         TWS=A1S*TTP**2
         DWS=1.d0+TWS
         DQS=(1.+TA)**4
         TZ=TVS/DVS/DWS/DQS
         SML=dlog(1.+TZ)
         SM=Sth+SML
* Derivatives for UM:
** over t=ln(TTP) -
         dTQ3dt=TQ3*2.5d0
         dQdt=-Q*dTQ3dt*A6/(1.d0+TQ3)
         dUdt=CVth-U
         dUMHdt=dUdt/(1.d0+Q)-U/(1.d0+Q)**2*dQdt
         dTVdt=(Zb*1.5d0*SqrtTB+2.d0*A4*TB)*TB
         dDVdt=dTVdt+A3*TB
         dDWdt=2.d0*TW
         DML=dTVdt/TV-dDVdt/DV-.5d0*dDWdt/DW
         dUMLdt=UML*DML
         dUMdt=dUMLdt+dUMHdt
         CMU=dUMdt+UM
** over b=ln(BP) -
         dQdb=-dQdt ! because Q=Q(TTP/BP):=Q(t-b)
         dUMHdb=-U/(1.d0+Q)**2*dQdb
         dTVdb=dTVdt
         dDVdb=dDVdt ! because DV=DV(TTP*BP):=DV(t+b)
         DMLb=dTVdb/TV-dDVdb/DV ! because DW=DW(TTP) => dDW/db=0
         dUMLdb=UML*DMLb
         dUMdb=dUMLdb+dUMHdb
* Derivatives for SM:
** over t=ln(TTP) -
         dTVSdt=Zb*(SqrtTB*2.5+A4S*TB*2.d0)*TB
         dDWSdt=2.d0*TWS
         DlnDQSdt=TA*4.d0/(1.+TA)
         DTZ=dTVSdt/TVS-A3S*TB/DVS-dDWSdt/DWS-DlnDQSdt
         dTZdt=TZ*DTZ
         CVL=dTZdt/(1.d0+TZ)
         CMS=CVL+CVth
** over b=ln(BP) -
         dTVSdb=dTVSdt
         DlnDQSdb=-DlnDQSdt
         DTZb=dTVSdb/TVS-A3S*TB/DVS-DlnDQSdb ! because dDWS/db=0
         dTZdb=TZ*DTZb
         dSMdb=dTZdb/(1.d0+TZ) ! because dSth/db=0
* Second derivatives:
** for U:
*** over tt
         dQdtt=-dQdt*dTQ3dt*A6/(1.d0+TQ3)+2.5d0*dQdt+
     +     Q*A6*(dTQ3dt/(1.d0+TQ3))**2
         dTVdtt=(Zb*2.25d0*SqrtTB+4.d0*A4*TB)*TB
         dDVdtt=dTVdtt+A3*TB
         dDWdtt=2.d0*dDWdt
         dDMLdt=dTVdtt/TV-dDVdtt/DV-.5d0*dDWdtt/DW-
     -     (dTVdt/TV)**2+(dDVdt/DV)**2+.5d0*(dDWdt/DW)**2
         dUMLdtt=dUMLdt*DML+dDMLdt*UML
         dUdtt=CWth-dUdt
         dUMHdtt=dUdtt/(1.d0+Q)-2.d0*dUdt/(1.d0+Q)**2*dQdt-
     -     U/(1.d0+Q)**2*dQdtt+2.d0*U/(1.d0+Q)**3*dQdt**2
         dUMdtt=dUMLdtt+dUMHdtt
         dCMUdt=dUMdtt+dUMdt
*** over bb
         dQdbb=dQdtt ! Q(t-b)
         dTVdbb=dTVdtt
         dDVdbb=dDVdtt ! DV(t+b)
         dDMLbdb=dTVdbb/TV-dDVdbb/DV-(dTVdb/TV)**2+(dDVdb/DV)**2
         dUMLdbb=dUMLdb*DMLb+dDMLbdb*UML
         dUMHdbb=-U/(1.d0+Q)**2*dQdbb+2.d0*U/(1.d0+Q)**3*dQdb**2
         dUMdbb=dUMLdbb+dUMHdbb
*** over tb
         dQdtb=-dQdtt ! Q(t-b)
         dTVdtb=dTVdtt
         dDVdtb=dDVdtt ! DV(t+b)
         dDMLbdt=dTVdtb/TV-dDVdtb/DV-dTVdt*dTVdb/TV**2+dDVdt*dDVdb/DV**2
         dUMLdtb=dUMLdt*DMLb+dDMLbdt*UML
         dUMHdtb=-(dUdt*dQdb+U*dQdtb)/(1.d0+Q)**2+
     +     2.d0*U/(1.d0+Q)**3*dQdb*dQdt
         dUMdtb=dUMLdtb+dUMHdtb
** for S:
*** over tt
         dTVSdtt=Zb*(SqrtTB*3.75d0+A4S*TB*4.d0)*TB
         dDWSdtt=4.d0*TWS
         dDTZdt=dTVSdtt/TVS-A3S*TB/DVS-dDWSdtt/DWS-DlnDQSdt-
     -  (dTVSdt/TVS)**2+(A3S*TB/DVS)**2+(dDWSdt/DWS)**2+DlnDQSdt**2/4.d0
         dTZdtt=dTZdt*DTZ+TZ*dDTZdt
         dCVLdt=dTZdtt/(1.d0+TZ)-(dTZdt/(1.d0+TZ))**2
         dCMSdt=dCVLdt+CWth
*** over bb
         dTVSdbb=dTVSdtt
         dDTZbdb=dTVSdbb/TVS-A3S*TB/DVS+DlnDQSdb-
     -     (dTVSdb/TVS)**2+(A3S*TB/DVS)**2+DlnDQSdb**2/4.d0
         dTZdbb=dTZdb*DTZb+TZ*dDTZbdb
         dSMdbb=dTZdbb/(1.d0+TZ)-(dTZdb/(1.d0+TZ))**2
*** over tb
         dTVSdtb=dTVSdtt
         dDTZbdt=dTVSdtb/TVS-A3S*TB/DVS-DlnDQSdb-
     -     dTVSdb*dTVSdt/TVS**2+(A3S*TB/DVS)**2-DlnDQSdb**2/4.d0
         dTZdtb=dTZdt*DTZb+TZ*dDTZbdt
         dSMdtb=dTZdtb/(1.d0+TZ)-dTZdt*dTZdb/(1.d0+TZ)**2
* Corrected values:
         SMth=SM+CMS-CMU
         UMth=UM+CMS-CMU
         CMth=CMS+dCMSdt-dCMUdt
         PMth=.5d0*(CMS+dSMdb-dUMdt-dUMdb)
         PDTth=PMth+.5d0*(dCMSdt+dSMdtb-dUMdtt-dUMdtb)
         PDRth=PMth+.25d0*(dUMdtt+2.d0*dUMdtb+dUMdbb-
     -     dCMSdt-2.d0*dSMdtb-dSMdbb)
* First moment of phonon frequencies and its variation:
         BC=CU1*BP*dsqrt(dsqrt(dsqrt(BP)))
         U1up=U1+BC*U1inf
         U1dn=1.d0+BC
         U1m=U1up/U1dn
         dLNu1=BC*1.125d0*(U1inf/U1up-1.d0/U1dn)
         d2LNu1=BC*CU1*1.265625d0*(U1*U1inf/U1up**2-1.d0/U1dn**2)
         B32=BP*dsqrt(BP)
         DU1=.0064*B32**2*dsqrt(B32)/((1.d0+B32)**2*dsqrt(1.d0+B32))
      endif
      return
      end

      subroutine FHARMAG(GAMI,TPT,BP,
     *     FMHL,UMHL,PMHL,CVMHL,SMHL,PDTMHL,PDRMHL,EVAR)
* Thermodynamic functions of a harmonic crystal in magnetic field
* 
*                                                       Version 01.06.12
* Input: GAMI - ionic Gamma, TPT=T_{p,i}/T
* Output: FMHL=F/(N_i T), UMHL=U/(N_i T), PMHL=P/(n_i T),
* CVMHL=C_V/N_i, SMHL=S/N_i
* PDTMHL = PMHL + d PMHL / d ln T, PDRMHL = PMHL + d PMHL/d ln\rho
* EVAR - energy variation with changing field direction
*  NB: zero-point cyclotron oscillations are excluded by U1m definition,
*  i.e., E0=E0(total)-\hbar\omega_c/(2kT)
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter(CM=.895929256d0) 
      call HLMAG(TPT,BP,UMth,CVMHL,SMHL,PMth,PDTMHL,PDRth,
     *  U1m,dLNu1,d2LNu1,DU1,UM,SM,CMU,CMS)
      U0=-CM*GAMI ! perfect lattice
      E0=1.5d0*U1m*TPT ! zero-point energy
      P0=.5d0*E0*(1.d0-dLNu1) ! zero-point pressure
      PDR0=E0*(.75d0-dLNu1+.25d0*dLNu1**2+.25d0*d2LNu1)
      EVAR=1.5d0*DU1*TPT
      UMph=UMth+E0 ! phonons: thermal + zero-point
      FMph=UMph-SMHL
      PMph=PMth+P0
      PDRph=PDRth+PDR0
      UMHL=U0+UMph
      FMHL=U0+FMph
      PMHL=U0/3.d0+PMph
      PDRMHL=U0/2.25d0+PDRph
      return
      end

      subroutine SPINION(GYRO,ZETI,MULTI,Fspin,Uspin,CVspin)
*                                                       Version 28.04.12
* Spin-degeneracy part of the entropy.
* Input: GYRO - gyrofactor
*        ZETI=\hbar\omega_c/kT
*        MULTI - spin multiplicity:
*                if MULTI=2, then assume fermions, otherwise bosons.
* Output: Fspin - spin free energy contribution / NkT
*         Uspin - spin internal energy contribution / NkT
*         CVspin - spin specific heat contribution / Nk
      implicit double precision (A-H), double precision (O-Z)
      parameter(TINY=1.d-7)
      CZ4=GYRO*ZETI*.25d0
      if (MULTI.gt.1.and.dabs(GYRO).gt.TINY) then ! add spin terms
        if (MULTI.eq.2) then
           CH=dcosh(GZ4)
           Fspin=-dlog(2.d0*CH)
           Uspin=-GZ4*dtanh(GZ4)
           CVspin=(GZ4/CH)**2
        else
           GZM4=GZ4*MULTI
           Fspin=-dlog(dsinh(GZM4)/dsinh(GZ4))
           Uspin=GZ4/dtanh(GZ4)-GZM4/dtanh(GZM4)
           CVspin=(GZ4/dsinh(GZ4))**2-(GZM4/dsinh(GZM4))**2
        endif
      else
         Fspin=0.
         Uspin=0.
         CVspin=0.
      endif
      return
      end

      subroutine QDMAG(TPT,ZETI,FWK,UWK,PWK,CVWK,PDTWK,PDRWK)
*                                                       Version 24.05.12
* Wigner-Kirkwood (quantum diffraction) corrections in a magnetic field
* Ref.: Alastuey A. & Jancovici B. (1980) Physica A 102, 327.
* Input: TPT - quantum plasma parameter (T_p/T)
*        ZETI=\hbar\omega_c/kT - cyclotron-to-thermal energy ratio
* Output: FWK, UWK, PWK - Wigner-Kirkwood contributions to the free and
*           internal energies and pressure, divuded by (n k T)
*         CVWK - contribution to the specific heat divided by (N k)
*         PDTWK=(1/nkT)*(d PWK/d ln T)_V
*         PDRWK=(1/nkT)*(d PWK/d ln\rho)_T
      implicit double precision (A-H), double precision (O-Z)
      if (ZETI.lt.0.) stop'WK-mag: ZETI < 0'
      if (ZETI.lt.1.d-2) then
         F1=-ZETI**2/11.25d0 ! u*f'(u)
         F0=1.d0+F1/2.d0 ! f(u)
         F2=F1 ! u^2*f''(u)
      elseif (ZETI.lt.20.) then
         CH=dcosh(ZETI)
         SH=dsinh(ZETI)
         F0=2.d0/ZETI*(CH/SH-1.d0/ZETI)+1.d0/3.d0 ! f(u)
         F1=4.d0/ZETI**2-2.d0/SH*(1.d0/SH+CH/ZETI) ! u*f'(u)
         F2=4.d0/SH*(CH/ZETI+ZETI*CH/SH**2+1.d0/SH)-12.d0/ZETI**2
      else ! u > 20, CH=SH\to\infty, CH/SH=1.
         F0=1.d0/3.d0+2.d0/ZETI-2.d0/ZETI**2 ! f(u)
         F1=4.d0/ZETI**2-2.d0/ZETI ! u*f'(u)
         F2=4.d0/ZETI-12.d0/ZETI**2 ! u^2*f''(u)
      endif
      WK=TPT**2/24.d0
      FWK=WK*F0
      UWK=WK*(2.d0*F0+F1)
      CVWK=-WK*(2.d0*F0+4.*F1+F2)
      PWK=FWK
      PDRWK=2.d0*PWK
      PDTWK=-WK*(F0+F1)
      return
      end
