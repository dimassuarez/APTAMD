      program Gcav
      implicit double precision (a-h,o-z)
C
C     We use three methods: 
C     Scaled Particle Theory formulae fed
C     with effectiv Rs and Rv radii and the 
C
C     The revised Pierotti formula of HOFINGER (fed with Rv of V-SES)
C     (Eq 2 in Chem. Soc. Rev., 2005, 34, 1012-1020) 
C
C     The empirical method of Grigoriev et al (JPC-B 2007 111  13748)
C     in which TI data are fitted against V ( Nota that V-SAS is required here)
C
      parameter(maxatm=100000)
      parameter (pi=3.14159265358979323846)
      parameter (ev=0.0361,lambda=22.3)
      parameter (rk1=0.427,rk2=-1.594,rk3=1.183)
      parameter (rd1=0.012,rd2=0.044,rd3=0.046)
      dimension surf(maxatm),ratom(maxatm)
C
      READ(5,*) nat
      READ (5,*) vses,sesnum
      READ(5,*)  sesana
      Sum=0.0
      DO I=1,NAT
         READ(5,*) surf(i)
         Sum=Sum+Surf(i)
      ENDDO
      DO I=1,NAT
         READ(5,*) ratom(i)
      ENDDO
C
      RS=SQRT(sesana/(4.0*pi))
      RV=((3.0*vses)/(4.0*pi))**(1.0/3.0)
C
C     Grigoriev (An approximated VSAS is used) 
C
      CHI=1.0D0 - (RV/RS) 
      GGRI= (EV*(VSES + SESANA*1.4) + ( lambda * CHI) )
C
C     Hofinger (Parameters at 300K  from Phys. Chem. Chem. Phys., 2006, 8, 5515)
C
      GHOF=RK1 + RK2*RV + RK3*RV*RV
C
C     SPT 
C
      GS=PIERCLAV(RS)
      GV=PIERCLAV(RV)
      GC=0.0
      DO I=1,NAT
        Gat=PIERCLAV(Ratom(I))
        GC=GC+(Surf(I)/(4.0*pi*(Ratom(I)**2)))*Gat
      ENDDO
C
      WRITE(6,'(2X,''GCAV(kcal/mol)  SPT-S= '',F10.3,
     &  '' STP-V= '',F10.3,
     &  '' C-SPT= '',F10.3,
     &  '' G_GRI= '',F10.3,
     &  '' G_HOF= '',F10.3,
     &  ''    RV= '',F10.3,
     &  ''    RS= '',F10.3,
     &  ''   SES= '',F10.3,
     &  ''  VSES= '',F10.3)') GS,GV,GC,GGRI,GHOF,RV,RS,sesana,vses
      stop
      end
C
      FUNCTION PIERCLAV(R)
      implicit double precision (a-h,o-z)
C
      Parameter(RWAT=1.40)
      Parameter(chi=0.383) 
      Parameter (pi=3.14159265358979323846)
      Parameter (RT=0.596164)
C
C    Adapted from Tomasi  Chem. Rev. 2005, 105, 2999-3093 Eq 82
C
C    Parameters taken from J. Phys. Chem. B 2006, 110, 11421-11426
C
      G=0.0
      G=-LOG(1.0-chi)
      G=G+( (3*chi)/(1.0-chi) )*(R/RWAT)
      G=G+( (3*chi)/(1.0-chi) + (9.0/2.0)*
     &      ((chi/(1.0-chi))**2) ) * ((R/RWAT)**2)    
      G=RT*G
C
      PIERCLAV=G
C
      RETURN
      END
       
  


