      PROGRAM fixq_pdbqt
      implicit real*8 (a-h,o-z)
C
C     D. Suarez 
C
      parameter(maxatm=500000,maxres=50000)
C
C     ADCK Gasteiger charges for terminal DNA residues
C     are not properly handled by the AUTODOCKTOOLS in MGLTOOLS 
C     This program fixes those charges
C
C     It also takes care of Na+ and Water molecules in the RDPDB
C     subroutine 
C
C     Fixes also Residue  Names. prepare_receptor4/prepare_ligand4
C     can modify residue names (e.g., CYX-->CYS). We revert
C     this by reding also the PDB file that generated 
C     the pdbqt and using its resnames  
C
C     PDBQT 1
C
      dimension x1(maxatm),y1(maxatm),z1(maxatm)
      dimension q1(maxatm),rad1(maxatm),qdock1(maxatm),qw(maxatm)
      character*3  atype1(maxatm) 
      dimension nfa1(maxres),nla1(maxres)
      dimension Idatm1(maxres),Idres1(maxres)
      character resnam1(maxres)*3
      character ta1(maxatm)*4
C     PDB 2
      dimension x2(maxatm),y2(maxatm),z2(maxatm)
      character*3  atype2(maxatm) 
      dimension nfa2(maxres),nla2(maxres)
      dimension Idatm2(maxres),Idres2(maxres)
      character resnam2(maxres)*3
      character ta2(maxatm)*4
C
      character line*120
      character finp1*256,fout*256,fref*256
      logical eof
C
C    Beginning, Okay!  
C
      finp1=' '
      fout=' '
      fref=' '
      CALL GETARG(1,finp1)
      CALL GETARG(2,fout)
      CALL GETARG(3,fref)
      IF ((finp1 .eq.' ').or.(fout.eq.' ').or.(fref.eq.' '))THEN 
        STOP 'Usage: fixq_pdbqt pdbqtinp  pdbqtout pdbref'
      ENDIF
C
C     Reading File 1
C
      OPEN (10,file=finp1,status='unknown')
      CALL RDPDBQT (10,NAT1,NRES1,NFA1,NLA1,TA1,RESNAM1,
     &     IDATM1,IDRES1,X1,Y1,Z1,Q1,RAD1,QDOCK1,ATYPE1)
      qsum=0.0d0
      DO i=1,nat1
        qsum=qsum+qdock1(i)
      ENDDO
      WRITE(6,'(''Sum of charges in PDBQT input file ='',F8.3)') qsum
C
      NUMION=0
      NUMWAT=0
      DO IRES=1,NRES1
        IF ( INDEX(RESNAM1(IRES),'WAT').ne.0 ) NUMWAT=NUMWAT+1 
        IF ( INDEX(RESNAM1(IRES),'Na+').ne.0 ) NUMION=NUMION+1 
      ENDDO
C
      IF (NUMION.gt.0) 
     &WRITE(6,'(''PDBQT data for '',I2,'' Na+ ions fixed'')') NUMION
      IF (NUMWAT.gt.0) 
     &WRITE(6,'(''PDBQT data for '',I3,'' waters fixed'')') NUMWAT
C
C
C     Fixing charges (DG Guanine residues 
C
      IDG=0
      DO IRES=1,NRES1
         qg0=0.0d0
         qg=0.0d0
         IF ( INDEX(RESNAM1(IRES),'DG').ne.0 ) THEN
            IDG=IDG+1
            DO IAT=NFA1(IRES),NLA1(IRES)
               qg0=qg0+qdock1(iat)
            ENDDO
            IF ( INDEX(RESNAM1(IRES),'DG5') .ne. 0 ) THEN
                JAT=NFA1(IRES)+6
            ELSE
                JAT=NFA1(IRES)+8
            ENDIF
            KAT=JAT+13
            DO IAT=JAT,KAT
               qa0=QDOCK1(IAT)
               if (INDEX(TA1(IAT),'N9').ne.0) QDOCK1(IAT)=-0.29377d0
               if (INDEX(TA1(IAT),'C8').ne.0) QDOCK1(IAT)= 0.22928d0
               if (INDEX(TA1(IAT),'N7').ne.0) QDOCK1(IAT)=-0.27585d0
               if (INDEX(TA1(IAT),'C5').ne.0) QDOCK1(IAT)= 0.20062d0
               if (INDEX(TA1(IAT),'C6').ne.0) QDOCK1(IAT)= 0.32362d0
               if (INDEX(TA1(IAT),'O6').ne.0) QDOCK1(IAT)=-0.32123d0
               if (INDEX(TA1(IAT),'N1').ne.0) QDOCK1(IAT)=-0.30093d0
               if (INDEX(TA1(IAT),'H1').ne.0) QDOCK1(IAT)= 0.18868d0
               if (INDEX(TA1(IAT),'C2').ne.0) QDOCK1(IAT)= 0.21973d0
               if (INDEX(TA1(IAT),'N2').ne.0) QDOCK1(IAT)=-0.39408d0
               if (INDEX(TA1(IAT),'H2').ne.0) QDOCK1(IAT)= 0.17316d0
               if (INDEX(TA1(IAT),'N3').ne.0) QDOCK1(IAT)=-0.23883d0
               if (INDEX(TA1(IAT),'C4').ne.0) QDOCK1(IAT)= 0.19346d0
c              WRITE(6,'(A3,'' q0='',F8.3,'' q='',F8.3)')
c    &         TA1(IAT),qa0,QDOCK1(IAT)
            ENDDO
            DO IAT=NFA1(IRES),NLA1(IRES)
               qg=qg+qdock1(iat)
            ENDDO
            WRITE(6,'( A5,1X,I3,'', q0='',F8.3,'' q_fixed='',F8.3)')
     &      RESNAM1(IRES),IRES,qg0,qg
         ENDIF
      ENDDO
C
      IF  ( IDG .eq. 0 ) THEN
         WRITE(6,'(''No guanine residues found -> No fix applied'')')
      ELSE
         WRITE(6,'(I3,'' guanine residues found. Fix applied'')') IDG
         qsum=0.0d0
         DO i=1,nat1
           qsum=qsum+qdock1(i)
         ENDDO
         WRITE(6,'(''Sum of fixed charges ='',F8.3)') qsum
      ENDIF
C
C     Identifying Z ADCK TYPE ATOMS: 
C
      NDZ=0
      qZ=0.0d0
      DO IRES=1,NRES1
         DO IAT=NFA1(IRES),NLA1(IRES)
            IF ( INDEX(ATYPE1(IAT),' Z ') .ne. 0 ) THEN
              NDZ=NDZ+1
              WRITE(6,'( A5,1X,I3,1X,A4 ,
     &        '' Z covalent atom -> charged nullified '')') 
     &        RESNAM1(IRES),IRES,TA1(IAT)
              qz=qz+QDOCK1(IAT)
              QDOCK1(IAT)=0.0d0
           ENDIF
         ENDDO
      ENDDO
C
C     Original CHARGE of Z ADCK atoms is redistributed 
C     to preserve charge integrity 
C
      IF ( NDZ .gt. 0 )  THEN
         qZ=qZ/DFLOAT(NAT1-NDZ) 
         DO IAT=1,NAT1
            IF ( INDEX(ATYPE1(IAT),' Z ') .eq. 0 ) 
     &      QDOCK1(IAT)=QDOCK1(IAT)+qZ 
         ENDDO
      ENDIF
C
C     Rounding off : an elaborate protocol is needed to
C     minimize truncation errors due to F10.3 format of charges
C
      DO IAT=1,NAT1
         QDOCK1(IAT)=dfloat(nint(QDOCK1(IAT)*1000d0))/1000.d0
      ENDDO
      qsum=0.0d0
      DO i=1,nat1
           qsum=qsum+qdock1(i)
      ENDDO
      WRITE(6,'(''Sum of charges with F8.3 format='',F8.3)') qsum
      qmax=maxval(abs(QDOCK1(1:NAT1)))
      qmean=sum(abs(QDOCK1(1:NAT1)))/dfloat(NAT1)
      q2mean=sum(abs(QDOCK1(1:NAT1)**2))/dfloat(NAT1)
      qsig2=q2mean - qmean*qmean 
      qsig2=qsig2*0.001d0
C  qdiff correction is weighed by a very narrow normal distribution of
C  of absolute value charges centered at the maximum absolute value
      qw=0.0d0
      DO IAT=1,NAT1
         x = - ( ( abs(QDOCK1(IAT)) - qmax )**2)/(2.0d0*qsig2)
         if ( x .lt. -50.d0 ) then
            qw(IAT)=0.0d0
         else
         qw(IAT)= (1.0d0/sqrt( 2.0d0 * 3.14159265359d0 * qsig2) )  * 
     &   exp ( x )
         endif
      ENDDO
      qw=qw/sum(qw) 
      qfsum= sum( dfloat ( nint ( QDOCK1(1:NAT1)*1000d0))  /1000.d0 )
      iq=NINT(qfsum)
      qdiff= qfsum - dfloat(iq)
      DO IAT=1,NAT1
         QDOCK1(IAT)=dfloat(nint(QDOCK1(IAT)*1000d0))/1000.d0
     &                    - qw(IAT)*qdiff
      ENDDO
      qfsum= sum( dfloat(nint(QDOCK1(1:NAT1)*1000d0))/1000.d0 )
      WRITE(6,
     &'(''Sum of charges in output PDBQT file='',F8.3)')
     &          qfsum

C
C    Reading PDB reference for checking residue names
C
      OPEN (13,file=fref,status='unknown')
      CALL RDPDBXYZ (13,NAT2,NRES2,NFA2,NLA2,TA2,RESNAM2,
     &     IDATM2,IDRES2,X2,Y2,Z2)
      CLOSE(13)
      IF ( NRES2 .ne. NRES1 ) THEN
         WRITE(6,'(''Reference PDB='',A,'' NRES='',I6,'' not equal'',
     &'' to NRES='',I6,'' in input PDBQT='',A)') fref,NRES2,NRES1,finp
         WRITE(6,'(''RESIDUE NAME CHECKING DEACTIVATED'')')
      ELSE
         WRITE(6,'(''Taking residue names from '',A)') fref
         DO IRES=1,NRES1
         resnam1(ires)=resnam2(ires)
         ENDDO
      ENDIF
C
C     Printing out 
C
      WRITE(6,'(''Printing output PDBQT file '')')
      OPEN (11,file=fout,status='unknown')
C
      REWIND(10) 
      IAT=0
      DO WHILE ( .NOT. EOF(line,10)  )
        IF ((LINE(1:4).eq.'ATOM' ).or.(LINE(1:6).eq.'HETATM')) THEN
          IAT=IAT+1
          IRES=IWHICHRES(IAT,NRES1,NFA1,NLA1)
             write(11,'(''ATOM'',
     +       2X,I5,1X,A4,1X,A3,1X,I5,4X,3F8.3,2F6.2,F10.3,A3,''   '')')
     +       Idatm1(iat),ta1(iat),resnam1(ires),Idres1(ires),
     +       x1(iat),y1(iat),z1(iat),q1(iat),rad1(iat),qdock1(iat),
     +       atype1(iat)
          IATOM=1
        ELSE
          L=0
          inner: DO J=120,1,-1
            IF (LINE(J:J) .ne. '') THEN
                L=J
                EXIT  inner
            ENDIF
          END DO  inner
          WRITE(11,'(A)') LINE(1:L)
        ENDIF
      ENDDO
C
      CLOSE(10)
      CLOSE(11)
C
      STOP
C
      END
C
C
C"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
C
      SUBROUTINE RDPDBQT (IUNIT,NUMATM,NUMRES,NFA,NLA,TA,RESNAM,
     &                  IDATM,IDRES,XR,YR,ZR,Q,RAD,QDOCK,ATYPE)
C     =============================================================
C
C"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
      implicit real*8  (a-h,o-z)
      parameter(maxatm=500000,maxres=50000)
C
C     All cards  ATOM or HETATM in IUNIT are
C     read to fill the coordinates and data for the protein
C
C     NOTE THAT THIS ROUTINE ASSUMES THAT BIG PDB FILES
C     ARE WELL WRITTEN, THAT IS, THAT RESIDUE NUMBER
C     IS WRITTEN WITH FORMAT I5
C
      integer iunit, numatm, numres
      dimension xr(*),yr(*),zr(*)
      dimension q(*),rad(*),qdock(*)
      character*3 atype(*)
      dimension nfa(*),nla(*),idres(*),idatm(*)
      character*4 ta(*)
      character*3 resnam(*) 
      character*80 line
C
      character head*4,tdum*3,tadum*4
C
      numatm=0
      numres=0
      iat=0
      ires=0
      idum0=0
      jres=0
C
      rewind(iunit)
 10   read(iunit,'(A4)',end=20,err=666) head
C
      IF (( head .eq. 'ATOM') .OR. 
     +         ( head .eq. 'HETA'))  THEN 

          backspace(iunit)
          iat=iat+1
          IF ( iat .gt. maxatm ) THEN 
               WRITE(6,'(''Reading Unit'',I3)') IUNIT  
               READ (IUNIT,'(A)') line
               WRITE(6,'(''Last line:'',/,A)') line   
               WRITE(6,'(''MAXATM='',I6)') maxatm 
               STOP '*** Too many atoms in PDB'
          ENDIF 
          tadum=''
          tdum=''
          read(iunit,'(12X,A4,1X,A3)',err=666,end=666) tadum,tdum
          backspace(iunit)
          IF ((index(tdum,'Na+').ne.0).or.(index(tdum,'WAT').ne.0)) THEN

              read(iunit,'(6X,I5,1X,A4,1X,A3,1X,I5,4X,3F8.3)',
     &        err=666,end=666)
     &        jat,tadum,tdum,jres,xr(iat),yr(iat),zr(iat)
              q(iat)=1.000d0
              rad(iat)=0.000d0
              IF (index(tdum,'Na+').ne.0) THEN
                  qdock(iat)=1.000d0
                  atype(iat)=' SD'
              ELSEIF (index(tadum,'O').ne.0) THEN 
                  qdock(iat)=-0.834d0
                  atype(iat)=' OA'
              ELSE
                  qdock(iat)=0.417d0
                  atype(iat)=' HD'
              ENDIF

          ELSE

          read(iunit,'(6X,I5,1X,A4,1X,A3,1X,I5,4X,3F8.3,2F6.2,F10.3,A3)'
     &    ,err=666,end=666)
     &     jat,tadum,tdum,jres,xr(iat),yr(iat),zr(iat),
     &     q(iat),rad(iat),qdock(iat),atype(iat)

          ENDIF
          Idatm(iat)=jat
C Checking prime issue
          l=0
          do i=1,4
            if ( tadum(i:i) .ne. '') l=l+1
          enddo
          if (tadum(1:1) .eq. '''' )  then
c             print*,'Problem tadum=',tadum 
              do i=1,l-1
                tadum(i:i)=tadum(i+1:i+1)
              enddo
              tadum(l:l)=''''
c             print*,'Fixed tadum=',tadum 
          endif
          ta(iat)=tadum
          IF ( idum0 .ne. jres) THEN
            ires=ires+1
            nfa(ires)=iat
            resnam(ires)=tdum
            idres(ires)=jres
            idum0=jres
          ENDIF
C
      ENDIF
C
      GOTO 10
C
 20   CONTINUE
C
      numatm=iat
      numres=ires
      DO ires=1,numres-1
        nla(ires)=nfa(ires+1)-1
      ENDDO
      nla(numres)=numatm
C
      IF (NUMATM .eq. 0) THEN 
        WRITE(6,*) 'Zero Atoms in PDB file UNIT= ',Iunit
        STOP
      ENDIF
C
      RETURN
C
 666  WRITE(6,*) 'PROBLEMS READING PDB FILE IN UNIT=', Iunit
      STOP
      END
C"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
C
      SUBROUTINE RDPDBXYZ (IUNIT,NUMATM,NUMRES,NFA,NLA,TA,RESNAM,
     &                   IDATM,IDRES,XR,YR,ZR)
C     =============================================================
C
C"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
      implicit real*8  (a-h,o-z)
      parameter(maxatm=500000,maxres=200000)
C
C     All cards  ATOM or HETATM in IUNIT are
C     read to fill the coordinates and data for the protein
C
C     NOTE THAT THIS ROUTINE ASSUMES THAT BIG PDB FILES
C     ARE WELL WRITTEN, THAT IS, THAT RESIDUE NUMBER
C     IS WRITTEN WITH FORMAT I5
C
      integer iunit, numatm, numres
      dimension xr(*),yr(*),zr(*)
      dimension nfa(*),nla(*),idres(*),idatm(*)
      dimension indx(10),rx(10)
      character*4 ta(*)
      character*3 resnam(*) 
      character*80 line
C
      character head*4,tdum*3
C
      numatm=0
      numres=0
      iat=0
      ires=0
      idum0=0
      jres=0
C
      rewind(iunit)
 10   read(iunit,'(A4)',end=20,err=666) head
C
      IF (( head .eq. 'ATOM') .OR. 
     +         ( head .eq. 'HETA'))  THEN 

          backspace(iunit)
          iat=iat+1
          IF ( iat .gt. maxatm ) THEN 
               WRITE(6,'(''Reading Unit'',I3)') IUNIT  
               READ (IUNIT,'(A)') line
               WRITE(6,'(''Last line:'',/,A)') line   
               WRITE(6,'(''MAXATM='',I6)') maxatm 
               STOP '*** Too many atoms in PDB'
          ENDIF 
          read(iunit,'(6X,I5,1X,A4,1X,A3,1X,I5,4X,3F8.3)'
     &    ,err=666,end=666)
     & jat,ta(iat),tdum,jres,xr(iat),yr(iat),zr(iat)
          Idatm(iat)=jat
          IF ( idum0 .ne. jres) THEN
            ires=ires+1
            nfa(ires)=iat
            resnam(ires)=tdum
            idres(ires)=jres
            idum0=jres
          ENDIF
C
      ENDIF
C
      GOTO 10
C
 20   CONTINUE
C
      numatm=iat
      numres=ires
      DO ires=1,numres-1
        nla(ires)=nfa(ires+1)-1
      ENDDO
      nla(numres)=numatm
C
      IF (NUMATM .eq. 0) THEN 
        WRITE(6,*) 'Zero Atoms in PDB file UNIT= ',Iunit
        STOP
      ENDIF
C
      RETURN
C
 666  WRITE(6,*) 'PROBLEMS READING PDB FILE IN UNIT=', Iunit
      STOP
      END
C
C
C
C"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
      LOGICAL function EOF (LINE,IUNIT)
C     =================================      
C
C"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
C
      implicit real*8 (a-h,o-z)
      character*(*)  LINE
      READ (IUNIT,'(a)',ERR=100,END=100)  LINE
      EOF = .false.
      RETURN
 100  CONTINUE
      EOF = .true.
      LINE=' '
      RETURN
      END
C"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
C                                                                       
      SUBROUTINE UCASE (STRNG)                                          
C    ==========================                                         
C                                                                       
C"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
C     CONVERTS LOWER CASE CHARACTERS TO UPPER CASE                       
C                                                                       
      CHARACTER*(*) STRNG                                               
C                                                                       
C      write(7,'(A)') STRNG                                             
C                                                                       
      DO 10 I = 1,LEN(STRNG)                                            
         IF ( ( STRNG(I:I).GE. 'a' )  .AND.                             
     +        ( STRNG(I:I) .LE. 'z' ) )                                 
     +   STRNG(I:I) = CHAR(ICHAR(STRNG(I:I)) + ICHAR('A') - ICHAR('a')) 
 10   CONTINUE                                                          
      RETURN                                                            
      END                                                               
C
C"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
      function iwhichRES (IAT,NRES,NFA,NLA)
C     =====================================      
C
C"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
C
      implicit real*8 (a-h,o-z)
      dimension NFA(*),NLA(*)
      DO IRES=1,NRES
         IF (( NFA(IRES) .le. IAT) .and. (IAT .le. NLA(IRES))) THEN
            iwhichRES=IRES
            RETURN
         ENDIF
      ENDDO
      iwhichRES=0
      RETURN
      END
C"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

