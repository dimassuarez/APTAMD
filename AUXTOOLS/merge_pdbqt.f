      PROGRAM merge_pdbqt
      implicit real (a-h,o-z)
C
C     D. Suarez 
C
      parameter(maxatm=500000,maxres=200000)
C
C     PDB 1
C
      dimension x1(maxatm),y1(maxatm),z1(maxatm)
      dimension q1(maxatm),rad1(maxatm),qdock1(maxatm)
      character*3  atype1(maxatm) 
      dimension nfa1(maxres),nla1(maxres)
      dimension Idatm1(maxres),Idres1(maxres)
      character resnam1(maxres)*3
      character ta1(maxatm)*4
C
C     PDB 2
C
      dimension x2(maxatm),y2(maxatm),z2(maxatm)
      dimension nfa2(maxres),nla2(maxres)
      dimension Idatm2(maxres),Idres2(maxres)
      character resnam2(maxres)*3
      character ta2(maxatm)*4
C
      dimension imatch(maxatm)
C
      character*4 ta1_temp,ta1_check,ta2_temp,ta2_check
      character finp1*256,finp2*256,fout*256
      character line*80, head*6
      logical diffnat,eof 
C
C    Beginning, Okay!  
C
      diffnat=.false.
      finp1=' '
      finp2=' '
      CALL GETARG(1,finp1)
      CALL GETARG(2,finp2)
      CALL GETARG(3,fout)
      IF ((finp1 .eq.' ').or.(finp2.eq.' ').or.(fout.eq.' ')) THEN 
        STOP 'Usage: merge_pdbqt pdbqtinp1  pdbinp2  pdbqtout'
      ENDIF
C
C     Reading File 1
C
      OPEN (10,file=finp1,status='unknown')
      CALL RDPDB (10,NAT1,NRES1,NFA1,NLA1,TA1,RESNAM1,
     &     IDATM1,IDRES1,X1,Y1,Z1,Q1,RAD1,QDOCK1,ATYPE1)
      CLOSE(10)
      qsum=0.0
      DO i=1,nat1
        qsum=qsum+q1(i)
      ENDDO
      WRITE(6,'(''Sum of charges in PDB1 input file ='',F8.2)') qsum
C
C     Reading File 2
C
      OPEN (10,file=finp2,status='unknown')
      CALL RDPDBXYZ (10,NAT2,NRES2,NFA2,NLA2,TA2,RESNAM2,
     &                  IDATM2,IDRES2,X2,Y2,Z2)
      CLOSE(10)
C
      IF ( NAT1 .ne. NAT2 ) THEN 
         print*,'Number of atoms do not match !'
         print*,nat1,nat2
         print*,'We proceed anyway'
         diffnat=.true.
      ENDIF
C
      IF ( NRES1 .ne. NRES2 ) THEN 
         print*,'Number of residues do not match !',nres1,nres2
         print*,'We cannot proceed further '
         stop 
      ENDIF
C
C     Printing out 
C
C
      isum=0
      DO ires=1,nres1
            DO iat=nfa1(ires),nla1(ires)
            ta1_check=''
            ta1_temp=ta1(iat)
            l1=0
            do k1=1,4
               if ((ta1_temp(k1:k1) .ne. ' ') .and. 
     &             (ta1_temp(k1:k1) .ne. '')) then 
                  l1=l1+1
                  ta1_check(l1:l1)=ta1_temp(k1:k1)
               endif
            enddo
            icheck=0
C
C Fix atom names of H's with an function
C
            DO jat=nfa2(ires),nla2(ires)
                l2=0
                ta2_check=''
                ta2_temp=ta2(jat)
                do k2=1,4
                    if ((ta2_temp(k2:k2) .ne. ' ') .and. 
     &                  (ta2_temp(k2:k2) .ne. '')) then 
                      l2=l2+1
                      ta2_check(l2:l2)=ta2_temp(k2:k2)
                   endif
                enddo
                IF (ta1_check(1:l1) .eq. ta2_check(1:l2) ) icheck=jat
            ENDDO
            IF ( icheck .ne. 0) THEN 
               isum=isum+1
               imatch(iat)=icheck
            ELSE   
               imatch(iat)=0
            ENDIF
            ENDDO
      ENDDO
C
      OPEN (10,file=finp1,status='unknown')
      OPEN (12,file=fout,status='unknown')

      iat=0
      DO WHILE (.not. EOF(LINE,10) )
         HEAD=LINE(1:6)
         IF ((HEAD .eq. 'ATOM  ') .or. (HEAD .eq. 'HETATM')) THEN
            iat=iat+1
            ires=iwhichRES(iat,nres1,nfa1,nla1)
            IF ( imatch(iat) .gt. 0 ) THEN
              jat=imatch(iat)
              write(12,'(''ATOM'',
     +        2X,I5,1X,A4,1X,A3,1X,I5,4X,3F8.3,2F6.2,F10.3,A3)')
     +         Idatm1(iat),ta1(iat),resnam1(ires),Idres1(ires),
     +         x2(jat),y2(jat),z2(jat),q1(iat),rad1(iat),qdock1(iat),
     +         atype1(iat)
            ELSE
              write(12,'(''ATOM'',
     +        2X,I5,1X,A4,1X,A3,1X,I5,4X,3F8.3,2F6.2,F10.3,A3,'' x '')')
     +         Idatm1(iat),ta1(iat),resnam1(ires),Idres1(ires),
     +         x1(iat),y1(iat),z1(iat),q1(iat),rad1(iat),qdock1(iat),
     +         atype1(iat)
            ENDIF
         ELSE
             WRITE(12,'(A80)') LINE
         ENDIF
      ENDDO
      CLOSE(10)
      CLOSE(12)
C
      WRITE(6,'(''Replaced coord for # atom ='',I8)') isum
C
      STOP
C
      END
C
C
C"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
C
      SUBROUTINE RDPDB (IUNIT,NUMATM,NUMRES,NFA,NLA,TA,RESNAM,
     &                  IDATM,IDRES,XR,YR,ZR,Q,RAD,QDOCK,ATYPE)
C     =============================================================
C
C"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
      implicit real  (a-h,o-z)
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
          read(iunit,'(6X,I5,1X,A4,1X,A3,1X,I5,4X,3F8.3,2F6.2,F10.3,A3)'
     &    ,err=666,end=666)
     & jat,tadum,tdum,jres,xr(iat),yr(iat),zr(iat),
     & q(iat),rad(iat),qdock(iat),atype(iat)
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
C
C
C"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
C
      SUBROUTINE RDPDBXYZ (IUNIT,NUMATM,NUMRES,NFA,NLA,TA,RESNAM,
     &                   IDATM,IDRES,XR,YR,ZR)
C     =============================================================
C
C"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
      implicit real  (a-h,o-z)
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
C     Esta rutina no esta igual en todos los FORTRAN
C     Si hay un fallo en algun codigo en lectura de numeros
C     actualizarla por esta
C
      SUBROUTINE RDLINEA(var,INDX,NI,RX,NR,NX,IOP)
      implicit real (a-h,o-z)
C
C     Esta subrutina, asi en castellano de vez en cuando
C     esta muy guapo, intentara leer indices numericos
C     contenidos en un variable alfanumerica
C     Es un poco chapuza, pero esperamos que efectiva !
C
C     Si IOP = 1 Lectura de numeros enteros
C     Si IOP = 2 Lectura de numeros reales
C     Si IOP = 3 Lectura de numeros enteros y reales
C
C     Si IOP=-1 en la salida....chungo !
C
C     Lectura ordenada
C
      character var*(*)
      character*1 ca,cc,cb
      logical NUMI,NUMR
      dimension INDX(*),RX(*)
      dimension ibra(50),jbra(50)
C
      NUMI=.false.
      IF (( IOP .eq. 1) .or. (IOP .eq. 3)) NUMI=.true.
      NUMR=.false.
      IF (( IOP .eq. 2) .or. (IOP .eq. 3)) NUMR=.true.
C
      L=LEN(var)
      IF ( L .eq. 0 ) THEN
          IOP=-1
          RETURN
      ENDIF
      nbra=0
      nleft=0
      nright=0
      DO I=1,L
        ja=I-1
        jb=I
        jc=I+1
C
        IF (ja .gt. 0 ) THEN 
           ca=VAR(ja:ja)
        ELSE
           ca=' '
        ENDIF
        cb=VAR(jb:jb)
        IF (jc .lt.  L ) THEN 
           cc=VAR(jc:jc)
        ELSE
           cc=' '
        ENDIF
C
        IF (((( cb .eq. '-') .or. ( cb  .eq. '.')  .or.
     +  ((cb .ge. '0') .and. (cb .le. '9')))) 
     +  .and. (ca .eq. ' ') .and. 
C
     +  ( ((cc .ge. '0') .and. (cc .le.'9')) .or. 
     +    ((cc .eq. '.') .or.  (cc .eq.' ')) ) ) THEN 
           nleft=nleft+1
           ibra(nleft)=I
        ENDIF
C
        IF ( ( ((ca .ge. '0') .and. (ca .le. '9')) .or.
     +          (ca .eq. '-')  .or. (ca .eq. '.')  .or. 
     +          (ca .eq. ' ') )  .and. 
     +          ((cb .ge. '0') .and. (cb .le. '9')) 
     +          .and. ( cc .eq. ' ') )  THEN
           nright=nright+1
           jbra(nright)=I
        ENDIF
      ENDDO
C
C
      IF ( nright .ne. nleft) THEN
        print*, nright,nleft
        write(6,'(2X,''Failing while reading numbers from '',//,A256)')
     +  var      
        IOP=-1
        RETURN
      ENDIF
C
      nbra=nleft
      NX=0
      NR=0
      NI=0
C
      DO I=1,nbra
        ipos=ibra(I)
        jpos=jbra(I)
        IF ((INDEX(var(ipos:jpos),'.') .ne. 0 ) .and. (NUMR)) THEN
           NX=NX+1
           NR=NR+1
           READ(var(ipos:jpos),*,err=666,end=666) RX(NR)
        ELSE IF ( NUMI) THEN
           NX=NX+1
           NI=NI+1
           READ(var(ipos:jpos),*,err=666,end=666) INDX(NI)
        ENDIF
      ENDDO
C
      RETURN
 666  IOP=-1
      RETURN
C
      END 
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
