      program FIXQTOP 
c
c     FIXING Atomic Charges in top file 
c
      parameter (maxatm=250000,maxres=125000)
      parameter(camber=18.2223D0)
      parameter(mbond=500000) 
C
C     PDB 1 (interface arrays....)
C
      dimension x1(maxatm),y1(maxatm),z1(maxatm),q1(maxatm),r1(maxatm)
      dimension nfa1(maxres),nla1(maxres)
      dimension Idatm1(maxatm),Idres1(maxres)
      character resnam1(maxres)*3
      character ta1(maxatm)*4
c
      dimension idatmq(maxatm)
c
      character*4 ta(maxatm),symb
      character*4 resnam(maxres)
      character*256 arg,fpqr,topin,topout,inpfile
      dimension q(maxatm),radii(maxatm)
      dimension ires(maxres)
      dimension rdummy(maxatm)
      character*256 line
      logical EOF,inrad,OK,check
c
      nat=0
      nres=0
      inrad=.true.
      iarg=0
      do i=1,maxatm
         x1(i)=0.0
         y1(i)=0.0
         z1(i)=0.0
         q1(i)=0.0
         r1(i)=0.0
      enddo
c
      narg= iargc()
      if (narg.lt.6) THEN
        write(6,'(5x,a)')
     &  'Usage: fixpqrtop -pqr [PDB pqr file] -topin [ parm file ]'
        write(6,'(5x,a)')
     &  '                 -topout [ parm file ]'
        stop
      endif
C
      DO WHILE (iarg .lt. narg)
         iarg = iarg + 1
         call getarg(iarg,arg)
         if (arg .eq. '-pqr') then
            iarg = iarg + 1
            call getarg(iarg,fpqr)
         else if (arg .eq. '-topin') then
            iarg = iarg + 1
            call getarg(iarg,topin)
         else if (arg .eq. '-topout') then
            iarg = iarg + 1
            call getarg(iarg,topout)
         endif
      ENDDO
C
      write(6,*)
      write(6,'(5x,''FIXPQRTOP'')')
      write(6,'(5x,''=========='')')
c
      inpfile=fpqr
      OPEN(10,FILE=inpfile,STATUS='UNKNOWN',ERR=666)
      write(6,'(5x,''Opening PDB pqr file '',A)') fpqr
c
      CALL RDPDBQ (10,NAT1,NRES1,NFA1,NLA1,TA1,RESNAM1,
     &  IDATM1,IDRES1,X1,Y1,Z1,Q1,R1)
      CLOSE(10)
C
      Rtot=0.0
      DO I=1,NAT1
        Rtot=Rtot+ABS(R1(i))
      ENDDO
      IF ( Rtot .lt. 0.001 ) inrad=.false.
C
      write(6,'(5x,''Read data in pqr for'',I4,'' atoms'')') NAT1
      write(6,'(5x,''and '',I4,'' residues'')') NRES1
      IF (.not. inrad) write(6,'(5x,''Radii not found in '',A)') fpqr
C
c
      inpfile=topin
      OPEN(10,FILE=inpfile,STATUS='UNKNOWN',ERR=666)
      write(6,'(5x,''Opening topology file '',A)') topin
C
      DO WHILE (.not. EOF(line,10) )
         IF (INDEX(line,'FLAG POINTERS') .ne. 0) THEN 
            read(10,*)
            read(10,'(10I8)') NATOM, NTYPES, NBONH,  MBONA,  NTHETH, 
     c      MTHETA, NPHIH,    MPHIA,  NHPARM, NPARM,  NNB, NRES,
     c      NBONA,    NTHETA, NPHIA,  NUMBND, NUMANG, NPTRA,
     c      NATYP,    NPHB,   IFPERT, NBPER,  NGPER,  NDPER,
     c      MBPER,    MGPER,  MDPER,  IFBOX,  NMXRS,  IFCAP,
     c      NUMEXTRA, NCOPY
            nat=NATOM
         ELSE IF (INDEX(line,'FLAG ATOM_NAME') .ne. 0) THEN 
            read(10,*) 
            read(10,'(20A4)') (ta(i),i=1,nat)
         ELSE IF (INDEX(line,'FLAG RESIDUE_LABEL') .ne. 0) THEN 
            read(10,*) 
            read(10,'(20A4)') (resnam(i),i=1,nres)
         ELSE IF (INDEX(line,'FLAG RESIDUE_POINTER') .ne. 0) THEN 
            read(10,*) 
            read(10,'(10I8)') (ires(i),i=1,nres)
         ELSE IF (INDEX(line,'FLAG CHARGE') .ne. 0) THEN 
            read(10,*) 
            read(10,'(5E16.8)') (q(i),i=1,nat)
         ELSE IF (INDEX(line,'FLAG RADII') .ne. 0) THEN 
            read(10,*) 
            read(10,'(5E16.8)') (radii(i),i=1,nat)
         ENDIF
      ENDDO
      write(6,'(5x,''Read data in top for'',I6,'' atoms'')') NAT
      write(6,'(5x,''                    '',I6,'' res'')') NRES
C
C
C    En funcion del ID atomico...cambiamos la carga
C
      qtot=0.0
      DO I=1,NAT
        qtot=qtot+q(i)
      ENDDO
      qtot=qtot/camber
C
      q0zone=0.0
      q1zone=0.0
      DO I=1,NAT
        q0zone=q0zone+q( idatm1(I) )
        q1zone=q1zone+q1( I )
      ENDDO
      q0zone=q0zone/camber
C
      write(6,'(5x,''Total charge in top'',F10.4)') qtot
      write(6,'(5x,''Charge of selected atoms in top'',F10.4)') q0zone
      write(6,'(5x,''Charge of atoms in PDB pqr '',F10.4)') q1zone
C
C     Charge rearrangement
C
      IF ( ABS(q1zone - q0zone) .gt. 0.001 ) THEN 
C
C       Find out what atoms in residues contaning the selected zone can be fixed 
C
        nqfix=0
        DO I=1,NRES1 
           IF ( I .lt. nres) THEN
              nfa=ires( idres1(I) )
              nla=ires( idres1(I)+1 ) - 1
           ELSE
             nfa=ires( idres1(I) )
             nla=nat 
           ENDIF
           DO J=nfa,nla
              OK=.true.
              DO K=1,NAT1
                 IF (J .eq. IDATM1(K)) OK=.false.
              ENDDO
              IF ( OK) THEN
                  nqfix=nqfix+1
                  idatmQ(nqfix)=J
              ENDIF
           ENDDO
        ENDDO
C
        q0fix=0.0
        DO I=1,nqfix
           q0fix=q0fix+q( idatmQ(I) )
        ENDDO
        q0fix=q0fix/camber
C
        write(6,'(5x,''Charge difference'',F10.4,'' redistributed in'',
     &  I4,'' atoms'')') ABS ( q1zone - q0zone), nqfix
        write(6,'(5x,''having originally '',F10.4,'' charge'')') q0fix
C
        qchange= - ( (q1zone - q0zone ) / float (nqfix) ) * camber
        write(6,'(5x,''Atom charges modified by '',F10.4)') 
     &  qchange/camber
C
        DO I=1,nqfix
            q( idatmQ(I) )= q( idatmQ(I) ) + qchange
        ENDDO
C
      ENDIF
C
C     Assigning new charges
C
      DO I=1,nat1
            q( idatm1(I) )= q1(I)*camber
      ENDDO
C
C     Checking total charge
C
      qtot=0.0
      DO I=1,NAT
        qtot=qtot+q(i)
      ENDDO
      qtot=qtot/camber
C
      write(6,'(5x,''Total charge in top after changes'',F10.4)') qtot
C
C     Assigning new chargesa (if relevant) 
C
      IF ( INRAD ) THEN 
C
      DO I=1,nat1
            radii( idatm1(I) )= r1(I)
      ENDDO
C
      ENDIF
C
C     Print out new topology
C
      rewind(10)
C
      OPEN(12,FILE=topout,STATUS='UNKNOWN')
      write(6,'(5x,''Writing new topology file '',A)') topout
C
      DO WHILE (.not. EOF(line,10) )
         IF (INDEX(line,'FLAG CHARGE') .ne. 0) THEN 
            write(12,'(A80)') line
            read(10,'(A)')  line 
            write(12,'(A80)') line
            write(12,'(5E16.8)') (q(i),i=1,nat)
            read(10,'(5E16.8)') (rdummy(i),i=1,nat)
         ELSE IF ((INRAD) .and.(INDEX(line,'FLAG RADII') .ne. 0)) THEN 
            write(12,'(A80)') line
            read(10,'(A)')  line 
            write(12,'(A80)') line
            write(12,'(5E16.8)') (radii(i),i=1,nat)
            read(10,'(5E16.8)') (rdummy(i),i=1,nat)
         ELSE
            write(12,'(A80)') line
         ENDIF
      ENDDO
C
      CLOSE(10)
      CLOSE(12)
C
      stop
 666  backspace(5)
      print*,'Problems opening ',inpfile
      stop
      end
C
C"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
      LOGICAL function EOF (LINE,IUNIT)
C     =================================      
C
C"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
C
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
      SUBROUTINE RDPDBQ   (IUNIT,NUMATM,NUMRES,NFA,NLA,TA,RESNAM,
     &                     IDATM,IDRES,X,Y,Z,Q,R)
C     ===========================================================
C
C"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
      parameter(maxatm=10000,maxres=3000)
C
C     All cards  ATOM or HETATM in IUNIT are
C     read to fill the Q array and data for the protein
C
C     NOTE THAT THIS ROUTINE ASSUMES THAT BIG PDB FILES
C     ARE WELL WRITTEN, THAT IS, THAT RESIDUE NUMBER
C     IS WRITTEN WITH FORMAT I5
C
      integer iunit, numatm, numres
      dimension x(*),y(*),z(*),q(*),r(*)
      dimension nfa(*),nla(*),idres(*),idatm(*)
      character*4 ta(*)
      character*3 resnam(*) 
      character*256 line
      dimension idummy(10),rdummy(10)
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
          read(iunit,'(6X,I5,1X,A4,1X,A3,1X,I5,4X)'
     &    ,err=666,end=666)
     &    jat,ta(iat),tdum,jres
          Idatm(iat)=jat
C
          backspace(iunit)
          read(iunit,'(a)') line 
          NR=0
          L=len(line)
          DO ii=30,L
            IF ((line(ii:ii) .ge. 'a').and.(line(ii:ii) .le. 'z'))
     &      line(ii:ii)=' ' 
            IF ((line(ii:ii) .ge. 'A').and.(line(ii:ii) .le. 'Z'))
     &      line(ii:ii)=' ' 
          ENDDO
          call rdlinea (line(30:L),idummy,NI,rdummy,NR,NX,2)
          IF (NR .lt. 4) THEN 
              PRINT*,line
              PRINT*,rdummy(1)
              PRINT*,'PROBLEMS WHILE READING CHARGES'
              STOP
          ENDIF
          x(iat)=rdummy(1)
          y(iat)=rdummy(2)
          z(iat)=rdummy(3)
          q(iat)=rdummy(4)
          IF (nr .EQ. 5) r(iat)=rdummy(5)
C
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
      RETURN
C
 666  WRITE(6,*) 'PROBLEMS READING PDBQ FILE IN UNIT=', Iunit
      STOP
      END
C
C"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
C
      SUBROUTINE RDLINEA(line,INDX,NI,RX,NR,NX,IOP)
C     ============================================
C
C
C"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
      parameter (maxnum=250)
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
      character line*(*)
      character var*512
      character*1 ca,cc,cb
      logical NUMI,NUMR,CHECK
      dimension INDX(*),RX(*)
      dimension ibra(maxnum),jbra(maxnum),ifirst(maxnum)
      dimension isecond(maxnum),itemp(maxnum),mpos(maxnum)
C
C
      DO I=1,maxnum
         ibra(I)=0
         jbra(I)=0
         ifirst(I)=0
         isecond(I)=0
         mpos(I)=0
         itemp(I)=0
      ENDDO
C
      NINTER=0
C
      NUMI=.false.
      IF (( IOP .eq. 1) .or. (IOP .eq. 3)) NUMI=.true.
      NUMR=.false.
      IF (( IOP .eq. 2) .or. (IOP .eq. 3)) NUMR=.true.
C
      L=LEN(line)
      I=L
      CHECK=(line(I:I) .ge. '0' ) .and. (line(I:I) .le. '9') 
      DO WHILE ( ( I .gt. 1 ) .and. ( .NOT. CHECK ) )
         I=I-1
         CHECK=(line(I:I) .ge. '0' ) .and. (line(I:I) .le. '9') 
      ENDDO
      L=I+2
C
      var(1:1)=' '
      var(2:I+1)=line(1:I)
      var(L:L)=' '
c     print*,'L=',L
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
     +  .and. ((ca .eq. ',') .or. (ca .eq. ' ')) .and. 
C
     +  ( ((cc .ge. '0') .and. (cc .le.'9')) .or. 
     +  ((cc .eq. '.') .or.  (cc .eq.' ') .or. (cc .eq. '-')) )  ) THEN 
           nleft=nleft+1
           ibra(nleft)=I
        ENDIF
        IF ( ((cb .ge. '0') .and. (cb .le. '9')) 
     +      .and. ((cc .eq. ',') .or. ( cc .eq. ' '))  )  THEN
           nright=nright+1
           jbra(nright)=I
           VAR(jc:jc)=' '
        ENDIF
C
        IF ( NUMI .and. ((cb .ge. '0') .and. (cb .le. '9')) 
     +      .and. ((cc .eq. '-') .or. ( cc .eq. ':'))  )  THEN
           VAR(jc:jc)=' ' 
           nright=nright+1
           jbra(nright)=I
           ninter=ninter+1
           IFIRST(ninter)=nright
           ISECOND(ninter)=nright+1
c          print*,'Interval', ninter,nright
c          print*,'Linea',var(1:L)
        ENDIF
C
      ENDDO
C
      IF ( nright .ne. nleft) THEN
        print*, nright,nleft
        write(6,'(2X,''Failure while reading numbers from '',A256)')
     +  line     
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
           MPOS(NI)=I
           READ(var(ipos:jpos),*,err=666,end=666) INDX(NI)
        ENDIF
      ENDDO
C
      IF ( NINTER .eq. 0) RETURN
C
      DO I=1,NI
          ITEMP(I)=INDX(I)
      ENDDO
      NTEMP=NI
C
      NI=0
      DO J=1,NTEMP-1
C
      ICHECK=0
      DO I=1,NINTER
         IF ( (IFIRST(I)  .eq. MPOS(J)) .and.
     &        (ISECOND(I) .eq. MPOS(J+1)) )   ICHECK=I
      ENDDO
C
      IF ( ICHECK .ne. 0 ) THEN 
C
          DO K=ITEMP(J),ITEMP(J+1)-1
              NI=NI+1
              INDX(NI)=K
          ENDDO
C
      ELSE
C
         NI=NI+1
         INDX(NI)=ITEMP(J)
C
      ENDIF
C
      ENDDO
C
      NI=NI+1
      INDX(NI)=ITEMP(NTEMP)
C
      RETURN
 666  IOP=-1
      RETURN
C
      END 
