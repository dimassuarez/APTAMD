      PROGRAM pdbcrd  
C
      implicit real*8 (a-h,o-z)
C     D. Suarez 
C     PSU July/1999
C
      parameter(maxatm=500000)
      parameter(maxres=200000)
      parameter(f1=1.00,f2=0.00)
C
C----------------------------------------------------------
C    Generation of a PDB file without SOLVENT
C    from an AMBER PDB file 
C
C    Input: PDB File (Unit 5)
C          
C    Output: PDB File (Unit 6)
C
C----------------------------------------------------------
C     Variables for Atoms-Hetatm 
C
      dimension xr(maxatm),yr(maxatm),zr(maxatm)
      dimension box(6)
      dimension idummy(10),rdummy(10)
      character finp*40,fout*40
      character line*80,head*4,tdum*3
      character title*80
      logical check,Water ,Lbox 
C
C     Reading Input file 
C
      ititl=0
      iat=0
      Lbox=.false.
      box=0.d0
C
 10   read(5,'(A4)',end=20,err=666) head
C
C      Reading TITLE and ATOMS 
C      lines from Input PDB File 
C
      IF (( head .eq. 'TITL') .AND. ( ititl .ne. 1)) THEN
          backspace(5) 
          read(5,'(a80)') line
          ititl=1
          title=line
          write(6,'(A80)') title
      ELSE IF (( head .eq. 'BOX ') .OR. ( head .eq. 'CRYS'))  THEN 
          backspace(5)
          read(5,'(A)') line 
          call rdlinea(line,idummy,ni,rdummy,nr,nx,2) 
          IF ( nr .eq. 6 ) THEN
              DO i=1,nr
                 box(i)=rdummy(i)
              ENDDO
          ELSE
              print*,'Number of box params:',nr
              STOP 'Problem reading box info'
          ENDIF
          lbox=.true.
C         print*, box
      ELSE IF (( head .eq. 'ATOM') .OR. 
     +         ( head .eq. 'HETA'))  THEN 
          backspace(5)
          iat=iat+1
          read(5,'(30X,3F8.3)',err=666,end=777)
     +    xr(iat),yr(iat),zr(iat)
       ENDIF
       GOTO 10
C
 20    CONTINUE
C
      numatm=iat
C
      write(6,*)   
      write(6,'(I6)') numatm 
      write(6,'(6F12.7)') (xr(i),yr(i),zr(i),i=1,numatm)
      IF ( LBOX ) write(6,'(6F12.7)') box
      STOP
C
 666  backspace(5)
      read(5,'(a80)') line
      print*, line
      Stop 'Error reading this line !'
 777  Stop 'Unexpected End of File'
C
      END
C
C     Esta rutina no esta igual en todos los FORTRAN
C     Si hay un fallo en algun codigo en lectura de numeros
C     actualizarla por esta
C
      SUBROUTINE RDLINEA(var,INDX,NI,RX,NR,NX,IOP)
      implicit real*8 (a-h,o-z)
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
      integer Aleft
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
C     Active tags
      Aleft=0       
C
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
           if (Aleft .eq.0) then 
              nleft=nleft+1
              ibra(nleft)=I
              Aleft=1
           endif
        ENDIF
C
        IF ( ( ((ca .ge. '0') .and. (ca .le. '9')) .or.
     +          (ca .eq. '-')  .or. (ca .eq. '.')  .or. 
     +          (ca .eq. ' ') )  .and. 
     +          ((cb .ge. '0') .and. (cb .le. '9')) 
     +          .and. ( cc .eq. ' ') )  THEN
           if (Aleft .eq.1) then
              nright=nright+1
              jbra(nright)=I
              Aleft=0
           endif
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
