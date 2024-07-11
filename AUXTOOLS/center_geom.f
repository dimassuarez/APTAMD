      PROGRAM Center_Geom
      implicit double precision (a-h,o-z)
C
C     D. Suarez 
C
      parameter(maxatm=500000,maxres=200000,maxhet=10000)
      parameter(maxcen=500)
      parameter(f1=1.00,f2=0.00)
C
C----------------------------------------------------------
C    Generation of a PDB file with coordinates in
C    the COM reference system or  in a given center
C    from a complete PDB file 
C
C----------------------------------------------------------
C     Variables for Atoms-Hetatm 
C
      dimension xr(maxatm),yr(maxatm),zr(maxatm)
      dimension icenter(maxcen)
      dimension idummy(10),rdummy(10)
      character resnam(maxres)*3
      character ta(maxatm)*4,t2*2,t4*4
C
      character finp*256,fout*256,arg*256
      character line*256,head*4,tdum*3
      character title*80
      logical water,centro,xyzcen,noatm 
C
C    Beginning, Okay!  
C
      noatm=.false.
      xyzcen=.false.
      narg=iargc()
      if ( narg .le. 3 ) then
         print*,'Usage; center [options]  -i input_pdb  -o  output_pdb'
         print*,'   options: '
         print*,'   -xyzcen "0.0,0.0,0.0" new center coordinates '
         print*,'           if not provided, COM by default '
         print*,'   -noatm  Only Geom measures are printed out '
         stop
      endif
      do iarg=1,narg
         call getarg(iarg,arg)
         if ( arg .eq. '-noatm') then
            noatm=.true.
         elseif ( arg .eq. '-xyzcen') then
            xyzcen=.true.
            call getarg(iarg+1,line)
            nr=0
            ni=0
            call RDLINEA(line,idummy,NI,rdummy,NR,NX,2)
            if ( nr .ne. 3 ) then
               print*,'Error while reading -xyzcen data'
               stop
            endif
            xrinput=rdummy(1)
            yrinput=rdummy(2)
            zrinput=rdummy(3)
         elseif ( arg .eq. '-i') then
                 call getarg(iarg+1,finp)
         elseif ( arg .eq. '-o') then
                 call getarg(iarg+1,fout)
         endif
      enddo
C
      open(10,file=finp,status='unknown')
      open(12,file=fout,status='unknown')
C
C     Reading Input file 
C
      ititl=0
      ires=0
      iat=0
      idum0=0
      ncen=0
      centro=.false.
C
 10   read(10,'(A4)',end=20,err=666) head
C
C      Reading TITLE and ATOMS 
C      lines from Input PDB File 
C
      IF (( head .eq. 'TITL') .AND. ( ititl .ne. 1)) THEN
          backspace(10) 
          read(10,'(a80)') line
          ititl=1
          title=line
      ELSE IF ( head .eq. 'CENT') THEN
          backspace(10) 
          read(10,'(a80)') line
          call RDLINEA(line,idummy,NI,rdummy,NR,NX,1)
          IF ( NI .eq. 1 ) THEN 
            ncen=ncen+1
            icenter(ncen)=idummy(1)
          ELSE IF ( NI .eq. 2) THEN
             DO I=idummy(1),idummy(2)
               ncen=ncen+1
               icenter(ncen)=I
             ENDDO
          ELSE
              print*,LINE
              print*,NI   
              STOP 'PROBLEM READING CENTERs'
          ENDIF
          centro=.true.
      ELSE IF ( head .eq. 'XYZC') THEN
          backspace(10) 
          read(10,'(a80)') line
          read(line(7:256),*) xrinput,yrinput,zrinput
          xyzcen=.true.
      ELSE IF (( head .eq. 'ATOM') .OR. 
     +         ( head .eq. 'HETA'))  THEN 
          backspace(10)
          iat=iat+1
          read(10,'(12X,A4,1X,A3,2X,I4,4X,3F8.3)',err=666,end=777)
     +    ta(iat),tdum,jres,xr(iat),yr(iat),zr(iat)
       ENDIF
       GOTO 10
C
 20    CONTINUE
C
       IF (ititl .ne. 0 ) write(12,'(A80)') title
       numatm=iat
       numres=ires
C
C     Computing CM and Molecular Mass 
C
      CALL rescen(numatm,xr,yr,zr,ta,xrcm,yrcm,zrcm,Rmass)     
      write(12,'(''COMMENT  NRES = '',I6,'' NAT='',I6)') numres,numatm
      write(12,'(''COMMENT  COM  = '',3F12.6)') xrcm,yrcm,zrcm
      write(12,'(''COMMENT  MASS = '',F12.2)') Rmass
      IF (( .not. centro) .and. ( .not. xyzcen) ) THEN 
         xrc=xrcm
         yrc=yrcm
         zrc=zrcm
         write(12,'(''COMMENT  CENTERED ON COM'' )') 
      ELSE IF (centro) THEN
         xrc=0.0
         yrc=0.0
         zrc=0.0
         do i=1,ncen
c        write(12,'(''COMMENT  CENTERED ON ATOM = '',I6,1X,A3)') 
c    &       Icenter(i),ta(Icenter(i))
             xrc=xrc+xr(Icenter(i))
             yrc=yrc+yr(Icenter(i))
             zrc=zrc+zr(Icenter(i))
         enddo
         xrc=xrc/float(ncen)
         yrc=yrc/float(ncen)
         zrc=zrc/float(ncen)
         write(12,'(''COMMENT  CENTERED ON XYZ = '',3F12.6)') 
     &   xrc,yrc,zrc
         distcom=sqrt ((xrc-xrcm)**2+(yrc-yrcm)**2+(zrc-zrcm)**2 )
         write(12,'(''COMMENT  DIST TO COM     = '',3F12.6)') distcom
      ELSE IF (xyzcen) THEN
         xrc=-xrinput
         yrc=-yrinput
         zrc=-zrinput
         write(12,'(''COMMENT  CENTERED ON XYZ = '',3F12.6)') 
     &   xrinput,yrinput,zrinput
         distcom=sqrt ((xrc-xrcm)**2+(yrc-yrcm)**2+(zrc-zrcm)**2 )
         write(12,'(''COMMENT  DIST TO COM     = '',3F12.6)') distcom
      ENDIF
C
C     Transforming coordinates
C
      rmin=9999999.0
      imin=0
      radgyr=0.0 
      radmax=0.0 
      DO i=1,numatm

         rad=sqrt( (xr(i)-xrcm)**2 + 
     &   (yr(i)-yrcm)**2 + (zr(i)-zrcm)**2 )
         radgyr=radgyr+rad
         if ( rad .gt. radmax) radmax=rad

         xr(i)=xr(i)-xrc
         yr(i)=yr(i)-yrc
         zr(i)=zr(i)-zrc
         r=sqrt( xr(i)**2 + yr(i)**2 +zr(i)**2 )
         if ( r .lt. rmin) then
              rmin=r
              imin=i
         endif
      ENDDO
      radgyr=radgyr/float(numatm)
      write(12,'(''COMMENT  CLOSEST ATOM = '',I6,1X,A3)') imin,ta(imin) 
      write(12,'(''COMMENT  RADGYR  = '',F12.6)')  radgyr
      write(12,'(''COMMENT  RADMAX  = '',F12.6)')  radmax 
      write(12,'(''COMMENT  MAX XYZ = '',3F12.6)')
     &              maxval(xr),maxval(yr),maxval(zr)
      write(12,'(''COMMENT  MIN XYZ = '',3F12.6)')
     &              minval(xr),minval(yr),minval(zr)

      write(12,'(''COMMENT  XYZ DIM = '',3F12.6)')
     &maxval(xr)-minval(xr),maxval(yr)-minval(yr),maxval(zr)-minval(zr)
C
      if ( .not.  noatm ) then 
C
C     Printing Output file...
C
      rewind(10)
      iat=0
 30   read(10,'(A4)',end=40,err=666) head
C
C
      IF (( head .eq. 'ATOM') .OR. 
     +         ( head .eq. 'HETA'))  THEN 
          backspace(10)
          iat=iat+1
      read(10,'(6X,I5,1X,A4,1X,A3,2X,I4)',err=666,end=777)
     +idum,t4,tdum,jres
C
C
      write(12,'(''ATOM'',1X,I6,1X,A4,1X,A3,2X,I4,4X,3F8.3)')
     +idum,t4,tdum,jres,xr(iat),yr(iat),zr(iat)
C
       ENDIF
C
       GOTO 30
C
 40    CONTINUE
c
      endif
C
      close(10)
      close(12)
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
c
      function Water(lab)
      implicit double precision (a-h,o-z)
      logical Water
      character*3 lab,watlab(3)
      data watlab/'WAT','HOH','IP3'/
      Water=.false.
      DO I=1,3
        IF ( lab .eq. watlab(I)) THEN
         Water=.true.
        ENDIF
      ENDDO
      END
C
      SUBROUTINE rescen(nat,x,y,z,t,xrc,yrc,zrc,rm)     
      implicit double precision (a-h,o-z)
C
      dimension x(*),y(*),z(*)
      character*4 t(*)
C
      real pat(7)
      character symbol(7)*1,t1*1,t2*1,t4*4
      data pat/12.01,14.00,16.00,32.35,1.01,30.97,65.4/
      data symbol/'C','N','O','S','H','P','Z'/
C
      xrc=0.0
      yrc=0.0
      zrc=0.0
      rm=0.0
      DO i=1,nat
        t4=t(i)
        t1=t4(1:1)
        t2=t4(2:2)
        rmi=0.0
        DO j=1,6
          if ( t1 .eq. symbol(j)) then
            rmi=pat(j)
          endif
          if ( t2 .eq. symbol(j)) then
            rmi=pat(j)
          endif
        ENDDO
        rm=rm+rmi
        xrc=xrc+x(i)*rmi
        yrc=yrc+y(i)*rmi
        zrc=zrc+z(i)*rmi
      ENDDO
      xrc=xrc/rm
      yrc=yrc/rm
      zrc=zrc/rm
C
      RETURN
C
      END
C
      SUBROUTINE RDLINEA(var,INDX,NI,RX,NR,NX,IOP)
      implicit real*8(a-h,o-z)
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
      DO I=1,L
         IF (var(I:I) .eq. '_') var(I:I)=' '
      ENDDO
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
     +    ((cc .eq. '.') .and. (cb .ne.'.')) .or. 
     +     (cc .eq.' ')) )  THEN 
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
