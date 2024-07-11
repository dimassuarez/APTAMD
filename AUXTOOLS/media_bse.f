      program media 
      implicit double precision (a-h,o-z)
      parameter (mxpts=1000000,mxcol=20)
      dimension V(mxpts),X(mxpts,mxcol),XMED(mxcol)
      dimension VBSE(mxpts),ISBSE(mxpts)
      dimension idummy(mxcol),rdummy(mxcol)
      character*512 linea
C
      n=0
      l=0
      NR0=0
      NR=0
      NI=0
      DO I=1,256
        linea(I:I)=' '
      ENDDO
      IOP=3
C
 10   read(5,'(a)',end=20,err=20) linea 
      l=l+1
      if ( index(linea,'#') .ne. 0) goto 10
      CALL RDLINEA(linea,idummy,NI,rdummy,NR,NT,IOP)
      IF ( IOP .eq . -1) STOP 'PROBLEM READING DATA'
      IF (((NR .ne. 0 ) .or. (NI .ne.0)) .and. ( n .eq. 0 )) THEN 
       NR0=NR+NI
      ELSE IF ( NR0 .ne. NR+NI) THEN
       WRITE(6,'(2X,''Row '',I5,'' has '',I3,'' REAL numbers'')')l,NR
       WRITE(6,'(2X,''Row '',I5,'' has '',I3,'' REAL numbers'')')l-1,NR0
       IF ( n .gt. 1 ) GOTO 10
      ENDIF
C
      IF ( NR+NI .gt. mxcol) STOP 'Too many columns' 
C
      IF ( NR .ne. 0 ) THEN 
       n=n+1
       DO J=1,NR
          X(n,J)=rdummy(J)
       ENDDO
      ENDIF
      IF ( NI .ne. 0 )  THEN 
       IF ( NR .eq. 0) n=n+1
       DO J=1,NI
          X(n,NR+J)=float(idummy(J))
       ENDDO
      ENDIF
      goto 10
C
 20   CONTINUE 
      write(6,*)
      write(6,*) 'Leidas ',l,' lineas'
      write(6,*) 'Solo ',n,' contenian datos'
      write(6,*)
C
      DO J=1,NR0
C
      DO I=1,n
      V(I)=X(I,J)
      ENDDO
      call average(n,v,vmed,vsd,vmax,vmin,vrms)
C
      nbse=0
      call bse(n,nbse,v,vbse,isbse)
c
      DO I=1,nbse
          X(I,J)=VBSE(I)
      ENDDO
      print*,nbse
c
C
      write(6,*) '-----------------------------------------------------'
      write(6,'(20X,''COLUMN '',I2)') J
      write(6,*) '              MEDIA=',vmed
      write(6,*) 'DESVIACION STANDARD =',vsd   
      write(6,*) '     ERROR STANDARD = ',vsd/sqrt(float(n))
      write(6,*) 'LIMIT BLOCKED ERROR =',vbse(nbse)
      write(6,*) 'MAX=',vmax,' MIN=',vmin  
      write(6,*) 'RMS (media del cuadrado)=',vrms
C
      ENDDO
      write(6,*) '-----------------------------------------------------'
C
      STOP
      END
C
      SUBROUTINE AVERAGE(n,V,Vmed,Vsd,Vmax,Vmin,Vrms)
      implicit double precision (a-h,o-z)
      dimension V(*)
      Vmax=-99999999999999999999999999.0
      Vmin=99999999999999999999999999.0
      Vmed=0.0D0
      Vmed2=0.0D0
      DO i=1,n
         IF (v(i) .gt. vmax) vmax=v(i) 
         IF (v(i) .lt. vmin) vmin=v(i) 
         Vmed=Vmed+V(i)
         Vmed2=Vmed2+V(i)**2
      ENDDO
      Vmed=Vmed/float(n)
      Vmed2=Vmed2/float(n)
      Vsd=sqrt(Vmed2-Vmed**2)
      Vrms=sqrt(Vmed2)
      RETURN
      END                         
C
      SUBROUTINE BSE(n,nbse,v,vbse,isbse)
      implicit double precision (a-h,o-z)
      parameter (mxpts=250000)
      parameter (mxbl=4)
      dimension V(*),VBSE(*),ISBSE(*),BMED(mxpts)
C
C     Maximo tamaño de bloque
C
      mmax=n/mxbl
C    
C     Maximo numero de puntos intermedios (100)
C
      kh=mmax/100
      if  (kh .eq. 0) kh=1
C   
      jbse=0
      do kbse=1,mmax,kh
C 
C     promedio de cada bloque
C
C
       lbse=0
       i=1
       do j=kbse,n,kbse
          vmed=0.0
          do k=i,j
            vmed=vmed+v(k)
          enddo
          vmed=vmed/float(j-i+1)   
          i=j
          lbse=lbse+1
          bmed(lbse)=vmed
       enddo
C
       call average (lbse,bmed,adummy,bsd,bdummy,cdummy,eduumy)
       jbse=jbse+1
       vbse(jbse)=bsd/sqrt(float(n)/float(kbse))
       isbse(jbse)=kbse
C
      enddo
C
      nbse=jbse
      return
      end
C
      SUBROUTINE RDLINEA(var,INDX,NI,RX,NR,NX,IOP)
      implicit real*8(a-h,o-z)
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
        IF (((( cb .eq. '-') .or. ( cb  .eq. '.')  .or. ( cb .eq. '+') 
     +    .or. ((cb .ge. '0') .and. (cb .le. '9')))) 
     +  .and. (ca .eq. ' ') .and. 
C
     +  ( ((cc .ge. '0') .and. (cc .le.'9')) .or. 
     +    ((cc .eq. '.') .or.  (cc .eq.' ')) ) ) THEN 
           nleft=nleft+1
           ibra(nleft)=I
        ENDIF
C
        IF ( ( ((ca .ge. '0') .and. (ca .le. '9')) .or.
     +     (ca .eq. '-')  .or. (ca .eq. '.')  .or.  (ca .eq. '+') .or.
     +     (ca .eq. ' ') )  .and. 
     +     ((cb .ge. '0') .and. (cb .le. '9')) 
     +     .and. ( cc .eq. ' ') )  THEN
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
