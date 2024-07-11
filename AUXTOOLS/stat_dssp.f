      program stat
      implicit double precision (a-h,o-z)
      parameter (maxvar=100000)
      character*128 record, inputfile
      character*256 line
      character*24 mark(maxvar),cmark
      character*1  sig(maxvar),csig
      integer T(maxvar),C(maxvar),H(maxvar),B(maxvar)
      integer G(maxvar),E(maxvar),TCB(maxvar)
      integer S(maxvar),HI(maxvar),PP(maxvar)
      integer CODE(maxvar)
      integer Ihist(maxvar),Idummy(maxvar)  
      logical gnuplot
C
      gnuplot=.false.
      do i=1,maxvar
        T(I)=0
        C(I)=0
        H(I)=0
        B(I)=0
        S(I)=0
        E(I)=0
        HI(I)=0
        PP(I)=0
        TCB(I)=0
        CODE(I)=0
        Ihist(I)=0
        Idummy(I)=0
      enddo
c
      inquire(file='info_for_gnuplot.dat',exist=gnuplot)
      if ( gnuplot ) then
         open(30,file='info_for_gnuplot.dat',status='unknown')
         read(30,*) ifirst,ilast
         rewind(30)
         write(30,'(''# DSSP data from '',I4,'' to '',I4)') ifirst,ilast
      endif
C
      open(10,file='traj.dat',status='unknown') 
C
      read(5,*) NVAR
C
      ifile=0
 100  read(5,'(a)',err=200,end=200) record
      ifile=ifile+1
      open(20,file=record)
      read(20,'(a)') line
      do while ( index(line,'#  RESIDUE AA STRUCTURE') .eq. 0 )
         read(20,'(a)') line
         if ( index(line,'HISTOGRAMS').ne.0 ) then 
         read(20,'(30I3)') (Idummy(j),j=1,30)
         do j=1,30
           Ihist(j)=Ihist(j)+IDummy(j) 
         enddo
        endif 
      enddo
      i=0
      do while ( I .lt. nvar)
         read(20,'(A16,A1)') cmark,csig
         if ( index(cmark,' ! ').eq.0 )then
           i=i+1
           mark(i)=cmark
           sig(i)=csig
         endif
      enddo
      close(20)
      Tsnap=0.0
      Csnap=0.0
      Hsnap=0.0
      Bsnap=0.0
      Gsnap=0.0
      Esnap=0.0
      Ssnap=0.0
      HIsnap=0.0
      PPsnap=0.0
C
C     T+C+B ->  COIL as usual custom by DSSP users
C
C     CODE 0 --> COIL
C     CODE 1 --> BEND
C     CODE 2 --> HELIX
C     CODE 3 --> STRAND
C     CODE 4 --> PPII 
C
C
      do i=1,nvar
         if ( sig(I) .eq. 'T') then 
            T(I)=T(I)+1
            Tsnap=Tsnap+1.0 
            TCB(I)=TCB(I)+1
            TCBsnap=TCBsnap+1.0 
            CODE(I)=0 
         else if ( sig(I) .eq. ' ') then
            C(I)=C(I)+1
            Csnap=Csnap+1.0 
            TCB(I)=TCB(I)+1
            TCBsnap=TCBsnap+1.0 
            CODE(I)=0 
         else if ( sig(I) .eq. 'H') then
            H(I)=H(I)+1
            Hsnap=Hsnap+1.0 
            CODE(I)=2 
         else if ( sig(I) .eq. 'B') then
            B(I)=B(I)+1
            Bsnap=Bsnap+1.0 
            TCB(I)=TCB(I)+1
            TCBsnap=TCBsnap+1.0 
            CODE(I)=0 
         else if ( sig(I) .eq. 'G') then
            G(I)=G(I)+1
            Gsnap=Gsnap+1.0 
            CODE(I)=2 
         else if ( sig(I) .eq. 'E') then
            E(I)=E(I)+1
            Esnap=Wsnap+1.0 
            CODE(I)=3
         else if ( sig(I) .eq. 'S') then
            S(I)=S(I)+1
            Ssnap=Ssnap+1.0 
            CODE(I)=1
         else if ( sig(I) .eq. 'I') then
            HI(I)=HI(I)+1
            HIsnap=HIsnap+1.0 
            CODE(I)=2
         else if ( sig(I) .eq. 'P') then
            PP(I)=PP(I)+1
            PPsnap=PPsnap+1.0 
            CODE(I)=4
         endif
         if ( sig(I) .eq. ' ') sig(I)='0'
      enddo
      Tsnap=(Tsnap/float(nvar))*100.0
      Csnap=(Csnap/float(nvar))*100.0
      Hsnap=(Hsnap/float(nvar))*100.0
      Bsnap=(Bsnap/float(nvar))*100.0
      Gsnap=(Gsnap/float(nvar))*100.0
      Esnap=(Esnap/float(nvar))*100.0
      Ssnap=(Ssnap/float(nvar))*100.0
      HIsnap=(HIsnap/float(nvar))*100.0
      PPsnap=(PPsnap/float(nvar))*100.0
      TCBnap=(TCBsnap/float(nvar))*100.0
      
      write(10,'(A10,1X,50(A1,1X))') record,(sig(I),I=1,nvar) 
      write(6,'(A20,5X,''T='',F8.2,'' C='',F8.2,'' E='',F8.2,
     &   '' H='',F8.2,'' B='',F8.2,'' G='',F8.2,
     &   '' I='',F8.2,'' S='',F8.2,'' P='',F8.2,'' TCB='',F8.2)') 
     & record,Tsnap,Csnap,Esnap,Hsnap,Bsnap,Gsnap,HIsnap,Ssnap,PPsnap,
     & TCBsnap
      if ( gnuplot ) then
          write(30,'(I6,100(I4,I2))') ifile,(J,CODE(J),J=IFirst,Ilast)
      endif
      goto 100
c
 200  continue
      write(6,'(''NRES='',I6)') NVAR
      write(6,'(''AVERAGES'')') 
      do i=1,NVAR
         write(6,'(A50,1X,''T='',F8.2,'' C='',F8.2,'' E='',F8.2,
     &   '' H='',F8.2,'' B='',F8.2,'' G='',F8.2,
     &   '' I='',F8.2,'' S='',F8.2,'' P='',F8.2,'' TCB='',F8.2)') 
     &   mark(I),float(T(I))/float(IFILE)*100.0, 
     &   float(C(I))/float(IFILE)*100.0, 
     &   float(E(I))/float(IFILE)*100.0, 
     &   float(H(I))/float(IFILE)*100.0, 
     &   float(B(I))/float(IFILE)*100.0,
     &   float(G(I))/float(IFILE)*100.0,
     &   float(HI(I))/float(IFILE)*100.0,
     &   float(S(I))/float(IFILE)*100.0,
     &   float(PP(I))/float(IFILE)*100.0,
     &   float(TCB(I))/float(IFILE)*100.0
      enddo
c     write(6,'(''Average Histogram of #Residues per Alpha Helix'')')
c     do j=1,30
c       write(6,'(I10,F10.2)') j,float(Ihist(j))/float(IFILE) 
c     enddo
      close(10)
      if (gnuplot) close(30)
      stop
      end

     
          
