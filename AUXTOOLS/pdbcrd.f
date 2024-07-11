      PROGRAM Reorder 
C
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
      dimension nfa(maxres),nla(maxres),box(3)
      character*3 resnam(maxres)
      character*4 ta(maxatm)
      character finp*40,fout*40
      character line*80,head*4,tdum*3
      character title*80
      logical check,Lbox 
C
C     Reading Input file 
C
      ititl=0
      iat=0
      ires=0
      idum0=0
      Lbox=.false.
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
      ELSE IF ( head .eq. 'BOX ')  THEN 
          backspace(5)
          read(5,'(A)') line 
          read(line(4:80),*) box
          lbox=.true.
C         print*, box
      ELSE IF (( head .eq. 'ATOM') .OR. 
     +         ( head .eq. 'HETA'))  THEN 
          backspace(5)
          iat=iat+1
          read(5,'(12X,A4,1X,A3,1X,I5,4X,3F8.3)',err=666,end=777)
     +    ta(iat),tdum,jres,xr(iat),yr(iat),zr(iat)
          IF ( idum0 .ne. jres) THEN
            ires=ires+1
            nfa(ires)=iat
            resnam(ires)=tdum
            idum0=jres
          ENDIF
       ENDIF
       GOTO 10
C
 20    CONTINUE
C
      numatm=iat
      numres=ires
C
      DO ires=1,numres-1
        nla(ires)=nfa(ires+1)-1
      ENDDO
C
      nla(numres)=numatm
C
      write(6,*)   
      write(6,'(I6)') numatm 
      write(6,'(6F12.7)') (xr(i),yr(i),zr(i),i=1,numatm)
      IF ( LBOX ) write(6,'(6F12.7)') box
      STOP
C
 666  backspace(10)
      read(5,'(a80)') line
      print*, line
      Stop 'Error reading this line !'
 777  Stop 'Unexpected End of File'
C
      END
C
