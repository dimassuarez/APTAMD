      PROGRAM Reorder 
C
C     D. Suarez 
C     PSU July/1999
C
      parameter(maxatm=500000)
      parameter(maxres=200000)
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
      dimension nfa(maxres),nla(maxres),ncon(maxres)
      dimension f1(maxatm),f2(maxatm) 
      character*3  resnam(maxres)
      character*4  ta(maxatm)
      character*26 info(maxatm)
      character finp*40,fout*40
      character line*80,head*4,tdum*3
      character title*80,label*3
C
C     Reading Input file 
C
      ititl=0
      iat=0
      ires=0
      idum0=0
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
      ELSE IF ( head .eq. 'TER') THEN
          idum0=-99
      ELSE IF ( head .eq. 'ENDM') THEN
          idum0=-99
      ELSE IF ( head .eq. 'MODE') THEN
          idum0=-99
      ELSE IF (( head .eq. 'ATOM') .OR. 
     +         ( head .eq. 'HETA'))  THEN 
          backspace(5)
          iat=iat+1
          read(5,'(12X,A4,1X,A3,2X,I4,4X,3F8.3,A)',err=666,end=777)
     +    ta(iat),tdum,jres,xr(iat),yr(iat),zr(iat),info(iat)
          IF ( idum0 .ne. jres) THEN
            ires=ires+1
            nfa(ires)=iat
            resnam(ires)=tdum
            ncon(ires)=0
            idum0=jres
          ENDIF  
C
 15       CONTINUE
C
       ENDIF
       GOTO 10
C
 20    CONTINUE
C
      numatm=iat
      numres=ires
C     print*,numatm
C     print*,numres
C
      DO ires=1,numres-1
        nla(ires)=nfa(ires+1)-1
      ENDDO
C
C
      nla(numres)=numatm
C
C     Printing OUTPUT PDB FILE
C 
      zero=0.d0
      DO ires=1,numres
        DO I=nfa(ires),nla(ires)
        write(6,'(''ATOM'',1X,I6,1X,A4,1X,A3,1X,I5,4X,3F8.3,A26)')
     +  I,ta(I),resnam(ires),ires,xr(I),yr(I),zr(I),info(I)
        ENDDO
        IF (( ires .gt. 1 ) .and. (ires .lt. numres)) THEN 
         IF ((index(resnam(ires),'NME').ne.0) .and.
     &       (index(resnam(ires+1),'ACE').ne.0)) WRITE(6,'(''TER'')')   
        ENDIF
      ENDDO
      WRITE(6,'(''TER'')')
C
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
