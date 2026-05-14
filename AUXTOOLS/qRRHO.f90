program qRRHOvib

implicit none
!Define constants
real*8,parameter :: R=8.3144648D0, kb=1.3806503D-23 !Ideal gas constant (J/mol/K), Boltzmann constant (J/K)
real*8,parameter :: NA=6.02214179D23 !Avogadro constant
real*8,parameter :: au2eV=27.2113838D0,au2kcal_mol=627.51D0,au2kJ_mol=2625.5D0,au2J=4.359744575D-18,au2cm_1=219474.6363D0 !Hartree to various units
real*8,parameter :: cal2J=4.184D0
real*8,parameter :: wave2freq=2.99792458D10 !cm^-1 to s^-1 (Hz)
real*8,parameter :: h=6.62606896D-34 !Planck constant, in J*s
real*8,parameter :: amu2kg=1.66053878D-27
real*8,parameter :: pi=3.141592653589793D0
real*8,parameter :: b2a=0.52917720859D0 !Bohr to Angstrom
real*8,parameter :: atm2Pa=101325 !atm to Pa
integer,parameter :: maxfreq=25000

real*8 freq(maxfreq)  ! Vib freq
real*8 wavenum(maxfreq)   
real*8 svib(maxfreq)   ! RRHO vib entropy contributions
real*8 svib_q(maxfreq) ! qRRHO vib entropy contributions
real*8 tmpS,miu,miup,Bav,wei,term,prefac,T,Sfree
real*8 svib_tot,svib_q_tot

integer i,j,k,nfreq

! Input
read(5,*) nfreq
do i=1,nfreq
   read(5,*) wavenum(i)
   freq(i)=wavenum(i)*wave2freq
enddo

T=300.0   ! 300K

svib_tot=0.0d0
svib_q_tot=0.0d0


write(6,'(''  MODE     Freq(cm-1)   Svib_RRHO (eu)   Svib_qRRHO(eq) '')')
write(6,'(''========================================================'')')

do i=1,nfreq

tmpS=0
svib(i)=0.0d0
svib_q(i)=0.0d0
if (freq(i)<=0) cycle  
prefac=h*freq(i)/(kb*T)
term=exp(-h*freq(i)/(kb*T))
tmpS=R*(prefac*term/(1-term)-log(1-term)) !RRHO
svib(i)=tmpS
miu=h/(8*pi**2*freq(i))
Bav=1D-44 !kg*m^2
miup=miu*Bav/(miu+Bav)
Sfree=R*( 0.5D0+log(dsqrt(8*pi**3*miup*kb*T/h**2)) )
wei=1/(1+(100D0/wavenum(i))**4)
tmpS=wei*tmpS+(1-wei)*Sfree
svib_q(i)=tmpS

svib_tot=svib_tot+svib(i)
svib_q_tot=svib_q_tot+svib_q(i)

write(6,'(2X,I5,2X,F10.2,2X,F8.3,2X,F8.3)') i,wavenum(i),svib(i)/cal2J,svib_q(i)/cal2J 

enddo
write(6,'(''========================================================'')')
write(6,'(''  TOTAL                '',2(2X,F8.3))') svib_tot/cal2J,svib_q_tot/cal2J
write(6,'(''  CORRECTION           '',2X,F8.3)') (svib_q_tot-svib_tot)/cal2J

end



