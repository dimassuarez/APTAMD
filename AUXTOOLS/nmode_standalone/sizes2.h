
!  --- sizes for normal mode analysis program:

!    MAXATOM is maximum number of atoms
!    MAXINT  is maximum number of internal coordinates on which
!            projections are to be made
!    MAXVEC  is the maximum number of normal mode eigenvectors that
!            will be used
!    MEMDRV  is a general workspace; program will complain if this
!            is too small and tell you how large to make it.

parameter (memdrv=20000000)
parameter (maxatom=25000)
parameter (maxint=1000)
parameter (maxvec=1000)

!     ----- SET THE LIMITS OF SOME ARRAY BOUNDS -----

parameter (maxdih = 30000)
parameter (maxdia = 24000)
parameter (maxinb = 48000)
parameter (maxbon = 24000)
parameter (maxbnh = 24000)
parameter (maxang = 24000)
parameter (maxanh = 24000)

