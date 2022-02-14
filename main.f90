program dynamicsOfOQS
use f95_precision
use lapack95
use blas95
use types_m
use parameters_m
use constants_m
use overlap_m
use system_hamiltonian_m
use rdftensor_m
use functions_m
use rkf_m
use time_evolution_m
use verlet_m 

implicit none

real*8 :: InitSitesRadius(nr, nc), InitSitesVel(nr, nc) 
real*8 :: InitRad(nsites), InitVel(nsites), re, spec_ho, spec_reserv
real*8 :: mode, ren 


    invdelta = 1.d0
    krdelta = 0.d0
    do i = 1, d_el
        invdelta(i, i) = 0.d0 
        krdelta(i, i) = 1.d0
    enddo


    open (unit = 1000, file = "spec_func.dat")
    do re = 0.d0, 2000.d0, 4.d0
        spec_ho = SpecFunc(ctFreq_ho, ren_ho, re)
        write(1000, '(100F20.5)') 5.3078d0 * re, spec_ho
    enddo
    close(1000)



    if (Therm_on == .true. ) then
      call init_random_seed()
      call gauss_dist( 1, .false. , 0.25d0, raioZero, InitSitesRadius )
      call gauss_dist( 2, .true., 0.25d0, velZero, InitSitesVel )
  else
      open ( unit=301, file="therm_radius.int", status="old", access="sequential" )
      read(301, *) ( ( InitRad(i) ), i = 1, nsites )
      close(301)

      open ( unit=302, file="therm_vel.int"  ,  status="old", access="sequential" )
      read(302, *) ( ( InitVel(i) ), i = 1, nsites )
      close(302)

      k = 1
      do j = 1, nc
          do i = 1, nr
              InitSitesRadius(i, j) = InitRad(k) 
              InitSitesVel(i, j) = InitVel(k)
              k = k + 1
          enddo
      enddo
      InitSitesRadius = InitSitesRadius*1.d-9
  endif

  
  if (FixSites == .true.) then
    InitSitesRadius = raioZero
    InitSitesVel = 0.d0 
  endif
  


   call System_Dynamics(InitSitesRadius, InitSitesVel)






end program DynamicsOfOQS
