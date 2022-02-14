module constants_m

!	Declarar Unidades

    
    real*8, parameter  :: HALF = 0.5d0 
    real*8, parameter :: HB_ev = 6.58211915d-16 ! hcortado em eV.s
    real*8, parameter :: kboltz = 8.6173324d-5 !constante de boltzmann em ev/K
    real*8, parameter :: kboltz_j = 1.380649d-23  !constante de Boltzmann em J/K 
    real*8, parameter :: HB_ev_ps = 6.58211915d-4  !hcortado em eV.ps
    real*8, parameter :: hz_to_thz = 1.d-12 !constante para passar as frequencias de Hz para Thz
    real*8, parameter :: thz_to_hz = 1.d12 !constante para passar as frequencias de Thz para Hz
    real*8, parameter :: joule_to_ev = 6.241509343d+18
    real*8, parameter :: ev_to_joule = 1.602176565d-19
    
    complex*16, parameter :: zi = (0.0d0,1.d0)
    
    integer :: pl, grid_size, i, m, j, l, k, n, r, c, alpha, beta, gama, delta, lambda, mu, kappa, sn
    
    real*8 , parameter :: R_ZERO = 0.d0
    real*8 , parameter :: electronMass = 9.1093897d-31 !Kg
    real*8 , parameter :: hbar_ev = 6.5821192815d-16 !eV.s
    real*8 , parameter :: hbar = 1.05457172647d-34 !J.s 
    real*8 , parameter :: PI = 3.141592653589793238462643383279502884197d0
    real*8 , parameter :: SQRT2 = 1.4142135623730950488d0
    real*8, parameter :: SQRT2_INV = 1.d0 / SQRT2   
    real*8, parameter :: twohpi = 2.d0 / PI     
    


end module constants_m
