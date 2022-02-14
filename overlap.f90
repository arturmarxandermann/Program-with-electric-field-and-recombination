module overlap_m

use types_m
use constants_m
use parameters_m

contains


subroutine define_sitios(InitSitesRadius, InitSitesVel, site, site_point)
!subrotina que define todos os sitios da celula OPV
    implicit none

    ! args
    real*8, intent(in) :: InitSitesRadius(nr, nc), InitSitesVel(nr, nc)
    type(quantum_site), target, allocatable, intent(out) :: site(:,:)
    type(obj_pointer), allocatable, intent(out) :: site_point(:)
    
    ! local
    real*8, parameter :: scale_nanometro = 1.d-9
    integer           :: cs, rs
    
    
    allocate(site(nr, nc))
    allocate(site_point(nsites)) !parametrizando os sitios com um único índice

    if ( SisDes == .true. ) then
        !CALCULO DOS INDICES DA AREA DE DESORDEM
        if ( ncdes /= 2 .AND. ncdes /= 4 .AND. ncdes /= 6 ) then
            print*, "CÁLCULO ERRADO DA ÁREA DE DESORDEM"
            print*, "ncdes PRECISA SER 2 OU 4 OU 6"  
            stop
        endif


        if ( ncdes == 2 ) then
            desi = ndo - nano
            desf = ndo + 3 * nano
        endif
        
        if ( ncdes == 4 ) then
            desi = ndo - 2 * nano
            desf = ndo + 4 * nano
        endif

        if ( ncdes == 6 ) then
            desi = ndo - 3 * nano
            desf = ndo + 5 * nano
        endif
        
        !SE HOUVER DESORDEM NO SISTEMA ALOCO OS VETORES QUE CONTEM OS INDECES DOS SITIOS    
        allocate(desd_r(ndes)  , source = 0 )
        allocate(desd_c(ndes)  , source = 0 )
        allocate(desa_r(ndes)  , source = 0 )
        allocate(desa_c(ndes)  , source = 0 )

        !LINHAS E COLUNAS DOS DOADORES ACEITADORES MODIFICADOS  
        !-- ESTRUTURA PARA 84 SÍTIOS --
        !desd_r = [1, 2, 3, 4, 5, 6] !DOADOR -> ACEITADOR
        !desd_c = [5, 5, 5, 6, 6, 6] !DOADOR -> ACEITADOR
        !desa_r = [5, 6, 7, 2, 3, 4] !ACEITADOR -> DOADOR
        !desa_c = [7, 7, 7, 8, 8, 8] !ACEITADOR -> DOADOR
        !-----------------------------

    endif
  
    
    !CRIANDO OS SITIOS COMO SE NAO HOUVESSE DESORDEM
    do j = 1, nr
        do i = 1, nanol  !nanol = numero de anodos em uma linha   
            site(j, i)%ovlptype = 1
            site(j, i)%mass    = siteMass  
            site(j, i)%radius0 = raioZero
            site(j, i)%radius  = InitSitesRadius(j, i)
            site(j, i)%vel     = InitSitesVel(j, i)
            site(j, i)%omega   = ( 2.d0 * hbar / ( me * (site(j, i)%radius)**2.0 ) ) * hz_to_thz   !largura a meia altura
            site(j, i)%omega0  = site(j, i)%omega * 1.d12 
            site(j, i)%vhomo   = ehomoAN
            site(j, i)%vlumo   = elumoAN
            site(j, i)%dashtype = 1
        enddo
    
        do i = iidi, ihd  !iidi = nanol + 1, ihd = nanol + ndol 
            if (i <  ihd ) then
                site(j, i)%ovlptype = 1
                site(j, i)%mass    = siteMass  
                site(j, i)%radius0 = raioZero
                site(j, i)%radius  = InitSitesRadius(j, i)
                site(j, i)%vel     = InitSitesVel(j, i)  
                site(j, i)%omega   = ( 2.d0 * hbar / ( me * (site(j, i)%radius)**2.0 ) ) * hz_to_thz   !largura a meia altura
                site(j, i)%omega0  = site(j, i)%omega * 1.d12 
                site(j, i)%vhomo   = ehomoD 
                site(j, i)%vlumo   = elumoD
                site(j, i)%dashtype = 2
            endif 

            if (i == ihd ) then
                site(j, i)%ovlptype = 1
                site(j, i)%mass    = siteMass  
                site(j, i)%radius0 = raioZero
                site(j, i)%radius  = InitSitesRadius(j, i)
                site(j, i)%vel     = InitSitesVel(j, i)  
                site(j, i)%omega   = ( 2.d0 * hbar / ( me * (site(j, i)%radius)**2.0 ) ) * hz_to_thz   !largura a meia altura
                site(j, i)%omega0  = site(j, i)%omega * 1.d12 
                site(j, i)%vhomo   = ehomoD
                site(j, i)%vlumo   = elumoD
                site(j, i)%dashtype = 2
            endif 
        enddo

        do i = ihd, nc
          
            if ( i == iha ) then !iha = nanol + ndol + 1
                site(j, i)%ovlptype = 2
                site(j, i)%mass    = siteMass  
                site(j, i)%radius0 = raioZero
                site(j, i)%radius  = InitSitesRadius(j, i)
                site(j, i)%vel     = InitSitesVel(j, i)  
                site(j, i)%omega   = ( 2.d0 * hbar / ( me * (site(j, i)%radius)**2.0 ) ) * hz_to_thz   !largura a meia altura
                site(j, i)%omega0  = site(j, i)%omega * 1.d12 
                site(j, i)%vhomo   = ehomoA
                site(j, i)%vlumo   = elumoA
                site(j, i)%dashtype = 3
            endif

            if ( i > iha .AND. i <= iiaf ) then !iiaf = nanol + ndol + nacl
                site(j, i)%ovlptype = 2
                site(j, i)%mass    = siteMass  
                site(j, i)%radius0 = raioZero
                site(j, i)%radius  = InitSitesRadius(j, i)
                site(j, i)%vel     = InitSitesVel(j, i)  
                site(j, i)%omega   = ( 2.d0 * hbar / ( me * (site(j, i)%radius)**2.0 ) ) * hz_to_thz   !largura a meia altura
                site(j, i)%omega0  = site(j, i)%omega * 1.d12 
                site(j, i)%vhomo   = ehomoA
                site(j, i)%vlumo   = elumoA
                site(j, i)%dashtype = 3
            endif

           if ( i == nc ) then 
                site(j, i)%ovlptype = 2
                site(j, i)%mass    = siteMass  
                site(j, i)%radius0 = raioZero
                site(j, i)%radius  = InitSitesRadius(j, i)
                site(j, i)%vel     = InitSitesVel(j, i)
                site(j, i)%omega   = ( 2.d0 * hbar / ( me * (site(j, i)%radius)**2.0 ) ) * hz_to_thz   !largura a meia altura
                site(j, i)%omega0  = site(j, i)%omega * 1.d12 
                site(j, i)%vhomo   = ehomoCA
                site(j, i)%vlumo   = elumoCA
                site(j, i)%dashtype = 4
           endif

        enddo    
    enddo  
    
    !MODIFICANDO OS SITIOS CASO HAJA DESORDEM
    if (SisDes == .true.) then
        !MODIFICANDO OS DOADORES -> ACEITADORES
       do j = 1, ndes
              rs = desd_r(j)
              cs = desd_c(j)

              site(rs, cs)%ovlptype = 2
              site(rs, cs)%mass    = siteMass  
              site(rs, cs)%radius0 = raioZero
              site(rs, cs)%radius  = InitSitesRadius(rs, cs)
              site(rs, cs)%vel     = InitSitesVel(rs, cs)  
              site(rs, cs)%omega   = ( 2.d0 * hbar / ( me * (site(rs, cs)%radius)**2.0 ) ) * hz_to_thz   !largura a meia altura
              site(rs, cs)%omega0  = site(rs, cs)%omega * 1.d12 
              site(rs, cs)%vhomo   = ehomoA
              site(rs, cs)%vlumo   = elumoA
              site(rs, cs)%dashtype = 3
        enddo   
       
        !MODIFICANDO OS ACEITADORES -> DOADORES
        do j = 1, ndes
              rs = desa_r(j)
              cs = desa_c(j)

              site(rs, cs)%ovlptype = 1
              site(rs, cs)%mass    = siteMass  
              site(rs, cs)%radius0 = raioZero
              site(rs, cs)%radius  = InitSitesRadius(rs, cs)
              site(rs, cs)%vel     = InitSitesVel(rs, cs)  
              site(rs, cs)%omega   = ( 2.d0 * hbar / ( me * (site(rs, cs)%radius)**2.0 ) ) * hz_to_thz   !largura a meia altura
              site(rs, cs)%omega0  = site(rs, cs)%omega * 1.d12 
              site(rs, cs)%vhomo   = ehomoD
              site(rs, cs)%vlumo   = elumoD
              site(rs, cs)%dashtype = 2
        enddo   
    endif
  
  


    do j = 1, nc
        do i = 1, nr
            !PRIMEIRO SITIO LOCALIZADO EM 0,0
            if (i == 1 .AND. j == 1) then
                site(i, j)%xPos = 0.d0
                site(i, j)%yPos = 0.d0
            endif

            !PRIMEIRA COLUNA DE SITIOS
            if (i > 1 .AND. j == 1 ) then
                site(i, j)%xPos = 0.d0
                site(i, j)%yPos = site(i-1 , j)%yPos + 2.d0 * raioZero + distance_molecule 
            endif
              
            !PRIMEIRA LINHA DE SITIOS  
            if (i == 1 .AND. j > 1 ) then
                site(i, j)%xPos = site(i , j-1)%xPos + 2.d0 * raioZero + distance_molecule
                site(i, j)%yPos = 0.d0
            
              !primeiro sitio aceitador
                if ( j == iha ) then
                    site(i, j)%xPos = site(i , j-1)%xPos + 2.d0 * raioZero + distance_int  
                endif
            endif

            if (i > 1 .AND. j > 1) then
                site(i, j)%xPos = site(i , j-1)%xPos + 2.d0 * raioZero + distance_molecule 
                site(i, j)%yPos = site(i-1 , j)%yPos + 2.d0 * raioZero + distance_molecule 

                if ( j == iha ) then
                    site(i, j)%xPos = site(i , j-1)%xPos + 2.d0 * raioZero + distance_int 
                endif
            endif
        enddo
    enddo
    
   ! do j = 1, nc
   !   do i = 1, nr
   !     print*, "tamanho molecula", i, j, "em nm é", 2.d0 * site(i, j)%radius * 1.d9
   !     print*, "omega da molecula", i, j, "em Thz",  site(i, j)%omega
   !     print*, "HBAR * omega da molecula:", i, j, "em ev é",  HB_ev_ps*site(i, j)%omega
   !   enddo
   ! enddo
    
   ! print*, "distancia entre os sitios utilizada:", distance_molecule



    k = 1
    do j = 1, nc
      do i = 1, nr
        site_point(k)%np => site(i, j) !site_pointer
        k = k + 1
      enddo
    enddo

    open( unit=500, file='positionsx.dat', status='replace')
    open( unit=501, file='positionsy.dat', status='replace')
    do k = 1, nsites
        write(500, "(60F20.8)") site_point(k)%np%xPos * 1.d9
        write(501, "(60F20.8)") site_point(k)%np%yPos * 1.d9
    enddo
    close(500)
    close(501) 

    open( unit=601, file='dashtypes.dat', status='replace')
    do k = 1, nsites
        write(601, "(I1)") site_point(k)%np%dashtype
    enddo
    close(601) 


end subroutine define_sitios


subroutine Basis_Builder_DerMtx(disx, disy)
    implicit none

    !args
    real*8, intent(in) :: disx(nsites, nsites), disy(nsites, nsites)
    !local
    real*8 :: DMtx 
    integer             :: i, j, k


    do j = 1, nsites
     do i = 1, nsites
        call DerivativeOverlap(i, j, disx, disy, DMtx)
        basis(i, j)%DerMtx = DMtx
      enddo
    enddo



end subroutine Basis_Builder_DerMtx

subroutine sites_distances(disx, disy)
   implicit none
   !args
   real*8, allocatable, intent(out) :: disx(:,:), disy(:,:)
   !local
   integer :: s1, s2
   integer, allocatable :: upsites(:), downsites(:)
   real*8 :: idst, dcc, threshold_diag, threshold_hor

   allocate(upsites(nc), source = 0 )
   allocate(downsites(nc), source = 0 )
   allocate(disx(nsites, nsites), source = 0.d0 )
   allocate(disy(nsites, nsites), source = 0.d0 )
    allocate(sdist(nsites, nsites),    source = 0.d0 ) 

   idst = float(nr) * ( distance_molecule + ( 2.d0 * raioZero )  )

   do i = 1, nc
       upsites(i) = i * nr
       downsites(i) = 1 + ( (i-1) * nr )
   enddo 


   do s1 = 1, nsites
       do s2 = 1, nsites 
           if (s2 > s1 ) then
              disx(s2, s1) = abs(site_point(s2)%np%xPos - site_point(s1)%np%xPos)
              disy(s2, s1) = abs(site_point(s2)%np%yPos - site_point(s1)%np%yPos)
              if ( Bound_on == .true. ) then
                  do i = 1, nc
                      do j = 1, nc
                          if ( s1 == downsites(i) .and. s2 == upsites(j) ) then
                                disy(s2, s1) = abs(site_point(s2)%np%yPos - idst &
                                               - site_point(s1)%np%yPos)
                          endif 
                          if ( s1 == upsites(i) .and. s2 == downsites(j) ) then
                                disy(s2, s1) = abs(site_point(s2)%np%yPos + idst &
                                               - site_point(s1)%np%yPos)
                          endif 
                      enddo
                  enddo
              endif
           endif
       enddo
   enddo
            
   disy = disy + transpose(disy)
   disx = disx + transpose(disx)
   
   !do i = 1, nsites
   !     do j = 1, nsites
   !       print*, "distancia em x entre o sitio", i, "e", j, "é", disx(i, j)
   !     enddo
   !enddo
   !
   !do i = 1, nsites
   !     do j = 1, nsites
   !       print*, "distancia em y entre o sitio", i, "e", j, "é", disy(i, j)
   !     enddo
   !enddo

   !stop 

    do s1 = 1, nsites-1
        do s2 = s1 + 1, nsites
             sdist(s2, s1) = sqrt( disx(s2, s1)**2.0 + disy(s2, s1)**2.0 )
        enddo
    enddo


    sdist = sdist + transpose(sdist)

    sdist = sdist * 1.d9 


    ! DEFININDO OS PRIMEIROS VIZINHOS

    dcc = 2.d0 * raioZero + distance_molecule

    threshold_diag = sqrt(2.d0 * dcc * dcc) + 0.1d-9 !TOTAL DIAGONAL

    threshold_hor = dcc * 1.d9 !PARA SITIOS HORIZONTAIS

    threshold_diag = threshold_diag * 1.d9
    
    nmfn_t = 0 !total

    nmfn_d = 0 !diagonal

    do s1 = 1, nsites-1
        do s2 = s1 + 1, nsites

            if ( sdist(s2, s1) > threshold_hor .and. sdist(s2, s1) <= threshold_diag ) then
              firstneighbors_diag(s2, s1) = 1.d0
              nmfn_d = nmfn_d + 1
            endif

            if ( sdist(s2, s1) <= threshold_diag ) then
                firstneighbors(s2, s1) = 1.d0
                nmfn_t = nmfn_t + 1
            endif
        enddo
    enddo

    firstneighbors_diag = firstneighbors_diag + transpose(firstneighbors_diag)
    firstneighbors = firstneighbors + transpose(firstneighbors)
    firstneighbors_hor = firstneighbors - firstneighbors_diag

    nmfn_h = nmfn_t - nmfn_d
    !=========================================


deallocate(upsites, downsites)
end subroutine sites_distances


subroutine Overlap(s2, s1, disx, disy, SMtx)
    implicit none

    ! args
    integer,             intent(in)  :: s1, s2
    real*8,              intent(in)  :: disx(nsites, nsites), disy(nsites, nsites)
    real*8,              intent(out) :: SMtx !feita para apenas um estado 
    
    ! local
    real*8  :: x1, y1, w1, x2, y2, w2, x12, y12
    real*8  :: d, d2
    real*8  :: sqmhbar, mhbar
    real*8  :: wproduct, wproduct2
    real*8  :: wden1
    real*8  :: Omega1PlusOmega2
    real*8  :: EXPd2
    
    
    x1 = 0.d0 
    y1 = 0.d0

    x2 = disx(s2, s1)
    y2 = disy(s2, s1)

    !print*, "x", s2, s1, x2
    !print*, "y", s2, s1, y2 

    w1 = site_point(s1)%np%omega*thz_to_hz 
    w2 = site_point(s2)%np%omega*thz_to_hz 
    x12 = x2
    y12 = y2

        
    d   = Sqrt(x12*x12 + y12*y12)
    d2  = d * d
    
    sqmhbar = sqrt( (me/hbar) )
    mhbar   = sqmhbar * sqmhbar
    
    
    wproduct  = sqrt(w1*w2) 
    wproduct2 = wproduct  * wproduct 
    

    Omega1PlusOmega2 = w1 + w2


    wden1 = 1.d0 /  Omega1PlusOmega2
    
    
    EXPd2 = EXP( - HALF * d2 * mhbar * wden1 * wproduct2)
    
    SMtx= 2.d0*EXPd2*wden1*wproduct
    
end subroutine Overlap

subroutine DerivativeOverlap(s2, s1, disx, disy, DerOvlp)        
    implicit none
    !SUBROTINA PARA CALCULAR OS ELEMENTOS <sitio(l1, c1, x, y) | d/dR sitio(l2, c2, x-x0, y-y0) >

    ! args
    integer,             intent(in)  :: s1, s2
    real*8,              intent(in)  :: disx(nsites, nsites), disy(nsites, nsites)
    real*8,              intent(out) :: DerOvlp
    
    ! local
    real*8    :: x1, y1, w1, x2, y2, w2, x12, y12
    real*8    :: w1_2, w2_2 
    real*8    :: d, d2
    real*8    :: sqmhbar, mhbar
    real*8    :: wden1, wden2, wden3, wproduct2
    real*8    :: omega1PlusOmega2
    real*8    :: w1Sw2half, w1Sw2half2
    real*8    :: EXPd2, DerTerm
    real*8    :: wSquareDiff


    x1 = 0.d0
    y1 = 0.d0

    x2 = disx(s2, s1)
    y2 = disy(s2, s1)

    w1 = site_point(s1)%np%omega * thz_to_hz
    w2 = site_point(s2)%np%omega * thz_to_hz

    w1_2 = w1 * w1
    w2_2 = w2 * w2

    x12 = x2
    y12 = y2

    d = sqrt(x12 * x12 + y12 * y12)
    d2 = d*d

    DerTerm = - ((4.d0 * hbar) / (me * (site_point(s2)%np%radius)**3.0)  )

    sqmhbar = sqrt( (me/hbar) )
    mhbar = sqmhbar * sqmhbar

    omega1PlusOmega2  = w1 + w2


    wSquareDiff = w1_2 - w2_2 

    w1Sw2half2 = (w1/w2)
    w1Sw2half  =  sqrt(w1Sw2half2)

    wproduct2 = w1 * w2

    wden1 = 1.d0 / (omega1PlusOmega2) 
    wden2 = wden1 * wden1
    wden3 = wden2 * wden1


    EXPd2 = EXP( - HALF * d2 * mhbar * wden1 * wproduct2)

    DerOvlp = HALF * EXPd2 * w1Sw2half * wden3 * ( wSquareDiff - mhbar * w1_2 * w2 * d2)       
    DerOvlp = DerOvlp * DerTerm  
end subroutine DerivativeOverlap

subroutine print_mat2(aa, nn, mm)
 implicit none
 
     ! args
     integer,   intent(in) :: nn, mm
     real*8,    intent(in) :: aa(nn, mm)
 
     !local
     integer :: i, j
 
     do i=1,mm
       write(*,'(100g12.4)') ( aa(i,j), j=1,nn )
     enddo
 end subroutine print_mat2




  
end module overlap_m
