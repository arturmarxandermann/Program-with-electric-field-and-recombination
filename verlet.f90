module verlet_m
use f95_precision
use lapack95
use blas95
use parameters_m
use constants_m 
use types_m 
use system_hamiltonian_m
use functions_m 

contains 

subroutine particle_energy(h_mtx, rho, ParticleEnergy, endoac)
    implicit none

    ! args
    real*8,     intent(in) :: h_mtx(d_el, d_el)
    complex*16, intent(in) :: rho(d_el, d_el)
    real*8,     intent(out):: ParticleEnergy
    real*8,     intent(out):: endoac
   
    ! local  
    integer :: i 
    complex*16, allocatable :: MtxEnergy(:,:)

    allocate(MtxEnergy(d_el, d_el), source = (0.d0, 0.d0) ) 
    
    call gemm(h_mtx, rho, MtxEnergy)
    
    ParticleEnergy = 0.d0
    endoac = 0.d0 

    do i = 1, d_el
       ParticleEnergy = ParticleEnergy + real(MtxEnergy(i, i))
    enddo

    do i = nano + 1, nsites-nca
       endoac = endoac + real(MtxEnergy(i, i))
    enddo
    
    deallocate(MtxEnergy) 
end subroutine particle_energy        

subroutine calculate_mean_vij(hMtx, vij_h, vij_d)
    implicit none
    !args
    real*8, intent(in) :: hMtx(d_el, d_el)
    real*8, intent(out) :: vij_h, vij_d

    vij_h = 0.d0
    vij_d = 0.d0

    do j = 1, d_el-1
      do i = j + 1, d_el
          vij_h = vij_h + hMtx(i, j) * firstneighbors_hor(i, j)
          vij_d = vij_d + hMtx(i, j) * firstneighbors_diag(i, j)
      enddo
    enddo
 
    vij_h = vij_h/float(nmfn_h)

    if (nmfn_d /= 0) then
      vij_d = vij_d/float(nmfn_d)
    else
      vij_d = 0.d0
    endif


end subroutine calculate_mean_vij
    
    
subroutine kinect_energy(VelSite, K_Energy)
    implicit none

    ! args
    real*8, intent(in)  :: VelSite(nr, nc)
    real*8, intent(out) :: K_Energy

    ! local
    real*8,  allocatable :: K_EnergySites(:,:)
    
    allocate(K_EnergySites(nr, nc), source = 0.d0)

    K_EnergySites = HALF * site%mass * VelSite * VelSite
    
    K_Energy = sum(K_EnergySites(:,:)) * joule_to_ev

    !write( 500, "(60F20.6)", advance = 'no') ti  
    !write( 500, "(60F20.6)") ( (K_EnergySites(1, i) * joule_to_ev ), i = 1, nsites )
    
    deallocate(K_EnergySites)
end subroutine kinect_energy

 

subroutine spring_energy(RadiusSite, V_Energy)
    implicit none

    ! args
    real*8, intent(in)  :: RadiusSite(nr, nc)
    real*8, intent(out) :: V_Energy
    
    ! local  
    real*8, allocatable :: V_EnergySites(:,:)
    real*8, allocatable :: omega0Squared(:,:)
    real*8, allocatable :: raDiff(:,:)
    real*8, allocatable :: raDiffSquared(:,:)

    
    
    allocate( V_EnergySites(nr, nc),         source = 0.d0 ) 
    allocate( omega0Squared(nr, nc),         source = 0.d0 ) 
    allocate( raDiff(nr, nc),                source = 0.d0 ) 
    allocate( raDiffSquared(nr, nc),         source = 0.d0 ) 
    
    
    
    omega0Squared = freqMode * freqMode
    raDiff        = RadiusSite - site%radius0
    raDiffSquared      = raDiff * raDiff


    V_EnergySites = HALF * site%mass * omega0Squared * raDiffSquared
    
    V_Energy = sum(V_EnergySites(:, :)) * joule_to_ev

    !write( 501, "(60F20.6)", advance = 'no') ti  
    !write( 501, "(60F20.6)") ( (V_EnergySites(1, i) * joule_to_ev ), i = 1, nsites )
    
    deallocate(V_EnergySites, omega0Squared, raDiff, raDiffSquared)
end subroutine spring_energy        



subroutine spring_force(Vforce)
    implicit none

    ! args
    real*8, allocatable, intent(out) :: Vforce(:,:)
    
    ! local 
    integer :: i, j 
    real*8, allocatable :: omega0Squared(:,:), raDiff(:,:)
    
    
    ALLOCATE( Vforce(nr, nc),        source = 0.d0) !Potential Force
    ALLOCATE( omega0Squared(nr, nc), source = 0.d0)
    ALLOCATE( raDiff(nr, nc),        source = 0.d0) !Radius Difference

    omega0Squared = freqMode * freqMode
    raDiff        = site%radius0 - site%radius

    Vforce = site%mass * omega0Squared * raDiff

    !Vforce = 0.d0 

    deallocate(omega0Squared, raDiff)
end subroutine spring_force

subroutine eletric_force(pl, rhoElReal, rhoHlReal, ti, disx, disy, Eforce)
    implicit none

    ! args
    integer,             intent(in)  :: pl
    real*8,              intent(in)  :: rhoElReal(d_el, d_el), rhoHlReal(d_el, d_el)
    real*8,              intent(in)  :: ti
    real*8,              intent(in)  :: disx(nsites, nsites), disy(nsites, nsites)
    real*8, allocatable, intent(out) :: Eforce(:,:)

    ! local
    real*8,  allocatable :: MtxEforceDiagonal(:,:)
    real*8               :: ForceStates
    real*8,  allocatable :: derTerm(:,:)
    real*8,  allocatable :: DerMtx(:,:), RhoTimesDer_el(:,:), RhoTimesDer_hl(:,:)
    real*8,  allocatable :: MtxEforceND1(:,:), MtxEforceND1_el(:,:), MtxEforceND1_hl(:,:)
    real*8               :: writepoint


    allocate( MtxEforceND1(nr, nc)      , source = 0.d0 )
    allocate( MtxEforceND1_el(nr, nc)      , source = 0.d0 )
    allocate( MtxEforceND1_hl(nr, nc)      , source = 0.d0 )
    allocate( RhoTimesDer_el(d_el, d_el)       , source = 0.d0 )
    allocate( RhoTimesDer_hl(d_el, d_el)       , source = 0.d0 )
    allocate( Eforce(nr, nc)            , source = 0.d0 )
    allocate( MtxEforceDiagonal(nr, nc) , source = 0.d0 )
    allocate( derTerm(nr, nc)           , source = 0.d0 )

    !========== CALCULO DAS DERIVADAS DW/DR ===========================
    derTerm = - (( 4.d0 * hbar ) / ( me * (site%radius)**3.0 ))
    !==================================================================

    !=========== CALCULO DA MATRIZ DERIVA =========
    call build_derivative_matrix(disx, disy, DerMtx)
    !==============================================


    !if (pl == 1 .OR. pl == nm_divisoes/2 .OR. pl == nm_divisoes ) then
    !  print*, "Matriz derivada"
    !  call print_mat2(DerMtx, d_el, d_el)
    !endif

    !======= PARTE DIAGONAL DA FORCA - DEPENDENTE DAS POPULACOES =============
    ForceStates = 0.d0 
    l = 1
    do j = 1, nc
       do i = 1, nr

          ForceStates = (rhoHlReal(l, l) + rhoElReal(l, l)) * hbar 


          MtxEforceDiagonal(i, j) = derTerm(i, j) * ForceStates
          
       l = l + 1
       enddo
    enddo
    !=========================================================================


   !===== PRIMEIRA PARTE NAO DIAGONAL - DEPENDENTE APENAS DAS COERÊNCIAS =====

   RhoTimesDer_el = matmul(rhoElReal, DerMtx)
   RhoTimesDer_hl = matmul(rhoHlReal, DerMtx)


   l = 1
   do j = 1, nc
     do i = 1, nr

       MtxEforceND1_el(i, j) = RhoTimesDer_el(l, l)
       MtxEforceND1_hl(i, j) = RhoTimesDer_hl(l, l)

        l = l + 1

       enddo
    enddo

   !==========================================================================
   MtxEforceND1 =  2.d0 * ev_to_joule * scD * MtxEforceND1_el + &
                    2.d0 * ev_to_joule * scD * MtxEforceND1_hl


   

   ! ================== FORCA ELETRICA TOTAL =====================
    if (WithForces == .true. ) then
        if ( Therm_on == .false. ) then
            Eforce = - forcecons * (MtxEforceDiagonal + MtxEforceND1)
        else
            Eforce = 0.d0
        endif
    else
            Eforce = 0.d0
    endif
   !==============================================================
    
    writepoint = mod( float(pl), float(wrcons) )
    if (writepoint == 0.d0 ) then
        write(450, "(60F20.8)") ti, forcecons * 1.d10 * (sum(Eforce(:,:))/float(nsites)) 
        write(451, "(60F20.8)") ti, - forcecons * 1.d10 * (sum(MtxEforceDiagonal(:,:))/float(nsites)) 
        write(452, "(60F20.8)") ti, - forcecons * 1.d10 * (sum(MtxEforceND1(:,:))/float(nsites)) 
    endif
   !forces.dat
   !if ( writepoint == 0.d0 ) then
   !     write( 150, "(60F20.8)", advance = "no" ) ti
   !     write( 150, "(60F20.8)", advance = "no" ) (( - MtxEforceDiagonal(1, i) * 1.d10 ), i = 1, nsites )
   !     write( 150, "(60F20.8)" ) (( - MtxEforceND1(1, i) * 1.d10 ), i = 1, nsites )
   !endif

   deallocate(derTerm, MtxEforceDiagonal, RhoTimesDer_el, RhoTimesDer_hl, MtxEforceND1, MtxEforceND1_el, MtxEforceND1_hl)
end subroutine eletric_force


subroutine calculate_temperature(K_erg, temperature)
    implicit none
    !args
    real*8,                intent(in)  :: K_erg   !energia cinetica
    real*8,                intent(out) :: temperature !temperatura instantanea


    temperature = ( 2.d0 * K_erg) / (float(nsites) * kboltz)   !temp instantanea variavel global


end subroutine calculate_temperature



subroutine calculate_BRDscale(temperature, bathTemp, BRDscale) 
    implicit none

    !args
    real*8,                intent(in)  :: temperature
    real*8,                intent(in)  :: bathTemp
    real*8,                intent(out) :: BRDscale

    !local
    real*8                             :: timeScale, tempScale, BRDsq, temp


    temp = temperature

    if (temp == 0.d0 ) temp = 1.d0 

    tempScale = bathTemp/temp 
    timeScale = (dt*1.d-12)/bathCoup  


    BRDsq = ( (tempScale - 1.d0) * timeScale ) + 1.d0 

    BRDscale = sqrt(BRDsq)  

end subroutine calculate_BRDscale


subroutine velocity_verlet(pl, bathTemp, BWRadius, BWVel, BWEforce, BWVforce, ti, rhoElReal, rhoHlReal, &
                          disx, disy, FWRadius, FWVel, FWEforce, FWVforce, insTemp, V_Energy, K_Energy )
    implicit none

    ! args
    integer,               intent(in)  :: pl
    real*8,                intent(in)  :: bathTemp  
    real*8,                intent(in), dimension(nr, nc) :: BWRadius, BWVel, BWEforce, BWVforce
    real*8,                intent(in)  :: ti
    real*8,                intent(in)  :: rhoElReal(d_el, d_el), rhoHlReal(d_el, d_el)
    real*8,                intent(in)  :: disx(nsites, nsites), disy(nsites, nsites)
    real*8,  allocatable,  intent(out) :: FWRadius(:,:), FWVel(:,:), FWEforce(:,:), FWVforce(:,:)
    real*8,                intent(out) :: insTemp, V_Energy, K_Energy
    
    ! local  
    real*8, allocatable :: BWForce(:,:), FWForce(:,:)
    real*8, allocatable :: BWAcc(:,:), FWAcc(:,:)
    integer             :: k, i, j 
    real*8              :: dt_ps
    real*8              :: BRDscale
    real*8, allocatable :: FWVelTherm(:,:)
    !BWRadius -> BackWardRadius (raio anterior).    FWRadius -> ForWardRadius (raio posterior).
    
    allocate( FWRadius(nr, nc)    , source = 0.d0 ) 
    allocate( FWVel(nr, nc)       , source = 0.d0 )
    allocate( BWForce(nr, nc)     , source = 0.d0 )
    allocate( FWForce(nr, nc)     , source = 0.d0 )
    allocate( BWAcc(nr, nc)       , source = 0.d0 )
    allocate( FWAcc(nr, nc)       , source = 0.d0 ) 
    allocate( FWVelTherm(nr, nc)  , source = 0.d0 )

    !ALOCO FWEforce na subrotina eletric_force
    dt_ps = dt * 1.d-12   

    !====== FORCA TOTAL EM t   
    BWForce =  BWEforce + BWVforce
    
    !======= ACELERACAO COM A FORCA EM t
    BWAcc = BWForce / site%mass
    
    !====== CALCULO DA ENERGIA CINETICA E POTENCIAL COM A FORCA EM t 
    call kinect_energy(BWVel, K_Energy) 
    call spring_energy(BWRadius, V_Energy)

    
    !===== CALCULO A TEMPERATURA DO SISTEMA COM A NOVA ENERGIA CINETICA
    call calculate_temperature(K_energy, insTemp) 

    if ( FixTemp == .true. ) then
        insTemp = bathTemp 
    endif

    !===== BERENDSEN SCALE 
    call calculate_BRDscale(insTemp, bathTemp, BRDscale)

    if (BRD_on == .false. ) BRDscale = 1.d0
    
    !======= CALCULO DOS NOVOS RAIOS EM t
    FWRadius = BWRadius + (BWVel * dt_ps) + (HALF * BWAcc * dt_ps * dt_ps)

    !======= ATUALIZO O RAIO E OMEGA 
    site%radius = FWRadius
    site%omega = ( ( 2.d0 * hbar) / ( me  * site%radius * site%radius ) ) * hz_to_thz

    !======= CALCULO DAS NOVAS FORCAS EM t + delta t
    call eletric_force(pl, rhoElReal, rhoHlReal, ti, disx, disy, FWEforce)
    call spring_force(FWVforce) 
    
    !======= TOTAL FORWARD FORCE 
    FWForce = FWEforce + FWVforce

    !======== NOVA ACELERACAO
    FWAcc = FWForce / site%mass
    
    !======== NOVA VELOCIDADE RADIAL
    FWVel = BWVel + HALF * dt_ps * (BWAcc + FWAcc) 
    
    !======== BERENDSEN THERMOSTAT FOR DONNOR AND ACCEPTOR
    FWVel = FWVel * BRDscale 
    
    !======== ATUALIZO AS VELOCIDADES
    site%vel = FWVel
    
 

    deallocate(BWForce, FWForce, BWAcc, FWAcc, FWVelTherm) 
end subroutine velocity_verlet

end module verlet_m
