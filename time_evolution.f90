module time_evolution_m
use f95_precision
use blas95
use lapack95
use types_m
use parameters_m
use constants_m
use functions_m
use system_hamiltonian_m
use rdftensor_m
use rkf_m
use omp_lib 
use verlet_m 

private

public :: PropOfWavePacket, f0, System_Dynamics


contains


!se as matrizes são definidas no parameters.f90, não preciso criar elas como intent(in)
!pois todas as subrotinas utilizam parameters.f90


subroutine System_Dynamics(InitSitesRadius, InitSitesVel)
    implicit none

    ! args 
    real*8,     intent(in)      :: InitSitesRadius(nr, nc), InitSitesVel(nr, nc)

    ! local 
    real*8                      :: ti, tf 
    real*8                      :: energy_el, energy_hl                       !energia do eletron
    real*8                      :: insTemp                          !temperatura instantanea do sistema
    real*8,     allocatable     :: BWRadius(:, :), BWVel(:, :) !Backward variables
    real*8                      :: V_Energy, K_Energy              !energia da mola e energia cinetica
    real*8,     allocatable     :: BWEforce(:,:)                   !aloco na subrotina eletric_force
    real*8,     allocatable     :: BWVforce(:,:)                   !aloco na subrotina spring_force
    real*8,     allocatable     :: hMtx_el(:,:), hMtx_hl(:,:)                       !hamitoniano calculado em t + delta t
    real*8,     allocatable     :: FWRadius(:,:), FWVel(:,:)
    real*8,     allocatable     :: FWEforce(:,:), FWVforce(:,:)    !Forward variables
    real*8,     allocatable     :: rReal_el(:,:), rReal_hl(:,:)                    !parte real da matriz densidade
    complex*16, allocatable     :: rSb_el(:,:), rSb_hl(:,:)                !matriz densidade inicial
    complex*16, allocatable     :: rSa_el(:,:), rSa_hl(:,:)
    real*8                      :: pop_el(nr, nc), pop_hl(nr, nc)          
    integer                     :: image_number 
    real*8                      :: elcoup, hlcoup, endoac_el, endoac_hl, vij_h, vij_d
    real*8,     allocatable     :: disx(:,:), disy(:,:)

    allocate( BWRadius(nr, nc),             source = 0.d0        )
    allocate( BWVel(nr, nc),                source = 0.d0        )
    allocate( rSb_el(d_el, d_el),      source = (0.d0, 0.d0)     )
    allocate( rSb_hl(d_el, d_el),      source = (0.d0, 0.d0)     )
    allocate( rReal_el(d_el, d_el ),         source = 0.d0       ) 
    allocate( rReal_hl(d_el, d_el ),         source = 0.d0       ) 
    allocate( basis(nsites, nsites)                              ) 
 
    !============= DEFINO A CONFIGURAÇÃO INICIAL DOS SÍTIOS
    call define_sitios(InitSitesRadius, InitSitesVel, site, site_point) 
    !======================================================

    !============= MONTA AS MATRIZES AUXILIARES ==========
    call monta_aux

    !============== DEFINO A MATRIZ DE DISTANCIA
    call sites_distances(disx, disy)
    !=======================================

    !============= DEFINO O TERMO DE ENERGIA PARA O CAMPO ELETRICO
    call calculate_efield(disx)
    !========================================
   
 
    !======= SUBROTINA QUE ABRE OS ARQUIVOS DE ESCRITA =============
    call open_write_files
    !==============================================================

    !============ CONDICOES INICIAIS NA BASE CARTESIANA DOS SITIOS para elétron e buraco
    rSb_el(initState, initState) = 1.d0 + 0.d0 * zi ! = matmul(TMtx, matmul(rSPolar_in, TMtxDg) ) 

    rSb_hl(initState, initState) = 1.d0 + 0.d0 * zi ! = matmul(TMtx, matmul(rSPolar_in, TMtxDg) ) 
    !===================================================================================
    

    !============ VELOCIDADES E RAIOS INICIAIS
    BWRadius = site(:, :)%radius                          !raio inicial
    
    BWVel = site(:, :)%vel                                !velocidade inicial
    
    rReal_el = real(rSb_el)                           !parte real da matriz densidade para eletrons
    
    rReal_hl = real(rSb_hl)                           !parte real da matriz densidade para buracos
    !====================================
    
    
    !============ HAMILTONIANA E FORÇAS INICIAIS

    !============= DEFINO A MATRIZ DE RECOMBINACAO ============
    call build_hamiltonian(1, 0.d0, rReal_el, rReal_hl, disx, disy, hMtx_el, elcoup)
    
    call build_hamiltonian(2, 0.d0, rReal_el, rReal_hl, disx, disy, hMtx_hl, hlcoup)
    
    call define_rec_ham(rSb_el, rSb_hl, hMtx_el, hMtx_hl)

    call eletric_force(0, rReal_el, rReal_hl, 0.d0, disx, disy, BWEforce) !Backward Eforce
    
    call spring_force(BWVforce)                           !Backward SpringForce
    !================================================================
        
    


    !============ TEMPO INICIAL
    ti = 0.d0


    !=========== CÁLCULO DAS ENERGIAS DAS PARTICULAS NO INSTANTE INICIAL

    call particle_energy(hMtx_el, rSb_el, energy_el, endoac_el)


    call particle_energy(hMtx_hl, rSb_hl, energy_hl, endoac_hl)

    
    !==================================================================
    
    
    !=========== CÁLCULO DA ENERGIA CINETICA, TEMPERATURA E POTENCIAL NO INSTANTE INICIAL ========
    call kinect_energy(BWVel, K_Energy)

    call calculate_temperature(K_Energy, insTemp)

    if (FixTemp == .true.) then
        insTemp = bathTemp 
    endif

    call spring_energy(BWRadius, V_Energy) 
    !===================================================

    call calculate_mean_vij(hMtx_el, vij_h, vij_d) 

    !=========== ESCREVO OS RESULTADOS NO INSTANTE INICIAL

    !energia cinética, potencial, temperatura e raios
    call write_files(ti, energy_el, energy_hl, elcoup, hlcoup, K_Energy, V_Energy, insTemp, endoac_el, endoac_hl, vij_h, vij_d)


    !populações
    call printa_resultado(15, 1, ti, rSb_el)

    call printa_resultado(16, 2, ti, rSb_hl)
    !===================================================





    !calcula dt 
    dt = tmax/float(nm_divisoes) 




    !================== COMEÇO DA DINÂMICA ====================================
    do pl = 1, nm_divisoes 
        !if (pl == 1 .OR. pl == nm_divisoes/2 .OR. pl == nm_divisoes ) then
        !    print*, "==== PASSO",pl,"===="
        !endif


        writepoint = mod( float(pl), float(wrcons) ) 


        !============= DEFINE O TEMPO FINAL E O DELTA_t ===========
        tf = float(pl)*(tmax)/float(nm_divisoes)
        !==========================================================



        !============== EVOLUCAO PARA O ELETRON DE ti ATÉ tf =======
        call PropOfWavePacket(1, pl, ti, tf, insTemp, hMtx_el, rSb_el, rSa_el) !propago função de onda
        
        call PropOfWavePacket(2, pl, ti, tf, insTemp, hMtx_hl, rSb_hl, rSa_hl) !propago função de onda
        !===========================================================

      
        !============= ATUALIZO AS CONDICOES INICIAIS QUANTICAS ==============
        rSb_el = rSa_el !atualizo cond. inicial para o eletron
        rReal_el = real(rSb_el) 

        rSb_hl = rSa_hl !atualizo cond. inicial para o eletron
        rReal_hl = real(rSb_hl) 
        !===========================================================
    

        !============ CALCULO A EVOLUCAO CLASSICA COM O ALGORITMO VELOCITY VERLET ===========================
    
        call velocity_verlet(pl, bathTemp, BWRadius, BWVel, BWEforce, BWVforce, ti, rReal_el, rReal_hl, disx, disy, &
                             FWRadius, FWVel, FWEforce, FWVforce, insTemp, V_Energy, K_Energy) 
                  !calculamos a energia cinetica e da mola no algoritmo de verlet em um tempo ti para ser igual ao resto do programa
        !====================================================================================================
        
        


        !=========== ATUALIZO A PARTE CLASSICA =============================================================
        BWRadius(:, :) = FWRadius(:, :) !aloco em velocity_verlet
        BWVel(:, :)    = FWVel(:, :)    !aloco em velocity_verlet
        BWEforce(:, :) = FWEforce(:, :) !aloco em eletric_force
        BWVforce(:, :) = FWVForce(:, :) !aloco em spring_force
        deallocate(FWEforce, FWVel, FWRadius, FWVforce) 
        !===================================================================================================
 

        !=============== CALCULO DO HAMILTONIANO - 1 el - 2 hl  ================
        call build_hamiltonian(1, tf, rReal_el, rReal_hl, disx, disy, hMtx_el, elcoup)


        call build_hamiltonian(2, tf, rReal_el, rReal_hl, disx, disy, hMtx_hl, hlcoup)
        !==========================================================

        !============= DEFINO A MATRIZ DE RECOMBINACAO ============
        call define_rec_ham(rSb_el, rSb_hl, hMtx_el, hMtx_hl) 
      


        if ( writepoint == 0.d0 ) then

            call particle_energy(hMtx_el, rSb_el, energy_el, endoac_el)
 
            call particle_energy(hMtx_hl, rSb_hl, energy_hl, endoac_hl)

            call rhomtx_to_pop(rSb_el, pop_el)

            call rhomtx_to_pop(rSb_hl, pop_hl)

            call calculate_mean_vij(hMtx_el, vij_h, vij_d) 
     

            !====== ESCREVE AS POPULACOES =====================
            call printa_resultado(15, 1, tf, rSb_el)
            
            call printa_resultado(16, 2, tf, rSb_hl)
            !==================================================
        
            !============ ESCREVE AS ENERGIAS, RAIOS E FORÇAS DE SAÍDA ============================================
            call write_files(tf, energy_el, energy_hl, elcoup, hlcoup, K_Energy, V_Energy, insTemp, endoac_el, endoac_hl, vij_h, vij_d)
            !====================================================================




        endif 
    

        !========= ATUALIZO O TEMPO ===========
        ti = tf
        !=====================================

    enddo

        if ( Therm_on == .true. ) then 
            print*, "Termalização acabou, escrevendo os arquivos de saída"
          
            open( unit = 301, file = "therm_radius.int", status = "replace" )
            write(301, "(60F20.8)") ( ( site_point(i)%np%radius*1.d9 ),  i = 1, nsites ) 
            close(301) 


           open( unit = 302, file = "therm_vel.int", status = "replace" )
           write(302, "(60F20.8)") ( ( site_point(i)%np%vel ),  i = 1, nsites ) 
           close(302) 
        endif     



    !========= FECHO OS ARQUIVOS DE ESCRITA =====
    call close_write_files
    !============================================
    
    
    
      DEALLOCATE(rSb_el, rSb_hl, rReal_el, rReal_hl, basis, BWEforce, BWVforce, hMtx_el, hMtx_hl, BWRadius, BWVel)
end subroutine System_Dynamics



subroutine PropOfWavePacket(eoh, step, ti, tf, insTemp, hMtx, rSb, rSa)
    implicit none
   
    ! args 
    integer,                 intent(in)  :: eoh
    integer,                 intent(in)  :: step 
    real*8,                  intent(in)  :: ti 
    real*8,                  intent(in)  :: tf 
    real*8,                  intent(in)  :: insTemp
    real*8,                  intent(in)  :: hMtx(d_el, d_el) 
    complex*16,              intent(in)  :: rSb(d_el, d_el)
    complex*16, allocatable, intent(out) :: rSa(:,:)
    
    ! local 
    real*8,     allocatable  :: energias(:), y(:), yp(:), work(:)
    REAL*8,     allocatable  :: phi(:,:), phi_transpose(:,:)
    real*8,     allocatable  :: frequency_matrix(:, :)
    complex*16, allocatable  :: recHb(:, :)
    COMPLEX*16, allocatable  :: rHb(:,:), rHa(:,:)
    INTEGER                  :: number_file_ham
    real*8                   :: abserr
    !integer                  :: iflag
    integer                  :: flag
    integer                  :: iwork(5)
    real*8                   :: relerr
    real*8                   :: t !esse t é a variável independente do ODE
    real*8                   :: tout 
    real*8                   :: energ(d_el) 
  
    
    
    
     allocate( rSa(d_el, d_el)           , source = (0.d0, 0.d0) )
     allocate( y(neqn_el)                     , source = 0.d0         )
     allocate( yp(neqn_el)                    , source = 0.d0         )
     allocate( work(100+21*neqn_el))
     allocate( rHb(d_el, d_el)          , source = (0.d0, 0.d0) ) !populaca inicial na base da hMtx
     allocate( rHa(d_el, d_el)             , source = (0.d0, 0.d0) )
     allocate( recHb(d_el, d_el)        , source = (0.d0, 0.d0)       )
    
     abserr = 1.d-8
     relerr = 1.d-8
     !iflag = 1 
     flag = 1
     t = ti !construo o tempo inicial de uma forma que o ODE nao atualiza!
     tout = tf !TOUT É O TEMPO DE SAÍDA, OU SEJA, DELTA T = tout(2) - tout(1)
     number_file_ham = 14
    
 
    !==== CALCULA OS AUTOESTADOS E AUTOVETORES =======
    call calculate_eigenvectors(step, hMtx, energias, phi, phi_transpose, frequency_matrix) 
    !===============================================================================================
    
    !==== PRINTA A MATRIZ DENSIDADE, AUTOVETORES E HAMILTONIANO =====
    call print_matrices(eoh, step, rSb, phi, hMtx) 
    !===============================================================================================

    !==== CONSTRUO O HAMILTONIANO DE RECOMBINACAO ===
    if ( eoh == 1 ) call recsite_to_recham(phi, phi_transpose, recSb, recHb) 
    if ( eoh == 2 ) call recsite_to_recham(phi, phi_transpose, recSb_hl, recHb) 
    !================================================

    !=========== CONDICOES INICIAIS NA BASE EXC ======================
    call rhosite_TO_rhoham(phi, phi_transpose, rSb, rHb)
    !=================================================================
      
    !======== TRANSFORMO AS CONDICOES INICIAIS PARA OS VETORES LINHAS DO ODE ========
    call monta_y(rHb, y) !condição inicial para resolver as ODE'S
    !================================================================================

    !====== CRIO A MATRIZ RW QUE É UTILIZADA PARA CALCULAR O ODE =================
    call createRwMatrix(eoh, rSb, hMtx, insTemp, frequency_matrix, phi, recHb, rwMtx)
    !=============================================================================
    
    
    !y esta na forma (para d=2):  y(1) = r(1, 1)  y(3) = real(r(2, 1))
    !                             y(2) = r(2, 2)  y(4) = aimag(r(2, 1))
    !---------------------------------------------------------------------------------- 
    
    !======== SUBROTINA QUE CALCULA AS EQUACOES DIFERENCIAIS DE T ATÉ TOUT ========
    !call ode ( f0, neqn_el, y, t, tout, relerr, abserr, iflag, work, iwork )
    call r8_rkf45 ( f0, neqn_el, y, yp, t, tout, relerr, abserr, flag )

    if ( flag /= 2 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PropOfWavePacket - Fatal error!'
      write ( *, '(a,i8)' ) '  RKF returned FLAG = ', flag
      stop
    end if
    !==============================================================================
    
    
    
    !========== TRANSFORMO O RESULTADO DO VETOR LINHA Y DO ODE PARA AS MATRIZES r ===========
    call monta_rho(y, rHa)
    !==========================================================================================
    
    
    !============== ESCREVO O OPERADOR DENSIDADE NA BASE DO SITIO =============================
    call rhoham_to_rhosite(phi, phi_transpose, rHa, rSa)
    !==========================================================================================
    
     
        
    
    DEALLOCATE(phi, phi_transpose, rHb, rHa, recHb, rwMtx)
    13 format(3es14.3E3)
end subroutine PropOfWavePacket



subroutine f0 ( t, y, yp )
!não definimos as condicoes iniciais aqui, apenas a funcao que vamos calcular a derivada. A condição inicial vem antes do call
    implicit none

    real*8 ::  t !tempo que vai entrar na equação diferencial
    real*8 ::  y(neqn_el)  !funcao que estamos querendol, neste caso x
    real*8 ::  yp(neqn_el) !derivada que estamos querendol calcular, por exemplo dx/dt = -x => sol: x = exp(-t)
    
    call gemv(rwMtx, y, yp)  !gemv calcula o produto da matriz RWMatrix com
                                !o vetor y resultandol em yp => muito mais rapido
    
    13 format(3es14.3E3)
end subroutine f0



end module time_evolution_m
