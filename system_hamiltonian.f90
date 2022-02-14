module system_hamiltonian_m
use f95_precision
use lapack95
use blas95
use types_m
use parameters_m
use constants_m
use overlap_m
use functions_m

contains

subroutine print_mat3(aa, nn, mm)
    implicit none

    ! args
    integer, intent(in) :: nn, mm
    real*8 , intent(in) :: aa(nn, mm)

    !local 
    integer :: i, j

    do, i=1,mm
        write(*,'(100g12.4)') ( aa(i,j), j=1,nn )
    enddo
end subroutine print_mat3


subroutine build_hamiltonian(eoh, time, rReal_el, rReal_hl, disx, disy, hMtx, elhlen)
    implicit none

    !args
    integer, intent(in) :: eoh
    real*8, intent(in) :: time 
    real*8, intent(in) :: rReal_el(nsites, nsites), rReal_hl(nsites, nsites)
    real*8, intent(in) :: disx(nsites, nsites), disy(nsites, nsites)
    real*8, intent(out), allocatable :: hMtx(:,:)
    real*8, intent(out) :: elhlen 

    !local
    real*8              :: SMtx
    integer             :: i, j, k
    real*8, allocatable :: eeh_vec(:, :), energyeh_a(:) 

    integer             :: c1, r1, s1, s2, c2, r2

    
    integer :: row_min, row_max, col_min, col_max
    
    
    allocate(hMtx(d_el , d_el)        , source = 0.d0 )
    allocate(eeh_vec(nsites, nsites)            , source = 0.d0 ) 
    allocate(energyeh_a(nsites)           , source = 0.d0 ) 

    
    
    ! === CONSTRUO OS BLOCOS DA MATRIZ ====
    ! === TERMO DE LIGACAO EL-HL
    if ( ElHlInt == .true. ) then
      elhlen = 0.d0
      do j = 1, nsites
        do i = 1, nsites 


                if ( eoh == 1 ) then
                    eeh_vec(i, j) = (rReal_hl(i, i) * ehcoup) / (1.d0 + sdist(i, j)) !coupling val
                    elhlen = elhlen + eeh_vec(i, j) * rReal_el(j, j)                 !energy val
                endif
                

                if ( eoh == 2 ) then
                    eeh_vec(i, j) = (rReal_el(i, i)  * ehcoup) / (1.d0 + sdist(i, j)) !coupling val
                    elhlen = elhlen + eeh_vec(i, j) * rReal_hl(j, j)                  !energy val
                endif 


        enddo
        energyeh_a(j) = sum(eeh_vec(:, j))
      enddo
    endif
    
    if (ElHlInt == .false. ) then
        energyeh_a = 0.d0
        elhlen = 0.d0
    endif
    !====================================================================



    ! OFF DIAGONAL PART OF HAMILTONIAN
     do j = 1, nsites-1
        do i = j + 1, nsites
            call Overlap(i, j, disx, disy, SMtx)
            hMtx(i,j) = SMtx * scD 
        enddo
     enddo

    hMtx = hMtx + transpose(hMtx) 

    !================================



    ! DIAGONAL PART OF HAMILTONIAN
    if (eoh == 1) then 
      do i = 1, nsites
          hMtx(i, i) = site_point(i)%np%vlumo + HB_ev_ps * site_point(i)%np%omega - energyeh_a(i) &
                       - eterm(i)
      enddo

    endif

    if (eoh == 2) then
      do i = 1, nsites
          hMtx(i, i) = site_point(i)%np%vhomo + HB_ev_ps * site_point(i)%np%omega - energyeh_a(i) &
                       + eterm(i)
      enddo
    endif
    ! ====================================

    ! ONLY FIRST NEIGHBORS
    do i = 1, nsites
        do j = 1, nsites
            if ( abs(hMtx(i, j)) <= 1.d-4) then
                hMtx(i, j) = 0.d0
            endif
        enddo
    enddo
    !-------------------

deallocate( eeh_vec, energyeh_a )
end subroutine build_hamiltonian


subroutine calculate_efield(disx)
  implicit none
  !args
  real*8, intent(in) :: disx(nsites, nsites)

  allocate( eterm(nsites),  source = 0.d0 ) 


    if (EField == .true. ) then
      do i = 1, nsites
          eterm(i) = disx(1, i) * el_field * 1.d9
      enddo
    else
      eterm = 0.d0
    endif
  


end subroutine calculate_efield

subroutine define_rec_ham(rden_el, rden_hl, hel, hhl)
    implicit none
    !args
    complex*16, intent(in) :: rden_el(d_el, d_el)
    complex*16, intent(in) :: rden_hl(d_el, d_el)
    real*8,     intent(in) :: hel(d_el, d_el), hhl(d_el, d_el) 
    !local
    integer :: j, i, k
    real*8                 :: trel, trhl

    recSb = 0.d0
    recSb_hl = 0.d0
    trel = 0.d0
    trhl = 0.d0
    k = 1

    ! ================== MODELO NÃO DIAGONAL SEM OPERADOR ===================
    do j = 1, nsites-1
        do i = j + 1, nsites
            recSb(i, j) = - rec_rate * rden_hl(i, j) 
            recSb_hl(i, j) = - rec_rate * rden_el(i, j) 
        enddo
    enddo

    recSb = recSb + transpose(conjg(recSb))
    recSb_hl = recSb_hl + transpose(conjg(recSb_hl))

    do i = 1, nsites
        recSb(i, i) = - rec_rate * rden_hl(i, i)
        recSb_hl(i, i) = - rec_rate * rden_el(i, i)
    enddo

    ! =========================================================

    ! ================== MODELO NÃO DIAGONAL COM OPERADOR ===================
    !do j = 1, nsites-1
    !    do i = j + 1, nsites
    !        recSb(i, j) = - rec_rate * (rden_hl(i, j) * hel(i, j))/scD
    !        recSb_hl(i, j) = - rec_rate * (rden_el(i, j) * hhl(i, j))/scD
    !    enddo
    !enddo

    !recSb = recSb + transpose(conjg(recSb))
    !recSb_hl = recSb_hl + transpose(conjg(recSb_hl))

    !do i = 1, nsites
    !    recSb(i, i) = - rec_rate * rden_hl(i, i)
    !    recSb_hl(i, i) = - rec_rate * rden_el(i, i)
    !enddo

    ! =========================================================

    ! ================= MODELO DIAGONAL =======================
    !do j = 1, nc
    !    do i = 1, nr
    !        !if ( j == (nanol + ndol) .OR. j == (nanol + ndol + 1) )  recSb(k, k) = - 5.d0 * rReal_hl(k, k)
    !        !if ( j == (nanol + ndol) .OR. j == (nanol + ndol + 1) )  recSb_hl(k, k) = - 5.d0 * rReal_el(k, k)
    !        !recSb(k, k) = - 5.d0 * rReal_hl(k, k)
    !        !recSb_hl(k, k) = - 5.d0 * rReal_el(k, k)
    !        k = k + 1
    !    enddo
    !enddo
    ! =========================================================

    ! ================ MODELO GRAZI ==========================

    !do i = 1, nsites
    !    trel = trel + real(rden_el(i, i))
    !    trhl = trhl + real(rden_hl(i, i))
    !enddo

    !recSb = -rec_rate * trel * rden_hl
    !recSb_hl = -rec_rate * trhl * rden_el
 
    ! =======================================================

end subroutine define_rec_ham


subroutine calculate_eigenvectors(pl, hamiltoniana, energias, phi, phi_transpose, omega_matrix)
!calcula os autovetores (phi) e autovalores (energias) da hamiltoniana. omega_matriz é a matriz de frequencias angulares
    implicit none

    ! args
    integer             , intent(in)  :: pl
    real*8              , intent(in)  :: hamiltoniana(d_el,d_el)
    real*8 , allocatable, intent(out) :: energias(:)
    real*8 , allocatable, intent(out) :: phi(:,:), phi_transpose(:,:), omega_matrix(:,:)

    ! local 
    real*8, allocatable ::  H_temp(:,:)
    character(len=1)    :: jobz, uplo
    integer             :: info, is, js 

    allocate(phi(d_el, d_el)               , source = 0.d0)
    allocate(phi_transpose(d_el, d_el)     , source = 0.d0)
    allocate(omega_matrix(d_el, d_el)      , source = 0.d0)
    allocate(energias(d_el)                , source = 0.d0)
    allocate(H_temp(d_el, d_el)            , source = 0.d0)
   
    info = 0 
    H_temp = hamiltoniana
    
    call syevd(H_temp, energias, "V", "L", info)
    
    !If info = 0, contains the n orthonormal eigenvectors, stored by columns.
    !(The i-th column corresponds to the ith eigenvalue.)
    !Como o iézimo autovetor é disposto na iézima coluna, C_n^N = <n|N> é o coeficiente da expansão do autovetor
    
    if ( info /= 0 ) write(*,*) info, "Error to found eigenvalues"
    
    !if ( pl == 1 .OR. pl == nm_divisoes/2 .OR. pl == nm_divisoes ) then
      !do i = 1, d_el
      !print*, "ENERGIA", i, energias(i)
      !enddo
    !endif


    !========== CONSTRUCAO DA MATRIZ DE AUTOVETORES ========
    do n = 1, d_el
      phi(n,:) = H_temp(n,:)
    enddo
    !=======================================================
    
    
    phi_transpose = transpose(phi)
    
    !============== CONSTRUÇÃO DA MATRIZ DE FREQUÊNCIAS ========================
    do j = 1, d_el-1
      do i = j + 1, d_el
        omega_matrix(i, j) = (energias(i) - energias(j)) / HB_ev_ps 
      enddo
    enddo

    omega_matrix = omega_matrix - transpose(omega_matrix) 
    !============================================================================


    13 format(3es14.5E3)
    
    deallocate(H_temp) 

end subroutine calculate_eigenvectors


subroutine build_derivative_matrix(disx, disy, DerMtx)
    implicit none

    !args
    real*8, intent(in) :: disx(nsites, nsites) , disy(nsites, nsites)
    real*8, allocatable, intent(out) :: DerMtx(:,:)

    !local 
    real*8 :: DMtx
    integer :: i, j, k
    
    allocate( DerMtx(d_el, d_el), source = 0.d0 )



    ! === CONSTRUO OS BLOCOS DA MATRIZ =======

    do j = 1, nsites
     do i = 1, nsites
        call DerivativeOverlap(i, j, disx, disy, DMtx)
        DerMtx(i, j) = DMtx
      enddo
    enddo



end subroutine build_derivative_matrix



end module system_hamiltonian_m
