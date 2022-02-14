module rdftensor_m
use f95_precision
use lapack95
use blas95
use types_m
use parameters_m
use constants_m
use functions_m
use omp_lib 
contains

pure real function soma_elemento(i1, i2, mtx)
    !funcao para somar os elementos \sum_n Gamma(a, n, n, b)
    implicit none
    !args
    integer, INTENT(IN) :: i1, i2
    real*8, DIMENSION(d_el, d_el, d_el, d_el), INTENT(IN) :: mtx
    !local
    real*8 :: temporary
    integer :: i 



    temporary = 0.d0

    do i = 1, d_el
        temporary = temporary + mtx(i1, i, i, i2)
    enddo

    soma_elemento = temporary
end function soma_elemento


pure function SpecFunc(cutfreq, renergy, freq)
    implicit none
    !args
    real*8, intent(in) :: cutfreq, freq, renergy
    real*8 :: SpecFunc 
    !local
    real*8 :: omgRatio, ctFreqSq, freqSq, den, ctFreqthz
   
    
    !multiplico por 1d^-12 para levar para 1/ps
    ctFreqthz = cutFreq * 1.d-12
    ctFreqSq = ctFreqthz * ctFreqthz


    freqSq = freq * freq
    den = freqSq + ctFreqSq

    SpecFunc = 2.d0 * renergy * ( (freq * ctFreqthz)  / den ) !modelo de drude DEPOIS


end function SpecFunc



function BEDist(freq, insTemp)  
    implicit none
    !args
    real*8, intent(in) :: freq, insTemp
    real*8 :: BEDist
    !local
    real*8 :: Etemp, Efreq, Expratio, Eratio
    real*8 :: freqT


    freqT = freq
    !if(freq == 0.d0) freqT = 1.d-5

    Etemp = kboltz * insTemp
    Efreq = hb_ev_ps * freqT
    Eratio = Efreq/Etemp
    Expratio = exp(Efreq/Etemp) 

    BEDist = 1.d0/(Expratio - 1.d0) 


end function BEDist

subroutine create_rdftensor(eoh, rho, hMtx, insTemp, fMtx, phi, rdfTsr)
    implicit none
    !args
    integer, intent(in) :: eoh
    complex*16, intent(in) :: rho(d_el, d_el) 
    real*8, intent(in) :: hMtx(d_el, d_el) 
    real*8, intent(in) :: insTemp
    real*8, intent(in) :: fMtx(d_el, d_el) 
    real*8, intent(in) :: phi(d_el, d_el)
    real*8, allocatable, intent(out) :: rdfTsr(:, :, :, :) 
      
    !local
    real*8 :: hbarsq, Etemp, tsrsm, gfac, rspmul, ctFreqthz, rsp_ho_all
    real*8, dimension(:, :),        allocatable :: rsp_ho,  qmtx_ho !Radius system part
    real*8, dimension(:, :, :, :),  allocatable :: gpl, gmn, gplft_ho, gmnft_ho
    integer :: a, b, c, d
    real*8 :: energyratio_ho, hbwc_ho

    allocate( rsp_ho(d_el, d_el)                 , source = 0.d0 )
    allocate( gpl(d_el, d_el, d_el, d_el)     , source = 0.d0 )
    allocate( gmn(d_el, d_el, d_el, d_el)     , source = 0.d0 )
    allocate( gplft_ho(d_el, d_el, d_el, d_el)   , source = 0.d0 )
    allocate( gmnft_ho(d_el, d_el, d_el, d_el)   , source = 0.d0 )
    allocate( rdfTsr(d_el, d_el, d_el, d_el)  , source = 0.d0 )

    !hbarsq = hb_ev_ps * hb_ev_ps
    Etemp = kboltz * insTemp 
    hbwc_ho = hb_ev_ps * ctFreq_ho * 1.d-12
    energyratio_ho = Etemp/hbwc_ho





    !CALCULANDO A TRANSFORMADA DE FOURIER DOS GAMMAS - analítico se o limite superior é infinito
    do j = 1, d_el
        do i = 1, d_el
            if (i > j) then
                gplft_ho(:, :, i, j) = SpecFunc(ctFreq_ho, ren_ho, fMtx(i, j)) * BEDist(fMtx(i, j), insTemp)
                gmnft_ho(i, j, :, :) = SpecFunc(ctFreq_ho, ren_ho, fMtx(i, j)) * ( BEDist(fMtx(i, j), insTemp) + 1.d0 ) 
            endif    

            if (i < j) then 
                gplft_ho(:, :, i, j) = SpecFunc(ctFreq_ho, ren_ho, fMtx(j, i)) * ( BEDist(fMtx(j, i), insTemp) + 1.d0 ) !gmnft(j, i, :, :)
                gmnft_ho(i, j, :, :) = SpecFunc(ctFreq_ho, ren_ho, fMtx(j, i)) * BEDist(fMtx(j, i), insTemp) !gplft(:, :, j, i)
            endif


            if (i == j) then
                gplft_ho(:, :, i, i) = 2.d0 * energyratio_ho * ren_ho        !modelo de drude w* = 0
                gmnft_ho(i, i, :, :) = gplft_ho(:, :, i, i)                  !modelo de drude w* = 0
            endif 

        enddo
    enddo


    
    !CALCULANDO GAMMA MAIS (gpl) E GAMMA MENOS (gmn)  
    do d = 1, d_el
        do c = 1, d_el
            do b = 1, d_el
                do a = 1, d_el
                  rsp_ho_all = sum(phi(:, a) * phi(:, b) * phi(:, c) * phi(:, d))

                  gpl(a, b, c, d) = rsp_ho_all * gplft_ho(a, b, c, d) 

                  gmn(a, b, c, d) = rsp_ho_all * gmnft_ho(a, b, c, d) 
                enddo
            enddo
        enddo
    enddo

    gpl = gpl / hb_ev_ps
    gmn = gmn / hb_ev_ps

    !CALCULANDO OS ELEMENTOS DO TENSOR DE REDFIELD Rabcd
    do d = 1, d_el
        do c = 1, d_el
            do b = 1, d_el
                do a = 1, d_el
                    tsrsm = gpl(d, b, a, c) + gmn(d, b, a, c)
                    if (a /= c .AND. d /= b) then
                        rdfTsr(a, b, c, d) = tsrsm
                    endif
        
                    if (a /= c .AND.  d == b) then
                        rdfTsr(a, b, c, d) = tsrsm - soma_elemento(a, c, gpl)
                    endif

                    if (a == c .AND. d /= b) then
                        rdfTsr(a, b, c, d) = tsrsm - soma_elemento(d, b, gmn)
                    endif
       
                    if (a == c .AND. d == b) then
                        rdfTsr(a, b, c, d) = tsrsm - soma_elemento(a, c, gpl) - soma_elemento(d, b, gmn)
                    endif
    
                    
                enddo
            enddo
        enddo
    enddo  


deallocate(gpl, gmn, gplft_ho, gmnft_ho, rsp_ho) 
end subroutine create_rdftensor    


    
subroutine createRwMatrix(eoh, rho, hMtx, insTemp, fMtx, phi, recHb, rwMtx)
    !subrotina que cria a mtx Rw que relaciona os vetores rhodot e rho
    !na forma rho_dot = matmul(Rw, rho).
    !Rw é composta por 9 blocos de mtxes, que são numerados na seguinte ordem
    ! Rw = ( ( 1 ) ( 2 )  ( 3 )  )
    !      ( ( 4 ) ( 5 )  ( 6 )  )
    !      (  (7 ) ( 8 )  ( 9 )  )
    !Como o bloco 1, 4 e 7 relaciona elementos rhoreal_aa, rhoreal_ab, rhoimag_ab
    !com rhoreal_aa
    !seu tamanho é d_el. Os outros blocos tem tamanho nd pois relacionam os
    !elementos de antes com rhoreal_ab, rhoimag_ab


    implicit none
    !args
    integer, intent(in) :: eoh
    complex*16, intent(in) :: rho(d_el, d_el)
    real*8, intent(in) :: hMtx(d_el, d_el) 
    real*8, intent(in) :: insTemp
    REAL*8, INTENT(IN), DIMENSION(d_el, d_el) :: fMtx
    real*8, intent(in), dimension(d_el, d_el) :: phi
    complex*16, intent(in) :: recHb(d_el, d_el) 
    REAL*8, INTENT(OUT), ALLOCATABLE, DIMENSION(:, :) :: rwMtx
    !local
    INTEGER :: nd, d_elSq, aux1, i, j, k, l, n, m, p, a, b
    real*8, allocatable :: rdfTsr(:, :, :, :)
    TYPE(tensor_product_matrices), DIMENSION(9) :: Rw_blocks
    TYPE(tensor_product_matrices), DIMENSION(9) :: rw_rec
    !TRmtc são mtxes utilziadas no produto tensorial para montar Rw
    real*8    :: temp1, temp2


    nd = (d_el * ( d_el -1 ) ) / 2
    d_elSq = d_el * d_el
    aux1 = d_el + nd

        


    ALLOCATE(rwMtx(d_elSq, d_elSq), source = 0.d0 )


    ALLOCATE(Rw_blocks(1)%elements(d_el, d_el), source = 0.d0 )
    ALLOCATE(Rw_blocks(2)%elements(d_el, nd), source = 0.d0 )
    ALLOCATE(Rw_blocks(3)%elements(d_el, nd), source = 0.d0 )

    ALLOCATE(Rw_blocks(4)%elements(nd, d_el), source = 0.d0 )
    ALLOCATE(Rw_blocks(5)%elements(nd, nd), source = 0.d0 )
    ALLOCATE(Rw_blocks(6)%elements(nd, nd), source = 0.d0 )

    ALLOCATE(Rw_blocks(7)%elements(nd, d_el), source = 0.d0 )
    ALLOCATE(Rw_blocks(8)%elements(nd, nd), source = 0.d0 )
    ALLOCATE(Rw_blocks(9)%elements(nd, nd), source = 0.d0 )

    ALLOCATE(rw_rec(1)%elements(d_el, d_el), source = 0.d0 )
    ALLOCATE(rw_rec(2)%elements(d_el, nd), source = 0.d0 )
    ALLOCATE(rw_rec(3)%elements(d_el, nd), source = 0.d0 )

    ALLOCATE(rw_rec(4)%elements(nd, d_el), source = 0.d0 )
    ALLOCATE(rw_rec(5)%elements(nd, nd), source = 0.d0 )
    ALLOCATE(rw_rec(6)%elements(nd, nd), source = 0.d0 )

    ALLOCATE(rw_rec(7)%elements(nd, d_el), source = 0.d0 )
    ALLOCATE(rw_rec(8)%elements(nd, nd), source = 0.d0 )
    ALLOCATE(rw_rec(9)%elements(nd, nd), source = 0.d0 )


    if ( Therm_on == .true. .OR. QD_on == .false. ) then
      allocate (rdfTsr(d_el, d_el, d_el, d_el), source = 0.d0 ) 
    else
      ! ---------- CRIO OS ELEMENTOS DO TENSOR DE REDFIELD Rabcd -------
      call create_rdftensor(eoh, rho, hMtx, insTemp, fMtx, phi, rdfTsr)
      !-----------------------------------------------------------------
    endif


    !---- bloco 1 ------
    do j = 1, d_el
        do i = 1, d_el
            Rw_blocks(1)%elements(i, j) = rdfTsr(i, i, j, j) 
        enddo
    enddo
    !-----------------

    !--- bloco 2 ------
    n = 0
    do i = 1, d_el
        do j = 1, d_el
            if (j > i ) then
                n = n + 1
                do l = 1, d_el
                    Rw_blocks(2)%elements(l, n) = rdfTsr(l, l, i, j) + rdfTsr(l, l, j, i)
                enddo
            endif
        enddo
    enddo
    !----------------


    !--- bloco 3 -------
    !bloco 3 é composto de zeros, já definido com o source
    !-------------------

    !--- bloco 4 --------

    do n = 1, d_el
        k = 0
        do i = 1, d_el
            do j = 1, d_el
                if (j > i ) then
                    k = k + 1
                    Rw_blocks(4)%elements(k, n) = rdfTsr(j, i, n, n)
                endif
            enddo
        enddo
    enddo

    !--------------------


    !--- blocos 5, 9, 6, 8, bloco 7 é zero


    l = 0
    k = 0

    do i = 1, d_el
        do j = 1, d_el
            l = 0
            if (j > i) then
                k = k + 1

                do m = 1, d_el
                    do n = 1, d_el
                        if (n > m) then
                            l = l + 1
                            Rw_blocks(5)%elements(l, k) = rdfTsr(n, m, j, i) + rdfTsr(n, m, i, j)
                            Rw_blocks(9)%elements(l, k) = rdfTsr(n, m, j, i) - rdfTsr(n, m, i, j) 
                            if (l == k) then
                            !l == k => mtxes diagonais
                                Rw_blocks(6)%elements(l, l) = fMtx(n, m) 
                            endif
                        endif
                    enddo
                enddo
            endif
        enddo
    enddo


    Rw_blocks(8)%elements = - Rw_blocks(6)%elements 



    if ( RecombinationOn == .true. ) then 


        ! MONTANDO O PRIMEIRO BLOCO DA PARTE DE RECOMBINAÇÃO
        ! d/dt rho_aa = sum_(i = a ) 2 H_ii rho^re_ii 
        do i = 1, d_el
            rw_rec(1)%elements(i, i) = 2.d0 * recHb(i, i)
        enddo

        ! MONTANDO O SEGUNDO BLOCO DA PARTE DE RECOMBINAÇÃO
        ! d/dt rho^re_aa = sum_(i diff a) 2 H_ai rho^re_ai
        do a = 1, d_el
            do i = 1, d_el
                if ( i /= a ) then
                    rw_rec(2)%elements(a, mtxaux(i, a)) = 2.d0 * recHb(i, a)
                endif
            enddo
        enddo
       !---------------------------------------------------

       ! MONTANDO O QUARTO BLOCO DA PARTE DE RECOMBINAÇÃO
       ! d/dt rho^re_ab = sum_i ( H_ai delta_ib + H_ib delta_ai ) rho^re_ii
  
       k = 1
       do b = 1, d_el 
          do a = 1, d_el
              if ( a > b ) then
                  do i = 1, d_el 
                  temp1 = 0.d0
                  temp2 = 0.d0
                      if ( i == b ) then
                          temp1 = recHb(a, i) 
                      endif
                      if ( i == a ) then 
                          temp2 = recHb(b, i) 
                      endif
                      rw_rec(4)%elements(k, i) = temp1 + temp2
                  enddo
             k = k + 1
             endif
         enddo
       enddo


      ! MONTANDO O QUINTO E NONO BLOCO DA PARTE DE RECOMBINAÇÃO
      ! 5: d/dt rho^re_ab = sum_i (H_ai rho^re_ib + rho^re_ai H_ib )
      ! 9: d/dt rho^im_ab = sum_i (H_ai rho^im_ib + rho^im_ai H_ib )

      k = 1
      do b = 1, d_el
          do a = 1, d_el 
              if ( a > b ) then
                  do i = 1, d_el
                      if ( i > b ) then
                          rw_rec(5)%elements(k, mtxaux(b, i)) =  rw_rec(5)%elements(k, mtxaux(b, i)) + recHb(a, i)
                          rw_rec(9)%elements(k, mtxaux(b, i)) =  rw_rec(9)%elements(k, mtxaux(b, i)) + recHb(a, i)
                      endif
                    
                      if ( i < b ) then 
                          rw_rec(5)%elements(k, mtxaux(b, i)) =  rw_rec(5)%elements(k, mtxaux(b, i)) + recHb(a, i)
                          rw_rec(9)%elements(k, mtxaux(b, i)) =  rw_rec(9)%elements(k, mtxaux(b, i)) - recHb(a, i)
                      endif

                      if ( i < a ) then
                          rw_rec(5)%elements(k, mtxaux(a, i)) =  rw_rec(5)%elements(k, mtxaux(a, i)) + recHb(i, b)
                          rw_rec(9)%elements(k, mtxaux(a, i)) =  rw_rec(9)%elements(k, mtxaux(a, i)) + recHb(i, b)
                      endif

                      
                      if ( i > a ) then
                          rw_rec(5)%elements(k, mtxaux(a, i)) =  rw_rec(5)%elements(k, mtxaux(a, i)) + recHb(i, b)
                          rw_rec(9)%elements(k, mtxaux(a, i)) =  rw_rec(9)%elements(k, mtxaux(a, i)) - recHb(i, b)
                      endif
                  enddo
              k = k + 1 
              endif
          enddo
      enddo

    endif         


    !----- montandol a mtx total --------
    !rwMtx(1:d_el, 1:d_el) = rw_rec(1)%elements
    !rwMtx(1:d_el, d_el+1:aux1) = rw_rec(2)%elements
    !rwMtx(1:d_el, aux1+1:d_elSq) = rw_rec(3)%elements

    !rwMtx(d_el+1:aux1, 1:d_el) = rw_rec(4)%elements
    !rwMtx(d_el+1:aux1, d_el+1:aux1) = rw_rec(5)%elements
    !rwMtx(d_el+1:aux1, aux1+1:d_elSq) = rw_rec(6)%elements

    !rwMtx(aux1+1:d_elSq, 1:d_el) = rw_rec(7)%elements
    !rwMtx(aux1+1:d_elSq, d_el+1:aux1) = rw_rec(8)%elements
    !rwMtx(aux1+1:d_elSq, aux1+1:d_elSq) = rw_rec(9)%elements
    
    rwMtx(1:d_el, 1:d_el) = Rw_blocks(1)%elements + rw_rec(1)%elements
    rwMtx(1:d_el, d_el+1:aux1) = Rw_blocks(2)%elements + rw_rec(2)%elements
    rwMtx(1:d_el, aux1+1:d_elSq) = Rw_blocks(3)%elements + rw_rec(3)%elements

    rwMtx(d_el+1:aux1, 1:d_el) = Rw_blocks(4)%elements + rw_rec(4)%elements
    rwMtx(d_el+1:aux1, d_el+1:aux1) = Rw_blocks(5)%elements + rw_rec(5)%elements
    rwMtx(d_el+1:aux1, aux1+1:d_elSq) = Rw_blocks(6)%elements + rw_rec(6)%elements

    rwMtx(aux1+1:d_elSq, 1:d_el) = Rw_blocks(7)%elements + rw_rec(7)%elements
    rwMtx(aux1+1:d_elSq, d_el+1:aux1) = Rw_blocks(8)%elements + rw_rec(8)%elements
    rwMtx(aux1+1:d_elSq, aux1+1:d_elSq) = Rw_blocks(9)%elements + rw_rec(9)%elements


!print*, "inicio"
!call print_mat2(rwMtx, d_elSq, d_elSq)
!print*, "fim"
!stop 


deallocate(rdfTsr) 
end subroutine createRwMatrix




end module rdftensor_m
