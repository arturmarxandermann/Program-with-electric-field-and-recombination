module parameters_m
use types_m
use constants_m 
!modulo que tem todos os parametros que usamos apenas no main.f90
!TODOS OS TEMPOS ESTAO EM UNIDADES DE PS E AS FREQUENCIAS EM UNIDADE DE THZ

    real*8,  parameter :: tmin = 0.d0                   !tempo inicial
    real*8,  parameter :: tmax = 10.d0                   !tempo final
    integer, parameter :: nm_divisoes = 20000            !numero de passos utilizado no programa
    integer, parameter :: wrcons = 20                    !escrevo o resultado a cada wrcons pontos
    
    integer, parameter :: nsites = 40                      !numero de sitios
    integer, parameter :: nr = 2                           !numero de linhas 
    integer, parameter :: nc = 20                           !numero de colunas 
    
    integer, parameter :: nanol = 1
    integer, parameter :: ndol = 9                      !numero de sitios doadores em uma linha
    integer, parameter :: nacl = 9                      !numero de sitios aceitadores em uma linha  
    integer, parameter :: ncal = 1

    integer, parameter :: ndes = 0!6                      !numero de pares com desordem
    integer, parameter :: ncdes = 0!4                     !numero de colunas com desordem

    integer, parameter :: ndo = 18                       !numero de sitios doadores totais
    integer, parameter :: nac = 18                       !numero de sitios aceitadores totais
    integer, parameter :: nano = nr                     !numero de drenos (buraco) 
    integer, parameter :: nca = nr                      !numero de drenos (elétron) 

    integer, parameter :: initState = 10               !estado inicial 

    integer, parameter :: nmdoc = ndo + nano
    integer, parameter :: nmacc = ndo + 2*nano

    real*8,  parameter :: ren_ho = 0.05d0             !energia de reorganizacao total em eV com reorg interna

    real*8, parameter :: el_field = 0.d0                 !campo eletrico em V/nanometro

    real*8,  parameter :: distance_int = 0.3d-9            !distancia da interface
    real*8,  parameter :: distance_molecule = 0.3d-9       !distancia entre os sitios
    real*8,  parameter :: scD = -1.468513d0*0.8d0 !1.527d0 !* 0.74d0          
                                                           
    real*8, parameter :: rec_rate = 5.d0
                                                               

    real*8, parameter :: forcecons = 1.d0

    logical, parameter :: BRD_on =     .true.          !true para Berendsen ligado e false para desligado 
    logical, parameter :: WithForces = .true.          !true para a força eletrica ligada 
    logical, parameter :: QD_on =      .true.          !true para o dissipador quantico ligado 
    logical, parameter :: Bound_on =   .false.          !true para a condição de contorno ligada
    logical, parameter :: ElHlInt =    .true.          !true para a interaçao el-hl ligada
    logical, parameter :: EField =     .false.
    logical, parameter :: SisDes =     .false.          !true para desordem no sistema
    
    logical, parameter :: Therm_on =   .false.          !true para termalizacao ligada

    logical, parameter :: FixSites = .false.             !true para evolucao com ri = raioZero
    logical, parameter :: FixTemp =  .false.             !true para temperatura fixa em bathTemp 

    logical, parameter :: RecombinationOn = .true.

    real*8, parameter :: bathTemp = 300.d0 
    real*8, parameter :: bathCoup = 1.d-13             !constante de acoplamento Berendsen em segundols

    real*8,  parameter :: siteMass = 1.9936d-26        !massa dos sitios
    real*8,  parameter :: effectiveMassMult = 1.d0      !1.6d0  !fator de multiplicacao da massa efetiva do eletron
    real*8, parameter :: freqMode = 2.d0 * PI * 8.15435d12 !14.8697d12 !frequencia do modo vibracional em Hz 
    real*8, parameter :: raioZero = 0.5d-9
    real*8, parameter :: ehcoup = 0.2d0   


    integer, parameter :: iidi = nanol + 1             !indice da interface doadora inicial
    integer, parameter :: ihd  = nanol + ndol          !indice da heterojuncao doadora
    integer, parameter :: iha  = nanol + ndol + 1      !indice da heterojuncao aceitadora
    integer, parameter :: iiaf = nanol + ndol + nacl   !indice da interface aceitadora final

    integer, parameter :: d_el = nsites         !tamanho da matriz rho, hamiltoniana
    integer, parameter :: neqn_el = d_el * d_el        !numero de equacoes diferenciais acopladas
    integer, parameter :: dcoup = ( d_el * ( d_el - 1)  ) / 2
   
    real*8,  parameter :: me = effectiveMassMult * electronMass !massa efetiva do eletron
    real*8,  parameter :: velZero  = sqrt( (kboltz * bathTemp * ev_to_joule ) / siteMass ) !vel inicial
    real*8, parameter :: ctFreq_ho = 40.d0 * 2.d0 * pi * 1.d12 !40.d0
    real*8, parameter :: freqZero = (2.d0 * hbar) / ( me * (raioZero)**2.0 ) 
    real*8, parameter :: enctFreq = HB_ev * freqZero

    !ESTRUTURA PARA O ÂNODO
    real*8, parameter :: ehomoAN = 0.d0 - enctFreq    !HOMO
    real*8, parameter :: elumoAN = 3.5d0 - enctFreq   !LUMO

    !ESTRUTURA PARA O DOADOR
    real*8,  parameter :: ehomoD = 0.3d0 - enctFreq                  
    real*8,  parameter :: elumoD = 1.5d0 - enctFreq                 
    
    !ESTRUTURA PARA O ACEITADOR
    real*8,  parameter :: ehomoA = 1.8d0 - enctFreq                   
    real*8,  parameter :: elumoA = 1.2d0 - enctFreq !0.3d0 - enctFreq !0.5d0 - enctFreq                 

    !ESTRUTURA PARA O CÁTODO
    real*8, parameter :: ehomoCA = 1.8d0 - enctFreq
    real*8, parameter :: elumoCA = 0.9d0 - enctFreq !0.9d0 - enctFreq


    real*8                          :: dt, writepoint  
    integer                         :: neqn, nmfn_t, nmfn_d, nmfn_h
    real*8,             allocatable :: h_mtx(:,:), rwMtx(:, :) 
    real*8,             allocatable :: sdist(:,:)
    real*8                          :: invdelta(d_el, d_el), krdelta(d_el, d_el)  
    real*8                          :: firstneighbors_diag(nsites, nsites), firstneighbors_hor(nsites, nsites)
    real*8                          :: firstneighbors(nsites, nsites) 
    type(quantum_site), target, allocatable :: site(:,:)
    type(obj_pointer)         , allocatable :: site_point(:)
    type(BasisBuild)          , allocatable :: basis(:,:)

    integer                            :: desi, desf 
    integer,               allocatable :: desd_r(:), desd_c(:)  !desordem doador rows, desordem doador columns
    integer,               allocatable :: desa_r(:), desa_c(:)  !desordem aceitador rows, desordem aceitador columns

    complex*16 :: recSb(d_el, d_el), recSb_hl(d_el, d_el)  
    integer :: mtxaux(d_el, d_el)
    real*8,                allocatable :: eterm(:)

end module parameters_m
