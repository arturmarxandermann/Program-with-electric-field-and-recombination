module types_m

type quantum_site
    real*8 :: radius
    real*8 :: radius0
    real*8 :: vel 
    real*8 :: mass
    real*8 :: vhomo
    real*8 :: vlumo 
    real*8 :: reorg_lumo
    real*8 :: reorg_homo
    real*8 :: omega
    real*8 :: omega0 
    real*8 :: xPos
    real*8 :: yPos
    real*8 :: gap
    integer :: dashtype 
    integer :: ovlptype
end type quantum_site

type obj_pointer
  type(quantum_site), pointer :: np 
end type obj_pointer

type BasisBuild
    real*8  :: hMtx
    real*8  :: DerMtx 
end type BasisBuild

type operators
    real*8, ALLOCATABLE, dimension(:,:) :: elements
end type operators


type tensor_product_matrices
    real*8, ALLOCATABLE, DIMENSION(:,:) :: elements
end type

end module types_m
