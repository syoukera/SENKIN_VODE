    module chemkin
    
    implicit none
      
    ! define length of work array
    !integer, parameter, private :: leniwk = 14271
    !integer, parameter, private :: lenrwk = 10177
    !integer, parameter, private :: lencwk = 10000
    integer :: leniwk = 14271
    integer :: lenrwk = 10177
    integer :: lencwk = 10000
    
    ! define unit number
    integer, parameter, private :: linck  = 25
    integer, parameter, private :: lout   = 16
    
    integer num_elem
    integer num_spec
    integer num_reac
    integer num_coef
    logical ierr
    
    ! declear work array
    integer,      allocatable :: ickwrk(:)
    real(8),      allocatable :: rckwrk(:)
    character(6), allocatable :: cckwrk(:)*16
    character(6), allocatable :: name_spec(:)*16
    
    ! integer 
      
    contains
      
      subroutine initialize()
        ! open input output files
        open(linck, file='cklink', form='unformatted')
        open(lout, file='skout', form='formatted')
        
        ! allocate work array
        call cklen (linck, lout, leniwk, lenrwk, lencwk)
        allocate(ickwrk(leniwk))
        allocate(rckwrk(lenrwk))
        allocate(cckwrk(lencwk))
        
        ! initialize chemkin data structure
        call ckinit(leniwk, lenrwk, lencwk, linck, lout,  &
                    ickwrk, rckwrk, cckwrk)
        call ckindx(ickwrk, rckwrk, num_elem, num_spec,   &
                    num_reac, num_coef)
        
        ! load species names
        allocate(name_spec(num_spec))
        call cksyms(cckwrk, lout, name_spec, ierr)
        
      end subroutine initialize
    
    end module chemkin