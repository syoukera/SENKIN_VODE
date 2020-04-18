    module chemkin
    
    implicit none
      
    ! define length of work array
    integer :: leniwk
    integer :: lenrwk
    integer :: lencwk
    
    ! define unit number
    integer, parameter, private :: linck  = 25
    integer, parameter, private :: lout   = 16
    
    ! define total number of the values
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
      
      subroutine ck_initialize()
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
        
      end subroutine ck_initialize
    
    end module chemkin