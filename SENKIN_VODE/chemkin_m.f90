    module chemkin
    
    implicit none
      
    ! define length of work array
    integer, parameter, private :: leniwk = 14271
    integer, parameter, private :: lenrwk = 10177
    integer, parameter, private :: lencwk = 10000
    
    ! define unit number
    integer, parameter, private :: linck  = 25
    integer, parameter, private :: lout   = 16
    
    integer num_elem
    integer num_spec
    integer num_reac
    integer num_coef
    !integer i
    logical ierr
    
    ! declear work array
    integer ickwrk(leniwk)
    real(8) rckwrk(lenrwk)
    character(6) cckwrk(lencwk)
    character(6), allocatable :: name_spec(:)
    
    ! integer 
      
    contains
      
      subroutine initialize()
        ! open input output files
        open(linck, file='cklink', form='unformatted')
        open(lout, file='skout', form='formatted')
        
        ! initialize chemkin data structure
        call ckinit(leniwk, lenrwk, lencwk, linck, lout,  &
          ickwrk, rckwrk, cckwrk)
        call ckindx(ickwrk, rckwrk, num_elem, num_spec,   &
             num_reac, num_coef)
        
        allocate(name_spec(num_spec))
        
        call cksyms(cckwrk, lout, name_spec, ierr)
        
       write(lout,*) name_spec(:)
        
      end subroutine initialize
    
    end module chemkin