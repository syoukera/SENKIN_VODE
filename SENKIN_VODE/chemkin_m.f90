    module chemkin
    
    implicit none
      
    ! define length of work array
    integer :: leniwk
    integer :: lenrwk
    integer :: lencwk
    
    ! define unit number
    integer, parameter, private :: linck  = 25
    integer, parameter, private :: lout   = 16
    
    !total number of values
    integer num_elem
    integer num_spec
    integer num_reac
    integer num_eqns
    integer num_coef
    logical ierr
    
    ! work array
    integer,      allocatable :: ickwrk(:)
    real(8),      allocatable :: rckwrk(:)
    character(6), allocatable :: cckwrk(:)*16
    
    ! species value
    character(6), allocatable :: name_spec(:)*16
    real(8),      allocatable :: weight_spec(:)
        
    contains
      
      subroutine ck_initialize()
        ! open input output files
        open(linck, file='cklink', form='unformatted')
        open(lout, file='skout', form='formatted')
        
        ! initialize work array
        call ck_init_wrk()
        
        ! get parameters of species
        call ck_get_params()
        
      end subroutine ck_initialize
      
      ! ---
      
      subroutine ck_init_wrk()
      
        ! allocate work array
        call cklen (linck, lout, leniwk, lenrwk, lencwk)
        allocate(ickwrk(leniwk))
        allocate(rckwrk(lenrwk))
        allocate(cckwrk(lencwk))
        
        ! initialize chemkin data structure
        call ckinit(leniwk, lenrwk, lencwk, linck, lout,  &
                    ickwrk, rckwrk, cckwrk)        
      end subroutine ck_init_wrk
      
      ! ---
      
      subroutine ck_get_params()
        ! get total number of values
        call ckindx(ickwrk, rckwrk, num_elem, num_spec,   &
                    num_reac, num_coef)
        num_eqns = num_spec + 1
      
        ! get species names
        allocate(name_spec(num_spec))
        call cksyms(cckwrk, lout, name_spec, ierr)
        
        ! get molecular weights
        allocate(weight_spec(num_spec))
        call ckwt(ickwrk, rckwrk, weight_spec)
        
      end subroutine ck_get_params
      
      ! ---
      
!      subroutine ck_fex(time, z, zdot, ickwrk, rckwrk)
!        IPRCK  = IPAR(2)
!        IPRD   = IPAR(3)
!        IPWDOT = IPAR(7)
!        IPH    = IPAR(8)
!        IPICK  = IPAR(9)
!        LOUT   = IPAR(10)
!        II     = IPAR(11)
!      
!      end subroutine ck_fex
!      
!      SUBROUTINE RCONP_CHEM (TIME, Z, ZP, DELTA, IRES, RPAR, IPAR)
!C
!      COMMON /RES1/ P
!C
!C  Residual of differential equations for constant pressure case
!C
!C  Variables:
!C    Z(1)   = temperature (Kelvin)
!C    Z(K+1) = species mass fractions
!C    P      = pressure (dyne/cm2) - constant in time
!C    RHO    = density (gm/cm3)
!C    RPAR   = array of reaction pre-exponential constants
!C
!      KK     = IPAR(1)
!      IPRCK  = IPAR(2)
!      IPRD   = IPAR(3)
!      IPWT   = IPAR(6)
!      IPWDOT = IPAR(7)
!      IPH    = IPAR(8)
!      IPICK  = IPAR(9)
!      LOUT   = IPAR(10)
!      II     = IPAR(11)
!C
!C        MODIFY CHEMKIN WORK ARRAY FOR PRE-EXPONENTIAL
!C
!      DO 15 I = 1, II
!        CALL CKRDEX (-I, RPAR(IPRCK), RPAR(IPRD+I-1))
!15    CONTINUE
!C
!C         CALL CHEMKIN SUBROUTINES
!C
!      CALL CKRHOY (P, Z(1), Z(2), IPAR(IPICK), RPAR(IPRCK), RHO)
!      CALL CKCPBS (Z(1), Z(2), IPAR(IPICK), RPAR(IPRCK), CPB)
!      CALL CKWYP  (P, Z(1), Z(2), IPAR(IPICK), RPAR(IPRCK),
!     1             RPAR(IPWDOT))
!      CALL CKHMS  (Z(1), IPAR(IPICK), RPAR(IPRCK), RPAR(IPH))
!      IF (RHO .EQ. 0.0) THEN
!         WRITE (LOUT, '(/1X,A)') 'Stop, zero density in RCONP_ZDP.'
!         STOP
!      ENDIF
!      VOLSP = 1. / RHO
!C
!C         ENERGY EQUATION
!C
!      SUM = 0.
!      DO 100 K = 1, KK
!         K1 = K-1
!         SUM = SUM + RPAR(IPH+K1) * RPAR(IPWDOT+K1) * RPAR(IPWT+K1)
! 100  CONTINUE
!      DELTA(1) = ZP(1) + VOLSP *SUM /CPB
!C
!C         SPECIES EQUATIONS
!C
!      DO 200 K = 1, KK
!         K1 = K-1
!         DELTA(K+1) = ZP(K+1) - RPAR(IPWDOT+K1) *RPAR(IPWT+K1) *VOLSP
! 200  CONTINUE
!C
!      RETURN
!      END
      
      
      
    end module chemkin