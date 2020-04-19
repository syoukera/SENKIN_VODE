    module main_module
    
    use chemkin
    
    ! initial conditions 
    real(8), parameter, private :: temperature_ini = 800d0    ! [K]
    real(8), parameter, private :: pressure        = 101325d0 ! [Pa]
    real(8), allocatable        :: y_ini(:)                   ! mass fraction [-]
    
    ! simulation conditions
    real(8), parameter, private :: time_end = 1.0d0   ! [s]
    real(8), parameter, private :: dt       = 1.0d-3  ! [s]
    
    ! define unit number
    integer, parameter, private :: mout   = 16
    
    ! dependent variables for VODE
    real(8) time, time_out
    real(8), allocatable :: z(:)
    real(8), allocatable :: zdot(:)
    
    CONTAINS
    
      subroutine fex(neq, time, z, zdot)
      
      end subroutine

      !SUBROUTINE FEX(NEQ,T,Y,YDOT)
      !  IMPLICIT NONE
      !  INTEGER, INTENT (IN) :: NEQ
      !  DOUBLE PRECISION, INTENT (IN) :: T
      !  DOUBLE PRECISION, INTENT (IN) :: Y(NEQ)
      !  DOUBLE PRECISION, INTENT (OUT) :: YDOT(NEQ)
      !  YDOT(1) = -.04D0*Y(1) + 1.D4*Y(2)*Y(3)
      !  YDOT(3) = 3.E7*Y(2)*Y(2)
      !  YDOT(2) = -YDOT(1) - YDOT(3)
      !  RETURN
      !END SUBROUTINE FEX

      !SUBROUTINE JEX(NEQ,T,Y,ML,MU,PD,NRPD)
      !  IMPLICIT NONE
      !  INTEGER, INTENT (IN) :: NEQ, ML, MU, NRPD
      !  DOUBLE PRECISION, INTENT (IN) :: T
      !  DOUBLE PRECISION, INTENT (IN) :: Y(NEQ)
      !  DOUBLE PRECISION, INTENT (OUT) :: PD(NRPD,NEQ)
      !  PD(1,1) = -.04D0
      !  PD(1,2) = 1.D4*Y(3)
      !  PD(1,3) = 1.D4*Y(2)
      !  PD(2,1) = .04D0
      !  PD(2,3) = -PD(1,3)
      !  PD(3,2) = 6.E7*Y(2)
      !  PD(2,2) = -PD(1,2) - PD(3,2)
      !  RETURN
      !END SUBROUTINE JEX
      
      ! ---
      
      subroutine mn_initialize()
      
        open(mout, file='mnout', form='formatted')
        
        ! linitialize chemkin module
        call ck_initialize()
        
        ! allocate dependent variables
        allocate(z(num_eqns))
        allocate(zdot(num_eqns))
    
      end subroutine mn_initialize
      
      ! ---
        
      subroutine mn_solve()
      
        time = 0.0d0
        
        do while (time < time_end)
          time_out = time + dt
          
          call mn_timestep()
          
          time = time_out
        enddo
      
      end subroutine mn_solve
      
      ! ---
        
      subroutine mn_timestep()
        
        !call dvode_f90(fex, num_eqns, z, time, time_out,        &
        !               itask, istate, options, j_fcn=jex))
        write(mout, *) 'time = ', time
      
      end subroutine mn_timestep

    end module main_module