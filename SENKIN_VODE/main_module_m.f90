    module main_module
    
    use chemkin
    
    implicit none
    
    ! initial conditions 
    real(8), parameter, private :: temperature_ini = 800d0    ! [K]
    real(8), parameter, private :: pressure_Pa     = 101325d0 ! [Pa]
    integer, parameter, private :: num_ini = 3
    real(8)      :: y_ini(num_ini)        ! name
    character(6) :: name_ini(num_ini)*16  ! mass fraction [-]
    
    ! simulation conditions
    real(8), parameter, private :: time_end = 1.0d0   ! [s]
    real(8), parameter, private :: dt       = 1.0d-3  ! [s]
    
    ! define unit number
    integer, parameter, private :: mout   = 16
    
    ! dependent variables for VODE
    real(8) time, time_out, pressure_ck
    real(8), allocatable :: z(:)
    real(8), allocatable :: zdot(:)
    
    data y_ini /0.0511d0, 0.179d0, 0.769d0/
    data name_ini /'CH4', 'N2', 'O2'/
    
    CONTAINS
        
      subroutine mn_initialize()
      
        open(mout, file='mnout', form='formatted')
        
        ! linitialize chemkin module
        call ck_initialize()
        
        ! allocate dependent variables
        allocate(z(num_eqns))
        allocate(zdot(num_eqns))
        
        ! assign initial values
        call mn_assign_ini()
        
      end subroutine mn_initialize
      
      ! ---
      
      subroutine mn_assign_ini()
      
        integer i, idx_ith_spec
        character(6) name_ith_spec
        
        ! assign initial temperature
        z(1) = temperature_ini
        
        ! assign initial mass fraction
        do i = 1, num_ini
          name_ith_spec = name_ini(i)
          ! get index from species name
          call ckcomp(name_ith_spec, name_spec, num_spec, idx_ith_spec)
          z(idx_ith_spec+1) = y_ini(i)
        enddo
        
        ! convert pressure unit to dyne/cm2
        pressure_ck = 0.1d0*pressure_Pa
        
      end subroutine mn_assign_ini
      
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
      
      ! ---
      
      subroutine fex()
      
      real(8) rho      ! gm/cm**3
      real(8) c_pb     ! ergs/(gm*K) 
      real(8) volume   ! cm**3/gm
      real(8) :: wdot(num_spec)     ! moles/(cm**3*sec)
      real(8) :: enthalpy(num_spec) ! ergs/gm
      integer i, kk, iprck, iprd, ipwt, ipwdot, iph, ipick, lout, ii
      
        KK     = ickwrk(1) ! num_spec
        IPRCK  = ickwrk(2)
        IPRD   = ickwrk(3)
        IPWT   = ickwrk(6)
        IPWDOT = ickwrk(7)
        IPH    = ickwrk(8)
        IPICK  = ickwrk(9)
        LOUT   = ickwrk(10)
        II     = ickwrk(11) ! num_reac
        
        ! ---------- prepare phisical values ----------
        ! pertubation factor
        do i = 1, num_reac
          call ckrdex(-i, rckwrk, rckwrk(iprd+i-1))
        enddo
        
        ! get rho, c_pb, wdot, enthalpy
        call ckrhoy(pressure_ck, z(1), z(2:), ickwrk, rckwrk, rho) 
        call ckcpbs(z(1), z(2), ickwrk, rckwrk, c_pb)
        call ckwyp(pressure_ck, z(1), z(2), ickwrk, rckwrk, wdot)
        call ckhms(z(1), ickwrk, rckwrk, enthalpy)
        
        if (rho < 0.0) then
          write (mout,*) 'Stop, zero density.'
        endif
        volume = 1.0d0/rho
        
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

    end module main_module