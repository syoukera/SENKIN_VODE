    program senkin_vode

      use dvode_f90_m
      use main_module
      
      call mn_initialize()
      
      call mn_solve()
      

      !IMPLICIT NONE
      !DOUBLE PRECISION ATOL, RTOL, T, TOUT, Y, RSTATS
      !INTEGER NEQ, ITASK, ISTATE, ISTATS, IOUT, IERROR, I
      !DIMENSION Y(3), ATOL(3), RSTATS(22), ISTATS(31)
      !
      !TYPE (VODE_OPTS) :: OPTIONS
      !
      !call ck_initialize()

!      OPEN (UNIT=6,FILE='example1.dat')
!      IERROR = 0
!      NEQ = 3
!      Y(1) = 1.0D0
!      Y(2) = 0.0D0
!      Y(3) = 0.0D0
!      T = 0.0D0
!      TOUT = 0.4D0
!      RTOL = 1.D-4
!      ATOL(1) = 1.D-8
!      ATOL(2) = 1.D-14
!      ATOL(3) = 1.D-6
!      ITASK = 1
!      ISTATE = 1
!!     OPTIONS = SET_OPTS(DENSE_J=.TRUE.,ABSERR_VECTOR=ATOL,RELERR=RTOL, &
!!       USER_SUPPLIED_JACOBIAN=.TRUE.)
!      OPTIONS = SET_NORMAL_OPTS(DENSE_J=.TRUE.,ABSERR_VECTOR=ATOL,      &
!        RELERR=RTOL,USER_SUPPLIED_JACOBIAN=.TRUE.)
!      DO IOUT = 1, 12
!        CALL DVODE_F90(FEX,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS,J_FCN=JEX)
!        CALL GET_STATS(RSTATS,ISTATS)
!        WRITE (6,90003) T, Y(1), Y(2), Y(3)
!        DO I = 1, NEQ
!          IF (Y(I)<0.0D0) IERROR = 1
!        END DO
!        IF (ISTATE<0) THEN
!          WRITE (6,90004) ISTATE
!          STOP
!        END IF
!        TOUT = TOUT*10.0D0
!      END DO
!      WRITE (6,90000) ISTATS(11), ISTATS(12), ISTATS(13), ISTATS(19), &
!        ISTATS(20), ISTATS(21), ISTATS(22)
!      IF (IERROR==1) THEN
!        WRITE (6,90001)
!      ELSE
!        WRITE (6,90002)
!      END IF
      
90000 FORMAT (/'  No. steps =',I4,'   No. f-s =',I4,'  No. J-s =',I4, &
        '   No. LU-s =',I4/'  No. nonlinear iterations =', &
        I4/'  No. nonlinear convergence failures =', &
        I4/'  No. error test failures =',I4/)
90001 FORMAT (/' An error occurred.')
90002 FORMAT (/' No errors occurred.')
90003 FORMAT (' At t =',D12.4,'   y =',3D14.6)
90004 FORMAT (///' Error halt: ISTATE =',I3)
      STOP
    end program senkin_vode
