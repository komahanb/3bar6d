! =============================================================================
!
!                            Main driver program
!
! =============================================================================


! Three Bar Truss Optimization
! Author : Komahan Boopathy, University of Dayton, OH 24/9/2013
! Thanks to IPOPT


program Threebar
  !
  implicit none
  !
  !     include the Ipopt return codes
  !
  include 'IpReturnCodes.inc'
  !
  !     Size of the problem (number of variables and equality constraints)
  !
  integer     N,     M,     NELE_JAC,     NELE_HESS,      IDX_STY
  parameter  (N = 6, M = 8, NELE_JAC = 48, NELE_HESS = 21)
  parameter  (IDX_STY = 1 )
  !
  !     Space for multipliers and constraints
  !
  double precision LAM(M)
  double precision G(M)
  !
  !     Vector of variables
  !
  double precision X(N)
  !
  !     Vector of lower and upper bounds
  !
  double precision X_L(N), X_U(N), Z_L(N), Z_U(N)
  double precision G_L(M), G_U(M)
  !
  !     Private data for evaluation routines
  !     This could be used to pass double precision and integer arrays untouched
  !     to the evaluation subroutines EVAL_*
  !
  double precision DAT(20)
  integer IDAT(20)
  !
  !     Place for storing the Ipopt Problem Handle
  !
  integer*8 IPROBLEM
  integer*8 IPCREATE
  !
  integer IERR
  integer IPSOLVE, IPADDSTROPTION
  integer IPADDNUMOPTION, IPADDINTOPTION
  integer IPOPENOUTPUTFILE
  !
  double precision F,pi
  integer i

  double precision  infbound
  parameter        (infbound = 1.d+20)

  !
  !     The following are the Fortran routines for computing the model
  !     functions and their derivatives - their code can be found further
  !     down in this file.
  !
  external EV_F, EV_G, EV_GRAD_F, EV_JAC_G, EV_HESS
  !!
  !!     The next is an optional callback method.  It is called once per
  !!     iteration.
  !!
  external ITER_CB
  !
  !     Set initial point and bounds:
  !

  pi=4.0*atan(1.0) ! constant for later use (visible globally)

  !Problem data and other constants
  dat(1)=10.0 !height ref
  dat(2)=1.0e7 !E
  dat(3)=0.1 !gamma
  dat(4)=45.0*pi/180.0
  dat(5)=20000.0

!!! Max constraint values

  !Tensile
  dat(6)=5000.0    ! psi tensile_sigma1_max=dat(6)      
  dat(7)=20000.0    ! psi tensile_sigma2_max=dat(7)
  dat(8)=5000.0    ! psi tensile_sigma3_max=dat(8)
  !Compressive
  dat(9)=5000.0    ! psi comp_sigma1_max=dat(9)
  dat(10)=20000.0   ! psi comp_sigma2_max=dat(10)
  dat(11)=5000.0   ! psi comp_sigma3_max=dat(11)
  !Displacement
  dat(12)=0.005    ! in  max_u_disp=dat(12)
  dat(13)=0.005    ! in  max_v_disp=dat(12)
  dat(14)=1.0      ! Factor of safety

  ! Area design variables

  do i=1,N-3
     X(i)   = 1.0  
     X_L(i) = 0.25 
     X_U(i) = infbound 
  end do

  ! orientation design variables

  ! phi(1)
  X(4)   = 45.0*pi/180.0
  X_L(4) = 30.0*pi/180.0
  X_U(4) = 60.0*pi/180.0

  !phi(2)
  X(5)   = 90.0*pi/180.0
  X_L(5) = 60.0*pi/180.0
  X_U(5) = 120.0*pi/180.0

  !phi(3)
  X(6)   = (90.0+45.0)*pi/180.0
  X_L(6) = (90.0+30.0)*pi/180.0
  X_U(6) = (90.0+60.0)*pi/180.0


  ! If you want to fix the orientation (phi's) uncomment the following
!!$
!!$      ! phi(1)
!!$      X(4)   = 45.0*pi/180.0
!!$      X_L(4) = 45.0*pi/180.0
!!$      X_U(4) = 45.0*pi/180.0
!!$      
!!$      !phi(2)
!!$      X(5)   = 90.0*pi/180.0
!!$      X_L(5) = 90.0*pi/180.0
!!$      X_U(5) = 90.0*pi/180.0
!!$      
!!$      !phi(3)
!!$      X(6)   = (90.0+45.0)*pi/180.0
!!$      X_L(6) = (90.0+45.0)*pi/180.0
!!$      X_U(6) = (90.0+45.0)*pi/180.0


  !
  !     Set bounds for the constraints
  !
  do i=1,M
     G_L(i)=-infbound
     G_U(i)=0.d0
  end do
  !
  !     First create a handle for the Ipopt problem (and read the options
  !     file)
  !


  IPROBLEM = IPCREATE(N, X_L, X_U, M, G_L, G_U, NELE_JAC, NELE_HESS,IDX_STY, EV_F, EV_G, EV_GRAD_F, EV_JAC_G, EV_HESS)
  if (IPROBLEM.eq.0) then
     write(*,*) 'Error creating an Ipopt Problem handle.'
     stop
  endif

  !
  !     Open an output file
  !
  IERR = IPOPENOUTPUTFILE(IPROBLEM, 'IPOPT.OUT', 50)
  if (IERR.ne.0 ) then
     write(*,*) 'Error opening the Ipopt output file.'
     goto 9000
  endif


  !
  !!
  !!     Set a callback function to give you control once per iteration.
  !!     You can use it if you want to generate some output, or to stop
  !!     the optimization early.
  !!
  call IPSETCALLBACK(IPROBLEM, ITER_CB)
  !
  !     As a simple example, we pass the constants in the constraints to
  !     the EVAL_C routine via the "private" DAT array.
  !

  !      DAT(2) = 0.d0
  !
  !     Call optimization routine
  !
  IERR = IPSOLVE(IPROBLEM, X, G, F, LAM, Z_L, Z_U, IDAT, DAT)

  !
  !     Output:
  !
  if( IERR.eq.IP_SOLVE_SUCCEEDED ) then
     write(*,*)
     write(*,*) 'The solution was found.'
     write(*,*)
     write(*,*) 'The final value of the objective function is ',F
     write(*,*)
     write(*,*) 'The optimal values of X are:'
     write(*,*)
     do i = 1, N
        if(i.GT.3) then
           write(*,*) 'X  (',i,') = ',X(i)*180.0/pi,'deg'
        else
           write(*,*) 'X  (',i,') = ',X(i),'in^2'
        end if
     enddo
!!$         write(*,*)
!!$         write(*,*) 'The multipliers for the lower bounds are:'
!!$         write(*,*)
!!$         do i = 1, N
!!$            write(*,*) 'Z_L(',i,') = ',Z_L(i)
!!$         enddo
!!$         write(*,*)
!!$         write(*,*) 'The multipliers for the upper bounds are:'
!!$         write(*,*)
!!$         do i = 1, N
!!$            write(*,*) 'Z_U(',i,') = ',Z_U(i)
!!$         enddo
     write(*,*)
     write(*,*) 'The multipliers for the equality constraints are:'
     write(*,*)
     do i = 1, M
        write(*,*) 'LAM(',i,') = ',LAM(i)
     enddo

  else
     write(*,*)
     write(*,*) 'An error occoured.'
     write(*,*) 'The error code is ',IERR
     write(*,*)
  endif
  !
9000 continue
  !
  !     Clean up
  !
  call IPFREE(IPROBLEM)
  stop
  !
9990 continue
  write(*,*) 'Error setting an option'
  goto 9000
end program Threebar
!
! =============================================================================
!
!                    Computation of objective function
!
! =============================================================================
!
      subroutine EV_F(N, X, NEW_X, F, IDAT, DAT, IERR)
      implicit none
      integer N, NEW_X,I
      double precision F, X(N)
      double precision DAT(*)
      integer IDAT(*)
      integer IERR

      call threebarf(0,dat,x,F)
          print *,' >> Obj:'
          print*,F
     
      IERR = 0

      return
    end subroutine EV_F
!
! =============================================================================
!
!                Computation of gradient of objective function
!
! =============================================================================
!
      subroutine EV_GRAD_F(N, X, NEW_X, GRAD, IDAT, DAT, IERR)
      implicit none
      integer N, NEW_X,i
      double precision GRAD(N), X(N)
      double precision DAT(*)
      integer IDAT(*)
      integer IERR

      call threebardf(0,dat,x,grad)
      print *,' >> Gradient of obj:'

      do i=1,n
         print*,grad(i)

      end do

      
      IERR = 0

      return
    end subroutine EV_GRAD_F
!
! =============================================================================
!
!                     Computation of equality constraints
!
! =============================================================================
!
      subroutine EV_G(N, X, NEW_X, M, G, IDAT, DAT, IERR)
      implicit none
      integer N, NEW_X, M,i
      double precision G(M), X(N)
      double precision DAT(*)
      integer IDAT(*)
      integer IERR

      do i=1,M
         call threebarf(i,dat,x,g(i))
      end do

      print*,''
      write(*,'(4x,a)') '>>Normalized Constraint Values:'

      do i=1,8
         write(*,'(E13.2)'),g(i)
      end do

      IERR = 0

      return
    end subroutine EV_G
    !
! =============================================================================
!
!                Computation of Jacobian of equality constraints
!
! =============================================================================
!
    subroutine EV_JAC_G(TASK, N, X, NEW_X, M, NZ, ACON, AVAR, A,IDAT, DAT, IERR)
      implicit none
      integer TASK, N, NEW_X, M, NZ
      double precision X(N), A(NZ)
      integer ACON(NZ), AVAR(NZ), I
      double precision DAT(*),cgrad(m,n)
      integer IDAT(*)
      integer IERR

      if( TASK.eq.0 ) then 
         !
         !     structure of Jacobian:
         !     


         ACON(1)  =1
         AVAR(1)  =1
         ACON(2)  =1
         AVAR(2)  =2
         ACON(3)  =1
         AVAR(3)  =3
         ACON(4)  =1 
         AVAR(4)  =4
         ACON(5)  =1
         AVAR(5)  =5
         ACON(6)  =1
         AVAR(6)  =6





         ACON(7)  =2 
         AVAR(7)  =1
         ACON(8)  =2
         AVAR(8)  =2
         ACON(9)  =2
         AVAR(9)  =3
         ACON(10) =2
         AVAR(10) =4
         ACON(11) =2
         AVAR(11) =5
         ACON(12) =2
         AVAR(12) =6



         ACON(13) =3
         AVAR(13) =1
         ACON(14) =3
         AVAR(14) =2
         ACON(15) =3
         AVAR(15) =3
         ACON(16) =3
         AVAR(16) =4
         ACON(17) =3
         AVAR(17) =5
         ACON(18) =3
         AVAR(18) =6



         ACON(19) =4
         AVAR(19) =1
         ACON(20) =4
         AVAR(20) =2
         ACON(21) =4
         AVAR(21) =3
         ACON(22) =4
         AVAR(22) =4
         ACON(23) =4
         AVAR(23) =5
         ACON(24) =4
         AVAR(24) =6



         ACON(25) =5
         AVAR(25) =1
         ACON(26) =5
         AVAR(26) =2
         ACON(27) =5
         AVAR(27) =3
         ACON(28) =5
         AVAR(28) =4
         ACON(29) =5
         AVAR(29) =5
         ACON(30) =5
         AVAR(30) =6




         ACON(31) =6
         AVAR(31) =1
         ACON(32) =6
         AVAR(32) =2
         ACON(33) =6
         AVAR(33) =3
         ACON(34) =6
         AVAR(34) =4
         ACON(35) =6
         AVAR(35) =5
         ACON(36) =6
         AVAR(36) =6




         ACON(37) =7
         AVAR(37) =1
         ACON(38) =7
         AVAR(38) =2
         ACON(39) =7
         AVAR(39) =3
         ACON(40) =7
         AVAR(40) =4
         ACON(41) =7
         AVAR(41) =5
         ACON(42) =7
         AVAR(42) =6




         ACON(43) =8
         AVAR(43) =1
         ACON(44) =8
         AVAR(44) =2
         ACON(45) =8
         AVAR(45) =3
         ACON(46) =8
         AVAR(46) =4
         ACON(47) =8
         AVAR(47) =5
         ACON(48) =8
         AVAR(48) =6


      else

         
         cgrad(:,:)=0.0 !initialize the jacobian

         do i=1,M
            call threebardf(i,dat,x,cgrad(i,1:6))
         end do
         

         ! Assemble
         A(1)=cgrad(1,1)
         A(2)=cgrad(1,2)
         A(3)=cgrad(1,3)
         A(4)=cgrad(1,4)
         A(5)=cgrad(1,5)
         A(6)=cgrad(1,6)

         A(7)=cgrad(2,1)
         A(8)=cgrad(2,2)
         A(9)=cgrad(2,3)
         A(10)=cgrad(2,4)
         A(11)=cgrad(2,5)
         A(12)=cgrad(2,6)

         A(13)=cgrad(3,1)
         A(14)=cgrad(3,2)
         A(15)=cgrad(3,3)
         A(16)=cgrad(3,4)
         A(17)=cgrad(3,5)
         A(18)=cgrad(3,6)

         A(19)=cgrad(4,1)
         A(20)=cgrad(4,2)
         A(21)=cgrad(4,3)
         A(22)=cgrad(4,4)
         A(23)=cgrad(4,5)
         A(24)=cgrad(4,6)

         A(25)=cgrad(5,1)
         A(26)=cgrad(5,2)
         A(27)=cgrad(5,3)
         A(28)=cgrad(5,4)
         A(29)=cgrad(5,5)
         A(30)=cgrad(5,6)

         A(31)=cgrad(6,1)
         A(32)=cgrad(6,2)
         A(33)=cgrad(6,3)
         A(34)=cgrad(6,4)
         A(35)=cgrad(6,5)
         A(36)=cgrad(6,6)

         A(37)=cgrad(7,1)
         A(38)=cgrad(7,2)
         A(39)=cgrad(7,3)
         A(40)=cgrad(7,4)
         A(41)=cgrad(7,5)
         A(42)=cgrad(7,6)

         A(43)=cgrad(8,1)
         A(44)=cgrad(8,2)
         A(45)=cgrad(8,3)
         A(46)=cgrad(8,4)
         A(47)=cgrad(8,5)
         A(48)=cgrad(8,6)

      end if

      IERR = 0
      return
    end subroutine EV_JAC_G

!
! =============================================================================
!
!                Computation of Hessian of Lagrangian
!
! =============================================================================
!
      subroutine EV_HESS(TASK, N, X, NEW_X, OBJFACT, M, LAM, NEW_LAM,NNZH, IRNH, ICNH, HESS, IDAT, DAT, IERR)
      implicit none
      integer TASK, N, NEW_X, M, NEW_LAM, NNZH, i, ir,j
      double precision X(N), OBJFACT, LAM(M), HESS(NNZH),OBJHESS(NNZH),CONHESS(M,NNZH)

      integer IRNH(NNZH), ICNH(NNZH)
      double precision DAT(*)
      integer IDAT(*)
      integer IERR
      double precision :: hesstmp


      if( TASK.eq.0 ) then
!
!     structure of sparse Hessian (lower triangle):
!

         IRNH(1) = 1
         ICNH(1) = 1

         IRNH(2) = 2
         ICNH(2) = 2

         IRNH(3) = 3
         ICNH(3) = 3

         IRNH(4) = 4
         ICNH(4) = 4

         IRNH(5) = 5
         ICNH(5) = 5

         IRNH(6) = 6
         ICNH(6) = 6

!!!!!!! diag 2
         IRNH(7) = 2
         ICNH(7) = 1

         IRNH(8) = 3
         ICNH(8) = 2
        
         IRNH(9) = 4
         ICNH(9) = 3
        
         IRNH(10) = 5
         ICNH(10) = 4
        
         IRNH(11) = 6
         ICNH(11) = 5

!!!!!   diag 3
         IRNH(12) = 3
         ICNH(12) = 1
        
         IRNH(13) = 4
         ICNH(13) = 2
        
         IRNH(14) = 5
         ICNH(14) = 3
        
         IRNH(15) = 6
         ICNH(15) = 4
        

!!!!! diag4

         IRNH(16) = 4
         ICNH(16) = 1
        
         IRNH(17) = 5
         ICNH(17) = 2
        
         IRNH(18) = 6
         ICNH(18) = 3
        

!!! diag5

         IRNH(19) = 5
         ICNH(19) = 1
        
         IRNH(20) = 6
         ICNH(20) = 2

!!! diag 6

       
         IRNH(21) = 6
         ICNH(21) = 1
         
 
      else
        
      IERR = 0
      
   endif
   return
 end subroutine EV_HESS
 !
! =============================================================================
 !
 !                   Callback method called once per iteration
 !
 ! =============================================================================
!
      subroutine ITER_CB(ALG_MODE, ITER_COUNT,OBJVAL, INF_PR, INF_DU,MU, DNORM, REGU_SIZE, ALPHA_DU, ALPHA_PR, LS_TRIAL, IDAT,DAT, ISTOP)
      implicit none
      integer ALG_MODE, ITER_COUNT, LS_TRIAL
      double precision OBJVAL, INF_PR, INF_DU, MU, DNORM, REGU_SIZE
      double precision ALPHA_DU, ALPHA_PR
      double precision DAT(*)
      integer IDAT(*)
      integer ISTOP

      if (ITER_COUNT .eq.0) then
         write(*,*) 
         write(*,*) 'iter    objective      ||grad||        inf_pr          inf_du         lg(mu)'
      end if

      write(*,'(i5,5e15.7)') ITER_COUNT,OBJVAL,DNORM,INF_PR,INF_DU,MU

      return
    end subroutine ITER_CB
