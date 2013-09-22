! Copyright (C) 2002, 2010 Carnegie Mellon University and others.
! All Rights Reserved.
! This code is published under the Eclipse Public License.
!
!    $Id: hs071_f.f.in 1861 2010-12-21 21:34:47Z andreasw $
!
! =============================================================================
!
!
!  This file contains routines to define a small Rosenbrock test problem.
!
!  f=0 at the optimal solution for x_i=1 for all i!
!
! =============================================================================
!
!
! =============================================================================
!
!                            Main driver program
!
! =============================================================================
!
      program example
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
      parameter  (N = 3, M = 2, NELE_JAC = 6, NELE_HESS = 6)
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
      double precision DAT(2)
      integer IDAT(1)
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
      double precision F
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

      DAT(1) = 1.0 !FS

!!$      do i=1,2
!!$         X(i)   =  200.0   !500 mm
!!$         X_L(i) =  10.0 !0 mm
!!$         X_U(i) =  250.0 !max 1 metre =1000 mm
!!$      end do


      !define R,T
      
!      R=x(1)
!      T=x(2)
!      L=x(3)


     

      x(1)  =220.0
      X_L(1)=0.0
      x_U(1)=200.0

      x(2)  =100.0
      X_L(2)=0.0
      x_U(2)=100.0

      x(3)=2000.0
      X_L(3)=500.0
      x_U(3)=5000.0
      
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
            write(*,*) 'X  (',i,') = ',X(i)
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
         write(*,*) 'This problem has infinitely many solutions--Komahan'
      else
         write(*,*)
         write(*,*) 'An error occoured.'
         write(*,*) 'The error code is ',IERR
         write(*,*)
      endif
!
9000  continue
!
!     Clean up
!
      call IPFREE(IPROBLEM)
      stop
!
9990  continue
      write(*,*) 'Error setting an option'
      goto 9000
    end program example
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
      real*8 :: rho, sigmay, pi, Fs, p, E,sigma_allow

      p=10.0e6               !10 MN
      E=207000.0               !N/mm2
      rho=7.833e-6           !kg/m3
      sigma_allow=248.0        !N/mm2
      pi=4.0*atan(1.0)
      Fs=dat(1)
!!$
!!$      R=x(1)
!!$      T=x(2)
!!$      L=x(3)
!!$      

!      x(1)=10.0
!      x(2)=20.0
!      x(3)=30.0

      f = 2.0*rho*x(3)*pi*x(1)*x(2)
!      print*,"objective",x,f

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
      real*8 :: rho,  sigmay, pi, Fs, p, E,sigma_allow
      real*8:: tau_allow,BM,V

      p=10.0e6               !10 MN
      E=207000.0               !N/mm2
      rho=7.833e-6           !kg/m3
      sigma_allow=248.0        !N/mm2
      pi=4.0*atan(1.0)
      Fs=dat(1)


!      x(1)=10.0
!      x(2)=20.0
!      x(3)=30.0

      grad(1) = 2.0*pi*rho*x(2)*x(3)
      grad(2) = 2.0*pi*rho*x(1)*x(3)
      grad(3) = 2.0*pi*rho*x(1)*x(2)
      
!      print*,'grad:',x,grad


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
      integer N, NEW_X, M
      double precision G(M), X(N)
      double precision DAT(*)
      integer IDAT(*)
      integer IERR
      real*8 :: rho,sigmay, pi, Fs, p, E,sigma_allow
      real*8::R,T,L
      p=10.0e6               !10 MN
      E=207000.0               !N/mm2
      rho=7.833e-6           !kg/m3
      sigma_allow=248.0        !N/mm2
      pi=4.0*atan(1.0)
      Fs=dat(1)
      
!      x(:)=200.0
      
!      x(1)=10.0
!      x(2)=20.0
!      x(3)=30.0
!!$
!!$      R=X(1)
!!$      T=X(2)
!!$      L=X(3)
      
      G(1) = P*Fs / (2.0*pi*x(1)*x(2)*sigma_allow) - 1.0
!      G(2) = (4.0*p*Fs*(L**2) / (T*E*(pi*R)**3)) - 1.0

      G(2)=(4.0*fs*P*x(3)**2/((pi**3)*E*x(2)*x(1)**3)) -1.0

!     print*,'Constraints:',x,g

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
      integer TASK, N, NEW_X, M, NZ
      double precision X(N), A(NZ)
      integer ACON(NZ), AVAR(NZ), I
      double precision DAT(*),dc(M,N)
      integer IDAT(*)
      integer IERR
      real*8 :: rho,sigmay, pi, Fs, p, E, sigma_allow

      p=10.0e6               !10 MN
      E=207000.0               !N/mm2
      rho=7.833e-6           !kg/m3
      sigma_allow=248.0        !N/mm2
      pi=4.0*atan(1.0)
      Fs=dat(1)

      if( TASK.eq.0 ) then 
         !
         !     structure of Jacobian:
         !     
         ACON(1) = 1
         AVAR(1) = 1

         ACON(2) = 1
         AVAR(2) = 2

         ACON(3) = 1
         AVAR(3) = 3

         ACON(4) = 2
         AVAR(4) = 1

         ACON(5) = 2
         AVAR(5) = 2

         ACON(6) = 2
         AVAR(6) = 3

      else

!         x(1)=10.0
!         x(2)=20.0
!         x(3)=30.0

         !---- GRADIENT OF CONSTRAINTS

         dc(:,:) =  0.0

         dc(1,1) = -P*FS / (2.0*pi*sigma_allow*(x(1)**2)*x(2))
         dc(1,2) = -P*FS / (2.0*pi*sigma_allow*x(1)*x(2)**2)
         dc(1,3) =  0.0

!         print*, 'Gr Cons 1',dc(1,:)


         dc(2,1) = -12.0*fs*P*x(3)**2 / (E*x(2)*(pi**3)*(x(1)**4))
         dc(2,2) = -4.0*P*FS*x(3)**2 / (E*(x(2)**2)*(pi**3)*(x(1)**3))
         dc(2,3) =  8.0*P*x(3)*FS / (E*x(2)*(pi**3)*(x(1)**3))

 !        print*, 'Gr Cons 2',dc(2,:)

         A(1)=dc(1,1)
         A(2)=dc(1,2)
         A(3)=dc(1,3)

         A(4)=dc(2,1)
         A(5)=dc(2,2)
         A(6)=dc(2,3)

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
      real*8 :: rho,sigmay, pi, Fs, p, E,sigma_allow !Truss design parameters
!      real*8::R,L,T
      p=10.0e6               !10 MN
      E=207000.0               !N/mm2
      rho=7.833e-6           !kg/m3
      sigma_allow=248.0        !N/mm2
      pi=4.0*atan(1.0)
      Fs=dat(1)

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

         IRNH(4) = 2
         ICNH(4) = 1

         IRNH(5) = 3
         ICNH(5) = 2

         IRNH(6) = 3
         ICNH(6) = 1
        
      else

!
!     calculate Hessian:
!

         ! Objective constraint
         
         objhess(1)=0.0
         objhess(2)=0.0
         objhess(3)=0.0
 
         objhess(4)=2.0*pi*rho*x(3)
         objhess(5)=2.0*pi*rho*x(1)

         objhess(6)=2.0*pi*rho*x(2)

 !        print*,''
 !        print*,'OBJ HESS:',objhess(:)
 !        print*,''

         ! First  constraint

         conhess(1,1)= P*FS/(pi*(x(1)**3)*x(2)*sigma_allow)
         conhess(1,2)= P*FS/(pi*x(1)*(x(2)**3)*sigma_allow)
         conhess(1,3)= 0.0

         conhess(1,4)= P*Fs/(2.0*pi*X(1)**2*X(2)**2*sigma_allow)
         conhess(1,5)= 0.0

         conhess(1,6)= 0.0

!         print*,'CONS HESS 1:',conhess(1,:)
!         print*,''
         ! SecondConstaint

         conhess(2,1)= 48.0*FS*(x(3)**2)*P/(E*x(2)*(pi**3)*(x(1)**5))
         conhess(2,2)= 8.0*FS*(x(3)**2)*P/(E*(pi*x(1)*x(2))**3)
         conhess(2,3)= 8.0*P*FS/(E*x(2)*(pi*x(1))**3)

         conhess(2,4)= (12.0*X(3)**2*P*FS)/(E*pi**3*X(1)**4*X(2)**2)

         conhess(2,5)= -8.0*x(3)*P*FS/((pi**3)*E*(x(1)**3)*(x(2)**2))

         conhess(2,6)= -(24.0*X(3)*P*FS)/(E*pi**3*X(1)**4*X(2))


!         print*,'CONS HESS 2:',conhess(2,:)
 !        print*,''
         ! Assemble
         
         HESS(:)=0.0
         do i=1,NNZH
            hesstmp=0.0
            do j=1,m
               hesstmp=hesstmp+lam(j)*conhess(j,i)
            end do
            hess(i)=hesstmp+objhess(i)
         end do
         
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
!      if (ITER_COUNT .gt. 1 .and. DNORM.le.1D-10) ISTOP = 1

      return
    end subroutine ITER_CB
