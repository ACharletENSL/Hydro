module grid
  implicit none
  ! creates grid and keeps it in memory for use
  integer,parameter :: Ncell=20,Ndim=1,Nvar=3,Nsteps=100
  real(kind=8),parameter :: Xmin=0.d0,Xmax=1.d0,dx=(Xmax-Xmin)/Ncell
  real(kind=8),dimension(1:Nvar,1:Ncell),save :: primVar,consVar
  character(len=16),save :: filename
  character(len=64),save :: path='/home/acharlet/4A CRAL/Code/Shared/Hydro/Resultats/'
  contains
    subroutine init()
      integer :: i
      primVar(1,1:floor(Ncell/2.))=3.d0;
      primVar(1,floor(Ncell/2.):Ncell)=1.d0 ! initial densities
      primVar(2,:)=primVar(1,:) ! initial pressure
      primVar(3,:)=0.d0  ! initial velocities
      
      consVar(1,:) = primVar(1,:)   ! D = rho
      consVar(2,:) = 4.d0*primVar(1,:) ! U = rho*h-p , h=5, rho=p
      consVar(3,:) = 0.d0          ! S = rho*h*v , v=0

      write (filename, "('data000.dat')")
      open(10,file=path//filename,status='new')
100   format (3e20.10)
      do i=1,Ncell
         write(10,100) primVar(1,i), primVar(2,i), primVar(3,i)
      end do
      close(10)
    end subroutine init

    ! equations d'evolution
    function evol1D(V,N) result(dV)
      real(kind=8),dimension(Nvar) :: V_N,V_Nm,V_Np,dV! V = (D,U,S), m indice-1, p indice+1
      real(kind=8),dimension(Nvar),intent(in) :: V
      real(kind=8) :: h,dU,dS,dW
      integer,intent(in) :: N

      V_N = V
      ! outflow conditions
      if (N==1) then
         V_Nm=V
         V_Np= consVar(:,N+1)
      elseif (N==Ncell) then
          V_Nm= consVar(:,N-1)
          V_Np= V
      else
         V_Nm= consVar(:,N-1)
         V_Np= consVar(:,N+1)
      end if
    
      h  = 0.2*(1.d0+4.d0*(V_N(2)/V_N(1)))
      dW = ((V_Np(3)*V_Np(3)/(V_Np(1)*h)+V_Np(2)-V_Np(1)*h)-(V_Nm(3)*V_Nm(3)/(V_Nm(1)*h)+V_Nm(2)-V_Nm(1)*h))/(2.d0*dx)
      dU = (V_Np(2)-V_Nm(2))/(2.d0*dx)
      dS = (V_Np(3)-V_Nm(3))/(2.d0*dx)
    
      dV(1) = -dS/h
      dV(2) = -dS
      dV(3) = -dW
    end function evol1D    
end module grid

program simulation
  use grid
  implicit none
  integer :: i,j
  !real(dp),parameter :: gam=4/3        ! isentropic expansion factor gamma
  
  call init()
  do i=1,Nsteps
     write (filename, "('data',I3.3,'.dat')") i
     open(10,file=path//filename,status='new')
     call evolution
     call cons2prim
100  format (3e20.10e3) 
     do j=1,Ncell
        write(10,100) primVar(1,j), primVar(2,j), primVar(3,j)
     end do
  end do
  
contains
! Runge Kutta 4
  function RK4(V_i,N,dt) result(dV)
    integer,intent(in) :: N
    real(kind=8),intent(in) :: dt
    real(kind=8),dimension(Nvar),intent(in) :: V_i
    real(kind=8),dimension(Nvar) :: k1,k2,k3,k4,dV
    k1 = evol1D(V_i,N)
    k2 = evol1D(V_i+0.5d0*dt*k1,N)
    k3 = evol1D(V_i+0.5d0*dt*k2,N)
    k4 = evol1D(V_i+dt*k3,N)
    dV = V_i + (dt/6.d0)*(k1+2.d0*k2+2.d0*k3+k4)
  end function RK4

! variables primitives vers variables conservatives
  subroutine cons2prim
    real(kind=8) :: h
    integer :: i
    do i=1,Ncell
       h=0.2d0*(1.d0+4.d0*(consVar(2,i)/consVar(1,i)))
       primVar(1,i) = consVar(1,i)
       primVar(2,i) = consVar(2,i)-h
       primVar(3,i) = consVar(3,i)/(h*consVar(1,i))
    end do
  end subroutine cons2prim

! subroutine d'evolution
  subroutine evolution
    real(kind=8),dimension(Nvar) :: temp
    real(kind=8),parameter :: dt=0.01d0
    integer :: i
    do i = 1,Ncell
       temp(:) = RK4(consVar(:,i),i,dt)
       consVar(:,i)=temp(:)
    end do
  end subroutine evolution
end program simulation
