module grid
  implicit none
  ! creates grid, initialise it and keeps it in memory for use
  integer, parameter,save:: dp=kind(0.d0)                   ! double precision
  integer,parameter,save :: Ncell=10,Ndim=1
  !real(dp),parameter,save :: Xmin=0.0_dp,Xmax=1.0_dp,dx=(Xmax-Xmin)/Ncell
  integer,parameter,save :: Nvar=Ndim+2
  real,dimension(Nvar,Ncell),save :: primVar,consVar
  integer :: i
  character(len=32) :: filename

  primVar(1,1:50)=3.0_dp; primVar(1,51:100)=1.0_dp ! initial densities
  primVar(2,:)=primVar(1,:) ! initial pressure
  primVar(3,:)=0.dp0  ! initial velocities

  consVar(1,:) = primVar(1,:)   ! D = rho
  consVar(2,:) = 4*primVar(1,:) ! U = rho*h-p , h=5, rho=p
  consVar(3,:) = 0.dp0          ! S = rho*h*v , v=0

  write (filename, "data000.dat")
  open(10,file=filename,status='new')
  100 format (3e20.10)
  do i=1,Ncell
     write(10,100) primVar(1,i) primVar(2,i) primVar(3,i)
  end do
  close(10)
  
  private
  public,only:primVar,consVar,Nvar,Ndim,Ncell
end module grid

program simulation
  use grid,only:primVar,consVar,Nvar,Ndim,Ncell
  implicit none
  integer :: i,j
  !real(dp),parameter :: gam=4/3        ! isentropic expansion factor gamma
  
! Fonction pour alleger les calculs  
  function W1D(V,h) result(W)
    implicit none
    real,dimension(Nvar) :: V
    real :: W
    real,intent(in) :: h
    W = (V(3)*V(3)/(V(1)*h)+V(2)-V(1)*h)
  end function W1D

 ! calcul de h avec variables primitives
  function hPrim(V) result(h)
    implicit none
    real,dimension(Nvar) :: V
    real :: h
    h = 1 + 4*(V(2)/V(1))
  end function hPrim
    
! equations d'evolution
  function evol1D(N) result(dtV)
    implicit none
    real,dimension(Nvar) :: V_N,V_nm,V_np,dtV ! V = (D,U,S), m indice-1, p indice+1
    real :: h,dU,dS,dW
    integer :: N
    V_N = consVar(:,N)
    V_Nm= consvar(:,N-1)
    V_Np= consvar(:,N+1)
    
    h  = 0.2*(1+4*(V_N(2)/V_N(1)))
    dW = (W1D(V_Np,h)-W1D(V_N,h))/dx
    dU = (V_Np(2)-V_Nm(2))/dx
    dS = (V_Np(3)-V_Nm(3))/dx
    
    dtV(1) = -dS/h
    dtV(2) = -dW
    dtV(3) = -dS
  end function evol1D 
  
! Runge Kutta 4
  function RK4(V_i,dt) result(V_f)
    implicit none
    real,intent(in) :: dt
    real,dimension(Nvar) :: V_i, V_f, V_t,k1,k2,k3,k4
    k1 = evol1D(V_i)
    k2 = evol1D(V_i+0.5*dt*k1)
    k3 = evol1D(V_i+0.5*dt*k2)
    k4 = evol1D(v_i+dt*k3)
    V_f = V_i + (dt/6.0)*(k1+2*k2+2*k3+k4)
  end function RK4

! variables primitives vers variables conservatives
  subroutine cons2prim
    implicit none
    real,dimension(Ncell) :: h

    h(:) = 0.2*(1+4*(consVar(2,:)/consVar(1,:))
    primVar(1,:) = consVar(1,:) ! rho = D
    primVar(2,:) = consVar(2,:)-h
    
  end subroutine cons2prim

! subroutine d'evolution
  subroutine evolution
    implicit none
    real,dimension(Nvar) :: temp
    real,parameter :: dt=0.01
    integer :: i
    do i = 1,Ncell
       temp(:) = RK4(consVar(:,i),dt)
       consVar(:,i)=temp(:)
    end do
  end subroutine evolution

  do i=1,10
     write (filename, "('data',I3.3,'.dat')") i
     open(10,file=filename,status='new')
     call evolution
     call cons2prim
     do j=1,Ncell
        write(10,100) primVar(1,i) primVar(2,i) primVar(3,i)
     end do

end program simulation
