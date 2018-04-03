program simulation
  implicit none
  integer, parameter:: dp=kind(0.d0)                   ! double precision
  real(dp),parameter :: Xmin=0.0_dp,Xmax=1.0_dp
  integer,parameter :: Ncell=100,Ndim=1
  real(dp),parameter :: dx=(Xmax-Xmin)/Ncell
  integer :: N
  integer,parameter :: Nvar=Ndim+2
  real,dimension(Ncell,Nvar) :: grid_primvar,grid_consvar,grid_fluxes
  !real(dp),parameter :: gam=4/3        ! isentropic expansion factor gamma
! test
! initialisation
  subroutine init
    implicit none

  end subroutine init

! Fonction pour alleger les calculs  
  function W1D(V,h) result(W)
    implicit none
    real,dimension(Nvar) :: V
    real :: W
    real,intent(in) :: h
    W = (V(3)*V(3)/(V(1)*h)+V(2)-V(1)*h)
  end function W1D
    
! equations d'evolution
  function evol1D(N) result(dtV)
    implicit none
    real,dimension(Nvar) :: V_N,V_nm,V_np,dtV ! V = (D,U,S), m indice-1, p indice+1
    real :: h,dU,dS,dW
    V_N = grid_consvar(N,:)
    V_Nm= grid_consvar(N-1,:)
    V_Np= grid_consvar(N+1,:)
    
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
  subroutine prim2cons(N)
    implicit none
    integer,intent(in) :: N
    real :: h
    real :: rho,P
    real,dimension(Ndim) :: v
    real :: D,U,S
    
    rho = grid_primvar(N,1)
    P = grid_primvar(N,2)
    v = grid_primvar(N,3)
    h = 1 + 4*P/rho   ! gam = 4/3
    
    D = rho
    U = rho*h-P
    S = rho*h*v

    grid_consvar(N,1) = D
    grid_consvar(N,2) = U
    grid_consvar(N,3) = v
  end subroutine prim2cons
  
! evolution des variables
  subroutine evol1D(N)
    implicit none
    real :: D,U,S
    real :: rho,P
    real,dimension(Ndim) :: v
    real,dimension(Ndim*Ndim) :: W
    
    rho = grid_primvar(N,1)
    P = grid_primvar(N,2)
    v = grid_primvar(N,3)
    h = 1 + gam/(gam-1)*P/rho
    W = rho*h*v*v

    
  end subroutine evol1D
end program simulation
