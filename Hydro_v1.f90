module grid
  implicit none
  ! creates grid and keeps it in memory for use
  integer,parameter :: Ncell=20,Ndim=1,Nvar=3,Nsteps=200
  real(kind=8),parameter :: Xmin=0.d0,Xmax=1.d0,dx=(Xmax-Xmin)/Ncell,dt_i=0.01d0,T_tot=Nsteps*dt_i
  real(kind=8),save :: c,dt,t
  real(kind=8),dimension(1:Nvar,1:Ncell),save :: primVar,consVar
  character(len=16),save :: filename
  character(len=64),save :: path='/home/acharlet/4A CRAL/Code/Shared/Hydro/Resultats/'
  contains
    subroutine init()
      integer :: i
      c=1; dt=dt_i; t=0
      primVar(1,1:floor(Ncell/2.))=3.d0;
      primVar(1,floor(Ncell/2.):Ncell)=1.d0 ! initial densities
      primVar(2,:) = 0.d0  ! initial velocities
      primVar(3,:) = primVar(1,:) ! initial pressure     
      consVar(1,:) = primVar(1,:)   ! D = rho
      consVar(2,:) = 0.d0          ! S = rho*h*v , v=0
      consVar(3,:) = 4.d0*primVar(1,:) ! U = rho*h-p , h=5, rho=p

      write (filename, "('data0000.dat')")
      open(10,file=path//filename,status='new')
100   format (3e20.10)
      do i=1,Ncell
         write(10,100) primVar(1,i), primVar(2,i), primVar(3,i)
      end do
      close(10)
    end subroutine init

 ! equations d'evolution
    function timeDer(V,N) result(dV)
      real(kind=8),dimension(Nvar) :: V_N,V_Nm,V_Np,dV! V = (D,S,U), m indice-1, p indice+1
      real(kind=8),dimension(Nvar),intent(in) :: V
      real(kind=8) :: h,hm,hp,dD,dU,dS,dW
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
    
      h  = 0.2*(1.d0+4.d0*(V_N(3)/V_N(1)))
      hp = 0.2*(1.d0+4.d0*(V_Np(3)/V_Np(1)))
      hm = 0.2*(1.d0+4.d0*(V_Nm(3)/V_Nm(1)))
      dW = ((V_Np(2)**2/(V_Np(1)*hp)+V_Np(3)-V_Np(1)*hp)-(V_Nm(2)**2/(V_Nm(1)*hm)+V_Nm(3)-V_Nm(1)*hm))/(2.d0*dx)
      dD = (V_Np(2)/hp-V_Nm(2)/hm)/(2.d0*dx)
      dU = (V_Np(3)-V_Nm(3))/(2.d0*dx)
      dS = (V_Np(2)-V_Nm(2))/(2.d0*dx)
    
      dV(1) = -dD
      dV(2) = -dW
      dV(3) = -dS
    end function timeDer
      
end module grid

program simulation
  use grid
  implicit none
  integer :: i,j
  !real(dp),parameter :: gam=4/3        ! isentropic expansion factor gamma
  
  call init()
  i=1
  do while(t<=T_tot)
     print*, "Step ",i,", t=",t," completion: ",100*t/T_tot,"%"
     write (filename, "('data',I4.4,'.dat')") i
     open(10,file=path//filename,status='new')
     call evolTVDLF
     call cons2prim
100  format (3e20.10e3) 
     do j=1,Ncell
        write(10,100) primVar(1,j), primVar(2,j), primVar(3,j)
     end do
     i = i+1    
  end do
  
contains
! Runge Kutta 4
  function RK4(V_i,N) result(V_f)
    integer,intent(in) :: N
    real(kind=8),dimension(Nvar),intent(in) :: V_i
    real(kind=8),dimension(Nvar) :: k1,k2,k3,k4,V_f
    k1 = timeDer(V_i,N)
    k2 = timeDer(V_i+0.5d0*dt*k1,N)
    k3 = timeDer(V_i+0.5d0*dt*k2,N)
    k4 = timeDer(V_i+dt*k3,N)
    V_f = V_i + (dt/6.d0)*(k1+2.d0*k2+2.d0*k3+k4)
  end function RK4

! Finite Volume
  function fV1D(V_i,N) result(V_f)
    integer,intent(in) :: N
    real(kind=8),dimension(Nvar),intent(in) :: V_i
    real(kind=8),dimension(Nvar) :: V_f
    V_f = V_i - dt*timeDer(V_i,N)
  end function fV1D

! calcul des flux
  function flux(V) result(F)
    ! flux > 0 quand va vers i croissants
    real(kind=8) :: h
    real(kind=8),dimension(Nvar),intent(in) :: V
    real(kind=8),dimension(Nvar) :: F
    h =  0.2*(1.d0+4.d0*(V(3)/V(1)))
    F(1) = V(2)/h
    F(2) = V(2)**2/(h*V(1))+V(3)-h*V(1)
    F(3) = V(2)
  end function flux  

! variables primitives vers variables conservatives
  subroutine cons2prim
    real(kind=8) :: h
    integer :: i
    do i=1,Ncell
       h=0.2d0*(1.d0+4.d0*(consVar(3,i)/consVar(1,i)))
       primVar(1,i) = consVar(1,i)
       primVar(2,i) = consVar(2,i)/(h*consVar(1,i))
       primVar(3,i) = consVar(1,i)*h-consVar(3,i)
    end do
  end subroutine cons2prim

  subroutine evolTVDLF
    integer :: i,j,k
    real(kind=8) :: ctemp
    real(kind=8),dimension(Nvar) :: Dm,Dp,V,Vm,Vp,Vnew
    real(kind=8),dimension(1:Nvar,1:Ncell) :: dV,Vh
    real(kind=8),dimension(1:Nvar,1:Ncell+1) :: F,Vl,Vr
    
    ! vérification du pas de temps
    do i = 1,Ncell
       ctemp = sqrt(primVar(3,i)/primVar(1,i))
       if (ctemp>c) c=ctemp
    end do
    if (dt > dx/c) dt=0.8*dx/c
    if (dt < 0.8*dx/c) dt=0.8*dx/c
    
    ! calcul limiteur minmod
    do i=1,Ncell
       V=consVar(:,i)
       if (i==1) then
          Vm = V
          Vp = consVar(:,i+1)
       else if (i==Ncell) then
          Vm = consVar(:,i-1)
          Vp = V
       else
          Vm = consvar(:,i-1)
          Vp = consVar(:,i+1)
       end if
       Dm = V-Vm
       Dp = Vp-V
       do j=1,Nvar
          if (Dm(j)>=0 .and. Dp(j)>=0) then
             dV(j,i) = min(Dm(j),Dp(j))
          else if (Dm(j)<=0 .and. Dp(j)<=0) then
             dV(j,i) = -min(-Dm(j),-Dp(j))
          else
             dV(j,i)=0
          end if
       end do
    end do
    
     ! calcul "demi pas de temps"
     do i = 1,Ncell
        Vh(:,i) = V - dt/(2*dx)*(flux(V+0.5*dV(:,i))-(V-0.5*dV(:,i)))
     end do
       
     ! calcul valeur variables aux interfaces, k numéro interface entre i-1 et i
     Vl(:,1) = Vh(:,1)+0.5*dV(:,1) !case "0" = copie case 1
     Vr(:,Ncell+1) = Vh(:,Ncell)-0.5*dV(:,Ncell) !case "Ncell+1" = copie case Ncell
     do k=1,Ncell
        Vl(:,k+1) = Vh(:,k)+0.5*dV(:,k)
        Vr(:,k) = Vh(:,k+1)-0.5*dV(:,k+1)
     end do

     ! calcul flux aux interfaces
     do k = 1,Ncell+1
        F(:,k) = 0.5*(flux(Vl(:,k))+flux(Vr(:,k)))
     end do
       
     ! calcul nouvelle valeur
     do i=1,Ncell
        Vnew = V - (dt/dx)*(F(:,i+1)-F(:,i))
        consVar(:,i) = Vnew
     end do
     t=t+dt
   end subroutine evolTVDLF
             
       
! subroutine d'evolution, v1
  subroutine evol
    real(kind=8),dimension(1:Nvar,1:Ncell) :: temp
    real(kind=8) :: ctemp
    integer :: i,j
    ! check if condition dt<dx/c is verified
    do i = 1,Ncell
       ctemp = sqrt(primVar(3,i)/primVar(1,i))
       if (ctemp>c) c=ctemp
    end do
    if (dt > dx/c) dt=0.8*dx/c
    if (dt < 0.8*dx/c) dt=0.8*dx/c
    do i = 1,Ncell
       temp(:,i) = RK4(consVar(:,i),i) 
    end do
    do i =1, Ncell
       consVar(:,i)=temp(:,i)
    end do
    t=t+dt
  end subroutine evol
end program simulation
