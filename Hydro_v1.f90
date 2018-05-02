module grid
  implicit none
  ! creates grid
  integer,parameter :: Ncell=1000,Ndim=1,Nvar=3,Nsteps=1000
  real(kind=8),parameter :: dx=0.01d0,dt_i=0.001d0,T_tot=Nsteps*dt_i
  real(kind=8),save :: c,dt,t
  real(kind=8),dimension(1:Nvar,1:Ncell),save :: primVar,consVar
   real(kind=8),dimension(1:Nvar,1:Ncell+1),save :: gridFlux,gridDiss,consR,consL
  character(len=16),save :: filename1,filename2,filename3
  character(len=128),save :: path='/media/acharlet/Data/Arthur/Documents/Cours/4A_CRAL/Code/Perso/Hydro/Resultats/',FileName
  
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

      write (filename1, "('data0000.dat')")
      FileName=trim(adjustl(path))//trim(adjustl(filename1))
      open(10,file=FileName,status='new')
100   format (3e20.10)
      do i=1,Ncell
         !write(10,100) consVar(1,i), consVar(2,i), consVar(3,i)
         write(10,100) primVar(1,i), primVar(2,i), primVar(3,i)
      end do
      close(10)
    end subroutine init

  ! variables primitives vers variables conservatives
  function cons2prim(Vc) result(Vp)
    real(kind=8),dimension(Nvar),intent(in) :: Vc
    real(kind=8),dimension(Nvar) :: Vp
    Vp(1) = Vc(1)
    Vp(3) = (Vc(3)-Vc(1))/3.d0
    Vp(2) = Vc(2)/(Vc(1)+4.d0*Vp(3))
  end function cons2prim

  ! variables conservatives vers variables primitives
  function prim2cons(Vp) result(Vc)
    real(kind=8),dimension(Nvar),intent(in) :: Vp
    real(kind=8),dimension(Nvar) :: Vc
    Vc(1) = Vp(1)
    Vc(2) = Vp(2)*(Vp(1)+4*Vp(1))
    Vc(3) = Vp(1)+3*Vp(3)
  end function prim2cons
  
end module grid

program simulation
  use grid
  implicit none
  integer :: iStep,i,l  
  call init()
  do iStep=1,Nsteps
     !print*, "Step ",iStep,", t=",t," completion: ",100*t/T_tot,"%"
     write (filename1, "('data',I4.4,'.dat')") iStep
     FileName=trim(adjustl(path))//trim(adjustl(filename1))
     open(10,file=FileName,status='new')
     call evol2t2s
     do l = 1,Ncell
        primVar(:,l)=cons2prim(consVar(:,l))
     end do
100  format (3e20.10e3) 
     do i=1,Ncell
        write(10,100) primVar(1,i), primVar(2,i), primVar(3,i)
        !write(10,100) consVar(1,i), consVar(2,i), consVar(3,i)
     end do
     
     !write (filename2, "('flux',I3.3,'.dat')") iStep
     !FileName=trim(adjustl(path))//trim(adjustl(filename2))
     !open(11,file=FileName,status='new')  
     !do i=1,Ncell+1
     !   write(11,100) gridFlux(1,i), gridFlux(2,i), gridFlux(3,i)
     !end do
     !write (filename3, "('diss',I3.3,'.dat')") iStep
     !FileName=trim(adjustl(path))//trim(adjustl(filename3))
     !open(12,file=FileName,status='new')
     !do i=1,Ncell+1 
     !   write(12,100) gridDiss(1,i),  gridDiss(2,i), gridDiss(3,i)
     !end do
  end do
  
contains
! calcul des flux
  function flux(Vc) result(F)
    ! flux > 0 quand va vers i croissants
    real(kind=8),dimension(Nvar),intent(in) :: Vc
    real(kind=8),dimension(Nvar) ::Vp,F
    Vp = cons2prim(Vc)
    F(1) = Vc(1)*Vp(2)
    F(2) = Vc(2)*Vp(2)+Vp(3)
    F(3) = Vp(2)*(Vc(3)+Vp(3))
  end function flux  

  subroutine evol1t1s
    integer :: i,j,k,l
    real(kind=8) :: ctemp
    real(kind=8),dimension(Nvar) :: V,Vnew
    real(kind=8),dimension(1:Nvar,1:Ncell) :: dV
    real(kind=8),dimension(1:Nvar,1:Ncell+1) :: F,diss
    ! vérification du pas de temps
    do i = 1,Ncell
       ctemp = sqrt((4.d0/3.d0)*(primVar(3,i)/primVar(1,i)))+abs(primVar(2,i))
       if (ctemp>c) c=ctemp
    end do
    if (dt > 0.4*dx/c) dt=0.1*dt
    if (dt < 0.05*dx/c) dt=10*dt
    
    ! calcul flux aux interfaces
    do k = 2,Ncell
       F(:,k) = 0.5d0*(flux(consVar(:,k-1))+flux(consVar(:,k)))
    end do
    F(:,1) = F(:,2)
    F(:,Ncell+1) = F(:,Ncell)
    !gridFlux=F

    ! calcul terme dissipatif
    do k = 2,Ncell+1
       diss(:,k) = (c+abs(primVar(2,k)))*0.5d0*(consVar(:,k-1)-consVar(:,k))
    end do
    diss(:,1) = diss(:,2)
    diss(:,Ncell+1) = diss(:,Ncell)

    ! calcul nouvelle valeur
    do i=1,Ncell
       do l=1,Nvar
          Vnew(l) = consVar(l,i) - (dt/dx)*(F(l,i+1)-F(l,i)) - (dt/dx)*(diss(l,i+1)-diss(l,i))
          consVar(l,i) = Vnew(l)
       end do
    end do
    t=t+dt
  end subroutine evol1t1s

  subroutine evol1t2s
   integer :: i,j,k,l
    real(kind=8) :: ctemp
    real(kind=8),dimension(Nvar) :: Vnew,Vm,Vp,Dm,Dp
    real(kind=8),dimension(1:Nvar,1:Ncell) :: dV,dVc
    real(kind=8),dimension(1:Nvar,1:Ncell+1) :: F,diss,Vl,Vr
    ! vérification du pas de temps
    do i = 1,Ncell
       ctemp = sqrt((4.d0/3.d0)*(primVar(3,i)/primVar(1,i)))+abs(primVar(2,i))
       if (ctemp>c) c=ctemp
    end do
    if (dt > 0.4*dx/c) dt=0.1*dt
    if (dt < 0.05*dx/c) dt=10*dt

    ! calcul limiteur minmod
    do i=1,Ncell
       if (i==1) then
          Vm = primVar(:,1)
          Vp = primVar(:,2)
       else if (i==Ncell) then
          Vm = primVar(:,Ncell-1)
          Vp = primVar(:,Ncell)
       else
          Vm = primVar(:,i-1)
          Vp = primVar(:,i+1)
       end if
       Dm = primvar(:,i) - Vm
       Dp = Vp - primVar(:,i)
       do j=1,Nvar
          if (Dm(j)>=0.d0 .and. Dp(j)>=0.d0) then
             dV(j,i) = min(Dm(j),Dp(j))
          else if (Dm(j)<=0.d0 .and. Dp(j)<=0.d0) then
             dV(j,i) = -min(-Dm(j),-Dp(j))
          else
             dV(i,j) = 0.d0
          end if
       end do
    end do
    do i=1,Ncell
       dVc(:,i) = prim2cons(dV(:,i))
    end do
    
    ! calcul valeur variables aux interfaces, k numéro interface entre i-1 et i
    do k=2,Ncell+1
       Vl(:,k) = consVar(:,k-1)+dVc(:,k-1)
    end do
    Vl(:,1) = Vl(:,2)
    do k=1,Ncell
       Vr(:,k) = consVar(:,k)-dVc(:,k)
    end do
    Vr(:,Ncell+1) = Vr(:,Ncell)

    ! calcul flux aux interfaces
    do k = 1,Ncell+1
       F(:,k) = 0.5d0*(flux(Vr(:,k))+flux(Vl(:,k)))
    end do
     
    ! calcul terme dissipatif
    do k = 1,Ncell
       diss(:,k) = (c+abs(primVar(2,k)))*0.5d0*(Vl(:,k)-Vr(:,k))
    end do
    diss(:,Ncell+1) = diss(:,Ncell)
     
    ! calcul nouvelle valeur
    do i=1,Ncell
       do l=1,Nvar
          Vnew(l) = consVar(l,i) - (dt/dx)*(F(l,i+1)-F(l,i)) - (dt/dx)*(diss(l,i+1)-diss(l,i))
          consVar(l,i) = Vnew(l)
       end do
    end do
    t=t+dt
  end subroutine evol1t2s

  subroutine evol2t1s
    integer :: i,j,k,l
    real(kind=8) :: ctemp
    real(kind=8),dimension(Nvar) :: V,Vnew,Dm,Dp,Vm,Vp
    real(kind=8),dimension(1:Nvar,1:Ncell) :: dV,Vh
    real(kind=8),dimension(1:Nvar,1:Ncell+1) :: F,diss
    ! vérification du pas de temps
    do i = 1,Ncell
       ctemp = sqrt((4.d0/3.d0)*(primVar(3,i)/primVar(1,i)))+abs(primVar(2,i))
       if (ctemp>c) c=ctemp
    end do
    if (dt > 0.4*dx/c) dt=0.1*dt
    if (dt < 0.05*dx/c) dt=10*dt
    
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
        Vh(:,i) = consVar(:,i) - dt/(2*dx)*(flux(consVar(:,i)+0.5*dV(:,i))-flux(consVar(:,i)-0.5*dV(:,i)))
     end do
     
    ! calcul flux aux interfaces
    do k = 2,Ncell
       F(:,k) = 0.5d0*(flux(Vh(:,k-1))+flux(Vh(:,k)))
    end do
    F(:,1) = F(:,2)
    F(:,Ncell+1) = F(:,Ncell)
    !gridFlux=F

    ! calcul terme dissipatif
    do l=1,Nvar
       do k = 2,Ncell+1
          diss(l,k) = (c+abs(primVar(2,k)))*0.5d0*(Vh(l,k-1)-Vh(l,k))
       end do
       diss(l,1) = diss(l,2)
       diss(l,Ncell+1) = diss(l,Ncell)
    end do
    
    ! calcul nouvelle valeur
    write(66,*) dt/dx
    do i=1,Ncell
       do l=1,Nvar
          Vnew(l) = consVar(l,i) - (dt/dx)*(F(l,i+1)-F(l,i)) - (dt/dx)*(diss(l,i+1)-diss(l,i))
          consVar(l,i) = Vnew(l)
       end do
    end do
    t=t+dt
  end subroutine evol2t1s
  
    
  subroutine evol2t2s
    integer :: i,j,k
    real(kind=8) :: ctemp
    real(kind=8),dimension(Nvar) :: Dm,Dp,V,Vm,Vp,Vnew
    real(kind=8),dimension(1:Nvar,1:Ncell) :: dV,dVc,Vh
    real(kind=8),dimension(1:Nvar,1:Ncell+1) :: F,Vl,Vr,diss
    
    ! vérification du pas de temps
    do i = 1,Ncell
       ctemp = sqrt((4.d0/3.d0)*(primVar(3,i)/primVar(1,i)))+abs(primVar(2,i))
       if (ctemp>c) c=ctemp
    end do
    if (dt > 0.4*dx/c) dt=0.1*dt
    if (dt < 0.05*dx/c) dt=10*dt
    
    ! calcul limiteur minmod
    do i=1,Ncell
       if (i==1) then
          Vm = primVar(:,1)
          Vp = primVar(:,2)
       else if (i==Ncell) then
          Vm = primVar(:,Ncell-1)
          Vp = primVar(:,Ncell)
       else
          Vm = primVar(:,i-1)
          Vp = primVar(:,i+1)
       end if
       Dm = primvar(:,i) - Vm
       Dp = Vp - primVar(:,i)
       do j=1,Nvar
          if (Dm(j)>=0.d0 .and. Dp(j)>=0.d0) then
             dV(j,i) = min(Dm(j),Dp(j))
          else if (Dm(j)<=0.d0 .and. Dp(j)<=0.d0) then
             dV(j,i) = -min(-Dm(j),-Dp(j))
          else
             dV(i,j) = 0.d0
          end if
       end do
    end do
    do i=1,Ncell
       dVc(:,i) = prim2cons(dV(:,i))
    end do
    
    ! calcul "demi pas de temps"
    do i = 1,Ncell
       Vh(:,i) = consVar(:,i) - dt/(2*dx)*(flux(consVar(:,i)+0.5*dVc(:,i))-flux(consVar(:,i)-0.5*dVc(:,i)))
    end do
       
    ! calcul valeur variables aux interfaces, k numéro interface entre i-1 et i
    do k=2,Ncell+1
       Vl(:,k) = Vh(:,k-1)+dVc(:,k-1)
    end do
    Vl(:,1) = Vl(:,2)
    do k=1,Ncell
       Vr(:,k) = Vh(:,k)-dVc(:,k)
    end do
    Vr(:,Ncell+1) = Vr(:,Ncell)

    ! calcul flux aux interfaces
    do k = 1,Ncell+1
       F(:,k) = 0.5d0*(flux(Vr(:,k))+flux(Vl(:,k)))
    end do

    ! calcul terme dissipatif
    do k = 1,Ncell
       diss(:,k) = (c+abs(primVar(2,k)))*0.5d0*(Vl(:,k)-Vr(:,k))
    end do
    diss(:,Ncell+1) = diss(:,Ncell)
       
    ! calcul nouvelle valeur
    do i=1,Ncell
       Vnew(:) = consVar(:,i) - (dt/dx)*(F(:,i+1)-F(:,i)) - (dt/dx)*(diss(:,i+1)-diss(:,i))
       consVar(:,i) = Vnew(:)
    end do
    t=t+dt
  end subroutine evol2t2s
end program simulation
