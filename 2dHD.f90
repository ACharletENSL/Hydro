module grid
  implicit none
  ! creates grid
  integer,parameter :: Ncell=15,Ndim=2,Nvar=4,Nsteps=20
  real(kind=8),parameter :: ds=0.01d0,dx=ds,dy=ds,dt_i=0.0001d0,T_tot=Nsteps*dt_i
  real(kind=8),save :: c,dt,t
  real(kind=8),dimension(1:Nvar,1:Ncell,1:Ncell),save :: primVar,consVar
  real(kind=8),dimension(1:Nvar,1:Ncell+1,1:Ncell),save :: gridFx,gridDx
  real(kind=8),dimension(1:Nvar,1:Ncell,1:Ncell+1),save :: gridFy,gridDy
  character(len=16),save :: filename1,filename2,filename3
  character(len=128),save :: path='/media/acharlet/Data/Arthur/Documents/Cours/4A_CRAL/Code/Perso/Hydro/Resultats2D/',FileName
  
  contains
    subroutine init()
      integer :: i,j
      c=1; dt=dt_i; t=0
      !density
      primVar(1,:,:) = 1.0d0
      primVar(1,6:10,6:10)=3.d0;
      !velocity
      primVar(2,:,:) = 0.d0
      primVar(3,:,:) = 0.d0
      !pressure
      primVar(4,:,:) = 0.1d0*primVar(1,:,:)    
      consVar(1,:,:) = primVar(1,:,:)   ! D = rho
      consVar(2,:,:) = 0.d0          ! S = rho*v , v=0
      consVar(3,:,:) = 0.d0
      consVar(4,:,:) = primVar(1,:,:)-primVar(4,:,:) ! U = rho-p

      write (filename1, "('data0000.dat')")
      FileName=trim(adjustl(path))//trim(adjustl(filename1))
      open(10,file=FileName,status='new')
101   format(2(1x,i5),2x,4e20.10)
      do i=1,Ncell
         do j=1,Ncell
            write(10,101) i, j, primVar(1,i,j), primVar(2,i,j), primVar(3,i,j), primVar(4,i,j)
         end do     
      end do
      close(10)
    end subroutine init

    ! variables primitives vers variables conservatives
  function cons2prim(Vc) result(Vp)
    real(kind=8),dimension(Nvar),intent(in) :: Vc
    real(kind=8),dimension(Nvar) :: Vp
    Vp(1) = Vc(1)
    Vp(2) = Vc(2)/Vc(1)
    Vp(3) = Vc(3)/Vc(1)
    Vp(4) = Vc(1)-Vc(4)
  end function cons2prim
  
  ! variables conservatives vers variables primitives
  function prim2cons(Vp) result(Vc)
    real(kind=8),dimension(Nvar),intent(in) :: Vp
    real(kind=8),dimension(Nvar) :: Vc
    Vc(1) = Vp(1)
    Vc(2) = Vp(1)*Vp(2)
    Vc(3) = Vp(1)*Vp(3)
    Vc(4) = Vp(1)-Vp(4)
  end function prim2cons
end module grid

program simulation
  use grid
  implicit none
  integer :: iStep,i,j
  
  call init()
  do iStep=1,Nsteps
     write (filename1, "('data',I4.4,'.dat')") iStep
     FileName=trim(adjustl(path))//trim(adjustl(filename1))
     open(10,file=FileName,status='new')
     call evol1t1s
     do i=1,Ncell
        do j=1,Ncell
           primVar(:,i,j) = cons2prim(consVar(:,i,j))
        end do
     end do
101  format(2(1x,i5),2x,4e20.10)
     do i=1,Ncell
        do j=1,Ncell
            write(10,101) i, j, primVar(1,i,j), primVar(2,i,j), primVar(3,i,j), primVar(4,i,j)
        end do
     end do
     !write (filename2, "('flux',I4.4,'.dat')") iStep
     !FileName=trim(adjustl(path))//trim(adjustl(filename2))
     !open(11,file=FileName,status='new')
     !do i=1,Ncell
     !   do j=1,Ncell
     !       write(11,101) i, j, gridFx(1,i,j), gridFy(1,i,j), gridFx(4,i,j), gridFy(4,i,j)
     !   end do
     !end do
     !write (filename3, "('diss',I4.4,'.dat')") iStep
     !FileName=trim(adjustl(path))//trim(adjustl(filename3))
     !open(12,file=FileName,status='new')
     !do i=1,Ncell
     !   do j=1,Ncell
     !       write(12,101) i, j, gridDx(1,i,j), gridDy(1,i,j), gridDx(4,i,j), gridDy(4,i,j)
     !   end do
     !end do
  end do
     
contains
  function fluxx(Vc) result(Fx)
    ! flux > 0 quand va vers i croissants
    real(kind=8),dimension(Nvar),intent(in) :: Vc
    real(kind=8),dimension(Nvar) :: Vp
    real(kind=8),dimension(Nvar) :: Fx
    Vp = cons2prim(Vc)
    Fx(1) = Vc(2)
    Fx(2) = Vc(2)*Vp(2)+Vp(4) !flux Sx selon x
    Fx(3) = Vc(3)*Vp(2)+Vp(4) !flux Sy selon x
    Fx(4) = Vc(2)
  end function fluxx
  function fluxy(Vc) result(Fy)
    ! flux > 0 quand va vers i croissants
    real(kind=8),dimension(Nvar),intent(in) :: Vc
    real(kind=8),dimension(Nvar) :: Vp
    real(kind=8),dimension(Nvar) :: Fy
    Vp = cons2prim(Vc)
    Fy(1) = Vc(3)
    Fy(2) = Vc(2)*Vp(3)+Vp(4) !flux Sx selon y
    Fy(3) = Vc(3)*Vp(3)+Vp(4) !flux Sy selon y
    Fy(4) = Vc(3)
  end function fluxy

  subroutine evol1t1s
    integer :: i,j,k,l
    real(kind=8) :: ctemp
    real(kind=8),dimension(Nvar) :: Vnew
    real(kind=8),dimension(1:Nvar,1:Ncell+1,1:Ncell) :: Fx,dissx
    real(kind=8),dimension(1:Nvar,1:Ncell,1:Ncell+1) :: Fy,dissy
    ! vérification pas de temps
    do i = 1,Ncell
      do  j = 1,Ncell
         ctemp = sqrt((4.d0/3.d0)*(primVar(4,i,j)/primVar(1,i,j)))+sqrt(primVar(2,i,j)**2+primVar(3,i,j))
         if (ctemp>c) c=ctemp
      end do
    end do
    if (dt > 0.4*dx/c) dt=0.1d0*dt
    if (dt < 0.05*dx/c) dt=10d0*dt

    ! calcul flux aux interfaces
    do l = 1,Ncell
       do k = 2,Ncell
          Fx(:,k,l) = 0.5d0*(fluxx(consVar(:,k-1,l))+fluxx(consVar(:,k,l)))
       end do
       Fx(:,1,l) = Fx(:,2,l)
       Fx(:,Ncell+1,l) = Fx(:,Ncell,l)
    end do
    gridFx(:,:,:) = Fx(:,:,:)
    do k = 1,Ncell
       do l = 2,Ncell
          Fy(:,k,l) = 0.5d0*(fluxy(consVar(:,k,l-1))+fluxy(consVar(:,k,l)))
       end do
       Fy(:,k,1) = Fy(:,k,2)
       Fy(:,k,Ncell+1) = Fy(:,k,Ncell)
    end do
    gridFy(:,:,:)=Fy(:,:,:)

    ! calcul terme dissipatif
    do l = 1,Ncell
       do k = 2,Ncell
          dissx(:,k,l) = c*0.5d0*(consVar(:,k-1,l)-consVar(:,k,l))
       end do
       dissx(:,1,l) = dissx(:,2,l)
       dissx(:,Ncell+1,l) = dissx(:,Ncell,l)
    end do
    gridDx(:,:,:)=dissx(:,:,:)
    do k = 1,Ncell
       do l = 2,Ncell
          dissy(:,k,l) = c*0.5d0*(consVar(:,k,l-1)-consVar(:,k,l))
       end do
       dissy(:,k,1) = dissy(:,k,2)
       dissy(:,k,Ncell+1) = dissy(:,k,Ncell)
    end do
    gridDy(:,:,:)=dissy(:,:,:)
    ! calcul nouvelle valeur
    do i=1,Ncell
       do j=1,Ncell
          Vnew(:) = consVar(:,i,j) -(dt/dx)*(Fx(:,i+1,j)-Fx(:,i,j)) - (dt/dx)*(dissx(:,i+1,j)-dissx(:,i,j))&
               - (dt/dy)*(Fy(:,i,j+1)-Fy(:,i,j)) - (dt/dy)*(dissy(:,i,j+1)-dissy(:,i,j))
          consVar(:,i,j) = Vnew(:)
       end do
    end do
    t=t+dt
  end subroutine evol1t1s
end program simulation
