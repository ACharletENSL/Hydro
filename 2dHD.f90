module grid
  implicit none
  ! creates grid
  integer,parameter :: Ncell=10,Ndim=2,Nvar=4,Nsteps=10
  real(kind=8),parameter :: dx=0.01d0,dt_i=0.001d0,T_tot=Nsteps*dt_i
  real(kind=8),save :: c,dt,t
  real(kind=8),dimension(1:Nvar,1:Ncell,1:Ncell),save :: primVar,consVar
   real(kind=8),dimension(1:Nvar,1:Ncell+1,1:Ncell+1),save :: gridFlux,gridDiss,consR,consL
  character(len=16),save :: filename1,filename2,filename3
  character(len=128),save :: path='/media/acharlet/Data/Arthur/Documents/Cours/4A_CRAL/Code/Perso/Hydro/Resultats2D/',FileName
  
  contains
    subroutine init()
      integer :: i,j
      c=1; dt=dt_i; t=0
      primVar(1,1:floor(Ncell/2.),1:floor(Ncell/2.))=3.d0;
      primVar(1,floor(Ncell/2.):Ncell,floor(Ncell/2.):Ncell)=1.d0 ! initial densities
      primVar(2,:,:) = 0.d0  ! initial velocities
      primVar(3,:,:) = 0.d0
      primVar(4,:,:) = primVar(1,:) ! initial pressure     
      consVar(1,:,:) = primVar(1,:)   ! D = rho
      consVar(2,:,:) = 0.d0          ! S = rho*h*v , v=0
      consVar(3,:,:) = 0.d0
      consVar(4,:,:) = 4.d0*primVar(1,:) ! U = rho*h-p , h=5, rho=p

      write (filename1, "('data0000.dat')")
      FileName=trim(adjustl(path))//trim(adjustl(filename1))
      open(10,file=FileName,status='new')
100   format (3e20.10)
      do i=1,Ncell
         do j=1,Ncell
            !write(10,100) consVar(1,i,j), consVar(2,i,j), consVar(3,i,j), consVar(4,i,j)
            write(10,100) primVar(1,i,j), primVar(2,i,j), primVar(3,i,j), primVar(4,i)
         end do     
      end do
      close(10)
    end subroutine init

    ! variables primitives vers variables conservatives
  function cons2prim(Vc) result(Vp)
    real(kind=8),dimension(Nvar),intent(in) :: Vc
    real(kind=8),dimension(Nvar) :: Vp
    Vp(1) = Vc(1)
    Vp(4) = (Vc(4)-Vc(1))/3.d0
    Vp(2) = Vc(2)/(Vc(1)+4.d0*Vp(4))
    Vp(3) = Vc(3)/(Vc(1)+4.d0*Vp(4))
  end function cons2prim
  
  ! variables conservatives vers variables primitives
  function prim2cons(Vp) result(Vc)
    real(kind=8),dimension(Nvar),intent(in) :: Vp
    real(kind=8),dimension(Nvar) :: Vc
    Vc(1) = Vp(1)
    Vc(2) = Vp(2)*(Vp(1)+4.d0*Vp(4))
    Vc(3) = Vp(3)*(Vp(1)+4.d0*Vp(4))
    Vc(4) = Vp(1)+3.d0*Vp(4)
  end function prim2cons
end module grid

program simulation
  use grid
  implicit none
  integer :: iStep,i,j
  
  call init()
  do Step=1,Nsteps
     write (filename1, "('data',I4.4,'.dat')") iStep
     FileName=trim(adjustl(path))//trim(adjustl(filename1))
     open(10,file=FileName,status='new')
     call evol
     do i=1,Ncell
        do j=1,Ncell
           primVar(:,i,j) = cons2prim(consVar(:,i,j))
        end do
     end do
     
100  format (3e20.10e3)
     do i=1,Ncell
        do j=1,Ncell
           !write(10,100) consVar(1,i,j), consVar(2,i,j), consVar(3,i,j), consVar(4,i,j)
           write(10,100) primVar(1,i,j), primVar(2,i,j), primVar(3,i,j), primVar(4,i)
        end do
     end do
  end do
     
contains
  function flux(Vc) result(F)
    ! flux > 0 quand va vers i croissants
    real(kind=8),dimension(Nvar),intent(in) :: Vc
    real(kind=8),dimension(Nvar) ::Vp,F
    Vp = cons2prim(Vc)
    F(1) = Vc(1)*Vp(2)
    F(2) = Vc(2)*Vp(2)+Vp(4)
    F(3) = Vc(3)*Vp(3)+Vp(4)
    F(4) = Vp(2)*(Vc(3)+Vp(4))
  end function flux
        
