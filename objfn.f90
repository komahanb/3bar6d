subroutine objfn
!Author : Komahan Boopathy, University of Dayton, OH, 11/22/2013
  use dimthree
  implicit none
  
  objF= x(1)*gamma*L(1) + x(2)*gamma*L(2) +  x(3)*gamma*L(3)

  return
end subroutine objfn

!++++++++++++++++++++++++++++++++++

subroutine consfn
!Author : Komahan Boopathy, University of Dayton, OH, 11/22/2013
  use dimthree
  implicit none

! INPUT:
! x : Design
! E : Young's modulus dat(4)
! L : Length dat(1:3)
! P : Load 
! theta: Load orientation
!
! OUTPUT:
! Output the value of all the constraints at the current design G(M)
!
! Author  : Komahan Boopathy 04/09/2013

  real*8::u,v,sigma(5)

  pu=P*cos(theta)
  pv=P*sin(theta)
 
  u=(L(2)/E)*(x(1)*pu + 2*dsqrt(2.0)*x(2)*pu + x(3)*pu + x(3)*pv - x(1)*pv)/(x(1)*x(2) + dsqrt(2.0)*x(1)*x(3) + x(2)*x(3))
  
  v=(L(2)/E)*(-x(1)*pu + x(3)*pu + x(1)*pv + x(3)*pv)/(x(1)*x(2) + dsqrt(2.0)*x(1)*x(3) + x(2)*x(3))

!!$if (loadcase.eq.2) then

    sigma(1)=-1.0*(sqrt(2.0)*x(2)*pu + x(3)*pu  +x(3)*pv)/(x(1)*x(2)+sqrt(2.0)*x(1)*x(3)+x(2)*x(3))

!!$else

!!$   sigma(1)=(sqrt(2.0)*x(2)*pu + x(3)*pu  +x(3)*pv)/(x(1)*x(2)+sqrt(2.0)*x(1)*x(3)+x(2)*x(3))
   
!!$ end if
 
 sigma(2)= (-(x(1)-x(3))*pu+(x(1)+x(3))*pv)/(x(1)*x(2)+sqrt(2.0)*x(1)*x(3)+x(2)*x(3))
 
 sigma(3)=(-sqrt(2.0)*x(2)*pu -x(1)*pu  +x(1)*pv)/(x(1)*x(2)+sqrt(2.0)*x(1)*x(3)+x(2)*x(3))
 

! stress constraints (normalized)

G(1) = (sigma(1) - 5000.0)/5000.0
G(4) = (-1.0*sigma(1)/5000) -1.0

G(2) = (sigma(2) - 20000.0)/20000.0
G(5) = (-1.0*sigma(2)/20000.0) -1.0

G(3) = (sigma(3) - 5000.0)/5000.0
G(6) = (-1.0*sigma(3) / 5000.0) -1.0
 
! displacement constraints (normalized)

G(7) = (-1.0*u -0.005)/0.005
G(8) = (v -0.005)/0.005

return

end subroutine consfn
!+++++++++++++++++++++++++++++++
subroutine loaddata
use dimthree
implicit none


  ! CONSTANT VALUES
  
  ! E=1.0e7 psi
  ! gamma=0.1 lb/in^3
  ! theta1=90
  ! theta2=135

  pi=4.0*atan(1.0) ! constant for later use (visible globally)

  !Length of the circular bars 1,2,3

  L(1) = 10.0*sqrt(2.0) !in 
  L(2) = 10.0           !in
  L(3) = 10.0*sqrt(2.0) !in
    
  E=1.0e7 !psi (Young's modulus)
  gamma=0.1   !lb/in3 (Weight density)
  
  ! 2 Loading cases

!!$  loadcase=1
!!$
!!$  if (loadcase.eq.1) then
!!$
!!$     theta = 90.0*pi/180.0 ! in radians
!!$     P=30000.0
!!$
!!$  else if (loadcase.eq.2) then

     theta = 135.0*pi/180.0 ! in radians
     P=20000.0

!!$  end if

  !Initial design (in)

  x(1)=1.0
  x(2)=1.0
  x(3)=1.0

end subroutine loaddata
!+++++++++++++++++++++++++++++++++++++++++++
subroutine gradobj
  use dimthree
  implicit none 
 
  grad(1) = L(1)*gamma
  grad(2) = L(2)*gamma
  grad(3) = L(3)*gamma
return
end subroutine gradobj
!++++++++++++++++++++++++++++++++++++++++++
subroutine gradconst
  use dimthree
  implicit none
  
  !  cgrad(M,N)
  
  cgrad(1,1)=((x(2) + 2.0**(1.0/2.0)*x(3))*(x(3)*pu + x(3)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(5000.0*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

  cgrad(1,2)=((x(1) + x(3))*(x(3)*pu + x(3)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(5000.0*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (2.0**(1.0/2.0)*pu)/(5000.0*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3)))

  cgrad(1,3)=((x(2) + 2.0**(1.0/2.0)*x(1))*(x(3)*pu + x(3)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(5000.0*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (pu + pv)/(5000.0*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3)))    

  cgrad(2,1)=((pu*(x(1) - x(3)) - pv*(x(1) + x(3)))*(x(2) + 2.0**(1.0/2.0)*x(3)))/(20000.0*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (pu - pv)/(20000.0*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3)))

  cgrad(2,2)=((x(1) + x(3))*(pu*(x(1) - x(3)) - pv*(x(1) + x(3))))/(20000.0*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

  cgrad(2,3)=(pu + pv)/(20000.0*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))) + ((pu*(x(1) - x(3)) - pv*(x(1) + x(3)))*(x(2) + 2.0**(1.0/2.0)*x(1)))/(20000.0*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

  cgrad(3,1)=((x(2) + 2.0**(1.0/2.0)*x(3))*(x(1)*pu - x(1)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(5000.0*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (pu - pv)/(5000.0*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3)))

  cgrad(3,2)=((x(1) + x(3))*(x(1)*pu - x(1)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(5000.0*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (2.0**(1.0/2.0)*pu)/(5000.0*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3)))

  cgrad(3,3)=((x(2) + 2.0**(1.0/2.0)*x(1))*(x(1)*pu - x(1)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(5000.0*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

  cgrad(4,1)=-((x(2) + 2.0**(1.0/2.0)*x(3))*(x(3)*pu + x(3)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(5000.0*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

  cgrad(4,2)=(2.0**(1.0/2.0)*pu)/(5000.0*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))) - ((x(1) + x(3))*(x(3)*pu + x(3)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(5000.0*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

  cgrad(4,3)=(pu + pv)/(5000.0*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))) - ((x(2) + 2.0**(1.0/2.0)*x(1))*(x(3)*pu + x(3)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(5000.0*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

  cgrad(5,1)=(pu - pv)/(20000.0*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))) - ((pu*(x(1) - x(3)) - pv*(x(1) + x(3)))*(x(2) + 2.0**(1.0/2.0)*x(3)))/(20000.0*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

  cgrad(5,2)=-((x(1) + x(3))*(pu*(x(1) - x(3)) - pv*(x(1) + x(3))))/(20000.0*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

  cgrad(5,3)=- (pu + pv)/(20000.0*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))) - ((pu*(x(1) - x(3)) - pv*(x(1) + x(3)))*(x(2) + 2.0**(1.0/2.0)*x(1)))/(20000.0*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

  cgrad(6,1)=(pu - pv)/(5000.0*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))) - ((x(2) + 2.0**(1.0/2.0)*x(3))*(x(1)*pu - x(1)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(5000.0*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

  cgrad(6,2)=(2.0**(1.0/2.0)*pu)/(5000.0*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))) - ((x(1) + x(3))*(x(1)*pu - x(1)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(5000.0*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

  cgrad(6,3)=-((x(2) + 2.0**(1.0/2.0)*x(1))*(x(1)*pu - x(1)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(5000.0*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

  !% Displacement
  cgrad(7,1) = (200.0*L(2)*(x(2) + 2.0**(1.0/2.0)*x(3))*(x(1)*pu + x(3)*pu - x(1)*pv + x(3)*pv + 2.0*2.0**(1.0/2.0)*x(2)*pu))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (200.0*L(2)*(pu - pv))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3)))

  cgrad(7,2) = (200.0*L(2)*(x(1) + x(3))*(x(1)*pu + x(3)*pu - x(1)*pv + x(3)*pv + 2.0*2.0**(1.0/2.0)*x(2)*pu))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (400.0*2.0**(1.0/2.0)*L(2)*pu)/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3)))

  cgrad(7,3) = (200.0*L(2)*(x(2) + 2.0**(1.0/2.0)*x(1))*(x(1)*pu + x(3)*pu - x(1)*pv + x(3)*pv + 2.0*2.0**(1.0/2.0)*x(2)*pu))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (200.0*L(2)*(pu + pv))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3)))

  cgrad(8,1)= -(200.0*L(2)*(pu - pv))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))) - (200.0*L(2)*(x(2) + 2.0**(1.0/2.0)*x(3))*(x(3)*pu - x(1)*pu + x(1)*pv + x(3)*pv))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

  cgrad(8,2)= -(200.0*L(2)*(x(1) + x(3))*(x(3)*pu - x(1)*pu + x(1)*pv + x(3)*pv))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

  cgrad(8,3)= (200.0*L(2)*(pu + pv))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))) - (200.0*L(2)*(x(2) + 2.0**(1.0/2.0)*x(1))*(x(3)*pu - x(1)*pu + x(1)*pv + x(3)*pv))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

  return  
end subroutine gradconst
