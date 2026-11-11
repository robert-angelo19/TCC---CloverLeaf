!Crown Copyright 2012 AWE.
!
! This file is part of CloverLeaf.
!
! CloverLeaf is free software: you can redistribute it and/or modify it under 
! the terms of the GNU General Public License as published by the 
! Free Software Foundation, either version 3 of the License, or (at your option) 
! any later version.
!
! CloverLeaf is distributed in the hope that it will be useful, but 
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
! details.
!
! You should have received a copy of the GNU General Public License along with 
! CloverLeaf. If not, see http://www.gnu.org/licenses/.

!>  @brief Fortran viscosity kernel.
!>  @author Wayne Gaudin
!>  @details Calculates an artificial viscosity using the Wilkin's method to
!>  smooth out shock front and prevent oscillations around discontinuities.
!>  Only cells in compression will have a non-zero value.

MODULE viscosity_kernel_module
  USE precision_mod, ONLY: wp, dp

CONTAINS

  SUBROUTINE viscosity_kernel(x_min,x_max,y_min,y_max,    &
                              celldx,celldy,              &
                              density0,                   &
                              pressure,                   &
                              viscosity,                  &
                              xvel0,                      &
                              yvel0                       )

    IMPLICIT NONE

    INTEGER     :: x_min,x_max,y_min,y_max
    REAL(KIND=wp), DIMENSION(x_min-2:x_max+2)                     :: celldx
    REAL(KIND=wp), DIMENSION(y_min-2:y_max+2)                     :: celldy
    REAL(KIND=wp), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2)     :: density0
    REAL(KIND=wp), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2)     :: pressure
    REAL(KIND=wp), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2)     :: viscosity
    REAL(KIND=wp), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3)     :: xvel0,yvel0

    INTEGER       :: j,k
    REAL(KIND=dp) :: ugrad,vgrad,grad2,pgradx,pgrady,pgradx2,pgrady2,grad     &
      ,ygrad,pgrad,xgrad,div,strain2,limiter,dirx,diry

!$ACC DATA &
!$ACC PRESENT(density0,pressure,viscosity,xvel0,yvel0)

!$ACC KERNELS

!$ACC LOOP INDEPENDENT
    DO k=y_min,y_max
!$ACC LOOP INDEPENDENT PRIVATE(ugrad,vgrad,div,strain2,pgradx,pgrady,pgradx2,pgrady2,limiter,pgrad,xgrad,ygrad,grad,grad2,dirx,diry)
      DO j=x_min,x_max
        ugrad=(REAL(xvel0(j+1,k  ),KIND=dp)+REAL(xvel0(j+1,k+1),KIND=dp)) &
          -(REAL(xvel0(j  ,k  ),KIND=dp)+REAL(xvel0(j  ,k+1),KIND=dp))

        vgrad=(REAL(yvel0(j  ,k+1),KIND=dp)+REAL(yvel0(j+1,k+1),KIND=dp)) &
          -(REAL(yvel0(j  ,k  ),KIND=dp)+REAL(yvel0(j+1,k  ),KIND=dp))

        div = REAL(celldx(j),KIND=dp)*ugrad+REAL(celldy(k),KIND=dp)*vgrad

        strain2 = 0.5_dp*(REAL(xvel0(j,  k+1),KIND=dp)+REAL(xvel0(j+1,k+1),KIND=dp) &
          -REAL(xvel0(j  ,k  ),KIND=dp)-REAL(xvel0(j+1,k  ),KIND=dp))/REAL(celldy(k),KIND=dp) &
          +0.5_dp*(REAL(yvel0(j+1,k  ),KIND=dp)+REAL(yvel0(j+1,k+1),KIND=dp) &
          -REAL(yvel0(j  ,k  ),KIND=dp)-REAL(yvel0(j  ,k+1),KIND=dp))/REAL(celldx(j),KIND=dp)

        pgradx=(REAL(pressure(j+1,k),KIND=dp)-REAL(pressure(j-1,k),KIND=dp)) &
          /(REAL(celldx(j),KIND=dp)+REAL(celldx(j+1),KIND=dp))
        pgrady=(REAL(pressure(j,k+1),KIND=dp)-REAL(pressure(j,k-1),KIND=dp)) &
          /(REAL(celldy(k),KIND=dp)+REAL(celldy(k+1),KIND=dp))

        pgradx2 = pgradx*pgradx
        pgrady2 = pgrady*pgrady

        limiter = ((0.5_dp*ugrad/REAL(celldx(j),KIND=dp))*pgradx2 &
          +(0.5_dp*vgrad/REAL(celldy(k),KIND=dp))*pgrady2+strain2*pgradx*pgrady) &
          /MAX(pgradx2+pgrady2,1.0e-16_dp)

        IF ((limiter.GT.0.0_dp).OR.(div.GE.0.0_dp))THEN
          viscosity(j,k) = 0.0_wp
        ELSE
          dirx=1.0_dp
          IF(pgradx.LT.0.0_dp) dirx=-1.0_dp
          pgradx = dirx*MAX(1.0e-16_dp,ABS(pgradx))
          diry=1.0_dp
          IF(pgradx.LT.0.0_dp) diry=-1.0_dp
          pgrady = diry*MAX(1.0e-16_dp,ABS(pgrady))
          pgrad = SQRT(pgradx**2+pgrady**2)
          xgrad = ABS(REAL(celldx(j),KIND=dp)*pgrad/pgradx)
          ygrad = ABS(REAL(celldy(k),KIND=dp)*pgrad/pgrady)
          grad  = MIN(xgrad,ygrad)
          grad2 = grad*grad

          viscosity(j,k)=REAL(2.0_dp*REAL(density0(j,k),KIND=dp)*grad2*limiter*limiter,KIND=wp)
        ENDIF

      ENDDO
    ENDDO


!$ACC END KERNELS

!$ACC END DATA

  END SUBROUTINE viscosity_kernel

END MODULE viscosity_kernel_module
