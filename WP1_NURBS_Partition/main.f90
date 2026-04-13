!=======================================================================
!  File:        main.f90
!  Description: Basic demo for reading a CAD and performing some 
!               geometric queries
!-----------------------------------------------------------------------
!  Author:      Rubén Sevilla and Xi Zou
!  Institution: Swansea University
!  Created:     2026
!=======================================================================

program main
use OCC_READER
use NURBS_PARTITION
implicit none

!> All surfaces
type(nurbsSurface), allocatable :: surfaces(:)

!> All intersection curves
type(nurbsCurve), allocatable :: curves(:)

!> Local curve for partitioning test
type(nurbsCurve) :: localCurve

!> Local surface for partitioning test
type(nurbsSurface) :: localSurface

double precision :: u, v, uv(2), C(3), dC(3), S(3), dSu(3), dSv(3), x(3)
double precision :: globalPt(3), localPt(3), ptWeight
integer :: iCurve, iSurface, iProc, nP
integer :: iPu, iPv, nPu, nPv


!=======================================================================
! Read and build data structure
!=======================================================================
call occ_read(curves, surfaces)

! Print curve 1 metadata
print *, "Curve 1 degree: ", curves(1)%q
print *, "Curve 1 nOfKnots: ", curves(1)%nOfKnots
print *, "Curve 1 knot vector: ", curves(1)%U
print *, "Curve 1 nOfControlPoints: ", curves(1)%nOfControlPoints

!=======================================================================
! Simulate the partitioning amongst 2 processors
!=======================================================================
nP = 2

do iProc = 1, nP
    print *, "===== Partitioning curve 1 for processor", iProc, "====="
    call partitionNurbsCurve(localCurve, curves(1), nP, iProc)
    
    print *, "Local knot vector: ", localCurve%U
    print *, "Num local control points: ", localCurve%nOfControlPoints
    
    u = 0.5D0*(localCurve%U(1) + localCurve%U(localCurve%nOfKnots))
    call nurbsCurvePointDerivative(C, dC, localCurve, u)
    print *, "Local eval at u=", u, ": ", C
    ! ========================================================================
    ! Check that the evaluation of the coordinates from the local and global curve return the same physical point

    !Declaring any variables that i need
    !in: u, out: 3D coordinates
    !evaluate the global curve at u
    call nurbsCurvePoint(globalPt, ptWeight, curves(1), u)
    !evaluate the local curve at the same u
    call nurbsCurvePoint(localPt, ptWeight, localCurve, u)

    print *, "Global eval at u=", u, ": ", globalPt
    print *, "Local eval at u=", u, ": ", localPt
    if (all(abs(globalPt - localPt) < 1.0D-12)) then
        print *, "Success: local and global evaluations match!"
    else
        print *, "Error: local and global evaluations do not match."
    end if
    deallocate(localCurve%U, localCurve%Pw)
end do

!=======================================================================
! Simulate the partitioning of surface 1 on a 2 x 2 processor grid
!=======================================================================
nPu = 2
nPv = 2

do iPu = 1, nPu
    do iPv = 1, nPv
        print *, "===== Partitioning surface 1 for processor (", iPu, ",", iPv, ") ====="
        call partitionNurbsSurface(localSurface, surfaces(1), nPu, nPv, iPu, iPv)

        print *, "Local U knot vector: ", localSurface%U
        print *, "Local V knot vector: ", localSurface%V
        print *, "Num local control points U: ", localSurface%nOfControlPointsU
        print *, "Num local control points V: ", localSurface%nOfControlPointsV

        ! Evaluate at the midpoint of the local parametric patch
        u = 0.5D0*(localSurface%U(1) + localSurface%U(localSurface%nOfKnotsU))
        v = 0.5D0*(localSurface%V(1) + localSurface%V(localSurface%nOfKnotsV))

        call nurbsSurfacePointDerivative(S, dSu, dSv, localSurface, u, v)
        print *, "Local surface eval at (u,v)=(", u, ",", v, "): ", S

        ! Check equivalence of local and global point evaluations
        call nurbsSurfacePoint(globalPt, ptWeight, surfaces(1), u, v)
        call nurbsSurfacePoint(localPt,  ptWeight, localSurface, u, v)

        print *, "Global surface eval at (u,v)=(", u, ",", v, "): ", globalPt
        print *, "Local  surface eval at (u,v)=(", u, ",", v, "): ", localPt

        if (all(abs(globalPt - localPt) < 1.0D-12)) then
            print *, "Success: local and global surface evaluations match!"
        else
            print *, "Error: local and global surface evaluations do not match."
        end if

        deallocate(localSurface%U, localSurface%V, localSurface%Pw)
    end do
end do


!=======================================================================
! Testing by calling some routines from the NURBS library
!=======================================================================
print *, "Example of forward operations for curves-----------------------------"
iCurve = 1
u = 0.3D0
call nurbsCurvePointDerivative(C, dC, curves(iCurve), u)
print *, "Curve Point:      ", C
print *, "Curve Derivative: ", dC

print *, "Example of forward operations for surfaces---------------------------"
iSurface = 1
u = 0.3D0
v = 20.0D0
call nurbsSurfacePointDerivative(S, dSu, dSv, surfaces(iSurface), u, v)
print *, "Surface Point:        ", S
print *, "Surface Derivative u: ", dSu
print *, "Surface Derivative v: ", dSv

print *, "Example of inverse operations for curves-----------------------------"
iCurve = 1
x = (/ -12277.9D0, -174.5D0, 2321.5D0 /)
u = nurbsCurvePointProjection(curves(iCurve), x)
print *, "Curve Parameter: ", u

print *, "Example of inverse operations for surfaces---------------------------"
iSurface = 1
x = (/ -12131.9D0, -195.2D0, 2454.6D0 /)
uv = nurbsSurfacePointProjection(surfaces(iSurface), x)
print *, "Surface Parameters: ", uv

end program main