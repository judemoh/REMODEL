!=======================================================================
!  File:        nurbsPartition_JH.f90
!  Description: Parallel partitioning of NURBS knot vectors
!=======================================================================
!  Author:      Jude Hussain
!  Institution: Queen's University Belfast
!  Created:     2026
!=======================================================================
module NURBS_PARTITION

use NURBS_CURVES
use NURBS_SURFACES

contains

!========================================================================
! Partition a NURBS curve knot vector across nP processors
! Returns the local curve for processor iProc (1-indexed)
!========================================================================

subroutine partitionNurbsCurve(localCurve, globalCurve, nP, iProc)
implicit none
!intent are like python type hints but enforced, inent(in) means input, intent(out) means output, intent(inout) means both
type(nurbsCurve), intent(out) :: localCurve
type(nurbsCurve), intent(in)  :: globalCurve
integer, intent(in) :: nP      ! number of processors
integer, intent(in) :: iProc   ! this processor (1 to nP)

double precision, allocatable :: uUnique(:)
integer :: nSpans, spanStart, spanEnd
integer :: knotStart, knotEnd
integer :: ctrlStart, ctrlEnd
integer :: nLocalKnots, nLocalCtrl
integer :: i, q

q = globalCurve%q

! Get unique knot values (the breakpoints)
uUnique = uniqueReal(globalCurve%U)
nSpans = size(uUnique) - 1

! --- Step 1: naive split of knot spans across processors ---
! divide spans as evenly as possible
! this is all in integers - so n_proc isnt even, then load is inherently unbalanced
!processor count starts from 0 in C but from 1 in fortran, so we use iProc-1 in the calculations
spanStart = (iProc-1)*nSpans/nP + 1
spanEnd   = iProc*nSpans/nP

print *, "Proc", iProc, "owns knot spans", spanStart, "to", spanEnd

! --- Step 2: expand by q spans on each side ---
! this is the ghost layer at the NURBS level
spanStart = max(1, spanStart - q)
spanEnd   = min(nSpans, spanEnd + q)

print *, "Proc", iProc, "expanded to spans", spanStart, "to", spanEnd

! --- Step 3: find the corresponding knot indices ---
! in the full knot vector (with repetitions)
! knotStart is first knot of spanStart
! knotEnd is last knot of spanEnd+1
knotStart = 1
knotEnd   = globalCurve%nOfKnots

! find where spanStart begins in full knot vector
i = 1
do while(i <= globalCurve%nOfKnots)
    if(abs(globalCurve%U(i) - uUnique(spanStart)) < 1.0D-12) then
        knotStart = i
        exit
    end if
    i = i + 1
end do

! find where spanEnd ends in full knot vector
i = globalCurve%nOfKnots
do while(i >= 1)
    if(abs(globalCurve%U(i) - uUnique(spanEnd+1)) < 1.0D-12) then
        knotEnd = i
        exit
    end if
    i = i - 1
end do

! --- Step 4: select control points ---
! control point index = knot index - q
ctrlStart = knotStart
ctrlEnd   = knotEnd - q - 1

nLocalKnots = knotEnd - knotStart + 1
nLocalCtrl  = ctrlEnd - ctrlStart + 1

print *, "Proc", iProc, "knot indices", knotStart, "to", knotEnd
print *, "Proc", iProc, "ctrl indices", ctrlStart, "to", ctrlEnd

! --- Step 5: build local curve ---
localCurve%q = q
localCurve%nOfKnots = nLocalKnots
localCurve%nOfControlPoints = nLocalCtrl

allocate(localCurve%U(nLocalKnots))
!4 x nLocalCtrl array - actual coordinates are recovered via dividing by w component
allocate(localCurve%Pw(4, nLocalCtrl))

!this copies the local portion of the knot vector
localCurve%U  = globalCurve%U(knotStart:knotEnd)
!copies the corresponding control points including the weights
localCurve%Pw = globalCurve%Pw(:, ctrlStart:ctrlEnd)

end subroutine partitionNurbsCurve

!========================================================================
! Partition a NURBS surface knot vectors across an nPu x nPv processor grid
! Returns the local surface for processor (iPu, iPv), both 1-indexed
!========================================================================

subroutine partitionNurbsSurface(localSurface, globalSurface, nPu, nPv, iPu, iPv)
implicit none

type(nurbsSurface), intent(out) :: localSurface
type(nurbsSurface), intent(in)  :: globalSurface
integer, intent(in) :: nPu      ! number of partitions in u direction
integer, intent(in) :: nPv      ! number of partitions in v direction
integer, intent(in) :: iPu      ! this partition in u (1 to nPu)
integer, intent(in) :: iPv      ! this partition in v (1 to nPv)

double precision, allocatable :: uUnique(:), vUnique(:)
integer :: nSpansU, nSpansV
integer :: ownedStartU, ownedEndU, ownedStartV, ownedEndV
integer :: spanStartU, spanEndU, spanStartV, spanEndV
integer :: knotStartU, knotEndU, knotStartV, knotEndV
integer :: ctrlStartU, ctrlEndU, ctrlStartV, ctrlEndV
integer :: nLocalKnotsU, nLocalKnotsV
integer :: nLocalCtrlU, nLocalCtrlV, nLocalCtrl
integer :: i, qU, qV
integer :: iuL, ivL, iuG, ivG, idxL, idxG

qU = globalSurface%qU
qV = globalSurface%qV

! Get unique knot values (the breakpoints) in each direction
uUnique = uniqueReal(globalSurface%U)
vUnique = uniqueReal(globalSurface%V)

nSpansU = size(uUnique) - 1
nSpansV = size(vUnique) - 1

! --- Step 1: naive split of knot spans across processor grid ---
ownedStartU = (iPu-1)*nSpansU/nPu + 1
ownedEndU   = iPu*nSpansU/nPu

ownedStartV = (iPv-1)*nSpansV/nPv + 1
ownedEndV   = iPv*nSpansV/nPv

print *, "Proc (", iPu, ",", iPv, ") owns U spans", ownedStartU, "to", ownedEndU
print *, "Proc (", iPu, ",", iPv, ") owns V spans", ownedStartV, "to", ownedEndV

! --- Step 2: expand by qU and qV spans on each side ---
! this is the ghost layer at the NURBS level
spanStartU = max(1, ownedStartU - qU)
spanEndU = min(nSpansU, ownedEndU + qU)

spanStartV = max(1, ownedStartV - qV)
spanEndV   = min(nSpansV, ownedEndV + qV)

print *, "Proc (", iPu, ",", iPv, ") expanded U spans", spanStartU, "to", spanEndU
print *, "Proc (", iPu, ",", iPv, ") expanded V spans", spanStartV, "to", spanEndV

! --- Step 3: find corresponding knot indices in the full knot vectors ---
! knotStartU is first knot of spanStartU
! knotEndU   is last knot of spanEndU+1
knotStartU = 1
knotEndU   = globalSurface%nOfKnotsU

i = 1
do while(i <= globalSurface%nOfKnotsU)
    if(abs(globalSurface%U(i) - uUnique(spanStartU)) < 1.0D-12) then
        knotStartU = i
        exit
    end if
    i = i + 1
end do

i = globalSurface%nOfKnotsU
do while(i >= 1)
    if(abs(globalSurface%U(i) - uUnique(spanEndU+1)) < 1.0D-12) then
        knotEndU = i
        exit
    end if
    i = i - 1
end do

! knotStartV is first knot of spanStartV
! knotEndV   is last knot of spanEndV+1
knotStartV = 1
knotEndV   = globalSurface%nOfKnotsV

i = 1
do while(i <= globalSurface%nOfKnotsV)
    if(abs(globalSurface%V(i) - vUnique(spanStartV)) < 1.0D-12) then
        knotStartV = i
        exit
    end if
    i = i + 1
end do

i = globalSurface%nOfKnotsV
do while(i >= 1)
    if(abs(globalSurface%V(i) - vUnique(spanEndV+1)) < 1.0D-12) then
        knotEndV = i
        exit
    end if
    i = i - 1
end do

! --- Step 4: select control points in both directions ---
! same idea as for curves, but now we build a rectangular block
ctrlStartU = knotStartU
ctrlEndU   = knotEndU - qU - 1

ctrlStartV = knotStartV
ctrlEndV   = knotEndV - qV - 1

nLocalKnotsU = knotEndU - knotStartU + 1
nLocalKnotsV = knotEndV - knotStartV + 1

nLocalCtrlU = ctrlEndU - ctrlStartU + 1
nLocalCtrlV = ctrlEndV - ctrlStartV + 1
nLocalCtrl  = nLocalCtrlU*nLocalCtrlV

! print *, "Proc (", iPu, ",", iPv, ") U knot indices", knotStartU, "to", knotEndU
! print *, "Proc (", iPu, ",", iPv, ") V knot indices", knotStartV, "to", knotEndV
! print *, "Proc (", iPu, ",", iPv, ") U ctrl indices", ctrlStartU, "to", ctrlEndU
! print *, "Proc (", iPu, ",", iPv, ") V ctrl indices", ctrlStartV, "to", ctrlEndV

! --- Step 5: build local surface ---
localSurface%qU = qU
localSurface%qV = qV
localSurface%nOfKnotsU = nLocalKnotsU
localSurface%nOfKnotsV = nLocalKnotsV
localSurface%nOfControlPointsU = nLocalCtrlU
localSurface%nOfControlPointsV = nLocalCtrlV
localSurface%nOfControlPoints  = nLocalCtrl

allocate(localSurface%U(nLocalKnotsU))
allocate(localSurface%V(nLocalKnotsV))
allocate(localSurface%Pw(4, nLocalCtrl))

localSurface%U = globalSurface%U(knotStartU:knotEndU)
localSurface%V = globalSurface%V(knotStartV:knotEndV)

! trimming p-curves are not partitioned here
localSurface%nOfPcurves = 0


! --- Step 6: copy the rectangular block of control points ---
! Pw is stored in flattened form with v varying fastest:
! index = iv + (iu-1)*nOfControlPointsV
do iuL = 1, nLocalCtrlU
    iuG = ctrlStartU + iuL - 1
    do ivL = 1, nLocalCtrlV
        ivG = ctrlStartV + ivL - 1

        idxG = ivG + (iuG - 1)*globalSurface%nOfControlPointsV
        idxL = ivL + (iuL - 1)*nLocalCtrlV

        localSurface%Pw(:, idxL) = globalSurface%Pw(:, idxG)
    end do
end do


end subroutine partitionNurbsSurface

!========================================================================

end module NURBS_PARTITION