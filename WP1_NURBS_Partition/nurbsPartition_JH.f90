!=======================================================================
!File:        nurbsPartition_JH.f90
!Description: Parallel partitioning of NURBS knot vectors
!=======================================================================
!Author:      Jude Hussain
!Institution: Queen's University Belfast
!Created:     2026
!=======================================================================
!
!Public routines
!---------------
!partitionNurbs        — general d-dimensional partitioning.
!d=1 (curve), d=2 (surface), d=3 (volume) are special cases set by nDim.
!
!partitionNurbsCurve: wrapper for d=1 (nurbsCurve)
!partitionNurbsSurface: wrapper for d=2 (nurbsSurface)
!
!Private methods
!---------------
!partitionOneDirection — all span arithmetic for one parametric
!direction - dimension agnostic
!
!Control point storage convention (last index fastest, i.e., v-fastest for surfaces):
!d=1:  index = i1
!d=2:  index = i2 + (i1-1)*n2
!=======================================================================

module NURBS_PARTITION

use NURBS_CURVES
use NURBS_SURFACES

implicit none

private :: partitionOneDirection

contains

!=======================================================================
!private subroutine: partitionOneDirection
!
!A dimension-agnostic generalization  of partitionNurbsCurve - contains
!all span arithmetic for one direction 
!
!In:  U        — full knot vector
!nOfKnots — length of U
!q        — polynomial degree in this direction
!nP       — number of processors in this direction
!iProc    — this processor (1-indexed)
!
!Out: knotStart, knotEnd  — local knot index range in U
!ctrlStarts, ctrlEnd  — local control point index range
!=======================================================================
subroutine partitionOneDirection(knotStart, knotEnd, ctrlStart, ctrlEnd, U, nOfKnots, q, nP, iProc)
implicit none

double precision, intent(in) :: U(:)
integer, intent(in)  :: nOfKnots, q, nP, iProc
integer, intent(out) :: knotStart, knotEnd, ctrlStart, ctrlEnd

double precision, allocatable :: uUnique(:)
integer :: nSpans, spanStart, spanEnd, i

!Get breakpoints (unique knots) and corresponding nspans
uUnique = uniqueReal(U(1:nOfKnots))
nSpans  = size(uUnique) - 1

!Naive split of integer spans as per Swansea apporach
spanStart = (iProc - 1)*nSpans/nP + 1
spanEnd   =  iProc     *nSpans/nP

!Ghost layer: expand by q on each side
!(Cox-de Boor needs q spans of context on each side to evaluate Ni,p)
spanStart = max(1,      spanStart - q)
spanEnd   = min(nSpans, spanEnd   + q)

!Find knot indices for expanded span range
!Forward search for knotStart, backward search for knotEnd
knotStart = 1
knotEnd   = nOfKnots

!Finds first instance of the a knot value to maintain multiplity for Ni,p-1 evaluation
i = 1
do while (i <= nOfKnots)
    if (abs(U(i) - uUnique(spanStart)) < 1.0D-12) then
        knotStart = i
        exit
    end if
    i = i + 1
end do

!Finds last repeated instance of the a knot value to maintain multiplity for Ni+1,p-1 evaluation
i = nOfKnots
do while (i >= 1)
    if (abs(U(i) - uUnique(spanEnd + 1)) < 1.0D-12) then
        knotEnd = i
        exit
    end if
    i = i - 1
end do

!Obtain control point indices from knot indices
ctrlStart = knotStart
ctrlEnd   = knotEnd - q - 1

end subroutine partitionOneDirection


!=======================================================================
!public subroutine: partitionNurbs
!
!General d-dimensional NURBS partitioning for any parametric dimension
!
!Curves (d=1), surfaces (d=2) set by nDims. The code is extendable to
!volumes (d=3), which requires a nurbsVolume definition in occ_inteface.cpp.

!
!The wrappers below (partitionNurbsCurve, partitionNurbsSurface)
!unpack Swansea structs, call this routine, and repack
!Extendable to other dimensions by defining an additional wrapper 
!
!inputs:
!nDims                        —     number of parametric directions
!knots(maxKnots, nDims)       — knot vectors as columns; pad shorter ones
!with zeros (only nKnots(d) entries are read)
!nKnots(nDims)                —     actual length of each knot vector
!degrees(nDims)               —     polynomial degree per direction
!nProcs(nDims)                —     number of processors per direction
!iProcs(nDims)                —     this processor per direction (1-indexed)
!
!Out:
!knotStarts(nDims)           — first knot index of local patch per direction
!knotEnds(nDims)             — last  knot index of local patch per direction
!ctrlStarts(nDims)           — first ctrl pt index per direction
!ctrlEnds(nDims)             — last  ctrl pt index per direction
!=======================================================================
subroutine partitionNurbs(nDims, knots, nKnots, degrees, nProcs, iProcs, knotStarts, knotEnds, ctrlStarts, ctrlEnds)
implicit none

integer,          intent(in)  :: nDims
double precision, intent(in)  :: knots(:,:)
integer,          intent(in)  :: nKnots(nDims)
integer,          intent(in)  :: degrees(nDims)
integer,          intent(in)  :: nProcs(nDims)
integer,          intent(in)  :: iProcs(nDims)
integer,          intent(out) :: knotStarts(nDims), knotEnds(nDims)
integer,          intent(out) :: ctrlStarts(nDims), ctrlEnds(nDims)

integer :: d

!Partition each direction independently, i.e. as in nurbsCurve
!New wrapper must be written for higher nDims
do d = 1, nDims
    call partitionOneDirection( &
        knotStarts(d), knotEnds(d), ctrlStarts(d), ctrlEnds(d), &
        knots(:, d), nKnots(d), degrees(d), nProcs(d), iProcs(d))

    print *, "Dim", d, "| knots", knotStarts(d), "to", knotEnds(d), &
             "| ctrl",  ctrlStarts(d), "to", ctrlEnds(d)
end do

end subroutine partitionNurbs


!=======================================================================
!PUBLIC: partitionNurbsCurve  (d = 1)
!
!Wrapper using nurbsCurve - unpacks the struct, calls partitionNurbs, repacks result
!=======================================================================
subroutine partitionNurbsCurve(localCurve, globalCurve, nP, iProc)
implicit none

type(nurbsCurve), intent(out) :: localCurve
type(nurbsCurve), intent(in)  :: globalCurve
integer,          intent(in)  :: nP, iProc

integer          :: knotStarts(1), knotEnds(1), ctrlStarts(1), ctrlEnds(1)
integer          :: nLocalKnots, nLocalCtrl
double precision :: knots(globalCurve%nOfKnots, 1)


!Potential: add block to ensure that size(globalCurve%U) = globalCurve%nOfKnots to catch mismatches
if (size(globalCurve%U) /= globalCurve%nOfKnots) then
    error stop "partitionNurbsCurve: size(globalCurve%U) != globalCurve%nOfKnots"
end if

!Reshape for partitionNurbs: pack knot vector as a column.
!partitionNurbs is dimension-agnostic so it expects a 2D array of knots, even for d=1
knots(:, 1) = globalCurve%U

call partitionNurbs( &
    nDims      = 1,                        &
    knots      = knots,                    &
    nKnots     = [globalCurve%nOfKnots],   &  
    degrees    = [globalCurve%q],          &
    nProcs     = [nP],                     &
    iProcs     = [iProc],                  &
    knotStarts = knotStarts,               &
    knotEnds   = knotEnds,                 &
    ctrlStarts = ctrlStarts,               &
    ctrlEnds   = ctrlEnds)

nLocalKnots = knotEnds(1) - knotStarts(1) + 1
nLocalCtrl  = ctrlEnds(1) - ctrlStarts(1) + 1

localCurve%q                = globalCurve%q
localCurve%nOfKnots         = nLocalKnots
localCurve%nOfControlPoints = nLocalCtrl

!allocate space for copies of the local knot vector and control points
allocate(localCurve%U(nLocalKnots))
allocate(localCurve%Pw(4, nLocalCtrl))

!Copy global portion of knot vectors and control points into local curve
localCurve%U  = globalCurve%U(knotStarts(1):knotEnds(1))
localCurve%Pw = globalCurve%Pw(:, ctrlStarts(1):ctrlEnds(1))

end subroutine partitionNurbsCurve


!=======================================================================
!partitionNurbsSurface  (public wrapper, d = 2)
!
!Wrapper using nurbSurface - unpacks the struct, calls partitionNurbs, repacks result
!=======================================================================
subroutine partitionNurbsSurface(localSurface, globalSurface, nPu, nPv, iPu, iPv)
implicit none

type(nurbsSurface), intent(out) :: localSurface
type(nurbsSurface), intent(in)  :: globalSurface
integer,            intent(in)  :: nPu, nPv, iPu, iPv

integer          :: knotStarts(2), knotEnds(2), ctrlStarts(2), ctrlEnds(2)
integer          :: nLocalKnotsU, nLocalKnotsV, nLocalCtrlU, nLocalCtrlV
integer          :: maxKnots
integer          :: iuL, ivL, iuG, ivG, idxL, idxG
double precision, allocatable :: knots(:,:)

!Pack both knot vectors as columns, pad shorter one with zeros
!partitionOneDirection reads only nKnots(d) entries so padding is safe.
maxKnots = max(globalSurface%nOfKnotsU, globalSurface%nOfKnotsV)
allocate(knots(maxKnots, 2))
knots = 0.0D0
knots(1:globalSurface%nOfKnotsU, 1) = globalSurface%U
knots(1:globalSurface%nOfKnotsV, 2) = globalSurface%V

call partitionNurbs( &
    nDims      = 2,                                                  &
    knots      = knots,                                              &
    nKnots     = [globalSurface%nOfKnotsU, globalSurface%nOfKnotsV], &
    degrees    = [globalSurface%qU,        globalSurface%qV],        &
    nProcs     = [nPu, nPv],                                         &
    iProcs     = [iPu, iPv],                                         &
    knotStarts = knotStarts,                                         &
    knotEnds   = knotEnds,                                           &
    ctrlStarts = ctrlStarts,                                         &
    ctrlEnds   = ctrlEnds)

nLocalKnotsU = knotEnds(1) - knotStarts(1) + 1
nLocalKnotsV = knotEnds(2) - knotStarts(2) + 1
nLocalCtrlU  = ctrlEnds(1) - ctrlStarts(1) + 1
nLocalCtrlV  = ctrlEnds(2) - ctrlStarts(2) + 1

localSurface%qU                = globalSurface%qU
localSurface%qV                = globalSurface%qV
localSurface%nOfKnotsU         = nLocalKnotsU
localSurface%nOfKnotsV         = nLocalKnotsV
localSurface%nOfControlPointsU = nLocalCtrlU
localSurface%nOfControlPointsV = nLocalCtrlV
localSurface%nOfControlPoints  = nLocalCtrlU * nLocalCtrlV
localSurface%nOfPcurves        = 0

allocate(localSurface%U(nLocalKnotsU))
allocate(localSurface%V(nLocalKnotsV))
allocate(localSurface%Pw(4, nLocalCtrlU * nLocalCtrlV))

localSurface%U = globalSurface%U(knotStarts(1):knotEnds(1))
localSurface%V = globalSurface%V(knotStarts(2):knotEnds(2))

!Copy the rectangular block of control points.
!Pw is stored v-fastest: index = iv + (iu-1)*nOfControlPointsV
do iuL = 1, nLocalCtrlU
    iuG = ctrlStarts(1) + iuL - 1
    do ivL = 1, nLocalCtrlV
        ivG  = ctrlStarts(2) + ivL - 1
        idxG = ivG + (iuG - 1)*globalSurface%nOfControlPointsV
        idxL = ivL + (iuL - 1)*nLocalCtrlV
        localSurface%Pw(:, idxL) = globalSurface%Pw(:, idxG)
    end do
end do

end subroutine partitionNurbsSurface

!=======================================================================
end module NURBS_PARTITION