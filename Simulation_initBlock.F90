!!****if* source/Simulation/SimulationMain/Cellular/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer, INTENT(in)::blockID)
!!
!!
!! DESCRIPTION
!!
!!   Initialize solution data in one block for a planar detonation,
!!   perturbed with some noise.  This detonation should be unstable
!!   to becoming cellular at the detonation front.
!!
!! ARGUMENTS
!!
!!  blockID:        integer  the number of the block to initialize
!!
!!
!! PARAMETERS
!!
!!    xhe4               mass fraction of he4
!!    xc12               mass fraction of c12
!!    xo16               mass fraction of o16
!!    rhoAmbient         density of the cold upstream material 
!!    tempAmbient        temperature of the cold upstream material
!!    velxAmbient        x-velocity of the cold upstream material
!!    rhoPerturb         density of the post shock material
!!    tempPerturb        temperature of the post shock material
!!    velxPerturb        x-velocity of the post shock material
!!    radiusPerturb      distance below which the perturbation is applied
!!    xCenterPerturb     origin of the of the perturbation
!!    yCenterPerturb     origin of the of the perturbation
!!    zCenterPerturb     origin of the of the perturbation
!!    usePseudo1d        .true. for a 1d initial configuration, with the ??
!!                          copied along the y and z directions
!!                       .false. for a spherical configuration
!!    noiseAmplitude     amplitude of the white noise added to the perturbation
!!    noiseDistance      distances above and below radiusPerturb get noise added
!!    xmax               boundary of domain
!!    xmin               boundary of domain
!!    ymax               boundary of domain
!!    ymin               boundary of domain
!!    zmax               boundary of domain
!!    zmin               boundary of domain
!!
!!    smallx             smallest allowed abundance
!!    smlrho             smallest allowed density
!!
!!  NOTES
!!    Species used in this simulation are HE4 (helium-4), C12 (carbon-12), O16 (oxygen-16)
!!
!! NOTES
!!   See paper: Timmes, FX; Zingale, M; Olson, K; Fryxell, B; The Astrophysical
!!               Journal, Nov. 10, 2000 : 543: 938-954
!!  
!!***

subroutine Simulation_initBlock(blockID)
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords, Grid_putPointData
  use User_interface, ONLY : User_initBlobCell
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"

  integer,INTENT(in) ::  blockID

!  Local variables

  real :: xx, yy, zz
  logical, parameter :: useGuardCell = .TRUE.

  integer, dimension(2,MDIM), save :: blockRange, blockExtent

  real,allocatable,dimension(:) :: xCoordsCell, yCoordsCell, zCoordsCell
  integer,dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: sizeX, sizeY, sizeZ

  integer :: i, j, k, n
  integer, dimension(MDIM) :: iPosition   !for putting data with Grid_putData

  ! variables needed for the eos call
  real, dimension(NUNK_VARS) :: unkArr

! ----------------------------------------------------------------------------------------------

 ! determine size of blocks
  call Grid_getBlkIndexLimits(blockID,blockRange,blockExtent)

  ! Get the indices of the blocks
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
  allocate(xCoordsCell(sizeX))
  sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
  allocate(yCoordsCell(sizeY))
  sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
  allocate(zCoordsCell(sizeZ))

!  if (NDIM == 3)  &
   call Grid_getCellCoords(KAXIS, blockID, CENTER, useGuardCell, zCoordsCell, sizeZ)
!  if (NDIM >= 2)  &
   call Grid_getCellCoords(JAXIS, blockID, CENTER, useGuardCell, yCoordsCell, sizeY)
   call Grid_getCellCoords(IAXIS, blockID, CENTER, useGuardCell, xCoordsCell, sizeX)

  !! NOTE the incorrect syntax below --  this causes a crash for KAXIS when NDIM < 3, 
  !!    works OK if you substitute 1 at the MAXCELLS location
!  call Grid_getCellCoords(KAXIS, blockID,CENTER,useGuardCell,zCoordsCell,1)
!  call Grid_getCellCoords(JAXIS, blockID,CENTER,useGuardCell,yCoordsCell,MAXCELLS)
!  call Grid_getCellCoords(IAXIS, blockID,CENTER,useGuardCell,xCoordsCell,MAXCELLS)

  ! now fill the master arrays

  do k = blockRange(LOW,KAXIS), blockRange(HIGH,KAXIS)
     if (NDIM == 3) then
        iPosition(3) = k
        zz = zCoordsCell(k)
     endif

     do j = blockRange(LOW,JAXIS),blockRange(HIGH,JAXIS)
        if (NDIM >= 2) then
           iPosition(2) = j
           yy = yCoordsCell(j)
        endif

        do i = blockRange(LOW,IAXIS),blockRange(HIGH,IAXIS)
           iPosition(1) = i
           xx = xCoordsCell(i)

           call User_initBlobCell(xx, yy, zz, unkArr)

           ! store the values
           ! fill the flash arrays

           call Grid_putPointData(blockID,CENTER,TEMP_VAR,EXTERIOR,iPosition,unkArr(TEMP_VAR))
           call Grid_putPointData(blockID,CENTER,DENS_VAR,EXTERIOR,iPosition,unkArr(DENS_VAR))
           call Grid_putPointData(blockID,CENTER,PRES_VAR,EXTERIOR,iPosition,unkArr(PRES_VAR))
           call Grid_putPointData(blockID,CENTER,EINT_VAR,EXTERIOR,iPosition,unkArr(EINT_VAR))
           call Grid_putPointData(blockID,CENTER,ENER_VAR,EXTERIOR,iPosition,unkArr(ENER_VAR))
           call Grid_putPointData(blockID,CENTER,GAMC_VAR,EXTERIOR,iPosition,unkArr(GAMC_VAR))
           call Grid_putPointData(blockID,CENTER,GAME_VAR,EXTERIOR,iPosition,unkArr(GAME_VAR))
           call Grid_putPointData(blockID,CENTER,VELX_VAR,EXTERIOR,iPosition,unkArr(VELX_VAR))
           call Grid_putPointData(blockID,CENTER,VELY_VAR,EXTERIOR,iPosition,unkArr(VELY_VAR))
           call Grid_putPointData(blockID,CENTER,VELZ_VAR,EXTERIOR,iPosition,unkArr(VELZ_VAR))

           do n = SPECIES_BEGIN,SPECIES_END
              call Grid_putPointData(blockID,CENTER,n,EXTERIOR,iPosition,unkArr(n))
           enddo

           !..end of 3d loops
        enddo  ! end of k loop
     enddo     ! end of j loop
  enddo        ! end of i loop

  ! cleanup
  !deallocate(rvec)
  deallocate(xCoordsCell)
  deallocate(yCoordsCell)
  deallocate(zCoordsCell)


  return
end subroutine Simulation_initBlock
