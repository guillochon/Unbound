subroutine User_updateEnvelopeComposition(blockCount, blockList)
    use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getCellCoords, Grid_getBlkPtr, &
                               Grid_releaseBlkPtr
    use Simulation_data, ONLY : sim_xfHot, sim_yMin, sim_yMax
    use Grid_data, ONLY : gr_meshComm, gr_meshMe
    implicit none

#include "Flash.h"
#include "constants.h"
#include "Eos.h"
#include "Flash_mpi.h"

    integer, INTENT(in)                        :: blockCount
    integer, INTENT(in), DIMENSION(blockCount)  :: blockList

    integer :: thisBlock, blockID, xSizeCoord, ySizeCoord, zSizeCoord, ierr, i, j, k
    real, dimension(:), allocatable :: xCoord, yCoord, zCoord
    real, dimension(SPECIES_BEGIN:SPECIES_END) :: new_xn, gnew_xn
    integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
    logical :: getGuardCells = .true.
    real, pointer, dimension(:,:,:,:)            :: solnData

    new_xn = 0.d0

    do thisBlock = 1, blockCount
        blockID = blockList(thisBlock)
        ! get dimensions/limits and coordinates
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
        xSizeCoord = blkLimitsGC(HIGH,IAXIS)
        ySizeCoord = blkLimitsGC(HIGH,JAXIS)
        zSizeCoord = blkLimitsGC(HIGH,KAXIS)
        !! allocate space for dimensions
        allocate(xCoord(xSizeCoord))
        allocate(yCoord(ySizeCoord))
        allocate(zCoord(zSizeCoord))

        call Grid_getCellCoords(IAXIS,blockID,CENTER,getGuardCells,xCoord,xSizeCoord)
        call Grid_getCellCoords(JAXIS,blockID,CENTER,getGuardCells,yCoord,ySizeCoord)
        call Grid_getCellCoords(KAXIS,blockID,CENTER,getGuardCells,zCoord,zSizeCoord)

        ! Get a pointer to solution data 
        call Grid_getBlkPtr(blockID,solnData)
        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                 ! Hard-coded to 5% of height of box for now
                 if (yCoord(j) .gt. 0.05d0*(sim_yMax - sim_yMin)) cycle
                 new_xn = new_xn + (1.d0-solnData(CORE_MSCALAR,i,j,k))*solnData(DENS_VAR,i,j,k)*&
                          solnData(SPECIES_BEGIN:SPECIES_END,i,j,k)
              enddo
           enddo
        enddo

        call Grid_releaseBlkPtr(blockID,solnData)
        deallocate(xCoord)
        deallocate(yCoord)
        deallocate(zCoord)
    enddo

    call MPI_ALLREDUCE(new_xn, gnew_xn, NSPECIES, FLASH_REAL, MPI_SUM, gr_meshComm, ierr)
    sim_xfHot = gnew_xn / sum(gnew_xn)
    if (gr_meshMe .eq. MASTER_PE) then
        write(*, '(A)', advance='no') 'sim_xfHot:'
        do i = SPECIES_BEGIN, SPECIES_END
            write(*, '(F15.10)', advance='no') sim_xfHot(i)
        enddo
        write(*, *)
    endif
end subroutine User_updateEnvelopeComposition
