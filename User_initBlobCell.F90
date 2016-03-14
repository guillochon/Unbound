subroutine User_initBlobCell(xx, yy, zz, unkArr)
    use Simulation_data, ONLY: sim_smallRho, sim_smallx, sim_radiusPerturb, sim_usePseudo1d, &
       sim_rhoPerturb, sim_tempPerturb, sim_velPerturb, sim_rhoAmbient, sim_tempAmbient, &
       sim_xCenterPerturb, sim_yCenterPerturb, sim_zCenterPerturb, &
       pos_vec, dens_vec, npts, sim_velyPerturb, sim_nChunks, sim_cloudRadius, &
       sim_chunkEllipticity, sim_chunkAngle, sim_chunkSeparation, &
       temp_vec, sim_xf, sim_xfAmb, sim_xfHot, sim_virTemp, sim_virDens, &
       sim_noiseAmplitude, sim_ranSeed, sim_rhoCloud, sim_tempCloud, sim_xCenterCloud
    use Eos_interface, ONLY : Eos
    use Grid_data, ONLY: gr_meshMe
    implicit none

#include "Flash.h"
#include "constants.h"
#include "Eos.h"

    real, intent(IN) :: xx, yy, zz
    real, dimension(NUNK_VARS), intent(INOUT) :: unkArr

    logical :: inChunk
    integer, allocatable :: dists(:)
    integer :: i, l
    real :: chunkRad, chunkA, chunkB, chunkC, xOffset, yOffset, coeff, distx, &
            disty, distz, xprime, yprime, zprime, cloudDist
    real, dimension(EOS_NUM)  :: eosData
    real, dimension(3) :: rvec

    unkArr = 0.d0
    
    allocate(dists(sim_nChunks))

    inChunk = .false.
    if (sim_nChunks .eq. 1) then
       sim_nChunks = 1
       chunkRad = sim_radiusPerturb
    else
       chunkRad = sim_radiusPerturb / (dble(sim_nChunks)**(1.d0/3.d0))
    endif

    if (sim_chunkEllipticity .ne. 0.d0) then
       chunkA = chunkRad*(1.d0-sim_chunkEllipticity)**(-2.d0/3.d0)
       chunkB = chunkA*(1.d0-sim_chunkEllipticity)
       chunkC = chunkB
    else
       chunkA = chunkRad
       chunkB = chunkRad
       chunkC = chunkRad
    endif
    xOffset = -sim_chunkSeparation*chunkA*dcos(sim_chunkAngle)
    yOffset = -sim_chunkSeparation*chunkA*dsin(sim_chunkAngle)
    ! compute the distance from the center
    do i = 1, sim_nChunks
       if (inChunk) exit
       if (sim_chunkEllipticity .eq. 0.d0) then
           if (NDIM .EQ. 1) then
              dists(i) = xx - (sim_xCenterPerturb + (i-1)*xOffset)
           else if (NDIM .EQ. 2) then
              if (sim_usePseudo1d) then
                 dists(i) = xx - sim_xCenterPerturb
              else
                 dists(i) = sqrt((xx - (sim_xCenterPerturb + (i-1)*xOffset))**2 + & 
                      &          (yy - (sim_yCenterPerturb + (i-1)*yOffset))**2)
              endif
           elseif (NDIM .EQ. 3) then
              if (sim_usePseudo1d) then
                 dists(i) = xx - sim_xCenterPerturb
              else
                 dists(i) = sqrt((xx - (sim_xCenterPerturb + (i-1)*xOffset))**2 + & 
                      &          (yy - (sim_yCenterPerturb + (i-1)*yOffset))**2 + &
                      &          (zz - sim_zCenterPerturb)**2)
              endif
           endif

           ! set the temperature, density, and velocity
           if (dists(i) .LE. chunkRad) then
              inChunk = .true.
              unkArr(DENS_VAR) = sim_rhoPerturb
              unkArr(TEMP_VAR) = sim_tempPerturb
              unkArr(VELX_VAR) = sim_velPerturb
              unkArr(VELY_VAR) = 0.d0
              unkArr(VELZ_VAR) = 0.d0
              unkArr(SPECIES_BEGIN:SPECIES_END) = sim_xf
           endif
       else
           if (NDIM .NE. 3) then
              call Driver_abortFlash('Error: Chunk Ellipticity only supports 3 dimensions')
           endif
           distx = xx - (sim_xCenterPerturb + (i-1)*xOffset)
           disty = yy - (sim_yCenterPerturb + (i-1)*yOffset)
           distz = zz - sim_zCenterPerturb
            
           xprime = distx*dcos(sim_chunkAngle) + disty*dsin(sim_chunkAngle)
           yprime = disty*dcos(sim_chunkAngle) - distx*dsin(sim_chunkAngle)
           zprime = distz
           if ((xprime/chunkA)**2.d0 + &
               (yprime/chunkB)**2.d0 + &
               (zprime/chunkC)**2.d0 .lt. 1.d0) then
               inChunk = .true.
               unkArr(DENS_VAR) = sim_rhoPerturb
               unkArr(TEMP_VAR) = sim_tempPerturb
               unkArr(VELX_VAR) = sim_velPerturb
               unkArr(VELY_VAR) = 0.d0
               unkArr(VELZ_VAR) = 0.d0
               unkArr(SPECIES_BEGIN:SPECIES_END) = sim_xf
           endif
       endif
    enddo
    if (.not. inChunk) then
        if (NDIM .EQ. 1) then
           cloudDist = xx - sim_xCenterCloud
        else if (NDIM .EQ. 2) then
           cloudDist = sqrt((xx - sim_xCenterCloud)**2 + & 
                 &          (yy - sim_yCenterPerturb)**2)
        elseif (NDIM .EQ. 3) then
           cloudDist = sqrt((xx - sim_xCenterCloud)**2 + & 
                 &          (yy - sim_yCenterPerturb)**2 + &
                 &          (zz - sim_zCenterPerturb)**2)
        endif

        ! set the temperature, density, and x-velocity
        if (cloudDist .LE. sim_cloudRadius) then
             unkArr(DENS_VAR) = sim_rhoCloud
             unkArr(TEMP_VAR) = sim_tempCloud
        else
             unkArr(DENS_VAR) = sim_rhoAmbient
             unkArr(TEMP_VAR) = sim_tempAmbient
        endif
        unkArr(VELX_VAR) = 0.d0
        unkArr(VELY_VAR) = 0.d0
        unkArr(VELZ_VAR) = 0.d0
        unkArr(SPECIES_BEGIN:SPECIES_END) = sim_xfAmb
    endif

    !  Need input of density and temperature
    eosData(EOS_DENS) = unkArr(DENS_VAR)
    eosData(EOS_TEMP) = unkArr(TEMP_VAR)

    call Eos(MODE_DENS_TEMP,1,eosData,unkArr(SPECIES_BEGIN:SPECIES_END))

    unkArr(DENS_VAR) = eosData(EOS_DENS)
    unkArr(TEMP_VAR) = eosData(EOS_TEMP)
    unkArr(PRES_VAR) = eosData(EOS_PRES)
    unkArr(EINT_VAR) = eosData(EOS_EINT)
    unkArr(GAMC_VAR) = eosData(EOS_GAMC)
    unkArr(ENER_VAR) = unkArr(EINT_VAR) + 0.5*sum(unkArr(VELX_VAR:VELZ_VAR)**2)
    unkArr(GAME_VAR) = unkArr(PRES_VAR)/(unkArr(ENER_VAR)*unkArr(DENS_VAR)) + 1.0

    deallocate(dists)
    return
end subroutine User_initBlobCell

