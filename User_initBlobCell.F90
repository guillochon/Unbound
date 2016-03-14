subroutine User_initBlobCell(xx, yy, zz, unkArr)
    use Simulation_data, ONLY: sim_smallRho, sim_smallx, sim_radiusPerturb, sim_usePseudo1d, &
       sim_rotVel, sim_xAtm, &
       sim_rhoPerturb, sim_tempPerturb, sim_velPerturb, &
       sim_xCenterPerturb, sim_yCenterPerturb, sim_zCenterPerturb, &
       pos_vec, dens_vec, npts, sim_velyPerturb, sim_nChunks, &
       sim_chunkEllipticity, sim_chunkAngle, sim_chunkSeparation, &
       temp_vec, sim_xf, sim_xfAmb, sim_xfHot, pres_vec, sim_virTemp, sim_virDens, &
       sim_noiseAmplitude, sim_ranSeed
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
            disty, distz, xprime, yprime, zprime
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

           ! set the temperature, density, and x-velocity
           if (dists(i) .LE. chunkRad) then
              inChunk = .true.
              unkArr(TEMP_VAR) = sim_tempPerturb
              do l = 1, npts
                  if (xx .lt. pos_vec(l)) then
                      coeff = (xx - pos_vec(l-1))/(pos_vec(l) - pos_vec(l-1))
                      unkArr(PRES_VAR) = pres_vec(l-1) + coeff*(pres_vec(l) - pres_vec(l-1))
                      exit
                  endif
              enddo
              eosData(EOS_DENS) = sim_rhoPerturb !guess
              eosData(EOS_TEMP) = unkArr(TEMP_VAR)
              eosData(EOS_PRES) = unkArr(PRES_VAR)
              call Eos(MODE_PRES_TEMP,1,eosData,sim_xf)
              unkArr(DENS_VAR) = eosData(EOS_DENS)
              unkArr(VELX_VAR) = sim_velPerturb*dcos(sim_chunkAngle)
              unkArr(VELY_VAR) = sim_velPerturb*dsin(sim_chunkAngle)
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
               unkArr(TEMP_VAR) = sim_tempPerturb

               do l = 1, npts
                   if (xx .lt. pos_vec(l)) then
                       coeff = (xx - pos_vec(l-1))/(pos_vec(l) - pos_vec(l-1))
                       unkArr(PRES_VAR) = pres_vec(l-1) + coeff*(pres_vec(l) - pres_vec(l-1))
                       exit
                   endif
               enddo

               eosData(EOS_DENS) = sim_rhoPerturb !guess
               eosData(EOS_TEMP) = unkArr(TEMP_VAR)
               eosData(EOS_PRES) = unkArr(PRES_VAR)
               call Eos(MODE_PRES_TEMP,1,eosData,sim_xf)
               unkArr(DENS_VAR) = eosData(EOS_DENS)
               call sim_ranmar(sim_ranSeed, rvec, 3)
               unkArr(VELX_VAR) = sim_velPerturb*(dcos(sim_chunkAngle) + sim_noiseAmplitude*(1.0 - 2.0 * rvec(1)))
               unkArr(VELY_VAR) = sim_velPerturb*(dsin(sim_chunkAngle) + sim_noiseAmplitude*(1.0 - 2.0 * rvec(2)))
               unkArr(VELZ_VAR) = sim_velPerturb*(sim_noiseAmplitude*(1.0 - 2.0 * rvec(3)))
               unkArr(SPECIES_BEGIN:SPECIES_END) = sim_xf
           endif
       endif
    enddo
    if (.not. inChunk) then
        do l = 1, npts
            if (xx .lt. pos_vec(l)) then
                coeff = (xx - pos_vec(l-1))/(pos_vec(l) - pos_vec(l-1))
                unkArr(DENS_VAR) = (dens_vec(l-1) + coeff*(dens_vec(l) - dens_vec(l-1)))
                unkArr(TEMP_VAR) = (temp_vec(l-1) + coeff*(temp_vec(l) - temp_vec(l-1)))
                exit
            endif
        enddo
        if (xx .lt. sim_xAtm) then
            unkArr(SPECIES_BEGIN:SPECIES_END) = sim_xfHot
            call sim_ranmar(sim_ranSeed, rvec, 3)
            unkArr(VELX_VAR) = sim_rotVel*(sim_noiseAmplitude*(1.0 - 2.0 * rvec(1)))
            unkArr(VELY_VAR) = sim_rotVel*(1.d0 + sim_noiseAmplitude*(1.0 - 2.0 * rvec(2)))
            unkArr(VELZ_VAR) = sim_rotVel*(sim_noiseAmplitude*(1.0 - 2.0 * rvec(3)))
        else
            unkArr(VELX_VAR) = 0.d0
            unkArr(VELY_VAR) = sim_rotVel !0.d0
            unkArr(VELZ_VAR) = 0.d0
            unkArr(SPECIES_BEGIN:SPECIES_END) = sim_xfAmb
            unkArr(CORE_MSCALAR) = 1.d0
        endif
    endif

    !  Need input of density and temperature
    eosData(EOS_TEMP) = unkArr(TEMP_VAR)
    eosData(EOS_DENS) = unkArr(DENS_VAR)

    call Eos(MODE_DENS_TEMP,1,eosData,unkArr(SPECIES_BEGIN:SPECIES_END))

    unkArr(TEMP_VAR) = eosData(EOS_TEMP)
    unkArr(DENS_VAR) = eosData(EOS_DENS)
    unkArr(PRES_VAR) = eosData(EOS_PRES)
    unkArr(EINT_VAR) = eosData(EOS_EINT)
    unkArr(GAMC_VAR) = eosData(EOS_GAMC)
    unkArr(ENER_VAR) = unkArr(EINT_VAR) + 0.5*sum(unkArr(VELX_VAR:VELZ_VAR)**2)
    unkArr(GAME_VAR) = unkArr(PRES_VAR)/(unkArr(ENER_VAR)*unkArr(DENS_VAR)) + 1.0

    deallocate(dists)
    return
end subroutine User_initBlobCell

