subroutine User_calculateAtmosphere()
    use Simulation_data
    use RuntimeParameters_interface, ONLY : RuntimeParameters_get
    use Grid_data, ONLY : gr_nBlockX
    use Eos_interface, ONLY : Eos
    implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

    real :: dz, envEntr, gconst
    integer :: i, imarker
    real, dimension(EOS_NUM)  :: eosData

    call RuntimeParameters_get('gconst', gconst)

    dz = (sim_xmax - sim_xmin)*(gr_nBlockX + 1.d0)/(npts - 1.d0)

    pos_vec(1) = sim_xmin - 4.d0*dz
    dens_vec = sim_virDens
    temp_vec = sim_virTemp

    do i = 2, npts
        pos_vec(i) = pos_vec(i-1) + dz
    enddo

    if (sim_xAtm .gt. sim_xmax) then
        call Driver_abortFlash('Error, xAtm greater than xmax')
    endif

    eosData(EOS_DENS) = sim_virDens
    eosData(EOS_TEMP) = sim_virTemp
    call Eos(MODE_DENS_TEMP,1,eosData,sim_xfHot)
    pres_vec = eosData(EOS_PRES)
    envEntr = eosData(EOS_ENTR)

    do i = 1, npts
        if (pos_vec(i) .ge. sim_xAtm) then
            imarker = i
            exit
        endif
    enddo

    do i = imarker-1, 1, -1
        pres_vec(i) = pres_vec(i+1) - gconst*dz*dens_vec(i+1) * (1.d0 - (sim_rotVel/sim_vCirc)**2)
        eosData(EOS_DENS) = dens_vec(i+1) !guess
        eosData(EOS_PRES) = pres_vec(i)
        eosData(EOS_ENTR) = envEntr
        call Eos(MODE_PRES_ENTR,1,eosData,sim_xfHot)
        dens_vec(i) = eosData(EOS_DENS)
        pres_vec(i) = eosData(EOS_PRES)
        temp_vec(i) = eosData(EOS_TEMP)
    enddo

    do i = imarker, npts
        pres_vec(i) = pres_vec(i-1) + gconst*dz*dens_vec(i-1) * (1.d0 - (sim_rotVel/sim_vCirc)**2)
        eosData(EOS_DENS) = dens_vec(i-1) !guess
        eosData(EOS_PRES) = pres_vec(i)
        eosData(EOS_TEMP) = sim_tempAmbient
        call Eos(MODE_PRES_TEMP,1,eosData,sim_xfAmb)
        dens_vec(i) = eosData(EOS_DENS)
        pres_vec(i) = eosData(EOS_PRES)
        temp_vec(i) = eosData(EOS_TEMP)
    enddo

    do i = 1, npts
        dens_vec(i) = max(dens_vec(i), sim_smallRho)
    enddo
end subroutine
