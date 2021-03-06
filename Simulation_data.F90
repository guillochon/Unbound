!!****if* source/Simulation/SimulationMain/Unbound/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data
!!
!! DESCRIPTION
!!  Store the simulation data for Cellular setup
!!   
!! ARGUMENTS
!!
!! PARAMETERS   
!!
!!    xhe4               mass fraction of he4
!!    xc12               mass fraction of c12
!!    xo16               mass fraction of o16
!!    rhoAmbient         density of the cold upstream material 
!!    tempAmbient        temperature of the cold upstream material
!!    rhoUDS         density of the post shock material
!!    tempUDS        temperature of the post shock material
!!    velxUDS        x-velocity of the post shock material
!!    radiusUDS      distance below which the perturbation is applied
!!    xCenterUDS     origin of the of the perturbation
!!    yCenterUDS     origin of the of the perturbation
!!    zCenterUDS     origin of the of the perturbation
!!                          copied along the y and z directions
!!                       .false. for a spherical configuration
!!    noiseAmplitude     amplitude of the white noise added to the perturbation
!!    noiseDistance      distances above and below radiusUDS get noise added
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
!!   No arguments.  All data passed by "use Simulation_data"
!! 
!! NOTES
!!   See paper: Timmes, FX; Zingale, M; Olson, K; Fryxell, B; The Astrophysical
!!               Journal, Nov. 10, 2000 : 543: 938-954
!!
!!***

module Simulation_data

  implicit none

#include "constants.h"
#include "Flash.h"



! ---- Runtime Parameters -------------------------

     real, save :: sim_smallRho, sim_smallx  ! minimum values of parameters
     real, save :: sim_radiusUDS     ! influence range of perturbation

     ! initial mass fractions
     real, save :: sim_xh1, sim_xhe4, sim_xc12, sim_xo16, sim_xne20, sim_xsi28, sim_xfe54
     real, save :: sim_ah1, sim_ahe4, sim_ac12, sim_ao16, sim_ane20, sim_asi28, sim_afe54
     real, save :: sim_hh1, sim_hhe4, sim_hc12, sim_ho16, sim_hne20, sim_hsi28, sim_hfe54

     ! ambient parameters
     real, save :: sim_rhoAmbient, sim_tempAmbient

     ! UDS parameters
     real, save :: sim_rhoUDS, sim_tempUDS, sim_velMedianUDS, sim_velSpreadUDS, &
                   sim_velExpansionUDS, sim_rhoSpreadUDS, sim_rhoScaleHeightsUDS, &
                   sim_xCenterUDS, sim_yCenterUDS, sim_zCenterUDS
     integer, save :: sim_nChunks
     real, save :: sim_chunkEllipticity, sim_chunkSeparation, sim_chunkAngle

     ! cloud parameters
     real, save :: sim_xCenterCloud, sim_yCenterCloud, sim_rhoCloud, sim_tempCloud, &
                   sim_cloudRadius, sim_cloudScaleHeights
     
     ! physical ranges of domain
     real, save :: sim_xmin, sim_xmax, sim_ymin, sim_ymax, sim_zmin, sim_zmax  

     real, dimension(SPECIES_BEGIN:SPECIES_END) :: sim_xfAmb, sim_xf, sim_xfHot

     integer, parameter :: npts = 10000
     real, dimension(npts) :: pos_vec, dens_vec, temp_vec
! -------------------------------------------------------------------------------------


     integer, save :: sim_meshMe
end module Simulation_data
