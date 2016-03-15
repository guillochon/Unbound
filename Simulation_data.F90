!!****if* source/Simulation/SimulationMain/Cellular/Simulation_data
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
     real, save :: sim_radiusPerturb     ! influence range of perturbation

     ! initial mass fractions
     real, save :: sim_xh1, sim_xhe4, sim_xc12, sim_xo16, sim_xne20, sim_xsi28, sim_xfe54
     real, save :: sim_ah1, sim_ahe4, sim_ac12, sim_ao16, sim_ane20, sim_asi28, sim_afe54
     real, save :: sim_hh1, sim_hhe4, sim_hc12, sim_ho16, sim_hne20, sim_hsi28, sim_hfe54

     ! ambient parameters
     real, save :: sim_rhoAmbient, sim_tempAmbient
     ! perturbed parameters
     real, save :: sim_rhoPerturb, sim_tempPerturb, sim_velPerturb, sim_velyPerturb

     real, save :: sim_xCenterCloud, sim_rhoCloud, sim_tempCloud, &
                   sim_cloudRadius, sim_cloudScaleHeights
     
     !  noise added to perturbation 
     real, save :: sim_noiseAmplitude, sim_noiseDistance

     logical, save :: sim_usePseudo1d    ! initial conditions depend only upon distance in x

     ! center of perturbation
     real, save :: sim_xCenterPerturb, sim_yCenterPerturb, sim_zCenterPerturb

     integer, save :: sim_nChunks

     real, save :: sim_chunkEllipticity, sim_chunkSeparation, sim_chunkAngle
     real, save :: sim_virTemp, sim_virDens, sim_vCirc

     real, save :: sim_rotVel

     integer, save :: sim_ranSeed

     ! physical ranges of domain
     real, save :: sim_xmin, sim_xmax, sim_ymin, sim_ymax, sim_zmin, sim_zmax  

     real, dimension(SPECIES_BEGIN:SPECIES_END) :: sim_xfAmb, sim_xf, sim_xfHot

     integer, parameter :: npts = 10000
     real, dimension(npts) :: pos_vec, dens_vec, temp_vec
! -------------------------------------------------------------------------------------


     integer, save :: sim_meshMe
end module Simulation_data
