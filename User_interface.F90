!!****h* source/User/User_interface
!!
!! This is the header file for the User module
!! that defines its public interfaces.
!!***
Module User_interface
  implicit none
#include "Flash.h"
#include "constants.h"
  interface
     subroutine User_calculateAtmosphere()
       implicit none
     end subroutine User_calculateAtmosphere

     subroutine User_updateEnvelopeComposition(blockCount, blockList)
       implicit none
       integer, intent(IN) :: blockCount
       integer, dimension(blockCount), intent(IN):: blockList
     end subroutine User_updateEnvelopeComposition

     subroutine User_initBlobCell(xx, yy, zz, unkArr)
       implicit none
       real, intent(IN) :: xx, yy, zz
       real, dimension(NUNK_VARS), intent(INOUT) :: unkArr
     end subroutine User_initBlobCell
  end interface

end Module User_interface

