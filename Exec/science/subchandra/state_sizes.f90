module state_sizes_module
   use network, only : nspec, naux
   implicit none
   integer, parameter :: nadv = 0
   integer, parameter :: NQAUX = 4
   integer, parameter :: NVAR = 7 + nadv + nspec + naux
   integer, parameter :: NGDNV = 6
   integer, parameter :: NQ = 8 + nadv + nspec + naux
   integer, parameter :: QVAR = 8 + nadv + nspec + naux
end module state_sizes_module
