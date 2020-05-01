module snapshot_m
  use precision_m
  implicit none
  public 
  integer(kind=4) :: nSnapshots
  type :: snapshot_t
    real(kind=fp_kind) :: coordinate ! coordinate (in a reduced dimension)
    real(kind=fp_kind) :: energyUnbiased ! Energy
  end type snapshot_t
end module snapshot_m


