module io_m
  implicit none
  public
  integer(kind=4),parameter :: id_meta_file=11
  integer(kind=4),parameter :: id_data_file=12
  integer(kind=4),parameter :: id_target_energy_file=13
  integer(kind=4),parameter :: id_target_weights_file=14
  integer(kind=4),parameter :: id_target_pmf_file=15
  character(len=80) :: metaFile
  character(len=80) :: targetHamiltonianEnergyFile
  character(len=80) :: targetHamiltonianPmfFile
  character(len=80) :: targetWeightsFile
end module io_m
