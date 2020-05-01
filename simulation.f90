module simulation_m
  use precision_m
  use snapshot_m
  private
  integer(kind=4), public :: nSimulation

  type :: simulation_t
    real(kind=fp_kind) :: beta  ! inverse temperature of this simulation
    integer(kind=4) :: nSnapshots ! Number of snapshots in this simulation
    type ( snapshot_t ), pointer :: snapshots(:) ! point to the snapshots
    contains
      procedure :: destroy
  end type simulation_t
  type ( simulation_t ), allocatable, public :: simulations(:) ! to save the information of all the simulations
  public :: readSimulationInfo, deleteSimulationInfo

contains

  subroutine readSimulationInfo(fid)
    use constant_m, only : kB
    implicit none
    integer(kind=4), intent(in) :: fid
    real(kind=fp_kind) :: simulationTemperature
    integer(kind=4) :: indexW, indexS
    allocate(simulations(nSimulation))
    do indexW = 1, nSimulation
      read(fid, *) simulationTemperature, simulations(indexW)%nSnapshots
      simulations(indexW)%beta = 1.d0 / ( kB * simulationTemperature )
      write(6,'(A,I3,A,F10.6,A,I10)') ' Simulation ', indexW, ':    beta = ', simulations(indexW)%beta, &
         & ', Number of snapshots:', simulations(indexW)%nSnapshots
      allocate(simulations(indexW)%snapshots(simulations(indexW)%nSnapshots))
    end do
  end subroutine readSimulationInfo

  subroutine destroy(this)
    implicit none
    class(simulation_t) :: this
    integer(kind=4) :: indexW
    do indexW = 1, nSimulation
      deallocate(this%snapshots)
    end do
  end subroutine destroy

  subroutine deleteSimulationInfo
    implicit none
    integer(kind=4) :: indexW
    do indexW = 1, nSimulation
      call simulations(indexW)%destroy()
    end do
    if(allocated(simulations))deallocate(simulations)
  end subroutine deleteSimulationInfo

end module simulation_m

