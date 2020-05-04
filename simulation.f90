module simulation_m
  use precision_m
  use constant_m
  use io_m
  use snapshot_m
  private
  integer(kind=4), public :: nSimulations, totalNumSnapshots

  type :: simulation_t
    real(kind=fp_kind) :: beta = 1.d0/(kB*298.0)  ! inverse temperature of this simulation
    real(kind=fp_kind) :: restraintStrength = 0.d0 ! strength of the restraint
    real(kind=fp_kind) :: restraintCenter ! center of the restraint
    integer(kind=4) :: nSnapshots ! Number of snapshots in this simulation
    type ( snapshot_t ), allocatable :: snapshots(:) ! point to the snapshots
    contains
      procedure :: destroy
  end type simulation_t
  type ( simulation_t ), allocatable, target, public :: simulations(:) ! to save the information of all the simulations
  public :: readSimulationInfo, deleteSimulationInfo
  public :: simulation_t
contains

  subroutine readSimulationInfo()
    use constant_m, only : kB
    use io_m
    implicit none
    real(kind=fp_kind) :: simulationTemperature
    integer(kind=4) :: indexW, indexS
    character(len=80) :: dataFile
    integer(kind=4) :: itmp
    if(idebug == 1) write(*,*) 'Entering readSimulationInfo'
    allocate(simulations(nSimulations))
    do indexW = 1, nSimulations
      read(id_meta_file, *) dataFile, simulations(indexW)%nSnapshots, simulationTemperature, &
                  & simulations(indexW)%restraintCenter, &
                  & simulations(indexW)%restraintStrength
      simulations(indexW)%beta = 1.d0 / ( kB * simulationTemperature )
      write(6,'(A,I3,A,F10.6,A,I10,A,A)') ' Simulation ', indexW, ':    beta = ', simulations(indexW)%beta, &
         & ', Number of snapshots:', simulations(indexW)%nSnapshots, &
         & ', Data file: ', dataFile
      allocate(simulations(indexW)%snapshots(simulations(indexW)%nSnapshots))
      open(id_data_file, file=dataFile, status = 'old')
      do indexS = 1, simulations(indexW)%nSnapshots
        read(id_data_file,*) itmp, simulations(indexW)%snapshots(indexS)%coordinate, &
                 & simulations(indexW)%snapshots(indexS)%energyUnbiased
      !  write(*,*)'   ', simulations(indexW)%snapshots(indexS)%coordinate, &
      !     & simulations(indexW)%snapshots(indexS)%energyUnbiased
      end do
      close(id_data_file)
    end do
    totalNumSnapshots = sum(simulations(:)%nSnapshots)
  end subroutine readSimulationInfo

  subroutine destroy(this)
    implicit none
    class(simulation_t) :: this
    integer(kind=4) :: indexW
    if(idebug == 1) write(*,*) 'Entering simulations%destroy'
    do indexW = 1, nSimulations
      if(allocated(this%snapshots))deallocate(this%snapshots)
    end do
  end subroutine destroy

  subroutine deleteSimulationInfo
    implicit none
    integer(kind=4) :: indexW
    if(idebug == 1) write(*,*) 'Entering deleteSimulationInfo'
    do indexW = 1, nSimulations
      call simulations(indexW)%destroy()
    end do
    if(allocated(simulations))deallocate(simulations)
  end subroutine deleteSimulationInfo

end module simulation_m

