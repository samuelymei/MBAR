program MBAR_caller
  use MBAR_m
  use io_m
  use simulation_m
  use reducedHamiltonian_m
  implicit none
  
  type (reducedHamiltonian_t), allocatable, target :: simulatedReducedHamiltonian(:)
  type (reducedHamiltonian_t), target :: targetReducedHamiltonian

  real(kind=fp_kind) :: targetBeta=1/(298*kB)

  real(kind=fp_kind), allocatable :: reducedEnergies(:,:)
  integer(kind=4), allocatable :: nSnapshotsInSimulation(:)
  real(kind=fp_kind), allocatable :: freeEnergies(:)
  real(kind=fp_kind), allocatable :: targetReducedEnergies(:)
  real(kind=fp_kind), allocatable :: weights(:)

  integer(kind=4) :: IndexW, IndexS
  integer(kind=4) :: JndexS
  nSimulations = 73

  open(id_meta_file,file= 'meta.dat')
  call readSimulationInfo()
  close(id_meta_file)
  allocate(simulatedReducedHamiltonian(nSimulations))
  do IndexW = 1, nSimulations
    simulatedReducedHamiltonian(IndexW)%ownSimulation => simulations(IndexW)
    call simulatedReducedHamiltonian(IndexW)%init( &
          &  simulatedReducedHamiltonian(IndexW)%ownSimulation%beta, &
          &  totalNumSnapshots,& 
          &  simulatedReducedHamiltonian(IndexW)%ownSimulation%nSnapshots)
  end do
  do IndexW = 1, nSimulations
    call simulatedReducedHamiltonian(IndexW)%processTrajectories
  end do

  allocate(reducedEnergies(totalNumSnapshots,nSimulations))
  allocate(nSnapshotsInSimulation(nSimulations))
  allocate(freeEnergies(nSimulations))
  allocate(targetReducedEnergies(totalNumSnapshots))
  allocate(weights(totalNumSnapshots))
  do IndexW = 1, nSimulations
    do JndexS = 1, totalNumSnapshots
      reducedEnergies(JndexS,IndexW) = simulatedReducedHamiltonian(IndexW)%reducedEnergies(jndexS)
    end do
    nSnapshotsInSimulation(IndexW) = simulatedReducedHamiltonian(IndexW)%nOwnSnapshots
  end do
  write(66,*)reducedEnergies
  call MBAR(50,1.D-7,nSimulations,totalNumSnapshots,reducedEnergies,nSnapshotsInSimulation,freeEnergies)
  do IndexW = 1, nSimulations
    simulatedReducedHamiltonian(IndexW)%freeEnergy = freeEnergies(IndexW)
  end do  

  call targetReducedHamiltonian%init(targetBeta,totalNumSnapshots,0)
  jndexS = 0
  do indexW = 1, nSimulations
    do indexS = 1, simulations(indexW)%nSnapshots
      jndexS = jndexS + 1
      targetReducedHamiltonian%reducedEnergies(jndexS) = &
            & simulations(indexW)%snapshots(indexS)%energyUnbiased * &
            & targetReducedHamiltonian%beta
    end do
  end do
  targetReducedEnergies = targetReducedHamiltonian%reducedEnergies
  call MBAR_weight(nSimulations,totalNumSnapshots,reducedEnergies,nSnapshotsInSimulation,freeEnergies,targetReducedEnergies,weights)

  call deleteSimulationInfo()
  deallocate(simulatedReducedHamiltonian)
  deallocate(reducedEnergies)
  deallocate(nSnapshotsInSimulation)
  deallocate(freeEnergies)
  deallocate(targetReducedEnergies)
  deallocate(weights)
end program MBAR_caller
