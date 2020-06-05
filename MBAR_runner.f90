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
  real(kind=fp_kind), allocatable :: weights(:,:)
  real(kind=fp_kind), allocatable :: extWeights(:,:)
  real(kind=fp_kind), allocatable :: extTheta(:,:)
  integer(kind=4), allocatable :: extNSnapshotsInSimulation(:)

  integer(kind=4) :: nbins
  real(kind=fp_kind) :: binmin,binmax,binwidth
  real(kind=fp_kind), allocatable :: bincenters(:)
  real(kind=fp_kind), allocatable :: pmf(:), pmfSE(:)
  real(kind=fp_kind), allocatable :: reweighting_entropy(:)
  integer(kind=4), allocatable :: nSnapshotsInBin(:)
  real(kind=fp_kind), allocatable :: sumOfWeightsInBin(:)

  integer(kind=4) :: IndexW, IndexS, IndexB
  integer(kind=4) :: JndexS

  character(len=80) :: metaFile
  character(len=80) :: targetHamiltonianEnergyFile
  character(len=80) :: targetHamiltonianPmfFile
  character(len=80) :: targetWeightsFile

  write(6,'(1X,A)')'Number of simulations:'
  read*,nSimulations
  write(6,'(1X,A)')'Name of the meta file:'
  read*,metaFile
  write(6,'(1X,A)')'Number of bins for output:'
  read*,nbins
  write(6,'(1X,A)')'Minimum and maximum of the bins:'
  read*,binmin,binmax
  write(6,'(1X,A)')'Name of the file containing the energies of the target Hamiltonian:'
  read*,targetHamiltonianEnergyFile
  write(6,'(1X,A)')'Name of the out file containing the PMF of the target Hamiltonian:'
  read*,targetHamiltonianPmfFile
  write(6,'(1X,A)')'Name of the out file containing the weights of the target Hamiltonian:'
  read*, targetWeightsFile

  binwidth = (binmax - binmin)/nbins

  open(id_meta_file,file= trim(metaFile))
  call readSimulationInfo()
  write(6,'(1X,A,I4,A,I8)')'Number of simulations:', nSimulations, '   Total Number of snapshots:', totalNumSnapshots
  close(id_meta_file)

  allocate(simulatedReducedHamiltonian(nSimulations))
  do IndexW = 1, nSimulations
    simulatedReducedHamiltonian(IndexW)%ownSimulation => simulations(IndexW)
    call simulatedReducedHamiltonian(IndexW)%init( &
          &  simulatedReducedHamiltonian(IndexW)%ownSimulation%beta, &
          &  totalNumSnapshots, &
          &  simulatedReducedHamiltonian(IndexW)%ownSimulation%nSnapshots)
  end do
  do IndexW = 1, nSimulations
    call simulatedReducedHamiltonian(IndexW)%processTrajectories
  end do

! allocate data space for MBAR
  allocate(reducedEnergies(totalNumSnapshots,nSimulations))
  allocate(nSnapshotsInSimulation(nSimulations))
  allocate(freeEnergies(nSimulations))
  do IndexW = 1, nSimulations
    do JndexS = 1, totalNumSnapshots
      reducedEnergies(JndexS,IndexW) = simulatedReducedHamiltonian(IndexW)%reducedEnergies(jndexS)
    end do
    nSnapshotsInSimulation(IndexW) = simulatedReducedHamiltonian(IndexW)%nOwnSnapshots
  end do
  allocate(weights(totalNumSnapshots,nSimulations))
  allocate(simulationsTheta(nSimulations,nSimulations))
! The main purpose of MBAR is to obtain the free energies for the simulated
! (reduced)Hamiltonians and associated uncertainties
  call MBAR(50,1.D-7,nSimulations,totalNumSnapshots,reducedEnergies,nSnapshotsInSimulation,freeEnergies, &
             & weights,simulationsTheta,'NRM')
  do IndexW = 1, nSimulations
    simulatedReducedHamiltonian(IndexW)%freeEnergy = &
       &  freeEnergies(IndexW)/simulatedReducedHamiltonian(IndexW)%beta
    simulatedReducedHamiltonian(IndexW)%freeEnergySD2 = &
       & simulationsTheta(IndexW,IndexW) - 2*simulationsTheta(IndexW,1) + &
       & simulationsTheta(1,1)
    simulatedReducedHamiltonian(IndexW)%freeEnergySD2 = &
       & simulatedReducedHamiltonian(IndexW)%freeEnergySD2 / &
       & simulatedReducedHamiltonian(IndexW)%beta
    write(6,'(A,I3,A,F10.3,A,G13.5)')' Free energy of Hamiltonian ', IndexW,  &
       &  ' is: ', simulatedReducedHamiltonian(IndexW)%freeEnergy, '+-', &
       &  sqrt(simulatedReducedHamiltonian(IndexW)%freeEnergySD2)
  end do  
  do IndexW = 1, nSimulations
    simulatedReducedHamiltonian(IndexW)%weights(:) = weights(:, IndexW)
  end do
  deallocate(weights)

  allocate(weights(totalNumSnapshots,1))
  allocate(targetReducedEnergies(totalNumSnapshots))


! Compute info for a target Hamiltonian
  call targetReducedHamiltonian%init(targetBeta,totalNumSnapshots,0)
  call targetReducedHamiltonian%processTrajectories(coordonly=.true.)
  open(id_target_energy_file, file=trim(targetHamiltonianEnergyFile))
  call targetReducedHamiltonian%readInEnergy(id_target_energy_file)
  close(id_target_energy_file)

  targetReducedEnergies = targetReducedHamiltonian%reducedEnergies
  call MBAR_weight(nSimulations,totalNumSnapshots,reducedEnergies,nSnapshotsInSimulation,freeEnergies, &
           & targetReducedEnergies,weights,targetReducedHamiltonian%freeenergy)
  targetReducedHamiltonian%weights(:)=weights(:,1)
  open(id_target_weights_file,file=targetWeightsFile)
  write(id_target_weights_file,'(I,E12.5,1X,F10.5)')(IndexS, targetReducedHamiltonian%weights(IndexS), &
        &  targetReducedHamiltonian%snapshots(IndexS)%coordinate, IndexS = 1, totalNumSnapshots)
  close(id_target_weights_file)


! Compute PMF for the target Hamiltonian
  allocate(bincenters(nbins))
  allocate(pmf(nbins))
  allocate(pmfSE(nbins))
  allocate(reweighting_entropy(nbins))
  allocate(nSnapshotsInBin(nbins))
  allocate(sumOfWeightsInBin(nbins))
  do IndexB = 1, nbins
    bincenters(IndexB) = binmin + (IndexB-0.5)*binwidth
  end do 
  pmf = 0.d0
  pmfSE = 0.d0
  reweighting_entropy = 0.d0
  nSnapshotsInBin = 0 
  sumOfWeightsInBin = 0.d0
  do IndexS = 1, totalNumSnapshots
    IndexB = int(( targetReducedHamiltonian%snapshots(IndexS)%coordinate - binmin )/binwidth) + 1
    if(IndexB > nbins .or. IndexB < 1) cycle
    pmf(IndexB) = pmf(IndexB) + targetReducedHamiltonian%weights(IndexS)
    reweighting_entropy(IndexB) = reweighting_entropy(IndexB) + &
      targetReducedHamiltonian%weights(IndexS)*log(targetReducedHamiltonian%weights(IndexS))
    nSnapshotsInBin(IndexB) = nSnapshotsInBin(IndexB) + 1
    sumOfWeightsInBin(IndexB) = sumOfWeightsInBin(IndexB) + &
      targetReducedHamiltonian%weights(IndexS)
  end do
  pmf = -log(pmf) / targetReducedHamiltonian%beta
  pmf = pmf - minval(pmf)
  reweighting_entropy(:) = -(reweighting_entropy(:)/sumOfWeightsInBin(:)-log(sumOfWeightsInBin(:))) &
       & /log(dble(nSnapshotsInBin(:)))

  call targetReducedHamiltonian%bootstrap(100,nbins,binmin,binwidth,pmf,pmfSE)
  open(id_target_pmf_file , file = targetHamiltonianPmfFile)
  do IndexB = 1, nbins
    write(id_target_pmf_file,'(F8.3,3F9.3)')bincenters(IndexB), pmf(IndexB), pmfSE(IndexB), reweighting_entropy(IndexB)
  end do
  close(id_target_pmf_file)

! delete data space
  do IndexW = 1, nSimulations
    call simulatedReducedHamiltonian(IndexW)%destroy()
  end do
  call targetReducedHamiltonian%destroy()
  call deleteSimulationInfo()

  deallocate(simulatedReducedHamiltonian)
  deallocate(reducedEnergies)
  deallocate(nSnapshotsInSimulation)
  deallocate(freeEnergies)
  deallocate(targetReducedEnergies)
  deallocate(weights)
  deallocate(bincenters)
  deallocate(pmf)
  deallocate(pmfSE)
  deallocate(reweighting_entropy)
  deallocate(nSnapshotsInBin)
  deallocate(sumOfWeightsInBin)
end program MBAR_caller
