program MBAR_caller
  use MBAR_m
  use io_m
  use simulation_m
  use reducedHamiltonian_m
  implicit none
  
  type (reducedHamiltonian_t), allocatable, target :: simulatedReducedHamiltonian(:)
  type (reducedHamiltonian_t), target :: targetReducedHamiltonian

  real(kind=fp_kind) :: targetBeta=1.d0/(298*kB)

  real(kind=fp_kind), allocatable :: reducedEnergies(:,:)
  integer(kind=4), allocatable :: nSnapshotsInSimulation(:)
  real(kind=fp_kind), allocatable :: freeEnergies(:)
  real(kind=fp_kind), allocatable :: targetReducedEnergies(:)
  real(kind=fp_kind), allocatable :: weights(:,:), weightsStdDev(:,:)
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
  write(6,'(1X,A)')'Target temperature:'
  read*,targetBeta
  targetBeta = 1.d0/ (targetBeta*kB)
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
  forall(IndexW = 1:nSimulations, JndexS = 1: totalNumSnapshots) &
    & reducedEnergies(JndexS,IndexW) = simulatedReducedHamiltonian(IndexW)%reducedEnergies(jndexS)
  forall(IndexW = 1:nSimulations) nSnapshotsInSimulation(IndexW) = simulatedReducedHamiltonian(IndexW)%nOwnSnapshots

  allocate(weights(totalNumSnapshots,nSimulations))
  allocate(covFreeEnergies(nSimulations,nSimulations))
! The main purpose of MBAR is to obtain the free energies for the simulated
! (reduced)Hamiltonians and associated uncertainties
  call MBAR(50,1.D-7,nSimulations,totalNumSnapshots,reducedEnergies,nSnapshotsInSimulation,freeEnergies, &
             & weights,covFreeEnergies,'NRM')
 
  forall (IndexW = 1 : nSimulations) 
    simulatedReducedHamiltonian(IndexW)%freeEnergy = freeEnergies(IndexW)
    simulatedReducedHamiltonian(IndexW)%freeEnergySD = sqrt( covFreeEnergies(IndexW,IndexW) - &
       & 2*covFreeEnergies(IndexW,1) + covFreeEnergies(1,1) )
  end forall
  do IndexW = 1, nSimulations
    write(6,'(A,I3,A,F10.3,A,G13.5)')' Free energy (kcal/mol) of Hamiltonian ', IndexW, ' is: ', &
       & simulatedReducedHamiltonian(IndexW)%freeEnergy/simulatedReducedHamiltonian(IndexW)%beta, & 
       & ' +- ', &
       & simulatedReducedHamiltonian(IndexW)%freeEnergySD/simulatedReducedHamiltonian(IndexW)%beta
  end do  

  forall(IndexW = 1 : nSimulations)simulatedReducedHamiltonian(IndexW)%weights(:) = weights(:, IndexW)
  deallocate(weights)

  allocate(weights(totalNumSnapshots,1))
  allocate(targetReducedEnergies(totalNumSnapshots))

! gather data for a target Hamiltonian
  call targetReducedHamiltonian%init(targetBeta, totalNumSnapshots, 0)
  call targetReducedHamiltonian%processTrajectories(coordonly=.true.)
  open(id_target_energy_file, file = trim(targetHamiltonianEnergyFile))
  call targetReducedHamiltonian%readInEnergy(id_target_energy_file)
  close(id_target_energy_file)

! compute weight for each sample under the target Hamiltonian
  targetReducedEnergies = targetReducedHamiltonian%reducedEnergies
  call MBAR_weight(nSimulations,totalNumSnapshots,reducedEnergies,nSnapshotsInSimulation,freeEnergies, &
           &  targetReducedEnergies, weights,targetReducedHamiltonian%freeenergy)
  write(6,'(A,F10.3)')' Free energy (kcal/mol) of the target Hamiltonian is ', &
      targetReducedHamiltonian%freeenergy/targetReducedHamiltonian%beta

  targetReducedHamiltonian%weights(:) = weights(:,1)
  open(id_target_weights_file, file = targetWeightsFile)
  write(id_target_weights_file,'(I,E12.5,1X,F10.5)')(IndexS, targetReducedHamiltonian%weights(IndexS), &
        &  targetReducedHamiltonian%snapshots(IndexS)%coordinate, IndexS = 1, totalNumSnapshots)
  close(id_target_weights_file)

  deallocate(weights)
  deallocate(covFreeEnergies)
  deallocate(nSnapshotsInSimulation)

! Compute PMF for the target Hamiltonian
  allocate(bincenters(nbins))
  allocate(pmf(nbins))
  allocate(pmfSE(nbins))
  allocate(reweighting_entropy(nbins))
  allocate(nSnapshotsInBin(nbins))
  allocate(sumOfWeightsInBin(nbins))

  allocate(weights(totalNumSnapshots,nSimulations+nbins))
  allocate(covFreeEnergies(nSimulations+nbins,nSimulations+nbins))
  allocate(nSnapshotsInSimulation(nSimulations+nbins))

  weights = 0.d0
  forall(IndexW = 1 : nSimulations)weights(:, IndexW) = simulatedReducedHamiltonian(IndexW)%weights(:)
  nSnapshotsInSimulation = 0
  forall(IndexW = 1 : nSimulations)nSnapshotsInSimulation(IndexW) = simulatedReducedHamiltonian(IndexW)%nOwnSnapshots
  forall(IndexB = 1:nbins) bincenters(IndexB) = binmin + (IndexB-0.5)*binwidth
  pmf = 0.d0
  pmfSE = 0.d0
  reweighting_entropy = 0.d0
  nSnapshotsInBin = 0 
  sumOfWeightsInBin = 0.d0
  do IndexS = 1, totalNumSnapshots
    IndexB = int(( targetReducedHamiltonian%snapshots(IndexS)%coordinate - binmin )/binwidth) + 1
    if(IndexB > nbins .or. IndexB < 1) cycle
    weights(IndexS,nSimulations+IndexB) = targetReducedHamiltonian%weights(IndexS)
    pmf(IndexB) = pmf(IndexB) + targetReducedHamiltonian%weights(IndexS)
    reweighting_entropy(IndexB) = reweighting_entropy(IndexB) + &
      & targetReducedHamiltonian%weights(IndexS)*log(targetReducedHamiltonian%weights(IndexS))
    nSnapshotsInBin(IndexB) = nSnapshotsInBin(IndexB) + 1
    sumOfWeightsInBin(IndexB) = sumOfWeightsInBin(IndexB) + &
      & targetReducedHamiltonian%weights(IndexS)
  end do

  forall(IndexB = 1 : nbins) weights(:,nSimulations+IndexB) = weights(:,nSimulations+IndexB)/sum(weights(:,nSimulations+IndexB))

  call ComputCovMatFromWeights(totalNumSnapshots,nSimulations+nbins,nSnapshotsInSimulation,weights,covFreeEnergies)

  pmf = -log(pmf)
  pmf = pmf - pmf(1)
  forall(IndexB=1:nbins) pmfSE(IndexB) = &
      & sqrt(  covFreeEnergies(nSimulations+IndexB,nSimulations+IndexB) &
          &  + covFreeEnergies(nSimulations+1,nSimulations+1) &
        &    - 2*covFreeEnergies(nSimulations+1,nSimulations+IndexB) )
  reweighting_entropy(:) = -(reweighting_entropy(:)/sumOfWeightsInBin(:)-log(sumOfWeightsInBin(:))) &
       & /log(dble(nSnapshotsInBin(:)))
  open(id_target_pmf_file , file = targetHamiltonianPmfFile)
  do IndexB = 1, nbins
    write(id_target_pmf_file,'(F8.3,3F9.3)')bincenters(IndexB), pmf(IndexB)/targetReducedHamiltonian%beta, &
         & pmfSE(IndexB)/targetReducedHamiltonian%beta, reweighting_entropy(IndexB)
  end do
  close(id_target_pmf_file)

! bootstrapping for the calculations of STD Err for PMF
  call targetReducedHamiltonian%bootstrap(100,nbins,binmin,binwidth,pmf,pmfSE)
  do IndexB = 1, nbins
    write(90,'(F8.3,3F9.3)')bincenters(IndexB), pmf(IndexB)/targetReducedHamiltonian%beta, &
        & pmfSE(IndexB)/targetReducedHamiltonian%beta, reweighting_entropy(IndexB)
  end do

! delete data space
  do IndexW = 1, nSimulations
    call simulatedReducedHamiltonian(IndexW)%destroy()
  end do
  call targetReducedHamiltonian%destroy()
  call deleteSimulationInfo()

  if(allocated(simulatedReducedHamiltonian))deallocate(simulatedReducedHamiltonian)
  if(allocated(reducedEnergies))deallocate(reducedEnergies)
  if(allocated(nSnapshotsInSimulation))deallocate(nSnapshotsInSimulation)
  if(allocated(freeEnergies))deallocate(freeEnergies)
  if(allocated(targetReducedEnergies))deallocate(targetReducedEnergies)
  if(allocated(weights))deallocate(weights)
  if(allocated(bincenters))deallocate(bincenters)
  if(allocated(pmf))deallocate(pmf)
  if(allocated(pmfSE))deallocate(pmfSE)
  if(allocated(reweighting_entropy))deallocate(reweighting_entropy)
  if(allocated(nSnapshotsInBin))deallocate(nSnapshotsInBin)
  if(allocated(sumOfWeightsInBin))deallocate(sumOfWeightsInBin)
end program MBAR_caller
