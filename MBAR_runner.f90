program MBAR_runner
  use MBAR_m
  use io_m
  use bin_m
  use simulation_m
  use reducedHamiltonian_m
  implicit none
  
  real(kind=fp_kind) :: targetBeta=1.d0/(298*kB)

  integer(kind=4) :: nbins
  real(kind=fp_kind) :: binmin, binmax
  integer(kind=4) :: iGaussSmooth

  real(kind=fp_kind), allocatable :: reducedEnergies(:,:)
  real(kind=fp_kind), allocatable :: freeEnergies(:)
  real(kind=fp_kind), allocatable :: targetReducedEnergies(:)

  real(kind=fp_kind), allocatable :: weights(:,:)
  real(kind=fp_kind), allocatable :: covFreeEnergies(:,:)
  integer(kind=4), allocatable :: nSnapshotsInSimulation(:)

  real(kind=fp_kind), allocatable :: overlap(:,:)
  real(kind=fp_kind), allocatable :: crossCorr(:,:)

  integer(kind=4) :: IndexW, IndexS, IndexB
  integer(kind=4) :: JndexW, JndexS

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
  write(6,'(1X,A)')'Perform Gaussian smoothing on the data? 0: No. 1: Yes.'
  read*,iGaussSmooth

  open( id_meta_file, file = trim(metaFile), status = 'old')
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
  allocate(crossCorr(nSimulations,nSimulations))
  allocate(overlap(nSimulations,nSimulations))
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
    write(6,'(A,I3,A,F8.3,A,F8.3)')' Free energy (kcal/mol) of Hamiltonian ', IndexW, ' is: ', &
       & simulatedReducedHamiltonian(IndexW)%freeEnergy/simulatedReducedHamiltonian(IndexW)%beta, & 
       & ' +- ', &
       & simulatedReducedHamiltonian(IndexW)%freeEnergySD/simulatedReducedHamiltonian(IndexW)%beta
  end do  

  forall(IndexW = 1 : nSimulations)simulatedReducedHamiltonian(IndexW)%weights(:) = weights(:, IndexW)

  do IndexW = 1, nSimulations
    do JndexW = IndexW, nSimulations
      crossCorr(IndexW, JndexW) = crossCorrelationBetweenHs(totalNumSnapshots, & 
                 & simulatedReducedHamiltonian(IndexW)%weights, &
                 & simulatedReducedHamiltonian(JndexW)%weights)
      crossCorr(JndexW, IndexW) = crossCorr(IndexW, JndexW)
      overlap(IndexW, JndexW) = MobleyOverlap(totalNumSnapshots, simulatedReducedHamiltonian(IndexW)%ownSimulation%nSnapshots, &
                 & simulatedReducedHamiltonian(IndexW)%weights, &
                 & simulatedReducedHamiltonian(JndexW)%weights)
      overlap(JndexW, IndexW) = MobleyOverlap(totalNumSnapshots, simulatedReducedHamiltonian(JndexW)%ownSimulation%nSnapshots, &
                 & simulatedReducedHamiltonian(JndexW)%weights, &
                 & simulatedReducedHamiltonian(IndexW)%weights)
    end do
  end do
  write(6,'(A)') 'Mobley overlap matrix:'
  do IndexW = 1, nSimulations
    write(6,'(10F8.3)') overlap(IndexW,:)
  end do
  write(6,'(A)') 'Cross correlation between simulations:'
  do IndexW = 1, nSimulations
    write(6,'(10F8.3)') crossCorr(IndexW,:)
  end do
  deallocate(weights)
  deallocate(covFreeEnergies)
  deallocate(crossCorr)
  deallocate(overlap)

  allocate(weights(totalNumSnapshots,nSimulations+1))
  allocate(targetReducedEnergies(totalNumSnapshots))
  allocate(crossCorr(nSimulations,1))
  allocate(overlap(nSimulations,1))
  
! gather data for a target Hamiltonian
  call targetReducedHamiltonian%init(targetBeta, totalNumSnapshots, 0)
  call targetReducedHamiltonian%processTrajectories(coordonly=.true.)
  open( id_target_energy_file, file = trim(targetHamiltonianEnergyFile), status = 'old')
  call targetReducedHamiltonian%readInEnergy(id_target_energy_file)
  close(id_target_energy_file)

! compute weight for each sample under the target Hamiltonian
  targetReducedEnergies = targetReducedHamiltonian%reducedEnergies
  forall(IndexW = 1 : nSimulations) weights(:, IndexW) = simulatedReducedHamiltonian(IndexW)%weights(:)
  write(6,'(A)') 'Computing weights for all the configurations under the target Hamiltonian'
  call MBAR_weight(nSimulations, totalNumSnapshots, reducedEnergies, nSnapshotsInSimulation, &
           & freeEnergies, targetReducedEnergies, &
           & weights, targetReducedHamiltonian%freeenergy, targetReducedHamiltonian%freeenergySD)
  write(6,'(A,F8.3,1X,A,F8.3)')' Free energy (kcal/mol) of the target Hamiltonian is ', &
           & targetReducedHamiltonian%freeenergy/targetReducedHamiltonian%beta, '+/-', &
           & targetReducedHamiltonian%freeenergySD/targetReducedHamiltonian%beta

  targetReducedHamiltonian%weights(:) = weights(:,nSimulations+1)

  open(id_target_weights_file, file = targetWeightsFile)
  write(id_target_weights_file,'(I8,E12.5,1X,F10.5)')(IndexS, &
        &  targetReducedHamiltonian%weights(IndexS)/sum(targetReducedHamiltonian%weights(:)), &
        &  targetReducedHamiltonian%snapshots(IndexS)%coordinate, IndexS = 1, totalNumSnapshots)
  close(id_target_weights_file)

  do IndexW = 1, nSimulations
    overlap(IndexW, 1) = MobleyOverlap(totalNumSnapshots, simulatedReducedHamiltonian(IndexW)%ownSimulation%nSnapshots, &
        & simulatedReducedHamiltonian(IndexW)%weights, targetReducedHamiltonian%weights)
    crossCorr(IndexW, 1) = crossCorrelationBetweenHs(totalNumSnapshots, &
        & simulatedReducedHamiltonian(IndexW)%weights, targetReducedHamiltonian%weights)
  end do
  write(6,'(A)') 'Mobley overlap between the target Hamiltonian and the simulated Hamiltonian'
  write(6,'(10F8.3)') overlap(:, 1)
  write(6,'(A)') 'Cross correlation between the target Hamiltonian and the simulated Hamiltonian'
  write(6,'(10F8.3)') crossCorr(:, 1)

  deallocate(weights)
  deallocate(nSnapshotsInSimulation)
  deallocate(overlap)
  deallocate(crossCorr)

! Compute PMF for the target Hamiltonian
  call targetReducedHamiltonian%computePMF(nbins, binmin, binmax, iGaussSmooth)

! bootstrapping for the calculations of STD Err for PMF
  call targetReducedHamiltonian%bootstrap(nbins, binmin, binmax, 100)

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
end program MBAR_runner
