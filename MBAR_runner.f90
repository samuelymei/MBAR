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
  real(kind=fp_kind), allocatable :: theta(:,:)
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
!  nSimulations = 73
!  nbins = 35
!  binmin = 1.5d0
!  binmax = 5.0d0
!  targetHamiltonianEnergyFile = 'pm6_ascii_6h2o_energy.out'
  binwidth = (binmax - binmin)/nbins

  open(id_meta_file,file= trim(metaFile))
  call readSimulationInfo()
  write(6,'(1X,A,I,A,I)')'Number of simulations:', nSimulations, ' Total Number of snapshots:', totalNumSnapshots
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
  do IndexW = 1, nSimulations
    do JndexS = 1, totalNumSnapshots
      reducedEnergies(JndexS,IndexW) = simulatedReducedHamiltonian(IndexW)%reducedEnergies(jndexS)
    end do
    nSnapshotsInSimulation(IndexW) = simulatedReducedHamiltonian(IndexW)%nOwnSnapshots
  end do
  allocate(weights(totalNumSnapshots,nSimulations))
  allocate(theta(nSimulations,nSimulations))
  call MBAR(50,1.D-7,nSimulations,totalNumSnapshots,reducedEnergies,nSnapshotsInSimulation,freeEnergies, &
             & weights,theta,'NRM')
  do IndexW = 1, nSimulations
    simulatedReducedHamiltonian(IndexW)%freeEnergy = freeEnergies(IndexW)
  end do  
  do IndexW = 1, nSimulations
    simulatedReducedHamiltonian(IndexW)%weights(:) = weights(:, IndexW)
  end do
!  allocate(extWeights(totalNumSnapshots,nSimulations+1))
!  allocate(extTheta(nSimulations+1,nSimulations+1))
!  extWeights(1:totalNumSnapshots,1:nSimulations)=weights(1:totalNumSnapshots,1:nSimulations)
  deallocate(weights)
  deallocate(theta)

  allocate(weights(totalNumSnapshots,1))
  allocate(targetReducedEnergies(totalNumSnapshots))

!  do IndexW = 1, nSimulations
!    targetReducedEnergies(:) = simulatedReducedHamiltonian(IndexW)%reducedEnergies(:)
!    call MBAR_weight(nSimulations,totalNumSnapshots,reducedEnergies,nSnapshotsInSimulation,freeEnergies, &
!            & targetReducedEnergies,simulatedReducedHamiltonian(IndexW)%freeenergy,weights)
!    simulatedReducedHamiltonian(IndexW)%weights(:)=weights(:,1)
!  end do

  call targetReducedHamiltonian%init(targetBeta,totalNumSnapshots,0)
  call targetReducedHamiltonian%processTrajectories(coordonly=.true.)
  open(id_target_energy_file, file=trim(targetHamiltonianEnergyFile))
  call targetReducedHamiltonian%readInEnergy(id_target_energy_file)
  close(id_target_energy_file)
!  jndexS = 0
!  do indexW = 1, nSimulations
!    do indexS = 1, simulations(indexW)%nSnapshots
!      jndexS = jndexS + 1
!      targetReducedHamiltonian%reducedEnergies(jndexS) = &
!            & simulations(indexW)%snapshots(indexS)%energyUnbiased * &
!            & targetReducedHamiltonian%beta
!    end do
!  end do
  targetReducedEnergies = targetReducedHamiltonian%reducedEnergies
  call MBAR_weight(nSimulations,totalNumSnapshots,reducedEnergies,nSnapshotsInSimulation,freeEnergies, &
           & targetReducedEnergies,weights,targetReducedHamiltonian%freeenergy)
  targetReducedHamiltonian%weights(:)=weights(:,1)
  open(id_target_weights_file,file=targetWeightsFile)
  write(id_target_weights_file,'(I,E12.5,1X,F10.5)')(IndexS, targetReducedHamiltonian%weights(IndexS), &
        &  targetReducedHamiltonian%snapshots(IndexS)%coordinate, IndexS = 1, totalNumSnapshots)
  close(id_target_weights_file)

!  extWeights(1:totalNumSnapshots,nSimulations+1)=weights(1:totalNumSnapshots,1)
!  allocate(extNSnapshotsInSimulation(nSimulations+1))
!  extNSnapshotsInSimulation(1:nSimulations)=nSnapshotsInSimulation(1:nSimulations)
!  extNSnapshotsInSimulation(nSimulations+1)=0
!  call weight2theta(totalNumSnapshots,nSimulations+1,extNSnapshotsInSimulation,extWeights,extTheta)
!  write(876,'(8E10.3)')extTheta
!  targetReducedHamiltonian%weightsSE=0.d0
!  do IndexS = 1, totalNumSnapshots
!    do IndexW = 1, nSimulations
!      targetReducedHamiltonian%weightsSE(IndexS)=targetReducedHamiltonian%weightsSE(IndexS) + &
!        & (extWeights(IndexS,nSimulations+1)*extWeights(IndexS,IndexW)*nSnapshotsInSimulation(IndexW))**2* &
!        & (targetReducedHamiltonian%freeenergy-simulatedReducedHamiltonian(IndexW)%freeEnergy)**2 * &
!        & (extTheta(nSimulations+1,nSimulations+1)+extTheta(IndexW,IndexW)-2*extTheta(nSimulations+1,IndexW))
!    end do
!    targetReducedHamiltonian%weightsSE(IndexS)=sqrt(targetReducedHamiltonian%weightsSE(IndexS))
!    write(74,*)targetReducedHamiltonian%weightsSE(IndexS)
!  end do
!  deallocate(extNSnapshotsInSimulation)

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
!    pmfSE(IndexB) = pmfSE(IndexB) + targetReducedHamiltonian%weightsSE(IndexS)**2
    reweighting_entropy(IndexB) = reweighting_entropy(IndexB) + &
      targetReducedHamiltonian%weights(IndexS)*log(targetReducedHamiltonian%weights(IndexS))
    nSnapshotsInBin(IndexB) = nSnapshotsInBin(IndexB) + 1
    sumOfWeightsInBin(IndexB) = sumOfWeightsInBin(IndexB) + &
      targetReducedHamiltonian%weights(IndexS)
  end do
!  write(987,'(G12.5,I6)')(sumOfWeightsInBin(IndexB),nSnapshotsInBin(IndexB), IndexB = 1, nbins)
  pmf = -log(pmf) / targetReducedHamiltonian%beta
!  pmfSE = sqrt(pmfSE) * exp(pmf)
  pmf = pmf - minval(pmf)
  reweighting_entropy(:) = -(reweighting_entropy(:)/sumOfWeightsInBin(:)-log(sumOfWeightsInBin(:))) &
       & /log(dble(nSnapshotsInBin(:)))
!  reweighting_entropy(:) = -reweighting_entropy(:)/log(dble(nSnapshotsInBin(:)))

  call targetReducedHamiltonian%bootstrap(100,nbins,binmin,binwidth,pmf,pmfSE)
  open(id_target_pmf_file , file = targetHamiltonianPmfFile)
  do IndexB = 1, nbins
    write(id_target_pmf_file,'(F8.3,3F9.3)')bincenters(IndexB), pmf(IndexB), pmfSE(IndexB), reweighting_entropy(IndexB)
  end do
  close(id_target_pmf_file)
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
