module reducedHamiltonian_m
  use precision_m
  use constant_m
  use io_m
  use snapshot_m
  use simulation_m
  implicit none
  private
  type :: reducedHamiltonian_t
    real(kind=fp_kind) :: beta
    integer(kind=4) :: totalNumSnapshots
    integer(kind=4) :: nOwnSnapshots
    type (snapshot_t), allocatable :: snapshots(:)
    real(kind=fp_kind), allocatable :: reducedEnergies(:)
    real(kind=fp_kind), allocatable :: weights(:)
    real(kind=fp_kind), allocatable :: weightsSE(:)
    real(kind=fp_kind) :: freeEnergy, freeEnergySD
    type (simulation_t), pointer :: ownSimulation =>Null()
    contains
      procedure :: destroy
      procedure :: init
      procedure :: processTrajectories
      procedure :: readInEnergy
      procedure :: computePMF
      procedure :: bootstrap
  end type reducedHamiltonian_t

  type (reducedHamiltonian_t), allocatable, public, target :: simulatedReducedHamiltonian(:)
  type (reducedHamiltonian_t), public, target :: targetReducedHamiltonian


  public :: reducedHamiltonian_t
  public :: computePMF
  public :: MobleyOverlap
  public :: crossCorrelationBetweenHs
  
  contains
!    function constructor(beta,TotalNumSnapshots,nOwnSnapshots)
!      implicit none
!      type (reducedHamiltonian_t) :: constructor
!      real(kind=fp_kind), intent(in) :: beta
!      integer(kind=4), intent(in) :: TotalNumSnapshots
!      integer(kind=4), intent(in) :: nOwnSnapshots
!      call constructor%init(beta,TotalNumSnapshots,nOwnSnapshots)
!    end function constructor

    subroutine init(this,beta,TotalNumSnapshots,nOwnSnapshots)
      class(reducedHamiltonian_t) :: this
      real(kind=fp_kind), intent(in) :: beta
      integer(kind=4), intent(in) :: TotalNumSnapshots
      integer(kind=4), intent(in) :: nOwnSnapshots
      if(idebug == 1) write(*,*) 'Entering reducedHamiltonian%init'
      this%beta = beta
      this%TotalNumSnapshots = TotalNumSnapshots
      this%nOwnSnapshots = nOwnSnapshots
      allocate(this%snapshots(this%TotalNumSnapshots))
      allocate(this%reducedEnergies(this%TotalNumSnapshots))
      allocate(this%weights(this%TotalNumSnapshots))
      allocate(this%weightsSE(this%TotalNumSnapshots))
    end subroutine init

    subroutine processTrajectories(this,coordonly)
      implicit none
      class(reducedHamiltonian_t) :: this
      logical, intent(in), optional :: coordonly
      logical :: crdonly=.false.
      integer(kind=4) :: indexW, indexS
      integer(kind=4) :: jndexS
      if(idebug == 1) write(*,*) 'Entering reducedHamiltonian%processTrajectories'

      if(present(coordonly))crdonly=coordonly
      jndexS = 0
      do indexW = 1, nSimulations
        do indexS = 1, simulations(indexW)%nSnapshots
          jndexS = jndexS + 1
          this%snapshots(jndexS)%coordinate = &
               & simulations(indexW)%snapshots(indexS)%coordinate
          if(.not.crdonly)then
            this%reducedEnergies(jndexS) = &
                 & simulations(indexW)%snapshots(indexS)%energyUnbiased &
                 & + 0.5*this%ownSimulation%restraintStrength &
                 & * (this%snapshots(jndexS)%coordinate - this%ownSimulation%restraintCenter)**2
            this%reducedEnergies(jndexS) = this%reducedEnergies(jndexS) * this%beta
          end if
        end do
      end do
    end subroutine processTrajectories

    subroutine readInEnergy(this,fileid)
      implicit none
      class (reducedHamiltonian_t) :: this
      integer(kind=4), intent(in) :: fileid
      real(kind=fp_kind) :: avgEner
      integer(kind=4) :: IndexS
      do IndexS = 1, this%totalNumSnapshots
        read(fileid,*) this%reducedEnergies(IndexS)
        this%reducedEnergies(IndexS) = this%reducedEnergies(IndexS) * this%beta
      end do
      avgEner = sum(this%reducedEnergies(:))/totalNumSnapshots
      this%reducedEnergies(:) = this%reducedEnergies(:) - avgEner
    end subroutine readInEnergy

    subroutine destroy(this)
      class(reducedHamiltonian_t) :: this
      if(idebug == 1) write(*,*) 'Entering reducedHamiltonian%destroy'
      if(allocated(this%snapshots))deallocate(this%snapshots)
      if(allocated(this%reducedEnergies))deallocate(this%reducedEnergies)
      if(allocated(this%weights))deallocate(this%weights)
      if(allocated(this%weightsSE))deallocate(this%weightsSE)
    end subroutine destroy

    subroutine computePMF(this)
      use io_m
      use precision_m
      use bin_m
      use MBAR_m
      use simulation_m
      implicit none
      class(reducedHamiltonian_t) :: this
      real(kind=fp_kind), allocatable :: extWeights(:,:)
      real(kind=fp_kind), allocatable :: extCovFreeEnergies(:,:)
      integer(kind=4), allocatable :: extNSnapshotsInSimulation(:)
      integer(kind=4) :: IndexW, IndexS, IndexB
      integer(kind=4) :: JndexS

      allocate(extWeights(totalNumSnapshots,nSimulations+nbins))
      allocate(extCovFreeEnergies(nSimulations+nbins,nSimulations+nbins))
      allocate(extNSnapshotsInSimulation(nSimulations+nbins))

      extWeights = 0.d0
      forall(IndexW = 1 : nSimulations)extWeights(:, IndexW) = simulatedReducedHamiltonian(IndexW)%weights(:)
      extNSnapshotsInSimulation = 0
      forall(IndexW = 1 : nSimulations)extNSnapshotsInSimulation(IndexW) = simulatedReducedHamiltonian(IndexW)%nOwnSnapshots
      bins(:)%pmf = 0.d0
      bins(:)%pmfSE = 0.d0
      bins(:)%reweightingEntropy = 0.d0
      bins(:)%nSnapshotsInBin = 0 
      bins(:)%sumOfWeightsInBin = 0.d0
      do IndexS = 1, totalNumSnapshots
        IndexB = int(( this%snapshots(IndexS)%coordinate - binmin )/binwidth) + 1
        if(IndexB > nbins .or. IndexB < 1) cycle
        extWeights(IndexS,nSimulations+IndexB) = this%weights(IndexS)
        bins(IndexB)%pmf = bins(IndexB)%pmf + this%weights(IndexS)
        bins(IndexB)%reweightingEntropy = bins(IndexB)%reweightingEntropy + &
          & this%weights(IndexS)*log(this%weights(IndexS))
        bins(IndexB)%nSnapshotsInBin = bins(IndexB)%nSnapshotsInBin + 1
        bins(IndexB)%sumOfWeightsInBin = bins(IndexB)%sumOfWeightsInBin + &
          & this%weights(IndexS)
      end do
  
      forall(IndexB = 1 : nbins) extWeights(:,nSimulations+IndexB) = &
            & extWeights(:,nSimulations+IndexB)/sum(extWeights(:,nSimulations+IndexB))
  
      call ComputCovMatFromWeights(totalNumSnapshots,nSimulations+nbins,extNSnapshotsInSimulation,extWeights,extCovFreeEnergies)
  
      bins(:)%pmf = -log(bins(:)%pmf)
      bins(:)%pmf = bins(:)%pmf - bins(1)%pmf
      forall(IndexB=1:nbins) bins(IndexB)%pmfSE = &
          & sqrt(  extCovFreeEnergies(nSimulations+IndexB,nSimulations+IndexB) &
              &  + extCovFreeEnergies(nSimulations+1,nSimulations+1) &
            &    - 2*extCovFreeEnergies(nSimulations+1,nSimulations+IndexB) )
      bins(:)%reweightingEntropy = -(bins(:)%reweightingEntropy/bins(:)%sumOfWeightsInBin-log(bins(:)%sumOfWeightsInBin)) &
           & /log(dble(bins(:)%nSnapshotsInBin))
      open(id_target_pmf_file , file = targetHamiltonianPmfFile)
      write(6,'(1X,A)')'Potential of mean force under the target Hamiltonian (kcal/mol)'
      do IndexB = 1, nbins
        write(id_target_pmf_file,'(F8.3,3F9.3)')bins(IndexB)%bincenter, bins(IndexB)%pmf/this%beta, &
             & bins(IndexB)%pmfSE/this%beta, bins(IndexB)%reweightingEntropy
        write(6,'(F8.3,3F9.3)')bins(IndexB)%bincenter, bins(IndexB)%pmf/this%beta, &
             & bins(IndexB)%pmfSE/this%beta, bins(IndexB)%reweightingEntropy
      end do
      close(id_target_pmf_file)
      deallocate(extWeights)
      deallocate(extNSnapshotsInSimulation)
      deallocate(extCovFreeEnergies)

    end subroutine computePMF

    subroutine bootstrap(this,nresamples)
      use random_m
      use bin_m
      implicit none
      class (reducedHamiltonian_t) :: this
      integer(kind=4), intent(in) :: nresamples
      real(kind=fp_kind) :: pmf(nbins), pmfSE(nbins)
 
      real(kind=fp_kind), allocatable :: accumulatedWeights(:)
      real(kind=fp_kind), allocatable :: pmfResampled(:,:)
      real(kind=fp_kind), allocatable :: weightsResampled(:,:)
      integer(kind=4), allocatable :: idBinOrigin(:)
      integer(kind=4), allocatable :: idBinResampled(:,:)

      real(kind=fp_kind) :: rand
      integer(kind=4) :: irand

      integer(kind=4) :: IndexB, IndexS, IndexR
      integer(kind=4) :: JndexS

      allocate(accumulatedWeights(this%TotalNumSnapshots))
      allocate(idBinOrigin(this%TotalNumSnapshots))
      allocate(idBinResampled(this%TotalNumSnapshots,nresamples))
      allocate(weightsResampled(this%TotalNumSnapshots,nresamples))
      allocate(pmfResampled(nbins,nresamples))

      forall(IndexS = 1 : this%TotalNumSnapshots) accumulatedWeights(IndexS) = sum(this%weights(1:IndexS))
      accumulatedWeights(this%TotalNumSnapshots) = 1.d0

      forall(IndexB = 1 : nbins) bins(IndexB)%bincenter = binmin + (IndexB-0.5)*binwidth
      forall(IndexS = 1 : this%TotalNumSnapshots) &
         & idBinOrigin(IndexS) = int((this%snapshots(IndexS)%coordinate - binmin )/binwidth) + 1

      pmfResampled = 0.d0
      do IndexR = 1, nresamples
        write(6,'(1X,A,I5)')'Bootstrapping cycle:',IndexR
        do IndexS = 1, this%TotalNumSnapshots
          rand = MyUniformRand()
          irand = rand*this%TotalNumSnapshots + 1
          do JndexS = 1, this%TotalNumSnapshots
            if(accumulatedWeights(JndexS) >= rand)exit
          end do
          idBinResampled(IndexS,IndexR) = idBinOrigin(irand)
          weightsResampled(IndexS,IndexR) = this%weights(irand)
        end do
        do IndexS = 1, this%TotalNumSnapshots
          if(idBinResampled(IndexS,IndexR) >= 1 .and. idBinResampled(IndexS,IndexR) <= nbins)then
            pmfResampled(idBinResampled(IndexS,IndexR),IndexR) = & 
               & pmfResampled(idBinResampled(IndexS,IndexR),IndexR) + weightsResampled(IndexS,IndexR)
          end if
        end do
      end do
      pmfResampled = -log(pmfResampled)
      forall( IndexR = 1 : nresamples) &
        & pmfResampled(:,IndexR) = pmfResampled(:,IndexR) - pmfResampled(1,IndexR)
      do IndexB = 1, nbins
        pmf(IndexB) = sum(pmfResampled(IndexB,:))/nresamples
        pmfSE(IndexB) = sqrt(sum((pmfResampled(IndexB,:)-pmf(IndexB))**2)/nresamples)
      end do
      pmf = pmf - pmf(1)

      write(6,'(1X,A)')'Potential of mean force under the target Hamiltonian from bootstrapping (kcal/mol)'
      do IndexB = 1, nbins
        write(6,'(F8.3,2F9.3)')bins(IndexB)%bincenter, pmf(IndexB)/targetReducedHamiltonian%beta, &
            & pmfSE(IndexB)/targetReducedHamiltonian%beta
      end do
      deallocate(accumulatedWeights)
      deallocate(idBinOrigin)
      deallocate(idBinResampled)
      deallocate(weightsResampled)
      deallocate(pmfResampled)
    end subroutine bootstrap

    function MobleyOverlap(n,n1,weights1,weights2) result(overlap)
      implicit none
      integer(kind=4), intent(in) :: n, n1
      real(kind=fp_kind), intent(in) :: weights1(n), weights2(n)
      real(kind=fp_kind) :: workWeights1(n), workWeights2(n)
      real(kind=fp_kind) :: overlap
      workWeights1(:) = weights1(:)
      workWeights2(:) = weights2(:)
      workWeights1(:) = workWeights1(:) / sum(workWeights1(:))
      workWeights2(:) = workWeights2(:) / sum(workWeights2(:))
      overlap = n1 * sum(workWeights1(:) * workWeights2(:))
    end function MobleyOverlap

    function crossCorrelationBetweenHs(n,weights1,weights2) result(cc)
      implicit none
      integer(kind=4), intent(in) :: n
      real(kind=fp_kind), intent(in) :: weights1(n), weights2(n)
      real(kind=fp_kind) :: workWeights1(n), workWeights2(n)
      real(kind=fp_kind) :: cc
      real(kind=fp_kind) :: cross_correlation
      workWeights1(:) = weights1(:)
      workWeights2(:) = weights2(:)
      workWeights1(:) = workWeights1(:) / sum(workWeights1(:))
      workWeights2(:) = workWeights2(:) / sum(workWeights2(:))
      cc = cross_correlation(n,workWeights1,workWeights2)
    end function crossCorrelationBetweenHs
end module reducedHamiltonian_m
