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

    subroutine computePMF(this, nbins, binmin, binmax, iGaussSmooth)
      use io_m
      use precision_m
      use bin_m
      use MBAR_m
      use simulation_m
      implicit none
      class(reducedHamiltonian_t) :: this
      integer(kind=4), intent(in) :: nbins
      real(kind=fp_kind), intent(in) :: binmin, binmax
      integer(kind=4), intent(in), optional :: iGaussSmooth

      type (binCollection_t) :: coordBins
      real(kind=fp_kind), allocatable :: extWeights(:,:)
      real(kind=fp_kind), allocatable :: extCovFreeEnergies(:,:)
      integer(kind=4), allocatable :: extNSnapshotsInSimulation(:)
      real(kind=fp_kind), allocatable :: weights4smoothing(:)
      integer(kind=4), allocatable :: samplesInThisBin(:)

      real(kind=fp_kind), allocatable :: properties(:)
      real(kind=fp_kind), allocatable :: avgProp(:), stderrProp(:)
      real(kind=fp_kind), allocatable :: weightSum(:)

      integer(kind=4) :: iDumpHistogram
      integer(kind=4) :: iComputeAvg
     
      integer(kind=4) :: IdxMin

      integer(kind=4) :: IndexW, IndexS, IndexB
      integer(kind=4) :: JndexS
     
      call coordBins%initbins(binmin, binmax, nbins)

      allocate(extWeights(totalNumSnapshots,nSimulations+coordBins%nbins))
      allocate(extCovFreeEnergies(nSimulations+nbins,nSimulations+coordBins%nbins))
      allocate(extNSnapshotsInSimulation(nSimulations+coordBins%nbins))
      allocate(weights4smoothing(totalNumSnapshots))
      allocate(samplesInThisBin(totalNumSnapshots))
      extWeights = 0.d0
      forall(IndexW = 1 : nSimulations)extWeights(:, IndexW) = simulatedReducedHamiltonian(IndexW)%weights(:)
      extNSnapshotsInSimulation = 0
      forall(IndexW = 1 : nSimulations)extNSnapshotsInSimulation(IndexW) = simulatedReducedHamiltonian(IndexW)%nOwnSnapshots

! run Gaussian smoothing to the distribution of the energies 
! On the right-hand side of the MBAR working equation, 
! 1/(\sum_{k=1}^K N_k exp(f_k-U_k(x))) works similar to the density-of-states
! we can Gaussian-smooth the distribution of density-of-states to reduce the
! probabilities of some low energy configurations
      if(present(iGaussSmooth) .and. iGaussSmooth > 0) then
        write(6,'(A)')'Perform Gaussian smoothing on the energy distribution in each bin'
        do IndexB = 1, coordBins%nbins
          weights4smoothing = 0.d0
          samplesInThisBin = 0
          do IndexS = 1, totalNumSnapshots
            if(int(( this%snapshots(IndexS)%coordinate - coordBins%binmin )/coordBins%binwidth) + 1 /= IndexB) cycle
            weights4smoothing(IndexS) = this%weights(IndexS)/exp(-this%reducedEnergies(IndexS))
            samplesInThisBin(IndexS) = 1
          end do
          iDumpHistogram = 0
          if(IndexB == 2)then
            iDumpHistogram = 1
          end if
          write(6,'(1X,A,I3,A,I)')'Number of samples in Bin', IndexB,' :', sum(samplesInThisBin)
          call GaussianSmoothing(totalNumSnapshots,samplesInThisBin,this%reducedEnergies,weights4smoothing,iDumpHistogram)
          do IndexS = 1, totalNumSnapshots
            if(int(( this%snapshots(IndexS)%coordinate - coordBins%binmin )/coordBins%binwidth) + 1 /= IndexB) cycle
            this%weights(IndexS) = weights4smoothing(IndexS) * exp(-this%reducedEnergies(IndexS))
          end do
        end do
      end if

! Compute the PMF
      coordBins%bins(:)%reweightingEntropy = 0.d0
      coordBins%bins(:)%nSnapshotsInBin = 0 
      coordBins%bins(:)%sumOfWeightsInBin = 0.d0
      do IndexS = 1, totalNumSnapshots
        IndexB = int(( this%snapshots(IndexS)%coordinate - coordBins%binmin )/coordBins%binwidth) + 1
        if(IndexB > coordBins%nbins .or. IndexB < 1) cycle
        extWeights(IndexS,nSimulations+IndexB) = this%weights(IndexS)
        coordBins%bins(IndexB)%reweightingEntropy = coordBins%bins(IndexB)%reweightingEntropy + &
          & this%weights(IndexS)*log(min(this%weights(IndexS)+1.D-300, 1.00d0))
        coordBins%bins(IndexB)%nSnapshotsInBin = coordBins%bins(IndexB)%nSnapshotsInBin + 1
        coordBins%bins(IndexB)%sumOfWeightsInBin = coordBins%bins(IndexB)%sumOfWeightsInBin + &
          & this%weights(IndexS)
      end do
      coordBins%bins(:)%pmf = coordBins%bins(:)%sumOfWeightsInBin / coordBins%bins(:)%binwidth
  
! Compute the variances of the PMF
      forall(IndexB = 1 : coordBins%nbins) extWeights(:,nSimulations+IndexB) = &
            & extWeights(:,nSimulations+IndexB)/sum(extWeights(:,nSimulations+IndexB))
  
      call ComputCovMatFromWeights(totalNumSnapshots,nSimulations+coordBins%nbins,extNSnapshotsInSimulation,extWeights, &
            & extCovFreeEnergies)
  
      coordBins%bins(:)%pmf = -log(coordBins%bins(:)%pmf)

      if(present(iGaussSmooth) .and. iGaussSmooth > 0) then
!       IdxMin = minloc(coordBins%bins(:)%pmf,1)
!       coordBins%bins(:)%pmf = coordBins%bins(:)%pmf - coordBins%bins(IdxMin)%pmf
!       forall(IndexB = 1 : coordBins%nbins) coordBins%bins(IndexB)%pmfSE = &
!           & sqrt(  extCovFreeEnergies(nSimulations+IndexB,nSimulations+IndexB) &
!               &  + extCovFreeEnergies(nSimulations+IdxMin,nSimulations+IdxMin) &
!               &  - 2*extCovFreeEnergies(nSimulations+IdxMin,nSimulations+IndexB) )

!       coordBins%bins(:)%pmf = coordBins%bins(:)%pmf - coordBins%bins(coordBins%nbins)%pmf
!       forall(IndexB = 1 : coordBins%nbins) coordBins%bins(IndexB)%pmfSE = &
!           & sqrt(  extCovFreeEnergies(nSimulations+IndexB,nSimulations+IndexB) &
!               &  + extCovFreeEnergies(nSimulations+coordBins%nbins,nSimulations+coordBins%nbins) &
!               &  - 2*extCovFreeEnergies(nSimulations+coordBins%nbins,nSimulations+IndexB) )

        coordBins%bins(:)%pmf = coordBins%bins(:)%pmf - coordBins%bins(1)%pmf
        forall(IndexB = 1 : coordBins%nbins) coordBins%bins(IndexB)%pmfSE = &
           & sqrt(  extCovFreeEnergies(nSimulations+IndexB,nSimulations+IndexB) &
               &  + extCovFreeEnergies(nSimulations+1,nSimulations+1) &
               &  - 2*extCovFreeEnergies(nSimulations+1,nSimulations+IndexB) )
      else
        coordBins%bins(:)%pmf = coordBins%bins(:)%pmf - coordBins%bins(1)%pmf
        forall(IndexB = 1 : coordBins%nbins) coordBins%bins(IndexB)%pmfSE = &
            & sqrt(  extCovFreeEnergies(nSimulations+IndexB,nSimulations+IndexB) &
                &  + extCovFreeEnergies(nSimulations+1,nSimulations+1) &
              &    - 2*extCovFreeEnergies(nSimulations+1,nSimulations+IndexB) )
      end if
      coordBins%bins(:)%reweightingEntropy = -(coordBins%bins(:)%reweightingEntropy/coordBins%bins(:)%sumOfWeightsInBin &
            & -log(coordBins%bins(:)%sumOfWeightsInBin)) /log(dble(coordBins%bins(:)%nSnapshotsInBin))

      open(id_target_pmf_file , file = targetHamiltonianPmfFile)
      write(6,'(1X,A)')'Potential of mean force under the target Hamiltonian (kcal/mol)'
      write(6,'(1X,A)')'    RC       f        df        RE'
      write(id_target_pmf_file,'(1X,A)')'#   RC       f        df        RE'
      do IndexB = 1, coordBins%nbins
        write(id_target_pmf_file,'(F8.3,3F9.3)')coordBins%bins(IndexB)%bincenter, coordBins%bins(IndexB)%pmf/this%beta, &
             & coordBins%bins(IndexB)%pmfSE/this%beta, coordBins%bins(IndexB)%reweightingEntropy
        write(6,'(F8.3,3F9.3)')coordBins%bins(IndexB)%bincenter, coordBins%bins(IndexB)%pmf/this%beta, &
             & coordBins%bins(IndexB)%pmfSE/this%beta, coordBins%bins(IndexB)%reweightingEntropy
      end do
      close(id_target_pmf_file)

! Compute ensemble average
      write(6,'(1X,A)')'Compute ensemble average or nor? 0. No. 1. Yes'
      read*,iComputeAvg
      if(iComputeAvg==1)then
        allocate(properties(totalNumSnapshots))
        allocate(avgProp(coordBins%nbins))
        allocate(stderrProp(coordBins%nbins))
        do IndexS = 1, totalNumSnapshots
          read(104,*)properties(IndexS)
        end do
        avgProp = 0.d0
        do IndexS = 1, totalNumSnapshots
          IndexB = int(( this%snapshots(IndexS)%coordinate - coordBins%binmin )/coordBins%binwidth) + 1
          if(IndexB > coordBins%nbins .or. IndexB < 1) cycle
          avgProp(IndexB) = avgProp(IndexB) + this%weights(IndexS)*properties(IndexS)
        end do
        avgProp = avgProp / max(coordBins%bins(:)%sumOfWeightsInBin,1.D-300)

        stderrProp = 0.d0
        do IndexS = 1, totalNumSnapshots
          IndexB = int(( this%snapshots(IndexS)%coordinate - coordBins%binmin )/coordBins%binwidth) + 1
          if(IndexB > coordBins%nbins .or. IndexB < 1) cycle
          stderrProp(IndexB) = stderrProp(IndexB) + this%weights(IndexS)*(properties(IndexS)-avgProp(IndexB))**2
        end do
        stderrProp(:) = sqrt(stderrProp(:)/coordBins%bins(:)%sumOfWeightsInBin)
        do IndexB = 1, nbins
          write(105,'(3F10.5)')coordBins%bins(IndexB)%bincenter, avgProp(IndexB), stderrProp(IndexB)
        end do
        deallocate(properties)
        deallocate(avgProp)
        deallocate(stderrProp)
      end if
      deallocate(extWeights)
      deallocate(extNSnapshotsInSimulation)
      deallocate(extCovFreeEnergies)
      deallocate(weights4smoothing)
      deallocate(samplesInThisBin)
      call coordBins%deleteBinInfo()
    end subroutine computePMF

    subroutine bootstrap(this, nbins, binmin, binmax, nresamples)
      use random_m
      use bin_m
      implicit none
      class (reducedHamiltonian_t) :: this
      integer(kind=4), intent(in) :: nbins
      real(kind=fp_kind), intent(in) :: binmin, binmax
      integer(kind=4), intent(in) :: nresamples
      real(kind=fp_kind) :: pmf(nbins), pmfSE(nbins)
 
      type (binCollection_t) :: coordBins
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

      call coordBins%initbins(binmin, binmax, nbins)
!      forall(IndexB = 1 : nbins) coordBins%bins(IndexB)%bincenter = binmin + (IndexB-0.5)*binwidth
      forall(IndexS = 1 : this%TotalNumSnapshots) &
         & idBinOrigin(IndexS) = int((this%snapshots(IndexS)%coordinate - coordBins%binmin )/coordBins%binwidth) + 1

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
          if(idBinResampled(IndexS,IndexR) >= 1 .and. idBinResampled(IndexS,IndexR) <= coordBins%nbins)then
            pmfResampled(idBinResampled(IndexS,IndexR),IndexR) = & 
               & pmfResampled(idBinResampled(IndexS,IndexR),IndexR) + weightsResampled(IndexS,IndexR)
          end if
        end do
      end do
      pmfResampled = -log(pmfResampled)
      forall( IndexR = 1 : nresamples) &
        & pmfResampled(:,IndexR) = pmfResampled(:,IndexR) - pmfResampled(1,IndexR)
      do IndexB = 1, coordBins%nbins
        pmf(IndexB) = sum(pmfResampled(IndexB,:))/nresamples
        pmfSE(IndexB) = sqrt(sum((pmfResampled(IndexB,:)-pmf(IndexB))**2)/nresamples)
      end do
      pmf = pmf - pmf(1)
!      pmf = pmf - minval(pmf(:))
!      pmf = pmf - pmf(coordBins%nbins)

      write(6,'(1X,A)')'Potential of mean force under the target Hamiltonian from bootstrapping (kcal/mol)'
      write(6,'(1X,A)')'    RC       f        df'
      do IndexB = 1, coordBins%nbins
        write(6,'(F8.3,2F9.3)')coordBins%bins(IndexB)%bincenter, pmf(IndexB)/targetReducedHamiltonian%beta, &
            & pmfSE(IndexB)/targetReducedHamiltonian%beta
      end do
      deallocate(accumulatedWeights)
      deallocate(idBinOrigin)
      deallocate(idBinResampled)
      deallocate(weightsResampled)
      deallocate(pmfResampled)
      call coordBins%deleteBinInfo()
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

    subroutine GaussianSmoothing(n,sampleInBin,energies,weights,iDumpHistogram)
      use bin_m
      implicit none
      integer(kind=4), intent(in) :: n
      integer(kind=4), intent(in) :: sampleInBin(n)
      real(kind=fp_kind), intent(in) :: energies(n)
      real(kind=fp_kind), intent(in out) :: weights(n)
      integer(kind=4), intent(in) :: iDumpHistogram

      type (binCollection_t) :: energyBins

      real(kind=fp_kind), allocatable :: energiesInBin(:)
      integer(kind=4), allocatable :: binID(:)
      real(kind=fp_kind) :: minE, maxE, width
      real(kind=fp_kind) :: meanE, sigmaE
      real(kind=fp_kind), allocatable :: scaleFactor(:)
      real(kind=fp_kind), allocatable :: gaussianCumulative(:)
      integer(kind=4) :: IndexB, IndexS
      integer(kind=4) :: JndexB, JndexS
! initialize bins 
      allocate(energiesInBin(sum(sampleInBin)))
      JndexS = 0
      do IndexS = 1, n
        if(sampleInBin(IndexS) == 1)then
          JndexS = JndexS +1
          energiesInBin(JndexS) = energies(IndexS)
        end if
      end do

!      minE = minval(energies(:))-1.0E-10
!      maxE = maxval(energies(:))+1.0E-10
      minE = minval(energiesInBin(:))-1.0E-10
      maxE = maxval(energiesInBin(:))+1.0E-10
      width = 0.2d0
      call energyBins%initbins2(minE, maxE, width)
      allocate(binID(n))
      allocate(scaleFactor(energyBins%nbins))
      allocate(gaussianCumulative(energyBins%nbins))

! compute histogram from fractional weights
      energyBins%bins(:)%sumOfWeightsInBin = 0.d0
      do IndexS = 1, n
        if(sampleInBin(IndexS) /= 1)cycle
        binID(IndexS) = int(( energies(IndexS) - energyBins%binmin )/energyBins%binwidth) + 1
        energyBins%bins(binID(IndexS))%sumOfWeightsInBin = energyBins%bins(binID(IndexS))%sumOfWeightsInBin &
           & + weights(IndexS)
      end do
! normalization, because cumulative_gaussian is normalized
      energyBins%bins(:)%sumOfWeightsInBin = &
         & energyBins%bins(:)%sumOfWeightsInBin / sum(energyBins%bins(:)%sumOfWeightsInBin)

      call Weighted_Mean_and_StandardDev(n, weights, energies, meanE, sigmaE)
      write(*,'(A,1X,E20.10,1X,A,1X,E20.10)')'Fitted Gaussian mean=', meanE, 'sigma=', sigmaE

      forall(IndexB = 1:energyBins%nbins) gaussianCumulative(IndexB) = &
         &   cumulative_gaussian(meanE, sigmaE, & 
                               & energyBins%bins(IndexB)%bincenter + energyBins%bins(IndexB)%binwidth/2) &
         & - cumulative_gaussian(meanE, sigmaE, & 
                               & energyBins%bins(IndexB)%bincenter - energyBins%bins(IndexB)%binwidth/2)

      if(iDumpHistogram==1) then
         do IndexB = 1, energyBins%nbins
            write(101,'(F10.6,4E20.10)')energyBins%bins(IndexB)%bincenter, &
               & energyBins%bins(IndexB)%sumOfWeightsInBin, &
               & gaussianCumulative(IndexB),&
!               & energyBins%bins(IndexB)%sumOfWeightsInBin * exp(-energyBins%bins(IndexB)%bincenter), &
               & gaussianCumulative(IndexB) * exp(-energyBins%bins(IndexB)%bincenter)
         end do
         do IndexS = 1, n
           if(sampleInBin(IndexS) /= 1)cycle
           write(100,*)IndexS,energies(IndexS)
         end do
      end if

      scaleFactor = 1.d0
      do IndexB = 1, energyBins%nbins
        if(energyBins%bins(IndexB)%sumOfWeightsInBin > 0.d0) then
          scaleFactor(IndexB) = gaussianCumulative(IndexB) / energyBins%bins(IndexB)%sumOfWeightsInBin
        end if
      end do

      do IndexS = 1, n
        if(sampleInBin(IndexS) /= 1)cycle
        weights(IndexS) = weights(IndexS) * scaleFactor(binID(IndexS))
      end do

      call energyBins%deleteBinInfo()
      deallocate(binID)
      deallocate(scaleFactor)
      deallocate(gaussianCumulative)
      deallocate(energiesInBin)
      contains
        function gaussian(x0,sigma,x)
          use precision_m
          implicit none
          real(kind=fp_kind) :: gaussian
          real(kind=fp_kind), intent(in) :: x0,sigma,x
          real(kind=fp_kind) :: pi
          pi = atan(1.0d0)*4.0d0
          gaussian = 1.d0/(sqrt(2*pi)*sigma)*exp(-(x-x0)**2/(2*sigma**2))
        end function gaussian

        pure function cumulative_gaussian(x0,sigma,x)
          implicit none
          real(kind=fp_kind) :: cumulative_gaussian
          real(kind=fp_kind), intent(in) :: x0,sigma,x
          real(kind=fp_kind) :: pi
          cumulative_gaussian = 0.5d0*(1+erf((x-x0)/(sqrt(2.d0)*sigma)))
        end function cumulative_gaussian
    end subroutine GaussianSmoothing
end module reducedHamiltonian_m
