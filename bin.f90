module bin_m
  use precision_m
  implicit none
  private
  integer(kind=4), public :: nbins
  real(kind=fp_kind), public :: binmin, binmax, binwidth
  type, public :: bin_t
    real(kind=fp_kind) :: bincenter
    real(kind=fp_kind) :: binwidth
    integer(kind=4) :: nSnapshotsInBin
    real(kind=fp_kind) :: sumOfWeightsInBin
    real(kind=fp_kind) :: reweightingEntropy
    real(kind=fp_kind) :: pmf
    real(kind=fp_kind)  :: pmfSE
    contains

  end type bin_t
  type (bin_t), allocatable, public :: bins(:)
  public :: deleteBinInfo
  public :: initbins
  public :: computePMF
contains
    subroutine initbins()
      implicit none
      integer(kind=4) :: IndexB
      binwidth = (binmax - binmin)/nbins
      allocate(bins(nbins))
      forall(IndexB = 1:nbins) 
         bins(IndexB)%bincenter = binmin + (IndexB-0.5)*binwidth
         bins(IndexB)%binwidth = binwidth
         bins(IndexB)%nSnapshotsInBin = 0
         bins(IndexB)%sumOfWeightsInBin = 0.d0
         bins(IndexB)%reweightingEntropy = 0.d0
         bins(IndexB)%pmf = 0.d0
         bins(IndexB)%pmfSE = 0.d0
      end forall
    end subroutine initbins

    subroutine deleteBinInfo()
      if(allocated(bins))deallocate(bins)
    end subroutine deleteBinInfo

    subroutine computePMF()
      use io_m
      use precision_m
      use MBAR_m
      use simulation_m
      use reducedHamiltonian_m
      implicit none

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
        IndexB = int(( targetReducedHamiltonian%snapshots(IndexS)%coordinate - binmin )/binwidth) + 1
        if(IndexB > nbins .or. IndexB < 1) cycle
        extWeights(IndexS,nSimulations+IndexB) = targetReducedHamiltonian%weights(IndexS)
        bins(IndexB)%pmf = bins(IndexB)%pmf + targetReducedHamiltonian%weights(IndexS)
        bins(IndexB)%reweightingEntropy = bins(IndexB)%reweightingEntropy + &
          & targetReducedHamiltonian%weights(IndexS)*log(targetReducedHamiltonian%weights(IndexS))
        bins(IndexB)%nSnapshotsInBin = bins(IndexB)%nSnapshotsInBin + 1
        bins(IndexB)%sumOfWeightsInBin = bins(IndexB)%sumOfWeightsInBin + &
          & targetReducedHamiltonian%weights(IndexS)
      end do
  
      forall(IndexB = 1 : nbins) extWeights(:,nSimulations+IndexB) = extWeights(:,nSimulations+IndexB)/sum(extWeights(:,nSimulations+IndexB))
  
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
        write(id_target_pmf_file,'(F8.3,3F9.3)')bins(IndexB)%bincenter, bins(IndexB)%pmf/targetReducedHamiltonian%beta, &
             & bins(IndexB)%pmfSE/targetReducedHamiltonian%beta, bins(IndexB)%reweightingEntropy
        write(6,'(F8.3,3F9.3)')bins(IndexB)%bincenter, bins(IndexB)%pmf/targetReducedHamiltonian%beta, &
             & bins(IndexB)%pmfSE/targetReducedHamiltonian%beta, bins(IndexB)%reweightingEntropy
      end do
      close(id_target_pmf_file)
      deallocate(extWeights)
      deallocate(extNSnapshotsInSimulation)
      deallocate(extCovFreeEnergies)

    end subroutine computePMF

end module bin_m
