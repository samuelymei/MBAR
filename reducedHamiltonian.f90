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
    real(kind=fp_kind) :: freeEnergy
    type (simulation_t), pointer :: ownSimulation 
    contains
      procedure :: destroy
      procedure :: init
      procedure :: processTrajectories
      procedure :: readInEnergy
      procedure :: bootstrap
  end type reducedHamiltonian_t

  ! type (reducedHamiltonian_t), allocatable, public :: simulatedReducedHamiltonian(:)


  public :: reducedHamiltonian_t
  
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

    subroutine bootstrap(this,nresamples,nbins,binmin,binwidth,pmf,pmfSE)
      use random_m
      implicit none
      class (reducedHamiltonian_t) :: this
      integer(kind=4), intent(in) :: nresamples
      integer(kind=4), intent(in) :: nbins
      real(kind=fp_kind), intent(in) :: binmin, binwidth
      real(kind=fp_kind), intent(out) :: pmf(nbins), pmfSE(nbins)
 
      real(kind=fp_kind), allocatable :: accumulatedWeights(:)
      real(kind=fp_kind), allocatable :: pmfResampled(:,:)
      real(kind=fp_kind), allocatable :: weightsResampled(:,:)
      integer(kind=4), allocatable :: idBinOrigin(:)
      integer(kind=4), allocatable :: idBinResampled(:,:)

      real(kind=fp_kind) :: rand
      integer(kind=4) :: irand

      real(kind=fp_kind), allocatable :: bincenters(:)

      integer(kind=4) :: IndexB, IndexS, IndexR
      integer(kind=4) :: JndexS

      allocate(accumulatedWeights(this%TotalNumSnapshots))
      allocate(idBinOrigin(this%TotalNumSnapshots))
      allocate(bincenters(nbins))
      allocate(idBinResampled(this%TotalNumSnapshots,nresamples))
      allocate(weightsResampled(this%TotalNumSnapshots,nresamples))
      allocate(pmfResampled(nbins,nresamples))

      accumulatedWeights(1) = this%weights(1)
      do IndexS = 2, this%TotalNumSnapshots
        accumulatedWeights(IndexS) = accumulatedWeights(IndexS-1) + this%weights(IndexS)
      end do
      accumulatedWeights(this%TotalNumSnapshots) = 1.d0

      do IndexB = 1, nbins
        bincenters(IndexB) = binmin + (IndexB-0.5)*binwidth
      end do 
      do IndexS = 1, this%TotalNumSnapshots
        idBinOrigin(IndexS) = int(( this%snapshots(IndexS)%coordinate - binmin )/binwidth) + 1
      end do

      pmfResampled = 0.d0
      do IndexR = 1, nresamples
        write(6,'(1X,A,I)')'Bootstrapping cycle:',IndexR
        do IndexS = 1, this%TotalNumSnapshots
          rand = MyUniformRand()
          irand = rand*this%TotalNumSnapshots + 1
          do JndexS = 1, this%TotalNumSnapshots
            if(accumulatedWeights(JndexS) >= rand)exit
          end do
          idBinResampled(IndexS,IndexR)=idBinOrigin(irand)
          weightsResampled(IndexS,IndexR)=this%weights(irand)
        end do
        weightsResampled(:,IndexR) = weightsResampled(:,IndexR) / sum(weightsResampled(:,IndexR))
        do IndexS = 1, this%TotalNumSnapshots
          if(idBinResampled(IndexS,IndexR) >= 1 .and. idBinResampled(IndexS,IndexR) <= nbins)then
            pmfResampled(idBinResampled(IndexS,IndexR),IndexR) = & 
               & pmfResampled(idBinResampled(IndexS,IndexR),IndexR) + weightsResampled(IndexS,IndexR)
               
          end if
        end do
      end do
      pmfResampled = -log(pmfResampled)/this%beta
      do IndexR = 1, nresamples
        pmfResampled(:,IndexR) = pmfResampled(:,IndexR) - minval(pmfResampled(:,IndexR))
      end do
      do IndexB = 1, nbins
        pmf(IndexB) = sum(pmfResampled(IndexB,:))/nresamples
        pmfSE(IndexB) = sqrt(sum((pmfResampled(IndexB,:)-pmf(IndexB))**2)/nresamples)
      end do
      pmf = pmf - minval(pmf)

      deallocate(accumulatedWeights)
      deallocate(idBinOrigin)
      deallocate(bincenters)
      deallocate(idBinResampled)
      deallocate(weightsResampled)
      deallocate(pmfResampled)
    end subroutine bootstrap

end module reducedHamiltonian_m
