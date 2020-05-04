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
    integer(kind=4) :: TotalNumSnapshots
    integer(kind=4) :: nOwnSnapshots
    type (snapshot_t), allocatable :: snapshots(:)
    real(kind=fp_kind), allocatable :: reducedEnergies(:)
    real(kind=fp_kind), allocatable :: weights(:)
    real(kind=fp_kind) :: freeEnergy
    type (simulation_t), pointer :: ownSimulation 
    contains
      procedure :: destroy
      procedure :: init
      procedure :: processTrajectories
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
    end subroutine init

    subroutine processTrajectories(this)
      implicit none
      class(reducedHamiltonian_t) :: this
      integer(kind=4) :: indexW, indexS
      integer(kind=4) :: jndexS
      if(idebug == 1) write(*,*) 'Entering reducedHamiltonian%processTrajectories'
      jndexS = 0
      do indexW = 1, nSimulations
        do indexS = 1, simulations(indexW)%nSnapshots
          jndexS = jndexS + 1
          this%snapshots(jndexS)%coordinate = &
               & simulations(indexW)%snapshots(indexS)%coordinate
          this%reducedEnergies(jndexS) = &
               & simulations(indexW)%snapshots(indexS)%energyUnbiased &
               & + 0.5*this%ownSimulation%restraintStrength &
               & * (this%snapshots(jndexS)%coordinate - this%ownSimulation%restraintCenter)**2
          this%reducedEnergies(jndexS) = this%reducedEnergies(jndexS) * this%beta
        end do
      end do
    end subroutine processTrajectories

    subroutine destroy(this)
      class(reducedHamiltonian_t) :: this
      if(idebug == 1) write(*,*) 'Entering reducedHamiltonian%destroy'
      if(allocated(this%snapshots))deallocate(this%snapshots)
      if(allocated(this%reducedEnergies))deallocate(this%reducedEnergies)
      if(allocated(this%weights))deallocate(this%weights)
    end subroutine destroy
end module reducedHamiltonian_m
