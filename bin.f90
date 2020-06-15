module bin_m
  use precision_m
  implicit none
  private
  type :: bin_t
    real(kind=fp_kind) :: bincenter
    real(kind=fp_kind) :: binwidth
    integer(kind=4) :: nSnapshotsInBin
    real(kind=fp_kind) :: sumOfWeightsInBin
    real(kind=fp_kind) :: reweightingEntropy
    real(kind=fp_kind) :: pmf
    real(kind=fp_kind)  :: pmfSE
  end type bin_t

  type, public :: binCollection_t
    integer(kind=4) :: nbins
    real(kind=fp_kind) :: binmin, binmax, binwidth
    type (bin_t), allocatable :: bins(:)
    contains
      procedure :: initbins
      procedure :: initbins2
      procedure :: deleteBinInfo
  end type binCollection_t

contains
    subroutine initbins(this, binmin, binmax, nbins)
      use precision_m
      implicit none
      class(binCollection_t) :: this   
      real(kind=fp_kind), intent(in) :: binmin, binmax
      integer(kind=4), intent(in) :: nbins
      integer(kind=4) :: IndexB
      this%binmin = binmin
      this%binmax = binmax
      this%nbins = nbins
      this%binwidth = (this%binmax - this%binmin)/this%nbins
      allocate(this%bins(nbins))
      forall(IndexB = 1:this%nbins) 
         this%bins(IndexB)%bincenter = this%binmin + (IndexB-0.5)*this%binwidth
         this%bins(IndexB)%binwidth = this%binwidth
         this%bins(IndexB)%nSnapshotsInBin = 0
         this%bins(IndexB)%sumOfWeightsInBin = 0.d0
         this%bins(IndexB)%reweightingEntropy = 0.d0
         this%bins(IndexB)%pmf = 0.d0
         this%bins(IndexB)%pmfSE = 0.d0
      end forall
    end subroutine initbins

    subroutine initbins2(this, binmin, binmax, binwidth)
      use precision_m
      implicit none
      class(binCollection_t) :: this
      real(kind=fp_kind), intent(in) :: binmin, binmax, binwidth
      integer(kind=4) :: IndexB
      this%binmin = binmin
      this%binmax = binmax
      this%binwidth = binwidth
      this%nbins = int( (this%binmax - this%binmin)/this%binwidth ) + 1
      allocate(this%bins(this%nbins))
      forall(IndexB = 1:this%nbins)
         this%bins(IndexB)%bincenter = this%binmin + (IndexB-0.5)*this%binwidth
         this%bins(IndexB)%binwidth = this%binwidth
         this%bins(IndexB)%nSnapshotsInBin = 0
         this%bins(IndexB)%sumOfWeightsInBin = 0.d0
         this%bins(IndexB)%reweightingEntropy = 0.d0
         this%bins(IndexB)%pmf = 0.d0
         this%bins(IndexB)%pmfSE = 0.d0
      end forall
    end subroutine initbins2

    subroutine deleteBinInfo(this)
      class(binCollection_t) :: this
      if(allocated(this%bins))deallocate(this%bins)
    end subroutine deleteBinInfo

end module bin_m
