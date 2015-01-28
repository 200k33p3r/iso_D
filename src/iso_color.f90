module iso_color

  !MESA modules
  use const_def, only: sp
  use utils_lib, only: alloc_iounit, free_iounit

  !other modules
  use color_def
  use color_lib

  !local modules
  use iso_eep_support

  implicit none

  type(bc_table), allocatable :: bc(:)
  character(len=file_path) :: color_input = 'input.cmd', bc_table_list, color_suffix

contains

  subroutine read_color_input(ierr)
    integer, intent(out) :: ierr
    integer :: io
    io=alloc_iounit(ierr)
    if(ierr/=0) return
    open(io,file=trim(color_input),action='read',status='old',iostat=ierr)
    if(ierr/=0) return
    read(io,'(a)') bc_table_list
    read(io,'(a)') color_suffix
    close(io)
    call free_iounit(io)
    call color_init(bc_table_list,bc,ierr)
  end subroutine read_color_input

  subroutine get_mags(iso,iT,ig,iL)
    type(isochrone), intent(inout) :: iso
    integer :: i, ierr
    real(sp), allocatable :: res(:)
    integer :: iT, ig, iL
    real(sp) :: logT, logg, logL
    iso% nfil = bc(1)% num_filter
    allocate(iso% mags(iso% nfil, iso% neep),res(iso% nfil))
    allocate(iso% labels(iso% nfil))
    iso% labels = bc(1)% labels
    res = 0.
    iso% mags = 0.0         
    do i=1,iso% neep
       logT = real(iso% data(iT,i),kind=sp)
       logg = real(iso% data(ig ,i),kind=sp)
       logL = real(iso% data(iL ,i),kind=sp)
       call color_get(bc(1), logT, logg, 1, 1, res, ierr)
       iso% mags(:,i) = SolBol - 2.5*logL - res
    enddo
    deallocate(res)
  end subroutine get_mags

end module iso_color