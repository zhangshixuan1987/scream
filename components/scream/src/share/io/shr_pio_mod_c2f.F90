module shr_pio_mod_c2f
  use iso_c_binding, only : c_int

  ! This modules contains a few C-interoperable routines, which allow C code
  ! to query global PIO information from the shr_pio_mod. This is needed b/c
  ! E3SM inits a pio subsystem for the atmosphere, and we need to get the
  ! specs of that subsystem in order to be able to access it.

  public :: get_io_sys, get_io_type, get_io_rearranger, get_io_format
contains
  function get_io_sys (atm_id) result(sys_id) bind(c)
    use shr_pio_mod, only : shr_pio_getiosys
    integer (kind=c_int), intent(in), value :: atm_id

    integer (kind=c_int) :: sys_id

    sys_id = shr_pio_getiosys(atm_id)%iosysid
  end function get_io_sys
  function get_io_type (atm_id) result(io_type) bind(c)
    use shr_pio_mod, only : shr_pio_getiotype
    integer (kind=c_int), intent(in), value :: atm_id

    integer (kind=c_int) :: io_type

    io_type= shr_pio_getiotype(atm_id)
  end function get_io_type
  function get_io_rearranger (atm_id) result(rearr_id) bind(c)
    use shr_pio_mod, only : shr_pio_getiorearranger
    integer (kind=c_int), intent(in), value :: atm_id

    integer (kind=c_int) :: rearr_id

    rearr_id = shr_pio_getiorearranger(atm_id)
  end function get_io_rearranger
  function get_io_format (atm_id) result(format_id) bind(c)
    use shr_pio_mod, only : shr_pio_getioformat
    integer (kind=c_int), intent(in), value :: atm_id

    integer (kind=c_int) :: format_id

    format_id = shr_pio_getioformat(atm_id)
  end function get_io_format
end module shr_pio_mod_c2f
