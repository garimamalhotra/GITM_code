! Reads in a netcdf file with lat-lon distribution of Kzz

subroutine read_eddy_map

  use netcdf
  use ModSizeGITM, only: nLons, nLats, nBlocks
  use ModGITM,     only: Longitude, Latitude
  use ModInputs,   only: EddyMapFile, EddyDiffusionCoef2D

  implicit none

  real, parameter :: rad2deg = 57.2957795130823209
  integer :: ncid, varid, did_lat, did_lon
  integer :: nLatFile, nLonFile
  integer :: iBlock, iLon, iLat
  real    :: lonDeg, latDeg
  real, allocatable :: latF(:), lonF(:), KzzF(:, :)

  call nc_check(nf90_open(trim(EddyMapFile), nf90_nowrite, ncid), &
                "open " // trim(EddyMapFile))

  call nc_check(nf90_inq_dimid(ncid, "lat", did_lat), "inq dimid lat")
  call nc_check(nf90_inquire_dimension(ncid, did_lat, len=nLatFile), &
                "inq dim lat")
  call nc_check(nf90_inq_dimid(ncid, "lon", did_lon), "inq dimid lon")
  call nc_check(nf90_inquire_dimension(ncid, did_lon, len=nLonFile), &
                "inq dim lon")

  allocate(latF(nLatFile), lonF(nLonFile), KzzF(nLonFile, nLatFile))

  call nc_check(nf90_inq_varid(ncid, "lat", varid),    "inq varid lat")
  call nc_check(nf90_get_var  (ncid, varid, latF),     "get var lat")
  call nc_check(nf90_inq_varid(ncid, "lon", varid),    "inq varid lon")
  call nc_check(nf90_get_var  (ncid, varid, lonF),     "get var lon")
  call nc_check(nf90_inq_varid(ncid, "Kzz_P0", varid), "inq varid Kzz_P0")
  call nc_check(nf90_get_var  (ncid, varid, KzzF),     "get var Kzz_P0")

  call nc_check(nf90_close(ncid), "close")

  do iBlock = 1, nBlocks
    do iLat = 1, nLats
      do iLon = 1, nLons
        lonDeg = Longitude(iLon, iBlock) * rad2deg
        latDeg = Latitude (iLat, iBlock) * rad2deg
        EddyDiffusionCoef2D(iLon, iLat, iBlock) = &
          bilinear(lonF, latF, KzzF, lonDeg, latDeg)
      end do
    end do
  end do

  deallocate(latF, lonF, KzzF)

contains

  subroutine nc_check(istat, where)
    integer,          intent(in) :: istat
    character(len=*), intent(in) :: where
    if (istat /= nf90_noerr) then
      write(*, *) "read_eddy_map: NetCDF error at ", trim(where), ": ", &
                  trim(nf90_strerror(istat))
      call stop_gitm("read_eddy_map failed")
    end if
  end subroutine nc_check

  function bilinear(lonF, latF, KzzF, lonDeg, latDeg) result(val)
    real, intent(in) :: lonF(:), latF(:), KzzF(:, :), lonDeg, latDeg
    real :: val, lonW, latW, lonQ, lonSpan
    integer :: nLonF, nLatF, i0, i1, j0, j1

    nLonF = size(lonF)
    nLatF = size(latF)

    lonSpan = 360.0
    lonQ = lonDeg
    do while (lonQ <  lonF(1));            lonQ = lonQ + lonSpan; end do
    do while (lonQ >= lonF(1) + lonSpan);  lonQ = lonQ - lonSpan; end do

    call locate(lonF, lonQ, i0, i1, lonW, periodic=.true., span=lonSpan)
    call locate(latF, latDeg, j0, j1, latW, periodic=.false., span=0.0)

    val = (1.0 - lonW) * (1.0 - latW) * KzzF(i0, j0) &
        +        lonW  * (1.0 - latW) * KzzF(i1, j0) &
        + (1.0 - lonW) *        latW  * KzzF(i0, j1) &
        +        lonW  *        latW  * KzzF(i1, j1)
  end function bilinear

  subroutine locate(x, q, i0, i1, w, periodic, span)
    real,    intent(in)  :: x(:), q, span
    logical, intent(in)  :: periodic
    integer, intent(out) :: i0, i1
    real,    intent(out) :: w
    integer :: n, lo, hi, mid

    n = size(x)
    if (q <= x(1)) then
      if (periodic) then
        i0 = n;     i1 = 1
        w  = (q - (x(n) - span)) / (x(1) - (x(n) - span))
      else
        i0 = 1;     i1 = 1;     w = 0.0
      end if
      return
    else if (q >= x(n)) then
      if (periodic) then
        i0 = n;     i1 = 1
        w  = (q - x(n)) / ((x(1) + span) - x(n))
      else
        i0 = n;     i1 = n;     w = 0.0
      end if
      return
    end if

    lo = 1;  hi = n
    do while (hi - lo > 1)
      mid = (lo + hi) / 2
      if (q >= x(mid)) then
        lo = mid
      else
        hi = mid
      end if
    end do
    i0 = lo;  i1 = hi
    w  = (q - x(i0)) / (x(i1) - x(i0))
  end subroutine locate

end subroutine read_eddy_map
