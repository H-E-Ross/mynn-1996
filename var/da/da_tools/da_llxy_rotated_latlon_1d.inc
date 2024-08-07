subroutine da_llxy_rotated_latlon_1d(lat,lon, proj, x, y)                          

   !----------------------------------------------------------------------- 
   ! Purpose: Compute the x/y location of a lat/lon on a rotated LATLON grid.
   ! Author :  Syed RH Rizvi,     MMM/NCAR
   !           06/01/2008
   !---------------------------------------------------------------------------

   implicit none

   real, intent(in)             :: lat(:)
   real, intent(in)             :: lon(:)
   type(proj_info), intent(in)  :: proj
   real, intent(out)            :: x(:)
   real, intent(out)            :: y(:)

   real, allocatable            :: rot_lat(:), rot_lon(:), deltalat(:), deltalon(:), lon360(:)
   real, allocatable            :: xlat(:), xlon(:)
   real                         :: cen_lat, cen_lon, latinc, loninc 
   integer                      :: n

   if (trace_use_frequent) call da_trace_entry("da_llxy_rotated_latlon_1d")

   n = size(lat)
   allocate ( rot_lat(n), rot_lon(n), deltalat(n), deltalon(n), lon360(n), xlat(n), xlon(n) )

   ! To account for issues around the dateline, convert the incoming
   ! longitudes to be 0->360.0
   where (lon < 0)
      lon360 = lon + 360.0 
   elsewhere
      lon360 = lon
   end where

   xlat = deg_to_rad*lat
   xlon = deg_to_rad*lon360
   cen_lat = deg_to_rad*proj%lat1
   cen_lon = deg_to_rad*proj%lon1
   if (cen_lon < 0.) cen_lon = cen_lon + 360.
  
   latinc = proj%latinc 
   loninc = proj%loninc  

   rot_lon = rad_to_deg*atan( cos(xlat) * sin(xlon-cen_lon)/ &
             (cos(cen_lat)*cos(xlat)*cos(xlon-cen_lon) + sin(cen_lat)*sin(xlat)))
   rot_lat = rad_to_deg*asin(  cos(cen_lat)*sin(xlat) - sin(cen_lat)*cos(xlat)*cos(xlon-cen_lon))


   deltalat = rot_lat 
   deltalon = rot_lon 
    
   ! Compute x/y
   x = proj%knowni + deltalon/loninc + 1.0
   y = proj%knownj + deltalat/latinc + 1.0

   deallocate ( rot_lat, rot_lon, deltalat, deltalon, lon360, xlat, xlon )
  
   if (trace_use_frequent) call da_trace_exit("da_llxy_rotated_latlon_1d")

end subroutine da_llxy_rotated_latlon_1d
