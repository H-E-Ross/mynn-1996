subroutine da_write_iv_rad_ascii (it, ob, iv )

   !---------------------------------------------------------------------------
   ! Purpose: write out innovation vector structure for radiance data.
   !---------------------------------------------------------------------------

   implicit none

   integer      ,     intent(in)  :: it       ! outer loop count
   type (y_type),     intent(in)  :: ob       ! Observation structure.
   type (iv_type),    intent(in)  :: iv       ! O-B structure.

   integer                        :: n        ! Loop counter.
   integer                        :: i, k, l  ! Index dimension.
   integer                        :: nlevelss ! Number of obs levels.

   integer            :: ios, innov_rad_unit
   character(len=filename_len)  :: filename
   character(len=7)   :: surftype
   integer            :: ndomain
   logical            :: amsr2, ahi, abi
   real               :: cip ! to output cloud-ice path
   integer            :: cloudflag ! to output cloudflag
   integer, dimension(1) :: maxl
   integer, allocatable :: level_max(:)

   real, allocatable :: dtransmt(:,:), transmt_jac(:,:), transmt(:,:), lod(:,:), lod_jac(:,:)

   if (trace_use) call da_trace_entry("da_write_iv_rad_ascii")

   write(unit=message(1),fmt='(A)') 'Writing radiance OMB ascii file'
   call da_message(message(1:1))

   do i = 1, iv%num_inst
      if (iv%instid(i)%num_rad < 1) cycle

      ! count number of obs within the loc%proc_domain
      ! ---------------------------------------------
      ndomain = 0
      do n =1,iv%instid(i)%num_rad
         if (iv%instid(i)%info%proc_domain(1,n)) then
            ndomain = ndomain + 1
         end if
      end do
      if (ndomain < 1) cycle

      if (rtm_option==rtm_option_crtm .and. write_jacobian ) then
         allocate ( dtransmt(iv%instid(i)%nchan,iv%instid(i)%nlevels) )
         allocate ( transmt_jac(iv%instid(i)%nchan,iv%instid(i)%nlevels) )
         allocate ( transmt(iv%instid(i)%nchan,iv%instid(i)%nlevels) )
         allocate ( lod(iv%instid(i)%nchan,iv%instid(i)%nlevels) )
         allocate ( lod_jac(iv%instid(i)%nchan,iv%instid(i)%nlevels) )
      end if

      ! Get variables for maximum level of weighting function
      !  Only filled if calc_weightfunc .and. use_crtm_kmatrix both are true
      if ( rtm_option == rtm_option_crtm .and. calc_weightfunc .and. use_crtm_kmatrix ) then
	 allocate(level_max(iv%instid(i)%nchan))
      endif

      amsr2 = index(iv%instid(i)%rttovid_string,'amsr2') > 0
      abi   = index(iv%instid(i)%rttovid_string,'abi') > 0
      ahi   = index(iv%instid(i)%rttovid_string,'ahi') > 0

      write(unit=filename, fmt='(i2.2,a,i4.4)') it,'_inv_'//trim(iv%instid(i)%rttovid_string)//'.', myproc

      call da_get_unit(innov_rad_unit)
      open(unit=innov_rad_unit,file=trim(filename),form='formatted',iostat=ios)
      if (ios /= 0 ) then
         call da_error(__FILE__,__LINE__, &
            (/"Cannot open innovation radiance file"//filename/))
      Endif
      write(unit=innov_rad_unit,fmt='(a,a,i7,a,i5,a)') trim(iv%instid(i)%rttovid_string), &
                        ' number-of-pixels : ', ndomain, &
                        ' channel-number-of-each-pixel : ', iv%instid(i)%nchan, &
                        ' index-of-channels : '
      write(unit=innov_rad_unit,fmt='(10i5)') iv%instid(i)%ichan
      if ( amsr2 ) then
         write(unit=innov_rad_unit,fmt='(a)') ' pixel-info : i date scanpos landsea_mask  elv lat lon  satzen satazi solzen solazi clw'
      else
         write(unit=innov_rad_unit,fmt='(a)') ' pixel-info : i date scanpos landsea_mask  elv lat lon  satzen satazi solzen solazi'
      end if
      write(unit=innov_rad_unit,fmt='(a)') ' grid%xb-surf-info : i t2m mr2m(ppmv) u10 v10 ps ts smois tslb snowh isflg &
                    & soiltyp vegtyp vegfra elev clwp cloud_frac cip cloudflag'
      ndomain = 0
      do n =1,iv%instid(i)%num_rad
         if (iv%instid(i)%info%proc_domain(1,n)) then
            ndomain=ndomain+1
            if ( amsr2 ) then ! write out clw
               write(unit=innov_rad_unit,fmt='(a,i7,2x,a,i6,i3,f6.0,6f8.2,f8.3)') 'INFO : ', ndomain, &
                                iv%instid(i)%info%date_char(n), &
                                iv%instid(i)%scanpos(n),        &
                                iv%instid(i)%landsea_mask(n),   &
                                iv%instid(i)%info%elv(n),       &
                                iv%instid(i)%info%lat(1,n),     &
                                iv%instid(i)%info%lon(1,n),     &
                                iv%instid(i)%satzen(n),         &
                                iv%instid(i)%satazi(n),         &
                                iv%instid(i)%solzen(n),         &
                                iv%instid(i)%solazi(n),         &
                                iv%instid(i)%clw(n)
            else ! no clw info
               write(unit=innov_rad_unit,fmt='(a,i7,2x,a,i6,i3,f6.0,6f8.2)') 'INFO : ', ndomain, &
                                iv%instid(i)%info%date_char(n), &
                                iv%instid(i)%scanpos(n),        &
                                iv%instid(i)%landsea_mask(n),   &
                                iv%instid(i)%info%elv(n),       &
                                iv%instid(i)%info%lat(1,n),     &
                                iv%instid(i)%info%lon(1,n),     &
                                iv%instid(i)%satzen(n),         &
                                iv%instid(i)%satazi(n),         &
                                iv%instid(i)%solzen(n),         &
                                iv%instid(i)%solazi(n)
            end if
            select case (iv%instid(i)%isflg(n))
            case (0) ;
               surftype = ' SEA : '
            case (1) ;
               surftype = ' ICE : '
            case (2) ;
               surftype = 'LAND : '
            case (3) ;
               surftype = 'SNOW : '
            case (4) ;
               surftype = 'MSEA : '
            case (5) ;
               surftype = 'MICE : '
            case (6) ;
               surftype = 'MLND : '
            case (7) ;
               surftype = 'MSNO : '
            end select

            ! Output cloud-ice path, cloudflag, pressure of weighting function peak
            if (rtm_option==rtm_option_crtm .and. crtm_cloud ) then
	       cip = iv%instid(i)%cip(n)
	    else
	       cip = -9999.0
	    endif
            if ( ahi ) then
               cloudflag = iv%instid(i)%cloudflag(n)
            else
               cloudflag = -999
            endif

	    ! %nlevels is unstaggered, so subtract 1 to get number of mass points
	    ! Order is top down in CRTM
            if ( rtm_option == rtm_option_crtm .and. calc_weightfunc .and. use_crtm_kmatrix ) then
	       level_max(:) = -1
	       do l=1,iv%instid(i)%nchan
		  maxl(:) = maxloc(iv%instid(i)%der_trans(l,:,n)) ! returns index of maximum value; if a tie, returns first occurrence of maximum
	          level_max(l) = maxl(1) ! model level of weighting function peak for this pixel and channel, going from top-->down
	       enddo
	      !level_max(:) = ( iv%instid(i)%nlevels - 1 ) - level_max(:) + 1 ! make order bottom-up
	    endif

            write(unit=innov_rad_unit,fmt='(a,i7,9f10.2,3i3,f8.3,f10.2,f8.3,f15.5,f15.5,i7)') surftype, n, &
                             iv%instid(i)%t2m(n), &
                             iv%instid(i)%mr2m(n),   &
                             iv%instid(i)%u10(n), &
                             iv%instid(i)%v10(n),  &
                             iv%instid(i)%ps(n),  &
                             iv%instid(i)%ts(n),  &
                             iv%instid(i)%smois(n),  &
                             iv%instid(i)%tslb(n),  &
                             iv%instid(i)%snowh(n), &
                             iv%instid(i)%isflg(n), &
                             nint(iv%instid(i)%soiltyp(n)), &
                             nint(iv%instid(i)%vegtyp(n)), &
                             iv%instid(i)%vegfra(n), &
                             iv%instid(i)%elevation(n), &
                             iv%instid(i)%clwp(n), &
                             iv%instid(i)%cloud_frac(n), &
                             cip, &
                             cloudflag

            write(unit=innov_rad_unit,fmt='(a)') 'OBS  : '
            write(unit=innov_rad_unit,fmt='(10f11.2)') ob%instid(i)%tb(:,n)
            write(unit=innov_rad_unit,fmt='(a)') 'BAK  : '
            write(unit=innov_rad_unit,fmt='(10f11.2)') iv%instid(i)%tb_xb(:,n)
            if (rtm_option==rtm_option_crtm .and. crtm_cloud .and. (amsr2 .or. ahi .or. abi) ) then
               write(unit=innov_rad_unit,fmt='(a)') 'BAK_clr : '
               write(unit=innov_rad_unit,fmt='(10f11.2)') iv%instid(i)%tb_xb_clr(:,n)
            endif
            if ( rtm_option == rtm_option_crtm .and. calc_weightfunc .and. use_crtm_kmatrix ) then
               write(unit=innov_rad_unit,fmt='(a)') 'WEIGHTFUNC_PEAK : '
               write(unit=innov_rad_unit,fmt='(10f11.2)') iv%instid(i)%pm( (/level_max/), n )
            endif
            write(unit=innov_rad_unit,fmt='(a)') 'IVBC : '
            write(unit=innov_rad_unit,fmt='(10f11.2)')  iv%instid(i)%tb_inv(:,n)
            write(unit=innov_rad_unit,fmt='(a)') 'EMS  : '
            write(unit=innov_rad_unit,fmt='(10f11.2)')  iv%instid(i)%emiss(1:iv%instid(i)%nchan,n)
            if (rtm_option==rtm_option_crtm .and. write_jacobian) then
               write(unit=innov_rad_unit,fmt='(a)') 'EMS_JACOBIAN : '
               write(unit=innov_rad_unit,fmt='(10f10.3)') iv%instid(i)%emiss_jacobian(1:iv%instid(i)%nchan,n)
            end if
            write(unit=innov_rad_unit,fmt='(a)') 'ERR  : '
            write(unit=innov_rad_unit,fmt='(10f11.2)') iv%instid(i)%tb_error(:,n)
            write(unit=innov_rad_unit,fmt='(a)') 'QC   : '
            write(unit=innov_rad_unit,fmt='(10i11)') iv%instid(i)%tb_qc(:,n)
            if ( abi .and. crtm_cloud ) then ! write out cloud_mod, cloud_obs
               write(unit=innov_rad_unit,fmt='(a)') 'CMOD   : '
               write(unit=innov_rad_unit,fmt='(10f11.2)') iv%instid(i)%cloud_mod(:,n)
               write(unit=innov_rad_unit,fmt='(a)') 'COBS   : '
               write(unit=innov_rad_unit,fmt='(10f11.2)') iv%instid(i)%cloud_obs(:,n)
	       write(unit=innov_rad_unit,fmt='(a)') 'CLOUD    : '
               write(unit=innov_rad_unit,fmt='(10i11)') iv%instid(i)%cloud_flag(:,n)
            end if

            if (write_profile) then
               nlevelss  = iv%instid(i)%nlevels
               if ( rtm_option == rtm_option_rttov ) then
#ifdef RTTOV
                  ! first, write RTTOV levels
                  write(unit=innov_rad_unit,fmt='(a)') 'RTM_level pres(mb) T(k) Q(ppmv)'
                  do k = 1, nlevelss
                     write(unit=innov_rad_unit,fmt='(i3,f10.2,f8.2,e11.4)') &
                        k, &                             ! RTTOV levels
                        coefs(i) % coef % ref_prfl_p(k) , &
                        iv%instid(i)%t(k,n) , &
                        iv%instid(i)%mr(k,n)
                  end do  ! end loop RTTOV level
                  ! second, write WRF model levels
                  write(unit=innov_rad_unit,fmt='(a)') &
                     'WRF_level pres(mb) T(k) q(g/kg) clw(g/kg) rain(g/kg)'
                  do k=kts,kte
                     write(unit=innov_rad_unit,fmt='(i3,f10.2,f8.2,3e11.4)') &
                        k,  &                     ! WRF model levels
                        iv%instid(i)%pm(k,n) , &
                        iv%instid(i)%tm(k,n) , &
                        iv%instid(i)%qm(k,n)*1000 , &    
                        iv%instid(i)%qcw(k,n)*1000.0, &
                        iv%instid(i)%qrn(k,n)*1000.0
                  end do ! end loop model level
#endif
               end if ! end if rtm_option_rttov

               if ( rtm_option == rtm_option_crtm ) then
#ifdef CRTM
                  write(unit=innov_rad_unit,fmt='(a)') &
                     'level fullp(mb) halfp(mb) t(k) q(g/kg) water(mm) ice(mm) rain(mm) snow(mm) graupel(mm) hail(mm)'
                  if (crtm_cloud) then
                     do k=1,iv%instid(i)%nlevels-1
                        write(unit=innov_rad_unit,fmt='(i3,2f10.2,f8.2,13f8.3)') &
                           k,  &
                           iv%instid(i)%pf(k,n), &
                           iv%instid(i)%pm(k,n), &
                           iv%instid(i)%tm(k,n), &
                           iv%instid(i)%qm(k,n), &
                           iv%instid(i)%qcw(k,n), &
                           iv%instid(i)%qci(k,n), &
                           iv%instid(i)%qrn(k,n), &
                           iv%instid(i)%qsn(k,n), &
                           iv%instid(i)%qgr(k,n), &
                           iv%instid(i)%qhl(k,n), &
                           iv%instid(i)%rcw(k,n), &
                           iv%instid(i)%rci(k,n), &
                           iv%instid(i)%rrn(k,n), &
                           iv%instid(i)%rsn(k,n), &
                           iv%instid(i)%rgr(k,n), &
                           iv%instid(i)%rhl(k,n)
                     end do ! end loop profile
                  else  ! no cloud
                     do k=1,iv%instid(i)%nlevels-1
                        write(unit=innov_rad_unit,fmt='(i3,2f10.2,f8.2,7f8.3)') &
                           k,  &
                           iv%instid(i)%pf(k,n), &
                           iv%instid(i)%pm(k,n), &
                           iv%instid(i)%tm(k,n), &
                           iv%instid(i)%qm(k,n), &
                           0.0, &
                           0.0, &
                           0.0, &
                           0.0, &
                           0.0, &
                           0.0
                     end do ! end loop profile
                  end if  ! end if crtm_cloud
#endif
               end if ! end if rtm_option_crtm

            end if  ! end if write_profile

            if ( rtm_option == rtm_option_crtm .and. write_jacobian) then
#ifdef CRTM

               if ( calc_weightfunc ) then
                  dtransmt(:,:)    = iv%instid(i)%der_trans(:,:,n)
                  transmt(:,:)     = iv%instid(i)%trans(:,:,n)
                  transmt_jac(:,:) = iv%instid(i)%trans_jacobian(:,:,n)
                  lod(:,:)         = iv%instid(i)%lod(:,:,n)
                  lod_jac(:,:)     = iv%instid(i)%lod_jacobian(:,:,n)
               else
                  dtransmt(:,:)    = 0.0
                  transmt(:,:)     = 0.0
                  transmt_jac(:,:) = 0.0
                  lod(:,:)         = 0.0
                  lod_jac(:,:)     = 0.0
               end if

               write(unit=innov_rad_unit,fmt='(a)') &
                  'channel level halfp(mb) t(k) q(g/kg) der_trans trans_jac trans lod_jac lod water(mm) ice(mm) rain(mm) snow(mm) graupel(mm) hail(mm)'
               if (crtm_cloud) then
                  do l=1,iv%instid(i)%nchan
                     do k=1,iv%instid(i)%nlevels-1
                        write(unit=innov_rad_unit,fmt='(i5,i3,f10.2,13f14.7,6f14.7)') &
                           l, k,  &
                           iv%instid(i)%pm(k,n), &
                           iv%instid(i)%t_jacobian(l,k,n), &
                           iv%instid(i)%q_jacobian(l,k,n), &
                           dtransmt(l,k),&
                           transmt_jac(l,k),&
                           transmt(l,k),&
                           lod_jac(l,k),&
                           lod(l,k),&
                           iv%instid(i)%water_jacobian(l,k,n), &
                           iv%instid(i)%ice_jacobian(l,k,n), & 
                           iv%instid(i)%rain_jacobian(l,k,n), &
                           iv%instid(i)%snow_jacobian(l,k,n), &
                           iv%instid(i)%graupel_jacobian(l,k,n), &
                           iv%instid(i)%hail_jacobian(l,k,n), &
                           iv%instid(i)%water_r_jacobian(l,k,n), &
                           iv%instid(i)%ice_r_jacobian(l,k,n), & 
                           iv%instid(i)%rain_r_jacobian(l,k,n), &
                           iv%instid(i)%snow_r_jacobian(l,k,n), &
                           iv%instid(i)%graupel_r_jacobian(l,k,n), &
                           iv%instid(i)%hail_r_jacobian(l,k,n)
                    end do ! end loop profile
                 end do ! end loop channels
               else  ! no cloud
                  do l=1,iv%instid(i)%nchan
                     do k=1,iv%instid(i)%nlevels-1
                        write(unit=innov_rad_unit,fmt='(i5,i3,f10.2,13f14.7,6f14.7)') &
                           l, k,  &
                           iv%instid(i)%pm(k,n), &
                           iv%instid(i)%t_jacobian(l,k,n), &
                           iv%instid(i)%q_jacobian(l,k,n), &
                           dtransmt(l,k),&
                           transmt_jac(l,k),&
                           transmt(l,k),&
                           lod_jac(l,k),&
                           lod(l,k),&
                           0., &
                           0., &
                           0., &
                           0., &
                           0., &
                           0., &
                           0., &
                           0., &
                           0., &
                           0., &
                           0., &
                           0.
                     end do ! end loop profile
                  end do ! end loop channels
               end if  ! end if crtm_cloud
#endif
            end if !  end if write_jacobian

         end if ! end if proc_domain
      end do ! end do pixels
      if (rtm_option==rtm_option_crtm .and. write_jacobian ) then
         deallocate ( dtransmt )
         deallocate ( transmt_jac )
         deallocate ( transmt )
         deallocate ( lod )
         deallocate ( lod_jac )
      end if
      if ( rtm_option == rtm_option_crtm .and. calc_weightfunc .and. use_crtm_kmatrix ) then
	 deallocate(level_max)
      endif
      close(unit=innov_rad_unit)
      call da_free_unit(innov_rad_unit)
   end do ! end do instruments

   if (trace_use) call da_trace_exit("da_write_iv_rad_ascii")

end subroutine da_write_iv_rad_ascii

