!=================================================================================================================
 module module_wildf_plumerise
 use mpas_log
 use mpas_kind_types,only: kind_phys => RKIND

 implicit none
 private
 public:: smk_pr_driver

 contains


!=================================================================================================================
 subroutine smk_pr_driver(&
                            ids,ide, jds,jde, kds,kde,                                   &
                            ims,ime, jms,jme, kms,kme,                                   &
                            its,ite, jts,jte, kts,kte,                                   &
                            zz, z, zmid, dz,                                             & 
                            t, th,                                                       &
                            u, v,                                                        &
                            rho, pres,                                                   & 
                            qv,                                                          & 
                            ivgtyp,                                                      & 
                            pm25_wild_emis_flux,                                         &
                            pm25wf_PRMsource                                             &
                            )


   

   !-------------------------------------------------------------------------------
    implicit none
   !-------------------------------------------------------------------------------
   !-- zz,                    d(zeta)/dz, vertical metric term [-]
   !-- z,                     height of layer  [m]
   !-- zmid,                  height of middle of layer [m]
   !-- dz,                    layer thickness  [m]
   !-- t,                     temperature  [K]
   !-- th,                    potential temperature (k)
   !-- u,                     u-velocity interpolated to theta points (m/s)
   !-- v,                     v-velocity interpolated to theta points (m/s)
   !-- rho,                   air density [kg/m3]
   !-- pres,                  pressure [Pa]
   !-- qv,                    water vapor mixing ratio (kg/kg)
   !-- ivgtyp,                dominant vegetation category [-]
   !-- pm25_wild_emis_flux    2d pm25 wildfire emissions at surface   [kg.m-2.s-1]
   !-- pm25wf_PRMsource       3d pm25 wildfire emissions after PRM calculations    [kg.m-2.s-1]

    integer,  intent(in   )   ::         ids,ide, jds,jde, kds,kde,                &
                                         ims,ime, jms,jme, kms,kme,                &
                                         its,ite, jts,jte, kts,kte

    real(kind=kind_phys),     dimension( ims:ime, kms:kme, jms:jme )             , &
             intent(in   )   ::                                                zz, &
                                                                                z, &
                                                                             zmid, &
                                                                               dz, &
                                                                                t, &
                                                                               th, &
                                                                                u, &
                                                                                v, &
                                                                              rho, &
                                                                             pres, &
                                                                               qv

    integer, dimension( ims:ime, jms:jme )                                       , &
             intent(in   )   ::                                            ivgtyp   

    real(kind=kind_phys),     dimension( ims:ime, jms:jme )                      , &
             intent(in   )   ::                                pm25_wild_emis_flux      

    real(kind=kind_phys),     dimension( ims:ime, kms:kme, jms:jme )             , &
             intent(out   )   ::                                  pm25wf_PRMsource

   





!------------------------------- !test routine! ----------------------------------------------------------
    !local
    integer ::  i,j,k
    real(kind=kind_phys) :: flame_fraction


   call mpas_log_write(' ')
   call mpas_log_write('--- enter subroutine smk_pr_driver:')


   !-----------------------------------------------------------------------------------------------
   !Flaming Consumption Freitas et al. (2006)  | LANDUSE.TBL = MODIFIED IGBP MODIS NOAH Land Cover! 
   !-----------------------------------------------------------------------------------------------
   ! 45,   ! Tropical Forest        |  igbp 1 thru 5 (Evergreen|Deciduos Needleleaf|Broadleaf Forest and Mixed)
   ! 75,   ! cerrado/woody savanna  |  igbp 6 thru 9 and 11 (Shrublands, savannas, Permanent wetlands(pantanal)?)
   ! 97    ! Grassland/cropland     |  igbp 10, 12 and 14 (grassland, cropland, natural cropland)
   !  0    ! case default           !  Tundra, barren vegetated, water, urban, ice
   !----------------------------------------------------------------------------------------------

   
   do j = jts,jte
      do i = its, ite


       ! Determine the vegetation type for the current (i, j) grid
       select case(ivgtyp(i, j))
         case(1:5)
            flame_fraction = 0.45
         case(6:9,11)
            flame_fraction = 0.75
         case(10,12,14)
            flame_fraction = 0.97
         case default
            flame_fraction = 0.0  ! Default value if not in the range
       end select
      

       ! Emissions at surface
       do k = kts,kts
            pm25wf_PRMsource(i,k,j) = ((1-flame_fraction) * pm25_wild_emis_flux(i,j)) &
                                      / ( (rho(i, k, j) * dz(i,k,j) / zz(i,k,j)) )
       end do
       

       ! Emissions at PRM levels
       do k = kts+20,kts+30
            pm25wf_PRMsource(i,k,j) = (flame_fraction * pm25_wild_emis_flux(i,j)) &
                                      / ( (rho(i, k, j) * dz(i,k,j)) / zz(i,k,j) )
       end do


     end do
  end do


   
   call mpas_log_write('--- exit subroutine smk_pr_driver:')
   call mpas_log_write(' ')

   end subroutine smk_pr_driver
!----------------------------------------------------------------------------------------------------------



!=================================================================================================================
 end module module_wildf_plumerise
!=================================================================================================================
