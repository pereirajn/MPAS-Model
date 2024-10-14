!==================================================================================================================
 module module_wildf_plumerise
 use mpas_log
 use mpas_derived_types !for messageType=MPAS_LOG_ERR
 use mpas_kind_types,only: kind_phys => RKIND

 use plumerise, only: plumerise_run
 
 implicit none
 private
 public:: smk_pr_driver

 contains

!=================================================================================================================
 subroutine smk_pr_driver(&
                            its,ite, jts,jte, kts,kte,                                   &
                            zz, z, zmid, dz,                                             & 
                            t, th,                                                       &
                            u, v,                                                        &
                            rho, rho_zz,                                                 & 
                            pres,                                                        & 
                            qv,                                                          & 
                            ivgtyp,                                                      & 
                            pm25_wild_emis_flux,                                         &
                            pm25wf_PRMsource                                             &
                            )
  

   !*******************************************************************************
    implicit none
   !*******************************************************************************

   !-------------------------------------------------------------------------------
   ! Variables and arrays from mpas_atmphys_driver_pbl
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
   !-- rho_zz,                dry air density divided by d(zeta)/dz  [kg/m3]
   !-- pres,                  pressure [Pa]
   !-- qv,                    water vapor mixing ratio (kg/kg)
   !-- ivgtyp,                dominant vegetation category [-]
   !-- pm25_wild_emis_flux    2d pm25 wildfire emissions at surface   [kg.m-2.s-1]
   !-- pm25wf_PRMsource       3d pm25 wildfire emissions after PRM calculations    [kg.m-2.s-1]

    integer,  intent(in   )   ::         its,ite, jts,jte, kts,kte

    real(kind=kind_phys),     dimension( its:ite, kts:kte, jts:jte )             , &
             intent(in   )   ::                                                zz, &
                                                                                z, &
                                                                             zmid, &
                                                                               dz, &
                                                                                t, &
                                                                               th, &
                                                                                u, &
                                                                                v, &
                                                                              rho, &
                                                                           rho_zz, &
                                                                             pres, &
                                                                               qv

    integer, dimension( its:ite, jts:jte )                                       , &
             intent(in   )   ::                                            ivgtyp   

    real(kind=kind_phys),     dimension( its:ite, jts:jte )                      , &
             intent(in   )   ::                                pm25_wild_emis_flux      

    real(kind=kind_phys),     dimension( its:ite, kts:kte, jts:jte )             , &
             intent(out   )   ::                                  pm25wf_PRMsource


   !-------------------------------------------------------------------------------
   ! Local Variables and Arrays to/from plumerise_run
   !-------------------------------------------------------------------------------
    
    real(kind=kind_phys):: &
       zbot_out,               & !
       ztop_out                  !

    integer, dimension( its:ite, jts:jte ) :: iveg_ag   

   !-------------------------------------------------------------------------------
   ! Local Variables and Arrays 
   !-------------------------------------------------------------------------------
    integer ::  i,j,k
    real(kind=kind_phys) :: flame_fraction
    
   
   call mpas_log_write(' ')
   call mpas_log_write('--- enter subroutine smk_pr_driver:')
    

   !*******************************************************************************
   ! Determines the vegetation combination from MPAS to PRM categories
   !*******************************************************************************
     
     !call mpas_log_write('--- vegetation category combination', messageType=MPAS_LOG_ERR)
     
     do j = jts,jte
         do i = its, ite
            select case(ivgtyp(i, j))
            case(2)
               iveg_ag(i,j) = 1 !Tropical Forest        !igbp = Evergreen Broadleaf Forest
            case(1,3:5)
               iveg_ag(i,j) = 2 !Boreal Forest          !igbp = Evergreen|Deciduos Needleleaf, Deciduous Broadleaf, and Mixed
            case(6:9,11)
               iveg_ag(i,j) = 3 !Cerrado/woody savanna  !igbp = Shrublands, savannas, Permanent wetlands(pantanal)?
            case(10,12,14)
               iveg_ag(i,j) = 4 !Grassland/cropland     !igbp = grassland, cropland, natural cropland
            case default
               iveg_ag(i,j) = 0  ! Default value if not in the range
            end select

            !call mpas_log_write('ivgtyp(i,j)=$i iveg_ag(i,j)=$i', intArgs=(/ivgtyp(i,j), iveg_ag(i,j)/), messageType=MPAS_LOG_ERR)

         enddo
      enddo


   !*******************************************************************************
   ! Call subroutine plumerise_run in the plumerise.F
   !*******************************************************************************
    
   !call mpas_log_write('kte=$i', intArgs=(/kte/))
   !call mpas_log_write('qv(its+300,kts+20,jte)=$r', realArgs=(/qv(its+300,kts+20,jte)/))
   !call mpas_log_write('t(its+300,kts+20,jte)=$r', realArgs=(/t(its+300,kts+20,jte)/))
   !call mpas_log_write('z(its+300,kts+20,jte)=$r', realArgs=(/z(its+300,kts+20,jte)/))
   

   call plumerise_run( &
         its = its, ite = ite, jts = jts, jte = jte, kts = kts, kte = kte                , &
         temp_in       = t         , pres_in       = pres      , u_in       = u          , &
         v_in          = v         , hml_in        = zmid      , qv_in      = qv         , &
         z_in          = z                                                               , &  
         zbot_out      = zbot_out  , ztop_out      = ztop_out                            , & 
         iveg_ag_in    = iveg_ag                                                         , & 
         pm25_wild_emis_flux_in = pm25_wild_emis_flux                                      &  
         )



   !******************************************************************************* 
   ! test routine
   !*******************************************************************************


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


       ! Se a classificacao de vegetacao permanecer com iveg_ag, mudar essa rotina para uma look up table 
       !Determine the flame fraction according to vegetation type for the current (i, j) grid
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
                                      / ( (rho_zz(i, k, j) * dz(i,k,j)) )
       end do
       

       ! Emissions at PRM levels
       do k = kts+20,kts+30
            pm25wf_PRMsource(i,k,j) = (flame_fraction * pm25_wild_emis_flux(i,j)) &
                                      / ( (rho_zz(i, k, j) * dz(i,k,j)) )
       end do


     end do
  end do


   
   call mpas_log_write('--- exit subroutine smk_pr_driver:')
   call mpas_log_write(' ')

   end subroutine smk_pr_driver
!=================================================================================================================



!=================================================================================================================
 end module module_wildf_plumerise
!=================================================================================================================