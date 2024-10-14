module plumerise
use mpas_log
use mpas_derived_types !for messageType=MPAS_LOG_ERR
use mpas_kind_types,only: kind_phys => RKIND


implicit none
 private
 public:: plumerise_run


 !=================================================================================================================
 ! Global configurations

     logical :: do_windeffect = .true.   !- turn ON (true), OFF (false) the env wind effect on plume rise
     integer :: sound_type = 0           !- set it to zero for the coupling with 3d host model. set 1 for reading a precreated sound already defined at PRM vertical grid
     integer, parameter :: heat_flag = 1 !- set 1 for use prescribed heat flux values. set 2 to use heat flux from FRP


 ! Global physical constants

     real(kind=kind_phys), parameter :: pi    = 3.14159265358979323846264338327950288
     

 ! Global variables and parameters    
     
     integer, parameter :: nveg_ag =   4   !cw nkp war 200
     integer, parameter ::   nkp   = 200   !cw nkp war 200
     integer, parameter ::   ntime = 200   !cw nkp war 200

     real(kind=kind_phys), dimension(nkp) ::  &
       w,      & !
       t,      & !
       theq,   & !
       qv,     & !
       qc,     & !
       qh,     & !
       qi,     & !
       sc,     & !  
       vth,    & !
       vti,    & !
       rho,    & !
       txs,    & !
       est,    & !
       dqsdz,  & !
       qsat,   & !
       qpas,   & !
       qtotal, & !
       td,     & !
       vel_p,  & !
       rad_p,  & !
       alpha
  real(kind=kind_phys), dimension(nkp) ::  &
       wc,    & !
       wt,    & ! cache variable initially set to zero by zero_plumegen_coms
       tt,    & !
       qvt,   & !
       qct,   & !
       qht,   & !
       qit,   & !
       sct,   & !
       vel_t, & !
       rad_t    !
  real(kind=kind_phys), dimension(nkp) ::  &
       dzm,    & ! used as inverse of full level vertical grid (inverse of interface of the layer)
       dzt,    & ! used as inverse of half level vertical grid (inverse of midle of the layer)
       vctr1,  & !
       vctr2,  & !
       vt3dc,  & !
       vt3df,  & !
       vt3dk,  & !
       vt3dg,  & !
       scr1,   & !
       vt3dj,  & !
       vt3dn,  & !
       vt3do,  & !
       vt3da,  & !
       scr2,   & ! 
       vt3db,  & !
       vt3dd,  & !
       vt3dl,  & !
       vt3dm,  & !
       vt3di,  & !
       vt3dh,  & !
       vt3de,  & !
       rbuoy,  & !
       dwdt_entr !
  real(kind=kind_phys), dimension(nkp) ::  & ! environment at plume grid
       the,  & !
       thee, & !
       rhe,  & !
       sce,  & !
       tde     ! 
  real(kind=kind_phys), dimension(nkp) :: &
       te,    & !
       pe,    & !
       upe,   & !
       vpe,   & !
       qvenv, & ! water vapor mixing ratio at PRM grid (kg h20v/kg dry air)
       vel_e, & !
       zt,    & ! thermo and water levels in the PRM grid (midle of the layer)
       zm       ! dynamical levels in the PRM model (interface of the layer)
  real(kind=kind_phys), dimension(1200) :: & ! environment at RAMS  grid
       wcon,   & !
       thtcon, & !
       zzcon,  & !
       scon      !
  real(kind=kind_phys) :: &
       dz,        & !
       visc(nkp), & !
       viscosity, & !
       tstpf        !   
  integer :: &
       n,   & !
       nm1, & !
       l      !counter
  real(kind=kind_phys) :: & ! advection
       advw,      & !
       advt,      & !
       advv,      & !
       advc,      & !
       advh,      & !
       advi,      & !
       cvh(nkp),  & !
       cvi(nkp),  & !
       adiabat,   & !
       wbar,      & !
       vhrel,     & !
       virel        !  
  real(kind=kind_phys) :: &
       zsurf, & ! set to 0 in zero_plumegen_coms and never changed afterwards
       zbase, & !
       ztop     !
  integer :: LBASE
  real(kind=kind_phys) :: & ! entrain
       area,     & !
       rsurf,    & ! entrainment surface radius (m) set in get_fire_properties
       radius(nkp) !  
  real(kind=kind_phys) :: & ! heating
       heating(ntime), & ! heating rates computed in get_fire_properties
       fmoist,         & ! fuel moisture fraction set to 1/10 in get_fire_properties
       bload,          & ! total loading set to 10 in get_fire_properties
       heat_fluxW        ! cache variable for heatflux in get_fire_properties
  real(kind=kind_phys) :: &
       dt,   & !
       time, & !
       tdur    ! number of seconds in the burn set in get_fire_properties to 60*mdur
  integer :: &
       mintime, & ! elapsed time, minutes - makeplume
       mdur,    & ! duration of burn, minutes 
       maxtime    ! set in in get_fire_properties
  real(kind=kind_phys) :: ztop_(ntime) !
 
  integer :: &
       kmt,     & ! will be set to nkp in the subroutine get_env_condition
       imm        ! counter for lower and upper boundary of the heat flux and for the plume
 !=================================================================================================================



 contains

!=================================================================================================================
 subroutine plumerise_run( &
               its,ite, jts,jte, kts,kte                , &
               temp_in      , pres_in      , u_in       , &
               v_in         , hml_in       , qv_in      , &
               z_in                                     , &
               zbot_out     , ztop_out                  , &
               iveg_ag_in                               , & 
               pm25_wild_emis_flux_in                     &  
               )

   !*******************************************************************************
    implicit none
   !*******************************************************************************

   !---------------------------------------------------------------
   ! Variables and arrays from/to module_wildf_plumerise 
   !---------------------------------------------------------------

    integer,  intent(in   )   ::  &
        its,ite,          & !index for the host model horizontal grid (its=1, ite = nCellsSolve)
        jts,jte,          & !index from WRF model (jts=jte=1)
        kts,kte             !index for the host model vertical grid (kts=1, kte = nVertLevels)

    real(kind=kind_phys), dimension( its:ite, kts:kte, jts:jte ), intent(in   )   :: &
        temp_in,       & ! Temperature of the host model (k)interpolated to theta points
        pres_in,       & ! Pressure of the host model (Pa) interpolated to theta points
        u_in,          & ! zonal wind of the host model (m/s) interpolated to theta points 
        v_in,          & ! meridional wind of the host model (m/s) interpolated to theta points 
        hml_in,        & ! geometric height of middle of layer of the host model (m MSL = meters above mean sea level)
        qv_in,         & ! water vapor mixing ratio of the host model (kg h20vapor/kg dry air)
        z_in             ! geometric height of layer interfaces of the host model [m MSL] located at w points


      
    real(kind=kind_phys), intent(inout   ) :: &
       zbot_out, & ! lowest height for smoke emissions (m)
       ztop_out    ! highest height for smoke emissions (m)


    integer, dimension( its:ite, jts:jte ),  intent(in   )   :: &
       iveg_ag_in                  ! dominant combined vegetation category   
       

   real(kind=kind_phys),dimension( its:ite, jts:jte ),  intent(in   )   :: &
       pm25_wild_emis_flux_in      ! pm25 emission flux (kg.m-2.s-1) read in by the host model

   !---------------------------------------------------------------
   ! Local Variables and arrays  
   !---------------------------------------------------------------
    
    
    integer ::  &
        m1,        & ! index of the upper most layer (numbered from bottom to top)
         i,        & ! horizontal index
         k,        & ! vertical index
         j,        & ! j index 
       imm           ! counter for lower and upper boundary of the heat flux and for the plume. !pq define aqui se e uma variavel global?
   
    real(kind=kind_phys), pointer :: &
        ucon(:),   & ! zonal wind (m/s)
        zcon(:),   & ! cartesian height (meters)
        vcon(:),   & ! meridional wind (m/s)
        tmpcon(:), & ! temperature
        prcon(:),  & ! pressure
        rvcon(:)     ! water vapor mixing ratio (kg h20v/kg dry air)

    integer, pointer :: iveg_ag

     real(kind=kind_phys) :: &
         & heat_flux_min, & ! minimal heat flux 
         & heat_flux_max, & ! maximum heat flux 
         & burnt_area       ! burnt area, fire size (m2)

   !********************************************************************************
   ! Initial settings
   !********************************************************************************
    
    call mpas_log_write(' ')
    call mpas_log_write('--- enter subroutine plumerise_run:')
    !call mpas_log_write('kte=$i', intArgs=(/kte/))
    !call mpas_log_write('qv_in(its+300,kts+20,jte)=$r', realArgs=(/qv_in(its+300,kts+20,jte)/))
    !call mpas_log_write('temp_in(its+300,kts+20,jte)=$r', realArgs=(/temp_in(its+300,kts+20,jte)/))
    !call mpas_log_write('z_in(its+300,kts+20,jte)=$r', realArgs=(/z_in(its+300,kts+20,jte)/))
    !call mpas_log_write('iveg_ag_in(its+300,jte)=$i', intArgs=(/iveg_ag_in(its+300,jte)/))
        
    !-allocate the local arrays to MPAS atmospheric conditions at 1D
    allocate(ucon(kte)  )
    allocate(zcon(kte)  )
    allocate(vcon(kte)  )
    allocate(tmpcon(kte)  )
    allocate(prcon(kte)  )
    allocate(rvcon(kte)  )
    allocate(iveg_ag)

    !Define local index of the MPAS upper most layer
    m1=kte 
    call mpas_log_write('MPAS upper most layer (m1)=$i', intArgs=(/m1/))

    !-initialize several parameters (only need to be done at the 1st time)
    call zero_plumegen_coms(ucon,zcon,vcon,tmpcon,prcon,rvcon)

    
   !********************************************************************************
   ! Loop on host model 
   !********************************************************************************
    
    do j = jts,jte
          call mpas_log_write('subdomain', messageType=MPAS_LOG_ERR)
          call mpas_log_write('pm25_wild_emis_flux_in min=$r pm25_wild_emis_flux_in max=$r', realArgs=(/minval(pm25_wild_emis_flux_in),maxval(pm25_wild_emis_flux_in)/), messageType=MPAS_LOG_ERR)
     
     do i = its, ite
     
          !---------------------------------------------------------------
          !- if the value of emission is close to zero, 
          !  there is not emission with plume rise => cycle  
          !---------------------------------------------------------------
            if (pm25_wild_emis_flux_in(i, j) .lt. 1.e-12) then
             
               !call mpas_log_write('nCellsSolve=$i', intArgs=(/i/))
               !call mpas_log_write('pm25_wild_emis_flux_in(i, j)=$r', realArgs=(/pm25_wild_emis_flux_in(i, j)/))
               !call mpas_log_write('No Fire: PRM aborted')
              
               cycle

               !add here o ztop e zbot =-9999 for the cycle case, so the model dont assign arbitrary values
               !check if this minimum value it is ok 

             else
               call mpas_log_write('-----------------------------', messageType=MPAS_LOG_ERR)
               call mpas_log_write('Fire detected: run PRM ', messageType=MPAS_LOG_ERR)
               call mpas_log_write('nCellsSolve=$i', intArgs=(/i/), messageType=MPAS_LOG_ERR)
               call mpas_log_write('pm25_wild_emis_flux_in(i, j)=$r', realArgs=(/pm25_wild_emis_flux_in(i, j)/),messageType=MPAS_LOG_ERR)
               
            endif
            

          !---------------------------------------------------------------
          ! read variables of host model  
          !---------------------------------------------------------------
           
           call mpas_log_write(' ', messageType=MPAS_LOG_ERR)
           call mpas_log_write('--- read variables of host model:', messageType=MPAS_LOG_ERR)

           !- fill local variables with the host model variables (MPAS grid)
            do k = kts,m1                                    ! z-direction of the host model 
               rvcon (k) = qv_in(i,k,j)                      ! fill with the water vapor mixing ratio (kg/kg) of the host model 
               zcon  (k) = hml_in(i,k,j) - z_in(i,kts,j)     ! fill with the cartesian height (meters) of the termodynamic grid point  
                                                        !(half-levels) above local surface--> muss full level sein, siehe aufruf htint
               tmpcon(k) = temp_in(i,k,j)                    !
               prcon (k) = pres_in(i,k,j)                    !
               ucon  (k) = u_in(i,k,j)                       !fill with the zonal wind (m/s) of the host model
               vcon  (k) = v_in(i,k,j)                       !fill with the meridional wind (m/s) of the host model
            enddo
        
           !- fill with the vegetation type combinated in the module_wildf_plumerise
           iveg_ag = iveg_ag_in(i,j)  

          call mpas_log_write(' qv_in(i,kts,j)=$r rvcon (kts)=$r', realArgs=(/qv_in(i,kts,j),rvcon(kts)/), messageType=MPAS_LOG_ERR)
          call mpas_log_write(' temp_in(i,kts,j)=$r tmpcon(kts)=$r', realArgs=(/temp_in(i,kts,j),tmpcon(kts)/), messageType=MPAS_LOG_ERR)
          call mpas_log_write(' hml_in(i,kts,j)=$r z_in(i,kts,j)=$r zcon(kts)=$r ', realArgs=(/hml_in(i,kts,j), z_in(i,kts,j), zcon(kts)/), messageType=MPAS_LOG_ERR)
          call mpas_log_write('iveg_ag_in(i,j)=$i iveg_ag=$i ', intArgs=(/iveg_ag_in(i,j), iveg_ag/), messageType=MPAS_LOG_ERR)

          
          !------------------------------------------------------------------------------
          ! get envinronmental state in the PRM grid (temp, water vapor mix ratio, ...)  
          !----------------------------------------------------------------------------
           
           call get_env_condition(kte,kmt,ucon,zcon,vcon,tmpcon,prcon,rvcon)


          !--------------------------------------------------------------------------------
          
          !uncomment the loop in iveg_ag in case of the host model has more than 1 ivgtyp in the same cell 
          !do iveg_ag = 1, nveg_ag !loop on vegetation type
          
           burnt_area=50.*1e4 !(m^2)
           
              do imm=1,2         !loop on heat_flux min and max
          
               !------------------------------------------------------------------------------
               ! get fire properties (burned area, plume radius, heating rates ...)  
               !----------------------------------------------------------------------------
                call get_fire_properties(burnt_area,heat_flux_max,heat_flux_min,iveg_ag, imm)
              
     
               !------------------------------------------------------------------------------
               ! generate plume rise 
               !----------------------------------------------------------------------------


              enddo !end loop on imm (heat_flux min and max)
          
          !enddo !end loop on vegetation type. uncomment when has more than one ivgtyp per cell
      enddo !end loop on i
     enddo !end loop on j 

     


   !********************************************************************************
   ! Final settings
   !********************************************************************************
    
    !-deallocate the local arrays that received MPAS atmospheric conditions at 1D
    deallocate(ucon  )
    deallocate(zcon )
    deallocate(vcon  )
    deallocate(tmpcon )
    deallocate(prcon )
    deallocate(rvcon )
    deallocate(iveg_ag)
   
    call mpas_log_write(' ')
    call mpas_log_write('--- exit subroutine plumerise_run:')

  end subroutine plumerise_run
!=================================================================================================================

!=================================================================================================================
subroutine zero_plumegen_coms(ucon,zcon,vcon,tmpcon,prcon,rvcon)

    implicit none
    real(kind=kind_phys), intent(inout) :: ucon(:),zcon(:),vcon(:),tmpcon(:),prcon(:),rvcon(:)

    call mpas_log_write('    --- enter subroutine zero_plumegen_coms:')

    w=0. ;t=0. ;td=0. ;theq=0. 
    qv=0. ;qc=0. ;qh=0. ;qi=0. ;sc=0. ;vel_p=0. ;rad_p=0. ;rad_t=0. 
    vth=0. ;vti=0. ;rho=0. ;txs=0. 
    vt3dc=0. ;vt3df=0. ;vt3dk=0. ;vt3dg=0. ;scr1=0. ;vt3dj=0. ;vt3dn=0. ;vt3do=0. 
    vt3da=0. ;scr2=0. ;vt3db=0. ;vt3dd=0. ;vt3dl=0. ;vt3dm=0. ;vt3di=0. ;vt3dh=0. ;vt3de=0. 
    wc=0. ;wt=0. ;tt=0. ;qvt=0. ;qct=0. ;qht=0. ;qit=0. ;sct=0. ;vel_t=0.
    est=0. ;qsat=0. ;qpas=0. ;qtotal=0. 
    dzm=0. ;dzt=0. ;zm=0. ;zt=0. ;vctr1=0. ;vctr2=0. 
    the=0. ;thee=0. ;rhe=0. ;sce=0. 
    ucon=0.;vcon=0.;wcon=0.;thtcon =0.;rvcon=0.;tmpcon=0.;prcon=0. 
    zcon=0.;zzcon=0.;scon=0. 
    dz=0.;dqsdz=0.;visc=0.;viscosity=0.;tstpf=0.
    advw=0.;advt=0.;advv=0.;advc=0.;advh=0.;advi=0.;cvh=0.;cvi=0.;adiabat=0.
    wbar=0.;vhrel=0.;virel=0.  
    zsurf=0.;zbase=0.;ztop=0.;area=0.;rsurf=0.;radius=0.;heating=0.
    fmoist=0.;bload=0.;dt=0.;time=0.;tdur=0.
    ztop_=0.
    n=0;nm1=0;l=0;lbase=0;mintime=0;mdur=0;maxtime=0
    dwdt_entr=0.
    kmt = 0
    vel_e(:)=0.;pe(:)=0.;te(:)=0.;qvenv(:)=0.;upe(:)=0.;vpe(:)=0.
    

    
    call mpas_log_write('w=$r', realArgs=(/w/))
    call mpas_log_write('ucon=$r', realArgs=(/ucon/))
    !call mpas_log_write('ztop=$r', realArgs=(/ztop/))

  end subroutine zero_plumegen_coms
!=================================================================================================================


!=================================================================================================================
subroutine get_env_condition(kte,kmt,ucon,zcon,vcon,tmpcon,prcon,rvcon)

   !*******************************************************************************
    implicit none
   !*******************************************************************************

    real(kind=kind_phys), intent (in) :: & 
         ucon(:),   & 
         zcon(:),   & 
         vcon(:),   & 
         tmpcon(:), & 
         prcon(:),  & 
         rvcon(:)     
    integer, intent (in) ::  kte  
    integer, intent (out) :: kmt
    integer :: & 
         k,    & !
         nk      !
         
    real(kind=kind_phys) :: &
         znz  !
         

    call mpas_log_write('--- enter subroutine get_env_condition:', messageType=MPAS_LOG_ERR)

   !*******************************************************************************
   ! Define vertical grid of plume model
   !*******************************************************************************
     call set_grid()
     
     !what this part is checking?
     znz=zcon(kte)  
     do k=nkp,1,-1 
       if(zt(k).LT.znz) then
          kmt=k
          exit
       endif
     enddo
    

   !*******************************************************************************
   ! Interpolates variables from MPAS grid to PRM grid
   !*******************************************************************************
     kmt=nkp 
     nk = kte
     
     if(sound_type .ne. 1) then 
      call htint(nk,   ucon, zcon, kmt,   upe,  zt)
      call htint(nk,   vcon, zcon, kmt,   vpe,  zt)
      call htint(nk,  rvcon, zcon, kmt, qvenv,  zt)
      call htint(nk, tmpcon, zcon, kmt,    te,  zt)
      call htint(nk,  prcon, zcon, kmt,    pe,  zt)
    else !Essa parte permanece mesmo que nao vamos usar radiossonda aqui? 
       te (:)=tmpcon(:)
       upe(:)=ucon(:)
       vpe(:)=vcon(:)
       pe (:)=prcon(:)
       qvenv(:)=rvcon(:)
    endif

   !*******************************************************************************
   ! Additional configurations: units and wind_effect
   !*******************************************************************************
     
     do k=1,kmt
       qvenv(k)=MAX(qvenv(k),1e-8)
     enddo

     ! environmental horizontal wind computed from u and v
     do k=1,kmt 
       vel_e(k) = sqrt(upe(k)**2+vpe(k)**2)
     enddo

     !ewe - env wind effect
     if(.not. do_windeffect)  vel_e = 0.

     !converts pressure from Pa to kPa
     do k=1,kmt 
      pe(k) = pe(k)*1.e-3
     enddo
     
end subroutine get_env_condition
!=================================================================================================================


!=================================================================================================================
subroutine set_grid()

   !*******************************************************************************
    implicit none
   !*******************************************************************************
     integer :: &
         k,   & ! loop index
         mzp    ! loop upper boundary

   call mpas_log_write('--- enter subroutine set_grid:', messageType=MPAS_LOG_ERR)

   !*******************************************************************************
   ! PRM grid configuration
   !*******************************************************************************
  
    dz=100. ! set constant grid spacing of plume grid model(meters)

    mzp=nkp 
    zt(1) = zsurf 
    zm(1) = zsurf 
    zt(2) = zt(1) + 0.5*dz
    zm(2) = zm(1) + dz
    do k=3,mzp
       zt(k) = zt(k-1) + dz ! thermo and water levels
       zm(k) = zm(k-1) + dz ! dynamical levels    
    enddo
       
    do k = 1,mzp-1
       dzm(k) = 1. / (zt(k+1) - zt(k))
    enddo
    dzm(mzp)=dzm(mzp-1) 
    
    do k = 2,mzp
       dzt(k) = 1. / (zm(k) - zm(k-1))
    enddo
    dzt(1) = dzt(2) * dzt(2) / dzt(3) 
    
    !   dzm(1) = 0.5/dz
    !   dzm(2:mzp) = 1./dz


   !call mpas_log_write('zsurf=$r', realArgs=(/zsurf/), messageType=MPAS_LOG_ERR)  
   !do k=1,mzp
   !  call mpas_log_write('zt(k)=$r', realArgs=(/zt(k)/), messageType=MPAS_LOG_ERR)
   !end do

end subroutine set_grid
!=================================================================================================================


!=================================================================================================================
subroutine htint (nzz1, vctra, eleva, nzz2, vctrb, elevb)

   !*******************************************************************************
    implicit none
   !*******************************************************************************
   
    integer,              intent(in ) :: nzz1        ! number of level from the host model (nlev)
    integer,              intent(in ) :: nzz2        ! number of level set by the global parameter nkp
    real(kind=kind_phys), intent(in ) :: vctra(nzz1) ! condition of environment from the host model 
    real(kind=kind_phys), intent(out) :: vctrb(nzz2) ! environmental conditions in the plume -  initially set to zero by zero_plumegen_coms
    real(kind=kind_phys), intent(in ) :: eleva(nzz1) ! vertical grid of the host model
    real(kind=kind_phys), intent(in ) :: elevb(nzz2) ! vertical grid set by set_grid()
    
    integer  :: l  ! counter  !pq declara aqui se e variavel global?
    integer  :: k  ! counter
    integer  :: kk ! counter for print(log.write) 
    real(kind=kind_phys) :: wt !pq declara aqui se e variavel global?

    call mpas_log_write('--- enter subroutine htint:', messageType=MPAS_LOG_ERR)

   !*******************************************************************************
   !Grid interpolation
   !*******************************************************************************
   
     l=1

     do k=1,nzz2
       do
          if ( (elevb(k) <  eleva(1)) .or. &
               ((elevb(k) >= eleva(l)) .and. (elevb(k) <= eleva(l+1))) ) then
             wt       = (elevb(k)-eleva(l))/(eleva(l+1)-eleva(l))
             vctrb(k) = vctra(l)+(vctra(l+1)-vctra(l))*wt
             exit
          else if ( elevb(k) >  eleva(nzz1))  then
             wt       = (elevb(k)-eleva(nzz1))/(eleva(nzz1-1)-eleva(nzz1))
             vctrb(k) = vctra(nzz1)+(vctra(nzz1-1)-vctra(nzz1))*wt
             exit
          end if
          
          !Essa parte aqui permance? Se sim vou trocar o print por log.write
          l=l+1
          if(l == nzz1) then
             ! PRINT *,'htint:nzz1',nzz1 !commented out for vectorisation 
             ! DO kk=1,l
             !   PRINT*,'kk,eleva(kk),elevb(kk)',eleva(kk),elevb(kk)
             ! END DO
             ! STOP 'htint'
             EXIT ! NP: Exit instead of Stop
          end if
       end do
     end do


end subroutine htint
!=================================================================================================================


!=================================================================================================================
subroutine get_fire_properties(burnt_area,heat_flux_max,heat_flux_min,iveg_ag, imm)

   !*******************************************************************************
    implicit none
   !*******************************************************************************
   
    real(kind=kind_phys), intent(in) :: &
         burnt_area,    & ! burnt area read in by the host model
         heat_flux_min, & ! minimal heatflux read in by the host model
         heat_flux_max    ! maximum heatflux read in by the host model

    integer, intent(in) :: iveg_ag, imm  
    
    integer :: &
         moist,   & ! fuel moisture
         i,       & ! loop counter ntime
         icount     ! loop counter number of seconds in the burn
    
    real(kind=kind_phys) :: &
         bfract,        & ! combustion factor
         effload,       & ! patchy burning NP: not used at all
         heat,          & ! joules/kg - tropical forest (mt) NP: there are other alternatives in the comments but only this one is used
         hinc             ! heat increase for a 5 minute ramp up in 4 steps
    
    real(kind=kind_phys), parameter :: beta = 0.88  !ref.: Paugam  et al., 2015
     
    real(kind=kind_phys),  dimension(2,4) :: heat_flux
         data heat_flux/  &
         !---------------------------------------------------------------------
         !  heat flux         !IGBP Land Cover      ! 
         !   kW/m^2           !Legend and           ! reference
         ! min   max          !description          ! 
         !--------------------------------------------------------------------
           30.0,  80.0,   &! iveg_ag=1:Tropical Forest         ! igbp 2
           30.0,  80.0,   &! iveg_ag=2:Boreal forest           ! igbp 1 & 3,4,5
            4.4,  23.0,   &! iveg_ag=3:cerrado/woody savanna   | igbp  6 thru 9, and 11
            3.3,   3.3    /! iveg_ag=4:Grassland/cropland      ! igbp 10, 12, 14
         !--------------------------------------------------------------------

   call mpas_log_write('', messageType=MPAS_LOG_ERR)
   call mpas_log_write('--- enter subroutine get_fire_properties:', messageType=MPAS_LOG_ERR)

   !*******************************************************************************
   ! Fire data inputs
   !*******************************************************************************

    !-- fire area at the surface
    area = burnt_area   ! area of burn, m^2
    

    !-- heat flux for the vegetation type
    if ( heat_flag == 1) then
         heat_fluxW = heat_flux(imm,iveg_ag) * 1000. !take the correspondent value in the table and convert to W/m^2
      
      !uncomment when have FRP
      !elseif ( heat_flag == 2) then
         ! "beta" factor converts FRP to convective energy
         !heat_fluxW = beta*(FRP/area)/0.55 ! in W/m^2
    endif

    

    !-- specific heat joules per kg of fuel consumed  
    !Why use the same value for the all vegetation since we have differents values?
    !heat = 21.e6     !-                  
    !heat = 15.5e6    !- joules/kg - cerrado/savannah
    heat = 19.3e6     !- joules/kg - tropical forest (mt)


    !-- others fires parameters
    mdur = 53        ! duration of burn, minutes
    !bload = 10.      ! total loading, kg/m**2 : not used at all 
    moist = 10       ! fuel moisture, %. average fuel moisture,percent dry
    maxtime =mdur-1  ! model time, min


   !*******************************************************************************
   ! Fix up inputs
   !*******************************************************************************

    if (MOD (maxtime, 2) .ne.0) maxtime = maxtime+1  !make maxtime even
    maxtime = maxtime * 60    ! and put in seconds

    rsurf = sqrt (area / pi)  !- entrainment surface radius (m)

    fmoist   = moist / 100.   !- fuel moisture fraction (%)


    !------------------------------------------------------------------------------
    ! calculate the energy flux and water content at lboundary.
    !----------------------------------------------------------------------------

     !fills heating() on a minute basis. 
     !could ask for a file at this po in the program. 
     !whatever is input has to be adjusted to a one minute timescale.

     do i = 1, ntime            !- make sure of energy release
       heating (i) = 0.0001     !- avoid possible divide by 0
     enddo

     tdur = mdur * 60.          !- number of seconds in the burn

     !bfract = 1.                !- combustion factor: not used at all

     !effload = bload * bfract   !- patchy burning: not used at all


    !------------------------------------------------------------------------------
    ! spread the burning evenly over the interval of the burn duration
    !----------------------------------------------------------------------------
     ! for the first few minutes (5min) the convective energy is lower for stability

     icount = 1  
     if(mdur > ntime) STOP 'Increase time duration (ntime) in min - see file "plumerise_mod.f90"'

     do while (icount.le.mdur)                             
          heating (icount) = heat_fluxW  * 0.55     ! W/m**2 (0.55 converts to convective energy)
          icount = icount + 1  
     enddo
    


    !BRAMS has a different procedure for imm=2. Check with Saulo
    !Ramp for 5 minutes
    !The convective energy increases linearly from 0.1 to heat_fluxW*0.55 in the first 5min. 
    hinc = heating (1) / 4.  
    heating (1) = 0.1              !convective energy = 0.1
    heating (2) = hinc             !convective energy = (heat_fluxW  * 0.55)*0.25
    heating (3) = 2. * hinc        !convective energy = (heat_fluxW  * 0.55)*0.50
    heating (4) = 3. * hinc        !convective energy = (heat_fluxW  * 0.55)*0.75


    call mpas_log_write('area of burn (m2)=$r', realArgs=(/area/),messageType=MPAS_LOG_ERR)
    call mpas_log_write('entrainment surface radius (m)=$r', realArgs=(/rsurf/),messageType=MPAS_LOG_ERR)  
    call mpas_log_write('duration of burn (min)=$i', intArgs=(/mdur/),messageType=MPAS_LOG_ERR)
    call mpas_log_write('specific heat (J/kg)=$r', realArgs=(/heat/),messageType=MPAS_LOG_ERR)
    call mpas_log_write('fuel moisture (%)=$i', intArgs=(/moist/),messageType=MPAS_LOG_ERR)
    call mpas_log_write('imm=$i', intArgs=(/imm/),messageType=MPAS_LOG_ERR)
    call mpas_log_write('iveg_ag=$i', intArgs=(/iveg_ag/),messageType=MPAS_LOG_ERR)
    call mpas_log_write('heat_fluxW (W/m2)=$r', realArgs=(/heat_fluxW/),messageType=MPAS_LOG_ERR)
    call mpas_log_write('convective energy flux (W/m2)=$r', realArgs=(/heating(6)/),messageType=MPAS_LOG_ERR)
    call mpas_log_write('', messageType=MPAS_LOG_ERR)


 return !I need this here because of the loop in imm?
end subroutine get_fire_properties
!=================================================================================================================



end module plumerise