!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module inject
!
! None
!
! :References: None
!
! :Owner: David Liptai
!
! :Runtime parameters:
!   - mdot          : *mass injection rate in grams/second*
!   - mdot_type     : *injection rate (0=const, 1=cos(t), 2=r^(-2))*
!   - npartperorbit : *particle injection rate in particles/binary orbit*
!
! :Dependencies: binaryutils, externalforces, infile_utils, io, options,
!   part, partinject, physcon, random, units
!
 use io, only:error
 use physcon, only:pi
 implicit none
 character(len=*), parameter, public :: inject_type = 'asteroidwind'
 real, public          :: mdot        = 6.33E+17     ! mass injection rate in grams/second
 real,save    :: dndt_scaling             ! scaling to get ninject correct

 public :: init_inject,inject_particles,write_options_inject,read_options_inject

 private

 real         :: npartperorbit = 100.     ! particle injection rate in particles per orbit
 integer      :: mdot_type     = 0        ! injection rate (0=const, 1=cos(t), 2=r^(-2))
 logical,save :: scaling_set              ! has the scaling been set (initially false)
 integer      :: shuffle_count = 0

contains
!-----------------------------------------------------------------------
!+
!  Initialize global variables or arrays needed for injection routine
!+
!-----------------------------------------------------------------------
subroutine init_inject(ierr)
 integer, intent(inout) :: ierr

 scaling_set = .false.

 ierr = 0

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  Inject particles
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npartoftype,dtinject)
 use io,            only:fatal
 use part,          only:nptmass,massoftype,igas,hfact,ihsoft,shuffle_part,isdead_or_accreted,kill_particle
 use partinject,    only:add_or_update_particle
 use physcon,       only:twopi,gg,kboltz,mass_proton_cgs
 use random,        only:get_random_pos_on_sphere,ran2
 use units,         only:umass, utime
 use options,       only:iexternalforce
 use externalforces,only:mass1
 use binaryutils,   only:get_orbit_bits
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 real,    dimension(3)  :: xyz,vxyz,r1,r2,v2,v1
 integer :: i,ipart,npinject,seed,pt
 real    :: dmdt,radd,h,u,speed,inject_this_step
 real    :: m1,m2,r
 real    :: dt,phi
 real, save :: have_injected,t_old
 real, save :: semia
 real ::tprec,beta,inc,rx,ry,rz,vx,vy,vz
 if (nptmass < 2 .and. iexternalforce == 0) call fatal('inject_asteroidwind','not enough point masses for asteroid wind injection')
 if (nptmass > 2) call fatal('inject_asteroidwind','too many point masses for asteroid wind injection')

 if (nptmass == 2) then
    pt = 1
    r1 = xyzmh_ptmass(1:3,2)
    m1 = xyzmh_ptmass(4,2)
    v1 = vxyz_ptmass(1:3,2)
 else
    pt = 1
    r1 = 0.
    m1 = mass1
    v1 = 0.
 endif

 r2        = xyzmh_ptmass(1:3,pt)
! rasteroid = xyzmh_ptmass(ihsoft,pt)
! radd = 0.0024*10.0D0**0.3333 ! mp=0.01
 radd = 10.0
!radd=0.015
 m2        = xyzmh_ptmass(4,pt)
 v2        = vxyz_ptmass(1:3,pt)
 speed     = sqrt(dot_product(v2,v2))
! vhat      = v2/speed
 r         = sqrt(dot_product(r1-r2,r1-r2))

!
! Add any dependency on radius to mass injection rate (and convert to code units)
!
 dmdt      = mdot*mdot_func(r,semia)/(umass/utime) ! Use semi-major axis as r_ref

!-- How many particles do we need to inject?
!   (Seems to need at least eight gas particles to not crash) <-- This statement may or may not be true...
!
 if (npartoftype(igas)<8) then
    npinject = 8-npartoftype(igas)
 else
!    print*,massoftype(igas)
    ! Calculate how many extra particles from previous step to now
    dt = time - t_old
    inject_this_step = dt*mdot/massoftype(igas)/(umass/utime)
    npinject = max(0, int(0.5 + have_injected + inject_this_step - npartoftype(igas) ))

    ! Save for next step (faster than integrating the whole thing each time)
    t_old = time
    have_injected = have_injected + inject_this_step
 endif
!
! clean particles 
!
if (shuffle_count > 100) then 
  write(*,*) "SHUFFLE: time is ", time, "dtlast is", dtlast
  do i=1,npart
   if (isdead_or_accreted(xyzh(4,i))) then
    call kill_particle(i)
   endif
  enddo
  call shuffle_part(npart) 
  shuffle_count = 0 
 else 
  shuffle_count = shuffle_count + 1 
 endif 

!
!-- Randomly inject particles around the radius 'radd'
!
 tprec=-45749.0
 inc=0.0D0/180D0*pi
 beta=2.0*pi*time/tprec
 do i=1,npinject
    phi = 2.*pi*(ran2(seed) - 0.5)
    !    xyz       = r2 +radd*(/cos(phi),sin(phi),0./)
    !    xyz       = r2 +radd*(/cos(phi)*cos(10.0D0/180.D0*pi),sin(phi),cos(phi)*sin(10.0D0/180.0D0*pi)/)
    !    vxyz      = v2+sqrt(m2/radd)*(/-sin(phi)*cos(10.0D0/180D0*pi),cos(phi),-sin(phi)*sin(10.0D0/180.0D0*pi)/)
    rx = cos(beta)*cos(phi)*cos(inc)-sin(beta)*sin(phi)
    ry= sin(beta)*cos(phi)*cos(inc)+cos(beta)*sin(phi)
    rz= cos(phi)*sin(inc)
    vx=-cos(beta)*sin(phi)*cos(inc)-sin(beta)*cos(phi)
    vy=-sin(beta)*sin(phi)*cos(inc)+cos(beta)*cos(phi)
    vz= -sin(phi)*sin(inc)
    xyz       = r2 +radd*(/rx,ry,rz/)
    vxyz      = v2+sqrt(m2/radd)*(/vx,vy,vz/)
!    vxyz      = v2+sqrt(m2/radd)*(/-sin(phi),cos(phi),0./)
    u         = 0. ! setup is isothermal so utherm is not stored
    h         = hfact*(radd/2.)
    ipart     = npart + 1
    call add_or_update_particle(igas,xyz,vxyz,h,u,ipart,npart,npartoftype,xyzh,vxyzu)
  enddo

 !
 !-- no constraint on timestep
 !
 dtinject = huge(dtinject)

end subroutine inject_particles

!-----------------------------------------------------------------------
!+
!  Returns dndt(t) depending on which function is chosen
!  Note that time in this function is strictly the fraction
!  of the orbit, not absolute time
!+
!-----------------------------------------------------------------------

real function mdot_func(r,r_ref)
 real, intent(in) :: r,r_ref

 select case (mdot_type)
 case (2)
    mdot_func = (r_ref/r)**2
 case default
    mdot_func = 1.0
 end select

end function mdot_func

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file.
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(mdot         ,'mdot'         ,'mass injection rate in grams/second'              ,iunit)
 call write_inopt(npartperorbit,'npartperorbit','particle injection rate in particles/binary orbit',iunit)
 call write_inopt(mdot_type    ,'mdot_type'    ,'injection rate (0=const, 1=cos(t), 2=r^(-2))'     ,iunit)

end subroutine write_options_inject

!-----------------------------------------------------------------------
!+
!  Reads input options from the input file.
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(name,valstring,imatch,igotall,ierr)
 use io, only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_options_inject'

 imatch  = .true.
 select case(trim(name))
 case('mdot')
    read(valstring,*,iostat=ierr) mdot
    ngot = ngot + 1
    if (mdot  <  0.) call fatal(label,'mdot < 0 in input options')
 case('npartperorbit')
    read(valstring,*,iostat=ierr) npartperorbit
    ngot = ngot + 1
    if (npartperorbit < 0.) call fatal(label,'npartperorbit < 0 in input options')
 case('mdot_type')
    read(valstring,*,iostat=ierr) mdot_type
    ngot = ngot + 1
 case default
    imatch = .false.
 end select

 igotall = (ngot >= 1)

end subroutine read_options_inject

end module inject
