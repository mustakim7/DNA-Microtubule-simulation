!
!    Flexible polymer
!    With self-propulsion at the end of polymer with polarity
!    Considered multiple DNAs of the same size within a square with PBC in x and y
!    Use Gaussian vol rep for iter < iter_soft to avoid overlapping beads, then start LJ
!
!---------------------------------------------------------------------------------------------------------------------------------!

 module cpu_mod
 save

 real*8, parameter ::  &
 ak_fene =   30.0000000000000      , &
 r_fene_sq =  0.444444444444444      , &
 eps_fene =  -6.66666666666667      , &
 !
 eps_LJ =   2.00000000000000      , trunc_volsq = 1.d0, &
 eps_LJ4 =   8.00000000000000      , &
 eps_LJ24 =   48.0000000000000      , &
 sigma_LJ =  0.296966239380113      , &
 sigma_LJ6 =  6.858710562414263E-004 , &
 sigma_LJ12 =  4.704191057897298E-007 , &
 !
 bead_size =  0.333333333333333      , &
 !
 eps_vol_rep =   20.0000000000000      , &
 alp_vol_rep =   5.00000000000000      , &
 fac_vol_rep =   200.000000000000      , &
 !
 fac_lang1 =  3.703703703703703E-002 , &
 fac_lang2 =  0.272165526975909      , &
 fac_lang_orient =   1.41421356237310      , &
 !
 delt =  1.000000000000000E-004 , &
 del_t_sqrt =  1.000000000000000E-002
 
 real*8, parameter ::  &
 pi =   3.14159265358979      , &
 pi4 =   12.5663706143592      , &
 pi4by3 =   4.18879020478639      , &
 stretch_init =  0.533333333333333     

 end module cpu_mod

!----------------------------------------------------------------------------------------------------------------------------------

 include 'mkl_vsl.f90'
!----------------------------------------------------------------------------------------------------------------------------------

 program main
 use cpu_mod
 use MKL_VSL_TYPE
 use MKL_VSL

 implicit integer*4(i-o)
 implicit real*8(a-h,p-z) 
 integer*4, allocatable, dimension(:,:) :: nbox_DNA
 integer*4, allocatable, dimension(:) :: nbx_DNA_updt_tag, nbx_DNA_updt_frm, nbx_DNA_updt_to !, itype_DNA
 real*8, allocatable, dimension(:) :: x_DNA, y_DNA, x1_DNA, y1_DNA, eta, det_DNA_fx, det_DNA_fy, &
                                      velx_DNA, vely_DNA, aMSD_DNA, delx_DNA, dely_DNA, theta, phi                                  
                                      
 dimension i_neigh(2)

 TYPE (VSL_STREAM_STATE) :: stream_G!, stream_U
 integer i_brng_G,           i_seed, method_G,           ndim_eta, nreplica
 !integer i_brng_G, i_brng_U, i_seed, method_G, method_U, ndim_eta, nreplica

 i_brng_G = VSL_BRNG_MT19937 !NIEDERR!
 !i_brng_U = VSL_BRNG_MT19937 !NIEDERR!
 method_G = VSL_RNG_METHOD_GAUSSIAN_ICDF !VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2
 !method_U = VSL_RNG_METHOD_UNIFORM_STD
 i_seed = 786

 open(11,file='system_size.inp')
 open(12,file='replica.inp')
 open(14,file='iter.inp')
 open(15,file='write.inp')
 open(13,file='input_pos.inp')
 open(16,file='control_para.inp')

 !open(111,file='DNA_bead_type.dat')
 open(112,file='small_delt.dat')
 open(113,file='time_series.dat')
 open(114,file='early_MSD_DNA.dat')
 open(115,file='MSD_DNA.dat')

 read(11,*) nDNA, nbead_per_DNA, len_box
 read(12,*) nreplica
 read(14,*) iter_tot, iter_ther, iter_soft
 read(15,*) ifreq_snap, ifreq_time_ser, ifreq_snap1, istore_i 
 read(16,*) spp, aniso_fric_tensor 

 close(11) ; close(12) ; close(14) ; close(15) 

 do ii = 1, (nreplica-1)*547
   xxx = grnd() 
 enddo

 nbead_tot_DNA = nDNA * nbead_per_DNA
 nbox_dim_DNA = nbead_tot_DNA                      ! Carefully update

 ndim_sys = 2

 ndim_eta = ndim_sys*nbead_tot_DNA !+ nDNA          ! bcz., nDNA leading beads have self-propulsion dirn.

 ndim_MSD = 9 * (floor(log10(1.0*iter_tot)) + 1)

 alen_box = dfloat(len_box)
 alen_box_half = alen_box * 0.5d0
 alen_box2 = alen_box * 2.d0
 len_box_sq = len_box*len_box
 len_box2p = 2*len_box + 1

 allocate( x_DNA(nbead_tot_DNA), y_DNA(nbead_tot_DNA), nbox_DNA(len_box_sq, nbox_dim_DNA),           & 
           x1_DNA(nbead_tot_DNA), y1_DNA(nbead_tot_DNA), eta(ndim_eta),                              & !itype_DNA(nbead_tot_DNA), 
           det_DNA_fx(nbead_tot_DNA), det_DNA_fy(nbead_tot_DNA), nbx_DNA_updt_tag(nbead_tot_DNA),    &
           nbx_DNA_updt_frm(nbead_tot_DNA), nbx_DNA_updt_to(nbead_tot_DNA),                          & 
           velx_DNA(nbead_tot_DNA), vely_DNA(nbead_tot_DNA), theta(nDNA),phi(nDNA) )

 allocate( aMSD_DNA(ndim_MSD), delx_DNA(nbead_tot_DNA), dely_DNA(nbead_tot_DNA) )
!read here the dna_pos fil (ibead1,nbead_tot)

 do ibead = 1, nbead_tot_DNA
   read(13,*) X_DNA(ibead), Y_DNA(ibead)
 enddo
close(13)

 call initial_place_DNA(nbead_tot_DNA, nbead_per_DNA, len_box, len_box_sq, x_DNA, y_DNA, nbox_DNA, nbox_dim_DNA) 


 !call type_set_DNA(nbead_tot_DNA, den_TDP_type, den_OTHER_type, itype_DNA, nTDPtype, nOTHERtype)
 !do ibead = 1, nbead_tot_DNA
 !  write(111,*) ibead, itype_DNA(ibead)
 !enddo
 !close(111)

 i_seed = i_seed + nreplica - 1
 ierr_G = vslnewstream( stream_G, i_brng_G, i_seed )
 !ierr_U = vslnewstream( stream_U, i_brng_U, i_seed )

 n_small_delt = 0

 vel_max = -10.d0

 delx_DNA = 0.d0 ; dely_DNA = 0.d0   ! initialization for MSD                                                     
 aMSD_DNA = 0.d0 ; icount_MSD = 0    !                                      

 call CPU_TIME(t1)


! do ibead = 1, nbead_tot_DNA
!     theta(ibead/nbead_per_DNA)=(2*grnd()-1)*pi    !initial random theta value 
! end do


  ierr_G = vdrnggaussian( method_G, stream_G, nDNA, theta, 0.d0, 1.d0 )  !initial random theta

  do iter = 1, iter_tot                                                                             ! start of iter loop

  ierr_G = vdrnggaussian( method_G, stream_G, nDNA, phi, 0.d0, 0.00010d0 )  !initial random theta
   
   ene = 0.d0
   nbx_DNA_updt_tag = 0
   itag_fene_wall = 0

   ierr_G = vdrnggaussian( method_G, stream_G, ndim_eta, eta, 0.d0, 1.d0 )

   !----------------------------------------------------------------------------------------------------------------------------------

   !$OMP PARALLEL DEFAULT(SHARED)
   !$OMP DO PRIVATE(ixref,iyref,fvolx,fvoly,ixx,iyy,ix,iy,icell,nn,jbead,delx,dely,r_sq,ff,ee,rr)  &
   !$OMP& PRIVATE(ffenex,ffeney,i_neigh,iDNA, i_dummy,icell1)                                      &
   !$OMP& PRIVATE(ffx,ffy)                                                                         &
   !$OMP& REDUCTION(+:ene)
   
   do ibead = 1, nbead_tot_DNA                                                                     ! start of ibead loop 
     
     ixref = ceiling(x_DNA(ibead)) ; iyref = ceiling(y_DNA(ibead)) 
                                                                                                   
     fvolx = 0.d0; fvoly = 0.d0

     do ixx = ixref-1, ixref+1                                                                     !
       ix = ixx                                                                                    !
       if(ix < 1) ix = len_box                                                                     ! PBC
       if(ix > len_box) ix = 1                                                                     !
                                                                                                   !
       do iyy = iyref-1, iyref+1                                                                   !
         iy = iyy                                                                                  !
         if(iy < 1) iy = len_box                                                                   ! PBC
         if(iy > len_box) iy = 1                                                                   !
                                                                                                   !
         icell = (iy-1)*len_box + ix                                                               !
                                                                                                   !
         do nn = 1, nbox_DNA(icell,1)                                                              !
           jbead = nbox_DNA(icell,nn+1)                                                            !
           if(jbead /= ibead) then                                                                 !
             delx = x_DNA(ibead) - x_DNA(jbead)                                                    !
             dely = y_DNA(ibead) - y_DNA(jbead)                                                    !
                                                                                                   !
             if(delx >  alen_box_half) delx = delx - alen_box                                      !
             if(dely >  alen_box_half) dely = dely - alen_box                                      !
             if(delx < -alen_box_half) delx = alen_box + delx                                      ! PBC
             if(dely < -alen_box_half) dely = alen_box + dely                                      !
                                                                                                   !
             r_sq = delx*delx + dely*dely                                                          !
                                                                                                   !
             if(r_sq < trunc_volsq) then                                                           !
                                                                                                   !
               if(iter < iter_soft) then                                                           !
                 ff = force_vol_rep_DNA(r_sq)                                                      !
                 ! ff = (1.0*iter)/iter_soft * force_vol_rep_DNA(r_sq)                                                      !
                 fvolx = fvolx + delx*ff                                                           !
                 fvoly = fvoly + dely*ff                                                           ! Soft potential to avoid overlap -> otherwise LJ would diverge
               !  fvolz = fvolz + delz*ff                                                           !
                 ee = energy_vol_rep_DNA(r_sq) * 0.5d0                                             !
                 ene = ene + ee                                                                    !
               else                                                                                !
                 ff = force_LJ(r_sq)                                                               !
                 fvolx = fvolx + delx*ff                                                           ! LJ: have both repul. & attr.
                 fvoly = fvoly + dely*ff                                                           !
                 ee = energy_LJ(r_sq) * 0.5d0                                                      !
                 ene = ene + ee                                                                    !
               endif                                                                               !
                                                                                                   !
             endif                                                                                 !
           endif                                                                                   !
         enddo                                                                                     !
       enddo                                                                                       !
     enddo                                                                                         !

     ffenex = 0.d0; ffeney = 0.d0                                                                  !
     i_neigh(1) = ibead - 1                                                                        !
     i_neigh(2) = ibead + 1                                                                        !
     if(mod(ibead, nbead_per_DNA) == 1) i_neigh(1) = ibead                                         !
     if(mod(ibead, nbead_per_DNA) == 0) i_neigh(2) = ibead                                         !
     do ii = 1, 2                                                                                  !
       delx = x_DNA(ibead) - x_DNA(i_neigh(ii))                                                    !
       dely = y_DNA(ibead) - y_DNA(i_neigh(ii))                                                    ! chain connectivity among DNA beads
                                                                                                   !
       if(delx >  alen_box_half) delx = delx - alen_box                                            !
       if(dely >  alen_box_half) dely = dely - alen_box                                            ! PBC
       if(delx < -alen_box_half) delx = alen_box + delx                                            !
       if(dely < -alen_box_half) dely = alen_box + dely                                            !
                                                                                                   !
       r_sq = delx*delx + dely*dely                                                                !
       ff = force_fene(r_sq)                                                                       !
       ffenex = ffenex + delx*ff                                                                   !
       ffeney = ffeney + dely*ff                                                                   !
     enddo                                                                                         !
     ee = energy_fene(r_sq)                                                                        !
     ene = ene + ee                                                                                !
     
    
     det_DNA_fx(ibead) = fac_lang1*(fvolx+ffenex)                                                  !
     det_DNA_fy(ibead) = fac_lang1*(fvoly+ffeney)                                                  ! overdamped Langevin update
     !! Last bead 
      !  write(222,*) mod(ibead,nbead_per_DNA)                                                                                      !
     if (mod(ibead,nbead_per_DNA)==0) then  
!       delx = x_DNA(ibead) - x_DNA(ibead-1)                                                    !
!       dely = y_DNA(ibead) - y_DNA(ibead-1)                                                    ! chain connectivity in last two DNA beads(due to SPP)
                                                                                                   !
!       if(delx >  alen_box_half) delx = delx - alen_box                                            !
!       if(dely >  alen_box_half) dely = dely - alen_box                                            ! PBC
!       if(delx < -alen_box_half) delx = alen_box + delx                                            !
!       if(dely < -alen_box_half) dely = alen_box + dely                                            !
                                                                                                   !
!       r_mod = sqrt(delx*delx + dely*dely)
       
!       x1_DNA(ibead) = x1_DNA(ibead) + delt*1.00d0*delx/r_mod              
!       y1_DNA(ibead) = y1_DNA(ibead) + delt*1.00d0*dely/r_mod
        F=spp
        rr=aniso_fric_tensor   
        xx1=F*cos(theta(ibead/nbead_per_DNA)) + det_DNA_fx(ibead)
        yy1=F*sin(theta(ibead/nbead_per_DNA)) + det_DNA_fy(ibead)
          rr1= (rr+1.0d0)/2.0d0    !rr=psi_parallel/psi_perpendicular
          rr2= (rr-1.0d0)/2.0d0  
       x1_DNA(ibead) = x1_DNA(ibead) + (rr1-rr2*cos(2.0d0*theta(ibead/nbead_per_DNA)))*delt*xx1 + &
                        rr2*sin(2.0d0*theta(ibead/nbead_per_DNA))*delt*yy1         
       y1_DNA(ibead) = y1_DNA(ibead) +  rr2*sin(2.0d0*theta(ibead/nbead_per_DNA))*delt*xx1 + &
                         (rr1+rr2*cos(2.0d0*theta(ibead/nbead_per_DNA)))*delt*yy1          

!       x1_DNA(ibead) = x1_DNA(ibead) + delt*1.00d0*cos(theta)     !         
!       y1_DNA(ibead) = y1_DNA(ibead) + delt*1.00d0*sin(theta)      !         
!        xx=del_t_sqrt*0.000010d0*eta(i_dummy) 
!        write(300,*)xx
!        write(301,*)xx1
!        write(302,*)theta(ibead/nbead_per_DNA)
!update theta here, initialize  theta randomly, 
!      theta(ibead/nbead_per_DNA)=theta(ibead/nbead_per_DNA) +xx ! del_t_sqrt*1.0d0*eta (i_dummy) 
      theta(ibead/nbead_per_DNA)=theta(ibead/nbead_per_DNA) + phi(ibead/nbead_per_DNA)* del_t_sqrt
!       theta = theta + del_t_sqrt*1.00d0* eta(i_dummy) 

!    write(223,*) iter, r_mod                                                         !
     
    !! The other beads
    else
     i_dummy = (ibead - 1)*ndim_sys                                                                !
     x1_DNA(ibead) = x_DNA(ibead) + delt*det_DNA_fx(ibead) + del_t_sqrt*fac_lang2*eta(i_dummy)     !         
                                                                                                   ! 
     i_dummy = i_dummy + 1                                                                         ! 
     y1_DNA(ibead) = y_DNA(ibead) + delt*det_DNA_fy(ibead) + del_t_sqrt*fac_lang2*eta(i_dummy)     !         
 
    endif
    !  write(221,*) mod(ibead,nbead_per_DNA)
!-----------------------------------  
!      190 format(' ',2F25.12,2I10)
!-------------------------------------
!     i_dummy = (ibead - 1)*ndim_sys                                                                !
!     det_DNA_fx(ibead) = fac_lang1*(fvolx+ffenex)                                                  !
!     if (mod(nbead_tot_DNA,nbead_per_DNA)==0) then  
!     r_mod = sqrt ((x_DNA(ibead)-x_DNA(ibead-1))**2 + (y_DNA(ibead)-y_DNA(ibead-1))**2) 
!     x1_DNA(ibead) = x_DNA(ibead) + delt*det_DNA_fx(ibead) + del_t_sqrt*fac_lang2*eta(i_dummy) + &
!                                                 0.00001d0*(x_DNA(ibead)-x_DNA(ibead-1))/r_mod     ! 
!     else     
!     x1_DNA(ibead) = x_DNA(ibead) + delt*det_DNA_fx(ibead) + del_t_sqrt*fac_lang2*eta(i_dummy)     !         
!     endif

!     i_dummy = i_dummy + 1                                                                         ! 
!     det_DNA_fy(ibead) = fac_lang1*(fvoly+ffeney)                                                  ! overdamped Langevin update
!     if (mod(nbead_tot_DNA,nbead_per_DNA)==0) then 
!     r_mod = sqrt ((x_DNA(ibead)-x_DNA(ibead-1))**2 + (y_DNA(ibead)-y_DNA(ibead-1))**2) 
!     y1_DNA(ibead) = y_DNA(ibead) + delt*det_DNA_fy(ibead) + del_t_sqrt*fac_lang2*eta(i_dummy) + &
!                                                 0.00001d0*(y_DNA(ibead)-y_DNA(ibead-1))/r_mod     ! 
!     else     
!     y1_DNA(ibead) = y_DNA(ibead) + delt*det_DNA_fy(ibead) + del_t_sqrt*fac_lang2*eta(i_dummy)     !         
!     endif
                                                                                                   !
     if(x1_DNA(ibead) > alen_box) x1_DNA(ibead) = x1_DNA(ibead) - alen_box                         ! 
     if(y1_DNA(ibead) > alen_box) y1_DNA(ibead) = y1_DNA(ibead) - alen_box                         ! PBC
     if(x1_DNA(ibead) < 0.d0)     x1_DNA(ibead) = x1_DNA(ibead) + alen_box                         !
     if(y1_DNA(ibead) < 0.d0)     y1_DNA(ibead) = y1_DNA(ibead) + alen_box                         !

     ix = ceiling(x_DNA(ibead)); ix1 = ceiling(x1_DNA(ibead))                                      !
     iy = ceiling(y_DNA(ibead)); iy1 = ceiling(y1_DNA(ibead))                                      !
     icell  = (iy-1) *len_box + ix                                                                 !
     icell1 = (iy1-1)*len_box + ix1                                                                !
                                                                                                   ! tagging nbox_DNA to be updated
     if(icell /= icell1) then                                                                      !
       nbx_DNA_updt_tag(ibead) = 1                                                                 !
       nbx_DNA_updt_frm(ibead) = icell                                                             !
       nbx_DNA_updt_to(ibead)  = icell1                                                            !
     endif                                                                                         !

   enddo                                                                                           ! end of ibead loop
   !$OMP END DO
   !$OMP END PARALLEL

   !----------------------------------------------------------------------------------------------------------------------------------

   !$OMP PARALLEL DEFAULT(SHARED)                                                                  !
   !$OMP DO PRIVATE(delx,dely,jbead)                                                               !
   do ibead = 1, nbead_tot_DNA-1                                                                   !
     jbead = ibead + 1                                                                             !
     if(mod(ibead, nbead_per_DNA) == 0) jbead = ibead                                              !
     delx = abs(x1_DNA(ibead) - x1_DNA(jbead))                                                     !
     dely = abs(y1_DNA(ibead) - y1_DNA(jbead))                                                     !
                                                                                                   ! no need of sign of del
     if(delx > alen_box_half) delx = alen_box - delx                                               !
     if(dely > alen_box_half) dely = alen_box - dely                                               ! PBC
                                                                                                   !
     if(delx*delx + dely*dely > r_fene_sq) then                                                    !
       !$OMP ATOMIC WRITE                                                                          !
       itag_fene_wall = 1                                                                          ! check FENE
       !write(20000+iter,173) ibead, delx, dely,       sqrt(delx*delx + dely*dely            )     !
     endif                                                                                         !
   enddo                                                                                           !
   !$OMP END DO                                                                                    !
   !$OMP END PARALLEL                                                                              !

   !----------------------------------------------------------------------------------------------------------------------------------

   if(itag_fene_wall > 0) then                                                                            !
     nbx_DNA_updt_tag = 0                                                                                 !
     n_small_delt = n_small_delt + 1                                                                      !
                                                                                                          !
     !$OMP PARALLEL DEFAULT(SHARED)                                                                       !
     !$OMP DO PRIVATE(i_dummy,ix,iy,ix1,iy1,icell,icell1)                                                 !
     do ibead = 1, nbead_tot_DNA                                                                          !
       i_dummy = (ibead - 1)*ndim_sys + 1                                                                 !
       x1_DNA(ibead) = x_DNA(ibead) + delt*0.000001d0*det_DNA_fx(ibead) + del_t_sqrt*0.001d0*eta(i_dummy) !                       
                                                                                                          !
       i_dummy = i_dummy + 1                                                                              !
       y1_DNA(ibead) = y_DNA(ibead) + delt*0.000001d0*det_DNA_fy(ibead) + del_t_sqrt*0.001d0*eta(i_dummy) !                       
                                                                                                          !
     if (mod(ibead,nbead_per_DNA)==0) then  
       delx = x1_DNA(ibead) - x1_DNA(ibead-1)                                                    !
       dely = y1_DNA(ibead) - y1_DNA(ibead-1)                                                    ! chain connectivity in last two DNA beads due to SPP
                                                                                                   !
       if(delx >  alen_box_half) delx = delx - alen_box                                            !
       if(dely >  alen_box_half) dely = dely - alen_box                                            ! PBC
       if(delx < -alen_box_half) delx = alen_box + delx                                            !
       if(dely < -alen_box_half) dely = alen_box + dely                                            !
                                                                                                   !
       r_mod = sqrt(delx*delx + dely*dely)

!       x1_DNA(ibead) = x1_DNA(ibead) + delt*0.000001d0*1.00d0*delx/r_mod     !         
!       y1_DNA(ibead) = y1_DNA(ibead) + delt*0.000001d0*1.00d0*dely/r_mod      !         
       
       x1_DNA(ibead) = x1_DNA(ibead) + delt*0.000001d0*1.00d0*cos(theta(ibead/nbead_per_DNA))     !         
       y1_DNA(ibead) = y1_DNA(ibead) + delt*0.000001d0*1.00d0*sin(theta(ibead/nbead_per_DNA))      !         
       xx=del_t_sqrt*1.0d0*eta(i_dummy)
!       write(*,*) xx
       theta(ibead/nbead_per_DNA)=theta(ibead/nbead_per_DNA) +xx ! del_t_sqrt*1.0d0*eta (i_dummy) 
     
endif

       if(x1_DNA(ibead) > alen_box) x1_DNA(ibead) = x1_DNA(ibead) - alen_box   !                          ! repeat overdamped Langevin update
       if(y1_DNA(ibead) > alen_box) y1_DNA(ibead) = y1_DNA(ibead) - alen_box   !                          !
       if(x1_DNA(ibead) < 0.d0)     x1_DNA(ibead) = x1_DNA(ibead) + alen_box   !                          !
       if(y1_DNA(ibead) < 0.d0)     y1_DNA(ibead) = y1_DNA(ibead) + alen_box   !                          !
                                                                                                          !
       ix = ceiling(x_DNA(ibead)); ix1 = ceiling(x1_DNA(ibead))         !                                 !
       iy = ceiling(y_DNA(ibead)); iy1 = ceiling(y1_DNA(ibead))         !                                 !
       icell  = (iy-1) *len_box + ix                                    !                                 !
       icell1 = (iy1-1)*len_box + ix1                                   !                                 !
                                                                        ! tagging nbox_DNA to updated     !                 
       if(icell /= icell1) then                                         !                                 !
         nbx_DNA_updt_tag(ibead) = 1                                    !                                 !
         nbx_DNA_updt_frm(ibead) = icell                                !                                 !
         nbx_DNA_updt_to(ibead)  = icell1                               !                                 !
       endif                                                            !                                 !
     enddo                                                                                                !
     !$OMP END DO                                                                                         !
     !$OMP END PARALLEL                                                                                   !
     write(112,*) iter, itag_fene_wall                                                                    !
     write(*,*) iter, itag_fene_wall                                                                    !
   endif                                                                                                  !

   !----------------------------------------------------------------------------------------------------------------------------------

   if(iter/ifreq_snap*ifreq_snap == iter) then                                                            !
                                                                                                          !
     if(iter/ifreq_time_ser*ifreq_time_ser == iter) then                                                  !
       write(113,173) iter, ene                                                                           !
     endif                                                                                                !
                                                                                                          !
     !$OMP PARALLEL DEFAULT(SHARED)                                                                       !
     !$OMP DO                                                                                             !
     do ibead = 1, nbead_tot_DNA                                                                          !
       velx_DNA(ibead) = x1_DNA(ibead) - x_DNA(ibead)                                                     !
       vely_DNA(ibead) = y1_DNA(ibead) - y_DNA(ibead)                                                     !
                                                                                                          !
       if(velx_DNA(ibead) >  alen_box_half) velx_DNA(ibead) = velx_DNA(ibead) - alen_box                  !
       if(vely_DNA(ibead) >  alen_box_half) vely_DNA(ibead) = vely_DNA(ibead) - alen_box                  !
       if(velx_DNA(ibead) < -alen_box_half) velx_DNA(ibead) = alen_box + velx_DNA(ibead)                  !
       if(vely_DNA(ibead) < -alen_box_half) vely_DNA(ibead) = alen_box + vely_DNA(ibead)                  !
                                                                                                          !
       !vel = velx_DNA(ibead)*velx_DNA(ibead) + vely_DNA(ibead)*vely_DNA(ibead)                            !
       !if(vel > vel_max) then                                                                             !
       !  vel_max = vel                                                                                    ! Correct only for serial execution
       !  iter_max = iter                                                                                  !
       !  ibead_max = ibead                                                                                !
       !endif                                                                                              !
     enddo                                                                                                !
     !$OMP END DO                                                                                         !
     !$OMP END PARALLEL                                                                                   ! velocity calculation and snapshot storing
                                                                                                          !
     do ibead = 1, nbead_tot_DNA                                                                          !
       write(1000000000+iter,172) x_DNA(ibead),      y_DNA(ibead),                          &             ! 
                                  det_DNA_fx(ibead), det_DNA_fy(ibead),                     &             !
                                  velx_DNA(ibead),   vely_DNA(ibead),                       & 
                          cos(theta(ibead/nbead_per_DNA)+1),  sin(theta(ibead/nbead_per_DNA)+1)               !
     enddo                                                                                                !
     close(1000000000+iter)                                                                               !
                                                                                                          !
   else if(iter>istore_i .and. iter/ifreq_snap1*ifreq_snap1==iter) then                                   !
                                                                                                          !
     do ibead = 1, nbead_tot_DNA                                                                          !
       write(1000000000+iter,172) x_DNA(ibead), y_DNA(ibead)                                              !
     enddo                                                                                                !
     close(1000000000+iter)                                                                               !
   endif                                                                                                  !

   !----------------------------------------------------------------------------------------------------------------------------------

   if(iter < iter_ther) then                                                                                           ! if(iter < iter_ther) then
                                                                                                                       !
     !$OMP PARALLEL DEFAULT(SHARED)                                                                                    !
     !$OMP DO PRIVATE(delx,dely)                                                                                       !
     do ibead = 1, nbead_tot_DNA                                                                                       !
       delx = x1_DNA(ibead) - x_DNA(ibead) ; dely = y1_DNA(ibead) - y_DNA(ibead)                                       !
                                                                                                                       ! calculating early-time MSD
       if(delx < -alen_box_half) delx = delx + alen_box ; if(delx > alen_box_half) delx = delx - alen_box              !
       if(dely < -alen_box_half) dely = dely + alen_box ; if(dely > alen_box_half) dely = dely - alen_box              !
                                                                                                                       !
       delx_DNA(ibead) = delx_DNA(ibead) + delx                                                                        !
       dely_DNA(ibead) = dely_DNA(ibead) + dely                                                                        !
     enddo                                                                                                             !
     !$OMP END DO                                                                                                      !
     !$OMP END PARALLEL                                                                                                !
                                                                                                                       !
     iter_sep = iter ; ii = 10**floor(log10(float(iter_sep)))                                                          !     
     if(mod(iter_sep,ii)==0) then                                                                                      !
       icount_MSD = icount_MSD + 1                                                                                     !
       do ibead = 1, nbead_tot_DNA                                                                                     !
         aMSD_DNA(icount_MSD) = aMSD_DNA(icount_MSD) + delx_DNA(ibead)*delx_DNA(ibead) +              &                !
                                   dely_DNA(ibead)*dely_DNA(ibead)                                                     !
       enddo                                                                                                           !
     endif                                                                                                             !
                                                                                                                       !
   else if(iter == iter_ther) then                                                                                     ! else if(iter == iter_ther) then
                                                                                                                       ! 
     if(nbead_tot_DNA > 0) aMSD_DNA = aMSD_DNA / nbead_tot_DNA                                                         !
                                                                                                                       !
     do ii = 1, icount_MSD                                                                                             !
       if(ii/9*9 /=ii ) then                                                                                           !
         jj = 10**(ii/9) * mod(ii,9)                                                                                   !
       else                                                                                                            !
         jj = 10**(ii/9 - 1) * 9                                                                                       !
       endif                                                                                                           ! writing early-time MSD
       write(114,173) jj, aMSD_DNA(ii)                                                                                 !           &                             
     enddo                                                                                                             ! initialization for s.s. MSD 
     close(114)                                                                                                        !
                                                                                                                       !
     delx_DNA = 0.d0 ; dely_DNA = 0.d0                                                                                 !
     aMSD_DNA = 0.d0 ; icount_MSD = 0                                                                                  !
                                                                                                                       !
   else                                                                                                                !! else
                                                                                                                       !
     !$OMP PARALLEL DEFAULT(SHARED)                                                                                    !
     !$OMP DO PRIVATE(delx,dely)                                                                                       !
     do ibead = 1, nbead_tot_DNA                                                                                       !
       delx = x1_DNA(ibead) - x_DNA(ibead) ; dely = y1_DNA(ibead) - y_DNA(ibead)                                       !
                                                                                                                       ! calculating steady-state MSD
       if(delx < -alen_box_half) delx = delx + alen_box ; if(delx > alen_box_half) delx = delx - alen_box              !
       if(dely < -alen_box_half) dely = dely + alen_box ; if(dely > alen_box_half) dely = dely - alen_box              !
                                                                                                                       !
       delx_DNA(ibead) = delx_DNA(ibead) + delx                                                                        !
       dely_DNA(ibead) = dely_DNA(ibead) + dely                                                                        !
     enddo                                                                                                             !
     !$OMP END DO                                                                                                      !
     !$OMP END PARALLEL                                                                                                !
                                                                                                                       !
     iter_sep = iter - iter_ther ; ii = 10**floor(log10(float(iter_sep)))                                              ! 
     if(mod(iter_sep,ii)==0) then                                                                                      !
       icount_MSD = icount_MSD + 1                                                                                     !
       do ibead = 1, nbead_tot_DNA                                                                                     !
         aMSD_DNA(icount_MSD) = aMSD_DNA(icount_MSD) + delx_DNA(ibead)*delx_DNA(ibead) +              &                !
                                   dely_DNA(ibead)*dely_DNA(ibead)                                                     !
       enddo                                                                                                           !
     endif                                                                                                             !
                                                                                                                       !
   endif                                                                                                               !! endif 

   !----------------------------------------------------------------------------------------------------------------------------------

   do ibead = 1, nbead_tot_DNA                                                                            !
     if(nbx_DNA_updt_tag(ibead) > 0) then                                                                 !
       n_old = nbox_DNA(nbx_DNA_updt_frm(ibead),1) + 1                                                    !
       nbox_DNA(nbx_DNA_updt_frm(ibead), 1) = nbox_DNA(nbx_DNA_updt_frm(ibead), 1) - 1                    !
       nbox_DNA(nbx_DNA_updt_to(ibead) , 1) = nbox_DNA(nbx_DNA_updt_to(ibead) , 1) + 1                    !
       nbox_DNA(nbx_DNA_updt_to(ibead) , nbox_DNA(nbx_DNA_updt_to(ibead),1)+1) = ibead                    !
                                                                                                          !
       i_spy = 0                                                                                          !
       do jj = 2, n_old                                                                                   ! nbox_DNA and pos_DNA update
         if(i_spy == 1) nbox_DNA(nbx_DNA_updt_frm(ibead), jj-1) = nbox_DNA(nbx_DNA_updt_frm(ibead),jj)    !
                                                                                                          !
         if(nbox_DNA(nbx_DNA_updt_frm(ibead), jj) == ibead) i_spy = 1                                     !
       enddo                                                                                              !
       nbox_DNA(nbx_DNA_updt_frm(ibead), n_old) = 0                                                       !
     endif                                                                                                !
   enddo                                                                                                  !
                                                                                                          !
   x_DNA = x1_DNA ; y_DNA = y1_DNA                                                                        !

   !----------------------------------------------------------------------------------------------------------------------------------
                                                  
   !---------------------------------- MSD : not for PBC --------------------------------------------------!
   !if(iter < iter_ther) then                                                                              ! if(iter < iter_ther) then
   !                                                                                                       !
   !  iter_sep = iter                                                                                      !
   !  ii = 10**floor(log10(float(iter_sep)))                                                               !
   !  if(mod(iter_sep,ii)==0) then                                                                         !
   !    icount_MSD = icount_MSD + 1                                                                        !
   !                                                                                                       !
   !    do iTDP = 1, nTDP                                                                                  !
   !      delx = x_TDP(iTDP) - x_TDP_0(iTDP)                                                               !
   !      dely = y_TDP(iTDP) - y_TDP_0(iTDP)                                                               !
   !      delz = z_TDP(iTDP) - z_TDP_0(iTDP)                                                               !
   !      aMSD_TDP(icount_MSD) = aMSD_TDP(icount_MSD) + delx*delx + dely*dely + delz*delz                  !
   !    enddo                                                                                              !
   !                                                                                                       !
   !    !$OMP PARALLEL DEFAULT(SHARED)                                                                     !
   !    !$OMP DO PRIVATE(ibead,jbead,delx,dely,delz,del)  &                                                ! calculating early-time MSD
   !    !$OMP& REDUCTION(+:aMSD_DNA_CM)                                                                    !
   !    do iDNA = 1, nDNA                                                                                  !
   !      ibead = (iDNA-1)*nbead_per_DNA                                                                   !
   !      jbead = ibead + nbead_per_DNA - 1                                                                !
   !                                                                                                       !
   !      xCM_DNA(iDNA) = sum( x_DNA(ibead:jbead) ) / nbead_per_DNA                                        !
   !      yCM_DNA(iDNA) = sum( y_DNA(ibead:jbead) ) / nbead_per_DNA                                        !
   !      zCM_DNA(iDNA) = sum( z_DNA(ibead:jbead) ) / nbead_per_DNA                                        !
   !                                                                                                       !
   !      delx = xCM_DNA(iDNA) - xCM_DNA_0(iDNA)                                                           !
   !      dely = yCM_DNA(iDNA) - yCM_DNA_0(iDNA)                                                           !
   !      delz = zCM_DNA(iDNA) - zCM_DNA_0(iDNA)                                                           !
   !      del  = delx*delx + dely*dely + delz*delz                                                         !
   !      aMSD_DNA_CM(icount_MSD) = aMSD_DNA_CM(icount_MSD) + del                                          !
   !    enddo                                                                                              !
   !    !$OMP END DO                                                                                       !
   !    !$OMP END PARALLEL                                                                                 !
   !  endif                                                                                                !
   !                                                                                                       !
   !else if(iter == iter_ther) then                                                                        ! else if(iter == iter_ther) then
   !                                                                                                       !
   !  aMSD_TDP = aMSD_TDP / nTDP ; aMSD_DNA_CM = aMSD_DNA_CM / nDNA                                        !
   !                                                                                                       !
   !  do ii = 1, icount_MSD                                                                                !
   !    if(ii/9*9 /=ii ) then                                                                              !
   !      jj = 10**(ii/9) * mod(ii,9)                                                                      !
   !    else                                                                                               !
   !      jj = 10**(ii/9 - 1) * 9                                                                          !
   !    endif                                                                                              !
   !    write(114,173) jj, aMSD_DNA_CM(ii), aMSD_TDP(ii)                                                   !
   !  enddo                                                                                                !
   !  close(114)                                                                                           !
   !                                                                                                       !
   !  aMSD_TDP = 0.d0 ; aMSD_DNA_CM = 0.d0 ; icount_MSD = 0                                                !
   !                                                                                                       ! writing early-time MSD
   !  x_TDP_0 = x_TDP ; y_TDP_0 = y_TDP ; z_TDP_0 = z_TDP                                                  !           &
   !                                                                                                       ! setting 0-th coords for steady-state MSD
   !  !$OMP PARALLEL DEFAULT(SHARED)                                                                       !
   !  !$OMP DO PRIVATE(ibead,jbead)                                                                        !
   !  do iDNA = 1, nDNA                                                                                    !
   !    ibead = (iDNA-1)*nbead_per_DNA                                                                     !
   !    jbead = ibead + nbead_per_DNA - 1                                                                  !
   !                                                                                                       !
   !    xCM_DNA_0(iDNA) = sum( x_DNA(ibead:jbead) ) / nbead_per_DNA                                        !
   !    yCM_DNA_0(iDNA) = sum( y_DNA(ibead:jbead) ) / nbead_per_DNA                                        !
   !    zCM_DNA_0(iDNA) = sum( z_DNA(ibead:jbead) ) / nbead_per_DNA                                        !
   !  enddo                                                                                                !
   !  !$OMP END DO                                                                                         !
   !  !$OMP END PARALLEL                                                                                   !
   !                                                                                                       !
   !else                                                                                                   ! else
   !                                                                                                       !
   !  iter_sep = iter - iter_ther                                                                          !
   !  ii = 10**floor(log10(float(iter_sep)))                                                               !
   !  if(mod(iter_sep,ii)==0) then                                                                         !
   !    icount_MSD = icount_MSD + 1                                                                        !
   !                                                                                                       !
   !    do iTDP = 1, nTDP                                                                                  !
   !      delx = x_TDP(iTDP) - x_TDP_0(iTDP)                                                               !
   !      dely = y_TDP(iTDP) - y_TDP_0(iTDP)                                                               !
   !      delz = z_TDP(iTDP) - z_TDP_0(iTDP)                                                               !
   !      aMSD_TDP(icount_MSD) = aMSD_TDP(icount_MSD) + delx*delx + dely*dely + delz*delz                  !
   !    enddo                                                                                              !
   !                                                                                                       !
   !    !$OMP PARALLEL DEFAULT(SHARED)                                                                     ! calculating steady-state MSD
   !    !$OMP DO PRIVATE(ibead,jbead,delx,dely,delz,del)  &                                                !                  
   !    !$OMP& REDUCTION(+:aMSD_DNA_CM)                                                                    !
   !    do iDNA = 1, nDNA                                                                                  !
   !      ibead = (iDNA-1)*nbead_per_DNA                                                                   !
   !      jbead = ibead + nbead_per_DNA - 1                                                                !
   !                                                                                                       !
   !      xCM_DNA(iDNA) = sum( x_DNA(ibead:jbead) ) / nbead_per_DNA                                        !
   !      yCM_DNA(iDNA) = sum( y_DNA(ibead:jbead) ) / nbead_per_DNA                                        !
   !      zCM_DNA(iDNA) = sum( z_DNA(ibead:jbead) ) / nbead_per_DNA                                        !
   !                                                                                                       !
   !      delx = xCM_DNA(iDNA) - xCM_DNA_0(iDNA)                                                           !
   !      dely = yCM_DNA(iDNA) - yCM_DNA_0(iDNA)                                                           !
   !      delz = zCM_DNA(iDNA) - zCM_DNA_0(iDNA)                                                           !
   !      del  = delx*delx + dely*dely + delz*delz                                                         !
   !      aMSD_DNA_CM(icount_MSD) = aMSD_DNA_CM(icount_MSD) + del                                          !
   !    enddo                                                                                              !
   !    !$OMP END DO                                                                                       !
   !    !$OMP END PARALLEL                                                                                 !
   !  endif                                                                                                !
   !                                                                                                       !
   !endif                                                                                                  ! endif 
   !-------------------------------------------------------------------------------------------------------!

   !----------------------------------------------------------------------------------------------------------------------------------

 enddo                                                                                             ! end of iter loop

 call CPU_TIME(t2)
 print*, 'CPU_TIME :', t2-t1

 ierr_G = vsldeletestream( stream_G )
 !ierr_U = vsldeletestream( stream_U )

! write(*,*)  'delt, vel_max, iter_max, ibead_max'       ! Correct only for serial execution
! write(*,175) delt, sqrt(vel_max), iter_max, ibead_max  !


 if(nbead_tot_DNA > 0) aMSD_DNA = aMSD_DNA / nbead_tot_DNA        !
                                                                  !
 do ii = 1, icount_MSD                                            !
   if(ii/9*9 /=ii ) then                                          !
     jj = 10**(ii/9) * mod(ii,9)                                  !
   else                                                           !
     jj = 10**(ii/9 - 1) * 9                                      !
   endif                                                          ! writing steady-state MSD
   write(115,173) jj, aMSD_DNA(ii)                                !
 enddo                                                            !
 close(115)                                                       !

 close(111) ; close(112) ; close(113)

  171 format(' ',3F25.12,I10)
  172 format(' ',10F25.12)
  173 format(' ',I10,8F25.12)
  174 format(' ',4I10)
  175 format(' ',2F25.12,2I10)
  176 format(' ',3I10,2F25.12)
  177 format(' ',2I10,5F25.12)
 
 end program main

!----------------------------------------------------------------------------------------------------------------------------------

 real*8 function force_vol_rep_DNA(r_sq)
 use cpu_mod

 implicit real*8(a-h,p-z)
 implicit integer*4(i-o)

 force_vol_rep_DNA = fac_vol_rep*exp(-alp_vol_rep*r_sq)

 end function force_vol_rep_DNA

!----------------------------------------------------------------------------------------------------------------------------------

 real*8 function energy_vol_rep_DNA(r_sq)
 use cpu_mod

 implicit real*8(a-h,p-z)
 implicit integer*4(i-o)

 energy_vol_rep_DNA = eps_vol_rep*exp(-alp_vol_rep*r_sq)

 end function energy_vol_rep_DNA

!----------------------------------------------------------------------------------------------------------------------------------

 real*8 function force_vol_rep_TDP(r_sq)
 use cpu_mod

 implicit real*8(a-h,p-z)
 implicit integer*4(i-o)

 force_vol_rep_TDP = fac_vol_rep_TDP*exp(-alp_vol_rep_TDP*r_sq)

 end function force_vol_rep_TDP

!----------------------------------------------------------------------------------------------------------------------------------

 real*8 function energy_vol_rep_TDP(r_sq)
 use cpu_mod

 implicit real*8(a-h,p-z)
 implicit integer*4(i-o)

 energy_vol_rep_TDP = eps_vol_rep_TDP*exp(-alp_vol_rep_TDP*r_sq)

 end function energy_vol_rep_TDP

!----------------------------------------------------------------------------------------------------------------------------------

 real*8 function force_vol_rep_TDP_DNA(r_sq)
 use cpu_mod

 implicit real*8(a-h,p-z)
 implicit integer*4(i-o)

 force_vol_rep_TDP_DNA = fac_vol_rep_TDP_DNA*exp(-alp_vol_rep_TDP_DNA*r_sq)

 end function force_vol_rep_TDP_DNA

!----------------------------------------------------------------------------------------------------------------------------------

 real*8 function energy_vol_rep_TDP_DNA(r_sq)
 use cpu_mod

 implicit real*8(a-h,p-z)
 implicit integer*4(i-o)

 energy_vol_rep_TDP_DNA = eps_vol_rep_TDP_DNA*exp(-alp_vol_rep_TDP_DNA*r_sq)

 end function energy_vol_rep_TDP_DNA

!----------------------------------------------------------------------------------------------------------------------------------

 real*8 function force_attr_Other(r)
 use cpu_mod

 implicit real*8(a-h,p-z)
 implicit integer*4(i-o)

 r_comp = peak_attr - r
 force_attr_Other = fac_attr * exp(-alp_attr*r_comp*r_comp) * (1 + alp_attr*r*r_comp)

 end function force_attr_Other

!----------------------------------------------------------------------------------------------------------------------------------

 real*8 function energy_attr_Other(r)
 use cpu_mod

 implicit real*8(a-h,p-z)
 implicit integer*4(i-o)

 r_comp = peak_attr - r
 energy_attr_Other = eps_attr*r*r*exp(-alp_attr*r_comp*r_comp)

 end function energy_attr_Other

!----------------------------------------------------------------------------------------------------------------------------------

 real*8 function force_attr_TDP(r)
 !real*8 function force_attr_TDP(r_sq)
 use cpu_mod

 implicit real*8(a-h,p-z)
 implicit integer*4(i-o)

 r_comp = peak_attr_TDP - r                                                                      ! pot minimum not at zero sep
 force_attr_TDP = fac_attr_TDP * exp(-alp_attr_TDP*r_comp*r_comp) * (1 + alp_attr_TDP*r*r_comp)  !

 !force_attr_TDP = fac_attr_TDP*exp(-alp_attr_TDP*r_sq)                                            ! pot minimum     at zero sep

 end function force_attr_TDP

!----------------------------------------------------------------------------------------------------------------------------------

 real*8 function energy_attr_TDP(r)
 !real*8 function energy_attr_TDP(r_sq)
 use cpu_mod

 implicit real*8(a-h,p-z)
 implicit integer*4(i-o)

 r_comp = peak_attr_TDP - r                                                                      ! pot minimum not at zero sep
 energy_attr_TDP = eps_attr_TDP*r*r*exp(-alp_attr_TDP*r_comp*r_comp)                             !

 !energy_attr_TDP = eps_attr_TDP*exp(-alp_attr_TDP*r_sq)                                           ! pot minimum     at zero sep

 end function energy_attr_TDP

!----------------------------------------------------------------------------------------------------------------------------------

 real*8 function force_attr_TDP_DNA_0(r)
 !real*8 function force_attr_TDP_DNA(r_sq)
 use cpu_mod

 implicit real*8(a-h,p-z)
 implicit integer*4(i-o)

 r_comp = peak_attr_TDP_DNA - r                                                                                 !
 force_attr_TDP_DNA_0 = fac_attr_TDP_DNA_0 * exp(-alp_attr_TDP_DNA*r_comp*r_comp) * (1 + alp_attr_TDP_DNA*r*r_comp) ! pot minimum not at zero sep

 !force_attr_TDP_DNA = fac_attr_TDP_DNA*exp(-alp_attr_TDP_DNA*r_sq)                                               ! pot minimum     at zero sep

 end function force_attr_TDP_DNA_0

!----------------------------------------------------------------------------------------------------------------------------------

 real*8 function energy_attr_TDP_DNA_0(r)
 !real*8 function energy_attr_TDP_DNA(r_sq)
 use cpu_mod

 implicit real*8(a-h,p-z)
 implicit integer*4(i-o)

 r_comp = peak_attr_TDP_DNA - r                                                                                !
 energy_attr_TDP_DNA_0 = eps_attr_TDP_DNA_0*r*r*exp(-alp_attr_TDP_DNA*r_comp*r_comp)                               ! pot minimum not at zero sep

 !energy_attr_TDP_DNA = eps_attr_TDP_DNA*exp(-alp_attr_TDP_DNA*r_sq)                                             ! pot minimum     at zero sep

 end function energy_attr_TDP_DNA_0

!----------------------------------------------------------------------------------------------------------------------------------

 real*8 function force_attr_TDP_DNA(r)
 !real*8 function force_attr_TDP_DNA(r_sq)
 use cpu_mod

 implicit real*8(a-h,p-z)
 implicit integer*4(i-o)

 r_comp = peak_attr_TDP_DNA - r                                                                                 !
 force_attr_TDP_DNA = fac_attr_TDP_DNA * exp(-alp_attr_TDP_DNA*r_comp*r_comp) * (1 + alp_attr_TDP_DNA*r*r_comp) ! pot minimum not at zero sep

 !force_attr_TDP_DNA = fac_attr_TDP_DNA*exp(-alp_attr_TDP_DNA*r_sq)                                               ! pot minimum     at zero sep

 end function force_attr_TDP_DNA

!----------------------------------------------------------------------------------------------------------------------------------

 real*8 function energy_attr_TDP_DNA(r)
 !real*8 function energy_attr_TDP_DNA(r_sq)
 use cpu_mod

 implicit real*8(a-h,p-z)
 implicit integer*4(i-o)

 r_comp = peak_attr_TDP_DNA - r                                                                                !
 energy_attr_TDP_DNA = eps_attr_TDP_DNA*r*r*exp(-alp_attr_TDP_DNA*r_comp*r_comp)                               ! pot minimum not at zero sep

 !energy_attr_TDP_DNA = eps_attr_TDP_DNA*exp(-alp_attr_TDP_DNA*r_sq)                                             ! pot minimum     at zero sep

 end function energy_attr_TDP_DNA

!----------------------------------------------------------------------------------------------------------------------------------

 real*8 function force_fene(r_sq)
 use cpu_mod

 implicit real*8(a-h,p-z)
 implicit integer*4(i-o)

 force_fene = - ak_fene / (1 - r_sq/r_fene_sq)

 end function force_fene

!----------------------------------------------------------------------------------------------------------------------------------

 real*8 function energy_fene(r_sq)
 use cpu_mod

 implicit real*8(a-h,p-z)
 implicit integer*4(i-o)

 energy_fene = eps_fene * log(1 - r_sq/r_fene_sq)

 end function energy_fene

!----------------------------------------------------------------------------------------------------------------------------------

 real*8 function force_wall_DNA(r_sep_memb)
 use cpu_mod

 implicit real*8(a-h,p-z)
 implicit integer*4(i-o)

 if(r_sep_memb > rad_coro_DNA) then
   force_wall_DNA = fac_wall_DNA1 * xi_DNA_2 * exp(akappa_DNA_sq*(rad_coro_DNA_sq-r_sep_memb*r_sep_memb)) ! -ve is absorbed within substraction
 else
   force_wall_DNA = fac_wall_DNA1 * (rad_coro_DNA/r_sep_memb + fac_wall_DNA2 * r_sep_memb/rad_coro_DNA)
 endif

 end function force_wall_DNA

!----------------------------------------------------------------------------------------------------------------------------------

 real*8 function energy_wall_DNA(r_sep_memb)
 use cpu_mod

 implicit real*8(a-h,p-z)
 implicit integer*4(i-o)

 if(r_sep_memb > rad_coro_DNA) then
   energy_wall_DNA = fac_wall_DNA_ene*zeta_DNA_fac*erfc(akappa_DNA*r_sep_memb)
 else
   a = r_sep_memb/rad_coro_DNA
   a_sq = a*a
   energy_wall_DNA = fac_wall_DNA_ene*(-log(a)-(a_sq-1)*xi_DNA_mod+zeta_DNA)
 endif

 end function energy_wall_DNA

!----------------------------------------------------------------------------------------------------------------------------------

 real*8 function force_wall_TDP(r_sep_memb)
 use cpu_mod

 implicit real*8(a-h,p-z)
 implicit integer*4(i-o)

 if(r_sep_memb > rad_coro_TDP) then
   force_wall_TDP = fac_wall_TDP1 * xi_TDP_2 * exp(akappa_TDP_sq*(rad_coro_TDP_sq-r_sep_memb*r_sep_memb)) ! -ve is absorbed within substraction
 else
   force_wall_TDP = fac_wall_TDP1 * (rad_coro_TDP/r_sep_memb + fac_wall_TDP2 * r_sep_memb/rad_coro_TDP)
 endif

 end function force_wall_TDP

!----------------------------------------------------------------------------------------------------------------------------------

 real*8 function energy_wall_TDP(r_sep_memb)
 use cpu_mod

 implicit real*8(a-h,p-z)
 implicit integer*4(i-o)

 if(r_sep_memb > rad_coro_TDP) then
   energy_wall_TDP = fac_wall_TDP_ene*zeta_TDP_fac*erfc(akappa_TDP*r_sep_memb)
 else
   a = r_sep_memb/rad_coro_TDP
   a_sq = a*a
   energy_wall_TDP = fac_wall_TDP_ene*(-log(a)-(a_sq-1)*xi_TDP_mod+zeta_TDP)
 endif

 end function energy_wall_TDP

!----------------------------------------------------------------------------------------------------------------------------------

 real*8 function force_sigmoidal_DNA(z,len_box)
 use cpu_mod

 implicit real*8(a-h,p-z)
 implicit integer*4(i-o)

 dum = exp(- sl_sigmoid * (z - len_box) )
 dum2 = 1.d0 + dum
 dum2 = dum2 * dum2
 force_sigmoidal_DNA = - h_sigmoid*sl_sigmoid*dum / dum2 

 end function force_sigmoidal_DNA

!----------------------------------------------------------------------------------------------------------------------------------

 real*8 function energy_sigmoidal_DNA(z,len_box)
 use cpu_mod

 implicit real*8(a-h,p-z)
 implicit integer*4(i-o)

 energy_sigmoidal_DNA = h_sigmoid / ( 1.d0 + exp(- sl_sigmoid * (z - len_box) ) )

 end function energy_sigmoidal_DNA

!----------------------------------------------------------------------------------------------------------------------------------

 real*8 function force_LJ(r_sq)
 use cpu_mod

 implicit real*8(a-h,p-z)
 implicit integer*4(i-o)

 r6 = r_sq*r_sq*r_sq
 r12 = r6*r6
 force_LJ = eps_LJ24 * (2*sigma_LJ12/r12 - sigma_LJ6/r6) / r_sq

 end function force_LJ

!----------------------------------------------------------------------------------------------------------------------------------

 real*8 function energy_LJ(r_sq)
 use cpu_mod

 implicit real*8(a-h,p-z)
 implicit integer*4(i-o)

 r6 = r_sq*r_sq*r_sq
 r12 = r6*r6
 energy_LJ = eps_LJ4 * (sigma_LJ12/r12 - sigma_LJ6/r6) !+ eps_LJ

 end function energy_LJ

!----------------------------------------------------------------------------------------------------------------------------------

 subroutine initial_place_DNA(nbead_tot_DNA, nbead_per_DNA, len_box, len_box_sq, x_DNA, y_DNA, nbox_DNA, nbox_dim_DNA) 
 use cpu_mod

 implicit real*8(a-h,p-z)
 implicit integer*4(i-o)
 dimension x_DNA(nbead_tot_DNA), y_DNA(nbead_tot_DNA), nbox_DNA(len_box_sq, nbox_dim_DNA) 

! epsinit = 0.1       !for making close distance among DNA

 nDNA = nbead_tot_DNA / nbead_per_DNA  ! Assumed equal lengths for all DNAs
 alen_box = dfloat(len_box)
 alen_box_half = alen_box / 2

 nbox_DNA = 0 

 do iDNA = 1, nDNA
   ibead = (iDNA - 1) * nbead_per_DNA + 1

!   x_DNA(ibead) = len_box * grnd() ; ix = ceiling(x_DNA(ibead))
!   y_DNA(ibead) = len_box * grnd() ; iy = ceiling(y_DNA(ibead)) 

!   x_DNA(ibead) = 0.5 * len_box + epsinit * len_box * (-0.5+grnd()) ;
 ix = ceiling(x_DNA(ibead))
!   y_DNA(ibead) = 0.5 * len_box + epsinit * len_box * (-0.5+grnd()) ; 
iy = ceiling(y_DNA(ibead))  !close separation among DNAs 

   icell = (iy-1)*len_box + ix
   nbox_DNA(icell, 1) = nbox_DNA(icell, 1) + 1
   nbox_DNA(icell, nbox_DNA(icell, 1) + 1) = ibead 

   write(1000000000,171) x_DNA(ibead), y_DNA(ibead)

   do jbead = 2, nbead_per_DNA
     ibead = (iDNA - 1) * nbead_per_DNA + jbead
    
     r_stretch = bead_size + (stretch_init - bead_size) * grnd()
     theta = (2*grnd() - 1) * pi 

!     x_DNA(ibead) = x_DNA(ibead-1) + r_stretch*cos(theta)
!     y_DNA(ibead) = y_DNA(ibead-1) + r_stretch*sin(theta)

     if(x_DNA(ibead) > alen_box) x_DNA(ibead) = x_DNA(ibead) - alen_box   !
     if(y_DNA(ibead) > alen_box) y_DNA(ibead) = y_DNA(ibead) - alen_box   !
     if(x_DNA(ibead) < 0.d0)     x_DNA(ibead) = x_DNA(ibead) + alen_box   ! PBC
     if(y_DNA(ibead) < 0.d0)     y_DNA(ibead) = y_DNA(ibead) + alen_box   !

     ix = ceiling(x_DNA(ibead)) ; iy = ceiling(y_DNA(ibead)) 
     icell = (iy-1)*len_box + ix
     nbox_DNA(icell, 1) = nbox_DNA(icell, 1) + 1
     nbox_DNA(icell, nbox_DNA(icell, 1) + 1) = ibead 

     write(1000000000,171) x_DNA(ibead), y_DNA(ibead)
   enddo
 enddo

 close(1000000000)

  171 format(' ',3G25.17)

 return
 end subroutine initial_place_DNA

!----------------------------------------------------------------------------------------------------------------------------------

 subroutine type_set_DNA(nbead_tot_DNA, den_TDP_type, den_OTHER_type, itype_DNA, nTDPtype, nOTHERtype)

 implicit real*8(a-h,p-z)
 implicit integer*4(i-o)
 dimension itype_DNA(nbead_tot_DNA) 

 !itype_DNA = 0                             !
 !nTDPtype = 0 ; nOTHERtype = 0             !
 !do ibead = 1, nbead_tot_DNA               !
 !  if(grnd() > 0.5d0) then                 !
 !    if(grnd() < den_TDP_type) then        !
 !      itype_DNA(ibead) = 1                !
 !      nTDPtype = nTDPtype + 1             !
 !    endif                                 ! default is 0-type
 !  else                                    !
 !    if(grnd() < den_OTHER_type) then      !
 !      itype_DNA(ibead) = 2                !
 !      nOTHERtype = nOTHERtype + 1         !
 !    endif                                 !
 !  endif                                   !
 !enddo                                     !

 itype_DNA = 1                              !
 nOTHERtype = 0                             !
 do ibead = 1, nbead_tot_DNA                !
   if(grnd() < den_OTHER_type) then         !
     itype_DNA(ibead) = 2                   ! default is 1-type
     nOTHERtype = nOTHERtype + 1            !
   endif                                    !
 enddo                                      !
 nTDPtype = nbead_tot_DNA - nOTHERtype      !

 return
 end subroutine type_set_DNA

!---------------------------------------------------------------------------------------------------------------------------------!

      subroutine sgrnd(seed)
      implicit integer(a-z)
      parameter(N     =  624)
      dimension mt(0:N-1)
      common /block/mti,mt
      save   /block/
      mt(0)= iand(seed,-1)
      do 1000 mti=1,N-1
        mt(mti) = iand(69069 * mt(mti-1),-1)
 1000 continue
      return
      end

!---------------------------------------------------------------------------------------------------------------------------------!

      double precision function grnd()
      implicit integer(a-z)
      real*8 tiny
      parameter (tiny=4.450147717014403D-308)
      parameter(N     =  624)
      parameter(N1    =  N+1)
      parameter(M     =  397)
      parameter(MATA  = -1727483681)
      parameter(UMASK = -2117483648)
      parameter(LMASK =  2147483647)
      parameter(TMASKB= -1658038656)
      parameter(TMASKC= -272236544)
      dimension mt(0:N-1)
      common /block/mti,mt
      save   /block/
      data   mti/N1/
      dimension mag01(0:1)
      data mag01/0, MATA/
      save mag01
      TSHFTU(y)=ishft(y,-11)
      TSHFTS(y)=ishft(y,7)
      TSHFTT(y)=ishft(y,15)
      TSHFTL(y)=ishft(y,-18)
      if(mti.ge.N) then
       if(mti.eq.N+1) then
         call sgrnd(4357)
       endif
        do 1000 kk=0,N-M-1
            y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
            mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
 1000   continue
        do 1100 kk=N-M,N-2
            y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
            mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
 1100   continue
        y=ior(iand(mt(N-1),UMASK),iand(mt(0),LMASK))
        mt(N-1)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
        mti = 0
      endif
      y=mt(mti)
      mti=mti+1
      y=ieor(y,TSHFTU(y))
      y=ieor(y,iand(TSHFTS(y),TMASKB))
      y=ieor(y,iand(TSHFTT(y),TMASKC))
      y=ieor(y,TSHFTL(y))
      if(y.lt.0) then
        grnd=(dble(y)+2.0d0**32)/2.0d0**32+tiny
      else
        grnd=dble(y)/2.0d0**32+tiny
      endif
      return
      end

!---------------------------------------------------------------------------------------------------------------------------------!

