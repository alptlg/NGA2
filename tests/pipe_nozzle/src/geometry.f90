!> Various definitions and tools for initializing NGA2 config
module geometry
   use ibconfig_class, only: ibconfig
   use precision,      only: WP
   implicit none
   private
   
   !> Single config
   type(ibconfig), public :: cfg
   
   !> Pipe diameter
   real(WP), public :: Dpipe
   
   public :: geometry_init
   

contains
   
   
   !> Initialization of problem geometry
   subroutine geometry_init
      use sgrid_class, only: sgrid
      use param,       only: param_read
      implicit none
      type(sgrid) :: grid
      
      
      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         integer :: i,j,k,nx,ny,nz,no
         real(WP) :: Lx,Ly,Lz,dx
         real(WP), dimension(:), allocatable :: x,y,z
         
         ! Read in grid definition
         call param_read('Pipe length',Lx)
         call param_read('Pipe diameter',Dpipe)
         call param_read('ny',ny); allocate(y(ny+1))
         call param_read('nx',nx); allocate(x(nx+1))
         call param_read('nz',nz); allocate(z(nz+1))
         
         dx=Lx/real(nx,WP)
         no=6
         if (ny.gt.1) then
            Ly=Dpipe+real(2*no,WP)*Dpipe/real(ny-2*no,WP)
         else
            Ly=dx
         end if
         if (nz.gt.1) then
            Lz=Dpipe+real(2*no,WP)*Dpipe/real(nz-2*no,WP)
         else
            Lz=dx
         end if
         
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=2,x=x,y=y,z=z,xper=.false.,yper=.true.,zper=.true.,name='pipe')
         
      end block create_grid
         
      
      ! Create a config from that grid on our entire group
      create_cfg: block
         use parallel, only: group
         integer, dimension(3) :: partition
         ! Read in partition
         call param_read('Partition',partition,short='p')
         ! Create partitioned grid
         cfg=ibconfig(grp=group,decomp=partition,grid=grid)
      end block create_cfg
      
      
      ! Create IB walls for this config (smooth CD nozzle)
      create_walls: block
         use ibconfig_class, only: sharp
         use param,         only: param_read
         integer :: i,j,k
         real(WP) :: xm, ym, zm, r_cyl
         real(WP) :: x0, L_conv, L_div, R_pipe, R_throat, R_x
         real(WP) :: t_conv, t_div
         real(WP), parameter :: pi = 3.141592653589793_WP

         ! CD nozzle parameters (read with defaults)
         R_pipe = 0.5_WP * Dpipe
         call param_read('Throat center x',x0,default=2.0_WP)
         call param_read('Convergent length',L_conv,default=0.5_WP*Dpipe)
         call param_read('Divergent length',L_div,default=0.5_WP*Dpipe)
         call param_read('Throat radius',R_throat,default=0.25_WP*Dpipe)

         ! Build signed distance: Gib = radial distance from axis - local radius R(x)
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_

                  xm = cfg%xm(i)
                  ym = cfg%ym(j)
                  zm = cfg%zm(k)
                  r_cyl = sqrt(ym*ym + zm*zm)

                  ! Local radius R(x): smooth CD profile (cosine blend, zero slope at ends)
                  if (xm <= x0 - L_conv) then
                     R_x = R_pipe
                  else if (xm < x0) then
                     t_conv = (xm - (x0 - L_conv)) / L_conv
                     R_x = R_pipe - (R_pipe - R_throat) * 0.5_WP * (1.0_WP - cos(pi * t_conv))
                  else if (xm <= x0 + L_div) then
                     t_div = (xm - x0) / L_div
                     R_x = R_throat + (R_pipe - R_throat) * 0.5_WP * (1.0_WP - cos(pi * t_div))
                  else
                     R_x = R_pipe
                  end if

                  cfg%Gib(i,j,k) = r_cyl - R_x

               end do
            end do
         end do
         ! Get normal vector
         call cfg%calculate_normal()
         ! Get VF field
         call cfg%calculate_vf(method=sharp,allow_zero_vf=.false.)

                  
      end block create_walls
      
      
   end subroutine geometry_init
   
   
end module geometry