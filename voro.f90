module voro

  use VorCells_Polygons, only: Polygon,VoronoiCell,voronoi_cell,VorCellInterConvPol, &
       VorCellInterConvPol_Collinear,voronoi_cell_destroy,polygon_destroy

  implicit none
  
  type :: pdata_type
     integer :: npols
     real(kind=8) :: c(2),rmax,volA,xmin,xmax,ymin,ymax,drawscale
     integer, allocatable :: nvpols(:),poledges(:,:)
     real(kind=8), allocatable :: xpol(:,:),ypol(:,:)
  end type pdata_type
  
  public :: voronoi,drawvor

contains

  ! ******************************************************************
  ! ******************************************************************

  subroutine voronoi(nsites,sites,An,Ax,Ay,Aflag,nvmax,sstart,nv,vx,vy,vflag,istop)
	
    implicit none
  
    ! SCALAR ARGUMENTS
    integer, intent(in) :: An,nsites,nvmax
    integer, intent(out) :: istop,nv
    
    ! ARRAY ARGUMENTS
    integer, intent(in) :: Aflag(An+1)
    real(kind=8), intent(in) :: Ax(An+1),Ay(An+1),sites(2,nsites)
    integer, intent(out) :: sstart(nsites+1),vflag(nvmax)
    real(kind=8), intent(out) :: vx(nvmax),vy(nvmax)
    
    ! LOCAL SCALARS
    integer :: allocerr,i,ierror,ind,j,k,status,ntriag
    type(VoronoiCell) :: Vi
    type(Polygon) :: Wi,P
    
    ! LOCAL ARRAYS
    integer :: indices(nsites),sitetv(3*2*nsites),start(nsites+1),til(3,2*nsites),tnbr(3,2*nsites)
    real(kind=8) :: vorxy(2,2*nsites)
    
   istop = 0

    if ( nsites .lt. 3 ) then
      write(*,*) 'In voronoi: nsites must be at least 3!'
      stop
    end if

    ! Polygon representing the domain

    P%n = An + 1
    P%deg = .false.

    allocate(P%v(2,P%n),P%e(P%n-1),stat=allocerr)
    if ( allocerr .ne. 0 ) then
       write(*,*) 'In voronoi: Allocation error.'
       stop
    end if
    
    P%v(1,1:P%n) = Ax(1:An+1)
    P%v(2,1:P%n) = Ay(1:An+1)
    P%e(1:P%n-1) = Aflag(1:An+1)
          
    ! Compute Delaunay triangulation

    indices(1:nsites) = (/ (i,i=1,nsites) /)

    call dtris2(nsites,sites,indices,ntriag,til,tnbr,ierror)

    if ( ierror .ne. 0 ) then
      istop = -1
      return
    end if
    
    if ( ierror .ne. 0 .and. ierror .ne. 225 ) then
       write(*,*) 'In voronoi: dtris2 returned an error different from 0 and 225: ',ierror
       istop = -1
    end if
  
    call voronoi_cell_patch(nsites,ntriag,til,start,sitetv)
    
    ! Compute Voronoi diagram
          
    if ( ntriag .gt. 2 * nsites ) then
       write(*,*) 'In voronoi: There is something wrong. According to the dtris2, '
       write(*,*) 'documentation the number of triangles in the Delaunay trangulation '
       write(*,*) 'cannot be larger than twice the number of sites!'
       stop
    end if
  
    do k = 1,ntriag
       call triangle_circumcenter_2d(sites(1:2,til(1:3,k)),vorxy(1:2,k))
       if ( isnan(vorxy(1,k)) .or. isnan(vorxy(1,k)) ) then
          write(*,*) 'In voronoi: triangle_circumcenter_2d returned NaN.'
          write(*,*) 'Triangle vertices: ',sites(1:2,til(1:3,k))
          istop = -2
       end if
    end do

    sstart(1) = 1

    do i = 1,nsites
       call voronoi_cell(i,nsites,sites,start,sitetv,ntriag,til,tnbr,vorxy,Vi,status)
       if ( status .ne. 0 ) then
          write(*,*) 'In voronoi: Error in voronoi_cell, status = ',status
          stop
       end if

       ! Compute the intersection of the ith Voronoi cell and the domain

       if ( nsites .ge. 3 .and. ierror .eq. 0 ) then
          call VorCellInterConvPol(sites(1:2,i),P,Vi,Wi)
       else ! if ( ( nsites .ge. 3 .and. ierror .eq. 225 ) .or. nsites .eq. 2 ) then
          call VorCellInterConvPol_Collinear(i,nsites,sites,P,Wi)
       end if

       if ( Wi%deg ) then
          sstart(i+1) = sstart(i)

       else
          sstart(i+1) = sstart(i) + ( Wi%n - 1 )
          
          if ( sstart(i+1) .gt. nvmax + 1 ) then
             write(*,*) 'Increase nvmax.'
             stop
          end if
    
          do j = 1,Wi%n - 1
             ind = sstart(i) + j - 1
             vx(ind) = Wi%v(1,j)
             vy(ind) = Wi%v(2,j)
             vflag(ind) = -Wi%e(j)
          end do
       end if
       
       call polygon_destroy(Wi)

       call voronoi_cell_destroy(Vi)
    end do

    nv = sstart(nsites+1) - 1
    
    deallocate(P%v,P%e,stat=allocerr)
    if ( allocerr .ne. 0 ) then
       write(*,*) 'In voronoi: Deallocation error.'
       stop
    end if

  end subroutine voronoi

  ! ******************************************************************
  ! ******************************************************************

  subroutine drawvor(nsites,sites,colors,An,Ax,Ay,sstart,nv,vx,vy)
	
    implicit none
  
    ! SCALAR ARGUMENTS
    integer, intent(in) :: An,nsites,nv
    
    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: Ax(An+1),Ay(An+1),sites(2,nsites),vx(nv),vy(nv)
    integer, intent(in) :: sstart(nsites+1),colors(nsites)
    
    ! LOCAL SCALARS
    integer :: i,j,jnext,k,color,maxdiffi
    real(kind=8) :: centerxi,centeryi,colormin,colormax,maxdiff,diff
    
    open(unit=10,file='voronoi.mp')                           

    write(10,10)
    
    maxdiff = colors(2) - colors(1)
    maxdiffi = 1
    do i = 2,nsites - 1
      diff = colors(i+1) - colors(i)
      if ( diff .gt. maxdiff ) then
         maxdiff = diff
         maxdiffi = i
      end if
    end do

    do i = 1,nsites
      centerxi = 0.0d0
      centeryi = 0.0d0
      do j = sstart(i),sstart(i+1) - 1
         centerxi = centerxi + vx(j)
         centeryi = centeryi + vy(j)
         if ( j .eq. sstart(i) ) then 
            write(10,21) vx(j),vy(j)
         else if ( j .lt. sstart(i+1) - 1 ) then 
            write(10,22) vx(j),vy(j)
         else
            if ( maxdiff .le. 10.0d0 ) then
               colormin = minval( colors(1:nsites) )
               colormax = maxval( colors(1:nsites) )
               color = max( 0,  min( ceiling( 100 + ( 255 - 100 ) * ( colors(i) - colormin ) / ( colormax - colormin ) ) , 255 ) )
               write(10,23) vx(j),vy(j),'(  0/255,',color,'/255,230/255)'
            else
               if ( i .le. maxdiffi ) then
                  colormin = minval( colors(1:maxdiffi) )
                  colormax = maxval( colors(1:maxdiffi) )
                  color = max( 0, min( ceiling( 100 * ( colors(i) - colormin ) / ( colormax - colormin ) ) , 100 ) )
                  write(10,24) vx(j),vy(j),'(',color,'/255,240/255,192/255)'
               else
                  colormin = minval( colors(maxdiffi+1:nsites) )
                  colormax = maxval( colors(maxdiffi+1:nsites) )
                  color = max( 0, min( ceiling( 100 * ( colors(i) - colormin ) / ( colormax - colormin ) ) , 100 ) )
                  write(10,23) vx(j),vy(j),'(159/255,',color,'/255,175/255)'
               end if
            end if
         end if
      end do
      centerxi = centerxi / ( sstart(i+1) - sstart(i) )
      centeryi = centeryi / ( sstart(i+1) - sstart(i) )
      ! write(10,40) centerxi,centeryi
      write(10,41) colors(i),centerxi,centeryi
    end do

   do i = 1,nsites
      ! Draw the cell (intersected with the domain)
      do j = sstart(i),sstart(i+1) - 1
         if ( j .lt. sstart(i+1) - 1 ) then 
            jnext = j + 1
         else
            jnext = sstart(i)
         end if
         write(10,20) vx(j),vy(j),vx(jnext),vy(jnext),'black'
      end do

      ! Draw the ith site
      ! write(10,30) sites(1:2,i)
      ! write(10,31) i,sites(1:2,i)
   end do
      
    ! Draw the domain
    do k = 1,An - 1
       write(10,20) Ax(k),Ay(k),Ax(k+1),Ay(k+1),'black'
    end do

    ! Close figure
    
    write(10,50)
    
    close(10)

    ! NON-EXECUTABLE STATEMENTS
10  format('prologues := 3;',/, &
         'outputtemplate := "%j.mps";',/, &
         'input mpcolornames;',/, &
         'beginfig( 1 );',/, &
         'color skyblue, dollarbill, desire, jonquil, royalblue, teagreen, ', &
         'orange, plum, bluea, blueb, bluec, blued, bluee, corn, deepsaffron, ferrarired;',/, &
         'skyblue      = (135/255,206/255,235/255);',/, &
         'dollarbill   = (133/255,187/255,101/255);',/, &
         'teagreen     = (208/255,240/255,192/255);',/, &
         'desire       = (234/255, 60/255, 83/255);',/, &
         'jonquil      = (244/255,202/255, 22/255);',/, &
         'royalblue    = ( 65/255,105/255,225/255);',/, &
         'orange       = (255/255,165/255,  0/255);',/, &
         'plum         = (221/255,160/255,221/255);',/, &
         'bluea        = (  0/255,230/255,230/255);',/, &
         'blueb        = (  0/255,191/255,230/255);',/, &
         'bluec        = (  0/255,153/255,230/255);',/, &
         'blued        = (  0/255,115/255,230/255);',/, &
         'bluee        = (  0/255, 77/255,230/255);',/, &
         'corn         = (255/255,230/255,102/255);',/, &
         'deepsaffron  = (255/255,153/255, 51/255);',/, &
         'ferrarired   = (255/255, 42/255,  0/255);',/, &
         'u = 4.4cm;')
20  format('draw ( ',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u) withcolor',a12,';')
21  format('fill ( ( ',f20.10,'u,',f20.10,'u)--')
22  format('       ( ',f20.10,'u,',f20.10,'u)--')
23  format('       ( ',f20.10,'u,',f20.10,'u)--cycle ) withcolor',a9,i3,a13,';')
24  format('       ( ',f20.10,'u,',f20.10,'u)--cycle ) withcolor',a1,i3,a21,';')
30  format('drawdot(',f20.10,'u,',f20.10,'u);')
31  format('dotlabel.bot(btex ',i2,' etex, (',f20.10,'u,',f20.10,'u));')
40  format('drawdot(',f20.10,'u,',f20.10,'u);')
41  format('label(btex ',i3,' etex, (',f20.10,'u,',f20.10,'u));')
50  format('endfig;',/,'end;')

  end subroutine drawvor
  
end module voro
