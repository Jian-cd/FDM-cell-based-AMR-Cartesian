module mesh
    use Global_Var
    use neighbor
    use frame
    use outfile
    implicit none
    contains

subroutine creat_mesh
    implicit none
    integer i,ni,nj
    type(grid),pointer::c
    real( prec ) tx,ty

    m=(Xmax-Xmin)/(x_factor*h);  n=(Ymax-Ymin)/h;  total=m*n
    allocate(cell(total)) 
    
    print*,"Creat Background mesh... "
    do i=1,total
        ni=int((i-1)/m)+1
        nj=i-m*(ni-1)
    
        cell(i)%num(1)=ni; cell(i)%num(2)=nj
        cell(i)%lvl=0
        cell(i)%center%x=0.5d0*(2.d0*nj-1.d0)*x_factor*h+Xmin
        cell(i)%center%y=(ni-0.5d0)*h+Ymin
        cell(i)%sp=0
    
        tx=cell(i)%center%x
        ty=cell(i)%center%y 
        cell(i)%sort=gridsort(tx,ty)
        
        cell(i)%father=>NULL()
        c=>cell(i)
        cell(i)%cross=gridcross(c)
        nullify(cell(i)%son1,cell(i)%son2,cell(i)%son3,cell(i)%son4)
    end do
    
    print*,"Background mesh done"
    
    if (lvlblock>0)    then
        print*,"Start Block fine...  "
        call  Block_refine
    end if   
    
return
end subroutine creat_mesh    

subroutine Block_refine
    implicit none
    integer i,j
    type(grid),pointer::c
    real( prec ) tStart,tEnd
    
    call CPU_TIME(tStart)
    
    do i=1,total
        c=>cell(i)   
        call Block_refi(c) 
    end do  
    
    do i=1,total
        c=>cell(i)   
        call BlockDis_refi(c) 
    end do  
    
    do j=1,lvlblock-1
        do i=1,total
            c=>cell(i) 
            if(associated(c%son1))then  
                call gridmodify(c)         
            end if
        end do
    end do      
    
    call CPU_TIME(tEnd)
    write(*,'(1X,A26,F10.2)') "Block refine time:        ", tEnd-tStart
return
end subroutine Block_refine  

integer function gridsort(xx,yy)         
    implicit none
    real( prec ):: xx,yy
    
    gridsort=0
return
end function gridsort  

integer function gridcross(c) 
    implicit none
    type(grid),pointer::c
    real( prec )::cx,cy   
    real( prec )::xx1,yy1,xx2,yy2,xx3,yy3,xx4,yy4
    integer       t0,t1,t2,t3,t4,lvl
        
    if( c%num(1)< 4 .and. c%num(2)< 4) then
        gridcross=BCcorner  
    else if( c%num(1)< 4 .and. c%num(2)> m-3) then
        gridcross=BCcorner   
    else if( c%num(1)> n-3 .and. c%num(2)> m-3) then
        gridcross=BCcorner    
    else if( c%num(1)> n-3 .and. c%num(2)< 4) then
        gridcross=BCcorner             
    else if( c%num(1)< 4 ) then
        gridcross=-1       
    else if( c%num(1)> n-3 ) then
        gridcross=-3            
    else if( c%num(2)> m-3 ) then
        gridcross=-2        
    else if( c%num(2)< 4 ) then
        gridcross=-4        
    else      
        lvl=c%lvl; cx=c%center%x; cy=c%center%y
        xx1=cx-x_factor*h/(2**(lvl+1))
        yy1=cy-h/(2**(lvl+1))
        xx2=cx+x_factor*h/(2**(lvl+1))
        yy2=cy-h/(2**(lvl+1))
        xx3=cx-x_factor*h/(2**(lvl+1))
        yy3=cy+h/(2**(lvl+1))
        xx4=cx+x_factor*h/(2**(lvl+1))
        yy4=cy+h/(2**(lvl+1))  
        
        t0=gridsort(cx,cy)
        t1=gridsort(xx1,yy1)   
        t2=gridsort(xx2,yy2)
        t3=gridsort(xx3,yy3)
        t4=gridsort(xx4,yy4)
        
        if(  t1+t2+t3+t4==4 )then     
            gridcross=3       
        else if( t0==1 )then   
            gridcross=2          
        else 
            if( t1+t2+t3+t4==0 )then  
                gridcross=0 
            else                                    
                gridcross=1
            end if
        end if
    end if   
return
end function gridcross  

recursive subroutine BlockDis_refi(c)
    implicit none
    type(grid),pointer  ::c
    real( prec )        ::dis
    integer i
    
    if( c%lvl>=lvlblock ) return
    if(associated(c%son1))then
        call BlockDis_refi(c%son1)
        call BlockDis_refi(c%son2)
        call BlockDis_refi(c%son3)
        call BlockDis_refi(c%son4)
        return
    endif    
    if( c%sort==1 ) return    
    
    dis=miblock(c)
    do  i=1, lvlblock
        if( dis<=refinesize(i).and.c%sort==0.and.c%lvl<i )then     
            call Refine_cell(c)                        
            call BlockDis_refi(c%son1)
            call BlockDis_refi(c%son2)
            call BlockDis_refi(c%son3)
            call BlockDis_refi(c%son4)
            exit
        end if
    end do
return
end subroutine BlockDis_refi

recursive subroutine Block_refi(c)
    implicit none
    type(grid),pointer::c
    
    if( c%lvl>=lvlblock ) return
    if(associated(c%son1))then
        call Block_refi(c%son1)
        call Block_refi(c%son2)
        call Block_refi(c%son3)
        call Block_refi(c%son4)
        return
    endif    
    if( c%sort==1 ) return
    
    if (c%center%x<block_size(1).or.c%center%x>block_size(2).or.&
        c%center%y<block_size(3).or.c%center%y>block_size(4) ) return
       
    call Refine_cell(c) 
    call Block_refi(c%son1)
    call Block_refi(c%son2)
    call Block_refi(c%son3)
    call Block_refi(c%son4)    
return
end subroutine Block_refi

subroutine Refine_cell(c)
    implicit none
    type(grid),pointer::c
    
    allocate(c%son1,c%son2,c%son3,c%son4)  
    c%son1%father=>c; c%son1%num=c%num; c%son1%lvl=c%lvl+1; c%son1%sp=1
    c%son1%center%x=c%center%x-x_factor*h/(2**(c%lvl+2))
    c%son1%center%y=c%center%y-h/(2**(c%lvl+2))
    c%son1%sort=gridsort(c%son1%center%x,c%son1%center%y)
    c%son1%cross=gridcross(c%son1)
    
    c%son1%U=c%U
    c%son1%Ut=c%son1%U; c%son1%Un=c%son1%U
    c%son1%rot=c%rot
    c%son1%div=c%div
    
    c%son1%eru=c%eru; c%son1%erv=c%erv; c%son1%erp=c%erp
    nullify(c%son1%son1,c%son1%son2,c%son1%son3,c%son1%son4)   

    c%son2%father=>c; c%son2%num=c%num; c%son2%lvl=c%lvl+1; c%son2%sp=2
    c%son2%center%x=c%center%x+x_factor*h/(2**(c%lvl+2))
    c%son2%center%y=c%center%y-h/(2**(c%lvl+2))
    c%son2%sort=gridsort(c%son2%center%x,c%son2%center%y)
    c%son2%cross=gridcross(c%son2)
    
    c%son2%U=c%U
    c%son2%Ut=c%U; c%son2%Un=c%son2%U
    c%son2%rot=c%rot
    c%son2%div=c%div

    c%son2%eru=c%eru; c%son2%erv=c%erv; c%son2%erp=c%erp
    nullify(c%son2%son1,c%son2%son2,c%son2%son3,c%son2%son4)
      
    c%son3%father=>c; c%son3%num=c%num; c%son3%lvl=c%lvl+1; c%son3%sp=3
    c%son3%center%x=c%center%x+x_factor*h/(2**(c%lvl+2))
    c%son3%center%y=c%center%y+h/(2**(c%lvl+2))
    c%son3%sort=gridsort(c%son3%center%x,c%son3%center%y)
    c%son3%cross=gridcross(c%son3)
    
    c%son3%U=c%U
    c%son3%Ut=c%U; c%son3%Un=c%son3%U
    c%son3%rot=c%rot
    c%son3%div=c%div

    c%son3%eru=c%eru; c%son3%erv=c%erv; c%son3%erp=c%erp
    nullify(c%son3%son1,c%son3%son2,c%son3%son3,c%son3%son4)
      
    c%son4%father=>c; c%son4%num=c%num; c%son4%lvl=c%lvl+1; c%son4%sp=4
    c%son4%center%x=c%center%x-x_factor*h/(2**(c%lvl+2))
    c%son4%center%y=c%center%y+h/(2**(c%lvl+2))
    c%son4%sort=gridsort(c%son4%center%x,c%son4%center%y)
    c%son4%cross=gridcross(c%son4)
          
    c%son4%U=c%U
    c%son4%Ut=c%U; c%son4%Un=c%son4%U
    c%son4%rot=c%rot
    c%son4%div=c%div

    c%son4%eru=c%eru; c%son4%erv=c%erv; c%son4%erp=c%erp
    nullify(c%son4%son1,c%son4%son2,c%son4%son3,c%son4%son4)
return
end subroutine Refine_cell

recursive subroutine gridmodify(c)
    implicit none
    type(grid),pointer::c,cjp,cjm,cim,cip

    if(associated(c%son1))then    
        call gridmodify(c%son1)
        call gridmodify(c%son2)
        call gridmodify(c%son3)
        call gridmodify(c%son4)
    else
        cjp=>up_neighbor(c)
        cjm=>down_neighbor(c)
        cim=>left_neighbor(c)
        cip=>right_neighbor(c)
        
        if(associated(cjp))then
            if((c%lvl-cjp%lvl)>1)then
                call Refine_cell(cjp)
            end if
        end if
  
        if(associated(cjm))then
            if((c%lvl-cjm%lvl)>1)then
                call Refine_cell(cjm)
            end if
        end if
  
        if(associated(cim))then  
            if((c%lvl-cim%lvl)>1)then
                call Refine_cell(cim)
            end if
        end if
  
        if(associated(cip))then  
            if((c%lvl-cip%lvl)>1)then
                call Refine_cell(cip)
            end if
        end if
    end if

return
end subroutine gridmodify

recursive subroutine gridmodify1(c)
    implicit none
    type(grid), pointer :: c, cjp, cjm, cim, cip
    logical :: left_has_son, right_has_son, up_has_son, down_has_son
    logical :: is_leaf

    if (associated(c%son1)) then
        call gridmodify1(c%son1)
        call gridmodify1(c%son2)
        call gridmodify1(c%son3)
        call gridmodify1(c%son4)
    end if

    if (c%cross /= 0) return

    cjp => up_neighbor(c)
    cjm => down_neighbor(c)
    cim => left_neighbor(c)
    cip => right_neighbor(c)

    up_has_son    = associated(cjp) .and. associated(cjp%son1)
    down_has_son  = associated(cjm) .and. associated(cjm%son1)
    left_has_son  = associated(cim) .and. associated(cim%son1)
    right_has_son = associated(cip) .and. associated(cip%son1)

    is_leaf = .not. associated(c%son1)

    if (is_leaf) then
        if ( (left_has_son .and. right_has_son) .or. &
             (up_has_son   .and. down_has_son) ) then
            call Refine_cell(c)
        end if
    else
        if ( .not. associated(c%son1%son1) .and. &
             .not. associated(c%son2%son1) .and. &
             .not. associated(c%son3%son1) .and. &
             .not. associated(c%son4%son1) ) then

            if ( (.not. left_has_son  .and. .not. right_has_son) .or. &
                 (.not. up_has_son    .and. .not. down_has_son) ) then

                call Coarse_cell(c)
                
                return
            end if
        end if
    end if

end subroutine gridmodify1

recursive subroutine gridmodify4(c)
    implicit none
    type(grid), pointer :: c, nbor, nnbor
    integer :: dir

    if (associated(c%son1)) then
        call gridmodify4(c%son1)
        call gridmodify4(c%son2)
        call gridmodify4(c%son3)
        call gridmodify4(c%son4)
        return
    end if

    do dir = 1, 4
        select case (dir)
        case (1); nbor => right_neighbor(c)
        case (2); nbor => left_neighbor(c)
        case (3); nbor => up_neighbor(c)
        case (4); nbor => down_neighbor(c)
        end select

        if (.not. associated(nbor)) cycle

        if (nbor%lvl <= c%lvl - 2) then
            call Refine_cell(nbor)
        end if

        select case (dir)
        case (1); nnbor => right_neighbor(nbor)
        case (2); nnbor => left_neighbor(nbor)
        case (3); nnbor => up_neighbor(nbor)
        case (4); nnbor => down_neighbor(nbor)
        end select

        if (associated(nnbor)) then
            if (nnbor%lvl <= c%lvl - 2) then
                call Refine_cell(nnbor)
            end if
        end if
    end do
end subroutine gridmodify4

subroutine Coarse_cell(c)
    implicit none
    type(grid),pointer::c,c1,c2,c3,c4
    real( prec )::var(Nvar)
    
    call interpA(c,var)
    c%U(1:7)=var(1:7); c%Un=c%U

    deallocate(c%son1,c%son2,c%son3,c%son4)
    nullify(c%son1,c%son2,c%son3,c%son4)
return
end subroutine Coarse_cell

function miwall(c) 
    implicit none
    type(grid),pointer::c
    real( prec )::miwall,xo,yo,x1,x2,y1,y2,a
    real( prec )::cos1,cos2
    integer i,j,dmax
    
    miwall=0.d0
    
return
end function miwall


function miblock(c) 
    implicit none
    type(grid),pointer::c
    real( prec )::miblock,xo,yo

    xo=c%center%x
    yo=c%center%y

    if( xo>block_size(2) )then
        if( yo>block_size(4) )then
            miblock=sqrt( (xo-block_size(2))*(xo-block_size(2))+&
                          (yo-block_size(4))*(yo-block_size(4)) )
        else if( yo<block_size(3) )then
            miblock=sqrt( (xo-block_size(2))*(xo-block_size(2))+&
                          (yo-block_size(3))*(yo-block_size(3)) )  
        else
            miblock=xo-block_size(2)
        end if
    else if( xo<block_size(1) )then
        if( yo>block_size(4) )then
            miblock=sqrt( (xo-block_size(1))*(xo-block_size(1))+&
                          (yo-block_size(4))*(yo-block_size(4)) )
        else if( yo<block_size(3) )then
            miblock=sqrt( (xo-block_size(1))*(xo-block_size(1))+&
                          (yo-block_size(3))*(yo-block_size(3)) )  
        else
            miblock=block_size(1)-xo
        end if 
    else
        if( yo>block_size(4) )then
            miblock=yo-block_size(4)
        else if( yo<block_size(3) )then
            miblock=block_size(3)-yo 
        else
            miblock=0.0
        end if         
    end if  
return
end function miblock

end module mesh
