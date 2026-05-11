Module solution
    use nndcenter
    use outfile
    use mesh
    implicit none
    procedure(Time_advance_RK), pointer     :: k_TimeMethod
    contains
  !======================================================================  
subroutine initial
    implicit none
    integer i
    type(grid),pointer::c    

    Neq=4
    
    do i=1,total            
        c=>cell(i) 
        call vortex(c)
    end do
        
    do i=1,total
        c=>cell(i)
        call initial_boundary(c)
    end do        
    call update
    
    print*, "Init Done!"
return
end subroutine initial

recursive subroutine vortex(c)
    implicit none
    type(grid),pointer::c
    real( prec )      ::xc,yc,rr,dT
    real( prec ),parameter:: cc=5.d0, Rc2=1.d0
    
    if(associated(c%son1))then
        call vortex(c%son1)
        call vortex(c%son2)
        call vortex(c%son3)
        call vortex(c%son4)
    else if(c%sort==0)then
        c%U=0.d0
        
        xc=0.5*(cell(1)%center%x+cell(m)%center%x)
        yc=0.5*(cell(m)%center%y+cell(total)%center%y)
        rr=(c%center%x-5.d0)**2+(c%center%y-5.d0)**2
        dT = -(gama-1.d0)*(cc*cc)*exp(1.d0-rr)/(8.d0*gama*Pi*Pi)
        
        c%U(1)=1.d0+(1.d0+dT)**(1.d0/(gama-1.d0))-1.d0
        c%U(2)=1.d0-(c%center%y-5.d0)*cc/(2.d0*Pi)*exp(0.5d0*(1.d0-rr))
        c%U(3)=1.d0+(c%center%x-5.d0)*cc/(2.d0*Pi)*exp(0.5d0*(1.d0-rr))
        c%U(5)=1.d0+(1.d0+dT)**(gama/(gama-1.d0))-1.d0
        
        c%U(4)=(c%U(5)*gama*Ma*Ma)/c%U(1) 
        
        c%Ut=c%U;   c%Un=c%U
        
        c%rot=0.d0;  c%div=0.d0 
        c%eru=0.d0; c%erv=0.d0; c%erp=0.d0
    else
        c%U=0.d0         
        c%Ut=c%U;   c%Un=c%U
        
        c%rot=0.d0;  c%div=0.d0        
        c%eru=0.d0; c%erv=0.d0; c%erp=0.d0
    end if
return
end subroutine vortex

recursive subroutine initial_boundary(c)
    implicit none    
    type(grid),pointer::c

    if(associated(c%son1))then
        call initial_boundary(c%son1)
        call initial_boundary(c%son2)
        call initial_boundary(c%son3)
        call initial_boundary(c%son4)       
    else if(c%sort==0)then     
        if( c%cross==-4 ) then 
            select case (kind_xmbc)
            case (BCcyc)    
                call BC_cyc(c,-4)
            end select
        
        else if( c%cross==-2 )   then 
            select case (kind_xpbc)
            case (BCcyc)    
                call BC_cyc(c,-2)
            end select       
            
        else if( c%cross==-3 )then  
            select case (kind_ypbc)
            case (BCcyc)    
                call BC_cyc(c,-3)
            end select

        else if( c%cross==-1 )then    
            select case (kind_ymbc)
            case (BCcyc)    
                call BC_cyc(c,-1) 
            end select 
        end if    
    end if
return
end subroutine initial_boundary 

subroutine BC_cyc(c,a)  
    implicit none 
    type(grid),pointer :: c, cn, cf
    real( prec ) :: cnvar(Nvar), xx, yy
    integer a, i, j, clvl
    integer, allocatable :: sp(:)
    
    if(a==-4)then
        i=c%num(1)*m+c%num(2)-6
        cn=>cell(i)
        xx = c%center%x + (Xmax - Xmin)
        yy = c%center%y
    else if(a==-2)then
        i=c%num(1)*m+c%num(2)-2*m+6
        cn=>cell(i)
        xx = c%center%x - (Xmax - Xmin)
        yy = c%center%y
    else if(a==-3)then
        i=c%num(1)*m+c%num(2)-(n-5)*m
        cn=>cell(i)
        xx = c%center%x
        yy = c%center%y  - (Ymax - Ymin)
    else if(a==-1)then
        i=c%num(1)*m+c%num(2)-(7-n)*m
        cn=>cell(i)
        xx = c%center%x
        yy = c%center%y  + (Ymax - Ymin)
    end if  
    
    if( c%lvl/=0 )then 
        clvl=c%lvl
        cf=>c
        
        allocate(sp(clvl))
        do j=1,clvl
            sp(j)=cf%sp
            cf=>cf%father
        end do
        
        do j=1,clvl
            if( .not.associated(cn%son1) ) exit
            if(sp(j)==1)then           
                cn=>cn%son1      
            else if(sp(j)==2)then   
                cn=>cn%son2 
            else if(sp(j)==3)then   
                cn=>cn%son3 
            else if(sp(j)==4)then   
                cn=>cn%son4       
            end if              
        end do
        deallocate(sp)
    end if
    
    call interpA(cn,cnvar)
    c%U(1:5)=cnvar(1:5)                                           
    
    c%U(6)=cnvar(6); c%U(7)=cnvar(7)
return
end subroutine BC_cyc 

!======================================================================
subroutine solute
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    implicit none
    type(grid),pointer::c
    integer i
    
    call IDgrid
    
    write(*,*) "Start the Calculation ..."
    print*, " --nstep-- ", " --nResgrid-- ", " --tstep-- ", " --t-- ", " --The R.M.S Residuals-- "
    
    do while (tt < t_end)
        it = it + 1
        
        l_max=0
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(c,i) SCHEDULE(DYNAMIC)     
        do i=1,nResgrid
            c=>cell_NS(i)%n
            lambda_x = max(lambda_x, ABS(c%U(2))+ SQRT(gama*c%U(5)/c%U(1)))
            lambda_y = max(lambda_y, ABS(c%U(3))+ SQRT(gama*c%U(5)/c%U(1)))    
            l_max = max(l_max, c%lvl)
        end do 
        !$OMP END PARALLEL DO
    
        if(iflag_localdt) then 
            tstep = cfl/(lambda_x/(x_factor*h/2**l_max)+lambda_y/(h/2**l_max))
        end if
        
        if( tt + tstep > t_end ) then
            tstep = t_end - tt
        end if
    
        call Time_advance_RK
        
        tt=tt+tstep
        
        call comput_error        
        
        do i=1,Neq 
            Res_rms(i)=sqrt(Res_rms(i)/(1.d0*nResgrid))     
            if (ieee_is_nan(Res_rms(i))) then
                write(*,*) 'ERROR: NaN after sqrt in equation ', i, ' at it =', it
                stop 'NaN detected in Res_rms after sqrt'
            end if
        enddo

        if(mod(it,Pstep_Res)==0 .or. tt==t_end)then       
            call output_Res   
        end if

        if (iflag_AMR)then
            call AMR
        
            !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(c,i) SCHEDULE(DYNAMIC)     
            do i=1,nBCgrid
                c=>cell_BC(i)%n
                call body_boundary_RK(c,1)
            end do 
            !$OMP END PARALLEL DO 
        end if
        
        if( tt==t_end )then
            call output_error
        endif          

    end do 
return
end subroutine solute


subroutine comput_error
    implicit none
    type(grid),pointer::c
    integer:: i
    real( prec )::uu,vv,pp,rhorho
    
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(c,i)  SCHEDULE(DYNAMIC)    
        do i=1,nResgrid 
            c=>cell_NS(i)%n
            
            rhorho = 1+0.2*dsin(Pi*(c%center%x+c%center%y-tt))
            
            uu=1.d0
            vv=1.d0
            pp=1.d0
        
            c%errho=c%U(1)-rhorho
            c%eru=(c%U(2)-uu)
            c%erv=(c%U(3)-vv)
            c%erp=(c%U(5)-pp)
        end do
    !$OMP END PARALLEL DO
return
end  subroutine comput_error

!======================================================================
subroutine Time_advance_RK
    implicit none
    type(grid),pointer::c
    integer:: i,j,KRK
    
    call getRKn
    j=3        
    
    do KRK=1,j 
        Res_max(:)=0.d0
        Res_rms(:)=0.d0
        
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(c,i) REDUCTION(MAX: Res_max) REDUCTION(+: Res_rms) SCHEDULE(DYNAMIC)            
        do i=1,nResgrid 
            c=>cell_NS(i)%n
            call solution_cell(c) 
        end do
        !$OMP END PARALLEL DO 

        call updateRK(KRK)
        
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(c,i) SCHEDULE(DYNAMIC)     
        do i=1,nBCgrid
            c=>cell_BC(i)%n
            call body_boundary_RK(c,j)
        end do 
        !$OMP END PARALLEL DO 
             
    end do
    
return
end  subroutine Time_advance_RK


!======================================================================
subroutine solution_cell(c)
    implicit none
    type(grid),pointer::c
    real( prec )::Uip(Nvar),Uim(Nvar),Ujp(Nvar),Ujm(Nvar)
     
    call Tip(c,Uip);  call Tim(c,Uim);   call Tjp(c,Ujp);   call Tjm(c,Ujm)
            
    call LXL(c,Uip,Uim,Ujp,Ujm)
    
return
end subroutine solution_cell

subroutine body_boundary_RK(c,j)
    implicit none    
    integer::j
    type(grid),pointer::c
    real( prec )::t
    
    if (j==1) then
        t=tt
    else if(j==2) then
        t=tt+tstep
    else if(j==3) then
        t=tt+0.5d0*tstep
    endif

    if(c%sort==0)then 
        if( c%cross==-4 ) then 
            select case (kind_xmbc)
            case (BCcyc) 
                call BC_cyc(c,-4)
            end select
        
        else if( c%cross==-2 )   then 
            select case (kind_xpbc)
            case (BCcyc)    
                call BC_cyc(c,-2)  
            end select     
                
        else if( c%cross==-3 )then  
            select case (kind_ypbc)
            case (BCcyc)    
                call BC_cyc(c,-3) 
            end select

        else if( c%cross==-1 )then  
            select case (kind_ymbc)
            case (BCcyc)    
                call BC_cyc(c,-1)
            end select
            
        end if    
    end if
return
end subroutine body_boundary_RK   
 
subroutine update
    implicit none
    type(grid),pointer::c
    integer i
    
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(c,i) SCHEDULE(DYNAMIC)    
    do i=1,nResgrid 
        c=>cell_NS(i)%n
        c%U(1:5)=c%Ut(1:5)
    end do
!$OMP END PARALLEL DO
return
end subroutine update

subroutine getRKn
    implicit none
    type(grid),pointer::c
    integer i
    
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(c,i) SCHEDULE(DYNAMIC)      
    do i=1,nResgrid 
        c=>cell_NS(i)%n
        c%Un(1:7)=c%U(1:7)
    end do    
!$OMP END PARALLEL DO    
return
end subroutine getRKn

subroutine updateRK(KRK)
    implicit none
    type(grid),pointer::c
    real( prec ):: Ralfa(3), Rbeta(3)
    integer i,j,KRK
    
    Ralfa(1)=0.d0 ;  Ralfa(2)=3.d0/4.d0 ; Ralfa(3)=1.d0/3.d0
    Rbeta(1)=1.d0 ;  Rbeta(2)=1.d0/4.d0 ; Rbeta(3)=2.d0/3.d0
    
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(c,i) SCHEDULE(DYNAMIC)      
    do i=1,nResgrid 
        c=>cell_NS(i)%n
        do j=1, Neq+1
            c%U(j)=Ralfa(KRK)*c%Un(j)+Rbeta(KRK)*c%Ut(j)
        end do
    end do
!$OMP END PARALLEL DO    
    return
end subroutine updateRK


!======================================================================
subroutine AMR
    implicit none    
    logical, save :: adapt_done = .false.
    integer i
    
    if( it>=Nstep_Adapt .and. mod(it,Pstep_Adapt)==1 )then  
        call rot_div
        call solute_ref
        call solute_gridmodify
        call idneighbor 
        call IDgrid
    endif
return
end subroutine AMR


subroutine rot_div
    implicit none
    type(grid),pointer::c
    integer i
    
        
    Srot=0
    Sdiv=0
    Ne=0  
    do i=1,total
        c=>cell(i)  
        if( c%cross==0 .or. c%cross==1 .or. c%cross==2 ) then
            call solute_div_rot(c)
        endif  
    end do
return
end subroutine rot_div

recursive subroutine solute_div_rot(c)
implicit none
type(grid),pointer::c
real( prec )::Uip(Nvar),Uim(Nvar),Ujp(Nvar),Ujm(Nvar)
real( prec )::Uipp(Nvar),Uimm(Nvar),Ujpp(Nvar),Ujmm(Nvar)
real( prec )::dx,dy

    if(.not.associated(c%son1).and.( c%cross==0 .or. c%cross==1 ) )then
        dy=h/2**c%lvl
        dx=x_factor*dy    
        call Tip(c,Uip);  call Tim(c,Uim); call Tjp(c,Ujp); call Tjm(c,Ujm)                                          
        call Urot(c,dx,dy,Uip,Uim,Ujp,Ujm) 
        call Udiv(c,dx,dy,Uip,Uim,Ujp,Ujm)
        
        Srot=Srot+c%rot**2
        Sdiv=Sdiv+c%div**2
        Ne=Ne+1
    elseif(associated(c%son1))then                                     
        call solute_div_rot(c%son1) 
        call solute_div_rot(c%son2)
        call solute_div_rot(c%son3)
        call solute_div_rot(c%son4)
    endif 
return
end subroutine solute_div_rot

subroutine solute_ref
    implicit none
    type(grid),pointer::c
    integer i
    real( prec )::adiv,arot
 
    arot=sqrt(Srot/Ne) 
    adiv=sqrt(Sdiv/Ne)
 
    do i=1,total
        c=>cell(i)       
        if(  c%cross==0 .or. c%cross==1  )then 
            call solute_ref_sub(c,arot,adiv)
        endif

        if( associated(c%son1) )then
            call solute_coarse_sub(c,arot,adiv)
        endif
    end do  
return
end subroutine solute_ref

recursive subroutine solute_coarse_sub(c,arot,adiv)
    implicit none
    type(grid),pointer::c   
    real( prec )::adiv,arot

    if(associated(c%son1).and.c%lvl<lvlAMR.and.&
    &.not.associated(c%son1%son1).and..not.associated(c%son2%son1).and.&
    &.not.associated(c%son3%son1).and..not.associated(c%son4%son1))then      
         if( ((c%son1%rot+c%son2%rot+c%son3%rot+c%son4%rot)/4<Coarse_goal_rot*arot .and.&
        & (c%son1%div+c%son2%div+c%son3%div+c%son4%div)/4<Coarse_goal_div*adiv    ))then
                call Coarse_cell(c) 
        end if    
     else if( associated(c%son1) )then
        if( associated(c%son1%son1) )then
            call solute_coarse_sub(c%son1,arot,adiv) 
        endif
        if( associated(c%son2%son1) )then
            call solute_coarse_sub(c%son2,arot,adiv)
        endif
        if( associated(c%son3%son1) )then  
            call solute_coarse_sub(c%son3,arot,adiv)
        endif
        if( associated(c%son4%son1) )then
            call solute_coarse_sub(c%son4,arot,adiv)
        end if
    endif 
return
end subroutine solute_coarse_sub

recursive subroutine solute_ref_sub(c,arot,adiv)
implicit none
type(grid),pointer::c
real( prec )::adiv,arot

    if( .not.associated(c%son1).and.c%lvl<lvlAMR )then 
        if( (c%rot>Refine_goal_rot*arot .or. c%div>Refine_goal_div*adiv ) ) then 
            call Refine_cell(c) 
        end if    
    elseif( associated(c%son1) .and. c%lvl<lvlAMR )then
        call solute_ref_sub(c%son1,arot,adiv) 
        call solute_ref_sub(c%son2,arot,adiv)
        call solute_ref_sub(c%son3,arot,adiv)
        call solute_ref_sub(c%son4,arot,adiv)
    endif 
return
end subroutine solute_ref_sub

subroutine solute_gridmodify
    implicit none
    type(grid),pointer::c
    integer i,j,k
    
    do k=1,lvlAMR
        do i=1,total
            c=>cell(i) 
            call gridmodify1(c)
        end do
        call idneighbor
    
        do i=1,total
            c=>cell(i) 
            call gridmodify4(c) 
        end do       
        call idneighbor
    end do
return
end subroutine solute_gridmodify

subroutine Urot(c,dx,dy,Uip,Uim,Ujp,Ujm)
    implicit none
    type(grid),pointer::c
    real( prec )::Uip(Nvar),Uim(Nvar),Ujp(Nvar),Ujm(Nvar)
    real( prec )::vx,uy,dx,dy

    vx=dux(Uip(3),Uim(3),dx); uy=dux(Ujp(2),Ujm(2),dy) 
    c%rot=abs(vx-uy)*dy**1.5   
return
end subroutine Urot

subroutine Udiv(c,dx,dy,Uip,Uim,Ujp,Ujm)
    implicit none  
    type(grid),pointer::c
    real( prec )::Uip(Nvar),Uim(Nvar),Ujp(Nvar),Ujm(Nvar)
    real( prec )::ux,vy,dx,dy
    real( prec )::qx2, qy2
    real(prec), parameter :: eps = 1.0e-15_prec
    
    ux=dux(Uip(1),Uim(1),dx); vy=dux(Ujp(1),Ujm(1),dy) 
    c%div=(abs(ux)+abs(vy))*dy**1.5
return
end subroutine Udiv


function dux(uip,uim,dx)
    implicit none
    real( prec )::dux,uip,uim,dx

    dux=(uip-uim)/(2.d0*dx)
return
end function dux


!========================================================
subroutine IDgrid
    implicit none
    type(grid),pointer::c
    integer i
    
    if(associated(cell_NS)) deallocate( cell_NS,cell_BC )  
    nResgrid=0; nBCgrid=0
    
    do i=1,total
        c=>cell(i)
        call IDNSgrid_num(c)
        call IDBCgrid_num(c)
    end do
    
    allocate(cell_NS(nResgrid),cell_BC(nBCgrid))
    nResgrid=0; nBCgrid=0
    do i=1,total
        c=>cell(i)
        call IDNSgrid(c)
        call IDBCgrid(c)
    end do    
return
end subroutine IDgrid


recursive subroutine IDNSgrid_num(c)
    implicit none
    type(grid),pointer::c
    
    if(associated(c%son1))then
        call IDNSgrid_num(c%son1)
        call IDNSgrid_num(c%son2)
        call IDNSgrid_num(c%son4)
        call IDNSgrid_num(c%son3)
    else if(c%sort==0)then
        if (c%cross<0)  return
        nResgrid=nResgrid+1
    end if    
return
end subroutine IDNSgrid_num

recursive subroutine IDBCgrid_num(c)
    implicit none    
    type(grid),pointer::c

    if(associated(c%son1))then
        call IDBCgrid_num(c%son1)
        call IDBCgrid_num(c%son2)
        call IDBCgrid_num(c%son3)
        call IDBCgrid_num(c%son4)       
    else if(c%sort==1)then
        if (c%cross==BCcorner)  return
        nBCgrid=nBCgrid+1
    else if(c%sort==0)then
        if (c%cross==BCcorner)  return
        if( c%cross<0 ) then
            nBCgrid=nBCgrid+1
        end if    
    end if
return
end subroutine IDBCgrid_num  

recursive subroutine IDNSgrid(c)
    implicit none
    type(grid),pointer::c
    
    if(associated(c%son1))then
        call IDNSgrid(c%son1)
        call IDNSgrid(c%son2)
        call IDNSgrid(c%son4)
        call IDNSgrid(c%son3)
    else if(c%sort==0)then   
        if (c%cross<0)  return
        nResgrid=nResgrid+1
        cell_NS(nResgrid)%n=>c
    end if    
return
end subroutine IDNSgrid

recursive subroutine IDBCgrid(c)
    implicit none    
    type(grid),pointer::c

    if(associated(c%son1))then
        call IDBCgrid(c%son1)
        call IDBCgrid(c%son2)
        call IDBCgrid(c%son3)
        call IDBCgrid(c%son4)       
    else if(c%sort==1)then
        if (c%cross==BCcorner)  return
        nBCgrid=nBCgrid+1
        cell_BC(nBCgrid)%n=>c
    else if(c%sort==0)then
        if (c%cross==BCcorner)  return
        if( c%cross<0 ) then
            nBCgrid=nBCgrid+1
            cell_BC(nBCgrid)%n=>c
        end if    
    end if
return
end subroutine IDBCgrid

end module solution
