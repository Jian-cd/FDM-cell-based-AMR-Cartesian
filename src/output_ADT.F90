module outfile
    use Global_Var
    implicit none
    real( prec ),allocatable::cVariables(:,:)
    contains

subroutine output_Res
    implicit none

    write(*, "(I8,I8,6E16.6)") it, nResgrid, tstep, tt, Res_rms(1:Neq)
    open(99,file='./Output/'//"Residual.dat",position="append")
    write(99,"(I8,I8,12E24.10)") it, nResgrid, tstep, tt, Res_max(1:Neq), Res_rms(1:Neq)
    close(99) 
return
end  subroutine output_Res    
!======================================================================

subroutine output_error  
    implicit none    
    integer             :: i, nResgrid_0, nResgrid_patch
    integer,save        :: count=1
    logical             :: ios
    type(grid),pointer  :: c
    real( prec )        :: L1_u, L2_u, Loo_u, L1_rho, L2_rho, Loo_rho
    real( prec )        :: P1_u, P2_u, Poo_u, P1_rho, P2_rho, Poo_rho
    real( prec )        :: L1_rho0, L2_rho0, Loo_rho0, L1_rho_patch, L2_rho_patch, Loo_rho_patch
    real( prec )        :: Unisum0, Unisum_patch, area, abs_rhorho
    real( prec )        :: L1_p, L2_p, Loo_p
    real( prec )        :: dx, dy, uu, pp, Unisum, Patchsum, xc, yc, rr, rhorho, d, dT
    real( prec ),parameter:: cc=0.02d0, Rc2=1.d0, cc1=5.d0
    
    open (83,file='./Output/accuracy_error.dat',position="append")
    
    if(count==1) write(83,*) 'nstep Ngrids L1_rho P1_rho L2_rho P2_rho Loo_rho Poo_rho '    
        
    xc=0.5*(cell(1)%center%x+cell(m)%center%x)
    yc=0.5*(cell(m)%center%y+cell(total)%center%y)      
        
    L1_rho = 0.d0; L2_rho = 0.d0; Loo_rho = 0.d0
    P1_rho = 0.d0; P2_rho = 0.d0; Poo_rho = 0.d0
    Unisum = 0.d0; Patchsum = 0.d0; ios = .false.
        
    do i=1,nResgrid
        c => cell_NS(i)%n
        dy = h/2**c%lvl
        dx = x_factor*dy  
        Unisum = Unisum+dx*dy
        rr = (c%center%x-5.d0)**2 + (c%center%y-5.d0)**2
        dT = -(gama-1.d0)*(cc1*cc1)*exp(1.d0-rr)/(8.d0*gama*Pi*Pi)
        
        rhorho = c%U(1) - (1.d0 + dT)**(1.d0/(gama - 1.d0))
            
        L1_rho = L1_rho + abs(rhorho)*dx*dy; L2_rho = L2_rho + rhorho*rhorho*dx*dy; Loo_rho = max(Loo_rho, abs(rhorho))
        if (c%lvl/=0) then
            ios=.true.
            Patchsum=Patchsum+dx*dy
            P1_rho = P1_rho + abs(rhorho)*dx*dy; P2_rho = P2_rho + rhorho*rhorho*dx*dy; Poo_rho = max(Poo_rho, abs(rhorho))
        end if
    end do
    
    L1_rho = L1_rho/Unisum; L2_rho = sqrt(L2_rho/Unisum)
    if (ios) then
        P1_rho = P1_rho/Patchsum; P2_rho = sqrt(P2_rho/Patchsum)
    end if
    write(83,"(2I8,6E14.6)") it, nResgrid, L1_rho, P1_rho, L2_rho, P2_rho, Loo_rho, Poo_rho
    
    write(83,"(I8, I8, 3E14.6, I8, 3E14.6, I8, 3E14.6)") it, nResgrid, &
        L1_rho, L2_rho, Loo_rho, &
        nResgrid_0, L1_rho0, L2_rho0, Loo_rho0, &
        nResgrid_patch, L1_rho_patch, L2_rho_patch, Loo_rho_patch
    
    close (83)  
    count=count+1
return
end subroutine output_error

end module outfile
