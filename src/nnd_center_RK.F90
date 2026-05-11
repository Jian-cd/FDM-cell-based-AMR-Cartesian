module nndcenter
    use frame
    implicit none
    contains

function minmod(x,y)
    implicit none
    real( prec )::x,y,minmod
    real a
    a=sign(1.d0,x)+sign(1.d0,y)
    minmod=0.5d0*a*min(abs(x),abs(y))
return
end function minmod

subroutine LXL(c,Uip,Uim,Ujp,Ujm)
    implicit none
    type(grid),pointer::c
    real( prec )::dx,dy
    real( prec )::Uip(Nvar),Uim(Nvar),Ujp(Nvar),Ujm(Nvar),fu(4),U(4),Ut(4)
    integer i
    
    dy=h/2**c%lvl
    dx=x_factor*dy

    U(1)=c%U(1)
    U(2)=c%U(1)*c%U(2)
    U(3)=c%U(1)*c%U(3)
    U(4)=c%U(5)/(gama-1.0d0)+0.5d0*c%U(1)*(c%U(2)*c%U(2)+c%U(3)*c%U(3))
    call Res_LXL(c,dx,dy,Uip,Uim,Ujp,Ujm,fu)  
    
    do i=1,4
        Ut(i)=U(i)+fu(i)*tstep
        fu(i)=fu(i)*dx*dy
        if(abs(fu(i)) .gt. Res_max(i))  Res_max(i)=abs(fu(i))  
        Res_rms(i)=Res_rms(i)+fu(i)*fu(i)   
    end do   
    
    c%Ut(1)=Ut(1)
    c%Ut(2)=Ut(2)/Ut(1)
    c%Ut(3)=Ut(3)/Ut(1)
    c%Ut(4)=(Ut(4)-0.5d0*Ut(1)*(c%Ut(2)**2+c%Ut(3)**2))/(cv*Ut(1))
    c%Ut(5)=(gama-1.0d0)*c%Ut(1)*(Ut(4)/c%Ut(1)-0.5d0*(c%Ut(2)*c%Ut(2)+c%Ut(3)*c%Ut(3)))
return
end subroutine LXL

subroutine Res_LXL(c,dx,dy,Uip,Uim,Ujp,Ujm,fu)
    implicit none
    type(grid),pointer::c, cup, cdown
    real( prec )::dx,dy
    real( prec )::Uip(Nvar),Uim(Nvar),Ujp(Nvar),Ujm(Nvar),Uipp(Nvar),Uimm(Nvar),Ujpp(Nvar),Ujmm(Nvar)
    real( prec )::Uippp(Nvar),Uimmm(Nvar),Ujppp(Nvar),Ujmmm(Nvar)
    real( prec )::fu(4),dE(4),dF(4),dEv(4),dFv(4)
    integer i
    
    call Tipp(c,Uipp); call Timm(c,Uimm); call Tjpp(c,Ujpp); call Tjmm(c,Ujmm)
    
    call LXL_X_wcns3(c,dx,Uip,Uim,Uipp,Uimm,dE)
    call LXL_Y_wcns3(c,dy,Ujp,Ujm,Ujpp,Ujmm,dF)
    
    do i=1,4
        fu(i)=-dE(i)-dF(i) 
    end do
return
end subroutine Res_LXL

subroutine LXL_X_wcns3(c, dx, Uip, Uim, Uipp, Uimm, dE)
    implicit none
    type(grid), pointer :: c
    real(prec), intent(in) :: dx
    real(prec), intent(in) :: Uip(Nvar), Uim(Nvar), Uipp(Nvar), Uimm(Nvar)
    real(prec), intent(out) :: dE(4)
    real(prec), dimension(4,5) :: U, Enew, w_i, g_i
    real(prec), dimension(5,5) :: U0
    real(prec), dimension(4,4) :: R_i, Rinv_i
    real(prec), dimension(4) :: w_L, w_R, g_L, g_R, flux_L, flux_R
    integer :: n
   
    U0(1:5,1) = Uimm(1:5)
    U0(1:5,2) = Uim(1:5)
    U0(1:5,3) = c%U(1:5)
    U0(1:5,4) = Uip(1:5)
    U0(1:5,5) = Uipp(1:5)

    call compute_E(U0, Enew, U)

    call REcalc((U(:,2) + U(:,3)) / 2.0, R_i, Rinv_i)
    w_i = matmul(Rinv_i, U)
    g_i = matmul(Rinv_i, Enew)

    call WCNS_Z3(dx, w_i(:,1), w_i(:,2), w_i(:,3), w_i(:,4), R_i, w_L)
    call WCNS_Z3(dx, w_i(:,4), w_i(:,3), w_i(:,2), w_i(:,1), R_i, w_R)
    call WCNS_Z3(dx, g_i(:,1), g_i(:,2), g_i(:,3), g_i(:,4), R_i, g_L)
    call WCNS_Z3(dx, g_i(:,4), g_i(:,3), g_i(:,2), g_i(:,1), R_i, g_R)

    flux_L = 0.5 * (g_L + g_R - lambda_x * (w_R - w_L))

    call REcalc((U(:,3) + U(:,4)) / 2.0, R_i, Rinv_i)
    w_i = matmul(Rinv_i, U)
    g_i = matmul(Rinv_i, Enew)

    call WCNS_Z3(dx, w_i(:,2), w_i(:,3), w_i(:,4), w_i(:,5), R_i, w_L)
    call WCNS_Z3(dx, w_i(:,5), w_i(:,4), w_i(:,3), w_i(:,2), R_i, w_R)
    call WCNS_Z3(dx, g_i(:,2), g_i(:,3), g_i(:,4), g_i(:,5), R_i, g_L)
    call WCNS_Z3(dx, g_i(:,5), g_i(:,4), g_i(:,3), g_i(:,2), R_i, g_R)

    flux_R = 0.5d0 * (g_L + g_R - lambda_x * (w_R - w_L))

    do n = 1, 4
        dE(n) = 4.0d0/(3.0d0*dx) * (flux_R(n) - flux_L(n)) &
              - 1.0d0/(6.0d0*dx) * (Enew(n,4) - Enew(n,2))
    end do
end subroutine LXL_X_wcns3

subroutine REcalc(U,R_i,Rinv_i)
    real( prec ),dimension(4,4) :: R_i,Rinv_i
    real( prec ),dimension(5) :: U0
    real( prec ),dimension(4) :: U
    real( prec ) :: a,q2
    
    U0(1)=U(1)
    U0(2)=U(2)/U(1)
    U0(3)=U(3)/U(1)
    U0(4)=(U(4)-0.5d0*U0(1)*(U0(2)**2+U0(3)**2))/(cv*U0(1))
    U0(5)=(gama-1.0d0)*U0(1)*(U(4)/U0(1)-0.5d0*(U0(2)*U0(2)+U0(3)*U0(3)))
    
    a = sqrt(gama*U0(5)/U0(1))
    
    q2 = 0.5*(U0(2)*U0(2)+U0(3)*U0(3))

    R_i(1,1) = 1.0
    R_i(2,1) = U0(2)-a
    R_i(3,1) = U0(3)
    R_i(4,1) = a*a/(gama - 1) + q2 - U0(2)*a

    R_i(1,2) = 1.0
    R_i(2,2) = U0(2)
    R_i(3,2) = U0(3)
    R_i(4,2) = q2

    R_i(1,3) = 1.0
    R_i(2,3) = U0(2) + a
    R_i(3,3) = U0(3)
    R_i(4,3) = a*a/(gama - 1) + q2 + U0(2)*a

    R_i(1,4) = 0.0
    R_i(2,4) = 0.0
    R_i(3,4) = -1.0
    R_i(4,4) = -U0(3)
    
    Rinv_i(1,1) = ((gama-1.0)*q2 + a*U0(2))/(2.0*a*a)
    Rinv_i(2,1) = (a*a - (gama-1.0)*q2)/(a*a)
    Rinv_i(3,1) = ((gama-1.0)*q2 - a*U0(2))/(2.0*a*a)
    Rinv_i(4,1) = U0(3)

    Rinv_i(1,2) = ((1.0-gama)*U0(2) - a)/(2.0*a*a)
    Rinv_i(2,2) = ((gama-1.0)*U0(2))/(a*a)
    Rinv_i(3,2) = ((1.0-gama)*U0(2) + a)/(2.0*a*a)
    Rinv_i(4,2) = 0.0

    Rinv_i(1,3) = ((1.0-gama)*U0(3))/(2.0*a*a)
    Rinv_i(2,3) = ((gama-1.0)*U0(3))/(a*a)
    Rinv_i(3,3) = ((1.0-gama)*U0(3))/(2.0*a*a)
    Rinv_i(4,3) = -1.0

    Rinv_i(1,4) = (gama-1.0)/(2.0*a*a)
    Rinv_i(2,4) = (1.0-gama)/(a*a)
    Rinv_i(3,4) = (gama-1.0)/(2.0*a*a)
    Rinv_i(4,4) = 0.0
end subroutine REcalc

subroutine compute_E(U0,Enew,U)
    implicit none
    real( prec ),dimension(5,5) :: U0
    real( prec ),dimension(4,5) :: Enew,U
    
    U(1,:) = U0(1,:)
    U(2,:) = U0(1,:)*U0(2,:)
    U(3,:) = U0(1,:)*U0(3,:)
    U(4,:) = U0(5,:)/(gama-1.0d0)+0.5d0*U0(1,:)*(U0(2,:)*U0(2,:)+U0(3,:)*U0(3,:))
    
    Enew(1,:) = U0(1,:)*U0(2,:)
    Enew(2,:) = U0(1,:)*U0(2,:)*U0(2,:)+U0(5,:)
    Enew(3,:) = U0(1,:)*U0(2,:)*U0(3,:)
    Enew(4,:) = U0(2,:)*(U0(5,:)/(gama-1.0d0)+0.5d0*U0(1,:)*(U0(2,:)*U0(2,:)+U0(3,:)*U0(3,:))+U0(5,:))
end subroutine compute_E

subroutine LXL_Y_wcns3(c,dy,Ujp,Ujm,Ujpp,Ujmm,dF)
    implicit none
    type(grid), pointer :: c
    real(prec), intent(in) :: dy
    real(prec), intent(in) :: Ujp(Nvar), Ujm(Nvar), Ujpp(Nvar), Ujmm(Nvar)
    real(prec), intent(out) :: dF(4)
    real(prec), dimension(4,5) :: U, Fnew, w_j, g_j
    real(prec), dimension(5,5) :: U0
    real(prec), dimension(4,4) :: R_j, Rinv_j
    real(prec), dimension(4) :: w_L, w_R, g_L, g_R, flux_L, flux_R
    integer :: n
    
    U0(1:5,1)=Ujmm(1:5)
    U0(1:5,2)=Ujm(1:5)
    U0(1:5,3)=c%U(1:5)
    U0(1:5,4)=Ujp(1:5)
    U0(1:5,5)=Ujpp(1:5)
    
    call compute_F(U0, Fnew, U)
       
    call RFcalc((U(:,2) + U(:,3)) / 2.0, R_j, Rinv_j)
    w_j = matmul(Rinv_j, U)
    g_j = matmul(Rinv_j, Fnew)
    
    call WCNS_Z3(dy, w_j(:,1), w_j(:,2), w_j(:,3), w_j(:,4), R_j, w_L)
    call WCNS_Z3(dy, w_j(:,4), w_j(:,3), w_j(:,2), w_j(:,1),R_j, w_R)
    call WCNS_Z3(dy, g_j(:,1), g_j(:,2), g_j(:,3), g_j(:,4), R_j, g_L)
    call WCNS_Z3(dy, g_j(:,4), g_j(:,3), g_j(:,2), g_j(:,1), R_j, g_R)
    
    flux_L = 0.5 * (g_L + g_R - lambda_y * (w_R - w_L))
    
    call RFcalc((U(:,3) + U(:,4)) / 2, R_j, Rinv_j)
    w_j = matmul(Rinv_j, U)
    g_j = matmul(Rinv_j, Fnew)
    
    call WCNS_Z3(dy, w_j(:,2), w_j(:,3), w_j(:,4), w_j(:,5), R_j, w_L)
    call WCNS_Z3(dy, w_j(:,5), w_j(:,4), w_j(:,3), w_j(:,2), R_j, w_R)
    call WCNS_Z3(dy, g_j(:,2), g_j(:,3), g_j(:,4), g_j(:,5), R_j, g_L)
    call WCNS_Z3(dy, g_j(:,5), g_j(:,4), g_j(:,3), g_j(:,2), R_j, g_R)
    
    flux_R = 0.5 * (g_L + g_R - lambda_y * (w_R - w_L))
    
    do n=1, 4
        dF(n)=4.d0/3.d0/dy*(flux_R(n) - flux_L(n)) - 1.d0/6.d0/dy*(Fnew(n,4) - Fnew(n,2))
    end do
    
    return
end subroutine LXL_Y_wcns3

subroutine RFcalc(U,R_j,Rinv_j)
    real( prec ),dimension(4,4) :: R_j,Rinv_j 
    real( prec ),dimension(5) :: U0
    real( prec ),dimension(4) :: U
    real( prec ) :: a,q2
    
    U0(1)=U(1)
    U0(2)=U(2)/U(1)
    U0(3)=U(3)/U(1)
    U0(4)=(U(4)-0.5d0*U0(1)*(U0(2)**2+U0(3)**2))/(cv*U0(1))
    U0(5)=(gama-1.0d0)*U0(1)*(U(4)/U0(1)-0.5d0*(U0(2)*U0(2)+U0(3)*U0(3)))
    
    a = sqrt(gama*U0(5)/U0(1))
    
    q2 = 0.5*(U0(2)*U0(2)+U0(3)*U0(3))

    R_j(1,1) = 1.0
    R_j(2,1) = U0(2)
    R_j(3,1) = U0(3)-a
    R_j(4,1) = a*a/(gama - 1) + q2 - U0(3)*a

    R_j(1,2) = 1.0
    R_j(2,2) = U0(2)
    R_j(3,2) = U0(3)
    R_j(4,2) = q2

    R_j(1,3) = 1.0
    R_j(2,3) = U0(2)
    R_j(3,3) = U0(3) + a
    R_j(4,3) = a*a/(gama - 1) + q2 + U0(3)*a

    R_j(1,4) = 0.0
    R_j(2,4) = 1.0
    R_j(3,4) = 0.0
    R_j(4,4) = U0(2)
    
    Rinv_j(1,1) = ((gama-1.0)*q2 + a*U0(3))/(2.0*a*a)
    Rinv_j(2,1) = (a*a - (gama-1.0)*q2)/(a*a)
    Rinv_j(3,1) = ((gama-1.0)*q2 - a*U0(3))/(2.0*a*a)
    Rinv_j(4,1) = -U0(2)

    Rinv_j(1,2) = ((1.0-gama)*U0(2))/(2.0*a*a)
    Rinv_j(2,2) = ((gama-1.0)*U0(2))/(a*a)
    Rinv_j(3,2) = ((1.0-gama)*U0(2))/(2.0*a*a)
    Rinv_j(4,2) = 1.0

    Rinv_j(1,3) = ((1.0-gama)*U0(3) - a)/(2.0*a*a)
    Rinv_j(2,3) = ((gama-1.0)*U0(3))/(a*a)
    Rinv_j(3,3) = ((1.0-gama)*U0(3) + a)/(2.0*a*a)
    Rinv_j(4,3) = 0.0

    Rinv_j(1,4) = (gama-1.0)/(2.0*a*a)
    Rinv_j(2,4) = (1.0-gama)/(a*a)
    Rinv_j(3,4) = (gama-1.0)/(2.0*a*a)
    Rinv_j(4,4) = 0.0
end subroutine RFcalc
    
subroutine compute_F(U0,Fnew,U)
    implicit none
    real( prec ),dimension(5,5) :: U0
    real( prec ),dimension(4,5) :: Fnew,U
    
    U(1,:) = U0(1,:)
    U(2,:) = U0(1,:)*U0(2,:)
    U(3,:) = U0(1,:)*U0(3,:)
    U(4,:) = U0(5,:)/(gama-1.0d0)+0.5d0*U0(1,:)*(U0(2,:)*U0(2,:)+U0(3,:)*U0(3,:))
    
    Fnew(1,:) = U0(1,:)*U0(3,:)
    Fnew(2,:) = U0(1,:)*U0(2,:)*U0(3,:)
    Fnew(3,:) = U0(1,:)*U0(3,:)*U0(3,:)+U0(5,:)
    Fnew(4,:) = U0(3,:)*(U0(5,:)/(gama-1.0d0)+0.5d0*U0(1,:)*(U0(2,:)*U0(2,:)+U0(3,:)*U0(3,:))+U0(5,:))
end subroutine compute_F  


subroutine WCNS_Z3(dx, u1, u2, u3, u4, R_i, hp_i)
    implicit none
    real( prec ), intent(in) :: dx 
    real( prec ), dimension(4), intent(in) :: u1, u2, u3, u4
    real( prec ), dimension(4,4), intent(in) :: R_i
    real( prec ), dimension(4), intent(out) :: hp_i
    real( prec ), dimension(4) :: IS1, IS2, IS3, tau3, a1, a2, w1, w2, p1, p2, uL
    real( prec ) :: d1, d2, b1, b2, xc1, xc2, a, b
    real( prec ), parameter :: c11 = 1.2, c21 = 0.1, c31 = 25., c12 = 1.2, c22 = 0.1, c32 = 35.
    real( prec ), parameter :: epsilon=1.e-40,epsilon2=1.e-6,epsilon1=1.e-12
    integer :: i
    
    d1 = 1./4.; d2 = 3./4.
    p1 = -1./2.*u1  + 3./2.*u2   
    p2 =  1./2.*u2  + 1./2.*u3
    
    IS1 = (-u1 + u2)**2
    IS2 = (-u2 + u3)**2
    IS3 = (u1 - 2.*u2 + u3)**2.0

    do i = 1, 4
        d1 = 1./11.; d2 = 10./11.
        p1(i) = u2(i)
        p2(i) = ((-u1(i) + 6.0d0 * u2(i) + 3.0d0 * u3(i)) / 8.0d0 - d1 * p1(i))/d2
            
        IS1(i) = min(IS1(i), IS2(i))
        IS2(i) = IS3(i) + u1(i)**2 + u2(i)**2 + u3(i)**2 - u1(i)*u2(i) - u2(i)*u3(i) - u1(i)*u3(i)
            
        tau3(i) = (10.*abs(IS2(i) - 3.*IS1(i)))**2
            
        a1(i) = d1 * (1.d0 + tau3(i)/(epsilon + IS1(i)))
        a2(i) = d2 * (1.d0 + tau3(i)/(epsilon + IS2(i)))
            
        w1(i) = a1(i) / (a1(i) + a2(i))
        w2(i) = a2(i) / (a1(i) + a2(i))
                
        UL(i) = w1(i) * p1(i) + w2(i) * p2(i)
    end do
        
    hp_i(1) = R_i(1,1) * uL(1) + R_i(1,2) * uL(2) + R_i(1,3) * uL(3) + R_i(1,4) * uL(4)
    hp_i(2) = R_i(2,1) * uL(1) + R_i(2,2) * uL(2) + R_i(2,3) * uL(3) + R_i(2,4) * uL(4)
    hp_i(3) = R_i(3,1) * uL(1) + R_i(3,2) * uL(2) + R_i(3,3) * uL(3) + R_i(3,4) * uL(4)
    hp_i(4) = R_i(4,1) * uL(1) + R_i(4,2) * uL(2) + R_i(4,3) * uL(3) + R_i(4,4) * uL(4)
end subroutine WCNS_Z3
   
end module nndcenter
