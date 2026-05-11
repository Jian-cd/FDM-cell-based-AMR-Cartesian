module frame
    use Global_Var 
    implicit none

    interface operator(.i.) 
        module procedure brinv 
    end interface

contains

    subroutine add_to_neighbor_list(new_neighbor, neighbors)
        type(p_array), allocatable, intent(inout) :: neighbors(:)
        type(grid), pointer, intent(in) :: new_neighbor
        type(p_array), allocatable :: temp(:)
        integer :: n
    
        if (.not. allocated(neighbors)) then
          allocate(neighbors(1))
          neighbors(1)%n => new_neighbor
        else
            do n = 1, size(neighbors)
                if (associated(neighbors(n)%n, new_neighbor)) return
            end do
      
            allocate(temp(size(neighbors)+1))
            temp(1:size(neighbors)) = neighbors
            temp(size(neighbors)+1)%n => new_neighbor
            call move_alloc(temp, neighbors)
        end if
    end subroutine add_to_neighbor_list


    function get_direction(xx, yy, x, y) result(dir)
        real(prec) :: xx, yy, x, y
        real(prec) :: dx, dy, angle
        integer :: dir
    
        dx = xx - x
        dy = yy - y
        angle = (atan2(dy, dx) + 1e-10) * 180.0 / pi 
    
        dir = mod(floor((angle + 22.5)/45.0 + 8), 8) + 1
    end function

    subroutine parent_value_weno(c, var)
        type(grid), pointer, intent(in) :: c
        real(prec), intent(out) :: var(nvar)  
        type(p_array), allocatable :: neighbors1(:),neighbors2(:)
        type(grid), pointer :: neighbor
        real(prec) :: x, y, xmin, xmax, ymin, ymax
        integer :: i, j, m, m1, m2
        logical :: dir_used(8)
    
        dir_used = .false.
    
        x = c%center%x
        y = c%center%y
    
        call search_cell_5(c%son1, x, y, neighbors1, neighbors2, dir_used)
        call search_cell_5(c%son2, x, y, neighbors1, neighbors2, dir_used)
        call search_cell_5(c%son3, x, y, neighbors1, neighbors2, dir_used)
        call search_cell_5(c%son4, x, y, neighbors1, neighbors2, dir_used)
    
        m1 = size(neighbors1); m2 = size(neighbors2)
        m = m1 + m2
        if (m < 8) then
            do i = 1, 2
                do j = 1, 2
                    if (i == 1 .and. j == 1) then
                        neighbor => c%leftnei
                    else if (i == 1 .and. j == 2) then
                        neighbor => c%rightnei
                    else if (i == 2 .and. j == 1) then
                        neighbor => c%downnei
                    else if (i == 2 .and. j == 2) then
                        neighbor => c%upnei
                    end if
                    if (associated(neighbor)) then
                        call search_cell_5(neighbor, x, y, neighbors1, neighbors2, dir_used)
                    end if
                end do
            end do
        end if
    
        m1 = size(neighbors1); m2 = size(neighbors2)
        m = m1 + m2
        if ( m < 7 ) then
            call interpA(c, var)
        else    
            call WLSI_US_WENO(c, x, y, neighbors1, neighbors2, m, m1, m2, var)
        end if
        if (allocated(neighbors1)) deallocate(neighbors1)
        if (allocated(neighbors2)) deallocate(neighbors2)
    end subroutine parent_value_weno

    subroutine child_value_weno(xx, yy, c, var)
        real(prec), intent(in) :: xx, yy
        type(grid), pointer, intent(in) :: c
        real(prec), intent(out) :: var(nvar)
        type(p_array), allocatable :: neighbors1(:), neighbors2(:)
        type(grid), pointer :: neighbor
        real(prec) :: x, y, xmin, xmax, ymin, ymax
        integer :: i, j, m, m1, m2
        logical :: dir_used(8)
    
        dir_used = .false.
    
        x = c%center%x
        y = c%center%y
    
        call add_to_neighbor_list(c, neighbors1)
    
        do i = -1, 1
            do j = -1, 1
                if (i == 0 .and. j == 0) cycle 
            
                if (i == -1 .and. j == -1) then
                    neighbor => c%downleftnei
                else if (i == -1 .and. j == 0) then
                    neighbor => c%leftnei
                else if (i == -1 .and. j == 1) then
                    neighbor => c%upleftnei
                else if (i == 0 .and. j == -1) then
                    neighbor => c%downnei
                else if (i == 0 .and. j == 1) then
                    neighbor => c%upnei
                else if (i == 1 .and. j == -1) then
                    neighbor => c%downrightnei
                else if (i == 1 .and. j == 0) then
                    neighbor => c%rightnei
                else if (i == 1 .and. j == 1) then
                    neighbor => c%uprightnei
                end if
                if (associated(neighbor)) then
                    call search_cell_5(neighbor, x, y, neighbors1, neighbors2, dir_used)
                end if
            end do
        end do
    
        m1 = size(neighbors1); m2 = size(neighbors2)
        m = m1 + m2
    
        if ( m<7 ) then
            var(1:7)=c%U(1:7)
        else
            call WLSI_US_WENO(c, xx, yy, neighbors1, neighbors2, m, m1, m2, var)
        end if
        if (allocated(neighbors1)) deallocate(neighbors1)
        if (allocated(neighbors2)) deallocate(neighbors2)
    end subroutine child_value_weno

    recursive subroutine search_cell_5(c, x, y, neighbors1, neighbors2, dir_used)
        type(grid), pointer :: c
        real(prec) :: x, y
        integer :: i, dir
        type(p_array), allocatable, intent(inout) :: neighbors1(:), neighbors2(:)
        logical :: dir_used(8)
    
        if (.not. associated(c%son1)) then 
            dir = get_direction(c%center%x, c%center%y, x, y)
            if (.not. dir_used(dir)) then 
                if (any((/1,5,3,7/) == dir)) then
                    call add_to_neighbor_list(c, neighbors1)
                else
                    call add_to_neighbor_list(c, neighbors2)
                end if
                dir_used(dir) = .true.
            end if
        else  
            do i = 1, 4
                call search_cell_5(c%son1, x, y, neighbors1, neighbors2, dir_used)
                call search_cell_5(c%son2, x, y, neighbors1, neighbors2, dir_used)
                call search_cell_5(c%son3, x, y, neighbors1, neighbors2, dir_used)
                call search_cell_5(c%son4, x, y, neighbors1, neighbors2, dir_used)
            end do
        end if
    end subroutine

    subroutine WLSI_US_WENO(c, xx, yy, neighbors1, neighbors2, m, m1, m2, var)
        implicit none
        integer, intent(in) :: m, m1, m2
        type(p_array), intent(in) :: neighbors1(:), neighbors2(:)
        type(grid), pointer :: c
        real(prec) :: xx, yy
        real(prec) :: var(Nvar)
    
        real(prec) :: cvar(m, Nvar), x(m), y(m), d(m), w(m), Phi(m,1), W_Matrix(m,m)
        real(prec), parameter :: zero = 0.d0, one = 1.d0, epsilon = 1.d-6
        integer :: i, k
    
        real(prec) :: P_large(m, 6), WP_large(m, 6), PTWP_large(6, 6), WPhi_large(m, 1), PTWPhi_large(6, 1)
        real(prec) :: arr_large(6, 1), var_large(Nvar), inv_PTWP_large(6, 6)
    
        real(prec) :: P_small1(m1, 3), WP_small1(m1, 3), PTWP_small1(3, 3), WPhi_small1(m1, 1), PTWPhi_small1(3, 1)
        real(prec) :: arr_small1(3, 1), inv_PTWP_small1(3, 3), var_small(2, Nvar), W_small1(m1,m1), Phi_small1(m1, 1)
    
        real(prec) :: P_small2(m2, 3), WP_small2(m2, 3), PTWP_small2(3, 3), WPhi_small2(m2, 1), PTWPhi_small2(3, 1)
        real(prec) :: arr_small2(3, 1), inv_PTWP_small2(3, 3), W_small2(m2,m2), Phi_small2(m2, 1) 
    
        real(prec) :: beta(3), alpha(3), omega(3)
        real(prec) :: smoothness_large(Nvar), smoothness_small(2, Nvar), tau
        real(prec), parameter :: linear_weights(3) = [0.98_prec, 0.01_prec, 0.01_prec]  
    
        integer :: stencil_indices1(m1), stencil_indices2(m2)
    
        do i = 1, m1
            cvar(i, 1:7) = neighbors1(i)%n%U(1:7)
            x(i) = neighbors1(i)%n%center%x
            y(i) = neighbors1(i)%n%center%y
        
            d(i) = sqrt((x(i) - xx)**2 + (y(i) - yy)**2)
            w(i) = 1.d0 / (d(i) + epsilon)
        end do
    
        do i = 1, m2
            cvar(m1+i, 1:7) = neighbors2(i)%n%U(1:7)
            x(m1+i) = neighbors2(i)%n%center%x
            y(m1+i) = neighbors2(i)%n%center%y
        
            d(m1+i) = sqrt((x(m1+i) - xx)**2 + (y(m1+i) - yy)**2)
            w(m1+i) = 1.d0 / (d(m1+i) + epsilon)
        end do
    
        P_large(:, 1) = one
        P_large(:, 2) = x(:)
        P_large(:, 3) = y(:)
        P_large(:, 4) = x(:) * y(:)
        P_large(:, 5) = x(:) * x(:)
        P_large(:, 6) = y(:) * y(:)
    
        W_Matrix = zero
        forall(i=1:m) W_Matrix(i,i) = w(i)
        WP_large = matmul(W_Matrix, P_large)
        PTWP_large = matmul(transpose(P_large), WP_large)
        inv_PTWP_large = .i.PTWP_large
    
        P_small1(:, 1) = one
        P_small1(:, 2) = x(1:m1)
        P_small1(:, 3) = y(1:m1)
    
        W_small1 = zero
        forall(i=1:m1) W_small1(i,i) = w(i)
        WP_small1 = matmul(W_small1, P_small1)
        PTWP_small1 = matmul(transpose(P_small1), WP_small1)
        inv_PTWP_small1 = .i.PTWP_small1
        
        P_small2(:, 1) = one
        P_small2(:, 2) = x(m1+1:m)
        P_small2(:, 3) = y(m1+1:m)
    
        W_small2 = zero
        forall(i=1:m2) W_small2(i,i) = w(m1+i)
        WP_small2 = matmul(W_small2, P_small2)
        PTWP_small2 = matmul(transpose(P_small2), WP_small2)
        inv_PTWP_small2 = .i.PTWP_small2
    
        do k = 1, Nvar
            if (k == 4) cycle 
            select case(k)
            case(1)
                Phi(:,1) = cvar(:,1)
            case(2)
                Phi(:,1) = cvar(:,2)
            case(3)
                Phi(:,1) = cvar(:,3)
            case(5)
                Phi(:,1) = cvar(:,5)
            case default
                cycle
            end select
        
            WPhi_large = matmul(W_Matrix, Phi)
            PTWPhi_large = matmul(transpose(P_large), WPhi_large)
            arr_large = matmul(inv_PTWP_large, PTWPhi_large)
        
            var_large(k) = arr_large(1, 1) + xx * arr_large(2, 1) + yy * arr_large(3, 1) + &
                           xx * yy * arr_large(4, 1) + xx*xx * arr_large(5, 1) + yy*yy * arr_large(6, 1)
            smoothness_large(k) = compute_smoothness(c, arr_large, 2)
        
            Phi_small1(:, 1) = Phi(1:m1, 1)
            WPhi_small1 = matmul(W_small1, Phi_small1)
            PTWPhi_small1 = matmul(transpose(P_small1), WPhi_small1)
            arr_small1 = matmul(inv_PTWP_small1, PTWPhi_small1)
        
            Phi_small2(:,1) = Phi(m1+1:m,1)
            WPhi_small2 = matmul(W_small2, Phi_small2)
            PTWPhi_small2 = matmul(transpose(P_small2), WPhi_small2)
            arr_small2 = matmul(inv_PTWP_small2, PTWPhi_small2)
        
            var_small(1, k) = arr_small1(1, 1) + xx * arr_small1(2, 1) + yy * arr_small1(3, 1)
            var_small(2, k) = arr_small2(1, 1) + xx * arr_small2(2, 1) + yy * arr_small2(3, 1)
            
            smoothness_small(1, k) = compute_smoothness(c, arr_small1, 1)
            smoothness_small(2, k) = compute_smoothness(c, arr_small2, 1)
        
            beta(1) = smoothness_large(k)
            beta(2:3) = smoothness_small(:,k)
        
            tau = (sum(abs(beta(1) - beta(2:3))) / 2.0)**2
        
            alpha = linear_weights * (1.0d0 + tau / (epsilon + beta))
            omega = alpha / sum(alpha)
    
            var(k) = omega(1) * (1.0d0/ linear_weights(1) * var_large(k) - &
                (linear_weights(2) / linear_weights(1) * var_small(1, k) + linear_weights(3) / linear_weights(1) * var_small(2, k))) + &
                (omega(2) * var_small(1, k) + omega(3) * var_small(2, k))
        end do    
    contains
        function compute_smoothness(c, coeff, order) result(beta)
            type(grid), pointer :: c
            real(prec) :: coeff(:,:)   
            integer :: order     
            real(prec) :: beta

            real(prec) :: x_min, x_max, y_min, y_max, hx, hy
            real(prec) :: a0, a1, a2, a3, a4, a5, b0, b1, b2

            hx = x_factor * h / 2.0_prec**(c%lvl + 1)
            hy = h / 2.0_prec**(c%lvl + 1)
            x_min = c%center%x - hx
            x_max = c%center%x + hx
            y_min = c%center%y - hy
            y_max = c%center%y + hy

            if (order == 2) then
                a0 = coeff(1,1); a1 = coeff(2,1); a2 = coeff(3,1)
                a3 = coeff(4,1); a4 = coeff(5,1); a5 = coeff(6,1)

                beta = 0.0_prec
                beta = beta + integrate_square_poly([a1, 2*a3, a4], [0,1,0], x_min, x_max, y_min, y_max)
                beta = beta + integrate_square_poly([a2, a4, 2*a5], [0,0,1], x_min, x_max, y_min, y_max)
                beta = beta + (2*a3)**2 * (x_max - x_min) * (y_max - y_min)
                beta = beta + a4**2 * (x_max - x_min) * (y_max - y_min)
                beta = beta + (2*a5)**2 * (x_max - x_min) * (y_max - y_min)

            else if (order == 1) then
                b0 = coeff(1,1); b1 = coeff(2,1); b2 = coeff(3,1)
                beta = (b1**2 + b2**2) * (x_max - x_min) * (y_max - y_min)
            end if
        end function

        function integrate_square_poly(k, xy_order, x_min, x_max, y_min, y_max) result(val)
            real(prec) :: k(3)      
            integer :: xy_order(3)  
            real(prec) :: x_min, x_max, y_min, y_max, val
            real(prec) :: term1, term2, term3

            term1 = k(1)**2 * (x_max - x_min) * (y_max - y_min)
            term2 = 2*k(1)*k(2) * (x_max**2 - x_min**2)/2 * (y_max - y_min) + &
                    2*k(1)*k(3) * (x_max - x_min) * (y_max**2 - y_min**2)/2
            term3 = k(2)**2 * (x_max**3 - x_min**3)/3 * (y_max - y_min) + &
                    2*k(2)*k(3) * (x_max**2 - x_min**2)/2 * (y_max**2 - y_min**2)/2 + &
                    k(3)**2 * (x_max - x_min) * (y_max**3 - y_min**3)/3
            val = term1 + term2 + term3
        end function
    end subroutine WLSI_US_WENO


    recursive subroutine interpA(c,var)  
        implicit none
        type(grid),pointer::c
        real( prec )::var(Nvar),son1var(Nvar),son2var(Nvar),son3var(Nvar),son4var(Nvar)
        integer i
  
        if(.not.associated(c%son1))then 
            var(1:7)=c%U(1:7)
        return
        end if
    
        if(associated(c%son1%son1))then 
            call interpA(c%son1,son1var)
        else
            son1var(1:7)=c%son1%U(1:7)
        end if
    
        if(associated(c%son2%son1))then
            call interpA(c%son2,son2var)
        else
            son2var(1:7)=c%son2%U(1:7)
        end if
    
        if(associated(c%son3%son1))then
            call interpA(c%son3,son3var)
        else
            son3var(1:7)=c%son3%U(1:7)  
        end if
    
        if(associated(c%son4%son1))then
            call interpA(c%son4,son4var)
        else
            son4var(1:7)=c%son4%U(1:7)
        end if
    
        do i=1, Nvar  
            var(i)=0.25d0*(son1var(i)+son2var(i)+son3var(i)+son4var(i))
        end do
    return
    end subroutine interpA


    subroutine Tip(c,var) 
        implicit none
        type(grid),pointer  :: c,cn,ct1,ct2,ct3,ct4,ct5,ct6,ct7,ct8,ct9,ct10,ct11,ct12
        real( prec )        :: var(Nvar), result_values(Nvar,1)
        real( prec )        :: xx,yy

        cn=>c%rightnei
        if(cn%lvl==c%lvl)then
            if(.not.associated(cn%son1))then
                var(1:7)=cn%U(1:7)
            else
                call parent_value_weno(cn, var)
            end if
        else
            xx=c%center%x+x_factor*h/2**(c%lvl);  yy=c%center%y
            call child_value_weno(xx, yy, cn, var)
        end if
    return 
    end subroutine Tip


    subroutine Tipp(c,var)      
        implicit none
        type(grid),pointer  :: c,cn,cnn,ct1,ct2,ct3,ct4,ct5,ct6,ct7,ct8,ct9,ct10,ct11,ct12
        real( prec )        :: var(Nvar), result_values(Nvar,1)
        real( prec )        :: xx,yy
    
        cn=>c%rightnei
        if(cn%lvl==c%lvl)then
            cnn=>cn%rightnei
            if(cnn%lvl==cn%lvl)then
                if(.not.associated(cnn%son1))then
                    var(1:7)=cnn%U(1:7)
                else
                    call parent_value_weno(cnn, var)
                end if
            else
                xx=c%center%x+2.0*x_factor*h/2**(c%lvl);  yy=c%center%y
                call child_value_weno(xx, yy, cnn, var)
            end if
        else
            xx=c%center%x+2.0*x_factor*h/2**(c%lvl);  yy=c%center%y
            call child_value_weno(xx, yy, cn, var) 
        end if  
    return
    end subroutine Tipp


    subroutine Tim(c,var)
        implicit none
        type(grid),pointer  :: c,cn,ct1,ct2,ct3,ct4,ct5,ct6,ct7,ct8,ct9,ct10,ct11,ct12
        real( prec )        :: var(Nvar), result_values(Nvar,1)
        real( prec )        :: xx,yy
   
        cn=>c%leftnei   
        if(cn%lvl==c%lvl)then  
            if(.not.associated(cn%son1))then
                var(1:7)=cn%U(1:7)
            else
                call parent_value_weno(cn, var)
            end if
        else 
            xx=c%center%x-x_factor*h/2**(c%lvl);  yy=c%center%y
            call child_value_weno(xx, yy, cn, var)
        end if
    return 
    end subroutine Tim


    subroutine Timm(c,var)               
        implicit none
        type(grid),pointer  :: c,cn,cnn,ct1,ct2,ct3,ct4,ct5,ct6,ct7,ct8,ct9,ct10,ct11,ct12
        real( prec )        :: var(Nvar), result_values(Nvar,1)
        real( prec )        :: xx,yy
  
        cn=>c%leftnei  
        if(cn%lvl==c%lvl)then
            cnn=>cn%leftnei 
            if(cnn%lvl==cn%lvl)then
                if(.not.associated(cnn%son1))then
                    var(1:7)=cnn%U(1:7)
                else
                    call parent_value_weno(cnn, var)
                end if
            else
                xx=c%center%x-2.0*x_factor*h/2**(c%lvl);  yy=c%center%y
                call child_value_weno(xx, yy, cnn, var)
            end if
        else
            xx=c%center%x-2.0*x_factor*h/2**(c%lvl);  yy=c%center%y
            call child_value_weno(xx, yy, cn, var)     
        end if
    return
    end subroutine Timm


    subroutine Tjp(c,var)
        implicit none
        type(grid),pointer  :: c,cn,ct1,ct2,ct3,ct4,ct5,ct6,ct7,ct8,ct9,ct10,ct11,ct12
        real( prec )        :: var(Nvar), result_values(Nvar,1)
        real( prec )        :: xx,yy
   
        cn=>c%upnei
        if(cn%lvl==c%lvl)then  
            if(.not.associated(cn%son1))then
                var(1:7)=cn%U(1:7)
            else
                call parent_value_weno(cn, var)
            end if  
        else
            xx=c%center%x;  yy=c%center%y+h/2**(c%lvl)
            call child_value_weno(xx, yy, cn, var)
        end if
    return 
    end subroutine Tjp


    subroutine Tjpp(c,var)      
        implicit none
        type(grid),pointer  :: c,cn,cnn,ct1,ct2,ct3,ct4,ct5,ct6,ct7,ct8,ct9,ct10,ct11,ct12
        real( prec )        :: var(Nvar), result_values(Nvar,1)
        real( prec )        :: xx,yy
    
        cn=>c%upnei 
        if(cn%lvl==c%lvl)then
            cnn=>cn%upnei 
            if(cnn%lvl==cn%lvl)then
                if(.not.associated(cnn%son1))then
                    var(1:7)=cnn%U(1:7)
                else
                    call parent_value_weno(cnn, var)
                end if
            else
                xx=c%center%x;  yy=c%center%y+2.0*h/2**(c%lvl)
                call child_value_weno(xx, yy, cnn, var)
            end if
        else
            xx=c%center%x;  yy=c%center%y+2.0*h/2**(c%lvl)
            call child_value_weno(xx, yy, cn, var)
        end if
    return
    end subroutine Tjpp


    subroutine Tjm(c,var)
        implicit none
        type(grid),pointer  :: c,cn,ct1,ct2,ct3,ct4,ct5,ct6,ct7,ct8,ct9,ct10,ct11,ct12
        real( prec )        :: var(Nvar), result_values(Nvar,1)
        real( prec )        :: xx,yy
    
        cn=>c%downnei 
        if(cn%lvl==c%lvl)then  
            if(.not.associated(cn%son1))then
                var(1:7)=cn%U(1:7)
            else
                call parent_value_weno(cn, var)
            end if   
        else
            xx=c%center%x;  yy=c%center%y-h/2**(c%lvl)
            call child_value_weno(xx, yy, cn, var)
        end if
    return 
    end subroutine Tjm


    subroutine Tjmm(c,var)               
        implicit none
        type(grid),pointer  :: c,cn,cnn,ct1,ct2,ct3,ct4,ct5,ct6,ct7,ct8,ct9,ct10,ct11,ct12
        real( prec )        :: var(Nvar), result_values(Nvar,1)
        real( prec )        :: xx,yy
    
        cn=>c%downnei 
        if(cn%lvl==c%lvl)then
            cnn=>cn%downnei 
            if(cnn%lvl==cn%lvl)then
                if(.not.associated(cnn%son1))then
                    var(1:7)=cnn%U(1:7)
                else
                    call parent_value_weno(cnn, var)
                end if
            else
                xx=c%center%x;  yy=c%center%y-2.0*h/2**(c%lvl)
                call child_value_weno(xx, yy, cnn, var)
            end if
        else
            xx=c%center%x;  yy=c%center%y-2.0*h/2**(c%lvl)
            call child_value_weno(xx, yy, cn, var)
        end if
    return
    end subroutine Tjmm
  

    function brinv(re) result(r)
        real( prec ),intent(in):: re(:,:)
        real( prec ):: r(size(re,1),size(re,1))
        integer     :: flag,n,i,j,k 
        real( prec ):: t,d 
        integer   is(size(re,1)),js(size(re,1))
    
        integer :: status

        n=size(re,1)
        r=re
        flag=1
        do k=1,n
            d=0.d0
            do i=k,n
                do j=k,n
                    if (abs(r(i,j)).gt.d) then
                        d=abs(r(i,j))
                        is(k)=i
                        js(k)=j
                    end if
                end do
            end do
            if (d+1.0.eq.1.0) then
                flag=0
                write(*,*) "flag=0,real matrix singularity!" 
                return
            end if
            do j=1,n
                t=r(k,j)
                r(k,j)=r(is(k),j)
                r(is(k),j)=t
            end do
            do i=1,n
                t=r(i,k)
                r(i,k)=r(i,js(k))
                r(i,js(k))=t
            end do
            r(k,k)=1/r(k,k)
            do j=1,n
                if (j.ne.k) r(k,j)=r(k,j)*r(k,k)
            end do
            do i=1,n
                if (i.ne.k) then
                    do j=1,n
                        if (j.ne.k) r(i,j)=r(i,j)-r(i,k)*r(k,j)
                    end do
                end if
            end do
            do i=1,n
                if (i.ne.k) r(i,k)=-r(i,k)*r(k,k)
            end do
        end do
    
        do k=n,1,-1
            do j=1,n
                t=r(k,j)
                r(k,j)=r(js(k),j)
                r(js(k),j)=t
            end do
            do i=1,n
                t=r(i,k)
                r(i,k)=r(i,is(k))
                r(i,is(k))=t
            end do
        end do
    
    end function

end module frame
