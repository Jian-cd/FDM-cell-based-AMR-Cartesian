module const_var
    implicit none

    integer,parameter       :: single_prec = KIND(1.0e0)
    integer,parameter       :: double_prec = KIND(1.0d0)
    integer,parameter       ::        prec = double_prec
    integer,parameter       :: BCcyc=6, BCcorner=-5
    integer,parameter       :: Nvar=9
    real( prec ),parameter  :: Pi=3.1415926535897932d0
end module const_var

    
module Global_Var
    use const_var
    implicit none
    real( prec )            :: Xmin, Xmax, Ymin, Ymax
    real( prec )            :: h, x_factor = 1.0
    integer                 :: m, n, total
    
    type point   
        real( prec )        :: x,y
    end type

    type grid
        integer             :: num(2)
        integer             :: lvl,sort,cross
        type(point)         :: center
        type(grid),pointer  :: father,son1,son2,son3,son4
        type(grid),pointer  :: upnei,downnei,leftnei,rightnei
        type(grid),pointer  :: upleftnei,downrightnei,downleftnei,uprightnei
        integer             :: sp
        real( prec )        :: U(7),Ut(7),Un(7)
        real( prec )        :: rot,div
        real( prec )        :: eru,erv,erp,errho 
    end type
       
    type p_array
        type(grid),pointer  :: n
    end type p_array
    
    type(grid),pointer      :: cell(:)
    type(p_array),pointer   :: cell_NS(:),cell_BC(:)
    type(point),allocatable :: Node(:)
    integer                 :: nNode,nElment
    real( prec )            :: Srot,Sdiv
    integer                 :: Ne
    real( prec )            :: Res_max(6),Res_rms(6)
    integer                 :: Neq,nResgrid,nBCgrid
    
    logical                 :: iflag_AMR
    integer                 :: BufferLayerLVL, lvlmax, lvlAMR, lvlblock
    real( prec )            :: Refine_goal_rot, Refine_goal_div, Coarse_goal_rot, Coarse_goal_div, Block_size(4)
    real( prec ),allocatable:: refinesize(:)
     
    real( prec )            :: Ma, gama, cv
    logical                 :: iflag_localdt
    integer                 :: kind_xpbc, kind_xmbc, kind_ypbc, kind_ymbc

    real( prec )            :: tstep
    real( prec )            :: CFL
    integer                 :: it, Nstep_Adapt
    integer                 :: Pstep_Res, Pstep_field, Pstep_Adapt
    
    real( prec )            :: lambda_x, lambda_y
    
    real( prec )            :: tt, t_end
    integer                 :: l_max
    integer                 :: Omp_num_threads  = 1
end module Global_Var
