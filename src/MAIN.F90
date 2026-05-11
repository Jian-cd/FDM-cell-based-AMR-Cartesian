program main      
    use omp_lib
    use mesh
    use solution
    implicit none    
    integer :: total_ompnum
    
    call  set_parameter 
    
    !$  call omp_set_num_threads(Omp_num_threads)
    !$	total_ompnum=OMP_GET_NUM_PROCS() 
    !$	write(*,*)  'OMP total threads number:',total_ompnum
    !$	write(*,*)  'OMP used  threads number:',Omp_num_threads
    !$OMP PARALLEL
    print*,  "                -----run ...-----                          "
    !$OMP END PARALLEL
    
    call  creat_mesh 
    call  idneighbor
    
    call  initial    
    call  solute
        
stop  
end program

    
!======================================================================    
subroutine set_parameter
    use Global_Var
    implicit none
    integer :: i
    
    Omp_num_threads = 4
    Xmin = 0.0; Xmax = 10.0
    Ymin = 0.0; Ymax = 10.0
    h = 0.5_prec
    gama = 1.4
    
    iflag_AMR = 0
    Refine_goal_rot = 1.4; Refine_goal_div = 1.4; Coarse_goal_rot = 0.4; Coarse_goal_div = 0.4
    BufferLayerLVL = 0; 
    lvlAMR = 0 
    lvlmax = 0
    lvlBlock = 0;
    Block_size(1) = 2.5; Block_size(2) = 7.5; Block_size(3) = 2.5; Block_size(4) = 7.5
    
    kind_xpbc = 6; kind_xmbc = 6; kind_ypbc = 6; kind_ymbc = 6
    
    t_end = 10.0; iflag_localdt = 0; tstep = 1e-4; CFL = 0.45
    
    Nstep_Adapt = 0; Pstep_Adapt = 30; Pstep_Res = 100
    
    !==========================================================================================
    Xmin = Xmin - 3.d0*x_factor*h; Xmax = Xmax + 3.d0*x_factor*h + 1.0e-10_prec
    Ymin = Ymin - 3.d0*h; Ymax = Ymax + 3.d0*h + 1.0e-10_prec    
    Block_size(2) = Block_size(2) + 1.0e-10_prec; Block_size(4) = Block_size(4) + 1.0e-10_prec
    
    it = 0.0; tt = 0.0
    
    Ma = 1.0; cv = 1.d0/(gama*(gama-1.d0)*Ma*Ma) 
    
    allocate( refinesize(lvlmax) )
    do i=1,lvlmax
        refinesize(i)=BufferLayerLVL*h/2**(i-1)
    end do         
    
    lambda_x=0;lambda_y=0
    l_max=0
return
end subroutine set_parameter    
