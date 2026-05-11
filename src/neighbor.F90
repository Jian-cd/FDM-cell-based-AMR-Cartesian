module neighbor  
    use Global_Var    
    implicit none
    contains
    
function initial_up_neighbor(c)
    implicit none
    type(grid),pointer::initial_up_neighbor,c
    integer i
  
    if( c%num(1)==n )then  
        initial_up_neighbor=>null()
    else
        i=c%num(1)*m+c%num(2)
        initial_up_neighbor=>cell(i)  
    end if  
return  
end function initial_up_neighbor


function initial_down_neighbor(c)
    implicit none
    type(grid),pointer::initial_down_neighbor,c           
    integer i
  
    if( c%num(1)==1 )then              
        initial_down_neighbor=>null()
    else
        i=c%num(1)*m-2*m+c%num(2)
        initial_down_neighbor=>cell(i)  
    end if  
return  
end function initial_down_neighbor


function initial_left_neighbor(c)
    implicit none
    type(grid),pointer::initial_left_neighbor,c
    integer i
  
    if( c%num(2)==1 )then 
        initial_left_neighbor=>null()
    else
        i=c%num(1)*m-m+c%num(2)-1
        initial_left_neighbor=>cell(i)  
    end if  
return  
end function initial_left_neighbor


function initial_right_neighbor(c)
    implicit none
    type(grid),pointer::initial_right_neighbor,c
    integer i
  
    if( c%num(2)==m )then 
        initial_right_neighbor=>null()
    else
        i=c%num(1)*m-m+c%num(2)+1
        initial_right_neighbor=>cell(i)  
    end if  
return  
end function initial_right_neighbor

recursive function up_neighbor(cell)    
    implicit none    
    type(grid),pointer::cellf,cellfn
    type(grid),pointer::cell,celln,up_neighbor  
  
    celln=>initial_up_neighbor(cell)  
    up_neighbor=>null() 
    cellf=>cell%father
    
    if(.not.associated(cell%father))then    
        up_neighbor=>celln 
        return                                
    else if(cell%sp==1)then                      
        up_neighbor=>cellf%son4
        return
    else if(cell%sp==2)then
        up_neighbor=>cellf%son3
        return  
    end if  
  
    cellfn=>up_neighbor(cellf) 
    if(associated(cellfn))then
        if(.not.associated(cellfn%son1))then    
            up_neighbor=>cellfn 
        else if(cell%sp==4)then
            up_neighbor=>cellfn%son1
        else if(cell%sp==3)then
            up_neighbor=>cellfn%son2
        end if    
    endif 
return 
end function up_neighbor

recursive function down_neighbor(cell)    
    implicit none    
    type(grid),pointer::cellf,cellfn
    type(grid),pointer::cell,celln,down_neighbor
  
    celln=>initial_down_neighbor(cell)  
    down_neighbor=>null()
    cellf=>cell%father
    
    if(.not.associated(cell%father))then  
        down_neighbor=>celln 
        return            
    else if(cell%sp==4)then
        down_neighbor=>cellf%son1
        return
    else if(cell%sp==3)then
        down_neighbor=>cellf%son2
        return  
    end if  
  
    cellfn=>down_neighbor(cellf)  
    if(associated(cellfn))then                
        if(.not.associated(cellfn%son1))then   
            down_neighbor=>cellfn 
        else if(cell%sp==1)then
            down_neighbor=>cellfn%son4
        else if(cell%sp==2)then
            down_neighbor=>cellfn%son3
        end if
    endif       
return 
end function down_neighbor

recursive function left_neighbor(cell)    
    implicit none    
    type(grid),pointer::cellf,cellfn
    type(grid),pointer::cell,celln,left_neighbor
  
    celln=>initial_left_neighbor(cell)
    left_neighbor=>null()
    cellf=>cell%father
    
    if(.not.associated(cell%father))then  
        left_neighbor=>celln 
        return            
    else if(cell%sp==2)then
        left_neighbor=>cellf%son1
        return
    else if(cell%sp==3)then
        left_neighbor=>cellf%son4
        return  
    end if  
  
    cellfn=>left_neighbor(cellf)  
    if(associated(cellfn))then                 
        if(.not.associated(cellfn%son1))then 
            left_neighbor=>cellfn 
        else if(cell%sp==4)then
            left_neighbor=>cellfn%son3
        else if(cell%sp==1)then
            left_neighbor=>cellfn%son2
        end if
    endif     
return 
end function left_neighbor

recursive function right_neighbor(c)    
    implicit none    
    type(grid),pointer::cellf,cellfn
    type(grid),pointer::c,celln,right_neighbor
  
    celln=>initial_right_neighbor(c)
    right_neighbor=>null() 
    cellf=>c%father
    
    if(.not.associated(c%father))then  
        right_neighbor=>celln 
        return            
    else if(c%sp==1)then
        right_neighbor=>cellf%son2
        return
    else if(c%sp==4)then
        right_neighbor=>cellf%son3
        return  
    end if  
  
    cellfn=>right_neighbor(cellf)  
  
    if(associated(cellfn))then                
        if(.not.associated(cellfn%son1))then 
            right_neighbor=>cellfn 
        else if(c%sp==2)then
            right_neighbor=>cellfn%son1
        else if(c%sp==3)then
            right_neighbor=>cellfn%son4
        end if
    endif
return 
end function right_neighbor


recursive function upleft_neighbor(cell)    
    implicit none    
    type(grid),pointer::cellf,cellfn
    type(grid),pointer::cell,celln,upleft_neighbor
  
    celln=>left_neighbor(cell)
    cellfn=>up_neighbor(cell)
    upleft_neighbor=>null()
    cellf=>cell%father
    
    if(.not.associated(cell%father))then  
        if((.not.associated(celln)).or.(.not.associated(cellfn)))then
            upleft_neighbor=>null()
        else
            upleft_neighbor=>up_neighbor(celln)
        endif
        return            
    else if(cell%sp==2)then
        upleft_neighbor=>cellf%son4
        return
    else if(cell%sp==4)then
        if((.not.associated(celln)).or.(.not.associated(cellfn)))then
            upleft_neighbor=>null()
        else
            upleft_neighbor=>up_neighbor(celln)
            do while((upleft_neighbor%lvl<cell%lvl).and.(associated(upleft_neighbor%son2)))        
                upleft_neighbor=>upleft_neighbor%son2
            end do
        endif
        return
    else if(cell%sp==1)then
        if((.not.associated(celln)).or.(celln%lvl<cell%lvl))then
            upleft_neighbor=>null()
        else
            upleft_neighbor=>up_neighbor(celln)
        endif
        return
    else if(cell%sp==3)then
        if((.not.associated(cellfn)).or.(cellfn%lvl<cell%lvl))then
            upleft_neighbor=>null()  
        else
            upleft_neighbor=>left_neighbor(cellfn)            
        endif
        return
    end if  
           
return 
end function upleft_neighbor


recursive function upright_neighbor(cell)    
    implicit none    
    type(grid),pointer::cellf,cellfn
    type(grid),pointer::cell,celln,upright_neighbor
  
    celln=>right_neighbor(cell)
    cellfn=>up_neighbor(cell)
    upright_neighbor=>null()
    cellf=>cell%father
    
    if(.not.associated(cell%father))then 
        if((.not.associated(celln)).or.(.not.associated(cellfn)))then
            upright_neighbor=>null()
        else
            upright_neighbor=>up_neighbor(celln)
        endif
        return            
    else if(cell%sp==1)then
        upright_neighbor=>cellf%son3
        return
    else if(cell%sp==3)then
        if((.not.associated(celln)).or.(.not.associated(cellfn)))then
            upright_neighbor=>null()
        else
            upright_neighbor=>up_neighbor(celln)
            do while((upright_neighbor%lvl<cell%lvl).and.(associated(upright_neighbor%son1)))
                upright_neighbor=>upright_neighbor%son1
            end do
        endif
        return
    else if(cell%sp==2)then
        if((.not.associated(celln)).or.(celln%lvl<cell%lvl))then
            upright_neighbor=>null()
        else
            upright_neighbor=>up_neighbor(celln)
        endif
        return
    else if(cell%sp==4)then
        if((.not.associated(cellfn)).or.(cellfn%lvl<cell%lvl))then
            upright_neighbor=>null()
        else
            upright_neighbor=>right_neighbor(cellfn)
        endif
        return
    end if  
           
return 
end function upright_neighbor


recursive function downright_neighbor(cell)    
    implicit none    
    type(grid),pointer::cellf,cellfn
    type(grid),pointer::cell,celln,downright_neighbor
  
    celln=>right_neighbor(cell)
    cellfn=>down_neighbor(cell)
    downright_neighbor=>null()
    cellf=>cell%father
    
    if(.not.associated(cell%father))then 
        if((.not.associated(celln)).or.(.not.associated(cellfn)))then
            downright_neighbor=>null()
        else
            downright_neighbor=>down_neighbor(celln)
        endif
        return              
    else if(cell%sp==4)then
        downright_neighbor=>cellf%son2
        return
    else if(cell%sp==2)then
        if((.not.associated(celln)).or.(.not.associated(cellfn)))then
            downright_neighbor=>null()
        else
            downright_neighbor=>down_neighbor(celln)
            do while((downright_neighbor%lvl<cell%lvl).and.(associated(downright_neighbor%son4)))
               downright_neighbor=>downright_neighbor%son4
            end do
        endif
        return
    else if(cell%sp==3)then
        if((.not.associated(celln)).or.(celln%lvl<cell%lvl))then
            downright_neighbor=>null()
        else
            downright_neighbor=>down_neighbor(celln)
        endif
        return
    else if(cell%sp==1)then
        if((.not.associated(cellfn)).or.(cellfn%lvl<cell%lvl))then
            downright_neighbor=>null()
        else
            downright_neighbor=>right_neighbor(cellfn)
        endif
        return
    end if  
           
return 
end function downright_neighbor


recursive function downleft_neighbor(cell)    
    implicit none    
    type(grid),pointer::cellf,cellfn
    type(grid),pointer::cell,celln,downleft_neighbor
  
    celln=>left_neighbor(cell)
    cellfn=>down_neighbor(cell)
    downleft_neighbor=>null()
    cellf=>cell%father
    
    if(.not.associated(cell%father))then 
        if((.not.associated(celln)).or.(.not.associated(cellfn)))then
            downleft_neighbor=>null()
        else
            downleft_neighbor=>down_neighbor(celln) 
        endif
        return            
    else if(cell%sp==3)then
        downleft_neighbor=>cellf%son1
        return
    else if(cell%sp==1)then
         if((.not.associated(celln)).or.(.not.associated(cellfn)))then
            downleft_neighbor=>null()
        else
            downleft_neighbor=>down_neighbor(celln)
            do while((downleft_neighbor%lvl<cell%lvl).and.(associated(downleft_neighbor%son3)))
                 downleft_neighbor=>downleft_neighbor%son3
            end do
        endif
        return
    else if(cell%sp==4)then
        if((.not.associated(celln)).or.(celln%lvl<cell%lvl))then
            downleft_neighbor=>null()
        else
            downleft_neighbor=>down_neighbor(celln)
        endif
        return
    else if(cell%sp==2)then
        if((.not.associated(cellfn)).or.(cellfn%lvl<cell%lvl))then
            downleft_neighbor=>null()
        else
            downleft_neighbor=>left_neighbor(cellfn)
        endif
        return
    end if  
           
return 
end function downleft_neighbor


subroutine idneighbor
    implicit none
    type(grid),pointer::c
    integer i
    
    do i=1,total
        c=>cell(i)   
        call idnei(c)
    end do 
return
end subroutine idneighbor


recursive subroutine idnei(c)
implicit none
type(grid),pointer::c,ct1,ct2,ct3,ct4,ct5,ct6,ct7,ct8
type(grid),pointer::ct9,ct10,ct11,ct12

    ct1=>up_neighbor(c)
    ct2=>down_neighbor(c)
    ct3=>right_neighbor(c)
    ct4=>left_neighbor(c)
    if(c%cross>=0)then
        ct5=>upleft_neighbor(c)  
        ct6=>upright_neighbor(c) 
        ct7=>downleft_neighbor(c)
        ct8=>downright_neighbor(c)
    else
        ct5=>null()  
        ct6=>null() 
        ct7=>null()
        ct8=>null()
    endif
  
    c%upnei=>ct1
    c%downnei=>ct2
    c%rightnei=>ct3
    c%leftnei=>ct4
    c%upleftnei=>ct5
    c%uprightnei=>ct6
    c%downleftnei=>ct7
    c%downrightnei=>ct8

    if(associated(c%son1))then
        call idnei(c%son1)
        call idnei(c%son2)
        call idnei(c%son3)
        call idnei(c%son4)
    end if
return
end subroutine idnei

end module neighbor
