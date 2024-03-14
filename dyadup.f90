    subroutine dyadup(matrix_out,matrix_in,out_size1,out_size2,in_size1,in_size2,type_)
    ! Variables
    integer*4 out_size1,out_size2,in_size1,in_size2
    real*8 matrix_out(out_size1,out_size2),matrix_in(in_size1,in_size2)
    character*3 type_   
    
    integer*4 i
    
    ! Body of dyadup
    matrix_out = 0.0
    
    
    if (type_ == 'row' ) then
        do i = 1,in_size1
            matrix_out(2*i-1,:) = matrix_in(i,:)
        end do                
    else if (type_ == 'col') then
        do i = 1,in_size2
            matrix_out(:,2*i-1) = matrix_in(:,i)
        end do       
        
    end if
    
    
    
    
    
    end subroutine dyadup