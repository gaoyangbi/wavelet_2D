    subroutine wextend(mode,y,x0,ysize_1,ysize_2,size_EXT,row_x0,col_x0)
    !   对输入数据矩阵进行反向扩充
    ! Variables
    implicit none
    integer*4 row_x0,col_x0,ysize_1,ysize_2
    integer*4 size_EXT,i,j
    real*8 y(ysize_1,ysize_2)
    real*8 x0(row_x0,col_x0)
    character*100 mode
    
    
    ! Body of wextend
    

!***************************************************行数不变，扩充列数。
    if (trim(mode) == 'addcol') then
        y(:,size_EXT+1:size_EXT+col_x0) = x0
    
        if (size_EXT .ge. col_x0) then
            print *, 'the col of data is not enough!!'
            return
        end if
        
        do j = 1,size_EXT
            y(:,size_EXT+1-j) = x0(:,j)
            y(:,size_EXT+col_x0+j) = x0(:,col_x0+1-j)
        end do
!****************************************************列数不变，扩充行数。   
    else if (trim(mode) == 'addrow') then
        y(size_EXT+1:size_EXT+row_x0,:) = x0
        
        if (size_EXT .ge. row_x0) then
            print *, 'the row of data is not enough!!'
            return
        end if
        
        do j = 1,size_EXT
            y(size_EXT+1-j,:) = x0(j,:)
            y(size_EXT+row_x0+j,:) = x0(row_x0+1-j,:)
        end do

    end if
    
10  format(A2 I2.2 x A14)    

    
    end subroutine wextend