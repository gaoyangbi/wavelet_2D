    subroutine wavedec2(c,s,data_ceshi,n,method,data_row,data_col,sum,wavelet_row)
    ! Variables
    implicit none
    integer*4 data_row,data_col
    integer*4 n,sum,row_,col_,wavelet_row
    integer*4 i,c_end,var
    real*8 data_ceshi(data_row,data_col)
    real*8 c(sum)
    integer*4 s(n+2,2)
    character*100 method
    
    real*8 ,allocatable::x0(:,:),x(:,:),h(:,:),v(:,:),d(:,:)
    
    
    


    row_ = data_row
    col_ = data_col
    
    s(n+2,1) = row_
    s(n+2,2) = col_
    
    ! Body of wavedec2
    do i = 1,n
        
        
        allocate(x(floor((row_+wavelet_row-1)/2.0),floor((col_+wavelet_row-1)/2.0)))
        allocate(h(floor((row_+wavelet_row-1)/2.0),floor((col_+wavelet_row-1)/2.0)))
        allocate(v(floor((row_+wavelet_row-1)/2.0),floor((col_+wavelet_row-1)/2.0)))
        allocate(d(floor((row_+wavelet_row-1)/2.0),floor((col_+wavelet_row-1)/2.0)))
        
        x = 0
        h = 0
        v = 0
        d = 0
        
        
        if (i == 1) then
            allocate(x0(row_,col_))
            x0 = data_ceshi
            c_end = sum
        end if

        
        call dwt2(x0,x,h,v,d,c,s,method,sum,row_,col_,n,i,c_end,wavelet_row)
        deallocate(x0)
        
        allocate(x0(floor((row_+wavelet_row-1)/2.0),floor((col_+wavelet_row-1)/2.0)))
        x0 = x
        deallocate(x,h,v,d)
        
        row_ = floor((row_+wavelet_row-1)/2.0)
        col_ = floor((col_+wavelet_row-1)/2.0)
        
        
        

    end do    
    c(:size(x0)) = reshape(x0,(/size(x0)/))
    s(1,1) = size(x0,1)
    s(1,2) = size(x0,2)

    
101 format(A100)    
    end subroutine wavedec2