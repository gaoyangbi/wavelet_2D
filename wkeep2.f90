    subroutine wkeep2(y,x,y_size1,y_size2,x_size1,x_size2)
    ! Variables
    integer*4 y_size1,y_size2,x_size1,x_size2
    real*8 y(y_size1,y_size2),x(x_size1,x_size2)
    
    integer*4 sx(1,2),siz(1,2),first(1,2),last(1,2)
    integer*4 k
    real*8 d(1,2)
    ! Body of wkeep2
    
    first    = 0
    last     = 0
    sx(1,1)  = x_size1
    sx(1,2)  = x_size2
    siz(1,1) = y_size1
    siz(1,2) = y_size2

    d        = (sx-siz)/2.0
    
    do k = 1,2
        first(1,k) = 1+floor(d(1,k));
        last(1,k)  = sx(1,k)-ceiling(d(1,k));
    end do
    
    
    y = x(first(1,1):last(1,1),first(1,2):last(1,2))
    
    end subroutine wkeep2