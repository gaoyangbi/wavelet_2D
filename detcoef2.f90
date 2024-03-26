    subroutine detcoef2(h,v,d,c,s,hvd_size1,hvd_size2,c_size,s_size,degree)
    ! Variables
    ! 'l' = 'all' 
    !  针对细节分量 hvd 进行重构
    integer*4 hvd_size1,hvd_size2,c_size,s_size,degree
    real*8 c(c_size)
    integer*4 s(s_size,2)    
    real*8 h(hvd_size1,hvd_size2),v(hvd_size1,hvd_size2),d(hvd_size1,hvd_size2)   
    
    integer*4 nmax,k,first,add,last
    ! Body of detcoef2
    
    nmax = size(s,1)-2;
    k    = size(s,1)-degree;
    first= s(1,1)*s(1,2) + 3*sum(s(2:k-1,1)*s(2:k-1,2))+1;
    add  = s(k,1)*s(k,2)
    
        
    last  = first+add-1;
    h     = reshape(c(first:last),(/s(k,:)/))
        
    first = first+add
    last  = first+add-1;
    v     = reshape(c(first:last),(/s(k,:)/))
        
    first = first+add
    last  = first+add-1;
    d     = reshape(c(first:last),(/s(k,:)/))

  
    end subroutine detcoef2