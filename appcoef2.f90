    subroutine appcoef2(out_matrix,c,s,Lo_R,Hi_R,matrix_size1,matrix_size2,c_size,s_size,wavelet_row,degree)
    ! Variables
    integer*4 matrix_size1,matrix_size2,c_size,s_size,wavelet_row,degree
    real*8 out_matrix(matrix_size1,matrix_size2),c(c_size)
    integer*4 s(s_size,2)
    real*8 Lo_R(wavelet_row),Hi_R(wavelet_row)
    
    real*8 ,allocatable::h(:,:),v(:,:),d(:,:)
    real*8 ,allocatable::linshi_matrix_in(:,:),linshi_matrix_out(:,:)
    
    
    integer*4 rmax,nmax,nl,nc,rm,i
    real*8 a(s(1,1),s(1,2))    
    
    ! Body of appcoef2
    rmax = size(s,1)
    nmax = rmax - 2
    
    nl = s(1,1);
    nc = s(1,2);
    
    a = 0.0
    a = reshape(c(1:nl*nc),(/nl,nc/))
    
    rm = rmax+1
    
    if (nmax < degree+1) then
        out_matrix = a
    end if
      
    
    do i = nmax,(degree+1),-1
        
        allocate(h(s(s_size-i,1),s(s_size-i,2)))
        allocate(v(s(s_size-i,1),s(s_size-i,2)))
        allocate(d(s(s_size-i,1),s(s_size-i,2)))
        allocate(linshi_matrix_out(s(rm-i,1),s(rm-i,2)))
        h                 = 0
        v                 = 0
        d                 = 0
        linshi_matrix_out = 0
        
        call detcoef2(h,v,d,c,s,s(s_size-i,1),s(s_size-i,2),c_size,s_size,i)
        
        if (i == nmax) then
            allocate(linshi_matrix_in(1,1))
            linshi_matrix_in = 0
            call idwt2(linshi_matrix_out,a,h,v,d,Lo_R,Hi_R,s(rm-i,1),s(rm-i,2),size(h,1),size(h,2),wavelet_row)
            
        else 
            call idwt2(linshi_matrix_out,linshi_matrix_in,h,v,d,Lo_R,Hi_R,s(rm-i,1),s(rm-i,2),size(h,1),size(h,2),wavelet_row)
            
        end if   
        
        deallocate(linshi_matrix_in)        
        
        allocate(linshi_matrix_in(s(rm-i,1),s(rm-i,2)))
        linshi_matrix_in  =  linshi_matrix_out
        
        if ( i-1 < (degree+1)) then
            out_matrix = linshi_matrix_out
        end if
               
        deallocate(h,v,d,linshi_matrix_out)
    end do
    
    
    end subroutine appcoef2