    subroutine wrcoef2(type_,x,c,s,method,degree,xsize_1,xsize_2,c_size,s_size,wavelet_row)
    ! Variables
    integer*4 degree,c_size,s_size,xsize_1,xsize_2
    integer*4 rmax,nmax,nmin,wavelet_row,imin
    character*1 type_
    character*100 method,wave_filename
    
    real*8 ,allocatable::Lo_D(:),Hi_D(:),Lo_R(:),Hi_R(:)    
    real*8 c(c_size),x(xsize_1,xsize_2)
    integer*4 s(s_size,2)
    real*8 ,allocatable::linshi_matrix1(:,:),linshi_matrix2(:,:)
    real*8 ,allocatable::h(:,:),v(:,:),d(:,:)
    real*8 ,allocatable::F1(:),F2(:)    
    integer*4 p
    
    ! Body of 
    rmax = size(s,1)
    nmax = rmax-2
    if (type_ == 'a') then 
        nmin = 0
    else 
        nmin = 1
    end if
    
!--------------------------------    
    wave_filename = trim(method)//'.txt'
    
    open(unit = 100, file = wave_filename , status='old', position='rewind')
    
    allocate(Lo_D(wavelet_row))
    allocate(Hi_D(wavelet_row))
    allocate(Lo_R(wavelet_row))
    allocate(Hi_R(wavelet_row))

    
    do j = 1,wavelet_row,1       
        read(100, *) Lo_D(j),Hi_D(j),Lo_R(j),Hi_R(j)
    end do     
    
    close(100)
!---------------------------------------    
    allocate(linshi_matrix1(s(s_size-degree,1),s(s_size-degree,2)))
    allocate(h(s(s_size-degree,1),s(s_size-degree,2)))
    allocate(v(s(s_size-degree,1),s(s_size-degree,2)))
    allocate(d(s(s_size-degree,1),s(s_size-degree,2)))
    allocate(F1(wavelet_row))
    allocate(F2(wavelet_row))    
    
    linshi_matrix1 = 0
    h              = 0
    v              = 0
    d              = 0    
    
    call detcoef2(h,v,d,c,s,s(s_size-degree,1),s(s_size-degree,2),c_size,s_size,degree)    
    
    if (type_ == 'a') then
        call appcoef2(linshi_matrix1,c,s,Lo_R,Hi_R,s(s_size-degree,1),s(s_size-degree,2),c_size,s_size,wavelet_row,degree)
        F1  =  Lo_R
        F2  =  Lo_R
    else if (type_ == 'h') then
        linshi_matrix1  =  h
        F1  =  Hi_R
        F2  =  Lo_R        
    else if (type_ == 'v') then
        linshi_matrix1  =  v
        F1  =  Lo_R
        F2  =  Hi_R
    else if (type_ == 'd') then
        linshi_matrix1  =  d
        F1  =  Hi_R
        F2  =  Hi_R
    end if
    
    imin  =  rmax - degree
    allocate(linshi_matrix2(s(imin+1,1),s(imin+1,2)))
    linshi_matrix2 = 0
    call upsconv2(linshi_matrix2,linshi_matrix1,F1,F2,s(imin+1,1),s(imin+1,2),size(linshi_matrix1,1),size(linshi_matrix1,2),wavelet_row)
    
    do p = 2,degree
        deallocate(linshi_matrix1)
        allocate(linshi_matrix1(s(imin+p,1),s(imin+p,2)))
        linshi_matrix1 = 0
        call upsconv2(linshi_matrix1,linshi_matrix2,Lo_R ,Lo_R,size(linshi_matrix1,1),size(linshi_matrix1,2),size(linshi_matrix2,1),size(linshi_matrix2,2),wavelet_row)
        
        deallocate(linshi_matrix2)
        allocate(linshi_matrix2(s(imin+p,1),s(imin+p,2)))
        
        linshi_matrix2  =  linshi_matrix1       
    end do
       
    
    x = linshi_matrix2
    
    deallocate(Lo_D,Hi_D,Lo_R,Hi_R,linshi_matrix1,h,v,d)    
    deallocate(linshi_matrix2)
    deallocate(F1,F2)
    
    
    
    end subroutine 