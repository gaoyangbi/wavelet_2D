    subroutine convdown(out_,z,F,size_EXT,first,last,out_size_1,out_size_2,zsize_1,zsize_2,Fsize_)
    !   进行卷积解码的计算
    ! Variables
    integer*4 size_EXT,linshi_size1,linshi_size2
    integer*4 out_size_1,out_size_2,zsize_1,zsize_2,Fsize_
    integer*4 first(2),last(2),j
    real*8 out_(out_size_1,out_size_2),z(zsize_1,zsize_2),F(Fsize_)
    real*8 ,allocatable::linshi_matrix1(:,:),linshi_matrix2(:,:),conv_core(:,:)
    character*100 mode_tend,mode_conv
    
    
    ! Body of convdown
    
    allocate(linshi_matrix1(size(z(:,first(2):last(2):2),1),size(z(:,first(2):last(2):2),2)))
    allocate(conv_core(1,Fsize_))
    linshi_matrix1 = 0
    conv_core = 0
    
    do j = 1,Fsize_,1
        conv_core(1,j) = F(j)
    end do
    
    linshi_matrix1 = z(:,first(2):last(2):2)
    
    allocate(linshi_matrix2(size(linshi_matrix1,1)+2*size_EXT,size(linshi_matrix1,2)))
    linshi_matrix2 = 0
    mode_tend = 'addrow' 
    call wextend(mode_tend,linshi_matrix2,linshi_matrix1,size(linshi_matrix2,1),size(linshi_matrix2,2),size_EXT,size(linshi_matrix1,1),size(linshi_matrix1,2))    
    
    deallocate(linshi_matrix1)
    
    linshi_size1 = size(transpose(linshi_matrix2),1) 
    linshi_size2 = size(transpose(linshi_matrix2),2) - size(conv_core,2) + 1    
    allocate(linshi_matrix1(linshi_size1,linshi_size2))
    linshi_matrix1 = 0
    
    mode_conv = 'valid'
    call conv2(mode_conv,linshi_matrix1,transpose(linshi_matrix2),conv_core,size(linshi_matrix1,1),size(linshi_matrix1,2),&
        size(transpose(linshi_matrix2),1),size(transpose(linshi_matrix2),2),&
        size(conv_core,1),size(conv_core,2))
    
    
    deallocate(linshi_matrix2)
    
    linshi_size1 = size(transpose(linshi_matrix1),1) 
    linshi_size2 = size(transpose(linshi_matrix1),2)
    allocate(linshi_matrix2(linshi_size1,linshi_size2))
    linshi_matrix2 = transpose(linshi_matrix1)
    
    
    
    deallocate(linshi_matrix1)   

    
    
    linshi_size1 = size(linshi_matrix2(first(1):last(1):2,:),1) 
    linshi_size2 = size(linshi_matrix2,2)
    allocate(linshi_matrix1(linshi_size1,linshi_size2))
    linshi_matrix1 = linshi_matrix2(first(1):last(1):2,:)
    
    
    out_ = linshi_matrix1
    
    
    deallocate(linshi_matrix1,linshi_matrix2,conv_core)
    
    end subroutine convdown