    subroutine upsconv2(matrix_out,matrix_in,R1_wavelet,R2_wavelet,out_size1,out_size2,in_size1,in_size2,wavelet_row)
    ! Variables
    integer*4 out_size1,out_size2,in_size1,in_size2,wavelet_row
    real*8 R1_wavelet(wavelet_row),R2_wavelet(wavelet_row)
    real*8 matrix_out(out_size1,out_size2),matrix_in(in_size1,in_size2)
    
    real*8 ,allocatable::linshi_matrix1(:,:),linshi_matrix2(:,:)
    real*8 ,allocatable::conv_core(:,:)
    character*100 mode_conv
    integer*4 i
    
    ! Body of upsconv2
    
    allocate(linshi_matrix1((2*in_size1-1),in_size2))    
    linshi_matrix1 = 0
    call dyadup(linshi_matrix1,matrix_in,(2*in_size1-1),in_size2,in_size1,in_size2,'row')
    
    allocate(conv_core(wavelet_row,1))
    allocate(linshi_matrix2((2*in_size1+wavelet_row-2),in_size2))
    mode_conv = 'full_'    
    conv_core = reshape(R1_wavelet,(/size(R1_wavelet),1/))    
    linshi_matrix2 = 0
    call conv2(mode_conv,linshi_matrix2,linshi_matrix1,conv_core,(2*in_size1+wavelet_row-2),in_size2,(2*in_size1-1),in_size2,wavelet_row,1)
    deallocate(linshi_matrix1)
    
    allocate(linshi_matrix1(size(linshi_matrix2,1),size(linshi_matrix2,2)*2-1))    
    linshi_matrix1 = 0
    call dyadup(linshi_matrix1,linshi_matrix2,size(linshi_matrix1,1),size(linshi_matrix1,2),size(linshi_matrix2,1),size(linshi_matrix2,2),'col')
    deallocate(linshi_matrix2)
    
    allocate(linshi_matrix2(size(linshi_matrix1,1),size(linshi_matrix1,2)+size(R2_wavelet)-1))
    conv_core = reshape(R2_wavelet,(/1,size(R2_wavelet)/))    
    linshi_matrix2 = 0
    call conv2(mode_conv,linshi_matrix2,linshi_matrix1,conv_core,size(linshi_matrix2,1),size(linshi_matrix2,2),size(linshi_matrix1,1),size(linshi_matrix1,2),1,wavelet_row)
    deallocate(linshi_matrix1)
    
    allocate(linshi_matrix1(out_size1,out_size2))   
    linshi_matrix1 = 0
    call wkeep2(linshi_matrix1,linshi_matrix2,out_size1,out_size2,size(linshi_matrix2,1),size(linshi_matrix2,2))
    
    
    
    matrix_out  =  linshi_matrix1
    
    deallocate(linshi_matrix1,linshi_matrix2,conv_core)
    
    end subroutine upsconv2