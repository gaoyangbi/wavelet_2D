    subroutine idwt2(matrix_out,matrix_in,h,v,d,Lo_R,Hi_R,out_size1,out_size2,h_size1,h_size2,wavelet_row)
    ! Variables
    integer*4 out_size1,out_size2,h_size1,h_size2,wavelet_row
    real*8 matrix_out(out_size1,out_size2)
    real*8 matrix_in(h_size1,h_size2),h(h_size1,h_size2),v(h_size1,h_size2),d(h_size1,h_size2)
    real*8 Lo_R(wavelet_row),Hi_R(wavelet_row)
    
    real*8 a_out(out_size1,out_size2),h_out(out_size1,out_size2),v_out(out_size1,out_size2),d_out(out_size1,out_size2)

    
    ! Body of idwt2
    

    a_out = 0
    h_out = 0
    v_out = 0
    d_out = 0
    
    call upsconv2(a_out,matrix_in,Lo_R,Lo_R,out_size1,out_size2,h_size1,h_size2,wavelet_row)
    call upsconv2(h_out,h,Hi_R,Lo_R,out_size1,out_size2,h_size1,h_size2,wavelet_row)
    call upsconv2(v_out,v,Lo_R,Hi_R,out_size1,out_size2,h_size1,h_size2,wavelet_row)
    call upsconv2(d_out,d,Hi_R,Hi_R,out_size1,out_size2,h_size1,h_size2,wavelet_row)
    
    
    matrix_out = a_out + h_out + v_out + d_out
    
    end subroutine idwt2