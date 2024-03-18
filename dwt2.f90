    subroutine dwt2(x0,x,h,v,d,c,s,method,sum,row_,col_,n,i,c_end,wavelet_row)
    !  针对单次的小波递进分解计算
    ! Variables
    implicit none
    integer*4 n,sum,row_,col_,c_end,wavelet_row,j,i
    integer*4 size_EXT
    integer*4 first(2),last(2)
    character*100 method,wave_filename,string,mode_tend
    character*100 mode_conv
    real*8 x0(row_,col_)
    real*8 x(floor((row_+wavelet_row-1)/2.0),floor((col_+wavelet_row-1)/2.0))
    real*8 h(floor((row_+wavelet_row-1)/2.0),floor((col_+wavelet_row-1)/2.0))
    real*8 v(floor((row_+wavelet_row-1)/2.0),floor((col_+wavelet_row-1)/2.0))
    real*8 d(floor((row_+wavelet_row-1)/2.0),floor((col_+wavelet_row-1)/2.0))
    real*8 c(sum)
    integer*4 s(n+2,2)
    real*8 ,allocatable::Lo_D(:),Hi_D(:),Lo_R(:),Hi_R(:)                           !  小波基函数 四分量
    real*8 ,allocatable::y(:,:),conv_core(:,:)
    real*8 ,allocatable::z(:,:)
    
    
!***********************************************************读取小波基函数
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
    
    first = (/2,2/)
    
    size_EXT = wavelet_row - 1
    last(1) = size(x0,1) + wavelet_row - 1
    last(2) = size(x0,2) + wavelet_row - 1
    
!*******************************************************对输入数据进行矩阵扩展
    allocate(y(row_ , col_ + 2 * size_EXT))    
    y = 0    
    mode_tend = 'addcol'
    
    call wextend(mode_tend,y,x0,row_,col_ + 2 * size_EXT,size_EXT,row_,col_)
    
    
    allocate(conv_core(1,wavelet_row))
    allocate(z(row_,(col_ + size_EXT)))
    z = 0
    conv_core = 0
    do j = 1,wavelet_row,1
        conv_core(1,j) = Lo_D(j)
    end do
    
    mode_conv = 'valid'    
    
    call conv2(mode_conv,z,y,conv_core,size(z,1),size(z,2),size(y,1),size(y,2),size(conv_core,1),size(conv_core,2))     
    x = 0
    h = 0
    call convdown(x,z,Lo_D,size_EXT,first,last,size(x,1),size(x,2),size(z,1),size(z,2),size(Lo_D))
    call convdown(h,z,Hi_D,size_EXT,first,last,size(h,1),size(h,2),size(z,1),size(z,2),size(Hi_D))
    
    
    do j = 1,wavelet_row,1
        conv_core(1,j) = Hi_D(j)
    end do
    
    call conv2(mode_conv,z,y,conv_core,size(z,1),size(z,2),size(y,1),size(y,2),size(conv_core,1),size(conv_core,2))
    v = 0
    d = 0
    call convdown(v,z,Lo_D,size_EXT,first,last,size(v,1),size(v,2),size(z,1),size(z,2),size(Lo_D))
    call convdown(d,z,Hi_D,size_EXT,first,last,size(d,1),size(d,2),size(z,1),size(z,2),size(Hi_D))
    
   
    c(c_end -size(h)+1 : c_end) = reshape(d,(/size(d)/)) 
    c(c_end -size(h)-size(v)+1 : c_end -size(h)) = reshape(v,(/size(v)/)) 
    c(c_end -size(h)-size(v)-size(d)+1 : c_end -size(h)-size(v)) = reshape(h,(/size(h)/)) 
    c_end = c_end -size(h)-size(v)-size(d)
    s(n+2-i,1) = size(x,1)
    s(n+2-i,2) = size(x,2)
    


    close(100)    
    deallocate(Lo_D,Hi_D,Lo_R,Hi_R,y,z,conv_core)
    
101 format(A100)    
    end subroutine dwt2
    