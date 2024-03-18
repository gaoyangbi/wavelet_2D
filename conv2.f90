    subroutine conv2(mode_conv,z,y,conv_core,zsize_1,zsize_2,ysize_1,ysize_2,conv_core_size_1,conv_core_size_2)
    ! 进行二维卷积运算
    implicit none
    integer*4 ysize_1,ysize_2,conv_core_size_1,conv_core_size_2,zsize_1,zsize_2
    integer*4 i,j,m,n
    real*8 y(ysize_1,ysize_2)
    real*8 conv_core(conv_core_size_1,conv_core_size_2),conv_core180(conv_core_size_1,conv_core_size_2)
    character*100 mode_conv
    real*8 z(zsize_1,zsize_2)    
    real*8 ,allocatable::y_(:,:)
    
    
    conv_core180 = 0
    z            = 0
    
!********************************************************计算卷积核函数
    do m = 1,conv_core_size_1,1
        do n = 1,conv_core_size_2,1
            conv_core180(conv_core_size_1+1-m,conv_core_size_2+1-n) = conv_core(m,n)                         
        end do                                          
    end do 
    
!*********************************************************针对两种卷积方式进行卷积计算
    if (mode_conv == 'valid') then
         do i = 1,zsize_1,1             
             do j = 1,zsize_2,1

                 do m = 1,conv_core_size_1,1
                     do n = 1,conv_core_size_2,1
                         z(i,j) = z(i,j) + y(i+m-1,j+n-1) * conv_core180(m,n)                         
                     end do                                          
                 end do                                  
             end do             
         end do
    else if(mode_conv == 'full_')then
        allocate(y_(ysize_1+2*(conv_core_size_1-1),ysize_2+2*(conv_core_size_2-1)))
        y_ = 0
        y_(conv_core_size_1:(conv_core_size_1+ysize_1-1),conv_core_size_2:(conv_core_size_2+ysize_2-1)) = y 
        
         do i = 1,zsize_1,1             
             do j = 1,zsize_2,1
                 
                 do m = 1,conv_core_size_1,1
                     do n = 1,conv_core_size_2,1
                         z(i,j) = z(i,j) + y_(i+m-1,j+n-1) * conv_core180(m,n)                         
                     end do                                          
                 end do                                  
             end do             
         end do
        
         deallocate(y_)
        
    end if
     
    
    
    end subroutine conv2 