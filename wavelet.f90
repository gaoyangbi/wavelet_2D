!  wavelet.f90 
!
!  FUNCTIONS:
!  wavelet - 进行二维小波滤波计算
!

!****************************************************************************
!
!  PROGRAM: wavelet
!
!  PURPOSE:  函数使用范例程序
!
!****************************************************************************

    program wavelet

    implicit none
    ! Variables 
    
    integer*4 data_row,data_col,n                  ! data_row: 数据行数    data_col: 数据列数  n：分解阶数
    integer*4 i,m,kold,k                           ! i,m,k,kold : 循环变量      
    character*10000 cfile,method,ofile             ! cfile : 所有输入文件路径，以英文`，`分隔    ofile : 输出文件夹的路径名   method : 选择小波滤波基
    character*200, allocatable::in_file(:)         ! in_file : 将每个输入文件进行拆分，保存至该数组

    
    open(unit = 90 ,file = 'wavelet.txt', status='old')
    read(90,'(a)') cfile    !  读取 cfile
    read(90,'(a)') ofile    !  读取 ofile
    read(90,*) data_row     
    read(90,*) data_col
    read(90,*) method
    read(90,*) n    
    close(90)
    

!******************************************** 此处确定输入文件的个数 m
    m = 0

    do i = 1,len(trim(cfile))
        if(cfile(i:i) == ',') then
            m = m + 1
        end if
    end do
    
    m = m + 1        ! 
    
!********************************************对cfile进行拆分，将每一个文件路径名保存至in_file
    allocate(in_file(m))
    m = 0
    kold = 0
    
    do k = 1,len(trim(cfile))
        if(cfile(k:k) == ',') then
            m = m + 1
            in_file(m) = cfile(kold+1:k-1)
            kold = k
        end if
    end do
    m = m + 1
    in_file(m) = cfile(kold+1:k-1)   
    
!*********************************************对每个输入文件进行小波分解
    do i = 1,m        
        call wavelet_main(in_file(i),ofile,data_row,data_col,method,n,i)
    end do

    
    deallocate(in_file)

    end program wavelet
    
    
    
    
    subroutine wavelet_main(ifile,ofile,data_row,data_col,method,n,num)
    ! Variables
    ! ifile    :  输入文件路径 
    ! ofile    :  输出文件夹路径
    ! data_row :  数据行数
    ! data_col :  数据列数
    ! method   :  滤波方法
    ! n        :  滤波阶数
    ! num      :  文件的编号数
    
    real*8 ,allocatable::data_(:,:)                           !  data_ : 用于保存读取的数据
    real*8 ,allocatable::c(:)                                 !  c : 保存中间数据的过程矩阵
    integer*4 ,allocatable::s(:,:)                            !  s : 保存输出矩阵维度的过程矩阵
    real*8 ,allocatable::a(:,:),h(:,:),v(:,:),d(:,:)          !  a,h,v,d :  保存分解结果，对应各自的小波成分
    
    integer*4 data_row,data_col,n,sum,row_,col_,wavelet_row,var  !  sum, row_, col_ : 中间变量       wavelet_row : 小波基函数的数据个数  var : 文件读取判断变量
    integer*4 i,j,m,num                                          !  i,j,m :  循环变量
    real*8 ,allocatable::latitude(:,:),longitude(:,:)            !  latitude, longitude :  数据经纬度
    character*100 method,wave_filename,string                    !  wave_filename :  小波基文件名   string : 临时字符串中间变量
    character*200 ifile
    character*200 file_in                                        !  file_in :  输出文件名
    character*100 degree                                         !  degree :  分解阶数
    character*100 num_string                                     !  num_string :  输出文件编号数
    character*10000 ofile
        
    allocate(data_(data_row,data_col))
    allocate(longitude(data_row,data_col))
    allocate(latitude(data_row,data_col))
    
!****************************************************读取数据文件
    open(unit = 95 ,file = trim(ifile), status='old')    
    do i = data_row,1,-1
        do j = 1,data_col
            read(95,*) longitude(i,j),latitude(i,j),data_(i,j)            
        end do
    end do  
    close(95)   

!****************************************************确定所选小波基函数的数据个数
    sum = 0

    wave_filename = trim(method)//'.txt'
    open(unit = 100, file = wave_filename , status='old', position='rewind')
    wavelet_row = 0    
    
    do while(.true.)   
        read(unit = 100,fmt = 101, iostat = var) string
        if (var < 0) then
            exit
        end if      
            
        wavelet_row = wavelet_row + 1            
    end do  
    
    close(100)     

!***************************************************** 确实过程矩阵c的维度，并进行小波分解与重构
    row_ = data_row
    col_ = data_col
    
    do i = 1,n-1
        sum = sum + floor((row_+wavelet_row-1)/2.0) * floor((col_+wavelet_row-1)/2.0) * 3
        row_ = floor((row_+wavelet_row-1)/2.0)
        col_ = floor((col_+wavelet_row-1)/2.0)
    end do
    sum = sum + floor((row_+wavelet_row-1)/2.0) * floor((col_+wavelet_row-1)/2.0) * 4

    allocate(c(sum))
    allocate(s(n+2,2))
    allocate(a(data_row,data_col))
    allocate(h(data_row,data_col))
    allocate(v(data_row,data_col))
    allocate(d(data_row,data_col))

    call wavedec2(c,s,data_,n,method,data_row,data_col,sum,wavelet_row)        

    
    call wrcoef2('a',a,c,s,method,n,data_row,data_col,sum,n+2,wavelet_row)
    call wrcoef2('h',h,c,s,method,n,data_row,data_col,sum,n+2,wavelet_row)
    call wrcoef2('v',v,c,s,method,n,data_row,data_col,sum,n+2,wavelet_row)
    call wrcoef2('d',d,c,s,method,n,data_row,data_col,sum,n+2,wavelet_row)
    
!*******************************************************将分解结果进行文件输出
    do m = 1,n
        write(degree,*) m
        write(num_string,*) num
        file_in = trim(adjustl(ofile))//'out_a'//trim(adjustl(degree))//'_'//trim(adjustl(num_string))//'.txt'
        open(unit = 200 , file = file_in , status = 'unknown')
        
        do i = data_row,1,-1
            do j = 1,data_col
                write(200,200) longitude(i,j),latitude(i,j),a(i,j)            
            end do
        end do                  
        close(200)
        
        file_in = trim(adjustl(ofile))//'out_h'//trim(adjustl(degree))//'_'//trim(adjustl(num_string))//'.txt'
        open(unit = 200 , file = file_in , status = 'unknown')
        
        do i = data_row,1,-1
            do j = 1,data_col
                write(200,200) longitude(i,j),latitude(i,j),h(i,j)            
            end do
        end do                  
        close(200)
        
        file_in = trim(adjustl(ofile))//'out_v'//trim(adjustl(degree))//'_'//trim(adjustl(num_string))//'.txt'
        open(unit = 200 , file = file_in , status = 'unknown')
        
        do i = data_row,1,-1
            do j = 1,data_col
                write(200,200) longitude(i,j),latitude(i,j),v(i,j)            
            end do
        end do                  
        close(200)
        
        file_in = trim(adjustl(ofile))//'out_d'//trim(adjustl(degree))//'_'//trim(adjustl(num_string))//'.txt'
        open(unit = 200 , file = file_in , status = 'unknown')
        
        do i = data_row,1,-1
            do j = 1,data_col
                write(200,200) longitude(i,j),latitude(i,j),d(i,j)            
            end do
        end do                  
        close(200)   
        
        
    end do 
    deallocate(c,s,a,h,v,d)
200 format(f12.6 x f12.6 x f12.6) 
101 format(A100)    
    
    end subroutine wavelet_main
    

