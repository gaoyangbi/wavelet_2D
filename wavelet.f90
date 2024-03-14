!  wavelet.f90 
!
!  FUNCTIONS:
!  wavelet - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: wavelet
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program wavelet

    implicit none
    ! Variables
    
    integer*4 data_row,data_col,n
    integer*4 i,m,kold,k
    character*10000 cfile,method,ofile
    character*200, allocatable::in_file(:)

    
    open(unit = 90 ,file = 'wavelet.txt', status='old')
    read(90,'(a)') cfile    !
    read(90,'(a)') ofile
    read(90,*) data_row
    read(90,*) data_col
    read(90,*) method
    read(90,*) n    
    close(90)
    
    m = 0
    do i = 1,len(trim(cfile))
        if(cfile(i:i) == ',') then
            m = m + 1
        end if
    end do
    
    m = m + 1
    
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
    
    do i = 1,m        
        call wavelet_main(in_file(i),ofile,data_row,data_col,method,n,i)
    end do

    
    deallocate(in_file)

    end program wavelet
    
    
    
    
    subroutine wavelet_main(ifile,ofile,data_row,data_col,method,n,num)
    ! Variables
    
    real*8 ,allocatable::data_(:,:)
    real*8 ,allocatable::c(:)
    integer*4 ,allocatable::s(:,:)
    real*8 ,allocatable::a(:,:),h(:,:),v(:,:),d(:,:)
    
    integer*4 data_row,data_col,n,sum,row_,col_,wavelet_row,var
    integer*4 i,j,m,num
    real*8 ,allocatable::latitude(:,:),longitude(:,:)
    character*100 method,wave_filename,string
    character*200 ifile
    character*200 file_in
    character*100 degree
    character*100 num_string    
    character*10000 ofile
        
    allocate(data_(data_row,data_col))
    allocate(longitude(data_row,data_col))
    allocate(latitude(data_row,data_col))
    
    open(unit = 95 ,file = trim(ifile), status='old')
    
    do i = data_row,1,-1
        do j = 1,data_col
            read(95,*) longitude(i,j),latitude(i,j),data_(i,j)            
        end do
    end do  
    close(95)   
    
    sum = 0

    wave_filename = trim(method)//'.txt'
    open(unit = 100, file = wave_filename , status='old', position='rewind')
    wavelet_row = 0
    
    
    do while(1)   !确认文件中数据内容行数
        read(unit = 100,fmt = 101, iostat = var) string
        if (var < 0) then
            exit
        end if      
            
        wavelet_row = wavelet_row + 1
            
    end do  
    
    close(100)      
    
    !计算矩阵C的长度
    row_ = data_row
    col_ = data_col
    
    do i = 1,n-1
        sum = sum + floor((row_+wavelet_row-1)/2.0) * floor((col_+wavelet_row-1)/2.0) * 3
        row_ = floor((row_+wavelet_row-1)/2.0)
        col_ = floor((col_+wavelet_row-1)/2.0)
    end do
    sum = sum + floor((row_+wavelet_row-1)/2.0) * floor((col_+wavelet_row-1)/2.0) * 4
     
    !动态矩阵初始化
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
    
    ! Body of wavelet
    
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
    

