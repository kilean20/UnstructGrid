!----------------------------------------------------------------
! (c) Copyright, 2017 by the Regents of the University of California.
! UnstructGridClass : 
! Version: 
! Author: 
! Description: 
! Comments:
!----------------------------------------------------------------
program preprocess
  implicit none
  integer :: num_args, ix, fileID, nInit_x, nInit_y, nInit_z
  double precision :: xDim,yDim,zDim
  character(len=25) :: fname, str_fileID, str_nInit_x, str_nInit_y, str_nInit_z
  type UnstructGrid
    integer :: nOutNode, nEle, nInNode(3), nInitNode(3)
    integer, allocatable :: Ele(:,:), InitEle(:,:,:), Neigh(:,:)
    double precision, allocatable :: InNode(:,:,:,:), OutNode(:,:)
    double precision :: dNode(3), dInitNode(3)
  end type UnstructGrid
  type(UnstructGrid) :: this

  num_args = command_argument_count()
  if (num_args == 2) then
    call get_command_argument(1,fname)
    call get_command_argument(2,str_fileID)
  elseif (num_args == 5) then
    call get_command_argument(1,fname)
    call get_command_argument(2,str_fileID)  
    call get_command_argument(3,str_nInit_x)
    call get_command_argument(4,str_nInit_y)
    call get_command_argument(5,str_nInit_z)
  else 
    print *, 'usage : ./preprocssor filename fileID [nInit_x, nInit_y, nInit_z]'
    call exit()
  endif
  

  read(str_fileID,*) fileID
  
  call construct_UnstructGrid(this,fname)
  call cal_max_dimension(this,xDim,yDim,zDim)
  call assume_nInit(this,xDim,yDim,zDim)
  
  if(num_args == 5) then
    read(str_nInit_x,*) nInit_x
    read(str_nInit_y,*) nInit_y
    read(str_nInit_z,*) nInit_z
    if (this%nInitNode(1) > nInit_x) then
      this%nInitNode(1) = nInit_x
    else
      print *, 'nInit_x is larger then maximum default. Setting to default value, nInit_x =', this%nInitNode(1)
    endif
    if (this%nInitNode(2) > nInit_y) then
      this%nInitNode(2) = nInit_y
    else
      print *, 'nInit_y is larger then maximum default. Setting to default value, nInit_y =', this%nInitNode(2)
    endif
    if (this%nInitNode(3) > nInit_z) then
      this%nInitNode(3) = nInit_z
    else
      print *, 'nInit_z is larger then maximum default. Setting to default value, nInit_z =', this%nInitNode(3)
    endif        
  endif
    
  this%dInitNode(1) = xDim/DBLE(this%nInitNode(1))+EPSILON(1.0)
  this%dInitNode(2) = yDim/DBLE(this%nInitNode(2))+EPSILON(1.0)
  this%dInitNode(3) = zDim/DBLE(this%nInitNode(3))+EPSILON(1.0)
  
  call locateInitTet(this)
  call print_UnstructGrid(this,fileID)
  
  contains 
  
  subroutine cal_max_dimension(this,xDim,yDim,zDim)
    implicit none
    type(UnstructGrid), intent(in) :: this
    double precision, intent(out) :: xDim,yDim,zDim
    integer :: i,j,k
    xDim = MAXVAL(abs(this%OutNode(:,1)))
    yDim = MAXVAL(abs(this%OutNode(:,2)))
    zDim = MAXVAL(abs(this%OutNode(:,3)))
  end subroutine cal_max_dimension
  
  subroutine assume_nInit(this,xDim,yDim,zDim)
    implicit none
    type(UnstructGrid), intent(inout) :: this
    double precision, intent(in) :: xDim,yDim,zDim
    integer :: i,j,k
    double precision :: jNode(3), kNode(3), xTetMax, yTetMax, zTetMax
    xTetMax = 0
    yTetMax = 0
    zTetMax = 0
    do i=1,this%nEle
      do j=1,4
        jNode = this%OutNode(this%Ele(i,j),1:3)
        do k=j+1,4
          kNode = this%OutNode(this%Ele(i,k),1:3)
          if (abs(jNode(1)-kNode(1)) > xTetMax) xTetMax = abs(jNode(1)-kNode(1))
          if (abs(jNode(2)-kNode(2)) > yTetMax) yTetMax = abs(jNode(2)-kNode(2))
          if (abs(jNode(3)-kNode(3)) > zTetMax) zTetMax = abs(jNode(3)-kNode(3))
        enddo
      enddo
    enddo
    this%nInitNode(1) = floor( xDim/xTetMax )
    this%nInitNode(2) = floor( yDim/yTetMax )
    this%nInitNode(3) = floor( zDim/zTetMax )
  end subroutine assume_nInit
  
  !======== locateInitTet ==========
  subroutine locateInitTet(this)
  !==================================
    implicit none
    type(UnstructGrid) :: this
    integer :: iEleCenter(3)
    double precision Node(this%nOutNode,3), EleCenter(3), CenterLoc(3), tempR
    double precision InitR(2*this%nInitNode(1),2*this%nInitNode(2),2*this%nInitNode(3))
    integer :: i,j,k, iTempInit
    
    !print*, 'Locating/Assigning Initial tetrahedron to each boxed...'
    do i=1,3
      Node(:,i) = this%OutNode(:,i)/this%dInitNode(i)
    enddo
    
    allocate(this%InitEle(2*this%nInitNode(1),2*this%nInitNode(2),2*this%nInitNode(3)))
    this%InitEle = -1
    InitR = 100d0
    
    ! locate centroid of the tetrahedrons to the boxes
    do i=1,this%nEle
      EleCenter(1) = sum(Node(This%Ele(i,:),1))/4d0 !center of mass
      EleCenter(2) = sum(Node(This%Ele(i,:),2))/4d0
      EleCenter(3) = sum(Node(This%Ele(i,:),3))/4d0
      iEleCenter(:) = ceiling(EleCenter(:))+this%nInitNode(:)
      CenterLoc = dble(iEleCenter(:)-this%nInitNode)-(/ 0.5d0,0.5d0,0.5d0 /)
      tempR = NORM2(EleCenter(:)-CenterLoc)
      if (tempR < InitR(iEleCenter(1),iEleCenter(2),iEleCenter(3))) then
        InitR(iEleCenter(1),iEleCenter(2),iEleCenter(3)) = tempR
        this%InitEle(iEleCenter(1),iEleCenter(2),iEleCenter(3)) = i
      endif
    enddo

    ! locate each node of tetrahedrons to the boxes when tetrahedron centroid is not located
    do j=1,4
      do i=1,this%nEle
        EleCenter = Node(This%Ele(i,j),:)
        iEleCenter(:) = ceiling(EleCenter(:))+this%nInitNode(:)
        do k=1,3
          if(iEleCenter(k)<1) then
            print*, i, iEleCenter, EleCenter
          endif
          if(iEleCenter(k)>2*this%nInitNode(k)) then
            print*, i, iEleCenter, EleCenter
          endif
        enddo
        iTempInit = this%InitEle(iEleCenter(1),iEleCenter(2),iEleCenter(3))
        if (iTempInit == -1) then
          CenterLoc = dble(iEleCenter(:)-this%nInitNode)-(/ 0.5d0,0.5d0,0.5d0 /)
          tempR = NORM2(EleCenter(:)-CenterLoc)
          if(tempR < InitR(iEleCenter(1),iEleCenter(2),iEleCenter(3))) then
            InitR(iEleCenter(1),iEleCenter(2),iEleCenter(3)) = tempR
            this%InitEle(iEleCenter(1),iEleCenter(2),iEleCenter(3)) = i
          endif
        endif
      enddo
    enddo
    
  end subroutine locateInitTet
  !======== end of locateInitTet ==========
  
  
  subroutine construct_UnstructGrid(this,fname)
    implicit none
    type (UnstructGrid), intent(out) :: this
    character(len=*), intent(in)  :: fname
    !integer,intent(in) :: fileID
    
    integer :: n, i,j,k
    double precision :: dum(9)
    !character(len=25) :: string
    
    print*, 'Constructing UnstructGrid...'
    !load core domain nodes
    !write(string,'(A8,I4.4,A8)') "unstruct", fileID, '.in.node'
    open(777,file=trim(fname)//'.in.node',action='read') 
    read(777,*) n, this%nInNode, this%dNode
    allocate(this%InNode(-this%nInNode(1):this%nInNode(1),&
                         -this%nInNode(2):this%nInNode(2),&
                         -this%nInNode(3):this%nInNode(3),6))   
    do k=-this%nInNode(3),this%nInNode(3)
      do j=-this%nInNode(2),this%nInNode(2)
        do i=-this%nInNode(1),this%nInNode(1)
          read(777,*) dum
          this%InNode(i,j,k,1:6) = dum(4:9)
        enddo
      enddo
    enddo
    close(777)
    print*, 'Core domain constructed...'
    
    !load halo domain nodes
    !write(string,'(A8,I4.4,A9)') "unstruct", fileID, '.out.node'
    open(777,file=trim(fname)//'.out.node',action='read') 
    read(777,*) this%nOutNode !, this%nInitNode, this%dNode
    allocate(this%OutNode(this%nOutNode,9))   
    do i=1,this%nOutNode
      read(777,*) this%OutNode(i,1:9)
    enddo 
    close(777)
    print*, '...'
        
    !load this%Ele
    !write(string,'(A8,I4.4,A4)') "unstruct", fileID, '.ele'
    open(777,file=trim(fname)//'.ele',action='read') 
    read(777,*) this%nEle
    allocate(this%Ele(this%nEle,4))   
    do i=1,this%nEle
      read(777,*) this%Ele(i,1:4)
    enddo 
    close(777)   
    print*, '...'
    
    
    !load this%Neigh 
    !write(string,'(A8,I4.4,A6)') "unstruct", fileID, '.neigh'
    open(777,file=trim(fname)//'.neigh',action='read') 
    read(777,*) j
    if(j/=this%nEle) then
      print *, 'the number of element in .neigh and .ele do not consistant' 
      call abort
    endif
    allocate(this%Neigh(j,4))   
    do i=1,j
      read(777,*) this%Neigh(i,1:4)
    enddo 
    close(777)   
    print*, 'Halo domain constructed...'
!    !load this%InitEle 
!    write(string,'(A8,I4.4,A5)') "unstruct", fileID, '.init'
!    open(777,file=string,action='read')
!    read(777,*)
!    allocate(this%InitEle( 2*this%nInitNode(1), &
!                           2*this%nInitNode(2), &
!                           2*this%nInitNode(3) ))
!    do k=1,2*this%nInitNode(3)
!      do j=1,2*this%nInitNode(2)
!        do i=1,2*this%nInitNode(1)
!          read(777,*) this%InitEle(i,j,k)
!        enddo
!      enddo
!    enddo
    
  end subroutine construct_UnstructGrid

  subroutine print_UnstructGrid(this,fileID)
    implicit none
    type (UnstructGrid), intent(in) :: this
    integer,intent(in) :: fileID
    character(len=25) :: str
          
          
    write (str, *) fileID
    str = "unstruct"//adjustl(str)
    open(777,file=str,form="unformatted",action='write')   
    
    !========Load Core Doamin ==========
    !write nodes (carteisn mesh)
    write(777) this%nInNode, this%dNode  ! nInNode=(nx,ny,nz), the total node is (2nx+1)(2ny+1)(2nz+1).  dInNode(dx,dy,dz) is the spacing between nodes
    write(777) this%InNode
    
    
    !========Load Halo Doamin ==========
    !write nodes (tetrahedral mesh)
    write(777) this%nOutNode ! nOutNode = number of node of tetrahedral grid
    write(777) this%OutNode
        
    !write element
    write(777) this%nEle
    write(777) this%Ele
    !write neigbor ( each element has 4 or less neighboring elements )  
    write(777) this%Neigh
    !write initial cadidate tetrahedral element to begin the tetrahedral searching algorithm (see subroutine InterpolateUnstructGrid)
    write(777) this%nInitNode, this%dInitNode
    write(777) this%InitEle
    close(777)

!    print *, '=======InNode========='
!    print *, this%nInNode
!    print *, this%dNode
!    print *, this%InNode(this%nInNode(1),-this%nInNode(2),-this%nInNode(3),:)
!    print *, this%InNode(8,-8,80,:)    
!    print *, '=======OutNode========='
!    print *, this%nOutNode
!    print *, this%OutNode(99999,:)
!    print *, this%OutNode(this%nOutNode,:)
!    print *, '=======Elem========='
!    print *, this%nEle
!    print *, this%Ele(99999,:)
!    print *, this%Ele(this%nEle,:)
!    print *, '=======Neigh========='
!    print *, this%Neigh(99999,:)
!    print *, this%Neigh(this%nEle,:)      
!    print *, '=======Neigh========='
!    print *, this%nInitNode
!    print *, this%dNode
!    print *, this%InitEle(20,-1,55)
!    print *, this%InitEle(this%nInitNode(1), &
!                           11, &
!                           this%nInitNode(3))

  endsubroutine print_unstructgrid
  
!  recursive subroutine sort(data2sort, iCol4Sort, nRow, nCol, iRowStart, iRowEnd)
!  ! quicksort algorithm
!    implicit none
!    double precision data2sort(nRow, nCol), mid, temp(nCol)
!    integer iCol4Sort, nCol, nRow, iRowStart, iRowEnd
!    integer i, j
    
!    mid = data2sort((iRowStart+iRowEnd)/2, iCol4Sort)
!    i = iRowStart
!    j = iRowEnd
!    do
!       do while (data2sort(i, iCol4Sort) < mid)
!          i=i+1      
!       end do
!       do while (mid < data2sort(j, iCol4Sort))
!          j=j-1
!       end do
!       if (i >= j) exit
!       temp = data2sort(i,:)
!       data2sort(i,:) = data2sort(j,:)
!       data2sort(j,:) = temp
!       i=i+1
!       j=j-1
!    end do
!    if (iRowStart < i-1) call sort(data2sort, iCol4Sort, nRow, nCol, iRowStart, i-1)
!    if (j+1 < iRowEnd)  call sort(data2sort, iCol4Sort, nRow, nCol, j+1, iRowEnd)
!  end subroutine sort
 
end program preprocess
