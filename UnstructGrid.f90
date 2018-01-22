!----------------------------------------------------------------
! (c) Copyright, 2017 by the Regents of the University of California.
! UnstructGridClass : 
! Version: 
! Author: 
! Description: 
! Comments:
!----------------------------------------------------------------
      module UnstructGridClass
      
        type UnstructGrid
          integer :: Nseg, Mapstp, Itype, nEle, dataID
          double precision :: Length, zEdge, field_scale, w, theta0
        end type UnstructGrid

        type UnstructGridData
          integer :: nInNode(3), nInitNode(3), nOutNode, nEle, fileID
          double precision :: dInNode(3), dInitNode(3)
          integer, allocatable :: Ele(:,:), InitEle(:,:,:), Neigh(:,:)
          double precision, allocatable :: InNode(:,:,:,:), OutNode(:,:)          
        end type
        
                      
      contains
      
        subroutine construct_UnstructGrid(this,that,numseg,nmpstp,type,blength,zEdge,field_scale,freq,theta0,refptcl,fileID)
          implicit none
          include 'mpif.h'
          type (UnstructGrid), intent(out) :: this
          type (UnstructGridData), intent(inout) :: that(:)
          integer,intent(in) :: numseg,nmpstp,type,fileID
          double precision, intent(in) :: blength,zEdge,field_scale,freq,theta0,refptcl
          
          integer i, n_that, myrank, ierr
          logical flag_written
          
          flag_written = .False.
          n_that = size(that)
          
          do i=1,n_that
            if(fileID == that(i)%fileID) then
              this%dataID = i
              flag_written = .True.
              exit
            endif
            if(that(i)%fileID < 0) then
              call construct_UnstructGridData(that(i),fileID)
              this%dataID = i
              flag_written = .True.
              exit
            endif
          enddo
          if (.not. flag_written) then
            call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
            if(myrank==0)  print*, &
                            'error : maximum number of allowed for UnstructGridData exceeded. NUnstructGridDataMax = ',  n_that
            call abort()            
          endif
              
          this%Nseg = numseg
          this%Mapstp = nmpstp
          this%Itype = type
          this%Length = blength
          this%zEdge = zEdge
          this%field_scale = field_scale
          this%w = 2d0*3.141592653589793*freq
          this%theta0 = theta0*3.141592653589793/180.0d0 + refptcl
         

          
        end subroutine construct_UnstructGrid


        subroutine construct_UnstructGridData(this,fileID)
          implicit none
          include 'mpif.h'          
          integer, intent(in) :: fileID
          type (UnstructGridData), intent(out) :: this
          character(len=25) :: str
          this%fileID = fileID
          write (str, *) fileID
          str = "unstruct"//adjustl(str)
          open(777,file=str,form="unformatted",action='read')   
          
          !========Load Core Doamin ==========
          !load nodes (carteisn mesh)
          read(777) this%nInNode, this%dInNode  ! nInNode=(nx,ny,nz), the total node is (2nx+1)(2ny+1)(2nz+1).  dInNode(dx,dy,dz) is the spacing between nodes
          allocate(this%InNode(-this%nInNode(1):this%nInNode(1),&
                               -this%nInNode(2):this%nInNode(2),&
                               -this%nInNode(3):this%nInNode(3),6))  ! 6 elements : Ex,Ey,Ez,Bx,By,Bz
          read(777) this%InNode
          
          !========Load Halo Doamin ==========
          !load nodes (tetrahedral mesh)
          read(777) this%nOutNode ! nOutNode = number of node of tetrahedral grid        
          allocate(this%OutNode(this%nOutNode,9))  ! 9 elements : x,y,z,Ex,Ey,Ez,Bx,By,Bz  
          read(777) this%OutNode     
          
          !load element
          read(777) this%nEle
          allocate(this%Ele(this%nEle,4))  ! 4 elements are node id which is the array index of this%OutNode
          read(777) this%Ele          
          
          !load neigbor ( each element has 4 or less neighboring elements )
          allocate(this%Neigh(this%nEle,4))   
          read(777) this%Neigh 
          
          !load initial cadidate tetrahedral element to begin the tetrahedral searching algorithm (see subroutine InterpolateUnstructGrid)
          read(777) this%nInitNode, this%dInitNode
          allocate(this%InitEle( 2*this%nInitNode(1), &
                                 2*this%nInitNode(2), &
                                 2*this%nInitNode(3)) )
          read(777) this%InitEle
          close(777)
          
          !print*, 'initEle boundary', this%nInitNode*this%dInitNode
        end subroutine construct_UnstructGridData


        subroutine getparam_UnstructGrid(this,blength,bnseg,bmapstp,&
                                       btype)
          implicit none
          type (UnstructGrid), intent(in) :: this
          double precision, intent(out) :: blength
          integer, intent(out) :: bnseg,bmapstp,btype

          blength = this%Length
          bnseg = this%Nseg
          bmapstp = this%Mapstp
          btype = this%Itype

        endsubroutine getparam_UnstructGrid


       
        subroutine InterpolateUnstructGrid(this, that, inputCoordinate, pAssign, tmpfld)
          implicit none
          include 'mpif.h'
          type (UnstructGrid), intent(in) :: this
          type (UnstructGridData), intent(in) :: that
          double precision, intent(in) :: inputCoordinate(4)
          integer, intent(inout) :: pAssign
          double precision, intent(out) :: tmpfld(6)
          
          double precision fMod(3), A(3,3), V, x12(3), x13(3), x14(3), x23(3), x24(3), x34(3), TetDomain(4), rTest, field(4,6)
          integer nParticles, i, j, iEle, nMod(3), iMin, myid
          logical flagLost
          double precision :: pData(3), t
          
          call MPI_COMM_RANK(MPI_COMM_WORLD,myid,i)
          
          pData = inputCoordinate(1:3)
          pData(3) = pData(3) - this%zEdge  - 0.5d0*this%Length
          t = inputCoordinate(4)
          tmpfld = 0.0
          
!          print *, '=======InNode========='
!          print *, that%nInNode
!          print *, that%dInNode
!          print *, that%InNode(that%nInNode(1),-that%nInNode(2),-that%nInNode(3),:)
!          print *, that%InNode(8,-8,80,:)    
!          print *, '=======OutNode========='
!          print *, that%nOutNode
!          print *, that%OutNode(99999,:)
!          print *, that%OutNode(that%nOutNode,:)
!          print *, '=======Elem========='
!          print *, that%nEle
!          print *, that%Ele(99999,:)
!          print *, that%Ele(that%nEle,:)
!          print *, '=======Neigh========='
!          print *, that%Neigh(99999,:)
!          print *, that%Neigh(that%nEle,:)      
!          print *, '=======Init========='
!          print *, that%nInitNode
!          print *, that%dInitNode
!          print *, that%InitEle(20,-1,55)    
!          print *, that%InitEle(that%nInitNode(1), &
!                                 11, &
!                                 that%nInitNode(3))
          
          
          !if(pAssign == -1) return  ! commented because it is already handled before calling this routine
          if ( abs(pData(3)) <= that%nInNode(3)*that%dInNode(3) ) then
            if (abs(pData(1)) < that%nInitNode(1)*that%dInitNode(1) .and. abs(pData(2)) < that%nInitNode(2)*that%dInitNode(2) ) then
              !print*, i
              if (abs(pData(1)) > that%nInNode(1)*that%dInNode(1) .or. abs(pData(2)) > that%nInNode(2)*that%dInNode(2) ) then
                ! initialize assignment of particle to a tettrahedron that%Element
                if (pAssign == 0) then ! if not yet assigned
                  nMod = ceiling(pData/that%dInitNode) + that%nInitNode
                  pAssign = that%InitEle(nMod(1), nMod(2), nMod(3))
                  if(pAssign == -1) then
                    !print*, 'lost by init', pData, inputCoordinate(3)
                    !write(6907,*) pData
                    !call flush(6907)
                    return
                  endif
                endif
                ! tetrahedral searching algorithm
                do while(.true.)       
                  x12(1:3) = that%OutNode(that%Ele(pAssign,1),1:3) - that%OutNode(that%Ele(pAssign,2),1:3)
                  x13(1:3) = that%OutNode(that%Ele(pAssign,1),1:3) - that%OutNode(that%Ele(pAssign,3),1:3)
                  x14(1:3) = that%OutNode(that%Ele(pAssign,1),1:3) - that%OutNode(that%Ele(pAssign,4),1:3)
                  x23(1:3) = that%OutNode(that%Ele(pAssign,2),1:3) - that%OutNode(that%Ele(pAssign,3),1:3)
                  x24(1:3) = that%OutNode(that%Ele(pAssign,2),1:3) - that%OutNode(that%Ele(pAssign,4),1:3)
                  x34(1:3) = that%OutNode(that%Ele(pAssign,3),1:3) - that%OutNode(that%Ele(pAssign,4),1:3)
                  A(1,1) =-x14(3)*x34(2) +x34(3)*x14(2)
                  A(2,1) =-x14(3)*x12(2) +x12(3)*x14(2)
                  A(3,1) = x23(3)*x12(2) -x12(3)*x23(2)
                  A(1,2) =-x14(1)*x34(3) +x34(1)*x14(3)
                  A(2,2) =-x14(1)*x12(3) +x12(1)*x14(3)
                  A(3,2) = x23(1)*x12(3) -x12(1)*x23(3)
                  A(1,3) =-x14(2)*x34(1) +x34(2)*x14(1)
                  A(2,3) =-x14(2)*x12(1) +x12(2)*x14(1)
                  A(3,3) = x23(2)*x12(1) -x12(2)*x23(1)
                  V = -x12(1)*( x13(2)*x14(3)-x13(3)*x14(2)) &
                      -x13(1)*(-x12(2)*x14(3)+x12(3)*x14(2)) &
                      -x14(1)*( x12(2)*x13(3)-x12(3)*x13(2)) 

                  TetDomain(2:4) = pData - that%OutNode(that%Ele(pAssign,1),1:3)
                  TetDomain(2:4) = matmul(A,TetDomain(2:4) )/V               
                  TetDomain(1) = 1d0-sum(TetDomain(2:4))
                  rTest = minval(TetDomain)
                  if ( rTest < -1d-15 ) then
                    iMin = minloc(TetDomain,1)                    
                    pAssign = that%Neigh(pAssign, iMin)
                    if (pAssign==-1) then
                      !print*, 'lost by neighboring search-----', pData, inputCoordinate(3)
                      !write(6907,*) pData, inputCoordinate
                      !call flush(6907)                      
                      tmpfld = 0d0
                      return
                    endif
                  else
                    exit
                  endif
        !          if(j>100) then
        !            print *, TetDomain
        !            if(j>120) call abort()
        !          endif
                enddo ! end of searching algorithm
                iEle = pAssign
                field(1,1:6) = that%OutNode(that%Ele(iEle,1),4:9)
                do j=2,4
                  field(j,1:6) = that%OutNode(that%Ele(iEle,j),4:9) - field(1,1:6)
                enddo
                tmpfld = field(1,1:6) + field(2,1:6)*TetDomain(2) + field(3,1:6)*TetDomain(3) + field(4,1:6)*TetDomain(4)
              else
                if (pAssign > 0) pAssign = 0 !remove assignment
                ! linaer interpolation of hexahedron Node
                fMod = pData/that%dInNode
                nMod = floor(fMod)
                fMod = fMod - nMod
                tmpfld= that%InNode(nMod(1)  , nMod(2)  , nMod(3)  , :)*(1d0-fMod(1))*(1d0-fMod(2))*(1d0-fMod(3)) +&
                        that%InNode(nMod(1)+1, nMod(2)  , nMod(3)  , :)*(    fMod(1))*(1d0-fMod(2))*(1d0-fMod(3)) +&
                        that%InNode(nMod(1)+1, nMod(2)+1, nMod(3)  , :)*(    fMod(1))*(    fMod(2))*(1d0-fMod(3)) +&
                        that%InNode(nMod(1)+1, nMod(2)  , nMod(3)+1, :)*(    fMod(1))*(1d0-fMod(2))*(    fMod(3)) +&
                        that%InNode(nMod(1)+1, nMod(2)+1, nMod(3)+1, :)*(    fMod(1))*(    fMod(2))*(    fMod(3)) +&                            
                        that%InNode(nMod(1)  , nMod(2)+1, nMod(3)  , :)*(1d0-fMod(1))*(    fMod(2))*(1d0-fMod(3)) +&
                        that%InNode(nMod(1)  , nMod(2)+1, nMod(3)+1, :)*(1d0-fMod(1))*(    fMod(2))*(    fMod(3)) +&
                        that%InNode(nMod(1)  , nMod(2)  , nMod(3)+1, :)*(1d0-fMod(1))*(1d0-fMod(2))*(    fMod(3))  
              endif
            else
              pAssign=-1
              !print*, 'lost even before initialization-----', pData , inputCoordinate(3)
!              write(6907,*) pData, inputCoordinate
!              call flush(6907)              
            endif
          else 
            return
          endif
          if (this%w .ne. 0.0d0) then
            tmpfld(1:3) = tmpfld(1:3)*cos(this%w*t+this%theta0)
            tmpfld(4:6) = tmpfld(4:6)*sin(this%w*t+this%theta0)
          endif
          tmpfld = tmpfld*this%field_scale
          
          
        end subroutine InterpolateUnstructGrid
       
      end module UnstructGridClass
