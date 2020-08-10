Program genus_count
! count genus based on vtk file
! Using Euler characteristic X=vertices-edges+faces
! Computational Materials Science 161 (2019) 135–142
implicit none

real(8):: a,b,c,rnd_3,dseed_3

Character(50) :: dummy
integer :: first,second,min,max,n,w,i,j,p,np,nt,ne,cap,repeat,classify,non_zero_points_count,line_prop
real(8):: x1,x2,x3,face_count,point_count,line_count
Character(11) :: filein
Character(5) :: string
real(8),allocatable :: points(:,:)
integer, allocatable :: lines(:,:),lines_code(:,:)
integer :: element(3)

!cap points=-1,normal points=0,facepoints=1,edge points=2,corner points=3
!normal lines=0 face lines=1 edge lines=2

!non_zero_points_count=classify(DBLE(56.0),DBLE(13.0),DBLE(240.0),DBLE(250.0))
!print*,non_zero_points_count
!np=max(1,2)
!print*,np
!print*,abs(-250.0)
!stop
filein="surface.vtk"
open(unit=20,file=filein,form='formatted',status='unknown')
open(unit=60,file=filein,form='formatted',status='unknown')
open(unit=80,file=filein,form='formatted',status='unknown')

do i=1,4
    read(20,*)dummy
    read(60,*)dummy
enddo

read(20,*)dummy,np,dummy
print*,"In total ", np, "points"
allocate (points(0:np-1,0:3))
allocate (lines(0:np-1,0:30))
allocate (lines_code(0:np-1,0:30)) !0:total number 1:point number,2:properties (face,edge)

non_zero_points_count=0
do i=0,np-1
    read(20,*)x1,x2,x3
!    if (abs(x1)==250.0) print*,"found"
    points(i,0)=-1
    points(i,1)=x1
    points(i,2)=x2
    points(i,3)=x3
    if (classify(x1,x2,x3,DBLE(250.0)).ne.0) non_zero_points_count=non_zero_points_count+1
    lines(i,0)=0
    do j=1,30
        lines(i,j)=-1
    enddo
enddo
print*,non_zero_points_count
!stop

do i=0,150000000
    read(20,*)dummy
    if (dummy.eq.'CELLS') then
        exit
    endif
enddo
backspace(20)
read(20,*)dummy,ne,dummy

do i=0,150000000
    read(60,*)dummy
    if (dummy.eq.'CELL_DATA') then
        print*,"Cell_Data starts from",i
        exit
    endif
enddo
print*,dummy
backspace(60)
read(60,*)dummy,nt
print*,nt,ne
if(nt.eq.ne) print*,"Cells and Cell_DATA match"
read(60,*)dummy
read(60,*)dummy

face_count=0
line_count=0
point_count=0

do i=1,ne
    read(60,*)cap
    if (cap.eq.0) then
        face_count=face_count+1
        read(20,*)dummy,element(1),element(2),element(3)
        points(element(1),0)=classify(points(element(1),1),points(element(1),2),points(element(1),3),DBLE(250.0))
        !if (points(element(1),0)>=1) print*,"found"
        points(element(2),0)=classify(points(element(2),1),points(element(2),2),points(element(2),3),DBLE(250.0))
        points(element(3),0)=classify(points(element(3),1),points(element(3),2),points(element(3),3),DBLE(250.0))
!        points(element(1),0)=0
!        points(element(2),0)=0
!        points(element(3),0)=0

        !for points 1 and 2
        first=min(element(1),element(2))
        second=max(element(1),element(2))
        if (lines(first,0)==0) then
            lines(first,0)=1
            lines_code(first,0)=1
            lines(first,1)=second
            lines_code(first,1)=line_prop(int(points(first,0)),int(points(second,0)))
        else
            repeat=0
            do j=1,lines(first,0)
                if (lines(first,j)==second) then
                    repeat=1
                    exit
                endif
            enddo
            if (repeat==0) then
                lines(first,0)=lines(first,0)+1
                lines_code(first,0)=lines_code(first,0)+1
                lines(first,lines(first,0))=second
                lines_code(first,lines_code(first,0))=line_prop(int(points(first,0)),int(points(second,0)))
            endif
        endif

        !for points 1 and 3
        first=min(element(1),element(3))
        second=max(element(1),element(3))
        if (lines(first,0)==0) then
            lines(first,0)=1
            lines_code(first,0)=1
            lines(first,1)=second
            lines_code(first,1)=line_prop(int(points(first,0)),int(points(second,0)))
        else
            repeat=0
            do j=1,lines(first,0)
                if (lines(first,j)==second) then
                    repeat=1
                exit
                endif
            enddo
            if (repeat==0) then
                lines(first,0)=lines(first,0)+1
                lines(first,lines(first,0))=second
                lines_code(first,0)=lines_code(first,0)+1
                lines(first,lines(first,0))=second
                lines_code(first,lines_code(first,0))=line_prop(int(points(first,0)),int(points(second,0)))
            endif
        endif

        !for points 2 and 3
        first=min(element(2),element(3))
        second=max(element(2),element(3))
        if (lines(first,0)==0) then
            lines(first,0)=1
            lines_code(first,0)=1
            lines(first,1)=second
            lines_code(first,1)=line_prop(int(points(first,0)),int(points(second,0)))
        else
            repeat=0
            do j=1,lines(first,0)
                if (lines(first,j)==second) then
                    repeat=1
                    exit
                endif
            enddo
            if (repeat==0) then
                lines(first,0)=lines(first,0)+1
                lines_code(first,0)=lines_code(first,0)+1
                lines(first,lines(first,0))=second
                lines_code(first,lines_code(first,0))=line_prop(int(points(first,0)),int(points(second,0)))
            endif
        endif

    else
        !if(mod(i,5000)==0) print*,"cap nonzero,for these lines",i
        read(20,*)dummy
    endif
    if(mod(i,10000)==0) print*,"elements calculated:", face_count,"total element scaned:",i
enddo

point_count=0
do i=0,np-1
    if (points(i,0)==0) then
        point_count=point_count+1
    else if (points(i,0)==1) then
        point_count=point_count+0.5
    else if (points(i,0)==2) then
        point_count=point_count+0.25
    else if (points(i,0)==3) then
        point_count=point_count+0.125
    endif
enddo

line_count=0
do i=0,np-1
    do j=1,lines_code(i,0)
        if (lines_code(i,j)==0) then
            line_count=line_count+1
        else if (lines_code(i,j)==1) then
            line_count=line_count+0.5
        else if (lines_code(i,j)==2) then
            line_count=line_count+0.25
        endif
    enddo
enddo

print*,"total vertices : ",point_count
print*,"total lines : ",line_count
print*,"total faces :",face_count
close(20)
close(60)

print*, "Euler characteristic is :",point_count-line_count+face_count
print*, "genus is  :", 1-((point_count-line_count+face_count)/2.0)




end Program

integer function min(a,b)
integer :: a,b
if(a<b) then
    min=a
else
    min=b
endif

end function min

function max(a,b) result(m)
integer,intent(in) :: a,b

if(a>b) then
    m=a
else
    m=b
endif
!print*,a,b
end function max

integer function classify(a,b,c,d)
real(8) :: a,b,c,d
integer :: count
count=0
!print*,abs(a),abs(b),abs(c),d
!if(((d-delta)<abs(a)).and.(abs(a)<(d+delta))) then
if(abs(a)==d) then
    !print(abs(a)==d)
    count=count+1
endif
!if(((d-delta)<abs(b)).and.(abs(b)<(d+delta))) then
if(abs(b)==d) then
    !print(abs(b)==d)
    count=count+1
endif
!if(((d-delta)<abs(c)).and.(abs(c)<(d+delta))) then
if(abs(c)==d) then
    count=count+1
endif
classify=count
end function classify

integer function line_prop(a,b)
integer :: a,b
count=0
if(a>=1.and.b>=1) then
    count=count+1
endif
if(a>=2.and.b>=2) then
    count=count+1
endif
line_prop=count
end function line_prop




