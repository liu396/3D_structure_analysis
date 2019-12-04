!=============================================
!=================READ ME=====================
!This code is based on the paper: C. Liu and P. S. Branicio, Efficient generation of non-cubic stochastic periodic bicontinuous nanoporous structures, 
!Computational Materials Science, 2019, Volume 169, 109101 
!Written by Chang Liu
!The code is for the generation of stochastic bicontinuous nanoporous in parallelpiped with arbitary edge lengths and angles between edges. 
!This code will assign random field value f_r on every single atoms in the input cfg files, and create an output cfg file. 
!The output cfg file can be visulized in Ovito, and delete atoms belongs to porous region by setting level cut on random field value f_r. 
!The parameters should be set in a seperate file called parameters.in and there are some parameters to set:
!	1. n_waves: total number of waves that you initially want to distribute over a sphere. 
!	2. q_0: wave number that controls the ligament size
!   3. upper_limit and lower_limit: After adjust originally created waves using Fibonacci_grid, to select only waves with appropriate wave number. 
!	4. Up to 5 types of atom are allowed. Atom type and atomic mass should be specified in parameters.in
!   5. The total number of atoms should be divisible by the number of processes requested for MPI run. 

!=============================================

Module mymod
implicit none
Character(15) :: filein
Character(17) :: fileout
Character(5) :: filebin
Character(13) :: filepara
Character(13) :: v_map
Character(11) :: filevec
Character(13) ::file_vec_origin
Character(13) ::file_phase
integer :: id,inx,iny,inz,n,w,i,j,ii,n_waves,ratio,count,goal_int_1,goal_int_2,goal_int_3,repeated
integer,allocatable :: goal_int(:),global_atomtype(:)
real(8) :: q_0,pi,theta,rnd_3,r,dseed_3,upper_limit,lower_limit
real(8) :: M11,M12,M13,M22,M23,M33,alpha,beta,gamma,a,b,c,q_tmp_1,q_tmp_2,q_tmp_3
real(8) :: mass,mass_1,mass_2,mass_3,mass_4,mass_5
integer, allocatable :: bin_H(:)
real(8),allocatable :: q(:), phi(:), qq(:),cap_H(:),global_xpos(:),global_ypos(:),global_zpos(:),local_xpos(:),local_ypos(:),local_zpos(:),global_fr(:),local_fr(:)
Character(5) :: string
Character(2) :: element_1,element_2,element_3,element_4,element_5

contains

    integer function closest(number)
    real(8) :: number
    if (number>0) then
        if ((number-int(number)).gt.0.5) then
            closest=int(number)+1
        else
            closest=int(number)
        endif
    else
        if (abs(number-int(number)).gt.0.5) then
        closest=int(number)-1
        else
        closest=int(number)
        endif
    endif

    return
    end

end module

Program Fibonacci_generation

!Unit is Angstrom
use mymod
implicit none
include "mpif.h"

integer :: ierr,myid,nprocs


Character(50) :: dummy
integer :: nsd,is
real(8) :: t1,t2,x,y,z
real(8) :: h(3,3)




call MPI_INIT(ierr)
t1=MPI_WTIME()
call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)


filein="input.cfg"
fileout="output.cfg"
filepara='parameters.in'
file_phase="phase_info"

upper_limit=1000000
lower_limit=-1000000

element_1='Au'
element_2=''
element_3=''
element_4=''
element_5=''
mass_1=196.0
mass_2=0
mass_3=0
mass_4=0
mass_5=0

if (myid==0) then
	print*, "reading parameters...."
	
	open(unit=5,file=filepara,form="formatted",status="unknown")
	
	rewind(5)
	do i=1,100
		read(5,*,end=120)dummy
		if (dummy=='q_0') then
			backspace(5)
			read(5,*,err=120)dummy,q_0
			exit
		endif
		cycle
120     print*,'q_0 not defined'
		stop		
	enddo
	
	rewind(5)
	do i=1,100
		read(5,*,end=121)dummy
		if (dummy=='n_waves') then
			backspace(5)
			read(5,*,err=121)dummy,n_waves
			exit
		endif
		cycle
121     print*,'n_waves not defined'
		stop
	enddo
	
	rewind(5)
	do i=1,100
		read(5,*,end=122)dummy
		if (dummy=='upper_limit') then
			backspace(5)
			read(5,*,err=122)dummy,upper_limit
			exit
		endif
		cycle
122     print*,'upper_limit not defined. Taking 1000000'
		exit		
	enddo	
	
	rewind(5)
	do i=1,100
		read(5,*,end=123)dummy
		if (dummy=='lower_limit') then
			backspace(5)
			read(5,*,err=123)dummy,lower_limit
			exit
		endif
		cycle
123     print*,'lower_limit not defined. Taking 1000000'
		exit		
	enddo	
	
	rewind(5)
	do i=1,100
		read(5,*,end=124)dummy
		if (dummy=='element_1') then
			backspace(5)
			read(5,*,err=124)dummy,dummy
			element_1=trim(dummy)
			exit
		endif
		cycle
124     print*,'element_1 not defined'
		exit		
	enddo	
	
	rewind(5)
	do i=1,100
		read(5,*,end=125)dummy
		if (dummy=='element_2') then
			backspace(5)
			read(5,*,err=125)dummy,dummy
			element_2=trim(dummy)
			exit
		endif
		cycle
125     print*,'element_2 not defined'
		exit		
	enddo	
	
	rewind(5)
	do i=1,100
		read(5,*,end=126)dummy
		if (dummy=='element_3') then
			backspace(5)
			read(5,*,err=126)dummy,dummy
			element_3=trim(dummy)
			exit
		endif
		cycle
126     print*,'element_3 not defined'
		exit		
	enddo	
	print*,element_3
	
	rewind(5)
	do i=1,100
		read(5,*,end=127)dummy
		if (dummy=='element_4') then
			backspace(5)
			read(5,*,err=127)dummy,dummy
			element_4=trim(dummy)
			exit
		endif
		cycle
127     print*,'element_4 not defined'
		exit		
	enddo	
	
	rewind(5)
	do i=1,100
		read(5,*,end=128)dummy
		if (dummy=='element_5') then
			backspace(5)
			read(5,*,err=128)dummy,dummy
			element_5=trim(dummy)
			exit
		endif
		cycle
128     print*,'element_5 not defined'
		exit		
	enddo	
	
	rewind(5)
	do i=1,100
		read(5,*,end=129)dummy
		if (dummy=='mass_1') then
			backspace(5)
			read(5,*,err=129)dummy,mass_1
			exit
		endif
		cycle
129     print*,'mass_1 not defined, taking 0'
		exit		
	enddo	
	
	rewind(5)
	do i=1,100
		read(5,*,end=130)dummy
		if (dummy=='mass_2') then
			backspace(5)
			read(5,*,err=130)dummy,mass_2
			exit
		endif
		cycle
130     print*,'mass_2 not defined, taking 0'
		exit		
	enddo	
	
	rewind(5)
	do i=1,100
		read(5,*,end=131)dummy
		if (dummy=='mass_3') then
			backspace(5)
			read(5,*,err=131)dummy,mass_3
			exit
		endif
		cycle
131     print*,'mass_3 not defined, taking 0'
		exit		
	enddo	
	
	rewind(5)
	do i=1,100
		read(5,*,end=132)dummy
		if (dummy=='mass_4') then
			backspace(5)
			read(5,*,err=132)dummy,mass_4
			exit
		endif
		cycle
132     print*,'mass_4 not defined, taking 0'
		exit		
	enddo
	
	rewind(5)
	do i=1,100
		read(5,*,end=133)dummy
		if (dummy=='mass_5') then
			backspace(5)
			read(5,*,err=133)dummy,mass_5
			exit
		endif
		cycle
133     print*,'mass_5 not defined, taking 0'
		exit		
	enddo	
	
	
		
	print*,"q_0: ",q_0
	print*,"n_waves: ",n_waves
	print*,"upper_limit: ",upper_limit
	print*,"lower_limit: ",lower_limit

	print*,"element_1: ",element_1
	print*,"element_2: ",element_2
	print*,"element_3: ",element_3
	print*,"element_4: ",element_4
	print*,"element_5: ",element_5

	print*,"mass_1: ",mass_1
	print*,"mass_2: ",mass_2
	print*,"mass_3: ",mass_3
	print*,"mass_4: ",mass_4
	print*,"mass_5: ",mass_5
	
	close(5)
endif

call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call MPI_BCAST(q_0,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(n_waves,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(lower_limit,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(upper_limit,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(mass_1,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(mass_2,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(mass_3,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(mass_4,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(mass_5,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(element_1,2,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(element_2,2,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(element_3,2,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(element_4,2,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(element_5,2,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)


allocate(q(3*n_waves))
allocate(qq(3*n_waves))
allocate(phi(n_waves))
allocate(goal_int(3*n_waves))

pi=3.1415926535897932384626433


open(unit=50,file=filein,form="formatted",status="unknown")
rewind(50)
read(50,*)dummy,dummy,dummy,dummy,n
read(50,*)dummy
do i=1,3
    do j=1,3
    read(50,*)dummy,dummy,dummy,h(i,j),dummy
    enddo
enddo

allocate(global_atomtype(n))
allocate(global_xpos(n))
allocate(global_ypos(n))
allocate(global_zpos(n))
allocate(global_fr(n))
nsd=n/nprocs
print*,n,nprocs,nsd
allocate(local_xpos(nsd))
allocate(local_ypos(nsd))
allocate(local_zpos(nsd))
allocate(local_fr(nsd))

if (myid==0) then
	
	print*, "calculating edges and angles"
    a=sqrt(h(1,1)*h(1,1)+h(1,2)*h(1,2)+h(1,3)*h(1,3))
    b=sqrt(h(2,1)*h(2,1)+h(2,2)*h(2,2)+h(2,3)*h(2,3))
    c=sqrt(h(3,1)*h(3,1)+h(3,2)*h(3,2)+h(3,3)*h(3,3))

    alpha=acos((h(2,1)*h(3,1)+h(2,2)*h(3,2)+h(2,3)*h(3,3))/b/c)

    beta=acos((h(1,1)*h(3,1)+h(1,2)*h(3,2)+h(1,3)*h(3,3))/a/c)

    gamma=acos((h(1,1)*h(2,1)+h(1,2)*h(2,2)+h(1,3)*h(2,3))/a/b)

    print*,"a,b,c equal to",a,b,c
    print*,"alpha,beta,gamma equal to", alpha/pi*180,beta/pi*180,gamma/pi*180

! Another way to generate random phases

!    call random_seed()
!    do i=1,n_waves
!        call random_number(rnd_3)
!        phi(i)=rnd_3*2*pi
!    enddo

! Another way to generate random phases
	open(unit=11,file=file_phase,form="formatted",status="unknown")
    dseed_3=43572d0
    i=1
    do while(i.lt.n_waves)
        call myrnd(rnd_3,dseed_3)
        phi(i)=rnd_3*2*pi
		write(11,'(f10.7)')phi(i)
        i=i+1
    enddo


    file_vec_origin="vec_ori.txt"
    open(unit=15,file=file_vec_origin,form="formatted",status="unknown")
    do i=1,n_waves
        z=-1+(i-1)*2.0/(n_waves-1)
        r=sqrt(1-z*z)
        theta=i*137.508/180*pi
        x=r*cos(theta)
        y=r*sin(theta)
        q(3*i)=z*q_0
        q(3*i-1)=y*q_0
        q(3*i-2)=x*q_0
        write(15,'(f10.7,f10.7,f10.7)')q(3*i-2),q(3*i-1),q(3*i)
    enddo
    close(15)

    M11=a
    M12=b*cos(gamma)
    M13=c*cos(beta)
    M33=c*sqrt(1-(cos(alpha))**2-(cos(beta))**2-(cos(gamma))**2+2*cos(alpha)*cos(beta)*cos(gamma))/sin(gamma)
    M22=b*sin(gamma)
    M23=c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma)

    count=0
	print*,"generating waves"
    do i=1,n_waves
        if (i==1) then
            count=count+1
            goal_int_1=closest(q(3*i-2)*M11)
            goal_int(3*count-2)=goal_int_1
            qq(3*count-2)=goal_int_1/M11


            goal_int_2=closest(qq(3*i-2)*M12+q(3*i-1)*M22)
            goal_int(3*count-1)=goal_int_2
            qq(3*i-1)=(goal_int_2-qq(3*i-2)*M12)/M22


            goal_int_3=closest(qq(3*i-2)*M13+qq(3*i-1)*M23+q(3*i)*M33)
            goal_int(3*count)=goal_int_3
            qq(3*i)=(goal_int_3-qq(3*i-2)*M13-qq(3*i-1)*M23)/M33

        else
            goal_int_1=closest(q(3*i-2)*M11)
            q_tmp_1=goal_int_1/M11

            goal_int_2=closest(q_tmp_1*M12+q(3*i-1)*M22)
            q_tmp_2=(goal_int_2-q_tmp_1*M12)/M22

            goal_int_3=closest(q_tmp_1*M13+q_tmp_2*M23+q(3*i)*M33)
            q_tmp_3=(goal_int_3-q_tmp_1*M13-q_tmp_2*M23)/M33

            repeated=0
            do j=1,i-1
                if (((goal_int_1).eq.(goal_int(3*j-2))).and.((goal_int_2).eq.(goal_int(3*j-1))).and.((goal_int_3).eq.(goal_int(3*j)))) then
                    repeated=1
                    exit
                endif
            enddo

            if ((repeated.ne.1).and.(((q_tmp_1)**2+(q_tmp_2)**2+(q_tmp_3)**2).lt.upper_limit**2).and.(((q_tmp_1)**2+(q_tmp_2)**2+(q_tmp_3)**2).gt.lower_limit**2)) then
                count=count+1
                goal_int(3*count-2)=goal_int_1
                goal_int(3*count-1)=goal_int_2
                goal_int(3*count)=goal_int_3
                qq(3*count-2)=q_tmp_1
                qq(3*count-1)=q_tmp_2
                qq(3*count)=q_tmp_3

            endif
        endif
    enddo
    print*,"there are in total ",count," waves saved"
!!!===============================================

    call draw_waves_v2(count) !After adjusted the waves to be compatible with the periodicity, the waves are saved in a cfg file called direction_int 
    call draw_waves(n_waves) !Original Fibonacci Grid is saved in a cfg file called direction_sph
    call write_Hbin(count) !The distribution of wave-number-related value H is saved in H.bin
    call write_vectors(count,1) !All final wave vectors are saved in vectors.txt. The numbering of the vectors starts from the second argument of the subroutine.

	print*,"master process reading input file..."
	
    read(50,*)dummy
    read(50,*)dummy,dummy,i
    i=i-3
    do j=1,i
      read(50,*)dummy,dummy,dummy
    enddo
	
	is=0
	i=1
	do while (i<n+1)
		if (mod(i,10000)==0) then
			print*,i
		endif
		read(50,*,end=300,err=200) global_xpos(i),global_ypos(i),global_zpos(i)
		global_atomtype(i)=is 
		i=i+1
		cycle
200 	continue
		backspace(50)
		backspace(50)
		read(50,*)dummy
		read(50,*)dummy
		string=trim(dummy)
    	if(string==element_1) then
        	is=1
    	elseif(string==element_2) then
        	is=2
    	elseif(string==element_3) then
        	is=3
		elseif(string==element_4) then
			is=4
		elseif(string==element_5) then
			is=5
    	endif
		cycle
300 	print*,"End of file found at atom:",i
		exit
	enddo
	
	print*,"distributing information to other processes"
endif
!stop

call MPI_SCATTER(global_xpos,nsd,MPI_REAL8,local_xpos,nsd,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
call MPI_SCATTER(global_ypos,nsd,MPI_REAL8,local_ypos,nsd,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
call MPI_SCATTER(global_zpos,nsd,MPI_REAL8,local_zpos,nsd,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call MPI_BCAST(count,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(phi,n_waves,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(q,3*n_waves,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(qq,3*count,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

close(50)

if (myid==0) print*,"calculating..."
print*,nsd,count
do i=1,nsd
	if (myid==0) then
		if (mod(i,2000)==0) print*,i*nprocs
	endif
	local_fr(i)=0
	do j=1,count
		local_fr(i)=local_fr(i)+cos(2*pi*qq(j*3-2)*(local_xpos(i)*h(1,1)+local_ypos(i)*h(2,1)+local_zpos(i)*h(3,1))+2*pi*qq(j*3-1)*(local_xpos(i)*h(1,2)+local_ypos(i)*h(2,2)+local_zpos(i)*(h(3,2)))+2*pi*qq(j*3)*(local_xpos(i)*h(1,3)+local_ypos(i)*h(2,3)+local_zpos(i)*h(3,3))+phi(j))
	enddo
	local_fr(i)=sqrt(2.0/count)*local_fr(i)
enddo

call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call MPI_GATHER(local_fr,nsd,MPI_REAL8,global_fr,nsd,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

if(myid==0) then
	print*,"writing output file..."
	open(unit=60,file=fileout,form="formatted",status="unknown")

	rewind(60)
	write(60,'("Number of particles = ",i9)') n
	write(60,'("A = 1.0 Angstrom (basic length-scale)")')

	do i=1,3
    	do j=1,3
        	write(60,'("H0(",i1,",",i1,") = ",f10.3," A")')i,j,h(i,j)
    	enddo
	enddo


	write(60,'(".NO_VELOCITY.")')
	write(60,'("entry_count = ",i2 )') 6
	write(60,'("auxiliary[0] = type")')
	write(60,'("auxiliary[1] = id")')
	write(60,'("auxiliary[2] = f_r")')

    if (global_atomtype(1)==1) then
		mass=mass_1
		string=element_1
	else if (global_atomtype(1)==2) then
		mass=mass_2
		string=element_2
	else if (global_atomtype(1)==3) then
		mass=mass_3
		string=element_3
	else if (global_atomtype(1)==4) then
		mass=mass_4
		string=element_4
	else if (global_atomtype(1)==5) then
		mass=mass_5
		string=element_5
	endif

	write(60,'(f17.10)') mass
	write(60,*) string
	write(60,'(es14.6,es14.6,es14.6,i3,i10,f17.10)') global_xpos(1), global_ypos(1),global_zpos(1),global_atomtype(1),1,global_fr(1)
	
	do i=2,n
		if (mod(i,10000)==0) print*,i
		if (global_atomtype(i)==global_atomtype(i-1)) then
			write(60,'(es14.6,es14.6,es14.6,i3,i10,f17.10)') global_xpos(i), global_ypos(i),global_zpos(i),global_atomtype(i),i,global_fr(i)
		else
		    if (global_atomtype(i)==1) then
				mass=mass_1
				string=element_1
			else if (global_atomtype(i)==2) then
				mass=mass_2
				string=element_2
			else if (global_atomtype(i)==3) then
				mass=mass_3
				string=element_3
			else if (global_atomtype(i)==4) then
				mass=mass_4
				string=element_4
			else if (global_atomtype(i)==5) then
				mass=mass_5
				string=element_5
			endif
			write(60,'(f17.10)') mass
			write(60,*) string
			write(60,'(es14.6,es14.6,es14.6,i3,i10,f17.10)') global_xpos(i), global_ypos(i),global_zpos(i),global_atomtype(i),i,global_fr(i)
		endif
	enddo

endif



call MPI_BARRIER(MPI_COMM_WORLD,ierr)
t2=MPI_WTIME()

if (myid==0) print*, "total time is ", t2-t1
call MPI_FINALIZE(ierr)


end program


subroutine myrnd(rnd,dseed)
real(8) rnd,dseed
!--------------------------------------------------------------------------------------------------
!  Random-number generator.
!--------------------------------------------------------------------------------------------------
real(8) d2p31m,d2p31
data d2p31m/2147483647d0/
data d2p31 /2147483648d0/
dseed=mod(dble(16807d0*dseed),dble(d2p31m))
rnd=dseed/d2p31
!print*,rnd
return
end
!==========================================



subroutine write_Hbin(waves)
    use mymod
    integer :: waves,ww,int_min_of,int_max_of
    real(8) :: min_of,max_of


    allocate(cap_H(waves)) !This is for H.bins

    do ww=1,waves
        cap_H(ww)=(qq(3*ww-2)*a)**2+(qq(3*ww-1)*a)**2+(qq(3*ww)*a)**2
    enddo
    print*,"capital H calculated"

    min_of=cap_H(1)
    do ww=2,waves
        if (cap_H(ww)<min_of) min_of=cap_H(ww)
    enddo
    print*,"min of H found",min_of

    max_of=cap_H(1)
    do ww=2,waves
        if (cap_H(ww)>max_of) max_of=cap_H(ww)
    enddo
    print*,"max of H found",max_of

    int_min_of=int(min_of)
    int_max_of=int(max_of)
    allocate(bin_H(int_min_of:int_max_of))

    bin_H=0
    do ww=1,waves
        bin_H(int(cap_H(ww)))=bin_H(int(cap_H(ww)))+1
    enddo

    filebin="H.bin"
    open(unit=20,file=filebin,form="formatted",status="unknown")
    print*,int_min_of,int_max_of
    rewind(20)
    do ww=int_min_of,int_max_of
        write(20,'(i7,i7)')ww,bin_H(ww)
    enddo
    close(20)
    print*,"bin statistics are written"
    deallocate(cap_H)
    deallocate(bin_H)

end subroutine
!==========================================

subroutine draw_waves(waves)
    use mymod
    integer :: waves
    real(8) :: box(3,3)

    v_map="direction_sph"
    box(1,1)=50000*q_0
    box(1,2)=0
    box(1,3)=0
    box(2,1)=0
    box(2,2)=50000*q_0
    box(2,3)=0
    box(3,1)=0
    box(3,2)=0
    box(3,3)=50000*q_0
    open(unit=70,file=v_map,form="formatted",status="unknown")
    rewind(70)

    write(70,'("Number of particles = ",i9)') waves
    write(70,'("A = 1.0 Angstrom (basic length-scale)")')

    do ii=1,3
    do j=1,3
    write(70,'("H0(",i1,",",i1,") = ",f10.3," A")')j,ii,box(j,ii)
    enddo
    enddo

    write(70,'(".NO_VELOCITY.")')
    write(70,'("entry_count = ",i2 )') 4
    write(70,'("auxiliary[1] = id")')

    write(70,'("196.0")')
    write(70,'("Au")')

    do j=1,waves
        write(70,'(es14.6,es14.6,es14.6,i10,i7)') q(3*j-2)/2.0/q_0+0.5,q(3*j-1)/2.0/q_0+0.5,q(j*3)/2.0/q_0+0.5,j
    enddo
    close(70)
    print*,"spherical waves are drwan"
end subroutine
!==========================================

subroutine draw_waves_v2(waves)
use mymod
integer :: waves
real(8) :: box(3,3)

v_map="direction_int"
box(1,1)=50000*q_0
box(1,2)=0
box(1,3)=0
box(2,1)=0
box(2,2)=50000*q_0
box(2,3)=0
box(3,1)=0
box(3,2)=0
box(3,3)=50000*q_0
open(unit=70,file=v_map,form="formatted",status="unknown")
rewind(70)

write(70,'("Number of particles = ",i9)') waves
write(70,'("A = 1.0 Angstrom (basic length-scale)")')

do ii=1,3
do j=1,3
write(70,'("H0(",i1,",",i1,") = ",f10.3," A")')j,ii,box(j,ii)
enddo
enddo

write(70,'(".NO_VELOCITY.")')
write(70,'("entry_count = ",i2 )') 4
write(70,'("auxiliary[1] = id")')

write(70,'("196.0")')
write(70,'("Au")')

do j=1,waves
write(70,'(es14.6,es14.6,es14.6,i10,i7)') qq(3*j-2)/2.0/q_0+0.5,qq(3*j-1)/2.0/q_0+0.5,qq(j*3)/2.0/q_0+0.5,j
enddo
close(70)
print*,"Integer waves are drawn"
end subroutine
!==========================================


subroutine convert_2_integer_v1(waves)
use mymod
integer :: trailx,traily,trailz,minx,miny,minz
real(8) :: distance
distance=1000
do i=1,waves
    trainlx=int(abs(q(3*i-2)))
    trainly=int(abs(q(3*i-1)))
    trailz=int(abs(q(3*i)))
    if (distance>minx**2+25*(miny**2)+25*(minz**2)-(q_0**2)*25) then
        distance=minx**2+25*(miny**2)+25*(minz**2)
        minx=trailx
        miny=traily
        minz=trailz
    endif

    trainlx=int(abs(q(3*i-2)))+1
    trainly=int(abs(q(3*i-1)))
    trailz=int(abs(q(3*i)))
    if (distance>minx**2+25*(miny**2)+25*(minz**2)-(q_0**2)*25) then
        distance=minx**2+25*(miny**2)+25*(minz**2)
        minx=trailx
        miny=traily
        minz=trailz
    endif

    trainlx=int(abs(q(3*i-2)))
    trainly=int(abs(q(3*i-1)))+1
    trailz=int(abs(q(3*i)))
    if (distance>minx**2+25*(miny**2)+25*(minz**2)-(q_0**2)*25) then
        distance=minx**2+25*(miny**2)+25*(minz**2)
        minx=trailx
        miny=traily
        minz=trailz
    endif

    trainlx=int(abs(q(3*i-2)))
    trainly=int(abs(q(3*i-1)))
    trailz=int(abs(q(3*i)))+1
    if (distance>minx**2+25*(miny**2)+25*(minz**2)-(q_0**2)*25) then
    distance=minx**2+25*(miny**2)+25*(minz**2)
    minx=trailx
    miny=traily
    minz=trailz
    endif

    trainlx=int(abs(q(3*i-2)))+1
    trainly=int(abs(q(3*i-1)))+1
    trailz=int(abs(q(3*i)))
    if (distance>minx**2+25*(miny**2)+25*(minz**2)-(q_0**2)*25) then
    distance=minx**2+25*(miny**2)+25*(minz**2)
    minx=trailx
    miny=traily
    minz=trailz
    endif

    trainlx=int(abs(q(3*i-2)))+1
    trainly=int(abs(q(3*i-1)))
    trailz=int(abs(q(3*i)))+1
    if (distance>minx**2+25*(miny**2)+25*(minz**2)-(q_0**2)*25) then
        distance=minx**2+25*(miny**2)+25*(minz**2)
        minx=trailx
        miny=traily
        minz=trailz
    endif

    trainlx=int(abs(q(3*i-2)))
    trainly=int(abs(q(3*i-1)))+1
    trailz=int(abs(q(3*i)))+1
    if (distance>minx**2+25*(miny**2)+25*(minz**2)-(q_0**2)*25) then
        distance=minx**2+25*(miny**2)+25*(minz**2)
        minx=trailx
        miny=traily
        minz=trailz
    endif

    trainlx=int(abs(q(3*i-2)))+1
    trainly=int(abs(q(3*i-1)))+1
    trailz=int(abs(q(3*i)))+1
    if (distance>minx**2+25*(miny**2)+25*(minz**2)-(q_0**2)*25) then
        distance=minx**2+25*(miny**2)+25*(minz**2)
        minx=trailx
        miny=traily
        minz=trailz
    endif

    qq(3*i-2)=minx
    qq(3*i-1)=miny
    qq(3*i)=minz
enddo
print*,"converted to integer waves"
end subroutine
!==========================================

subroutine convert_2_integer_v2(waves)
use mymod
integer :: waves
do i=1,3*waves
    qq(i)=closest(q(i))
enddo
print*,"converted to integer waves"
end subroutine
!===============================================


subroutine delete_repeat(waves)
use mymod
integer,allocatable :: q_swap(:)
integer :: waves,flag
allocate(q_swap(3*waves))
count=0
count=count+1
q_swap(3*1-2)=qq(3*1-2)
q_swap(3*1-1)=qq(3*1-1)
q_swap(3*1)=qq(3*1)
do j=2,waves
    flag=0
    do i=1,count
        if((qq(3*j-2).eq.q_swap(3*i-2)).and.(qq(3*j-1).eq.q_swap(3*i-1)).and.(qq(3*j).eq.q_swap(3*i))) then
            flag=1
            exit
        endif
    enddo
    if (flag==0) then
        count=count+1
        q_swap(3*count-2)=qq(3*j-2)
        q_swap(3*count-1)=qq(3*j-1)
        q_swap(3*count)=qq(3*j)
    endif
enddo

deallocate(qq)
allocate(qq(3*count))

do j=1,count
qq(3*j-2)=q_swap(3*j-2)
qq(3*j-1)=q_swap(3*j-1)
qq(3*j)=q_swap(3*j)
!print*,q(j,1:3)
enddo

deallocate(q_swap)

end subroutine
!===========================================

subroutine write_vectors(waves,start)
use mymod
integer :: waves,start
filevec="vectors.txt"
open(unit=10,file=filevec,form="formatted",status="unknown")

do i=1,waves
    write(10,'("q(",i5,",1:3) = (/",f10.7,",",f10.7,",",f10.7,"/)")') i+start-1,qq(3*i-2),qq(3*i-1),qq(3*i)
!    write(10,'(f10.7,f10.7,f10.7)') qq(3*i-2),qq(3*i-1),qq(3*i)
enddo

close(10)


end subroutine
