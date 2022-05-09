program tulli12
use mod_global
implicit none
integer :: i,j,k,total_trajectories,number_of_cores,TT,t,loop,ham,knot
real*8 :: x1,dx1,start,finish
real*8, dimension(:),allocatable :: population_matrix
real*8, dimension(:),allocatable :: Time_matrix
!call setup_initial_values1
open(1, file ='input.txt')
read(1,*) number_of_cores,total_trajectories,iseed
close(1)

call setup_initial_values1
ntraj=int(total_trajectories/number_of_cores)


TT=int(total_time/dtc)
allocate(Time_matrix(TT))
allocate(population_matrix(TT))
do j=1,TT
population_matrix(TT)=0
enddo
do i=1,ntraj
call setup_initial_values2




call CPU_TIME(start)

!open(13,file='Energy-time.dat')
open(12,file='vmd_file.xyz')
call classical_evolution 
!call velocity_verlet
!call MD_dynamics 
close(12)
!close(13)

do k=1,TT
population_matrix(k)=population_matrix(k)+population_mat(k)
enddo


call CPU_time(finish)
write(119,*)finish-start
enddo

do t=1,TT
Time_matrix(t)=10*t
enddo

!write(93,*) H

write(77,*)   0.0000000000,     1.000000000
do loop=1,TT
write(77,*) Time_matrix(loop),1-population_matrix(loop)/ntraj
enddo

end program
!..............................................................................


