program tulli12
use mod_global
implicit none
integer :: i,j,total_trajectories,number_of_cores,TT
real*8 :: start,finish
real*8, dimension(:,:),allocatable :: population_matrix,corr_mat,therm_mat
!real*8, dimension(:),allocatable :: final_pop,final_imp,final_long

open(1, file ='input.txt')
read(1,*) number_of_cores,total_trajectories,iseed
close(1)

call setup_initial_values1
ntraj=int(total_trajectories/number_of_cores)


TT=int(total_time/dtc)

allocate(population_matrix(TT,ntraj))
allocate(corr_mat(TT,ntraj))
allocate(therm_mat(TT,ntraj))

!allocate(final_pop(TT))
!allocate(final_imp(TT))
!allocate(final_long(TT))

do traj_no=1,ntraj
call setup_initial_values2

!call draw_pes
call CPU_TIME(start)
!call classical_evolution 
call MD_dynamics

population_matrix(:,traj_no)=pop_mat
corr_mat(:,traj_no)=imp_mat
therm_mat(:,traj_no)=long_mat

call CPU_time(finish)
write(119,*)finish-start
enddo


call pop_averaging(population_matrix,150)!final_pop,150)
call pop_averaging(corr_mat,180)!final_imp,180)
call pop_averaging(therm_mat,140)!final_long,140)



!final_pop=0

!do j=1,ntraj
!   final_pop(:)=final_pop(:)+population_matrix(:,j)
!enddo


!do j=1,int(total_time/dtc)
!   if (mod(int(j*dtc),500).eq.0) then
!      write(21,*) j*dtc,final_pop(j)/real(ntraj)
 !  end if
!enddo
end program
!..............................................................................


