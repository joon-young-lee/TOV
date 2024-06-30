program main
    use TOV_solver
    implicit none

    real(8) :: p_i = 1.0e30_8, p_f = 1.0e31_8
    real(8), parameter :: del_r = 1.0e2_8
    integer, parameter :: maxit = 1e7
    integer, parameter :: step = 1000
    real(8) :: result(3)
    integer :: i, j
    
    ! Call the TOV function
    result = solve_TOV(p_i, del_r, maxit)
    ! result(1) = mass
    ! result(1) = radius
    ! result(1) = iteration
    
    ! Output the results
    print *, 'Initial pressure', p_i, 'dyen/cm^2'
    print *, 'Delta r', del_r, 'cm'
    print *, 'Mass:', result(1), 'M_0'
    print *, 'Radius:', result(2), 'km'
    print *, 'Iteration', result(3), 'iterations'
    
    ! Write arrays to file
    open(unit=10, file='R_M_gp.txt', status='old', action='write')
    write(10, *) 'Initial Pressure, ', 'Mass, ', 'Radius'
    do j = 0, 11
        
        do i = 1, step+1
            result = solve_TOV(p_i + (i-1) * (p_f - p_i)/step, del_r, maxit)
            write(10, *) p_i + (i-1) * (p_f - p_i)/step, result(1), result(2)
            if (mod(i, 100) == 0) then
                print *, '---------------------------------------------------'
                print *, 'Initial pressure', p_i + (i-1) * (p_f - p_i)/step, 'dyen/cm^2'
                print *, 'Delta r', del_r, 'cm'
                print *, 'Mass:', result(1), 'M_0'
                print *, 'Radius:', result(2), 'km'
                print *, 'Iteration', result(3), 'iterations'
            end if
        end do
        p_i = p_i * 10
        p_f = p_f * 10
    end do
end program main

