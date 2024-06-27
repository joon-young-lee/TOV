program main
    use TOV_solver
    implicit none

    real(8) :: p_i = 1.0e30_8, p_f = 1.0e31_8
    real(8) :: p_max = 3.5e35
    real(8), parameter :: del_r = 1.0e2_8
    integer, parameter :: maxit = 1e7
    integer, parameter :: step = 1000
    real(8) :: result(3)
    integer :: i, j
    
    ! Call the TOV function
    result = solve_TOV(p_max, del_r, maxit)
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
    do j = 0, 11
        open(unit=10, file='R_M.txt', status='old', access='append')
        do i = 1, step+1
            result = solve_TOV(p_i + (i-1) * (p_f - p_i)/step, del_r, maxit)
            write(10, *) p_i + (i-1) * (p_f - p_i)/step, result(1), result(2)
        end do
        p_i = p_i * 10
        p_f = p_f * 10
    end do
end program main

