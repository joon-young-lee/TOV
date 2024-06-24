program main
    use TOV_solver
    implicit none

    real(8), parameter :: p0 = 1.0e32_8
    real(8), parameter :: del_r = 1.0e1_8
    integer, parameter :: maxit = 1e6
    real(8) :: result(3)

    ! Call the TOV function
    result = solve_TOV(p0, del_r, maxit)

    ! Output the results
    print *, 'Initial pressure', p0
    print *, 'Delta r', del_r
    print *, "Final mass:", result(1), 'M_0'
    print *, "Final radius:", result(2), 'km'
    print *, 'Iteration', result(3), 'iterations'

    
end program main

