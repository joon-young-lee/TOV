program main
    implicit none
    
end program mainprogram main
    use TOV_solver
    implicit none

    real(8), parameter :: p0 = 1.0e32_8
    real(8), parameter :: del_r = 100.0_8
    integer, parameter :: maxit = 1000000
    real(8) :: result(2)

    ! Call the TOV function
    result = solve_TOV(p0, del_r, maxit)

    ! Output the results
    print *, "Final mass:", result(1)
    print *, "Final radius:", result(2)

end program main

