module TOV_solver
    implicit none

    ! Constants and initial conditions
    real(8), parameter :: gamma = 5.0_8 / 3.0_8
    real(8), parameter :: c = 2.99792458e10_8
    real(8), parameter :: c2 = c ** 2
    real(8), parameter :: K = 6.428e-26_8
    real(8), parameter :: pi = 3.1415926543_8
    real(8), parameter :: G = 6.67259e-8_8
    real(8), parameter :: M0 = 1.9884e33_8
contains

    function solve_TOV(p0, del_r, maxit) result(final_mass_radius_iteration)
        implicit none
        real(8), intent(in) :: p0, del_r
        integer, intent(in) :: maxit
        real(8) :: final_mass_radius_iteration(3)

        integer :: i, j
        real(8) :: p, m, r, p_new, m_new, r_new
        real(8) :: k1_p, k1_m, k2_p, k2_m, k3_p, k3_m, k4_p, k4_m
        real(8), allocatable :: r_arr(:), p_arr(:), m_arr(:)

        ! Allocate arrays
        allocate(r_arr(0:maxit+1), p_arr(0:maxit+1), m_arr(0:maxit+1))

        ! Initial values
        p_arr(0) = p0
        m_arr(0) = 1.0e-10_8
        r_arr(0) = 1.0e-10_8

        i = 0
        do while (i <= maxit .and. p_arr(i) > 0.0_8)
            r = r_arr(i)
            p = p_arr(i)
            m = m_arr(i)

            ! Calculate the Runge-Kutta coefficients for pressure and mass
            k1_p = dp_dr(r, p, m)
            k1_m = dm_dr(r, p, m)

            k2_p = dp_dr(r + del_r / 2.0_8, p + k1_p / 2.0_8, m + k1_m / 2.0_8)
            k2_m = dm_dr(r + del_r / 2.0_8, p + k1_p / 2.0_8, m + k1_p / 2.0_8)

            k3_p = dp_dr(r + del_r / 2.0_8, p + k2_p / 2.0_8, m + k2_m / 2.0_8)
            k3_m = dm_dr(r + del_r / 2.0_8, p + k2_p / 2.0_8, m + k2_p / 2.0_8)

            k4_p = dp_dr(r + del_r, p + k3_p, m + k3_m)
            k4_m = dm_dr(r + del_r, p + k3_p, m + k3_p)

            ! Update pressure and mass using the weighted sum of the Runge-Kutta coefficients
            p_new = p_arr(i) + (k1_p + 2.0_8 * k2_p + 2.0_8 * k3_p + k4_p) / 6.0_8 * del_r
            m_new = m_arr(i) + (k1_m + 2.0_8 * k2_m + 2.0_8 * k3_m + k4_m) / 6.0_8 * del_r
            r_new = r_arr(i) + del_r

            ! Append the new values to the arrays
            p_arr(i+1) = p_new
            m_arr(i+1) = m_new
            r_arr(i+1) = r_new

            ! Increment the iteration counter
            i = i + 1
        end do

        ! Output the final mass and radius
        final_mass_radius_iteration(1) = m_arr(i)/M0
        final_mass_radius_iteration(2) = r_arr(i)/1e5_8
        final_mass_radius_iteration(3) = i

    


        ! Deallocate arrays
        deallocate(r_arr, p_arr, m_arr)
    end function solve_TOV

    real(8) function dp_dr(r, p, m)
        implicit none
        real(8), intent(in) :: r, p, m
        ! Your derivative calculation here
        dp_dr = - G * (p/K) ** (1/gamma) * m * &
        (1 + p/((p/K) ** (1/gamma))) * &
        (1 + 4 * pi * p * r ** 3/(m * c2)) / &
        (1 - 2 * G * m / (c2 * r)) / ((c * r) ** 2)  ! Example calculation
    end function dp_dr

    real(8) function dm_dr(r, p, m)
        implicit none
        real(8), intent(in) :: r, p, m
        ! Your derivative calculation here
        dm_dr = 4 * pi * r**2 * (p/K) ** (1/gamma) / c2
    end function dm_dr

end module TOV_solver
