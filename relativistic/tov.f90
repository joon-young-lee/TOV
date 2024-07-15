module TOV_solver ! Relativistic, pure neutron
    ! no consideration of gravitational potential
    implicit none
    
    ! Constants and initial conditions
    real(8), parameter :: c = 2.99792e10_8
    real(8), parameter :: c2 = c ** 2
    real(8), parameter :: pi = 3.14159_8
    real(8), parameter :: G = 6.67259e-8_8
    real(8), parameter :: M0 = 1.98892e33_8
    real(8), parameter :: m_n = 1.67493e-24_8
    real(8), parameter :: hbar = 1.05457e-27_8
    real(8), parameter :: hbar3 = hbar ** 3
    real(8), parameter :: e_0 = m_n ** 4 * c ** 5/(pi ** 2 * hbar3)
contains

    function solve_TOV(p0, del_r, maxit) result(final_mass_radius_iteration)
        implicit none
        
        ! character(len=36), intent(in) :: file_name
        real(8), intent(in) :: p0, del_r
        integer, intent(in) :: maxit !k

        real(8) :: final_mass_radius_iteration(3)
        

        integer :: i, j
        real(8) :: p, m, r, r2, r3, r_4, p2, p3, p4, m2, m3, m4
        real(8) :: k1_p, k1_m, k2_p, k2_m, k3_p, k3_m, k4_p, k4_m
        

        ! Initial values
        p = p0
        m = 1.0e-20_8
        r = 1.0e-20_8
        

        i = 0
        do while (i <= maxit .and. p > 1.0e15_8)
            
            ! Calculate the Runge-Kutta coefficients for pressure and mass
            k1_p = dp_dr(r, p, m)
            k1_m = dm_dr(r, p, m)

            r2 = r + del_r/2.0_8
            p2 = p + del_r * k1_p / 2
            m2 = m + del_r * k1_m / 2
            k2_p = dp_dr(r2, p2, m2)
            k2_m = dm_dr(r2, p2, m2)


            r3 = r2
            p3 = p + del_r * k2_p / 2.0_8
            m3 = m + del_r * k2_m / 2.0_8
            k3_p = dp_dr(r3, p3, m3)
            k3_m = dm_dr(r3, p3, m3)

            r_4 = r2 + del_r/2.0_8
            p4 = p + del_r * k3_p
            m4 = m + del_r * k3_m
            k4_p = dp_dr(r_4, p4, m4)
            k4_m = dm_dr(r_4, p4, m4)

            ! Update pressure and mass using the weighted sum of the Runge-Kutta coefficients
            p = p + (k1_p + 2 * k2_p + 2 * k3_p + k4_p) / 6 * del_r
            m = m + (k1_m + 2 * k2_m + 2 * k3_m + k4_m) / 6 * del_r
            r = r + del_r


            ! Increment the iteration counter
            i = i + 1
        end do

        ! Output the final mass and radius
        final_mass_radius_iteration(1) = m/M0
        final_mass_radius_iteration(2) = r/1e5_8
        final_mass_radius_iteration(3) = i



        end function solve_TOV

    real(8) function dp_dr(r, p, m)
        implicit none
        real(8), intent(in) :: r, p, m
        
        dp_dr = - G * e(root(p)) * m * &
        (1 + p / e(root(p))) * &
        (1 + 4 * pi * p * r ** 3/(m * c2)) / &
        (1 - 2 * G * m / (c2 * r)) /&
        (c2 * r ** 2)
    end function dp_dr

    real(8) function dm_dr(r, p, m)
        implicit none
        real(8), intent(in) :: r, p, m
        ! TOV equation
        dm_dr = 4 * pi * r ** 2 * e(root(p)) / c2
        
    end function dm_dr

    real(8) function P(x) ! Function of pressure
        implicit none
        real(8), intent(in) :: x
        P = e_0/24 * ((2 * x ** 3 - 3 * x)*sqrt(x**2 + 1) + 3 * log(x + sqrt(x**2 + 1.0)))
        ! arcsinh = log(x + sqrt(x*x + 1.0))
    end function P
    
    real(8) function deriv(y) ! Derivative of p
        implicit none
        real(8), intent(in) :: y
        real(8), parameter :: h = 1.0e-6
        ! f'(x) = (f(x+h)-f(x))/h
        deriv = (P(y+h) - P(y-h))/(2 * h)
    end function deriv
    
    real(8) function root(p0)
    implicit none
    real(8), intent(in) :: p0
    real(8) :: x_old, x_new, tol
    integer :: j, maxit

    ! Parameters
    maxit = 1000
    tol = 1.0e-4_8

    ! Initial guess
    x_old = 1.0e-1_8

    j = 0
    do
        x_new = x_old - (P(x_old) - p0) / deriv(x_old)
        if (abs(x_new - x_old) < tol) then
            root = x_new
            return
        end if
        x_old = x_new
        j = j + 1
        if (j >= maxit) then
            root = x_new
            return
        end if
    end do
end function root

    real(8) function e(x) ! Function of energy density
        implicit none
        real(8), intent(in) :: x
        e = e_0/8 * ((2 * x ** 3 + x)*sqrt(x**2 + 1) - log(x + sqrt(x**2 + 1.0)))
        ! arcsinh = log(x + sqrt(x*x + 1.0))
    end function e

end module TOV_solver
