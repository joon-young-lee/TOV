module TOV_solver ! Relativistic, pure neutron
    implicit none
    ! p = e_0/8 * (2 * x ** 3 + x)(x**2 + 1)**(1/2) - arcsinh(x)
    ! e = e_0/24 * (2 * x ** 3 - 3 * x)(x**2 + 1)**(1/2) + 3 * arcsinh(x)
    ! Constants and initial conditions
    real(8), parameter :: c = 2.99792458e10_8
    real(8), parameter :: c2 = c ** 2
    real(8), parameter :: pi = 3.14159265359_8
    real(8), parameter :: G = 6.67259e-8_8
    real(8), parameter :: M0 = 1.9884e33_8
    real(8), parameter :: m_n = 1.6749e-24_8
    real(8), parameter :: hbar = 1.0546e-27_8
    real(8), parameter :: hbar3 = hbar ** 3
    real(8), parameter :: e_0 = m_n ** 4 * c ** 5/(pi ** 2 * hbar3)

    print *, e(root(1.0e30))
    
contains


    real(8) function P(x) ! Function of pressure
        implicit none
        real(8), intent(in) :: x
        P = e_0/8 * (2 * x ** 3 + x)*(x**2 + 1)**(1/2) - log(x + sqrt(x**2 + 1.0))
        ! arcsinh = log(x + sqrt(x*x + 1.0))
    end function P
    
    real(8) function deriv(y) ! Derivative of p
        implicit none
        real(8), intent(in) :: y
        real(8), parameter :: h = 1.0e-5
        ! f'(x) = (f(x+h)-f(x))/h
        deriv = (P(y+h) - P(y-h))/(2 * h)
    end function deriv
    
    real(8) function root(l)
    implicit none
    real(8), intent(in) :: l
    real(8) :: x_old, x_new, tol
    integer :: j, maxit

    ! Parameters
    maxit = 1000
    tol = 1.0e0_8

    ! Initial guess
    x_old = 1.0e4_8

    j = 0
    do
        x_new = x_old - (P(x_old) - l) / deriv(x_old)
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
        e = e_0/24 * (2 * x ** 3 - 3 * x)*(x**2 + 1)**(1/2) + 3 * log(x + sqrt(x**2 + 1.0))
        ! arcsinh = log(x + sqrt(x*x + 1.0))
    end function e

end module TOV_solver
