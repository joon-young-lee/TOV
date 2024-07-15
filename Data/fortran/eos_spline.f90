module cubic_spline
    implicit none
    private
    public :: compute_cubic_spline_coefficients, evaluate_cubic_spline

contains

    function compute_cubic_spline_coefficients(x, fx, coefs)
        ! Computes the coefficients for cubic spline interpolation
        implicit none
        real, dimension(:), intent(in) :: x, fx
        real, dimension(:,:), allocatable, intent(out) :: coefs
        integer :: n, i
        real, dimension(:), allocatable :: h, b, u, v, z

        n = size(x)
        allocate(coefs(n-1, 4))  ! Each interval has 4 coefficients
        allocate(h(n-1), b(n-1), u(n), v(n), z(n))

        ! Step 1: Calculate h and b
        do i = 1, n-1
            h(i) = x(i+1) - x(i)
            b(i) = (fx(i+1) - fx(i)) / h(i)
        end do

        ! Step 2: Solve the tridiagonal system
        u(1) = 2.0 * (h(1) + h(2))
        v(1) = 6.0 * (b(2) - b(1))
        do i = 2, n-2
            u(i) = 2.0 * (h(i) + h(i+1)) - h(i-1)**2 / u(i-1)
            v(i) = 6.0 * (b(i+1) - b(i)) - h(i-1) * v(i-1) / u(i-1)
        end do

        ! Back substitution
        z(n-1) = v(n-2) / u(n-2)
        do i = n-2, 2, -1
            z(i) = (v(i-1) - h(i) * z(i+1)) / u(i-1)
        end do
        z(1) = 0.0
        z(n) = 0.0

        ! Step 3: Compute the spline coefficients
        do i = 1, n-1
            coefs(i, 1) = fx(i)
            coefs(i, 2) = b(i) - h(i) * (2.0 * z(i) + z(i+1)) / 6.0
            coefs(i, 3) = z(i) / 2.0
            coefs(i, 4) = (z(i+1) - z(i)) / (6.0 * h(i))
        end do

        deallocate(h, b, u, v, z)
    end subroutine compute_cubic_spline_coefficients

    real function evaluate_cubic_spline(x, coefs, xi)
        ! Evaluates the cubic spline at the given point xi
        implicit none
        real, dimension(:), intent(in) :: x
        real, dimension(:,:), intent(in) :: coefs
        real, intent(in) :: xi
        integer :: i, n

        n = size(x)
        evaluate_cubic_spline = 0.0

        ! Find the right interval
        do i = 1, n-1
            if (xi >= x(i) .and. xi <= x(i+1)) then
                evaluate_cubic_spline = coefs(i, 1) + &
                    coefs(i, 2) * (xi - x(i)) + &
                    coefs(i, 3) * (xi - x(i))**2 + &
                    coefs(i, 4) * (xi - x(i))**3
                exit
            end if
        end do
    end function evaluate_cubic_spline

end module cubic_spline

