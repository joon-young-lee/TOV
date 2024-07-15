module TOV_solver ! Given EoS
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

    function solve_TOV(file_name, p0, del_r, maxit) result(final_mass_radius_iteration)
        implicit none
        
        character(len=36), intent(in) :: file_name
        real(8), intent(in) :: p0, del_r
        integer, intent(in) :: maxit, k

        real(8) :: final_mass_radius_iteration(3)
        

        integer :: i, j
        real(8) :: p, m, r, r2, r3, r4, p2, p3, p4, m2, m3, m4
        real(8) :: k1_p, k1_m, k2_p, k2_m, k3_p, k3_m, k4_p, k4_m
        

        ! Initial values
        p = p0
        m = 1.0e-10_8
        r = 1.0e-10_8
        

        i = 0
        do while (i <= maxit .and. p_arr(i) > 0.0_8)
            r = r_arr(i)
            p = p_arr(i)
            m = m_arr(i)
            
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

            r4 = r2 + del_r/2.0_8
            p4 = p + del_r * k3_p
            m4 = m + del_r * k3_m
            k4_p = dp_dr(r4, p4, m4)
            k4_m = dm_dr(r4, p4, m4)

            ! Update pressure and mass using the weighted sum of the Runge-Kutta coefficients
            p += (k1_p + 2 * k2_p + 2 * k3_p + k4_p) / 6 * del_r
            m += (k1_m + 2 * k2_m + 2 * k3_m + k4_m) / 6 * del_r
            r += del_r


            ! Increment the iteration counter
            i = i + 1
        end do

        ! Output the final mass and radius
        final_mass_radius_iteration(1) = m/M0
        final_mass_radius_iteration(2) = r/1e5_8
        final_mass_radius_iteration(3) = i



        
    end function solve_TOV


    function eos_arr(file_name, k) result(arrays)
        implicit none
        character(len=*), intent(in) :: file_name
        integer, intent(in) :: k
    
        ! Define a derived type to return multiple arrays
        type :: eos_arrays
            real(8), allocatable :: x(:)
            real(8), allocatable :: y(:)
        end type eos_arrays

        type(eos_arrays) :: arrays

        integer :: ios, num_lines, j
    

        ! Write the filename
        !write(file_name, '(A, I6.6, A)') 'Archive/eft_pnm32_', k, '_ldm_eos.dat'

        ! First pass: Count the number of lines
        open(unit=10, file=file_name, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, "Error opening file:", trim(filename)
            return
        endif

        num_lines = 0
        do
            read(10, '(A)', iostat=ios) line
            if (ios /= 0) exit
            num_lines = num_lines + 1
        end do
        close(10)

        ! Allocate arrays based on the number of lines
        if (num_lines > 0) then
            allocate(arrays%x(num_lines), arrays%y(num_lines))
        else
            print *, "File is empty:", trim(filename)
            return  ! Skip the file if it is empty
        endif

        ! Second pass: Read the data into arrays
        open(unit=10, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, "Error reopening file:", trim(filename)
            return
        endif

        do j = 1, num_lines
            read(10, *) arrays%x(j), arrays%y(j)
        end do
        close(10)

    end function eos_arr


    function evaluate_cubic_spline(p, e, x, coefs, xi)
        ! Evaluates the cubic spline at the given point xi
        implicit none
        real, dimension(:), intent(in) :: x
        real, dimension(:,:), intent(in) :: coefs = compute_cubic_spline_coefficients(p, e)
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

    
    
    real(8) function dp_dr(r, p, m, file_name, k)
        implicit none
        arrays = eos_arr(file_name, k)
        real(8), intent(in) :: r, p, m
        
        dp_dr = - G * e(root(p)) * m * &
        (1 + p / e(root(p))) * &
        (1 + 4 * pi * p * r ** 3/(m * c2)) / &
        (1 - 2 * G * m / (c2 * r)) /&
        (c2 * r ** 2)
    end function dp_dr

    real(8) function dm_dr(file_name, k, r, p, m)
        implicit none
        character(len=*), intent(in) :: file_name
        integer, intent(in) :: k
        real(8), intent(in) :: r, p, m
        
        
        ! TOV equation
        dm_dr = 4 * pi * r**2 * energt_dens / c2
        
    end function dm_dr

end module TOV_solver
