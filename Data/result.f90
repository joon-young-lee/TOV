function eos
    ! input (data of eos, del_r, maxit, p_i, p_f, step)
    ! output Mass, Radius
    implicit none

    ! Parameters
    integer, parameter :: num_files = 100

    ! Variables
    integer :: i, j, ios, num_lines
    character(len=36) :: filename
    real(8), allocatable :: x(:), y(:)
    character(len=100) :: line

    ! Loop over each file
    do i = 1, num_files
        write(filename, '(A, I6.6, A)') 'Archive/eft_pnm32_', i, '_ldm_eos.dat'

        ! First pass: Count the number of lines
        open(unit=10, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, "Error opening file:", filename
            cycle
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
            allocate(x(num_lines), y(num_lines))
        else
            cycle  ! Skip the file if it is empty
        endif

        ! Second pass: Read the data into arrays
        open(unit=10, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, "Error reopening file:", filename
            cycle
        endif

        do j = 1, num_lines
            read(10, *) x(num_lines - j + 1), y(num_lines - j + 1)
        end do
        close(10)

        
        print *, "File:", filename, "Points read:", num_lines
        

        ! Deallocate arrays to prepare for the next file
        deallocate(x, y)
    end do

end function eos

