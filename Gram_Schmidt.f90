module basic_types

    implicit none
	
    private
	
    public :: int4, real8
	
    ! define basic types
    integer, parameter :: int4	= selected_int_kind(8)
    integer, parameter :: real8 = selected_real_kind(p=15)
	
end module basic_types
	
module io_files
use basic_types
	
	implicit none
	
	private
	
	public :: output_unit
	
	! define basec units
	integer(int4), parameter :: output_unit = 6 ! sreen
	
end module io_files
	
program Gram_Schmidt
use io_files
use basic_types

    implicit none
    
	integer(int4) :: i, j, k, m, n
	real(real8), dimension(:,:), allocatable :: A, B, Q, R
    real(real8), dimension(:), allocatable :: u  
    
    m = 3; n = 3 ! Size of the sistem
	
	allocate(A(m,n),Q(m,n),R(m,n))
	allocate(u(m))
	
	! Initialize the matrix A
    !A = reshape((/ 1.0,  1.0, &
    !               3.0, -1.0  /), shape(A))
	
	A = reshape((/ 1.0, 1.0, 0.0, &
                   1.0, 0.0, 1.0, &
                   0.0, 1.0, 1.0 /), shape(A))
	
	! Perform Gram-Schmidt process
    do k = 1, n
        u = A(:, k) ! u <- A
        do j = 1, k-1
            R(j, k) = dot_product(Q(:, j), A(:, k))
            u = u - R(j, k) * Q(:, j) ! a2 - <a2 . e1> e1
        end do
        R(k, k) = sqrt(dot_product(u,u)) 
        Q(:, k) = u / R(k, k)
    end do
	
	write(output_unit,'(/5X,A16/)') "The matrix Q is:"
	do i = 1, n		
	    write(output_unit,'(5X,6F8.4)') Q(i,:)		
	end do
	
	write(output_unit,*)
	
	write(output_unit,'(/5X,A16/)') "The matrix R is:"
	do i = 1, n		
	    write(output_unit,'(5X,6F8.4)') R(i,:)		
	end do
	
	write(output_unit,*)
	
	write(output_unit,'(/5X,A21/)') "The matrix A = QR is:"
	B = matmul(Q,R)
	do i = 1, n		
	    write(output_unit,'(5X,6F8.4)') B(i,:)		
	end do 	
	
	write(output_unit,*)
	
end program Gram_Schmidt