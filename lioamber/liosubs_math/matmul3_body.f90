!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!function matmul3_xxx( Amat, Bmat, Cmat ) result( Dmat )
!  implicit none
!  xxxxxxxxxx, intent(in) :: Amat(:,:)
!  xxxxxxxxxx, intent(in) :: Bmat(:,:)
!  xxxxxxxxxx, intent(in) :: Cmat(:,:)
!  xxxxxxxxxx             :: Dmat(:,:)
   logical :: error_found
!------------------------------------------------------------------------------!
   error_found = .false.
   error_found = (error_found) .or. ( size(Amat,2) /= size(Bmat,1) )
   error_found = (error_found) .or. ( size(Cmat,1) /= size(Bmat,2) )

   if (error_found) then
      print*, 'ERROR INSIDE matmul3_xxx'
      print*, 'Wrong sizes of input/output'
      print*; stop
   endif

!------------------------------------------------------------------------------!
  allocate(Emat(size(Bmat,1),size(Bmat,2)))
	write(*,*) "flag 1 matmul3"
  Emat = matmul( Amat, Bmat )
	write(*,*) "flag 2 matmul3"
  Dmat = matmul( Emat, Cmat )
  deallocate(Emat)
!a  Dmat=Emat
	write(*,*) "flag 3 matmul3"

!end function matmul3_ddd
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
