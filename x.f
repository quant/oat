	program x
	real, allocatable :: a(:)
	call mtrace
!dir$ attributes c :: mtrace
	do i=1,1000000
	  allocate(a(i))
	  if (mod(loc(a),16).eq.15) print*,loc(a)
	  deallocate(a)
	enddo
	end
