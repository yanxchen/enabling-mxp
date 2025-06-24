c-----------------------------------------------------------------------
      subroutine write2file(var,size,filename)
      implicit none
      integer size
      real var(size)
      character*(*) filename
      integer :: i
      open(unit=10, file=filename, status='replace', action='write')

      do i=1,size
         write(10,*) var(i)
      enddo
      close(10)
      end
c-----------------------------------------------------------------------
      subroutine write2file_f(var,size,filename)
      implicit none
      integer size
      real*4 var(size)
      character*(*) filename
      integer :: i
      open(unit=10, file=filename, status='replace', action='write')

      do i=1,size
         write(10,*) var(i)
      enddo
      close(10)
      end
c-----------------------------------------------------------------------
