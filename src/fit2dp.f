      subroutine fit2dp(n, x, y, work, ipvt, work2, info, coeff)
      integer n, info, ipvt(n)
      real*8 x(n,3), y(n), work(3,3), work2(3,n), det(2), coeff(3)
      do i = 1, 3
         do j = 1, 3
            work(i,j) = 0
            do k = 1, n
               work(i,j) = work(i,j) + x(k,i)*x(k,j)
            enddo
         enddo
      enddo
      call dgefa(work,3,3,ipvt,info)
      call dgedi(work,3,3,ipvt,det,work2,1)
      if (info .eq. 0) then
         do i = 1, 3
            do j = 1, n
               work2(i,j) = 0
               do k = 1, 3
                  work2(i,j) = work2(i,j) + work(i,k)*x(j,k)
               enddo
            enddo
         enddo
         do i = 1, 3
            coeff(i) = 0
            do j = 1, n
               coeff(i) = coeff(i) + work2(i,j)*y(j)
            enddo
         enddo
      endif
      end
