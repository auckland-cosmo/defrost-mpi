      subroutine p3dfft_ftran_r2c(IN,OUT)
      
      use p3dfft

      real(mytype) IN(1,1,*)
      complex(mytype) OUT(1,1,*)
c      real(mytype) IN(nx_fft,jistart:jiend,kjstart:kjend)
c      complex(mytype) OUT(iistart:iiend,jjstart:jjend,nz_fft)
c      logical flg_inplace

      call ftran_r2c(IN,OUT)

      return
      end

      subroutine p3dfft_btran_c2r(IN,OUT)
      
      use p3dfft

      real(mytype) OUT(1,1,*)
      complex(mytype) IN(1,1,*)
c      real(mytype) OUT(nx_fft,jistart:jiend,kjstart:kjend)
c      complex(mytype) IN(iistart:iiend,jjstart:jjend,nz_fft)
c      logical flg_inplace

      call btran_c2r(IN,OUT)

      return
      end
