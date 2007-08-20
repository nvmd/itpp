c     zdotusub.f
c
c     The program is a fortran wrapper for zdotu.
c     Witten by Keita Teranishi.  2/11/1998
c     Interface modified by Adam Piatyszek to match ACML zdotusub_()
c
      subroutine zdotusub(dotu,n,x,incx,y,incy)
c
      external zdotu
      double complex zdotu,dotu
      integer n,incx,incy
      double complex x(*),y(*)
c
      dotu=zdotu(n,x,incx,y,incy)
      return
      end
