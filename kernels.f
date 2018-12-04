


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This file contains subroutines for evaluating the various
c       integral kernels associated with the Helmholtz equation.
c
c       These versions of the routines are not accelerated in any 
c       fashion.  See the file fastkernel.f for the accelerated
c       (but somewhat more limited) versions of these routines.
c
c   ksingle - evaluate the kernel of the acoustic single layer operator 
c        S_k with a user-specified wavenumber
c
c   kdouble - evaluate the kernel of the acoustic double layer operator 
c        D_k with a user-specified wavenumber
c
c   kcombined - evaluate the kernel of the acoustic cominbed layer 
c        operator S_k+D_k
c
c   ksingleprime - evaluate the kernel of the derivative of the 
c       acoustic single layer operator S_k'
c
c   kdoubleprime - evaluate the kernel of the derivative of the acoustic
c       double layer operator D_k'
c
c
c
c   ksingle0 - evaluate the kernel of the Laplace single layer operator 
c       S_0
c
c   kdouble0 - eveluate the kernel of the Laplace double layer operator
c       D_0
c
c   ksingle0prime - evaluate the kernel of the derivative of the Laplace
c       single layer operator S_0'
c
c   kdouble0prime - evaluate the kernel of the derivative of the Laplace
c       double layer operator D_0'
c
c
c
c   ksingle_single0 - evaluate the kernel of the operator S_k-S_0
c
c   kdouble_double0 - evaluate the kernel of the operator D_k-D_0
c
c   ksingleprime_single0prime - evaluate the kernel of the operator 
c       (S_k-S_0)'
c
c   kdoubleprime_double0prime - evaluate the kernel of the operator 
c       (D_k-D_0)'
c
c
c
c   ksingle_single - evaluate the kernel of the operator 
c       S_k1 - S_k2 for two different wavenumbers
c
c   kdouble_double - evaluate the kernel of the operator 
c       D_k1 - D_k2 for two different wavenumbers 
c
c   ksingleprime_singleprime - evaluate the kernel of the
c       operator S_k1' - S_k2'
c
c   kdoubleprime_doubleprime - evaluate the kernel of the
c       operator D_k1' - D_k2' 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine ksingle2(x,y,val,rk,par2,par3,par4)
        implicit double precision (a-h,o-z)
        dimension x(10),y(10)
        double complex val,rk,z,h0,h1,ima
        data ima   / (0.0d0,1.0d0) /
        data pi    / 3.14159265358979323846264338327950288d0 /
c
c       Evaluate the kernel of the single layer operator S_k:
c
c           i/4 * H_0(rk*|x-y|).
c
c
        whtx  = x(1)
        x1    = x(2)
        x2    = x(3)
        dx1   = x(4)
        dx2   = x(5)
c
        irx   = x(6)
        rx1   = x(7)
        rx2   = x(8)
c
        whty  = y(1)
        y1    = y(2)
        y2    = y(3)
        dy1   = y(4)
        dy2   = y(5)
c
        iry   = y(6)
        ry1   = y(7)
        ry2   = y(8)
c
        if (irx .ne. iry) then
        x1 = x1+rx1
        x2 = x2+rx2
        y1 = y1+ry1
        y2 = y2+ry2
        endif
c
        r  = sqrt((x1-y1)**2+(x2-y2)**2)
        z  = rk*r
        ifexpon = 1
        call hank103(z,h0,h1,ifexpon)
c
        val = ima/4*h0
        val = val*sqrt(whty)
c
        end


        subroutine kdouble2(x,y,val,rk,par2,par3,par4)
        implicit double precision (a-h,o-z)
        dimension x(10),y(10)
        double complex val,rk,z,h0,h1,ima
        data ima   / (0.0d0,1.0d0) /
        data pi    / 3.14159265358979323846264338327950288d0 /
c
c       Evaluate the kernel of the double layer operator D_k:
c
c                                          (x-y)\cdot n_y
c           i/4(H_1(rk*|x-y|) rk|x-y|)*    --------------
c                                             |x-y|^2
c
c
        whtx  = x(1)
        x1    = x(2)
        x2    = x(3)
        dx1   = x(4)
        dx2   = x(5)
c
        irx   = x(6)
        rx1   = x(7)
        rx2   = x(8)
c
        whty  = y(1)
        y1    = y(2)
        y2    = y(3)
        dy1   = y(4)
        dy2   = y(5)
c
        iry   = y(6)
        ry1   = y(7)
        ry2   = y(8)
c
        if (irx .ne. iry) then
        x1 = x1+rx1
        x2 = x2+rx2
        y1 = y1+ry1
        y2 = y2+ry2
        endif
c
        dry = sqrt(dy1**2+dy2**2)
        r   = sqrt((x1-y1)**2+(x2-y2)**2)
        z   = rk*r
        ifexpon = 1
        call hank103(z,h0,h1,ifexpon)
c
        val = ima/4*h1*rk
        val = val*( (x1-y1)*dy2 - (x2-y2)*dy1)/(r*dry) 
c
c        val = ima/4*h1*rk*z
c        val = val*( (x1-y1)*dy2 - (x2-y2)*dy1)/((r**2)*dry) 

c
        val = val*sqrt(whty)
c        
        end


        subroutine ksingle(x,y,val,rk,par2,par3,par4)
        implicit double precision (a-h,o-z)
        dimension x(10),y(10)
        double complex val,rk,z,h0,h1,ima
        data ima   / (0.0d0,1.0d0) /
        data pi    / 3.14159265358979323846264338327950288d0 /
c
c       Evaluate the kernel of the single layer operator S_k:
c
c           i/4 * H_0(rk*|x-y|).
c
c
        whtx  = x(1)
        x1    = x(2)
        x2    = x(3)
        dx1   = x(4)
        dx2   = x(5)
c
        irx   = x(6)
        rx1   = x(7)
        rx2   = x(8)
c
        whty  = y(1)
        y1    = y(2)
        y2    = y(3)
        dy1   = y(4)
        dy2   = y(5)
c
        iry   = y(6)
        ry1   = y(7)
        ry2   = y(8)
c
        if (irx .ne. iry) then
        x1 = x1+rx1
        x2 = x2+rx2
        y1 = y1+ry1
        y2 = y2+ry2
        endif
c
        r  = sqrt((x1-y1)**2+(x2-y2)**2)
        z  = rk*r
        ifexpon = 1
        call hank103(z,h0,h1,ifexpon)
c
        val = ima/4*h0
        val = val*sqrt(whtx)*sqrt(whty)
c
        end


        subroutine kdouble(x,y,val,rk,par2,par3,par4)
        implicit double precision (a-h,o-z)
        dimension x(10),y(10)
        double complex val,rk,z,h0,h1,ima
        data ima   / (0.0d0,1.0d0) /
        data pi    / 3.14159265358979323846264338327950288d0 /
c
c       Evaluate the kernel of the double layer operator D_k:
c
c                                          (x-y)\cdot n_y
c           i/4(H_1(rk*|x-y|) rk|x-y|)*    --------------
c                                             |x-y|^2
c
c
        whtx  = x(1)
        x1    = x(2)
        x2    = x(3)
        dx1   = x(4)
        dx2   = x(5)
c
        irx   = x(6)
        rx1   = x(7)
        rx2   = x(8)
c
        whty  = y(1)
        y1    = y(2)
        y2    = y(3)
        dy1   = y(4)
        dy2   = y(5)
c
        iry   = y(6)
        ry1   = y(7)
        ry2   = y(8)
c
        if (irx .ne. iry) then
        x1 = x1+rx1
        x2 = x2+rx2
        y1 = y1+ry1
        y2 = y2+ry2
        endif
c
        dry = sqrt(dy1**2+dy2**2)
        r   = sqrt((x1-y1)**2+(x2-y2)**2)
        z   = rk*r
        ifexpon = 1
        call hank103(z,h0,h1,ifexpon)
c
        val = ima/4*h1*rk
        val = val*( (x1-y1)*dy2 - (x2-y2)*dy1)/(r*dry) 
c
c        val = ima/4*h1*rk*z
c        val = val*( (x1-y1)*dy2 - (x2-y2)*dy1)/((r**2)*dry) 

c
        val = val*sqrt(whtx)*sqrt(whty)
c        
        end


        subroutine kcombined(x,y,val,rk,par2,par3,par4)
        implicit double precision (a-h,o-z)
        dimension x(1),y(1)
        double complex val,rk,z,h0,h1,ima
        data ima   / (0.0d0,1.0d0) /
        data pi    / 3.14159265358979323846264338327950288d0 /
c
c       Evaluate the kernel of the combined potential 
c
c               S_k + D_k
c
        whtx  = x(1)
        x1    = x(2)
        x2    = x(3)
        dx1   = x(4)
        dx2   = x(5)
c
        irx   = x(6)
        rx1   = x(7)
        rx2   = x(8)
c
        whty  = y(1)
        y1    = y(2)
        y2    = y(3)
        dy1   = y(4)
        dy2   = y(5)
c
        iry   = y(6)
        ry1   = y(7)
        ry2   = y(8)
c
        if (irx .ne. iry) then
        x1 = x1+rx1
        x2 = x2+rx2
        y1 = y1+ry1
        y2 = y2+ry2
        endif
c
        dry = sqrt(dy1**2+dy2**2)
        r   = sqrt((x1-y1)**2+(x2-y2)**2)
        z  = rk*r
        ifexpon = 1
        call hank103(z,h0,h1,ifexpon)
c
        val = ima/4.0d0*h1*rk*( (x1-y1)*dy2 - (x2-y2)*dy1)/(r*dry)
        val = val+h0
c
        val = val*sqrt(whtx)*sqrt(whty)
c
c
        end


        subroutine kcombined2(x,y,val,rk,par2,par3,par4)
        implicit double precision (a-h,o-z)
        dimension x(1),y(1)
        double complex val,rk,z,h0,h1,ima
        data ima   / (0.0d0,1.0d0) /
        data pi    / 3.14159265358979323846264338327950288d0 /
c
c       Evaluate the kernel of the combined potential 
c
c               S_k + D_k
c
        whtx  = x(1)
        x1    = x(2)
        x2    = x(3)
        dx1   = x(4)
        dx2   = x(5)
c
        irx   = x(6)
        rx1   = x(7)
        rx2   = x(8)
c
        whty  = y(1)
        y1    = y(2)
        y2    = y(3)
        dy1   = y(4)
        dy2   = y(5)
c
        iry   = y(6)
        ry1   = y(7)
        ry2   = y(8)
c
        if (irx .ne. iry) then
        x1 = x1+rx1
        x2 = x2+rx2
        y1 = y1+ry1
        y2 = y2+ry2
        endif
c
        dry = sqrt(dy1**2+dy2**2)
        r   = sqrt((x1-y1)**2+(x2-y2)**2)
        z  = rk*r
        ifexpon = 1
        call hank103(z,h0,h1,ifexpon)
c
        val = ima/4.0d0*h1*rk*( (x1-y1)*dy2 - (x2-y2)*dy1)/(r*dry)
        val = val+h0
c
        val = val*sqrt(whty)
c
c
        end


        subroutine ksingleprime(x,y,val,rk,par2,par3,par4)
        implicit double precision (a-h,o-z)
        dimension x(10),y(10)
        double complex val,rk,z,h0,h1,ima
        data ima   / (0.0d0,1.0d0) /
        data pi    / 3.14159265358979323846264338327950288d0 /
c
c       Evaluate the kernel of the derivative of the single layer
c       operator S'_k:
c
c                                          (y-x)\cdot n_x
c           i/4(H_1(rk*|x-y|) rk|x-y|)*    --------------
c                                             |x-y|^2
c
        whtx  = x(1)
        x1    = x(2)
        x2    = x(3)
        dx1   = x(4)
        dx2   = x(5)
c
        irx   = x(6)
        rx1   = x(7)
        rx2   = x(8)
c
        whty  = y(1)
        y1    = y(2)
        y2    = y(3)
        dy1   = y(4)
        dy2   = y(5)
c
        iry   = y(6)
        ry1   = y(7)
        ry2   = y(8)
c
        if (irx .ne. iry) then
        x1 = x1+rx1
        x2 = x2+rx2
        y1 = y1+ry1
        y2 = y2+ry2
        endif
c
        drx = sqrt(dx1**2+dx2**2)
c
        r  = sqrt((x1-y1)**2+(x2-y2)**2)
        z  = rk*r
        ifexpon = 1
        call hank103(z,h0,h1,ifexpon)
c
        val = ima/4*h1*rk
        val = val*( (y1-x1)*dx2 - (y2-x2)*dx1)/(r*drx)
c
        val = val*sqrt(whtx)*sqrt(whty)
c
        end


        subroutine kdoubleprime(x,y,val,rk,par2,par3,par4)
        implicit double precision (a-h,o-z)
        dimension x(5),y(5)
        double complex val,rk,z,h0,h1,ima
        data ima   / (0.0d0,1.0d0) /
        data pi    / 3.14159265358979323846264338327950288d0 /
c
c       Evaluate the kernel of the derivative of the double layer
c       operator D_k'
c
c
        whtx  = x(1)
        x1    = x(2)
        x2    = x(3)
        dx1   = x(4)
        dx2   = x(5)
c
        irx   = x(6)
        rx1   = x(7)
        rx2   = x(8)
c
        whty  = y(1)
        y1    = y(2)
        y2    = y(3)
        dy1   = y(4)
        dy2   = y(5)
c
        iry   = y(6)
        ry1   = y(7)
        ry2   = y(8)
c
        if (irx .ne. iry) then
        x1 = x1+rx1
        x2 = x2+rx2
        y1 = y1+ry1
        y2 = y2+ry2
        endif
c
c       Compute the values of the Hankel functions.
c
        dd  = sqrt((x1-y1)**2+(x2-y2)**2)
        z   = rk*dd
        ifexpon = 1
        call hank103(z,h0,h1,ifexpon)
c
        val = (ima*rk*((rk*(dx2*(x1 - y1) + dx1*(-x2 + y2))*
     -         (dy2*(x1 - y1) + dy1*(-x2 + y2))*
     -         h0)/
     -       ((x1 - y1)**2 + (x2 - y2)**2) + 
     -      ((dx1*(2*dy2*(x1 - y1)*(x2 - y2) + 
     -              dy1*(x1 + x2 - y1 - y2)*(x1 - x2 - y1 + y2)) + 
     -           dx2*(2*dy1*(x1 - y1)*(x2 - y2) - 
     -              dy2*(x1 + x2 - y1 - y2)*(x1 - x2 - y1 + y2)))*
     -         h1)/
     -       ((x1 - y1)**2 + (x2 - y2)**2)**1.5))/
     -  (4.0d0*Sqrt(dx1**2 + dx2**2)*Sqrt(dy1**2 + dy2**2))

c       Scale by the weights so as to make this an L^2(\Gamma)\to
c       L^2(\Gamma) discretization.
c
        val = val*sqrt(whty)*sqrt(whtx)
c                
        end



        subroutine ksingle0(x,y,val,par1,par2,par3,par4)
        implicit double precision (a-h,o-z)
        dimension x(5),y(5)
        double complex val
        data pi    / 3.14159265358979323846264338327950288d0 /
        data pi4   /12.56637061435917295385057353311801154d0 /
c
c       Evaluate the Laplace single layer kernel:
c
c             1
c        -  ----  log |x-y| 
c           2*pi
c
c
        whtx  = x(1)
        x1    = x(2)
        x2    = x(3)
        dx1   = x(4)
        dx2   = x(5)
c
        irx   = x(6)
        rx1   = x(7)
        rx2   = x(8)
c
        whty  = y(1)
        y1    = y(2)
        y2    = y(3)
        dy1   = y(4)
        dy2   = y(5)
c
        iry   = y(6)
        ry1   = y(7)
        ry2   = y(8)
c
        if (irx .ne. iry) then
        x1 = x1+rx1
        x2 = x2+rx2
        y1 = y1+ry1
        y2 = y2+ry2
        endif
c
        r2   = (x1-y1)**2+(x2-y2)**2
        val0 = log(r2)/(pi4)
        val0 = -val0*sqrt(whtx)*sqrt(whty)
c
        val = val0
ccccccccccccccccccccccccccccccccc
c
        end


        subroutine kdouble0(x,y,val,rk,par2,par3,par4)
        implicit double precision (a-h,o-z)
        dimension x(5),y(5)
        double complex val
        data ima   / (0.0d0,1.0d0) /
        data pi    / 3.14159265358979323846264338327950288d0 /
c
c       Evaluate the kernel of the double layer operator D_0:
c
        whtx  = x(1)
        x1    = x(2)
        x2    = x(3)
        dx1   = x(4)
        dx2   = x(5)
c
        irx   = x(6)
        rx1   = x(7)
        rx2   = x(8)
c
        whty  = y(1)
        y1    = y(2)
        y2    = y(3)
        dy1   = y(4)
        dy2   = y(5)
c
        iry   = y(6)
        ry1   = y(7)
        ry2   = y(8)
c
        if (irx .ne. iry) then
        x1 = x1+rx1
        x2 = x2+rx2
        y1 = y1+ry1
        y2 = y2+ry2
        endif
c
        dry = sqrt(dy1**2+dy2**2)
        dd  = (x1-y1)**2+(x2-y2)**2
c
        val0 = ( (x1-y1)*dy2 - (x2-y2)*dy1)/(dd*dry) 
        val0 = val0*sqrt(whtx)*sqrt(whty)
c
        val = val0/(2*pi)
c        
        end


        subroutine kdouble20(x,y,val,rk,par2,par3,par4)
        implicit double precision (a-h,o-z)
        dimension x(5),y(5)
        double complex val
        data ima   / (0.0d0,1.0d0) /
        data pi    / 3.14159265358979323846264338327950288d0 /
c
c       Evaluate the kernel of the double layer operator D_0:
c
        whtx  = x(1)
        x1    = x(2)
        x2    = x(3)
        dx1   = x(4)
        dx2   = x(5)
c
        irx   = x(6)
        rx1   = x(7)
        rx2   = x(8)
c
        whty  = y(1)
        y1    = y(2)
        y2    = y(3)
        dy1   = y(4)
        dy2   = y(5)
c
        iry   = y(6)
        ry1   = y(7)
        ry2   = y(8)
c
        if (irx .ne. iry) then
        x1 = x1+rx1
        x2 = x2+rx2
        y1 = y1+ry1
        y2 = y2+ry2
        endif
c
        dry = sqrt(dy1**2+dy2**2)
        dd  = (x1-y1)**2+(x2-y2)**2
c
        val0 = ( (x1-y1)*dy2 - (x2-y2)*dy1)/(dd*dry) 
        val0 = val0*sqrt(whty)
c
        val = val0/(2*pi)
c        
        end


        subroutine ksingle0prime(x,y,val,rk,par2,par3,par4)
        implicit double precision (a-h,o-z)
        dimension x(5),y(5)
        double complex val
        data ima   / (0.0d0,1.0d0) /
        data pi    / 3.14159265358979323846264338327950288d0 /
c
c       Evaluate the kernel of the double layer operator S_0':
c
        whtx  = x(1)
        x1    = x(2)
        x2    = x(3)
        dx1   = x(4)
        dx2   = x(5)
c
        irx   = x(6)
        rx1   = x(7)
        rx2   = x(8)
c
        whty  = y(1)
        y1    = y(2)
        y2    = y(3)
        dy1   = y(4)
        dy2   = y(5)
c
        iry   = y(6)
        ry1   = y(7)
        ry2   = y(8)
c
        if (irx .ne. iry) then
        x1 = x1+rx1
        x2 = x2+rx2
        y1 = y1+ry1
        y2 = y2+ry2
        endif
c
        drx = sqrt(dx1**2+dx2**2)
        dd  = (x1-y1)**2+(x2-y2)**2
c
        val0 = ( (y1-x1)*dx2 - (y2-x2)*dx1)/(dd*drx)
        val0 = val0/(2*pi)
        val0 = val0*sqrt(whtx)*sqrt(whty)
c
        val = val0
c        
        end




        subroutine kdoubleprime_double0prime(x,y,val,rk,par2,
     1    par3,par4)
        implicit double precision (a-h,o-z)
        dimension x(10),y(10)
        double complex val,ima,h0,h1,z,omega,f,g,df,dg
        data ima   / (0.0d0,1.0d0) /
        data pi /3.14159265358979323846264338327950288d0/
c
c
        whtx  = x(1)
        x1    = x(2)
        x2    = x(3)
        dx1   = x(4)
        dx2   = x(5)
c
        irx   = x(6)
        rx1   = x(7)
        rx2   = x(8)
c
        whty  = y(1)
        y1    = y(2)
        y2    = y(3)
        dy1   = y(4)
        dy2   = y(5)
c
        iry   = y(6)
        ry1   = y(7)
        ry2   = y(8)
c
        if (irx .ne. iry) then
        x1 = x1+rx1
        x2 = x2+rx2
        y1 = y1+ry1
        y2 = y2+ry2
        endif
c
c       Compute the values of the Hankel functions.
c
        r   = sqrt((x1-y1)**2+(x2-y2)**2)
        z   = rk*r
        ifexpon = 1
        call hank103(z,h0,h1,ifexpon)
c
        drx = sqrt(dx1**2+dx2**2)
        dry = sqrt(dy1**2+dy2**2)
c
        dx1 = dx1 / drx
        dx2 = dx2 / drx
        dy1 = dy1 / dry
        dy2 = dy2 / dry        
c
c       Compute f, g grad(f), grad(g)
c
        f = ((x1-y1)*dy2-(x2-y2)*dy1)/(r**2)
        f = f/(2*pi)
c
c       Compute g-1
c
        if (abs(z) .le. 1.0d-3) then
        z=rk*r
        call hank1(z,g)
        else
        g = pi*ima/2.0d0 * h1*rk*r-1.0d0
        endif
c
        dg = rk**2*pi*ima/2.0d0*( (x1-y1)*dx2-(x2-y2)*dx1)*h0
c
        df = (dx1*dy1+dx2*dy2)/r**2
        df = df-2*((x1-y1)*dx2-(x2-y2)*dx1)*((x1-y1)*dy2-(x2-y2)*dy1)/
     1    r**4
c
        df = df/(2*pi)
        val = (g)*df+dg*f
c
c       Scale by the weights so as to make this an L^2(\Gamma)\to
c       L^2(\Gamma) discretization.
c
        val = val*sqrt(whty)*sqrt(whtx)
c
        end


        subroutine hank1(z,val)
        implicit double precision (a-h,o-z)
        double complex z,val,val0
        double complex j1
c
c       Evaluate the function 
c
c              Pi*I/2 H_1(z)*z-1
c
c       for z near 0.
c
c
c
c       First compute J1.
c
        j1=z/2.0d0 - z**3/16.0d0 + z**5/384.0d0 - z**7/18432.0d0 + 
     1     z**9/1.47456d6
c
        val0 = -j1*z*log(z/2.0d0)
c
c       Now compute the analytic part
c
        val = (-0.038607832450766430303256045041201216d0,
     -    0.785398163397448309615660845819875721d0)*z**2 - 
     -  (0.0420490209436541962120929943698498481d0,
     -    0.0981747704246810387019576057274844651d0)*z**4 + 
     -  (0.00283711198376336928661498587652152145d0,
     -    0.00409061543436170994591490023864518605d0)*z**6 - 
     -  (0.0000749304290598850082859603539090131783d0,
     -    0.0000852211548825356238732270883051080426d0)*z**8 + 
     -  (1.08921825387356260357450442386266473d-6,
     -    1.06526443603169529841533860381385053d-6)*z**10
c
        val = val0+val
        end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       The egregious laziness below needs to be fixed.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine ksingle_single(x,y,val,rk1,rk2,par3,par4)
        implicit double precision (a-h,o-z)
        dimension x(1),y(1)
        double complex val1,val2,val
c
        call ksingle(x,y,val1,rk1,par2,par3,par4)
        call ksingle(x,y,val2,rk2,par2,par3,par4)
c
        val = val1-val2
        end


        subroutine kdouble_double(x,y,val,rk1,rk2,par3,par4)
        implicit double precision (a-h,o-z)
        dimension x(1),y(1)
        double complex val1,val2,val
c
        call kdouble(x,y,val1,rk1,par2,par3,par4)
        call kdouble(x,y,val2,rk2,par2,par3,par4)
c
        val = val1-val2
        end


        subroutine ksingleprime_singleprime(x,y,val,rk1,rk2,par3,par4)
        implicit double precision (a-h,o-z)
        dimension x(1),y(1)
        double complex val1,val2,val
c
        call ksingleprime(x,y,val1,rk1,par2,par3,par4)
        call ksingleprime(x,y,val2,rk2,par2,par3,par4)
c
        val = val1-val2
        end


        subroutine kdoubleprime_doubleprime(x,y,val,rk1,rk2,par3,par4)
        implicit double precision (a-h,o-z)
        dimension x(1),y(1)
        double complex val1,val2,val
c
        call kdoubleprime(x,y,val1,rk1,par2,par3,par4)
        call kdoubleprime(x,y,val2,rk2,par2,par3,par4)
c
        val = val1-val2
        end


