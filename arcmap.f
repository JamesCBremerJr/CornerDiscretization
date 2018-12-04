        implicit double precision (a-h,o-z)
c
        dimension ss(1 000 000),swhts(1 000 000)
        dimension ts(1 000 000),twhts(1 000 000)
c
        dimension w(10 000 000)
c
        external ellipse,snowcone,fundb
c
        pi     = acos(-1.0d0)        
c
c       Find the arclength of the ellipse.
c
        a  = -pi
        b  =  pi
        c  = a
c
        call arclen(a,b,ellipse,par1,par2,par3,par4,rl)
        call prin2("rl=*",rl,1)
c
c       Construct a quadrature in the arclength variable.
c
        epslen = .01d0
        nlege  = 30
        call mesh(epslen,0.0d0,rl,nquad,ss,swhts)
c
        call prinf("nquad=*",nquad,1)
c
c       Map the quadrature formula.
c
        call elapsed(t1)
        call arcmap(ier,c,nquad,ss,swhts,ts,twhts,w,
     1     ellipse,par1,par2,par3,par4)
        call elapsed(t2)
c        
       call prin2("ss=*",ss,nquad)
c        call prin2("ts=*",ts,nquad)
c
        call prin2("arcmap time = *",t2-t1,1)
c
c       Test the formula by integrating a double layer potential --- the
c       result should be 0.
c
        x0 = 10.0d0
        y0 = -3.0d0
c
        sum = 0

        do 1000 j=1,nquad
        t   = ts(j)
        wht = twhts(j)        
        call fundb(x0,y0,t,val,ellipse,par1,par2,par3,par4)
        sum = sum+val*wht
 1000 continue
c
        errabs = abs(sum)
        call prin2("ellipse errabs = *",errabs,1)
c
c       Find the extents of the snowcone curve.
c
        theta=pi/2
c
        a = -pi
        b = pi
        z = 0
c
        call arclen(z,a,snowcone,theta,par2,par3,par4,as)
        call arclen(z,b,snowcone,theta,par2,par3,par4,bs)
c
        call prin2("as=*",as,1)
        call prin2("bs=*",bs,1)
c
c       Build a graded quadrature on the snowcone.
c
        delta  = .1d0
        dsub   = 5
        epslen = .1d0
        epscut = 1.0d-30
c
        nlege  = 16
c
        call gradedmesh(as,bs,nlege,epslen,epscut,dsub,delta,
     1    nquad,ss,swhts)
c
c       Map the formula.
c
        c = 0
c
        call elapsed(t1)
        call arcmap(ier,c,nquad,ss,swhts,ts,twhts,w,
     1     snowcone,theta,par2,par3,par4)
        call elapsed(t2)
c
        call prinf("nquad=*",nquad,1)
        call prin2("ss=*",ss,nquad)
        call prin2("ts=*",ts,nquad)
c
        call prin2("arcmap time = *",t2-t1,1)
c
c       Test the formula by integrating a double layer potential --- the
c       result should be 0.
c
        x0 = 3.0d0
        y0 = 0.0d0
c
        sum = 0

        do 2000 j=1,nquad
        t   = ts(j)
        wht = twhts(j)
        call fundb(x0,y0,t,val,snowcone,theta,par2,par3,par4)
        sum = sum+val*wht
 2000 continue
c
        errabs = abs(sum)
        call prin2("snowcone errabs = *",errabs,1)
c
        end
c
c
c
        subroutine fundb(x0,y0,t,val,curve,par1,par2,par3,par4)
        implicit double precision (a-h,o-z)
        external curve        
c
        call curve(t,x,y,dx,dy,par1,par2,par3,par4)
c
        pi  = acos(-1.0d0)
c
        val = ( (x0-x)*dy-(y0-y)*dx ) / ((x0-x)**2+(y0-y)**2)
        val = val/(2*pi)
c
        end



        subroutine mesh(epslen,a,b,n,xs,whts)
        implicit double precision (a-h,o-z)
        dimension xs(1),whts(1),ab(2,10 000)
        dimension xslege(1000),whtslege(1000)
c
c       Fetch a Legendre quadrature.
c
        nlege = 30
        call legequad(nlege,xslege,whtslege)
c
c       Determine how long each interval should really be.
c
        dd = (b-a)/epslen
        nints = ceiling(dd)
c
        dd = nints
        dd = (b-a)/dd
c
c       Construct the list of intervals.
c
        t=a

        do 1000 j=1,nints
        ab(1,j)=t
        ab(2,j)=t+dd
        t=t+dd
 1000 continue
c
c       Construct the quadrature.
c
        n=0
c
        do 1100 j=1,nints
        aa = ab(1,j)
        bb = ab(2,j)
c
        alpha = (bb-aa)/2
        beta  = (bb+aa)/2
c
        do 1200 i=1,nlege
        n=n+1
c
        xs(n)=xslege(i)*alpha+beta
        whts(n)=whtslege(i)*alpha
c
 1200 continue
 1100 continue
c
        end



        subroutine gradedmesh(a,b,nlege,epslen,epscut,dsub,delta,
     1    nquad,xs,whts)
        implicit double precision (a-h,o-z)
        dimension xs(1),whts(1)
c
c       Automatically allocated variables.
c
        dimension xslege(nlege),whtslege(nlege)
        dimension ab(2,10 000)
c
c       Build a graded mesh on the interval [a,b].
c
        call legequad(nlege,xslege,whtslege)
c
        nquad=0
c
        if (-delta .ne. a) then
c      
c       Build the portion of the quadrature on [a,-delta].
c
        dd    = (-delta-a)/epslen
        nints = ceiling(dd)
        dd    = nints
        dd    = (-delta-a)/dd
c
        do 2000 j=1,nints
        aa = a+dd*(j-1)
        bb = a+dd*j
c
        if (j .eq. nints) bb = -delta
c
        alpha = (bb-aa)/2
        beta  = (bb+aa)/2
c
        do 2100 i=1,nlege
        x=xslege(i)*alpha+beta
        wht=whtslege(i)*alpha
        nquad=nquad+1
        xs(nquad)=x
        whts(nquad)=wht
 2100 continue
 2000 continue
c
        endif
c
c       Build the quadrature on [-delta,0]
c
        dd = -log((epscut/delta)**(1.0d0/dsub))/log(2.0d0)
        nn = ceiling(dd)
c
        do 1000 j=1,nn
c
        aa = -2.0d0**(-j+1)
        bb = -2.0d0**(-j)
c
        alpha = (bb-aa)/2
        beta  = (bb+aa)/2
c
        do 1100 i=1,nlege
        x           = xslege(i)*alpha+beta
        wht         = whtslege(i)*alpha
c
        nquad=nquad+1
        xs(nquad)   = delta*x**dsub
        whts(nquad) = delta*wht*x**(dsub-1)*dsub
 1100 continue
 1000 continue
c
c       Build the quadrature on [0,delta]
c
        do 1200 j=nn,1,-1
c
        aa = 2.0d0**(-j)
        bb = 2.0d0**(-j+1)
c
        alpha = (bb-aa)/2
        beta  = (bb+aa)/2
c
        do 1300 i=1,nlege
        x           = xslege(i)*alpha+beta
        wht         = whtslege(i)*alpha
c
        nquad=nquad+1
        xs(nquad)   = delta*x**dsub
        whts(nquad) = delta*wht*x**(dsub-1)*dsub
 1300 continue
 1200 continue
c
c
c       Build the portion of the quadrature on [delta,b].
c
        if (delta .ne. b) then
c
        dd    = (b-delta)/epslen
        nints = ceiling(dd)
        dd    = nints
        dd    = (b-delta)/dd
c
        do 2200 j=1,nints
        aa = delta+dd*(j-1)
        bb = delta+dd*j
c
        if (j .eq. nints) bb = b
c
        alpha = (bb-aa)/2
        beta  = (bb+aa)/2
c
        do 2300 i=1,nlege
        x=xslege(i)*alpha+beta
        wht=whtslege(i)*alpha
        nquad=nquad+1
        xs(nquad)=x
        whts(nquad)=wht
 2300 continue
 2200 continue
c
        endif
c
        end


        subroutine ellipse(t,x,y,dx,dy,par1,par2,par3,par4)
        implicit double precision (a-h,o-z)
        data pi /3.14159265358979323846264338327950288d0/
c
        a = 2.0d0
        b = 1.0d0
c
        x = a*cos(t)
        y = b*sin(t)
c
        dx = -a*sin(t)
        dy = b*cos(t)
c
        end


        subroutine snowcone(t,x,y,dx,dy,theta,par2,par3,par4)
        implicit double precision (a-h,o-z)
c
        beta = tan(theta/2)
c
        if (t .lt. 0) then
        x=-2*sin(t/2)
        y=-beta*sin(t)
c
        dx = -cos(t/2)
        dy =-beta*cos(t)
c
        else
        x= 2*sin(t/2)
        y=-beta*sin(t)
c
        dx = cos(t/2)
        dy =-beta*cos(t)
c
        endif
        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging code and beginning of the
c       code proper.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This file contains routines for constructing parameterization
c       independent quadratures on planar curve segments.  Given
c       a user-supplied parameterization r(t) on an interval [a,b],
c       the subroutine arcmap maps a quadrature formula given in the 
c       arclength distance variable
c
c                     t
c          s(t) =  int  |r'(u)| du,
c                     c
c       
c       where c is a user-specified point in the interval [a,b], into 
c       the original parameterization variable t.  
c
c       That is, given a quadrature rule {s_j,w_j} such that
c
c          \int f(x) ds(x) \approx  \sum  f(\gamma(s_j)) w_j            (1)
c
c       where \gamma is the parameterization of the curve segment in 
c       terms of the arclength variable s, arcmap will construct a 
c       quadrature rule {t_j,v_j} such that
c
c          \int f(x) ds(x) \approx \sum f(r(t_j)) |r'(t_j)| v_j.        (2)
c
c       The following subroutines are user-callable:
c
c    arclen - return the arclength of a planar curve segment
c
c    arcradius - find the intersections of a small segment of a planar
c       curve with a ball of given radius centered at a point on the
c       curve
c
c    arcmap - map a quadrature formula given in the arclength variable
c       s into a quadrature formula given in the original user-suppied
c       parameterization variable t
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



        subroutine arclen(a,b,curve,par1,par2,par3,par4,rl)
        implicit double precision (a-h,o-z)
        external curve
c
c       Return the arclength of a curve segment given a user-supplied 
c       parameterization of that curve.
c
c                          Input Parameters:
c
c   (a,b) - the parameterization interval
c   curve - an external subroutine specifying a parameterization of the
c       curve; the calling sequence for curve is as follows:
c
c       subroutine curve(t,x,y,dx,dy,par1,par2,par3,par4)
c
c       Return the coordinates (x(t),y(t)) and their derivatives with
c       respect to t.  The parameters par? are user-supplied and of
c       arbitrary type.
c
c                          Output Parameters:
c
c   rl - the arclength of the curve segment
c
        call get_arclength(a,b,rl,curve,par1,par2,par3,par4)
c
        end



        subroutine arcradius(ier,epsrad,curve,par1,par2,par3,par4,
     1     a,b)
        implicit double precision (a-h,o-z)
        external curve
c
c       Given a user-supplied parameterization r(t) defined over an
c       interval containing the point 0, this subroutine searches for
c       a and b such that
c
c                       |r(a)| = epsrad = |r(b)|.
c
c       If no such points exist, or there are more than two 
c       intersections, then this subroutine can fail spectactulary,
c       including returning false results. 
c
c                            Input Parameters:
c
c   eps - radius of the 
c   curve - an external subroutine specifying a parameterization of the
c       curve; the calling sequence for curve is as follows:
c
c       subroutine curve(t,x,y,dx,dy,par1,par2,par3,par4)
c
c       Return the coordinates (x(t),y(t)) and their derivatives with
c       respect to t.  The parameters par? are user-supplied and of
c       arbitrary type.
c
c   par? - arbitrarily typed user supplied parameters which are simply
c       passed on to the 
c
c                           Output Parameters:
c
c   ier - an error return code;
c       ier = 0    indicates successful execution
c       ier = 32   means that the Newton procedure failed to find th
c                 
c       ier = 64   means that the Newton procedure
c
c
        ier  = 0
        eps0 = 1.0d-15
c
c       Perform Newton iterations with t=epsrad as the initial guess.
c
        maxiters = 7
        maxsteps = 4
c
        t0 = epsrad/10.0d0
        call curve(t0,x,y,dx,dy,par1,par2,par3,par4)
        r0   = sqrt(x**2+y**2)
        dr0  = (x*dx+y*dy)/r0 
        df0  = r0-epsrad
c
        do 1000 iter=1,maxiters
c
        if (abs(df0/dr0) .lt. eps0) then
        b = t
        goto 1200
        endif
c
        do 1100 istep=1,maxsteps
        alpha   = 2.0d0**(istep-1)
        t       = t0-alpha*(df0/dr0)
c
        if (t .lt. 0) goto 1100
c
        call curve(t,x,y,dx,dy,par1,par2,par3,par4)
        r       = sqrt(x**2+y**2)
        dr      = (x*dx+y*dy)/r
        df      = r-epsrad
c
c       Accept the step if the error is sufficiently improved.
c
        if (abs(df) .lt. abs(df0)*.9d0) then
        r0   = r
        t0   = t
        dr0  = dr
        df0  = r0-epsrad
        goto 1000
        endif
 1100 continue
c
c       The maximum number of steps has been exceeded
c
        ier = 32
        return
c        
 1000 continue
c
c       The maximum number of iterations has been exceeded.
c
        ier = 64
        return
 1200 continue
c
        t0 = -epsrad/10.0d0
        call curve(t0,x,y,dx,dy,par1,par2,par3,par4)
        r0   = sqrt(x**2+y**2)
        dr0  = (x*dx+y*dy)/r0 
        df0  = r0-epsrad
c
        do 2000 iter=1,maxiters
c
        if (abs(df0/dr0) .lt. eps0) then
        a = t
        goto 2200
        endif
c
        do 2100 istep=1,maxsteps
        alpha   = 2.0d0**(istep-1)
        t       = t0-alpha*(df0/dr0)
c
        if (t .gt. 0) goto 2100
c
        call curve(t,x,y,dx,dy,par1,par2,par3,par4)
        r       = sqrt(x**2+y**2)
        dr      = (x*dx+y*dy)/r
c
c       Accept the step if the error is sufficiently improved.
c
        if (abs(df) .lt. abs(df0)*.9d0) then
        r0   = r
        t0   = t
        dr0  = dr
        df0  = r0-epsrad
        goto 2000
        endif
 2100 continue
c
c       The maximum number of steps has been exceeded
c
        ier = 32
        return
c        
 2000 continue
c
c       The maximum number of iterations has been exceeded.
c
        ier = 64
        return
 2200 continue
c
        end


        subroutine arcmap(ier,c,nquad,ss,swhts,ts,twhts,w,
     1     curve,par1,par2,par3,par4)
        implicit double precision (a-h,o-z)
        dimension ss(1),swhts(1),ts(1),twhts(1),w(nquad,1)
        external curve
c     
c       This subroutine takes as input a user-specified parameterization 
c       r(t) given over an interval [a,b] and a quadrature formula
c       specified in the arclength variable
c
c                    t
c         s(t) = \int  |r'(u)|du
c                    c
c
c       and maps the quadrature formula into the original 
c       parameterization variable t.
c
c                         Input Parameters:
c
c   c - point in the interval [a,b] which specified what point on the
c       curve with respect to which arclength distances will be
c       measures
c   nquad - the number of points in the input quadrature
c   ss - the nodes of the quadrature in the arc length coordinate s
c   swhts - the weights of the quadrature in the arclength coordinate
c
c   w - a work array which must be of length at least 3*nquad
c
c   curve - an external subroutine specifying a parameterization of the
c       curve; the calling sequence for curve is as follows:
c
c       subroutine curve(t,x,y,dx,dy,par1,par2,par3,par4)
c
c       Return the coordinates (x(t),y(t)) and their derivatives with
c       respect to t.  The parameters par? are user-supplied and of
c       arbitrary type.
c
c                        Output Parameters:
c
c   ier - an error return code;
c       ier = 0    indicates successful execution
c       ier = 16   means that the Newton procedure failed in 
c                  convergence in a reasonable number of iterations
c       ier = 128  means that the Newton step-length control procedure
c                  failed
c   ts - the nodes of the mapped quadrature formula
c   twhts - the weights of the mapped quadrature formula
c
c
        ier = 0
c 
c       Call the auxillary routine to shape the work array.
c    
        call arcmap0(ier,a,b,c,nquad,ss,swhts,ts,twhts,
     1     w(1,1),w(1,2),w(1,3),curve,par1,par2,par3,par4)
c
        end


        subroutine arcmap0(ier,a,b,c,nquad,ss,swhts,ts,twhts,
     1     ss0,ts0,idxes,curve,par1,par2,par3,par4)
        implicit double precision (a-h,o-z)
        dimension ss(1),swhts(1),ts(1),twhts(1),ss0(1),ts0(1),idxes(1)
        external curve
        ier = 0
c
c       Collect all of the points less than c.
c
        nn = 0
        do 1000 j=1,nquad
        if (ss(j) .lt. c) then
        nn=nn+1
        ss0(nn)=ss(j)
        idxes(nn)=j
        endif
 1000 continue
c
c       Sort the points less than 0 in descending order.
c
        call quicksort(nn,ss0,idxes)
        call reverse(nn,ss0,idxes)
c
c       Perform the mapping
c
        t0=c
        s0=0
c
        call arcmap00(ier,t0,s0,curve,par1,par2,par3,par4,nn,ts0,ss0)
        if (ier .ne. 0) return
c
c       Copy 'em into the output array.
c
        do 1100 j=1,nn
        idx=idxes(j)
        ts(idx)=ts0(j)
        twhts(idx)=swhts(idx)
 1100 continue
c
c       Now gather all the points > c.
c
        nn = 0
        do 2000 j=1,nquad
        if (ss(j) .ge. c) then
        nn=nn+1
        ss0(nn)=ss(j)
        idxes(nn)=j
        endif
 2000 continue
c
c       Sort 'em in ascending order.
c
        call quicksort(nn,ss0,idxes)
c
c       Perform the mapping
c
        t0=c
        s0=0
c
        call arcmap00(ier,t0,s0,curve,par1,par2,par3,par4,nn,ts0,ss0)
        if (ier .ne. 0) return
c
c       Copy 'em into the output array ... make sure to update the 
c       weights.
c
        do 2100 j=1,nn
        idx=idxes(j)
        ts(idx)=ts0(j)
        twhts(idx)=swhts(idx)
 2100 continue
c
c       Reweight the quadrature formula.
c
        do 3000 j=1,nquad
        t    = ts(j)
        twht = twhts(j)
c
        call curve(t,x,y,dx,dy,par1,par2,par3,par4)
        dr = sqrt(dx**2+dy**2)
        twhts(j) = twht/dr
 3000 continue
c
        end



        subroutine arcmap00(ier,t0,s0,curve,par1,par2,par3,par4,n,ts,ss)
        implicit double precision (a-h,o-z)
        dimension ts(1),twhts(1),ss(1),swhts(1)
        data isinit /-1/
        external curve
        save
c
        ier = 0
c
        if (isinit .eq. -1) then
        isinit=1
        call mach_zero(eps)
        eps=eps*10
        endif
c
        m=30
c
        t00 = t0
        s00 = s0
c
        do 1000 j=1,n
c
        s=ss(j)
c
        call arcmap_newt(ier,curve,par1,par2,par3,par4,eps,t00,s00,t,s)
c
        if (ier .ne. 0) return
c
        ts(j)=t
        t00 = t
        s00 = s
c
 1000 continue
c
        end
c
c
c
        subroutine arcmap_newt(ier,curve,par1,par2,par3,par4,
     1    eps,t0,s0,t,s)
        implicit double precision (a-h,o-z)
        external curve
c
        ier      = 0
        maxiters = 12
        maxsteps = 12
c
c       Do liner interpolation to get an initial guess.
c
        call curve(t0,x0,y0,dx0,dy0,par1,par2,par3,par4)
        dsdt=sqrt(dx0**2+dy0**2)
        t = t0+(s-s0)/dsdt
c
c       Perform Newton iterations.
c
        do 1000 iter=1,maxiters
        call curve(t,x,y,dx,dy,par1,par2,par3,par4)
        df=sqrt(dx**2+dy**2)
        call get_arclength(t0,t,f,curve,par1,par2,par3,par4)
        f=(f-(s-s0))
        delta = f/df
c
        if (abs(f) .lt. eps) goto 1200
c
        alpha=1
        do 1100 istep=1,maxsteps
        tnew=t-alpha*delta
        call get_arclength(t0,tnew,fnew,curve,par1,par2,par3,par4)
        fnew=(fnew-(s-s0))
        if (abs(fnew) .lt. abs(f)) goto 1300
        alpha=alpha/2.0d0
 1100 continue
        ier = 128
        return
 1300 continue
        t=tnew
 1000 continue
        ier = 16
 1200 continue
        end
c
c
c
        subroutine get_arclength(a,b,rl,curve,par1,par2,par3,par4)
        implicit double precision (a-h,o-z)
        external arcwrapper,curve
        data isinit /-1/
        save
c
c       Return the arclength of a portion of the curve.
c
        if (isinit .eq. -1) then
        m = 30
        call mach_zero(eps)
        eps=eps*10
        isinit=1
        endif
c
        call arcgauss(ier,a,b,eps,m,rl,arcwrapper,curve,par1,par2,par3,
     1      par4,par5,par6,par7)
c
        end
c
c
c
        subroutine arcwrapper(t,val,curve,par1,par2,par3,par4,par5,par6,
     1    par7)
        implicit double precision (a-h,o-z)
        external curve
        call curve(t,x,y,dx,dy,par1,par2,par3,par4)
c
        val=sqrt(dx**2+dy**2)
        end


        subroutine arcgauss(ier,a,b,eps,m,val,funuser,par1,par2,par3,
     1    par4,par5,par6,par7,par8)
        implicit double precision (a-h,o-z)
        save
        dimension w(1 000 000),xs(1000),whts(1000)
        data m0 /0/
        external funuser
c
c       This procedure uses adaptive Legendre integration to evaluate
c       an integral of the form
c
c           \int_a^b f(x) dx
c
c       where f(x) is a function supplied by the user via an external
c       subroutine.
c
c       NOTE: the failure of this procedure is usually caused by asking
c       for an unreasonable level of accuracy.
c
c                              Input Parameters:
c
c   (a,b) - the domain of integration
c   eps - precision for the computation 
c   m - order of the quadrature formula 
c   funuser - a user-supplied external subroutine with calling sequence
c
c          subroutine funuser(x,val,par1,par2,par3,par4,par5,par6,par7,
c            par8)
c
c          Return the value of the function f at the point x.
c
c   par? - user-supplied arbitrarily typed variables
c
c                              Output Parameters:
c
c   ier - error return code;
c       ier = 0    means the procedure was successful
c       ier = 4    means that the internal work array was insufficient
c       ier = 512  means that maximum recursion depth was exceeded                   
c
c
        ier = 0 
        done = 1
        dtwo = 2
        MAXDEPTH = 1000
c
c       Initialize the quadrature, if necessary.
c
        if (m0 .ne. m) then
           m0 = m
           call legequad(m,xs,whts)
c           itype=1
c           call legeexps(itype,m,xs,u,v,whts)
        endif
c
        val = 0
c
c       Initialize a stack at the end of the work array.
c
        istack0 = 1 000 000
        istack  = istack0
c
        call arcint(a,b,m,val0,funuser,par1,par2,par3,par4,par5,
     1    par6,par7,par8,xs,whts)
c
        ilevel=0
c
        w(istack)=a
        w(istack-1)=b
        w(istack-2)=val0
        w(istack-3)=ilevel
c
        istack=istack-4
c
 1000 continue
c
c       Pop an element off the stack, assuming it isn't empty.
c
        if (istack .eq. istack0) return
c
        istack=istack+4
c
        a0=w(istack)
        b0=w(istack-1)
        val0=w(istack-2)
        ilevel=w(istack-3)
c
        if (ilevel .gt. MAXDEPTH) then
           ier = 512
           return
        endif
c
c       Compute the integral over the left and right half of the
c       interval ...
c
        c0 = (a0+b0)/dtwo
c
        call arcint(a0,c0,m,vall,funuser,par1,par2,par3,par4,par5,
     1    par6,par7,par8,xs,whts)
c
        call arcint(c0,b0,m,valr,funuser,par1,par2,par3,par4,par5,
     1    par6,par7,par8,xs,whts)
c
c       Compare the two approximations we have obtained ... and
c       if we have integrated with sufficient accuracy, add the
c       contribution to the output value.
c
        val1 = valr+vall
        errabs = abs(val0-val1)
c
        if (errabs .le. eps) then
           val=val+val0
           goto 1000
        endif
c
c       Otherwise push the two children onto the stack and continue.
c
        if (istack .lt. 9) then
           ier = 4
           return
        endif
c
        w(istack)=a0
        w(istack-1)=c0
        w(istack-2)=vall
        w(istack-3)=ilevel+1
        istack=istack-4
c
        w(istack)=c0
        w(istack-1)=b0
        w(istack-2)=valr
        w(istack-3)=ilevel+1
        istack=istack-4
c
        goto 1000
c
        end


        subroutine arcint(a,b,m,val,funuser,par1,par2,par3,par4,par5,
     1    par6,par7,par8,xs,whts)
        implicit double precision (a-h,o-z)
        dimension xs(1),whts(1)
        external funuser
c
        dtwo=2
        alpha=(b-a)/dtwo
        beta=(b+a)/dtwo
c
        val=0
c
        do 1000 j=1,m
        x=alpha*xs(j)+beta
        wht=alpha*whts(j)        
        call funuser(x,f,par1,par2,par3,par4,par5,par6,par7,par8)
        val=val+f*wht
 1000 continue
c
c
        end
