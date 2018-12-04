cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This file contains subroutines which provide the user with
c       parameterizations of several simple planar curves.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        subroutine curves_set_parameters(ifarc,epsrad,epscut,dsub,
     -    istrat,epslen)
        implicit double precision (a-h,o-z)
c
c       Set reasonable defaults in case the user neglects to call the 
c       routine.
c     
        data ifarc0   / 1       /
        data epsrad0  / 0.1d0   /
        data epscut0  / 1.0d-15 /
        data dsub0    / 3.0d0   /
        data istrat0  / 1       /
        data epslen0  / 1.0d0   /
c
        save
c
c       Set the parameters for discretization.
c
        ifarc0  = ifarc
        epsrad0 = epsrad
        epscut0 = epscut
        dsub0   = dsub
        istrat0 = istrat
        epslen0 = epslen
        return


        entry curves_get_parameters(ifarc,epsrad,epscut,dsub,istrat,
     -    epslen)
        ifarc   = ifarc0
        epsrad  = epsrad0
        epscut  = epscut0
        dsub    = dsub0
        epslen  = epslen0
        istrat  = istrat0
        end


        subroutine ellipse(t,x,y,dx,dy,a,b,angle,par4)
        implicit double precision (a-h,o-z)
c
c       This routine supplies a parameterization of an ellipse with
c       axis-lengths a and b.  It is centered at (0,0) and a counter-
c       clockwise rotation of angle radians is applied.
c
        x = cos(angle)*cos(t)*a-sin(angle)*sin(t)*b
        y = sin(angle)*cos(t)*a+cos(angle)*sin(t)*b
c
        dx = cos(angle)*(-sin(t)*a)-sin(angle)*(cos(t)*b)
        dy = sin(angle)*(-sin(t)*a)+cos(angle)*(cos(t)*b)
c
        end



        subroutine clover(t,x,y,dx,dy,angle,par2,par3,par4)
        implicit double precision (a-h,o-z)
c
c       This routine supplies a parameterization of an 4-leaf clover.
c       It is centered at (0,0) and a counter-clockwise rotation of angle 
c       radians is applied.
c
        dleafs  = 4.0d0
        delta   = 2.0d0
c
        xx  = cos(t)*(delta+sin(dleafs*t))
        dxx = dleafs*cos(dleafs*t)*cos(t)-sin(t)*(delta+sin(dleafs*t))
c
        yy  = sin(t)*(delta+sin(dleafs*t))
        dyy = dleafs*cos(dleafs*t)*sin(t)+cos(t)*(delta+sin(dleafs*t))
c
        rot11 = cos(angle)
        rot12 = -sin(angle)
        rot21 = sin(angle)
        rot22 = cos(angle)
c
        x = rot11*xx+rot12*yy
        y = rot21*xx+rot22*yy
c
        dx = rot11*dxx+rot12*dyy
        dy = rot21*dxx+rot22*dyy
c
        end


        subroutine snowcone(t,x,y,dx,dy,theta,angle,par3,par4)
        implicit double precision (a-h,o-z)
        dimension rot(2,2)
c
c       Supply a parameterization of a snowcone curve with a corner of
c       angle theta at (0,0).  It is rotated about the point (-1,0) in a
c       counter-clockwise fashion by angle radians.
c
        beta = tan(theta/2)
c
        if (t .lt. 0) then
c
        xx  = -2*sin(t/2)
        yy  = -beta*sin(t)
c
        dxx =  -cos(t/2)
        dyy = -beta*cos(t)
c
        else
        xx  =  2*sin(t/2)
        yy  =  -beta*sin(t)
c
        dxx  =  cos(t/2)
        dyy  = -beta*cos(t)
c
        endif
c
        rot11 = cos(angle)
        rot12 = -sin(angle)
        rot21 = sin(angle)
        rot22 = cos(angle)
c
        x = rot11*xx + rot12*yy
        y = rot21*xx + rot22*yy
c
        dx = rot11*dxx + rot12*dyy
        dy = rot21*dxx + rot22*dyy
c
        end



        subroutine boomerang(t,x,y,dx,dy,theta,angle,par3,par4)
        implicit double precision (a-h,o-z)
c
        beta = tan(theta/2)
c
        if (t .lt. 0) then
        xx=2*sin(3*t/2)
        yy=beta*sin(t)
c
        dxx = 3*cos(3*t/2)
        dyy = beta*cos(t)
c
        else
        xx=-2*sin(3*t/2)
        yy=beta*sin(t)
c
        dxx = -3*cos(3*t/2)
        dyy = beta*cos(t)
c
        endif
c
        rot11 = cos(angle)
        rot12 = -sin(angle)
        rot21 = sin(angle)
        rot22 = cos(angle)
c
        x = rot11*xx + rot12*yy
        y = rot21*xx + rot22*yy
c
        dx = rot11*dxx + rot12*dyy
        dy = rot21*dxx + rot22*dyy
c
        end



        subroutine clover2(t,x,y,dx,dy,iregion,angle,par3,par4)
        implicit double precision (a-h,o-z)
        data s2  / 1.41421356237309504880168872420969808d0 /
        data pi  / 3.14159265358979323846264338327950288d0 /
        data s3  / 0.35355339059327376220042218105242452d0 /
        data eps / 1.0d-7 /
c
c       This subroutine supplies a parameterization of a four-leaf
c       clover with 4 corner points.  
c
c        if (iregion .eq. 1) then
c
        if (iregion .eq. 1) goto 1000
        if (iregion .eq. 2) goto 2000
        if (iregion .eq. 3) goto 3000
        if (iregion .eq. 4) goto 4000
        stop
c
c       REGION AROUND THE FIRST CORNER 
c
 1000 continue

        if (t .gt. -eps .AND. t .lt. 0) then

        x = (5*t)/(2.*s2) - (7*t**2)/(4.*s2) - 
     -  (29*t**3)/(12.*s2) + (79*t**4)/(48.*s2) + 
     -  (49*t**5)/(48.*s2)

        y = (3*t)/(2.*s2) + (9*t**2)/(4.*s2) - 
     -  (9*t**3)/(4.*s2) - (27*t**4)/(16.*s2) + 
     -  (81*t**5)/(80.*s2)

        dx = 5/(2.*s2) - (7*t)/(2.*s2) - 
     -  (29*t**2)/(4.*s2) + (79*t**3)/(12.*s2) + 
     -  (245*t**4)/(48.*s2) - (727*t**5)/(240.*s2)

        dy = 3/(2.*s2) + (9*t)/(2.*s2) - 
     -  (27*t**2)/(4.*s2) - (27*t**3)/(4.*s2) + 
     -  (81*t**4)/(16.*s2) + (243*t**5)/(80.*s2)
c
        goto 5000
        endif

        
        if ( t .gt. 0 .AND. t .lt. eps) then
        x = (-3*t)/(2.*s2) + (9*t**2)/(4.*s2) + 
     -  (9*t**3)/(4.*s2) - (27*t**4)/(16.*s2) - 
     -  (81*t**5)/(80.*s2)

        y =(-5*t)/(2.*s2) - (7*t**2)/(4.*s2) + 
     -  (29*t**3)/(12.*s2) + (79*t**4)/(48.*s2) - 
     -  (49*t**5)/(48.*s2)

        dx = -3/(2.*s2) + (9*t)/(2.*s2) + 
     -  (27*t**2)/(4.*s2) - (27*t**3)/(4.*s2) - 
     -  (81*t**4)/(16.*s2) + (243*t**5)/(80.*s2)

        dy = -5/(2.*s2) - (7*t)/(2.*s2) + 
     -  (29*t**2)/(4.*s2) + (79*t**3)/(12.*s2) - 
     -  (245*t**4)/(48.*s2) - (727*t**5)/(240.*s2)
        goto 5000
        endif

        if (t .lt. 0) then
        x = (s2+ Cos(Pi/4. + t)*(-2 + 8*Cos(t)*Sin(t)))/4.0d0
        y = (s2-2*Cos(Pi/4. + 3*t))/4.0d0
        dx = (2*Cos(t) + 3*Cos(3*t) + 2*Sin(t) - 3*Sin(3*t))*s3
        dy = (3*(Cos(3*t) + Sin(3*t)))*s3
        goto 5000
        endif

        if (t .gt. 0) then
        x = (s2 - 2*Sin(Pi/4. + 3*t))/4.0d0
        y = (s2 - 2*(1 + 4*Cos(t)*Sin(t))*Sin(Pi/4. + t))/4.
        dx = (-3*(Cos(3*t) - Sin(3*t)))*s3
        dy = -(2*Cos(t) + 3*Cos(3*t) - 2*Sin(t) + 3*Sin(3*t))*s3
        goto 5000
        endif

c
c       REGION AROUND THE SECOND CORNER
c
 2000 continue
c
        if (t .gt. -eps .AND. t .lt. 0) then
        x = (-3*t)/(2.*s2) - (9*t**2)/(4.*s2) + 
     -  (9*t**3)/(4.*s2) + (27*t**4)/(16.*s2) - 
     -  (81*t**5)/(80.*s2)

        y = (5*t)/(2.*s2) - (7*t**2)/(4.*s2) - 
     -  (29*t**3)/(12.*s2) + (79*t**4)/(48.*s2) + 
     -  (49*t**5)/(48.*s2)

        dx = -3/(2.*s2) - (9*t)/(2.*s2) + 
     -  (27*t**2)/(4.*s2) + (27*t**3)/(4.*s2) - 
     -  (81*t**4)/(16.*s2) - (243*t**5)/(80.*s2)

        dy = 5/(2.*s2) - (7*t)/(2.*s2) - 
     -  (29*t**2)/(4.*s2) + (79*t**3)/(12.*s2) + 
     -  (245*t**4)/(48.*s2) - (727*t**5)/(240.*s2)
        goto 5000
        endif

        if (t .gt. 0 .AND. t .lt. eps) then
c
        x = (5*t)/(2.*s2) + (7*t**2)/(4.*s2) - 
     -  (29*t**3)/(12.*s2) - (79*t**4)/(48.*s2) + 
     -  (49*t**5)/(48.*s2)

        y = (-3*t)/(2.*s2) + (9*t**2)/(4.*s2) + 
     -  (9*t**3)/(4.*s2) - (27*t**4)/(16.*s2) - 
     -  (81*t**5)/(80.*s2)

        dx = 5/(2.*s2) + (7*t)/(2.*s2) - 
     -  (29*t**2)/(4.*s2) - (79*t**3)/(12.*s2) + 
     -  (245*t**4)/(48.*s2) + (727*t**5)/(240.*s2)

        dy = -3/(2.*s2) + (9*t)/(2.*s2) + 
     -  (27*t**2)/(4.*s2) - (27*t**3)/(4.*s2) - 
     -  (81*t**4)/(16.*s2) + (243*t**5)/(80.*s2)
        goto 5000
        endif

        if (t .lt. 0) then
        x = (-1 + Cos(3*t) - Sin(3*t))*s3
        y = (s2+ Cos(Pi/4. + t)*(-2 + 4*Sin(2*t)))/4.0d0
        dx = (-3*(Cos(3*t) + Sin(3*t)))*s3
        dy = (2*Cos(t) + 3*Cos(3*t) + 2*Sin(t) - 3*Sin(3*t))*s3
        goto 5000
        endif
c
        if (t .gt. 0) then
        x = (-s2 + (2 + 4*Sin(2*t))*Sin(Pi/4. + t))/4.0d0
        y =-(-1 + Cos(3*t) + Sin(3*t))*s3
        dx = (2*Cos(t) + 3*Cos(3*t) - 2*Sin(t) + 3*Sin(3*t))*s3
        dy = (-3*(Cos(3*t) - Sin(3*t)))*s3
        goto 5000
        endif
        return
c
c       REGION AROUND THE THIRD CORNER
c
 3000 continue

        if (t .gt. -eps .AND. t .lt. 0) then
        x = (-5*t)/(2.*s2) + (7*t**2)/(4.*s2) + 
     -  (29*t**3)/(12.*s2) - (79*t**4)/(48.*s2) - 
     -  (49*t**5)/(48.*s2)

        y = (-3*t)/(2.*s2) - (9*t**2)/(4.*s2) + 
     -  (9*t**3)/(4.*s2) + (27*t**4)/(16.*s2) - 
     -  (81*t**5)/(80.*s2)

        dx = -5/(2.*s2) + (7*t)/(2.*s2) + 
     -  (29*t**2)/(4.*s2) - (79*t**3)/(12.*s2) - 
     -  (245*t**4)/(48.*s2) + (727*t**5)/(240.*s2)

        dy = -3/(2.*s2) - (9*t)/(2.*s2) + 
     -  (27*t**2)/(4.*s2) + (27*t**3)/(4.*s2) - 
     -  (81*t**4)/(16.*s2) - (243*t**5)/(80.*s2)
        goto 5000
        endif
c
        if (t .gt. 0 .AND. t .lt. eps) then

        x = (3*t)/(2.*s2) - (9*t**2)/(4.*s2) - 
     -  (9*t**3)/(4.*s2) + (27*t**4)/(16.*s2) + 
     -  (81*t**5)/(80.*s2)

        y = (5*t)/(2.*s2) + (7*t**2)/(4.*s2) - 
     -  (29*t**3)/(12.*s2) - (79*t**4)/(48.*s2) + 
     -  (49*t**5)/(48.*s2)

        dx = 3/(2.*s2) - (9*t)/(2.*s2) - 
     -  (27*t**2)/(4.*s2) + (27*t**3)/(4.*s2) + 
     -  (81*t**4)/(16.*s2) - (243*t**5)/(80.*s2)

        dy =5/(2.*s2) + (7*t)/(2.*s2) - 
     -  (29*t**2)/(4.*s2) - (79*t**3)/(12.*s2) + 
     -  (245*t**4)/(48.*s2) + (727*t**5)/(240.*s2)
        goto 5000
        endif

        if (t .lt. 0) then
        x = (-s2 + Cos(Pi/4. + t)*(2 - 8*Cos(t)*Sin(t)))/4.0d0
        y = (-s2 + 2*Cos(Pi/4. + 3*t))/4.0d0
        dx = -(2*Cos(t) + 3*Cos(3*t) + 2*Sin(t) - 3*Sin(3*t))*s3
        dy = (-3*(Cos(3*t) + Sin(3*t)))*s3
        goto 5000
        endif

        if (t .gt. 0) then
        x = (-s2 + 2*Sin(Pi/4. + 3*t))/4.0d0
        y = (-s2 + (2 + 8*Cos(t)*Sin(t))*Sin(Pi/4. + t))/4.0d0
        dx = (3*(Cos(3*t) - Sin(3*t)))*s3
        dy = (2*Cos(t) + 3*Cos(3*t) - 2*Sin(t) + 3*Sin(3*t))*s3
        goto 5000
        endif

        return
c
c       REGION AROUND THE FOURTH CORNER
c
 4000 continue

        if (t .gt. -eps .AND. t .lt. 0) then
        x = (3*t)/(2.*s2) + (9*t**2)/(4.*s2) - 
     -  (9*t**3)/(4.*s2) - (27*t**4)/(16.*s2) + 
     -  (81*t**5)/(80.*s2)

        y = (-5*t)/(2.*s2) + (7*t**2)/(4.*s2) + 
     -  (29*t**3)/(12.*s2) - (79*t**4)/(48.*s2) - 
     -  (49*t**5)/(48.*s2)

        dx = 3/(2.*s2) + (9*t)/(2.*s2) - 
     -  (27*t**2)/(4.*s2) - (27*t**3)/(4.*s2) + 
     -  (81*t**4)/(16.*s2) + (243*t**5)/(80.*s2)

        dy = -5/(2.*s2) + (7*t)/(2.*s2) + 
     -  (29*t**2)/(4.*s2) - (79*t**3)/(12.*s2) - 
     -  (245*t**4)/(48.*s2) + (727*t**5)/(240.*s2)
        goto 5000
        return
        endif
c
        if (t .gt. 0 .AND. t .lt. eps) then
        x = (-5*t)/(2.*s2) - (7*t**2)/(4.*s2) + 
     -  (29*t**3)/(12.*s2) + (79*t**4)/(48.*s2) - 
     -  (49*t**5)/(48.*s2)

        y = (3*t)/(2.*s2) - (9*t**2)/(4.*s2) - 
     -  (9*t**3)/(4.*s2) + (27*t**4)/(16.*s2) + 
     -  (81*t**5)/(80.*s2)

        dx = -5/(2.*s2) - (7*t)/(2.*s2) + 
     -  (29*t**2)/(4.*s2) + (79*t**3)/(12.*s2) - 
     -  (245*t**4)/(48.*s2) - (727*t**5)/(240.*s2)

        dy = 3/(2.*s2) - (9*t)/(2.*s2) - 
     -  (27*t**2)/(4.*s2) + (27*t**3)/(4.*s2) + 
     -  (81*t**4)/(16.*s2) - (243*t**5)/(80.*s2)
        goto 5000
        return
        endif

        if (t .lt. 0) then
        x = (1 - Cos(3*t) + Sin(3*t))*s3
        y = (-s2 + Cos(Pi/4. + t)*(2 - 8*Cos(t)*Sin(t)))/4.0d0

        dx = (3*(Cos(3*t) + Sin(3*t)))*s3
        dy = -(2*Cos(t) + 3*Cos(3*t) + 2*Sin(t) - 3*Sin(3*t))*s3
        goto 5000
        endif
c
        if (t .gt. 0) then
        x = (1 - (Cos(t) + Sin(t))*(1 + 2*Sin(2*t)))*s3
        y = (-s2 + 2*Sin(Pi/4. + 3*t))/4.0d0
        dx = -(2*Cos(t) + 3*Cos(3*t) - 2*Sin(t) + 3*Sin(3*t))*s3
        dy = (3*(Cos(3*t) - Sin(3*t)))*s3
        goto 5000
        endif

 5000 continue
        xx  = x
        yy  = y
        dxx = dx
        dyy = dy
c
        rot11 = cos(angle)
        rot12 = -sin(angle)
        rot21 = sin(angle)
        rot22 = cos(angle)
c
        x = rot11*xx+rot12*yy
        y = rot21*xx+rot22*yy
c
        dx = rot11*dxx+rot12*dyy
        dy = rot21*dxx+rot22*dyy
        
        end


        subroutine polygon(t,x,y,dx,dy,nverts,verts,ivert,angle)
        implicit double precision (a-h,o-z)
        dimension verts(2,1)
c        
        rot11 = cos(angle)
        rot12 =-sin(angle)
        rot21 = sin(angle)
        rot22 = cos(angle)
c
        x1 = verts(1,ivert)
        y1 = verts(2,ivert)
c
        if (t .gt. 0) then
        ivert2 = ivert+1
        if (ivert2 .gt. nverts) ivert2 = 1
c
        x2 = verts(1,ivert2)
        y2 = verts(2,ivert2)
c
        dxx = x2-x1
        dyy = y2-y1
c
        xx = dxx*t
        yy = dyy*t
c
        else
        ivert2 = ivert-1
        if (ivert2 .eq. 0 ) ivert2 = nverts
c
        x2 = verts(1,ivert2)
        y2 = verts(2,ivert2)
c
        dxx = x1-x2
        dyy = y1-y2
c
        xx = dxx*t
        yy = dyy*t
c
        endif
c
c
c
        x = rot11*xx+rot12*yy
        y = rot21*xx+rot22*yy
c
        dx = rot11*dxx+rot12*dyy
        dy = rot21*dxx+rot22*dyy
c
        end



        subroutine add_ellipse(disc,iregion,x0,y0,d1,d2,angle)
        implicit double precision (a-h,o-z)        
        dimension disc(1)
        external ellipse
c
        call disc_start_component(ier,disc)
c
        iregion = iregion + 1
        rx      = x0
        ry      = y0
c
        call curves_get_parameters(ifarc,epsrad,epscut,dsub,istrat,
     -    epslen)
        eps = epslen
c
c$$$        istrat  = 0
c$$$        eps     = .5d0
c
        pi = acos(-1.0d0)
        a  = -pi
        b  = pi
c
        call disc_add_smooth(ier,disc,iregion,rx,ry,istrat,eps,
     1    a,b,ellipse,d1,d2,angle,par4)
        if (ier .ne. 0) then
        call prinf("after disc_add, ier = *",ier,1)
        stop
        endif
c
        call disc_end_component(disc)
        end



        subroutine add_circle(disc,iregion,x0,y0,r)
        implicit double precision (a-h,o-z)        
        dimension disc(1)
        external ellipse
c
        call disc_start_component(ier,disc)
c
        iregion = iregion + 1
        rx      = x0
        ry      = y0
c
c        call curves_get_parameters(ifarc,epsrad,epscut,dsub,istrat,
c     -    epslen)
c
        istrat  = 0
        eps     = 1.5d0
c
        pi = acos(-1.0d0)
        a  = -pi
        b  = pi
c
        call disc_add_smooth(ier,disc,iregion,rx,ry,istrat,eps,
     1    a,b,ellipse,r,r,0.0d0,par4)
        if (ier .ne. 0) then
        call prinf("after disc_add, ier = *",ier,1)
        stop
        endif
c
        call disc_end_component(disc)
        end



        subroutine add_clover(disc,iregion,x0,y0,angle)
        implicit double precision (a-h,o-z)        
        dimension disc(1)
        external clover
c
        call disc_start_component(ier,disc)
c
        iregion = iregion + 1
        rx      = x0
        ry      = y0
c
        call curves_get_parameters(ifarc,epsrad,epscut,dsub,istrat,
     -    epslen)
c
        eps = epslen
c
        pi = acos(-1.0d0)
        a  = -pi
        b  =  pi
c
        call disc_add_smooth(ier,disc,iregion,rx,ry,istrat,eps,
     1    a,b,clover,angle,par2,par3,par4)
        if (ier .ne. 0) then
        call prinf("after disc_add, ier = *",ier,1)
        stop
        endif
c
        call disc_end_component(disc)
c
        end


        subroutine add_snowcone(disc,iregion,x0,y0,theta,angle)
        implicit double precision (a-h,o-z)
        external snowcone
c
        call disc_start_component(ier,disc)
        pi = acos(-1.0d0)
c
        iregion = iregion+1
c
        call curves_get_parameters(ifarc,epsrad,epscut,dsub,istrat,
     -    epslen)
c
c        
c        ifarc  = 0
c        epsrad = 0.2d0
c        epscut = 1.0d-15
c        dsub   = 3
c
        rot11 = cos(angle)
        rot12 =-sin(angle)
        rot21 = sin(angle)
        rot22 = cos(angle)
c
        rxx = -1.0d0
        ryy = 0.0d0
c
        rx = rot11*(rxx)+rot12*(ryy)+x0
        ry = rot21*(rxx)+rot22*(ryy)+y0
c
        call disc_add_corner(ier,disc,iregion,rx,ry,ifarc,
     1    a,b,epsrad,epscut,dsub,snowcone,theta,angle,par3,par4)
        if (ier .ne. 0) then
        call prinf("after disc_add, ier = *",ier,1)
        stop
        endif
c
c        istrat = 1
c        epslen = 0.5d0
c
        aa = -pi
        bb = a
c
        call disc_add_smooth(ier,disc,iregion,rx,ry,istrat,epslen,aa,bb,
     1    snowcone,theta,angle,par3,par4)
        if (ier .ne. 0) then
        call prinf("after disc_add, ier = *",ier,1)
        stop
        endif
c
        aa = b
        bb = pi
c
        call disc_add_smooth(ier,disc,iregion,rx,ry,istrat,epslen,aa,bb,
     1    snowcone,theta,angle,par3,par4)
        if (ier .ne. 0) then
        call prinf("after disc_add, ier = *",ier,1)
        stop

        endif
c
        call disc_end_component(disc)
c
        end



        subroutine add_boomerang(disc,iregion,x0,y0,theta,angle)
        implicit double precision (a-h,o-z)
        external boomerang
c
        call disc_start_component(ier,disc)
        iregion = iregion+1
c
        call curves_get_parameters(ifarc,epsrad,epscut,dsub,istrat,
     -    epslen)

c$$$        ifarc  = 0
c$$$        epsrad = 0.1d0
c$$$        epscut = 1.0d-20
c$$$        dsub   = 3
c
        rot11 = cos(angle)
        rot12 =-sin(angle)
        rot21 = sin(angle)
        rot22 = cos(angle)
c
        rxx = -1.0d0
        ryy = 0.0d0
c
        rx = rot11*rxx+rot12*ryy+x0
        ry = rot21*rxx+rot22*ryy+y0
c
        call disc_add_corner(ier,disc,iregion,rx,ry,ifarc,
     1    a,b,epsrad,epscut,dsub,boomerang,theta,angle,par3,par4)
        if (ier .ne. 0) then
        call prinf("after disc_add, ier = *",ier,1)
        stop
        endif
c
        pi = acos(-1.0d0)
c
c$$$        istrat = 3
c$$$        epslen    = 1.0d-12
c
        aa = -pi
        bb = a
c
        call disc_add_smooth(ier,disc,iregion,rx,ry,istrat,epslen,aa,bb,
     1    boomerang,theta,angle,par3,par4)
        if (ier .ne. 0) then
        call prinf("after disc_add, ier = *",ier,1)
        stop
        endif
c
        aa = b
        bb = pi
c
        call disc_add_smooth(ier,disc,iregion,rx,ry,istrat,epslen,aa,bb,
     1    boomerang,theta,angle,par3,par4)
        if (ier .ne. 0) then
        call prinf("after disc_add, ier = *",ier,1)
        stop
        endif
c
        call disc_end_component(disc)
c
        end


        subroutine add_clover2(disc,iregion,x0,y0,angle)
        implicit double precision (a-h,o-z)
        external clover2
c
c       Add a 4-leaf clover which has 4 corner points.
c
        call disc_start_component(ier,disc)
        
        pi = acos(-1.0d0)
c
        rot11 = cos(angle)
        rot12 =-sin(angle)
        rot21 = sin(angle)
        rot22 = cos(angle)
c
        call curves_get_parameters(ifarc,epsrad,epscut,dsub,istrat,
     -    epslen)

c        ifarc  = 1
c       epsrad = 0.2d0
c       epscut = 1.0d-24
c       dsub   = 3
c
c        istrat = 1
c        epslen = 0.125d0
c
        iregion = iregion+1
        rxx = -1.0d0/(2.0d0*Sqrt(2.0d0))
        ryy = -1.0d0/(2.0d0*Sqrt(2.0d0))
c
        rx = rot11*rxx+rot12*ryy+x0
        ry = rot21*rxx+rot22*ryy+y0
c
        ii = 1
c
        call disc_add_corner(ier,disc,iregion,rx,ry,ifarc,
     1    a,b,epsrad,epscut,dsub,clover2,ii,angle,par3,par4)
        if (ier .ne. 0) then
        call prinf("after disc_add, ier = *",ier,1)
        stop
        endif
c
        aa =  -pi/4
        bb =  a
c
        call disc_add_smooth(ier,disc,iregion,rx,ry,istrat,epslen,aa,bb,
     1    clover2,ii,angle,par3,par4)        
        if (ier .ne. 0) then
        call prinf("after disc_add, ier = *",ier,1)
        stop
        endif
c
        aa =  b
        bb =  pi/4
c
        call disc_add_smooth(ier,disc,iregion,rx,ry,istrat,epslen,aa,bb,
     1    clover2,ii,angle,par3,par4)        
        if (ier .ne. 0) then
        call prinf("after disc_add, ier = *",ier,1)
        stop
        endif
c
        iregion = iregion+1
        rxx =  1.0d0/(2.0d0*Sqrt(2.0d0))
        ryy =- 1.0d0/(2.0d0*Sqrt(2.0d0))
c
        rx = rot11*rxx+rot12*ryy+x0
        ry = rot21*rxx+rot22*ryy+y0
c
        ii = 2
c
        call disc_add_corner(ier,disc,iregion,rx,ry,ifarc,
     1    a,b,epsrad,epscut,dsub,clover2,ii,angle,par3,par4)
        if (ier .ne. 0) then
        call prinf("after disc_add, ier = *",ier,1)
        stop
        endif
c
        aa =  -pi/4
        bb =  a
c
        call disc_add_smooth(ier,disc,iregion,rx,ry,istrat,epslen,aa,bb,
     1    clover2,ii,angle,par3,par4)        
        if (ier .ne. 0) then
        call prinf("after disc_add, ier = *",ier,1)
        stop
        endif
c
        aa =  b
        bb =  pi/4
c
        call disc_add_smooth(ier,disc,iregion,rx,ry,istrat,epslen,aa,bb,
     1    clover2,ii,angle,par3,par4)        
        if (ier .ne. 0) then
        call prinf("after disc_add, ier = *",ier,1)
        stop
        endif
c
c
c
        iregion = iregion+1
        rxx = 1.0d0/(2.0d0*Sqrt(2.0d0))
        ryy = 1.0d0/(2.0d0*Sqrt(2.0d0))
c
        rx = rot11*rxx+rot12*ryy+x0
        ry = rot21*rxx+rot22*ryy+y0
c
        ii = 3
c
        call disc_add_corner(ier,disc,iregion,rx,ry,ifarc,
     1    a,b,epsrad,epscut,dsub,clover2,ii,angle,par3,par4)
        if (ier .ne. 0) then
        call prinf("after disc_add, ier = *",ier,1)
        stop
        endif
c
        aa =  -pi/4
        bb =  a
c
        call disc_add_smooth(ier,disc,iregion,rx,ry,istrat,epslen,aa,bb,
     1    clover2,ii,angle,par3,par4)        
        if (ier .ne. 0) then
        call prinf("after disc_add, ier = *",ier,1)
        stop
        endif
c
        aa =  b
        bb =  pi/4
c
        call disc_add_smooth(ier,disc,iregion,rx,ry,istrat,epslen,aa,bb,
     1    clover2,ii,angle,par3,par4)        
        if (ier .ne. 0) then
        call prinf("after disc_add, ier = *",ier,1)
        stop
        endif
c
c
c
        iregion = iregion+1
        rxx = -  1.0d0/(2.0d0*Sqrt(2.0d0))
        ryy =    1.0d0/(2.0d0*Sqrt(2.0d0))
c
        rx = rot11*rxx+rot12*ryy+x0
        ry = rot21*rxx+rot22*ryy+y0
c
        ii = 4
c
        call disc_add_corner(ier,disc,iregion,rx,ry,ifarc,
     1    a,b,epsrad,epscut,dsub,clover2,ii,angle,par3,par4)
        if (ier .ne. 0) then
        call prinf("after disc_add, ier = *",ier,1)
        stop
        endif
c
        aa =  -pi/4
        bb =  a
c
        call disc_add_smooth(ier,disc,iregion,rx,ry,istrat,epslen,aa,bb,
     1    clover2,ii,angle,par3,par4)        
        if (ier .ne. 0) then
        call prinf("after disc_add, ier = *",ier,1)
        stop
        endif
c
        aa =  b
        bb =  pi/4
c
        call disc_add_smooth(ier,disc,iregion,rx,ry,istrat,epslen,aa,bb,
     1    clover2,ii,angle,par3,par4)        
        if (ier .ne. 0) then
        call prinf("after disc_add, ier = *",ier,1)
        stop
        endif
c
        call disc_end_component(disc)
c        
        end


        subroutine add_polygon(disc,iregion,nverts,verts,x0,y0,angle)
        implicit double precision (a-h,o-z)
        dimension verts(2,nverts)
        external polygon
c
        call disc_start_component(ier,disc)
c
        rot11 = cos(angle)
        rot12 = -sin(angle)
        rot21 = sin(angle)
        rot22 = cos(angle)
c
c        ifarc  = 1
c        epsrad = 0.1d0
c        epscut = 1.0d-20
c        dsub   = 3
c
c        istrat = 1
c        eps    = 2.0d0
c
        call curves_get_parameters(ifarc,epsrad,epscut,dsub,istrat,
     -    epslen)
        eps = epslen
c
        do 1000 ivert=1,nverts
        iregion = iregion+1
        rxx = verts(1,ivert)
        ryy = verts(2,ivert)
c
        rx = rot11*rxx+rot12*ryy+x0
        ry = rot21*rxx+rot22*ryy+y0
c
        call disc_add_corner(ier,disc,iregion,rx,ry,ifarc,
     1    a,b,epsrad,epscut,dsub,polygon,nverts,verts,ivert,angle)
        if (ier .ne. 0) then
        call prinf("after disc_add, ier = *",ier,1)
        stop
        endif
c
        aa = -0.5d0
        bb =  a
c
        call disc_add_smooth(ier,disc,iregion,rx,ry,istrat,eps,aa,bb,
     -    polygon,nverts,verts,ivert,angle)
        if (ier .ne. 0) then
        call prinf("after disc_add, ier = *",ier,1)
        stop
        endif
c
        aa =  b
        bb =  0.5d0
c
        call disc_add_smooth(ier,disc,iregion,rx,ry,istrat,eps,aa,bb,
     -    polygon,nverts,verts,ivert,angle)
        if (ier .ne. 0) then
        call prinf("after disc_add, ier = *",ier,1)
        stop
        endif
c
 1000 continue
c
        call disc_end_component(disc)
c
        end


     
        subroutine add_triangle(disc,iregion,x0,y0,angle)
        implicit double precision (a-h,o-z)
        parameter (nverts = 3)
        dimension verts(2,nverts)
        data verts / -1.25d0,-0.75d0, 
     -                0.75d0,-0.75d0, 
     -                0.75d0, 1.25d0 /
c
        call add_polygon(disc,iregion,nverts,verts,x0,y0,angle)
c
        end


        subroutine add_lshaped(disc,iregion,x0,y0,angle)
        implicit double precision (a-h,o-z)
        parameter ( nverts = 6 )
        dimension verts(2,nverts),verts0(2,nverts)
        data verts  /
     -   0.000000d0, 0.000000d0,
     -   2.000000d0, 0.000000d0,
     -   2.000000d0, 1.000000d0,
     -   1.000000d0, 1.000000d0,
     -   1.000000d0, 2.000000d0,
     -   0.000000d0, 2.000000d0 /
c
        do 1000 i=1,nverts
        verts0(1,i)=verts(1,i)
        verts0(2,i)=verts(2,i)
 1000 continue
c
        call add_polygon(disc,iregion,nverts,verts0,x0,y0,angle)

        end




        subroutine add_starburst(disc,iregion,x0,y0,angle)
        implicit double precision (a-h,o-z)
        parameter ( nverts = 38 )
        dimension verts(2,nverts),verts0(2,nverts)
        data verts  /
     -   0.23086D+01 , 0.16074D+00 , 0.99293D+00 , 0.25591D+00,
     -   0.18847D+01 , 0.84256D+00 , 0.11494D+01 , 0.72927D+00,
     -   0.20278D+01 , 0.19093D+01 , 0.12278D+01 , 0.15695D+01,
     -   0.14012D+01 , 0.25399D+01 , 0.47457D+00 , 0.14481D+01,
     -   0.38038D+00 , 0.21376D+01 , 0.00000D+00 , 0.10000D+01,
     -   -.38038D+00 , 0.21376D+01 , -.47457D+00 , 0.14481D+01,
     -   -.14012D+01 , 0.25399D+01 , -.12278D+01 , 0.15695D+01,
     -   -.20278D+01 , 0.19093D+01 , -.11494D+01 , 0.72927D+00,
     -   -.18847D+01 , 0.84256D+00 , -.99293D+00 , 0.25591D+00,
     -   -.23086D+01 , 0.16074D+00 , -.16818D+01 , -.11710D+00,
     -   -.28805D+01 , -.74240D+00 , -.17670D+01 , -.78996D+00,
     -   -.22281D+01 , -.14138D+01 , -.88446D+00 , -.83276D+00,
     -   -.12368D+01 , -.15810D+01 , -.53098D+00 , -.96248D+00,
     -   -.77108D+00 , -.23530D+01 , -.32039D+00 , -.18005D+01,
     -   0.00000D+00 , -.30000D+01 , 0.32039D+00 , -.18005D+01,
     -   0.77108D+00 , -.23530D+01 , 0.53098D+00 , -.96248D+00,
     -   0.12368D+01 , -.15810D+01 , 0.88446D+00 , -.83276D+00,
     -   0.22281D+01 , -.14138D+01 , 0.17670D+01 , -.78996D+00,
     -   0.28805D+01 , -.74240D+00 , 0.16818D+01 , -.11710D+00 /
c
c       Center the starburst on (0,0).
c
        do 1000 i=1,nverts
        verts0(1,i)=verts(1,i)
        verts0(2,i)=verts(2,i)+.23005d0
 1000 continue
c
        call add_polygon(disc,iregion,nverts,verts0,x0,y0,angle)
        end



        subroutine add_starburst0(disc,iregion,x0,y0,angle)
        implicit double precision (a-h,o-z)
        parameter ( nverts = 22 )
        dimension verts(2,nverts)
        data verts  /
     -   0.21944D+01 , 0.23389D+00 , 0.97840D+00 , 0.42446D+00,
     -   0.20121D+01 , 0.18703D+01 , 0.10433D+01 , 0.16765D+01,
     -   0.80020D+00 , 0.23486D+01 , 0.61232D-16 , 0.10000D+01,
     -   -.80020D+00 , 0.23486D+01 , -.10433D+01 , 0.16765D+01,
     -   -.20121D+01 , 0.18703D+01 , -.97840D+00 , 0.42446D+00,
     -   -.21944D+01 , 0.23389D+00 , -.17831D+01 , -.19005D+00,
     -   -.26911D+01 , -.11675D+01 , -.91767D+00 , -.85299D+00,
     -   -.10701D+01 , -.17196D+01 , -.48984D+00 , -.14377D+01,
     -   -.55109D-15 , -.30000D+01 , 0.48984D+00 , -.14377D+01,
     -   0.10701D+01 , -.17196D+01 , 0.91767D+00 , -.85299D+00,
     -   0.26911D+01 , -.11675D+01 , 0.17831D+01 , -.19005D+00/
c
        call add_polygon(disc,iregion,nverts,verts,x0,y0,angle)
        end




        subroutine add_tank(disc,iregion,x0,y0)
        implicit double precision (a-h,o-z)
        dimension verts(2,12),dknots(2,12)
        data nverts / 12 /
        data verts  / -4.0d0,0.0d0, 4.0d0,0.0d0, 5.0d0,1.5d0,
     1                 2.0d0,1.5d0, 2.1d0,2.0d0, 6.0d0,3.3d0,
     2                 5.8d0,3.5d0, 2.1d0,2.4d0, 2.0d0,3.0d0,
     3                -1.5d0,3.0d0,-1.5d0,1.5d0,-5.0d0,1.5d0 /
        data dknots /  0.0d0,0.0d0, 5.2d0,1.0d0, 0.0d0,0.0d0,
     1                 0.0d0,0.0d0, 0.0d0,0.0d0, 0.0d0,0.0d0,
     2                 0.0d0,0.0d0, 0.0d0,0.0d0, 0.5d0,3.2d0,
     3                -1.67d0,2.2d0, 0.0d0,0.0d0,-5.2d0,1.0d0 /
        external sscurve
        save
c
c
c       Build trivial knots.
c        
        do 0100 j=1,nverts-1
c
        if (abs(dknots(1,j))+abs(dknots(2,j)) .eq. 0.0) then
        dknots(1,j) = (verts(1,j)+verts(1,j+1))/2.0d0
        dknots(2,j) = (verts(2,j)+verts(2,j+1))/2.0d0
        endif
 0100 continue
c
        if (abs(dknots(1,nverts-1))+abs(dknots(1,1)) .eq. 0.0) then
        dknots(1,12) = (verts(1,1)+verts(1,12))/2.0d0
        dknots(2,12) = (verts(2,1)+verts(2,12))/2.0d0
        endif
c
        call sscurve_init(nverts,verts,dknots)
        call disc_start_component(ier,disc)
c
c        rot11 = cos(angle)
c        rot12 = -sin(angle)
c        rot21 = sin(angle)
c        rot22 = cos(angle)
c
c$$$        ifarc  = 1
c$$$        epsrad = 0.1d0
c$$$        epscut = 1.0d-19
c$$$        dsub   = 3
c$$$c
c$$$        istrat = 1
c$$$        eps    = 0.125d0
c
        call curves_get_parameters(ifarc,epsrad,epscut,dsub,istrat,
     -    epslen)
        eps=epslen
c
        do 1000 ivert=1,nverts
c
        call sscurve_delta(ivert,xx,yy)
        iregion = iregion+1
        rx      = x0+xx
        ry      = y0+yy
c
c
        call disc_add_corner(ier,disc,iregion,rx,ry,ifarc,
     1    a,b,epsrad,epscut,dsub,sscurve,ivert,angle,p3,p4)
        if (ier .ne. 0) then
        call prinf("after disc_add, ier = *",ier,1)
        stop
        endif
c
        aa =-0.5d0
        bb = a
c
        call disc_add_smooth(ier,disc,iregion,rx,ry,istrat,eps,aa,bb,
     -    sscurve,ivert,angle,p3,p4)
        if (ier .ne. 0) then
        call prinf("after disc_add, ier = *",ier,1)
        stop
        endif
c
c     
        aa =  b
        bb =  0.5d0
c
        call disc_add_smooth(ier,disc,iregion,rx,ry,istrat,eps,aa,bb,
     -    sscurve,ivert,angle,p3,p4)
        if (ier .ne. 0) then
        call prinf("after disc_add, ier = *",ier,1)
        stop
        endif
c
 1000 continue
c
        call disc_end_component(disc)
        return
        end


        subroutine add_sharkfin(disc,iregion,x0,y0,angle)
        implicit double precision (a-h,o-z)
        dimension verts(2,3),dknots(2,3)
        data nverts / 3 /
c
        data verts  /-2.5d0,-1.5d0,  1.5d0,-1.5d0,   1.5d0, 2.5d0 /
        data dknots /-0.5d0,-1.1d0,  1.1d0, 0.5d0,  -0.5d0, 1.5d0 /

        external sscurve
        save
c
        call sscurve_init(nverts,verts,dknots)
        call disc_start_component(ier,disc)
c
c        rot11 = cos(angle)
c        rot12 = -sin(angle)
c        rot21 = sin(angle)
c        rot22 = cos(angle)
c
        call curves_get_parameters(ifarc,epsrad,epscut,dsub,istrat,
     -    epslen)
        eps=epslen
c
        do 1000 ivert=1,nverts
c
        call sscurve_delta(ivert,xx,yy)
c
        iregion = iregion+1
        rx      = x0+xx
        ry      = y0+yy
c
        call disc_add_corner(ier,disc,iregion,rx,ry,ifarc,
     1    a,b,epsrad,epscut,dsub,sscurve,ivert,angle,p3,p4)
        if (ier .ne. 0) then
        call prinf("after disc_add, ier = *",ier,1)
        stop
        endif
c
        aa =-0.5d0
        bb = a
c
        call disc_add_smooth(ier,disc,iregion,rx,ry,istrat,eps,aa,bb,
     -    sscurve,ivert,angle,p3,p4)
        if (ier .ne. 0) then
        call prinf("after disc_add, ier = *",ier,1)
        stop
        endif
c     
        aa =  b
        bb =  0.5d0
c
        call disc_add_smooth(ier,disc,iregion,rx,ry,istrat,eps,aa,bb,
     -    sscurve,ivert,angle,p3,p4)
        if (ier .ne. 0) then
        call prinf("after disc_add, ier = *",ier,1)
        stop
        endif
c
 1000 continue
c
        call disc_end_component(disc)
        return
        end



        subroutine sscurve(t,x,y,dx,dy,ivert,angle,par3,par4)
        implicit double precision (a-h,o-z)
        dimension verts(2,1000),dknots(2,1000),coefs(6,1000)
        dimension verts0(2,1),dknots0(2,1),coefs2(6,1000)
        save
c
        if (t .gt. 0) then
        i1 = ivert
c
        a0 = coefs(1,i1)
        a1 = coefs(2,i1)
        a2 = coefs(3,i1)
        b0 = coefs(4,i1)
        b1 = coefs(5,i1)
        b2 = coefs(6,i1)
c
        x = a1*t+a2*t**2
        y = b1*t+b2*t**2
c
        dx = a1+2*a2*t
        dy = b1+2*b2*t
c
        return
        endif
c
        i1 = ivert
c
        a0 = coefs2(1,i1)
        a1 = coefs2(2,i1)
        a2 = coefs2(3,i1)
        b0 = coefs2(4,i1)
        b1 = coefs2(5,i1)
        b2 = coefs2(6,i1)
c
        x = a1*t+a2*t**2
        y = b1*t+b2*t**2
c
        dx = a1+2*a2*t
        dy = b1+2*b2*t
c
        return


        entry sscurve_init(nverts0,verts0,dknots0)
        nverts = nverts0
c
        do 1000 j=1,nverts
        verts(1,j)=verts0(1,j)
        verts(2,j)=verts0(2,j)
 1000 continue
c
        do 1100 j=1,nverts
        dknots(1,j)=dknots0(1,j)
        dknots(2,j)=dknots0(2,j)
 1100 continue
c
c       Build the coefficients.
c
        do 1200 j=1,nverts0
c
        t0 = 0.00d0
        t1 = 0.50d0
        t2 = 1.00d0
c
        j1 = j
        j2 = j1+1
        j0 = j1-1
c
        if (j2 .gt. nverts) j2 = 1
        if (j0 .lt. 1)      j0 = nverts
c
        x0 = verts(1,j1)
        x1 = dknots(1,j1)
        x2 = verts(1,j2)
c
        y0 = verts(2,j1)
        y1 = dknots(2,j1)
        y2 = verts(2,j2)
c
        call sspoly_coefs(t0,t1,t2,x0,x1,x2,coefs(1,j))
        call sspoly_coefs(t0,t1,t2,y0,y1,y2,coefs(4,j))
c
        t0 = 0.00d0
        t1 = -0.50d0
        t2 = -1.00d0

        x0 = verts(1,j1)
        x1 = dknots(1,j0)
        x2 = verts(1,j0)
c
        y0 = verts(2,j1)
        y1 = dknots(2,j0)
        y2 = verts(2,j0)
c
        call sspoly_coefs(t0,t1,t2,x0,x1,x2,coefs2(1,j))
        call sspoly_coefs(t0,t1,t2,y0,y1,y2,coefs2(4,j))
c
 1200 continue
c
        return


        entry sscurve_delta(ivert,xx,yy)
        xx = coefs(1,ivert)
        yy = coefs(4,ivert)
        end


c
        subroutine sspoly_coefs(x0,x1,x2,f0,f1,f2,b)
        implicit double precision (a-h,o-z)
        dimension amatr(3,3),b(3)
c
c       Compute the coefficients of the 2nd order polynomial which
c       goes through the three specified points.
c
        b(1) = f0
        b(2) = f1
        b(3) = f2
c
c       Form the vandermonde matrix.
c
        amatr(1,1)=1
        amatr(1,2)=x0
        amatr(1,3)=x0**2
c
        amatr(2,1)=1
        amatr(2,2)=x1
        amatr(2,3)=x1**2
c
        amatr(3,1)=1.0d0
        amatr(3,2)=x2
        amatr(3,3)=x2**2
c
c
        call qrsolv(amatr,3,b,rcond)
c
        return
        end
