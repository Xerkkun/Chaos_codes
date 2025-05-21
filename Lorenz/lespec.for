
        program ode
c
c       n=number of nonlinear odes
c       nn=n*(n+1)=total number of odes
c       nn1=nn+1 (for compilers that start arrays at element 0)
c
        parameter (n=3, nn=12, nn1=13)
c
        implicit double precision(a-h,o-z)
c
        external fcn
c
        dimension y(nn1),yprime(nn1),v(nn1),A(nn1),B(nn1),C(nn1)
        dimension D(nn1),cum(n),znorm(n),gsc(n)
c
c       NSTEP is the total number of reorthonormalization steps that
c           will be performed.  Set this to any very large value.
c       IRATE is the number of numerical integration time steps per
c          reorthonormalization.
c       STPSZE is the integration stepsize in seconds.  Choose a #
c          of "seconds" small compared to the "mean orbital period".
c       IO is the rate at which output is generated.  Rather than
c          outputing the spectrum for each of the NSTEP reorthonorms,
c          values are generated for every IOth reorthonorm. 
c
        write(*,111)
111     format(1x,'nstep, stpsze, irate, io : ')
        read(*,*)  nstep,stpsze,irate,io
c
c       initial conditions for nonlinear ODEs
c       (** Choose within the system's basin of attraction **)
c
        v(1)=1.0
        v(2)=1.0
        v(3)=1.0
        tme=0.0
c
c       initial conditions for linearized ODEs
c       (** Leave these #s alone! They are problem independent! **)
c
        do 10 i=n+1,nn
           v(i)=0.0
10      continue
        do 20 i=1,n
           v((n+1)*i)=1.0
           cum(i)=0.0
20      continue
c
        do 100 m=1,nstep
c
        do 25 j=1,irate
c
c       *****************************************************
c
        do 26 i=1,nn
           y(i)=v(i)
26      continue
        t=tme
        call fcn(t, y, yprime)
        do 27 i=1,nn
           A(i)=yprime(i)
27      continue
c
c       *****************************************************
c
        do 28 i=1,nn
           y(i)=v(i)+(stpsze*A(i))/2.0
28      continue
        t=tme+stpsze/2.0
        call fcn(t, y, yprime)
        do 29 i=1,nn
        B(i)=yprime(i)
29      continue
c
c       *****************************************************
c
        do 30 i=1,nn
           y(i)=v(i)+(stpsze*B(i))/2.0
30      continue
        t=tme+stpsze/2.0
        call fcn(t, y, yprime)
        do 31 i=1,nn
           C(i)=yprime(i)
31      continue
c
c       *****************************************************
c
        do 32 i=1,nn
           y(i)=v(i)+(stpsze*C(i))
32      continue
        t=tme+stpsze
        call fcn(t, y, yprime)
        do 33 i=1,nn
           D(i)=yprime(i)
33      continue
c
c       *****************************************************
c
        do 34 i=1,nn
           v(i)=v(i)+stpsze*(A(i)+D(i)+2.0*(B(i)+C(i)))/6.0
34      continue
        tme=tme+stpsze
c
c       *****************************************************
c
25      continue
c
c       construct new orthonormal basis by gram-schmidt
c
c       normalize first vector
c
        znorm(1)=0.0
        do 38 j=1,n
           znorm(1)=znorm(1)+v(n*j+1)**2
38      continue
        znorm(1)=sqrt(znorm(1))
        do 40 j=1,n
           v(n*j+1)=v(n*j+1)/znorm(1)
40      continue
c
c       generate new orthonormal set
c
c       make j-1 gsr coefficients
c
        do 80 j=2,n
c
        do 50 k=1,j-1
           gsc(k)=0.0
           do 50 l=1,n
              gsc(k)=gsc(k)+v(n*l+j)*v(n*l+k)
50      continue
c
c       construct a new vector
c
        do 60 k=1,n
           do 60 l=1,j-1
              v(n*k+j)=v(n*k+j)-gsc(l)*v(n*k+l)
60      continue
c
c       calculate the vector's norm
c
        znorm(j)=0.0
        do 70 k=1,n
           znorm(j)=znorm(j)+v(n*k+j)**2
70      continue
        znorm(j)=sqrt(znorm(j))
c
c       normalize the new vector
c
        do 80 k=1,n
           v(n*k+j)=v(n*k+j)/znorm(j)
80      continue
c
c       update running vector magnitudes
c
        do 90 k=1,n
           cum(k)=cum(k)+alog(znorm(k))/alog(2.0)
90      continue
c
c       normalize exponent and print every io iterations
c
        if(mod(m,io).ne.0) goto 100
        write(5,334) tme,(cum(k)/tme,k=1,n)
334     format(1x,f12.6,2x,e12.6,2x,e12.6,2x,e12.6)
c
100     continue
c
        end
c
c
***************************************************************
c
        subroutine fcn(t,y,yprime)
c
        implicit double precision(a-h,o-z)
c
        dimension y(13),yprime(13)
c
c       parameters for Lorenz ODEs
c
        b=4.0
        sg=16.0
        r=45.92
c
c       Nonlinear Lorenz equations
c
        yprime(1)=sg*(y(2)-y(1))
        yprime(2)=-y(1)*y(3)+r*y(1)-y(2)
        yprime(3)=y(1)*y(2)-b*y(3)
c
c       Linearized Lorenz equations
c
	do 10 i=0,2
           yprime(4+i)=sg*(y(7+i)-y(4+i))
           yprime(7+i)=(r-y(3))*y(4+i)-y(7+i)-y(1)*y(10+i)
           yprime(10+i)=y(2)*y(4+i)+y(1)*y(7+i)-b*y(10+i)
10      continue
c
        return
        end

