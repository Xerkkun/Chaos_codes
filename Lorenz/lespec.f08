program ode
    ! Definición de dp para double precision
    integer, parameter :: dp = kind(0.0d0)
    implicit none

    integer, parameter :: n=3, nn=12, nn1=13
    real(dp) :: y(nn1), yprime(nn1), v(nn1), A(nn1), B(nn1), C(nn1)
    real(dp) :: D(nn1), cum(n), znorm(n), gsc(n)
    real(dp) :: tme, stpsze, t
    integer :: nstep, irate, io, i, j, k, l, m

    interface
        subroutine fcn(t, y, yprime)
            real(dp), intent(in) :: t
            real(dp), intent(in) :: y(:)
            real(dp), intent(out) :: yprime(:)
        end subroutine fcn
    end interface

    ! Lectura de parámetros
    write(*,'(A)', advance='no') 'nstep, stpsze, irate, io : '
    read(*,*) nstep, stpsze, irate, io

    ! Condiciones iniciales
    v(1:3) = 1.0_dp
    tme = 0.0_dp
    v(n+1:nn) = 0.0_dp
    do i = 1, n
        v((n+1)*i) = 1.0_dp
        cum(i) = 0.0_dp
    end do

    ! Bucle principal
    do m = 1, nstep
        do j = 1, irate
            ! Cálculos
            y = v
            t = tme
            call fcn(t, y, yprime)
            A = yprime

            y = v + (stpsze*A)/2.0_dp
            t = tme + stpsze/2.0_dp
            call fcn(t, y, yprime)
            B = yprime

            y = v + (stpsze*B)/2.0_dp
            t = tme + stpsze/2.0_dp
            call fcn(t, y, yprime)
            C = yprime

            y = v + (stpsze*C)
            t = tme + stpsze
            call fcn(t, y, yprime)
            D = yprime

            v = v + stpsze*(A + D + 2.0_dp*(B + C))/6.0_dp
            tme = tme + stpsze

            ! Gram-Schmidt Orthonormalization
            ! [Código para la orthonormalización aquí]

            ! Impresión de resultados
            if (mod(m, io) == 0) then
                write(*,'(F12.6, 3(2X, E12.6))') tme, (cum(k)/tme, k=1,n)
            end if
        end do
    end do
end program ode

subroutine fcn(t, y, yprime)
    real(dp), intent(in) :: t
    real(dp), intent(in) :: y(:)
    real(dp), intent(out) :: yprime(:)
    real(dp), parameter :: b=4.0_dp, sg=16.0_dp, r=45.92_dp

    ! Ecuaciones de Lorenz
    yprime(1) = sg*(y(2) - y(1))
    yprime(2) = -y(1)*y(3) + r*y(1) - y(2)
    yprime(3) = y(1)*y(2) - b*y(3)

    ! Ecuaciones linealizadas de Lorenz
    do i = 0, 2
        yprime(4+i) = sg*(y(7+i) - y(4+i))
        yprime(7+i) = (r - y(3))*y(4+i) - y(7+i) - y(1)*y(10+i)
        yprime(10+i) = y(2)*y(4+i) + y(1)*y(7+i) - b*y(10+i)
    end do
end subroutine fcn

