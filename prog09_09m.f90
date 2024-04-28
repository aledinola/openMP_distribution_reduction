!##############################################################################
! MODULE globals
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
module globals

    use omp_lib ! OpenMP library
    implicit none

    integer, parameter :: par = 0 ! 1: parallel, 0: serial
    
    ! model parameters
    real*8, parameter :: gamma = 0.50d0
    real*8, parameter :: beta = 0.96d0
    real*8, parameter :: alpha = 0.36d0
    real*8, parameter :: delta = 0.08d0
    real*8, parameter :: rho = 0.6d0
    real*8, parameter :: sigma_eps = 0.04d0*(1d0-rho**2)

    ! numerical parameters
    real*8, parameter :: a_l = 0d0
    real*8, parameter :: a_u = 100d0
    real*8, parameter :: a_grow = 0.01d0
    real*8, parameter :: sig_in = 1d-7
    real*8, parameter :: sig_out = 1d-6
    integer, parameter :: itermax = 50000

    ! macroeconomic variables
    real*8 :: r, w, KK, AA, LL, YY, CC, II

    ! the shock process
    integer, parameter :: NS = 7
    real*8 :: pi(NS, NS), eta(NS), weights(NS)

    ! policy and value function
    integer, parameter :: NA = 1000
    real*8 :: a(0:NA), c(0:NA, NS), V(0:NA, NS)

    ! variables to numerically determine policy function
    real*8 :: c_new(0:NA, NS), RHS(0:NA, NS)
    logical :: check

    ! variables to communicate with function
    real*8 :: a_com
    integer :: is_com

    ! variables to numerically determine the steady state distribution
    real*8 :: phi(0:NA, NS), phi_new(0:NA, NS)

contains


    ! calculate the difference on the asset market given an interest rate
    function asset_market(r_input)

        implicit none
        real*8, intent(in) :: r_input
        real*8 :: asset_market, t1, t2

        ! set the interest rate
        r = r_input

        ! calculate the corresponding capital stock
        KK = (alpha/(r+delta))**(1d0/(1d0-alpha))*LL

        ! get wages and output
        w = (1d0-alpha)*(KK/LL)**alpha
        YY = KK**alpha*LL**(1d0-alpha)

        ! solve the household value function iteration problem
        call solve_household()

        ! calculate the invariant distribution directly
        t1 = omp_get_wtime()
        !call get_distribution()
        call get_distribution2()
        t2 = omp_get_wtime()
        write(*,'(A,G0)') "Time for distribution: ", t2-t1
        
        ! get aggregate assets
        AA = sum(a*sum(phi, 2))

        ! get the difference between asset demand and supply
        asset_market = (AA - KK)/AA

        write(*, '(a)')'       K/Y         r      diff_K'
        write(*,'(2f10.4,f12.8)')KK/YY*100d0, r*100d0, AA-KK

    end function asset_market


    ! determines the solution to the household optimization problem
    subroutine solve_household()

        use toolbox

        implicit none
        integer :: ia, is, iter
        real*8 :: con_lev, p

        ! set tolerance level for interior optimization
        call settol_root(sig_in)

        ! do a policy function iteration with rootfinding
        do iter = 1, itermax

            ! interpolate coefficients
            call interpolate()

            ! solve the household problem for all gridpoints
            do ia = 0, NA
                do is = 1, NS

                    ! initialize starting value and communication variables
                    p = c(ia, is)
                    a_com = a(ia)
                    is_com = is

                    ! find the optimal consumption level
                    call fzero(p, foc, check)
                    if(check)write(*,*)'ERROR IN ROOTFINDING PROCESS'

                    ! check for borrowing constraint
                    if(p > (1d0+r)*a_com+w*eta(is_com) - a_l) &
                       p = (1d0+r)*a_com+w*eta(is_com) - a_l

                    ! get optimal consumption and value function
                    c_new(ia, is) = p

                enddo
            enddo

            ! get convergence level
            con_lev = maxval(abs(c_new(:, :) - c(:, :))/max(abs(c(:, :)), 1d-10))

            ! check for convergence
            if(con_lev < sig_in)then
                return
            endif

            c = c_new
        enddo

        write(*,*)'Policy Function Iteration did not converge'

    end subroutine solve_household


    ! the first order condition
    function foc(x_in)

        use toolbox

        implicit none
        real*8, intent(in) :: x_in
        real*8 :: foc, aplus, varphi
        integer :: ial, iar

        ! future assets
        aplus = (1d0+r)*a_com + w*eta(is_com) - x_in

        ! calculate future expected marginal utility
        call linint_Grow(aplus, a_l, a_u, a_grow, NA, ial, iar, varphi)
        foc = varphi*RHS(ial,is_com) + (1d0-varphi)*RHS(iar,is_com)

        ! get first order condition
        foc = x_in - foc

    end function foc


    ! Set up data for interpolation
    subroutine interpolate()

        implicit none
        integer :: is, ia, is_p

        do is = 1, NS
            do ia = 0, NA

                ! calculate the RHS of the first order condition
                RHS(ia, is) = 0d0
                do is_p = 1, NS
                    RHS(ia, is) = RHS(ia, is) + pi(is, is_p)*c_new(ia, is_p)**(-1d0/gamma)
                enddo
                RHS(ia, is) = (beta*(1d0+r)*RHS(ia, is))**(-gamma)
            enddo
        enddo

    end subroutine interpolate


    ! calculates the invariant distribution of households
    subroutine get_distribution()

        use toolbox
        use omp_lib

        implicit none
        integer :: ia, is, iter, is_p
        real*8 :: aplus, con_lev,phi_val
        integer :: ial(0:NA, NS), iar(0:NA, NS)
        real*8 :: varphi(0:NA, NS)

        ! get the interpolation shares and points
        do ia = 0, NA
            do is = 1, NS

                ! calculate where this guy would go
                aplus = (1d0+r)*a(ia) + w*eta(is) - c(ia,is)
                aplus = min(max(aplus, a_l), a_u)

                ! determine the gridpoints in between this decision lies
                call linint_Grow(aplus, a_l, a_u, a_grow, NA, &
                    ial(ia, is), iar(ia, is), varphi(ia, is))
            enddo
        enddo

        ! iterate until the distribution function converges
        phi_val = 1d0/((NA+1)*NS) ! initial guess
        phi =phi_val ! set initial guess
        write(*,*) "sum(phi) = ", sum(phi)
        
        do iter = 1, itermax

            phi_new = 0d0

            !$omp parallel if (par==1) default(shared) private(ia,is,is_p) reduction(+: phi_new)
            !$omp do collapse(2)
            do ia = 0, NA
                do is = 1, NS
                    do is_p = 1, NS
                        phi_new(ial(ia,is), is_p) = phi_new(ial(ia,is), is_p) + &
                            pi(is,is_p)*varphi(ia,is)*phi(ia,is)
                        phi_new(iar(ia,is), is_p) = phi_new(iar(ia,is), is_p) + &
                            pi(is,is_p)*(1d0-varphi(ia,is))*phi(ia,is)
                    enddo
                enddo
            enddo
            !$omp enddo
            !$omp end parallel
            
            con_lev = maxval(abs(phi_new(:, :) - phi(:, :))/max(abs(phi(:, :)), 1d-10))

            ! update distribution
            phi = phi_new
            
            write(*,*) "iter = ", iter, "con_lev = ", con_lev

            ! check for convergence
            if(con_lev < sig_in)then
                return
            endif
        enddo !end iteration
        
        
        write(*,*)'Distribution Function Iteration did not converge'

    end subroutine get_distribution

    ! calculates the invariant distribution of households
    subroutine get_distribution2()

        use toolbox
        use omp_lib

        implicit none
        integer :: ia, is, iter, is_p
        real*8 :: aplus, con_lev
        integer :: ial(0:NA, NS), iar(0:NA, NS)
        real*8 :: varphi(0:NA, NS)

        ! get the interpolation shares and points
        do ia = 0, NA
            do is = 1, NS

                ! calculate where this guy would go
                aplus = (1d0+r)*a(ia) + w*eta(is) - c(ia,is)
                aplus = min(max(aplus, a_l), a_u)

                ! determine the gridpoints in between this decision lies
                call linint_Grow(aplus, a_l, a_u, a_grow, NA, &
                    ial(ia, is), iar(ia, is), varphi(ia, is))
            enddo
        enddo

        ! iterate until the distribution function converges
        phi =  1d0/((NA+1)*NS) ! set initial guess
        write(*,*) "sum(phi) = ", sum(phi)
        
        do iter = 1, itermax

            phi_new = 0d0

            ! First, calculate phi_new(a',s) using the law of motion for the 
            ! endogenous state variable a. 
            do ia = 0, NA
            do is = 1, NS
                phi_new(ial(ia,is),is)=phi_new(ial(ia,is),is)+varphi(ia,is)*phi(ia,is)
                phi_new(iar(ia,is),is)=phi_new(iar(ia,is),is)+(1d0-varphi(ia,is))*phi(ia,is)
            enddo
            enddo
            ! Second, calculate phi_new(a',s') using the transition matrix pi(s,s')
            ! phi_new(a',s') = phi_new(a',s)*pi(s,s'), where * is a matrix product
            phi_new = matmul(phi_new, pi)
            
            con_lev = maxval(abs(phi_new(:, :) - phi(:, :))/max(abs(phi(:, :)), 1d-10))

            ! update distribution
            phi = phi_new
            
            write(*,*) "iter = ", iter, "con_lev = ", con_lev

            ! check for convergence
            if(con_lev < sig_in)then
                return
            endif
        enddo !end iteration
        
        
        write(*,*)'Distribution Function Iteration did not converge'

    end subroutine get_distribution2

end module
