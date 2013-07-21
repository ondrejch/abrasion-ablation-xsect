module binomial_module

implicit none
contains
function binom(n,k) ! n!/((n–k)! k!)
  implicit none
  !real (selected_real_kind (15,307)) :: binom,xn,xk
  real (selected_real_kind (16)) :: binom,xn,xk
  integer, intent(in) :: n, k
  xn=dble(n); xk=dble(k)
  binom = gamma(xn+1)/(gamma(xn-xk+1)*gamma(xk+1))
end function binom

end module binomial_module

