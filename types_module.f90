module types_module

implicit none
integer, parameter :: sp = selected_real_kind(6,37)   ! single precision
integer, parameter :: dp = selected_real_kind(15,307) ! double precision
integer, parameter :: qp = selected_real_kind(16)     ! quad precision

end module types_module
