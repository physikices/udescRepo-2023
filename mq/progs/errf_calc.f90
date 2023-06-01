program calc_errf
  implicit none

  real :: x, sroot, result_1, result_2

  ! entrada de dados
  write(*,*) 'Digite o valor de x:'
  read(*,*) x
  
  sroot=sqrt(x)
  result_1=erf(x)
  result_2=erf(sroot)

  
  ! saida
  write(*,*) 'O valor de erf(',x,') é:', sroot, result_1, result_2
end program calc_errf
