module precision_mod
  use, intrinsic :: iso_fortran_env, only: real32, real64
  implicit none
  public

  ! Kinds portáveis (padrão Fortran)
  integer, parameter :: sp = real32   ! Single precision  (32 bits)
  integer, parameter :: dp = real64   ! Double precision  (64 bits)

  ! Kind específico do nvfortran para meia precisão (16 bits)
  ! ATENÇÃO: só funciona com o compilador nvfortran e GPU CC >= 6.0
  integer, parameter :: hp = 2        ! Half precision    (16 bits)

  ! >>> AQUI você escolhe a precisão de trabalho do experimento <<<
  ! Descomente APENAS uma das linhas abaixo:
  ! integer, parameter :: wp = sp
  ! integer, parameter :: wp = dp
  integer, parameter :: wp = hp

end module precision_mod