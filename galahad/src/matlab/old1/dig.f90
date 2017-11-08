program dig
  integer :: i, j, k
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
  real :: val
  REAL( kind = wp ) :: dval
  i = digits( i )
  write( 6, * ) i, ' integer digits '
  j = digits( val )
  write( 6, * ) j, ' real digits '
  k = digits( dval )
  write( 6, * ) k, ' double precision digits '
end program dig
