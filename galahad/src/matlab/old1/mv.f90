      SUBROUTINE mv( A, tot )
      USE GALAHAD_SMT_double
      TYPE ( SMT_type ) :: A
      REAL * 8 :: tot
        tot = sum( A%val( : A%ne ) ) / A%ne
      RETURN
      END SUBROUTINE mv

