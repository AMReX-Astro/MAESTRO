      COMMON /MCCONS/ VERS, PREC, KERR, LENI, LENR

      SAVE /MCCONS/

!$omp threadprivate(/ MCCONS /)
