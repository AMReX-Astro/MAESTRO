pgf90 -o tdaxpy_gpu \
  -fast -acc -ta=tesla -Mcuda=cuda7.0 -Minfo=accel \
  blas_mod.F90 tdaxpy.F90

pgf90 -o tdaxpy_cpu \
  -fast \
  blas_mod.F90 tdaxpy.F90
