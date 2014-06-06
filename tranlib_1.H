      COMMON /MCMCMC/ RU, PATMOS, SMALL, NKK, NO, NLITE, INLIN, IKTDIF,
     1                IPVT, NWT, NEPS, NSIG, NDIP, NPOL, NZROT, NLAM,
     2                NETA, NDIF, NTDIF, NXX, NVIS, NXI, NCP, NCROT,
     3                NCINT, NBIND, NEOK, NSGM, NAST, NBST, NCST, NXL,
     4                NR, NWRK, K3

      SAVE /MCMCMC/

!$omp threadprivate(/ MCMCMC /)
