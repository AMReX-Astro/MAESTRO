#!/bin/sh

MAESTRO_EXE=./main.Linux.Intel.debug.exe
COMPARE_EXE=fcompare.Linux.Intel.exe

${MAESTRO_EXE} inputs_2d_xm &> advect_2d_xm.out
${MAESTRO_EXE} inputs_2d_xp &> advect_2d_xp.out
${MAESTRO_EXE} inputs_2d_ym &> advect_2d_ym.out
${MAESTRO_EXE} inputs_2d_yp &> advect_2d_yp.out

echo " " >> advect_2d_report.out
echo "2-d advection -X direction error" >> advect_2d_report.out
${COMPARE_EXE} --infile1 dens_orig --infile2 dens_xm_final >> advect_2d_report.out

echo " " >> advect_2d_report.out
echo "2-d advection +X direction error" >> advect_2d_report.out
${COMPARE_EXE} --infile1 dens_orig --infile2 dens_xp_final >> advect_2d_report.out

echo " " >> advect_2d_report.out
echo "2-d advection -Y direction error" >> advect_2d_report.out
${COMPARE_EXE} --infile1 dens_orig --infile2 dens_ym_final >> advect_2d_report.out

echo " " >> advect_2d_report.out
echo "2-d advection +Y direction error" >> advect_2d_report.out
${COMPARE_EXE} --infile1 dens_orig --infile2 dens_yp_final >> advect_2d_report.out


