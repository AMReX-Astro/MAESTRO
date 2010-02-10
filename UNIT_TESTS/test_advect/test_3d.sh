#!/bin/sh

MAESTRO_EXE=./main.Linux.Intel.debug.exe
COMPARE_EXE=fcompare.Linux.Intel.exe

${MAESTRO_EXE} inputs_3d_xm &> advect_3d_xm.out
${MAESTRO_EXE} inputs_3d_xp &> advect_3d_xp.out
${MAESTRO_EXE} inputs_3d_ym &> advect_3d_ym.out
${MAESTRO_EXE} inputs_3d_yp &> advect_3d_yp.out
${MAESTRO_EXE} inputs_3d_zm &> advect_3d_zm.out
${MAESTRO_EXE} inputs_3d_zp &> advect_3d_zp.out

echo " " >> advect_3d_report.out
echo "3-d advection -X direction error" >> advect_3d_report.out
${COMPARE_EXE} --infile1 dens_orig --infile2 dens_xm_final >> advect_3d_report.out

echo " " >> advect_3d_report.out
echo "3-d advection +X direction error" >> advect_3d_report.out
${COMPARE_EXE} --infile1 dens_orig --infile2 dens_xp_final >> advect_3d_report.out

echo " " >> advect_3d_report.out
echo "3-d advection -Y direction error" >> advect_3d_report.out
${COMPARE_EXE} --infile1 dens_orig --infile2 dens_ym_final >> advect_3d_report.out

echo " " >> advect_3d_report.out
echo "3-d advection +Y direction error" >> advect_3d_report.out
${COMPARE_EXE} --infile1 dens_orig --infile2 dens_yp_final >> advect_3d_report.out

echo " " >> advect_3d_report.out
echo "3-d advection -Z direction error" >> advect_3d_report.out
${COMPARE_EXE} --infile1 dens_orig --infile2 dens_zm_final >> advect_3d_report.out

echo " " >> advect_3d_report.out
echo "3-d advection +Z direction error" >> advect_3d_report.out
${COMPARE_EXE} --infile1 dens_orig --infile2 dens_zp_final >> advect_3d_report.out


