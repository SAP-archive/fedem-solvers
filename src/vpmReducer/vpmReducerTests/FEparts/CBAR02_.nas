$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*                  I-DEAS 8 NASTRAN TRANSLATOR
$*                    FOR MSC/NASTRAN VERSION 70.7
$*
$*           MODEL FILE: /project/oeystein/FEDEMTest/test_i01.mf1
$*           INPUT FILE: /project/oeystein/FEDEMTest/CBAR02_.nas
$*             EXPORTED: AT 11:02:31 ON 30-Apr-01
$*                 PART: CBAR02_
$*                  FEM: Fem1
$*
$*                UNITS: SI-Meter (newton)
$*                      ... LENGTH : meter
$*                      ... TIME   : sec
$*                      ... MASS   : kilogram (kg)
$*                      ... FORCE  : newton(N)
$*                      ... TEMPERATURE : deg Celsius
$*
$*        SUBSET EXPORT: OFF
$*
$*     REAL DATA FILTER: ON (    .1000E-14)
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
$*  BULK DATA
$*
$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$*
BEGIN BULK
$*
$*
$*
$*  GRID CARDS
$*
GRID*                  1               0             0.0             0.0+
*                    0.0               0
GRID*                  2               0             0.0             0.0+
*           -1.000000000               0
GRID*                  3               0             0.0             0.0+
*           -0.125000000               0
GRID*                  4               0             0.0             0.0+
*           -0.250000000               0
GRID*                  5               0             0.0             0.0+
*           -0.375000000               0
GRID*                  6               0             0.0             0.0+
*           -0.500000000               0
GRID*                  7               0             0.0             0.0+
*           -0.625000000               0
GRID*                  8               0             0.0             0.0+
*           -0.750000000               0
GRID*                  9               0             0.0             0.0+
*           -0.875000000               0
$*
$*  ELEMENT CARDS
$*
CBAR*                  1               3               1               3+
*                    0.0     1.000000000             0.0
CBAR*                  2               3               3               4+
*                    0.0     1.000000000             0.0
CBAR*                  3               3               4               5+
*                    0.0     1.000000000             0.0
CBAR*                  4               3               5               6+
*                    0.0     1.000000000             0.0
CBAR*                  5               3               6               7+
*                    0.0     1.000000000             0.0
CBAR*                  6               3               7               8+
*                    0.0     1.000000000             0.0
CBAR*                  7               3               8               9+
*                    0.0     1.000000000             0.0
CBAR*                  8               3               9               2+
*                    0.0     1.000000000             0.0
$*
$*  MATERIAL CARDS
$*
$*
$*  I-DEAS Material: 1  name: GENERIC_ISOTROPIC_STEEL
MAT1*                  1 2.067999949E+11 8.015504179E+10     0.289999992+
*         7820.000000000 1.169999996E-05    21.850000381                +
*        1.500000000E+08 1.500000000E+08 6.800000000E+07
$*
$*  PROPERTY CARDS
$*
$*
$*  I-DEAS property: 3  name: LINEAR BEAM3
$*  Fore Section   : 2  name: RECTANGLE 0.2 X 0.1
PBAR*                  3               1     0.020000001 1.666666503E-05+
*        6.666666741E-05 4.577562868E-05             0.0                +
*            0.050000001     0.100000001     0.050000001    -0.100000001+
*           -0.050000001    -0.100000001    -0.050000001     0.100000001+
*            0.833333313     0.833333313
$*
ASET,1,1,2,1
ENDDATA
