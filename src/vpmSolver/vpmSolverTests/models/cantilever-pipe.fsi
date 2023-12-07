&HEADING
  modelFile = 'C:\Users\I331884\Fedem-test\CantileverCylinderBeam.fmm'
  version = 3.0
/

'This input file models a uniform cantilever beam with a circular pipe
'cross section, subjected to a constant shear load at the free end.

&MECHANISM
  id = 1
  extId = 1
  extDescr = 'Cantilever circular pipe with constant transverse tip load'
  weightTranslation =  1.000000000e-01
  weightRotation    =  1.000000000e+00
  weightGeneralized =  1.000000000e+00
/

&TRIAD
  id = 16
  extId = 1
  extDescr = 'Fixed end'
  nDOFs = 6
  ur  =   1.000000000e+00   0.000000000e+00   0.000000000e+00   0.000000000e+00
          0.000000000e+00   1.000000000e+00   0.000000000e+00   0.000000000e+00
          0.000000000e+00   0.000000000e+00   1.000000000e+00   0.000000000e+00
  BC = 0 0 0 0 0 0
/

&TRIAD
  id = 17
  extId = 2
  extDescr = 'Free end'
  nDOFs = 6
  ur  =   1.000000000e+00   0.000000000e+00   0.000000000e+00   0.000000000e+00
          0.000000000e+00   1.000000000e+00   0.000000000e+00   0.000000000e+00
          0.000000000e+00   0.000000000e+00   1.000000000e+00   1.000000000e+01
/

&TRIAD
  id = 22
  extId = 3
  nDOFs = 6
  ur  =   1.000000000e+00   0.000000000e+00   0.000000000e+00   0.000000000e+00
          0.000000000e+00   1.000000000e+00   0.000000000e+00   0.000000000e+00
          0.000000000e+00   0.000000000e+00   1.000000000e+00   2.000000000e+00
/

&TRIAD
  id = 23
  extId = 4
  nDOFs = 6
  ur  =   1.000000000e+00   0.000000000e+00   0.000000000e+00   0.000000000e+00
          0.000000000e+00   1.000000000e+00   0.000000000e+00   0.000000000e+00
          0.000000000e+00   0.000000000e+00   1.000000000e+00   4.000000000e+00
/

&TRIAD
  id = 25
  extId = 5
  nDOFs = 6
  ur  =   1.000000000e+00   0.000000000e+00   0.000000000e+00   0.000000000e+00
          0.000000000e+00   1.000000000e+00   0.000000000e+00   0.000000000e+00
          0.000000000e+00   0.000000000e+00   1.000000000e+00   6.000000000e+00
/

&TRIAD
  id = 27
  extId = 6
  nDOFs = 6
  ur  =   1.000000000e+00   0.000000000e+00   0.000000000e+00   0.000000000e+00
          0.000000000e+00   1.000000000e+00   0.000000000e+00   0.000000000e+00
          0.000000000e+00   0.000000000e+00   1.000000000e+00   8.000000000e+00
/

&SUP_EL
  id = 18
  extId = 1
  numTriads = 2
  triadIds = 16 22
  elPropId = 64
  shadowPosAlg = 1
  refTriad1Id = 16, offset1 =  0.000000000e+00   0.000000000e+00   0.000000000e+00
  refTriad2Id = 22, offset2 =  0.000000000e+00   0.000000000e+00   0.000000000e+00
  refTriad3Id = 16, offset3 =  0.000000000e+00   2.000000000e+00   0.000000000e+00
  massCorrFlag = 0
  stiffScale =  1.000000000e+00
  massScale  =  1.000000000e+00
  alpha1 =  0.000000000e+00,  alpha2 =  0.000000000e+00
  supPos =  0.000000000e+00  -0.000000000e+00  -1.000000000e+00   0.000000000e+00
            0.000000000e+00   1.000000000e+00   0.000000000e+00   0.000000000e+00
            1.000000000e+00  -0.000000000e+00   0.000000000e+00   0.000000000e+00
/
&TRIAD_UNDPOS
  supElId = 18
  triadId = 16
  undPosInSupElSystem =  0.000000000e+00   0.000000000e+00  -1.000000000e+00   0.000000000e+00
                         0.000000000e+00   1.000000000e+00   0.000000000e+00   0.000000000e+00
                         1.000000000e+00   0.000000000e+00   0.000000000e+00   0.000000000e+00
/
&TRIAD_UNDPOS
  supElId = 18
  triadId = 22
  undPosInSupElSystem =  0.000000000e+00   0.000000000e+00  -1.000000000e+00   2.000000000e+00
                         0.000000000e+00   1.000000000e+00   0.000000000e+00   0.000000000e+00
                         1.000000000e+00   0.000000000e+00   0.000000000e+00   0.000000000e+00
/

&SUP_EL
  id = 24
  extId = 2
  numTriads = 2
  triadIds = 22 23
  elPropId = 64
  shadowPosAlg = 1
  refTriad1Id = 22, offset1 =  0.000000000e+00   0.000000000e+00   0.000000000e+00
  refTriad2Id = 23, offset2 =  0.000000000e+00   0.000000000e+00   0.000000000e+00
  refTriad3Id = 22, offset3 =  0.000000000e+00   2.000000000e+00   0.000000000e+00
  massCorrFlag = 0
  stiffScale =  1.000000000e+00
  massScale  =  1.000000000e+00
  alpha1 =  0.000000000e+00,  alpha2 =  0.000000000e+00
  supPos =  0.000000000e+00  -0.000000000e+00  -1.000000000e+00   0.000000000e+00
            0.000000000e+00   1.000000000e+00   0.000000000e+00   0.000000000e+00
            1.000000000e+00  -0.000000000e+00   0.000000000e+00   2.000000000e+00
/
&TRIAD_UNDPOS
  supElId = 24
  triadId = 22
  undPosInSupElSystem =  0.000000000e+00   0.000000000e+00  -1.000000000e+00   0.000000000e+00
                         0.000000000e+00   1.000000000e+00   0.000000000e+00   0.000000000e+00
                         1.000000000e+00   0.000000000e+00   0.000000000e+00   0.000000000e+00
/
&TRIAD_UNDPOS
  supElId = 24
  triadId = 23
  undPosInSupElSystem =  0.000000000e+00   0.000000000e+00  -1.000000000e+00   2.000000000e+00
                         0.000000000e+00   1.000000000e+00   0.000000000e+00   0.000000000e+00
                         1.000000000e+00   0.000000000e+00   0.000000000e+00   0.000000000e+00
/

&SUP_EL
  id = 26
  extId = 3
  numTriads = 2
  triadIds = 23 25
  elPropId = 64
  shadowPosAlg = 1
  refTriad1Id = 23, offset1 =  0.000000000e+00   0.000000000e+00   0.000000000e+00
  refTriad2Id = 25, offset2 =  0.000000000e+00   0.000000000e+00   0.000000000e+00
  refTriad3Id = 23, offset3 =  0.000000000e+00   2.000000000e+00   0.000000000e+00
  massCorrFlag = 0
  stiffScale =  1.000000000e+00
  massScale  =  1.000000000e+00
  alpha1 =  0.000000000e+00,  alpha2 =  0.000000000e+00
  supPos =  0.000000000e+00  -0.000000000e+00  -1.000000000e+00   0.000000000e+00
            0.000000000e+00   1.000000000e+00   0.000000000e+00   0.000000000e+00
            1.000000000e+00  -0.000000000e+00   0.000000000e+00   4.000000000e+00
/
&TRIAD_UNDPOS
  supElId = 26
  triadId = 23
  undPosInSupElSystem =  0.000000000e+00   0.000000000e+00  -1.000000000e+00   0.000000000e+00
                         0.000000000e+00   1.000000000e+00   0.000000000e+00   0.000000000e+00
                         1.000000000e+00   0.000000000e+00   0.000000000e+00   0.000000000e+00
/
&TRIAD_UNDPOS
  supElId = 26
  triadId = 25
  undPosInSupElSystem =  0.000000000e+00   0.000000000e+00  -1.000000000e+00   2.000000000e+00
                         0.000000000e+00   1.000000000e+00   0.000000000e+00   0.000000000e+00
                         1.000000000e+00   0.000000000e+00   0.000000000e+00   0.000000000e+00
/

&SUP_EL
  id = 28
  extId = 4
  numTriads = 2
  triadIds = 25 27
  elPropId = 64
  shadowPosAlg = 1
  refTriad1Id = 25, offset1 =  0.000000000e+00   0.000000000e+00   0.000000000e+00
  refTriad2Id = 27, offset2 =  0.000000000e+00   0.000000000e+00   0.000000000e+00
  refTriad3Id = 25, offset3 =  0.000000000e+00   2.000000000e+00   0.000000000e+00
  massCorrFlag = 0
  stiffScale =  1.000000000e+00
  massScale  =  1.000000000e+00
  alpha1 =  0.000000000e+00,  alpha2 =  0.000000000e+00
  supPos =  0.000000000e+00  -0.000000000e+00  -1.000000000e+00   0.000000000e+00
            0.000000000e+00   1.000000000e+00   0.000000000e+00   0.000000000e+00
            1.000000000e+00  -0.000000000e+00   0.000000000e+00   6.000000000e+00
/
&TRIAD_UNDPOS
  supElId = 28
  triadId = 25
  undPosInSupElSystem =  0.000000000e+00   0.000000000e+00  -1.000000000e+00   0.000000000e+00
                         0.000000000e+00   1.000000000e+00   0.000000000e+00   0.000000000e+00
                         1.000000000e+00   0.000000000e+00   0.000000000e+00   0.000000000e+00
/
&TRIAD_UNDPOS
  supElId = 28
  triadId = 27
  undPosInSupElSystem =  0.000000000e+00   0.000000000e+00  -1.000000000e+00   2.000000000e+00
                         0.000000000e+00   1.000000000e+00   0.000000000e+00   0.000000000e+00
                         1.000000000e+00   0.000000000e+00   0.000000000e+00   0.000000000e+00
/

&SUP_EL
  id = 29
  extId = 5
  numTriads = 2
  triadIds = 27 17
  elPropId = 64
  shadowPosAlg = 1
  refTriad1Id = 27, offset1 =  0.000000000e+00   0.000000000e+00   0.000000000e+00
  refTriad2Id = 17, offset2 =  0.000000000e+00   0.000000000e+00   0.000000000e+00
  refTriad3Id = 27, offset3 =  0.000000000e+00   2.000000000e+00   0.000000000e+00
  massCorrFlag = 0
  stiffScale =  1.000000000e+00
  massScale  =  1.000000000e+00
  alpha1 =  0.000000000e+00,  alpha2 =  0.000000000e+00
  supPos =  0.000000000e+00  -0.000000000e+00  -1.000000000e+00   0.000000000e+00
            0.000000000e+00   1.000000000e+00   0.000000000e+00   0.000000000e+00
            1.000000000e+00  -0.000000000e+00   0.000000000e+00   8.000000000e+00
/
&TRIAD_UNDPOS
  supElId = 29
  triadId = 27
  undPosInSupElSystem =  0.000000000e+00   0.000000000e+00  -1.000000000e+00   0.000000000e+00
                         0.000000000e+00   1.000000000e+00   0.000000000e+00   0.000000000e+00
                         1.000000000e+00   0.000000000e+00   0.000000000e+00   0.000000000e+00
/
&TRIAD_UNDPOS
  supElId = 29
  triadId = 17
  undPosInSupElSystem =  0.000000000e+00   0.000000000e+00  -1.000000000e+00   2.000000000e+00
                         0.000000000e+00   1.000000000e+00   0.000000000e+00   0.000000000e+00
                         1.000000000e+00   0.000000000e+00   0.000000000e+00   0.000000000e+00
/

&ELEMENT_PROPERTY
  id = 64
  extId = 1
  geometry = 1.539380400e-02 4.621989652e-04 4.621989652e-04 9.243979304e-04 2.0 2.0 0.0 0.0
  material = 7.850000000e+03 2.100000000e+11 8.139534884e+10
/

&LOAD
  id = 19
  extId = 1
  triadId = 17
  lDof = 1
  f0 = 100000.0
/
