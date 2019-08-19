species(
    label = 'C=C([O])C1[CH]CC=CC=C1(2238)',
    structure = SMILES('C=C([O])C1[CH]CC=CC=C1'),
    E0 = (262.341,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.04339,0.10061,-8.92457e-05,4.06621e-08,-7.38355e-12,31742.4,20.7418], Tmin=(100,'K'), Tmax=(1326.26,'K')), NASAPolynomial(coeffs=[21.3094,0.0331933,-1.29969e-05,2.33409e-09,-1.58673e-13,25813.3,-93.409], Tmin=(1326.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(262.341,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3-Cycloheptadiene) + radical(Cs_S) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C7H8(11)',
    structure = SMILES('C1=CC=CCC=C1'),
    E0 = (165.503,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3872.77,'J/mol'), sigma=(6.33517,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=604.92 K, Pc=34.56 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.976321,0.0470646,2.43918e-05,-7.34711e-08,3.37246e-11,20032,17.1406], Tmin=(100,'K'), Tmax=(948.557,'K')), NASAPolynomial(coeffs=[19.6256,0.016921,-4.63468e-06,8.3158e-10,-6.48659e-14,14312.2,-83.3476], Tmin=(948.557,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(165.503,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene)"""),
)

species(
    label = 'CH2CO(43)',
    structure = SMILES('C=C=O'),
    E0 = (-60.8183,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2120,512.5,787.5],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13241,0.0181319,-1.74093e-05,9.35336e-09,-2.01725e-12,-7148.09,13.3808], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.75871,0.00635124,-2.25955e-06,3.62322e-10,-2.15856e-14,-8085.33,-4.9649], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-60.8183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CH2CO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH2]C1([O])C2C=CC=CCC21(3450)',
    structure = SMILES('[CH2]C1([O])C2C=CC=CCC21'),
    E0 = (400.583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.255213,0.11318,-0.000181617,1.97262e-07,-8.36439e-11,48313,11.8192], Tmin=(100,'K'), Tmax=(790.348,'K')), NASAPolynomial(coeffs=[-4.65938,0.0889336,-4.728e-05,9.44808e-09,-6.69903e-13,50462.6,41.2255], Tmin=(790.348,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(400.583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s2_3_7_ane) - ring(Cycloheptane) - ring(Cyclopropane) + ring(1,3-Cycloheptadiene) + ring(Cyclopropane) + radical(CJC(C)2O) + radical(CC(C)2OJ)"""),
)

species(
    label = 'H(19)',
    structure = SMILES('[H]'),
    E0 = (211.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (1.00794,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1205.6,'J/mol'), sigma=(2.05,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,25472.7,-0.459566], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,25472.7,-0.459566], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.792,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'C=C([O])C1C=CC=CCC=1(3515)',
    structure = SMILES('C=C([O])C1C=CC=CCC=1'),
    E0 = (160.288,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (133.167,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.58682,0.0839792,-3.59616e-05,-2.21476e-08,1.6973e-11,19458.5,27.2251], Tmin=(100,'K'), Tmax=(966.278,'K')), NASAPolynomial(coeffs=[22.6556,0.0272367,-9.15148e-06,1.63023e-09,-1.16433e-13,13124.1,-93.6442], Tmin=(966.278,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.288,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3,5-Cycloheptatriene) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])C1C=CC=CC=C1(3463)',
    structure = SMILES('C=C([O])C1C=CC=CC=C1'),
    E0 = (164.955,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (133.167,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.930506,0.0933263,-6.03638e-05,-2.68579e-10,1.02316e-11,20030.7,28.3633], Tmin=(100,'K'), Tmax=(960.374,'K')), NASAPolynomial(coeffs=[23.9698,0.024692,-7.95062e-06,1.37861e-09,-9.7214e-14,13630.4,-99.1814], Tmin=(960.374,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(164.955,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3,5-Cycloheptatriene) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2][C]=O(360)',
    structure = SMILES('[CH2][C]=O'),
    E0 = (160.185,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,539.612,539.669],'cm^-1')),
        HinderedRotor(inertia=(0.000578908,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39563,0.0101365,2.30741e-06,-8.97567e-09,3.68242e-12,19290.3,10.0703], Tmin=(100,'K'), Tmax=(1068.9,'K')), NASAPolynomial(coeffs=[6.35055,0.00638951,-2.69368e-06,5.4221e-10,-4.02477e-14,18240.9,-6.33602], Tmin=(1068.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CJC=O) + radical(CsCJ=O)"""),
)

species(
    label = '[CH]1C=C[CH]CC=C1(125)',
    structure = SMILES('[CH]1C=C[CH]CC=C1'),
    E0 = (333.186,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.70798,-0.040239,0.000339088,-4.54205e-07,1.84685e-10,40184,19.5412], Tmin=(100,'K'), Tmax=(904.636,'K')), NASAPolynomial(coeffs=[40.627,-0.029508,2.54911e-05,-5.11056e-09,3.33198e-13,26023.7,-199.942], Tmin=(904.636,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(333.186,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(Allyl_S) + radical(C=CCJC=C)"""),
)

species(
    label = 'C=C([O])[C]1C=CC=CCC1(3465)',
    structure = SMILES('C=C([O])C1=C[CH]C=CCC1'),
    E0 = (166.192,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.354442,0.0124675,0.000257393,-3.96227e-07,1.69594e-10,20183.2,31.964], Tmin=(100,'K'), Tmax=(901.444,'K')), NASAPolynomial(coeffs=[49.1487,-0.0259374,2.49216e-05,-5.11515e-09,3.38541e-13,4149.51,-238.517], Tmin=(901.444,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(166.192,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(C=C(C)OJ) + radical(C=CCJC=C)"""),
)

species(
    label = 'C=C([O])C1[CH]C=CC=CC1(2239)',
    structure = SMILES('C=C([O])C1C=C[CH]C=CC1'),
    E0 = (187.021,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,350,440,435,1725,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4629.19,'J/mol'), sigma=(7.55705,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=723.07 K, Pc=24.34 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.554994,0.0115722,0.000248281,-3.78125e-07,1.60819e-10,22677.8,32.2663], Tmin=(100,'K'), Tmax=(903.379,'K')), NASAPolynomial(coeffs=[45.8044,-0.0203488,2.1608e-05,-4.45426e-09,2.92929e-13,7629.38,-219.478], Tmin=(903.379,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(187.021,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(C=C(C)OJ) + radical(C=CCJC=C)"""),
)

species(
    label = 'C=C([O])C1[C]=CC=CCC1(3467)',
    structure = SMILES('C=C([O])C1[C]=CC=CCC1'),
    E0 = (305.64,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.19627,0.105658,-0.000100276,4.91801e-08,-9.58357e-12,36954.4,19.2041], Tmin=(100,'K'), Tmax=(1243.31,'K')), NASAPolynomial(coeffs=[21.3316,0.0331805,-1.28352e-05,2.29407e-09,-1.55895e-13,31352.6,-94.3856], Tmin=(1243.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(305.64,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3-Cycloheptadiene) + radical(C=C(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'C=C([O])C1C=CC=[C]CC1(3464)',
    structure = SMILES('C=C([O])C1C=CC=[C]CC1'),
    E0 = (305.64,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.19627,0.105658,-0.000100276,4.91801e-08,-9.58357e-12,36954.4,19.2041], Tmin=(100,'K'), Tmax=(1243.31,'K')), NASAPolynomial(coeffs=[21.3316,0.0331805,-1.28352e-05,2.29407e-09,-1.55895e-13,31352.6,-94.3856], Tmin=(1243.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(305.64,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3-Cycloheptadiene) + radical(Cds_S) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C(O)[C]1[CH]CC=CC=C1(3516)',
    structure = SMILES('C=C(O)C1[CH]CC=C[CH]C=1'),
    E0 = (169.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.189035,0.011177,0.000274061,-4.22189e-07,1.80836e-10,20591.7,29.9504], Tmin=(100,'K'), Tmax=(902.101,'K')), NASAPolynomial(coeffs=[52.9572,-0.0318283,2.80214e-05,-5.68835e-09,3.75586e-13,3300.64,-262.258], Tmin=(902.101,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(169.5,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(C=CCJC=C) + radical(Allyl_S)"""),
)

species(
    label = '[CH]=C(O)C1[CH]CC=CC=C1(3517)',
    structure = SMILES('[CH]=C(O)C1[CH]CC=CC=C1'),
    E0 = (371.632,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,350,440,435,1725,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.68754,0.110478,-0.000106539,5.13746e-08,-9.66786e-12,44914.1,21.9583], Tmin=(100,'K'), Tmax=(1300.51,'K')), NASAPolynomial(coeffs=[25.6475,0.0264027,-9.56672e-06,1.66425e-09,-1.11861e-13,37804.2,-117.1], Tmin=(1300.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(371.632,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3-Cycloheptadiene) + radical(Cs_S) + radical(Cds_P)"""),
)

species(
    label = 'C=C(O)C1[C]=CC=CC[CH]1(3518)',
    structure = SMILES('C=C(O)C1[C]=CC=CC[CH]1'),
    E0 = (362.378,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3615,1277.5,1000,350,440,435,1725,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.41087,0.108849,-0.000105547,5.18876e-08,-1.00491e-11,43787.3,21.0311], Tmin=(100,'K'), Tmax=(1256.42,'K')), NASAPolynomial(coeffs=[23.2716,0.0302683,-1.17312e-05,2.10815e-09,-1.43995e-13,37585,-103.682], Tmin=(1256.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(362.378,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3-Cycloheptadiene) + radical(Cs_S) + radical(Cds_S)"""),
)

species(
    label = 'C=C([O])C1C=[C]C=CCC1(3468)',
    structure = SMILES('C=C([O])C1C=[C]C=CCC1'),
    E0 = (266.794,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,350,440,435,1725,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.24934,0.107725,-0.000105173,5.38173e-08,-1.09522e-11,32283.4,18.5159], Tmin=(100,'K'), Tmax=(1193.57,'K')), NASAPolynomial(coeffs=[20.6714,0.0342624,-1.2849e-05,2.24972e-09,-1.50995e-13,27050.7,-91.1175], Tmin=(1193.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(266.794,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3-Cycloheptadiene) + radical(C=CJC=C) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])C1C=C[C]=CCC1(3466)',
    structure = SMILES('C=C([O])C1C=C[C]=CCC1'),
    E0 = (266.794,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.24934,0.107725,-0.000105173,5.38173e-08,-1.09522e-11,32283.4,18.5159], Tmin=(100,'K'), Tmax=(1193.57,'K')), NASAPolynomial(coeffs=[20.6714,0.0342624,-1.2849e-05,2.24972e-09,-1.50995e-13,27050.7,-91.1175], Tmin=(1193.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(266.794,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3-Cycloheptadiene) + radical(C=C(C)OJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=C([O])C1C=CC=CCC1(3469)',
    structure = SMILES('[CH]=C([O])C1C=CC=CCC1'),
    E0 = (314.895,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,350,440,435,1725,3120,650,792.5,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.46759,0.107233,-0.00010112,4.85174e-08,-9.15322e-12,38080.9,20.1114], Tmin=(100,'K'), Tmax=(1292.08,'K')), NASAPolynomial(coeffs=[23.6978,0.0293259,-1.06749e-05,1.85075e-09,-1.23787e-13,31577.8,-107.746], Tmin=(1292.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(314.895,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3-Cycloheptadiene) + radical(Cds_P) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C(O)C1[CH]C=C[CH]C=C1(3422)',
    structure = SMILES('C=C(O)C1[CH]C=C[CH]C=C1'),
    E0 = (190.329,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3615,1277.5,1000,350,440,435,1725,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.390068,0.0102755,0.000264974,-4.04124e-07,1.72079e-10,23086.2,30.251], Tmin=(100,'K'), Tmax=(903.921,'K')), NASAPolynomial(coeffs=[49.6123,-0.0262385,2.47071e-05,-5.02728e-09,3.29959e-13,6780.75,-243.216], Tmin=(903.921,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(190.329,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(Allyl_S) + radical(C=CCJC=C)"""),
)

species(
    label = 'C=C(O)C1[CH]CC=C[C]=C1(3519)',
    structure = SMILES('C=C(O)C1[CH]CC=C[C]=C1'),
    E0 = (323.531,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3615,1277.5,1000,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.45813,0.110856,-0.000110262,5.63294e-08,-1.13491e-11,39116.1,20.3215], Tmin=(100,'K'), Tmax=(1210.13,'K')), NASAPolynomial(coeffs=[22.5813,0.0313953,-1.17686e-05,2.06897e-09,-1.39499e-13,33297.9,-100.24], Tmin=(1210.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(323.531,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3-Cycloheptadiene) + radical(C=CJC=C) + radical(Cs_S)"""),
)

species(
    label = 'C=C(O)C1[CH]C[C]=CC=C1(3520)',
    structure = SMILES('C=C(O)C1[CH]C[C]=CC=C1'),
    E0 = (362.378,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3615,1277.5,1000,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.41087,0.108849,-0.000105547,5.18876e-08,-1.00491e-11,43787.3,21.0311], Tmin=(100,'K'), Tmax=(1256.42,'K')), NASAPolynomial(coeffs=[23.2716,0.0302683,-1.17312e-05,2.10815e-09,-1.43995e-13,37585,-103.682], Tmin=(1256.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(362.378,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3-Cycloheptadiene) + radical(Cs_S) + radical(Cds_S)"""),
)

species(
    label = 'C=C(O)C1[CH]CC=[C]C=C1(3521)',
    structure = SMILES('C=C(O)C1[CH]CC=[C]C=C1'),
    E0 = (323.531,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3615,1277.5,1000,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.45813,0.110856,-0.000110262,5.63294e-08,-1.13491e-11,39116.1,20.3215], Tmin=(100,'K'), Tmax=(1210.13,'K')), NASAPolynomial(coeffs=[22.5813,0.0313953,-1.17686e-05,2.06897e-09,-1.39499e-13,33297.9,-100.24], Tmin=(1210.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(323.531,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3-Cycloheptadiene) + radical(Cs_S) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]1CC=CC=CC1[C]1CO1(3522)',
    structure = SMILES('[CH]1CC=CC=CC1[C]1CO1'),
    E0 = (409.058,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.506127,0.0883096,-5.21041e-05,-1.69808e-09,9.84839e-12,49370.5,19.6746], Tmin=(100,'K'), Tmax=(921.133,'K')), NASAPolynomial(coeffs=[18.3072,0.0348247,-1.09486e-05,1.76533e-09,-1.1569e-13,44707.7,-76.0402], Tmin=(921.133,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(409.058,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(Ethylene_oxide) + ring(1,3-Cycloheptadiene) + radical(Cs_S) + radical(C2CsJO)"""),
)

species(
    label = 'C=C([O])C1[CH]C2C=CCC21(3523)',
    structure = SMILES('C=C([O])C1[CH]C2C=CCC21'),
    E0 = (308.147,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.649363,0.0520033,4.24263e-05,-8.8877e-08,3.64317e-11,37201.4,30.7217], Tmin=(100,'K'), Tmax=(990.148,'K')), NASAPolynomial(coeffs=[17.3889,0.0368407,-1.40793e-05,2.67925e-09,-1.95974e-13,31314.9,-62.8573], Tmin=(990.148,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.147,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + polycyclic(s2_4_5_ene_1) + radical(cyclobutane) + radical(C=C(C)OJ)"""),
)

species(
    label = 'O=C1CC2CC=C[CH][CH]C12(3125)',
    structure = SMILES('O=C1CC2C[CH][CH]C=CC12'),
    E0 = (143.144,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.114083,0.0936147,-8.5959e-05,4.92022e-08,-1.28079e-11,17349.2,-1.80427], Tmin=(100,'K'), Tmax=(879.249,'K')), NASAPolynomial(coeffs=[7.53504,0.059856,-2.83694e-05,5.53858e-09,-3.93488e-13,16044.2,-36.6515], Tmin=(879.249,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(143.144,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclobutane) + ring(Cycloheptene) + ring(Cyclobutane) + radical(RCCJCC) + radical(Allyl_S)"""),
)

species(
    label = 'C=C([O])C1C=CC2[CH]CC21(3524)',
    structure = SMILES('C=C([O])C1C=CC2[CH]CC21'),
    E0 = (307.599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.370148,0.0639798,6.62467e-07,-4.03932e-08,1.80553e-11,37140,29.3408], Tmin=(100,'K'), Tmax=(1039.36,'K')), NASAPolynomial(coeffs=[15.4979,0.0408366,-1.65599e-05,3.12415e-09,-2.22353e-13,32100.8,-53.3401], Tmin=(1039.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(307.599,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + polycyclic(s2_4_5_ene_1) + radical(cyclobutane) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C1OC2[CH]C1[CH]CC=C2(3360)',
    structure = SMILES('C=C1OC2[CH]C1[CH]CC=C2'),
    E0 = (254.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.854083,0.0835474,-1.45779e-05,-5.01687e-08,2.70348e-11,30798.7,4.98883], Tmin=(100,'K'), Tmax=(987.206,'K')), NASAPolynomial(coeffs=[26.4516,0.0267212,-9.99824e-06,1.95464e-09,-1.47805e-13,22785.2,-139.675], Tmin=(987.206,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(254.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Tetrahydrofuran) + ring(Cycloheptene) + ring(Tetrahydrofuran) + radical(CCJCO) + radical(Cs_S)"""),
)

species(
    label = 'C=C([O])C1C=C[CH]C2CC21(3525)',
    structure = SMILES('C=C([O])C1[CH]C=CC2CC21'),
    E0 = (240.056,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.587369,0.0492406,6.11674e-05,-1.15695e-07,4.77991e-11,29018.2,27.0518], Tmin=(100,'K'), Tmax=(968.569,'K')), NASAPolynomial(coeffs=[20.1493,0.0322275,-1.12499e-05,2.12928e-09,-1.59483e-13,22237.4,-82.1408], Tmin=(968.569,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(240.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + polycyclic(s2_3_6_ene_1) + radical(cyclohexene-allyl) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])C1C2[CH]C=CCC21(3526)',
    structure = SMILES('C=C([O])C1C2C=C[CH]CC21'),
    E0 = (240.056,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.587369,0.0492406,6.11674e-05,-1.15695e-07,4.77991e-11,29018.2,27.0518], Tmin=(100,'K'), Tmax=(968.569,'K')), NASAPolynomial(coeffs=[20.1493,0.0322275,-1.12499e-05,2.12928e-09,-1.59483e-13,22237.4,-82.1408], Tmin=(968.569,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(240.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + polycyclic(s2_3_6_ene_1) + radical(cyclohexene-allyl) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C1OC2[CH]C[CH]C1C=C2(3404)',
    structure = SMILES('C=C1OC2[CH]C[CH]C1C=C2'),
    E0 = (259.773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.76939,0.109231,-9.63843e-05,4.05251e-08,-6.63004e-12,31466,4.93504], Tmin=(100,'K'), Tmax=(1481.7,'K')), NASAPolynomial(coeffs=[29.4738,0.0248875,-1.0999e-05,2.1075e-09,-1.48058e-13,22207.3,-158.08], Tmin=(1481.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(259.773,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxane) + ring(oxepane) + ring(36dihydro2hpyran) + radical(CCJCO) + radical(Cs_S)"""),
)

species(
    label = 'C=C1OC2C=C[CH]C1[CH]C2(3310)',
    structure = SMILES('C=C1OC2[CH]C=CC1[CH]C2'),
    E0 = (177.971,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.72891,0.0781651,6.41453e-05,-1.80196e-07,8.50416e-11,21655.3,3.23895], Tmin=(100,'K'), Tmax=(935.703,'K')), NASAPolynomial(coeffs=[47.3434,-0.00956822,9.14277e-06,-1.61525e-09,8.57512e-14,7129.2,-258.795], Tmin=(935.703,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(177.971,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxane) + ring(Cycloheptane) + ring(Oxane) + radical(Cs_S) + radical(C=CCJCO)"""),
)

species(
    label = 'C=C1OC2[CH]C=CC[CH]C12(3375)',
    structure = SMILES('C=C1OC2[CH]C=CC[CH]C12'),
    E0 = (171.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.31767,0.0959308,-4.32168e-05,-2.21856e-08,1.72848e-11,20872.9,0.447323], Tmin=(100,'K'), Tmax=(1003.13,'K')), NASAPolynomial(coeffs=[27.0941,0.0285713,-1.11773e-05,2.16823e-09,-1.60754e-13,12861.7,-148.231], Tmin=(1003.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(171.801,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Oxetane) + ring(Cycloheptene) + ring(Oxetane) + radical(Cs_S) + radical(C=CCJCO)"""),
)

species(
    label = '[CH]=CC=CC(C=C)C(=C)[O](3527)',
    structure = SMILES('[CH]=CC=CC(C=C)C(=C)[O]'),
    E0 = (387.07,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,3120,650,792.5,1650,350,440,435,1725,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.627202,'amu*angstrom^2'), symmetry=1, barrier=(14.4206,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.627229,'amu*angstrom^2'), symmetry=1, barrier=(14.4212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.627022,'amu*angstrom^2'), symmetry=1, barrier=(14.4165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.626946,'amu*angstrom^2'), symmetry=1, barrier=(14.4147,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.15418,0.103631,-9.96427e-05,4.9657e-08,-9.7906e-12,46747.5,38.7428], Tmin=(100,'K'), Tmax=(1234.35,'K')), NASAPolynomial(coeffs=[21.2688,0.0309678,-1.13405e-05,1.96513e-09,-1.31247e-13,41212,-74.1561], Tmin=(1234.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(387.07,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=CC=CCC=CC(=C)[O](3528)',
    structure = SMILES('[CH]=CC=CCC=CC(=C)[O]'),
    E0 = (363.302,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.787332,'amu*angstrom^2'), symmetry=1, barrier=(18.1023,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.786551,'amu*angstrom^2'), symmetry=1, barrier=(18.0844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.786478,'amu*angstrom^2'), symmetry=1, barrier=(18.0827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.78689,'amu*angstrom^2'), symmetry=1, barrier=(18.0921,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.15307,0.0998955,-7.63033e-05,1.6492e-08,4.1078e-12,43892.7,37.7363], Tmin=(100,'K'), Tmax=(970.421,'K')), NASAPolynomial(coeffs=[23.7051,0.0269013,-9.02565e-06,1.566e-09,-1.08645e-13,37680.5,-88.5927], Tmin=(970.421,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(363.302,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = 'C=C(O)C1C=CC=CCC=1(3529)',
    structure = SMILES('C=C(O)C1C=CC=CCC=1'),
    E0 = (22.4829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.938062,0.087771,-2.7774e-05,-4.01783e-08,2.54432e-11,2900.8,26.985], Tmin=(100,'K'), Tmax=(952.459,'K')), NASAPolynomial(coeffs=[26.3041,0.0240293,-7.18183e-06,1.25876e-09,-9.26745e-14,-4586.78,-115.18], Tmin=(952.459,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(22.4829,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3,5-Cycloheptatriene)"""),
)

species(
    label = 'C=C(O)C1C=CC=CC=C1(3471)',
    structure = SMILES('C=C(O)C1C=CC=CC=C1'),
    E0 = (27.1507,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.28074,0.0971032,-5.21078e-05,-1.84154e-08,1.87653e-11,3472.97,28.1197], Tmin=(100,'K'), Tmax=(946.62,'K')), NASAPolynomial(coeffs=[27.6263,0.021472,-5.97396e-06,1.00555e-09,-7.33264e-14,-4083.99,-120.763], Tmin=(946.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(27.1507,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3,5-Cycloheptatriene)"""),
)

species(
    label = 'C=C([O])C1[CH]CCC=CC=1(3530)',
    structure = SMILES('C=C([O])C1[CH]C=CCCC=1'),
    E0 = (166.192,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.354442,0.0124675,0.000257393,-3.96227e-07,1.69594e-10,20183.2,31.964], Tmin=(100,'K'), Tmax=(901.444,'K')), NASAPolynomial(coeffs=[49.1487,-0.0259374,2.49216e-05,-5.11515e-09,3.38541e-13,4149.51,-238.517], Tmin=(901.444,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(166.192,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(C=C(C)OJ) + radical(C=CCJC=C)"""),
)

species(
    label = '[O]C12CC3C=CC1C2[CH]C3(3531)',
    structure = SMILES('[O]C12CC3C=CC1C2[CH]C3'),
    E0 = (393.342,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.877077,0.0593115,-1.20585e-05,-1.25475e-08,5.4027e-12,47428.1,23.5553], Tmin=(100,'K'), Tmax=(1252.57,'K')), NASAPolynomial(coeffs=[10.8055,0.046641,-1.96808e-05,3.64217e-09,-2.50208e-13,43447.6,-32.5403], Tmin=(1252.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(393.342,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s3_6_6_ene_1) + polycyclic(s2_3_6_ane) + polycyclic(s2_3_6_ene_1) - ring(Cyclohexane) - ring(Cyclopropane) - ring(Cyclohexene) + radical(bicyclo[4.1.0]heptane-C6-2) + radical(CC(C)2OJ)"""),
)

species(
    label = 'C=C([O])[CH]C1C=CC=CC1(3461)',
    structure = SMILES('[CH2]C([O])=CC1C=CC=CC1'),
    E0 = (206.212,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.227949,0.0766482,-2.09872e-05,-3.08845e-08,1.83348e-11,24968.3,30.8308], Tmin=(100,'K'), Tmax=(978.032,'K')), NASAPolynomial(coeffs=[19.6452,0.0329817,-1.17009e-05,2.1058e-09,-1.48907e-13,19282.2,-73.8005], Tmin=(978.032,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(206.212,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(C=C(C)OJ) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]C1C=CC=CC1C(=C)[O](3473)',
    structure = SMILES('[CH2]C1C=CC=CC1C(=C)[O]'),
    E0 = (256.651,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,350,440,435,1725,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.154912,0.0741253,-1.17509e-05,-4.28194e-08,2.34639e-11,31033.2,31.1609], Tmin=(100,'K'), Tmax=(956.119,'K')), NASAPolynomial(coeffs=[19.9857,0.03155,-1.03535e-05,1.8047e-09,-1.26992e-13,25276.5,-75.066], Tmin=(956.119,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(256.651,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3-Cyclohexadiene) + radical(C=C(C)OJ) + radical(Isobutyl)"""),
)

species(
    label = 'C=C1OC2CC=CC=CC12(3347)',
    structure = SMILES('C=C1OC2CC=CC=CC12'),
    E0 = (-44.3644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.23955,0.167464,-0.000321901,3.2722e-07,-1.24392e-10,-5140.36,-17.6623], Tmin=(100,'K'), Tmax=(842.373,'K')), NASAPolynomial(coeffs=[1.19937,0.08122,-4.38342e-05,8.62649e-09,-5.99106e-13,-3239.21,-18.9398], Tmin=(842.373,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-44.3644,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Oxetane) + ring(1,3-Cycloheptadiene) + ring(Oxetane)"""),
)

species(
    label = 'C=C([O])C1[CH]CC2C=CC21(3532)',
    structure = SMILES('C=C([O])C1[CH]CC2C=CC21'),
    E0 = (343.765,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.696157,0.0502744,4.73807e-05,-9.44977e-08,3.85568e-11,41484.2,29.9248], Tmin=(100,'K'), Tmax=(987.042,'K')), NASAPolynomial(coeffs=[17.6568,0.0359857,-1.36438e-05,2.60284e-09,-1.91302e-13,35483.9,-65.1139], Tmin=(987.042,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + Estimated bicyclic component: polycyclic(s2_4_5_ane) - ring(Cyclopentane) - ring(Cyclobutane) + ring(Cyclopentane) + ring(Cyclobutene) + radical(Cs_S) + radical(C=C(C)OJ)"""),
)

species(
    label = 'O(20)',
    structure = SMILES('[O]'),
    E0 = (243.005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (15.9994,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,29226.7,5.11107], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,29226.7,5.11107], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.005,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'C=[C]C1[CH]CC=CC=C1(3533)',
    structure = SMILES('C=[C]C1[CH]CC=CC=C1'),
    E0 = (576.936,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2950,3100,1380,975,1025,1650,1685,370,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (118.176,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0642286,0.086085,-6.89639e-05,2.88751e-08,-4.97435e-12,69538.2,16.1368], Tmin=(100,'K'), Tmax=(1353.77,'K')), NASAPolynomial(coeffs=[15.9241,0.0388438,-1.66195e-05,3.09783e-09,-2.14053e-13,65209.3,-65.8404], Tmin=(1353.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(576.936,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3-Cycloheptadiene) + radical(Cds_S) + radical(Cs_S)"""),
)

species(
    label = 'CC(=O)[C]1[CH]CC=CC=C1(3534)',
    structure = SMILES('CC([O])=C1[CH]CC=CC=C1'),
    E0 = (182.846,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0693066,0.0201347,0.000237042,-3.76519e-07,1.62762e-10,22195.2,30.0022], Tmin=(100,'K'), Tmax=(902.08,'K')), NASAPolynomial(coeffs=[50.079,-0.0269482,2.48875e-05,-5.08098e-09,3.35804e-13,6065.82,-245.503], Tmin=(902.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(182.846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(Cycloheptane) + radical(C=C(C)OJ) + radical(Allyl_S)"""),
)

species(
    label = 'CC(=O)C1[C]=CC=CC[CH]1(3535)',
    structure = SMILES('CC(=O)C1[C]=CC=CC[CH]1'),
    E0 = (342.249,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,375,552.5,462.5,1710,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.484367,0.0950961,-8.10422e-05,3.62578e-08,-6.63447e-12,41327.2,20.3161], Tmin=(100,'K'), Tmax=(1288.96,'K')), NASAPolynomial(coeffs=[17.2195,0.0401567,-1.71083e-05,3.19072e-09,-2.21018e-13,36763.2,-69.5889], Tmin=(1288.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(342.249,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(CCJCC=O) + radical(Cds_S)"""),
)

species(
    label = 'CC(=O)C1[CH]C=C[CH]C=C1(3439)',
    structure = SMILES('CC(=O)C1C=C[CH][CH]C=C1'),
    E0 = (201.824,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.806355,0.0153196,0.000212392,-3.1926e-07,1.33402e-10,24440,29.9546], Tmin=(100,'K'), Tmax=(914.384,'K')), NASAPolynomial(coeffs=[38.0359,-0.00495435,1.17416e-05,-2.42905e-09,1.50543e-13,11670.7,-178.919], Tmin=(914.384,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(201.824,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(Allyl_S) + radical(Allyl_S)"""),
)

species(
    label = 'CC(=O)C1[CH]CC=C[C]=C1(3536)',
    structure = SMILES('CC(=O)C1[CH]CC=C[C]=C1'),
    E0 = (303.402,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,375,552.5,462.5,1710,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.542428,0.0971969,-8.59642e-05,4.08495e-08,-7.9663e-12,36656.6,19.6479], Tmin=(100,'K'), Tmax=(1215.34,'K')), NASAPolynomial(coeffs=[16.3583,0.0415725,-1.73116e-05,3.19068e-09,-2.19765e-13,32548.6,-65.1842], Tmin=(1215.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.402,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(C=CJC=C) + radical(CCJCC=O)"""),
)

species(
    label = 'CC(=O)C1[CH]C[C]=CC=C1(3537)',
    structure = SMILES('CC(=O)C1[CH]C[C]=CC=C1'),
    E0 = (342.249,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,375,552.5,462.5,1710,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.484367,0.0950961,-8.10422e-05,3.62578e-08,-6.63447e-12,41327.2,20.3161], Tmin=(100,'K'), Tmax=(1288.96,'K')), NASAPolynomial(coeffs=[17.2195,0.0401567,-1.71083e-05,3.19072e-09,-2.21018e-13,36763.2,-69.5889], Tmin=(1288.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(342.249,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(CCJCC=O) + radical(Cds_S)"""),
)

species(
    label = 'CC(=O)C1[CH]CC=[C]C=C1(3538)',
    structure = SMILES('CC(=O)C1[CH]CC=[C]C=C1'),
    E0 = (303.402,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,375,552.5,462.5,1710,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.542428,0.0971969,-8.59642e-05,4.08495e-08,-7.9663e-12,36656.6,19.6479], Tmin=(100,'K'), Tmax=(1215.34,'K')), NASAPolynomial(coeffs=[16.3583,0.0415725,-1.73116e-05,3.19068e-09,-2.19765e-13,32548.6,-65.1842], Tmin=(1215.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.402,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(C=CJC=C) + radical(CCJCC=O)"""),
)

species(
    label = 'C=C1OC2CC=C[CH][CH]C12(3331)',
    structure = SMILES('C=C1OC2C[CH]C=C[CH]C12'),
    E0 = (142.567,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.871068,0.0824987,-5.42559e-06,-6.19525e-08,3.17462e-11,17344.4,-3.12538], Tmin=(100,'K'), Tmax=(977.858,'K')), NASAPolynomial(coeffs=[26.9131,0.0270758,-9.73235e-06,1.88103e-09,-1.4257e-13,9126.53,-150.781], Tmin=(977.858,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(142.567,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Oxetane) + ring(Cycloheptene) + ring(Oxetane) + radical(Allyl_S) + radical(Allyl_S)"""),
)

species(
    label = 'O=C1CC2[CH]C1[CH]CC=C2(3073)',
    structure = SMILES('O=C1CC2[CH]C1[CH]CC=C2'),
    E0 = (209.582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.684393,0.0936118,-7.25684e-05,2.8467e-08,-4.50742e-12,25383,5.27979], Tmin=(100,'K'), Tmax=(1485.35,'K')), NASAPolynomial(coeffs=[20.3131,0.0370662,-1.5465e-05,2.83735e-09,-1.93677e-13,19145.3,-104.328], Tmin=(1485.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(209.582,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclopentane) + ring(Cycloheptene) + ring(Cyclopentane) + radical(CCJCC=O) + radical(CCJCC=O)"""),
)

species(
    label = 'O=C1CC2[CH]C[CH]C1C=C2(3061)',
    structure = SMILES('O=C1CC2[CH]C[CH]C1C=C2'),
    E0 = (207.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.41434,0.101035,-0.000103767,6.02161e-08,-1.45469e-11,25053.6,4.58917], Tmin=(100,'K'), Tmax=(987.116,'K')), NASAPolynomial(coeffs=[13.1422,0.0461008,-2.02893e-05,3.83708e-09,-2.68029e-13,22377.3,-60.6375], Tmin=(987.116,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(207.013,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclohexane) + ring(Cycloheptene) + ring(Cyclohexene) + radical(CCJCC=O) + radical(Cs_S)"""),
)

species(
    label = 'O=C1CC2C=C[CH]C1[CH]C2(3111)',
    structure = SMILES('O=C1CC2[CH]C=CC1[CH]C2'),
    E0 = (148.588,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.287527,0.0924059,-7.07701e-05,2.78601e-08,-4.53602e-12,18026.3,-0.464518], Tmin=(100,'K'), Tmax=(1411.22,'K')), NASAPolynomial(coeffs=[16.6829,0.0443043,-1.96422e-05,3.70705e-09,-2.57253e-13,13236.5,-88.1827], Tmin=(1411.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(148.588,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclohexane) + ring(Cycloheptene) + ring(Cyclohexane) + radical(Allyl_S) + radical(CCJCC=O)"""),
)

species(
    label = 'O=C1CC2C=C[CH]C[CH]C12(3091)',
    structure = SMILES('O=C1CC2[CH]C=CC[CH]C12'),
    E0 = (149.723,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.580062,0.0927089,-6.96612e-05,2.63438e-08,-4.03895e-12,18178.6,1.26619], Tmin=(100,'K'), Tmax=(1519.71,'K')), NASAPolynomial(coeffs=[19.8458,0.0389469,-1.65971e-05,3.06589e-09,-2.09664e-13,11970.2,-105.826], Tmin=(1519.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(149.723,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclobutane) + ring(Cycloheptene) + ring(Cyclobutane) + radical(CCJCC=O) + radical(Allyl_S)"""),
)

species(
    label = 'CC(=O)C1C=CC=CCC=1(3539)',
    structure = SMILES('CC(=O)C1C=CC=CCC=1'),
    E0 = (6.39494,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.761656,0.0864355,-3.60545e-05,-2.10874e-08,1.54745e-11,956.7,28.2842], Tmin=(100,'K'), Tmax=(1004.18,'K')), NASAPolynomial(coeffs=[23.3479,0.0302397,-1.16244e-05,2.20353e-09,-1.60286e-13,-5894.13,-98.1332], Tmin=(1004.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.39494,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)HHH) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene)"""),
)

species(
    label = 'CC(=O)C1C=CC=CC=C1(3475)',
    structure = SMILES('CC(=O)C1C=CC=CC=C1'),
    E0 = (-5.24664,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.11893,0.0988779,-8.40889e-05,3.63453e-08,-6.23677e-12,-435.304,29.0429], Tmin=(100,'K'), Tmax=(1403.64,'K')), NASAPolynomial(coeffs=[22.5326,0.0314772,-1.20609e-05,2.13513e-09,-1.43629e-13,-7074.93,-93.0815], Tmin=(1403.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-5.24664,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene)"""),
)

species(
    label = '[CH2]C12OC3C=CC1C2[CH]C3(3540)',
    structure = SMILES('[CH2]C12OC3C=CC1C2[CH]C3'),
    E0 = (412.848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23212,0.0631928,-1.46953e-05,-1.15817e-08,4.19426e-12,49751.3,22.2582], Tmin=(100,'K'), Tmax=(1597.53,'K')), NASAPolynomial(coeffs=[22.6161,0.0415496,-2.43254e-05,4.93636e-09,-3.46741e-13,38848.5,-103.665], Tmin=(1597.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(412.848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s3_6_6_ene_1) + polycyclic(s2_3_6_ane) + polycyclic(s2_3_6_ene_1) - ring(Oxane) - ring(Cyclopropane) - ring(36dihydro2hpyran) + radical(CJC(C)OC) + radical(Cs_S)"""),
)

species(
    label = 'O=[C]CC1[CH]CC=CC=C1(2244)',
    structure = SMILES('O=[C]CC1[CH]CC=CC=C1'),
    E0 = (289.021,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4440.56,'J/mol'), sigma=(7.29278,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=693.61 K, Pc=25.98 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.570022,0.0938456,-7.73464e-05,3.29559e-08,-5.69593e-12,34931.1,22.2447], Tmin=(100,'K'), Tmax=(1367.25,'K')), NASAPolynomial(coeffs=[18.6404,0.037644,-1.56881e-05,2.89157e-09,-1.9871e-13,29678,-76.4436], Tmin=(1367.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(289.021,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3-Cycloheptadiene) + radical(CCCJ=O) + radical(Cs_S)"""),
)

species(
    label = 'O=C1CC2CC=CC=CC12(2251)',
    structure = SMILES('O=C1CC2CC=CC=CC12'),
    E0 = (-96.5854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.84629,0.170942,-0.00036477,3.95716e-07,-1.54816e-10,-11446.8,-15.5347], Tmin=(100,'K'), Tmax=(856.815,'K')), NASAPolynomial(coeffs=[-8.63575,0.0977841,-5.31287e-05,1.0408e-08,-7.17812e-13,-6434.52,38.6316], Tmin=(856.815,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-96.5854,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclobutane) + ring(1,3-Cycloheptadiene) + ring(Cyclobutane)"""),
)

species(
    label = '[CH2]C1=C[CH]CC=CC=CO1(3541)',
    structure = SMILES('C=C1C=CC[CH]C=C[CH]O1'),
    E0 = (251.425,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.861262,0.090351,-4.77407e-05,-6.62152e-09,9.68712e-12,30428.8,17.7242], Tmin=(100,'K'), Tmax=(1028.95,'K')), NASAPolynomial(coeffs=[22.3932,0.0332977,-1.31823e-05,2.48476e-09,-1.78039e-13,23877.9,-103.707], Tmin=(1028.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(251.425,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclononane) + radical(C=CCJ(O)C) + radical(Allyl_S)"""),
)

species(
    label = 'CH2(T)(33)',
    structure = SMILES('[CH2]'),
    E0 = (381.08,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([971.045,2816.03,3444.23],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.71758,0.00127391,2.17347e-06,-3.48858e-09,1.65209e-12,45872.4,1.75298], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.14632,0.00303671,-9.96474e-07,1.50484e-10,-8.57336e-15,46041.3,4.72342], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(381.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'O=[C]C1[CH]CC=CC=C1(3542)',
    structure = SMILES('O=[C]C1[CH]CC=CC=C1'),
    E0 = (317.966,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,1855,455,950,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0365209,0.0773619,-5.37797e-05,1.24525e-08,1.1532e-12,38397.6,17.8025], Tmin=(100,'K'), Tmax=(1076.62,'K')), NASAPolynomial(coeffs=[18.1407,0.0283326,-1.12518e-05,2.08318e-09,-1.46128e-13,33411.1,-76.2148], Tmin=(1076.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(317.966,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3-Cycloheptadiene) + radical(CCJCC=O) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'N2',
    structure = SMILES('N#N'),
    E0 = (-8.69489,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0135,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(810.913,'J/mol'), sigma=(3.621,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.61263,-0.00100893,2.49898e-06,-1.43376e-09,2.58636e-13,-1051.1,2.6527], Tmin=(100,'K'), Tmax=(1817.04,'K')), NASAPolynomial(coeffs=[2.9759,0.00164141,-7.19722e-07,1.25378e-10,-7.91526e-15,-1025.84,5.53757], Tmin=(1817.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.69489,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""N2""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'Ne',
    structure = SMILES('[Ne]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (20.1797,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1235.53,'J/mol'), sigma=(3.758e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,3.35532], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,3.35532], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ne""", comment="""Thermo library: primaryThermoLibrary"""),
)

transitionState(
    label = 'TS1',
    E0 = (262.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (402.203,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (385.845,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (397.905,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (342.422,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (315.706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (390.26,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (386.598,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (447.059,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (447.059,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (419.451,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (563.761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (406.686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (311.102,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (311.102,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (359.203,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (323.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (356.571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (514.037,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (443.084,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (493.371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (493.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (393.198,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (387.582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (393.198,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (343.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (263.596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (263.596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (358.691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (289.537,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (311.832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (428.115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (432,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (325.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (287.314,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (308.293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (667.813,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (422.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (507.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (270.625,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (384.754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (819.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (379.399,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (386.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (321.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (336.442,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (484.818,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (445.972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (387.688,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (343.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (358.691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (289.537,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (329.196,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (325.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (287.314,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (687.318,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (534.622,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (270.625,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (565.225,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (699.046,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    products = ['C7H8(11)', 'CH2CO(43)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    products = ['[CH2]C1([O])C2C=CC=CCC21(3450)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.73316e+11,'s^-1'), n=0.231216, Ea=(139.863,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_csHNd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(19)', 'C=C([O])C1C=CC=CCC=1(3515)'],
    products = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.45e+08,'cm^3/(mol*s)'), n=1.64, Ea=(13.7654,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2607 used for Cds-CdCd_Cds-CsH;HJ
Exact match found for rate rule [Cds-CdCd_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(19)', 'C=C([O])C1C=CC=CC=C1(3463)'],
    products = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.0545849,'m^3/(mol*s)'), n=2.81111, Ea=(21.1569,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 26 used for Cds-CdH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CdH_Cds-CsH;HJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C7H8(11)', '[CH2][C]=O(360)'],
    products = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.0147295,'m^3/(mol*s)'), n=2.37207, Ea=(16.7338,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CdH_Cds-CsH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]1C=C[CH]CC=C1(125)', 'CH2CO(43)'],
    products = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7051.36,'m^3/(mol*s)'), n=1.07323, Ea=(43.3386,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cdd_Od;CJ] + [Ck_O;YJ] for rate rule [Ck_O;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    products = ['C=C([O])[C]1C=CC=CCC1(3465)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.15233e+09,'s^-1'), n=1.014, Ea=(127.919,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    products = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.022e+10,'s^-1'), n=1.34, Ea=(199.577,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 9 used for R2H_S;C_rad_out_H/(Cd-Cd-Cd);Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/(Cd-Cd-Cd);Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=C([O])C1[C]=CC=CCC1(3467)'],
    products = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.32e+09,'s^-1'), n=0.99, Ea=(141.419,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 206 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=C([O])C1C=CC=[C]CC1(3464)'],
    products = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.32e+09,'s^-1'), n=0.99, Ea=(141.419,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 206 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    products = ['C=C(O)[C]1[CH]CC=CC=C1(3516)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.19923e+08,'s^-1'), n=1.46351, Ea=(157.11,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;XH_out] + [R3H_SS_2Cd;Y_rad_out;XH_out] for rate rule [R3H_SS_2Cd;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=C(O)C1[CH]CC=CC=C1(3517)'],
    products = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C(O)C1[C]=CC=CC[CH]1(3518)'],
    products = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;XH_out] for rate rule [R4H_SS(Cd)S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=C([O])C1C=[C]C=CCC1(3468)'],
    products = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=C([O])C1C=C[C]=CCC1(3466)'],
    products = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=C([O])C1C=CC=CCC1(3469)'],
    products = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    products = ['C=C(O)C1[CH]C=C[CH]C=C1(3422)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(61085.3,'s^-1'), n=1.8624, Ea=(60.8997,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R5HJ_3;O_rad_out;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=C(O)C1[CH]CC=C[C]=C1(3519)'],
    products = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_single;XH_out] for rate rule [R5H_DSSS;Cd_rad_out_singleDe_Cd;O_H_out]
Euclidian distance = 3.74165738677
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=C(O)C1[CH]C[C]=CC=C1(3520)'],
    products = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.81182e+08,'s^-1'), n=1.25566, Ea=(151.659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;XH_out] for rate rule [R6HJ_2;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C(O)C1[CH]CC=[C]C=C1(3521)'],
    products = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.62236e+06,'s^-1'), n=1.59783, Ea=(119.552,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_Cd;XH_out] + [R6H_SMSSR;Y_rad_out;XH_out] for rate rule [R6H_SMSSR;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 3.16227766017
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]1C=C[CH]CC=C1(125)', '[CH2][C]=O(360)'],
    products = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.47003e+07,'m^3/(mol*s)'), n=0.0284742, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    products = ['[CH]1CC=CC=CC1[C]1CO1(3522)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    products = ['C=C([O])C1[CH]C2C=CCC21(3523)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(130.857,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 116 used for Rn0c7_gamma_short_SSDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCs
Exact match found for rate rule [Rn0c7_gamma_short_SSDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    products = ['O=C1CC2CC=C[CH][CH]C12(3125)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.82767e+08,'s^-1'), n=1.02667, Ea=(125.241,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    products = ['C=C([O])C1C=CC2[CH]CC21(3524)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(130.857,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 116 used for Rn0c7_gamma_short_SSDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCs
Exact match found for rate rule [Rn0c7_gamma_short_SSDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    products = ['C=C1OC2[CH]C1[CH]CC=C2(3360)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(9.9e+10,'s^-1'), n=0.06, Ea=(80.7512,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn2cx_beta;doublebond_intra_pri_HCd;radadd_intra] for rate rule [Rn2c7_beta_short;doublebond_intra_pri_HCd;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    products = ['C=C([O])C1C=C[CH]C2CC21(3525)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.1e+12,'s^-1'), n=0.14, Ea=(1.2552,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    products = ['C=C([O])C1C2[CH]C=CCC21(3526)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.1e+12,'s^-1'), n=0.14, Ea=(1.2552,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    products = ['C=C1OC2[CH]C[CH]C1C=C2(3404)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4.54269e+17,'s^-1'), n=-0.500262, Ea=(96.3509,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_cyclic;doublebond_intra_pri_HCd;radadd_intra] for rate rule [Rn2c7_gamma_long;doublebond_intra_pri_HCd;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    products = ['C=C1OC2C=C[CH]C1[CH]C2(3310)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(110000,'s^-1'), n=1.18, Ea=(27.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn2c7_gamma_long;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    products = ['C=C1OC2[CH]C=CC[CH]C12(3375)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.99832e+10,'s^-1'), n=0.37247, Ea=(49.4918,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri;radadd_intra] for rate rule [Rn2c7_alpha_long;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 3.16227766017
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]=CC=CC(C=C)C(=C)[O](3527)'],
    products = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(4.625e+11,'s^-1'), n=0.16, Ea=(41.045,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;doublebond_intra_pri;radadd_intra_cdsingleH] for rate rule [R7_linear;doublebond_intra_pri_2H;radadd_intra_cdsingleH]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]=CC=CCC=CC(=C)[O](3528)'],
    products = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(4.58366e+14,'s^-1'), n=-0.170131, Ea=(68.698,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7;doublebond_intra_pri;radadd_intra_cdsingleH] + [R7;doublebond_intra_pri_HCd;radadd_intra] for rate rule [R7_linear;doublebond_intra_pri_HCd;radadd_intra_cdsingleH]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    products = ['C=C(O)C1C=CC=CCC=1(3529)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    products = ['C=C(O)C1C=CC=CC=C1(3471)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    products = ['C=C([O])C1[CH]CCC=CC=1(3530)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2.125e+09,'s^-1'), n=0.991, Ea=(45.9529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_3_pentadiene;CH_end;CdHC_2] for rate rule [1_3_pentadiene;CH(CJ)_1;CdHC_2]
Euclidian distance = 1.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[O]C12CC3C=CC1C2[CH]C3(3531)'],
    products = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.24e+10,'s^-1'), n=1.27, Ea=(274.47,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [cyclohexene] for rate rule [cyclohexene_1inring]
Euclidian distance = 1.0
family: Intra_Retro_Diels_alder_bicyclic"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    products = ['C=C([O])[CH]C1C=CC=CC1(3461)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    products = ['[CH2]C1C=CC=CC1C(=C)[O](3473)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HH)CJ;CsJ;C] for rate rule [cCs(-HH)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    products = ['C=C1OC2CC=CC=CC12(3347)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    products = ['C=C([O])C1[CH]CC2C=CC21(3532)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1,3-butadiene_backbone;C=C_1;C=C_2] for rate rule [1,3-butadiene_backbone;CdH(C)_1;CdH(C)_2]
Euclidian distance = 1.41421356237
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction42',
    reactants = ['O(20)', 'C=[C]C1[CH]CC=CC=C1(3533)'],
    products = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction43',
    reactants = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    products = ['CC(=O)[C]1[CH]CC=CC=C1(3534)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(2.30996,'s^-1'), n=3.5644, Ea=(117.059,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['CC(=O)C1[C]=CC=CC[CH]1(3535)'],
    products = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    products = ['CC(=O)C1[CH]C=C[CH]C=C1(3439)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(62296.1,'s^-1'), n=1.86, Ea=(59.4128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;Cs_H_out_H/Cd] for rate rule [R5HJ_3;C_rad_out_2H;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['CC(=O)C1[CH]CC=C[C]=C1(3536)'],
    products = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_single;Cs_H_out_2H] for rate rule [R5H_DSSS;Cd_rad_out_singleDe_Cd;Cs_H_out_2H]
Euclidian distance = 3.60555127546
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['CC(=O)C1[CH]C[C]=CC=C1(3537)'],
    products = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(5.59786e+07,'s^-1'), n=1.58088, Ea=(142.57,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_2H] for rate rule [R6HJ_2;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['CC(=O)C1[CH]CC=[C]C=C1(3538)'],
    products = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(5.59786e+07,'s^-1'), n=1.58088, Ea=(142.57,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_2H] for rate rule [R6H_SMSSR;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    products = ['C=C1OC2CC=C[CH][CH]C12(3331)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(4.7342e+08,'s^-1'), n=0.920995, Ea=(125.348,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_csHCs] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_csHCs]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction50',
    reactants = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    products = ['O=C1CC2[CH]C1[CH]CC=C2(3073)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(9.9e+10,'s^-1'), n=0.06, Ea=(80.7512,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn2cx_beta;doublebond_intra_pri_HCd;radadd_intra_cs2H] for rate rule [Rn2c7_beta_short;doublebond_intra_pri_HCd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction51',
    reactants = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    products = ['O=C1CC2[CH]C[CH]C1C=C2(3061)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(4.54269e+17,'s^-1'), n=-0.500262, Ea=(96.3509,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_cyclic;doublebond_intra_pri_HCd;radadd_intra_cs] for rate rule [Rn2c7_gamma_long;doublebond_intra_pri_HCd;radadd_intra_cs2H]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction52',
    reactants = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    products = ['O=C1CC2C=C[CH]C1[CH]C2(3111)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(110000,'s^-1'), n=1.18, Ea=(27.196,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R7;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H] for rate rule [Rn2c7_gamma_long;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction53',
    reactants = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    products = ['O=C1CC2C=C[CH]C[CH]C12(3091)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(2.54811e+07,'s^-1'), n=0.785425, Ea=(66.8554,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri;radadd_intra_cs2H] for rate rule [Rn2c7_alpha_long;doublebond_intra_pri;radadd_intra_cs2H]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction54',
    reactants = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    products = ['CC(=O)C1C=CC=CCC=1(3539)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction55',
    reactants = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    products = ['CC(=O)C1C=CC=CC=C1(3475)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction56',
    reactants = ['[CH2]C12OC3C=CC1C2[CH]C3(3540)'],
    products = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(1.24e+10,'s^-1'), n=1.27, Ea=(274.47,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [cyclohexene] for rate rule [cyclohexene_1inring]
Euclidian distance = 1.0
family: Intra_Retro_Diels_alder_bicyclic"""),
)

reaction(
    label = 'reaction38',
    reactants = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    products = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction58',
    reactants = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    products = ['O=C1CC2CC=CC=CC12(2251)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction59',
    reactants = ['[CH2]C1=C[CH]CC=CC=CO1(3541)'],
    products = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction60',
    reactants = ['CH2(T)(33)', 'O=[C]C1[CH]CC=CC=C1(3542)'],
    products = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(2.23625e+06,'m^3/(mol*s)'), n=0.36814, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

network(
    label = '491',
    isomers = [
        'C=C([O])C1[CH]CC=CC=C1(2238)',
    ],
    reactants = [
        ('C7H8(11)', 'CH2CO(43)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '491',
    Tmin = (300,'K'),
    Tmax = (2000,'K'),
    Tcount = 8,
    Tlist = ([302.47,323.145,369.86,455.987,609.649,885.262,1353.64,1896.74],'K'),
    Pmin = (0.01,'bar'),
    Pmax = (100,'bar'),
    Pcount = 5,
    Plist = ([0.0125282,0.0667467,1,14.982,79.8202],'bar'),
    maximumGrainSize = (0.5,'kcal/mol'),
    minimumGrainCount = 250,
    method = 'modified strong collision',
    interpolationModel = ('Chebyshev', 6, 4),
    activeKRotor = True,
    activeJRotor = True,
    rmgmode = True,
)

