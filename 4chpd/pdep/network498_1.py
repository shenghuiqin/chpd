species(
    label = 'O=[C]CC1[CH]C=CC=CC1(2245)',
    structure = SMILES('O=[C]CC1C=C[CH]C=CC1'),
    E0 = (213.702,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.819605,0.00725946,0.000251568,-3.74444e-07,1.57485e-10,25875.7,34.5197], Tmin=(100,'K'), Tmax=(906.275,'K')), NASAPolynomial(coeffs=[43.2142,-0.0160425,1.90032e-05,-3.91762e-09,2.54642e-13,11464.2,-202.952], Tmin=(906.275,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(213.702,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(1,4-Cycloheptadiene) + radical(CCCJ=O) + radical(C=CCJC=C)"""),
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
    label = 'O=C1CC2[CH]C1[CH]C=CC2(3112)',
    structure = SMILES('O=C1CC2[CH]C1C=C[CH]C2'),
    E0 = (148.588,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.287527,0.0924059,-7.07701e-05,2.78601e-08,-4.53602e-12,18026.3,-0.464518], Tmin=(100,'K'), Tmax=(1411.22,'K')), NASAPolynomial(coeffs=[16.6829,0.0443043,-1.96422e-05,3.70705e-09,-2.57253e-13,13236.5,-88.1827], Tmin=(1411.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(148.588,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclopentane) + ring(Cycloheptene) + ring(Cyclopentane) + radical(Allyl_S) + radical(CCJCC=O)"""),
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
    label = 'O=C=CC1C=C[CH]C=CC1(3113)',
    structure = SMILES('O=C=CC1C=C[CH]C=CC1'),
    E0 = (176.784,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2120,512.5,787.5,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (133.167,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.73319,0.00869962,0.000246418,-3.72057e-07,1.57574e-10,21439.1,31.1762], Tmin=(100,'K'), Tmax=(905.209,'K')), NASAPolynomial(coeffs=[44.7845,-0.0206239,2.1039e-05,-4.29813e-09,2.80744e-13,6690.26,-214.374], Tmin=(905.209,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(176.784,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCd(CCO)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + ring(1,4-Cycloheptadiene) + radical(C=CCJC=C)"""),
)

species(
    label = 'O=[C]CC1C=CC=C=CC1(3114)',
    structure = SMILES('O=[C]CC1C=CC=C=CC1'),
    E0 = (317.379,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (133.167,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.141159,0.0707136,-2.14547e-05,-2.00404e-08,1.18138e-11,38323,28.1297], Tmin=(100,'K'), Tmax=(1039.58,'K')), NASAPolynomial(coeffs=[16.6378,0.0365967,-1.45874e-05,2.7203e-09,-1.92287e-13,33306.7,-59.7274], Tmin=(1039.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(317.379,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C]CC1C=C=CC=CC1(3115)',
    structure = SMILES('O=[C]CC1C=C=CC=CC1'),
    E0 = (317.379,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (133.167,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.141159,0.0707136,-2.14547e-05,-2.00404e-08,1.18138e-11,38323,28.1297], Tmin=(100,'K'), Tmax=(1039.58,'K')), NASAPolynomial(coeffs=[16.6378,0.0365967,-1.45874e-05,2.7203e-09,-1.92287e-13,33306.7,-59.7274], Tmin=(1039.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(317.379,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(CCCJ=O)"""),
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
    label = 'O=C[CH]C1C=C[CH]C=CC1(3116)',
    structure = SMILES('[O]C=CC1C=C[CH]C=CC1'),
    E0 = (196.445,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.561624,0.00563692,0.000277166,-4.14845e-07,1.7527e-10,23816.4,31.9534], Tmin=(100,'K'), Tmax=(905.223,'K')), NASAPolynomial(coeffs=[48.881,-0.0247254,2.39867e-05,-4.87564e-09,3.18379e-13,7564.4,-237.797], Tmin=(905.223,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(196.445,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(1,4-Cycloheptadiene) + radical(C=CCJC=C) + radical(C=COJ)"""),
)

species(
    label = 'O=[C]CC1C=CC[C]=CC1(3101)',
    structure = SMILES('O=[C]CC1C=CC[C]=CC1'),
    E0 = (349.819,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.478119,0.0241397,0.000189453,-3.00692e-07,1.28701e-10,42250.2,36.5479], Tmin=(100,'K'), Tmax=(908.201,'K')), NASAPolynomial(coeffs=[39.7522,-0.00932495,1.43066e-05,-2.98617e-09,1.92773e-13,29362.9,-180.821], Tmin=(908.201,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(349.819,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(1,4-Cycloheptadiene) + radical(CCCJ=O) + radical(Cds_S)"""),
)

species(
    label = 'O=[C]CC1C=[C]CC=CC1(3102)',
    structure = SMILES('O=[C]CC1C=[C]CC=CC1'),
    E0 = (349.819,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2850,1437.5,1250,1305,750,350,1855,455,950,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.478119,0.0241397,0.000189453,-3.00692e-07,1.28701e-10,42250.2,36.5479], Tmin=(100,'K'), Tmax=(908.201,'K')), NASAPolynomial(coeffs=[39.7522,-0.00932495,1.43066e-05,-2.98617e-09,1.92773e-13,29362.9,-180.821], Tmin=(908.201,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(349.819,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(1,4-Cycloheptadiene) + radical(CCCJ=O) + radical(Cds_S)"""),
)

species(
    label = 'O=CC[C]1C=C[CH]C=CC1(3117)',
    structure = SMILES('O=CC[C]1C=C[CH]C=CC1'),
    E0 = (185.722,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1556,-0.00092565,0.000271432,-3.9059e-07,1.61962e-10,22499.1,31.4411], Tmin=(100,'K'), Tmax=(908.072,'K')), NASAPolynomial(coeffs=[41.2828,-0.0123916,1.73347e-05,-3.59073e-09,2.3109e-13,8396.48,-195.804], Tmin=(908.072,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(185.722,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(1,4-Cycloheptadiene) + radical(C=CCJC=C) + radical(Allyl_T)"""),
)

species(
    label = 'O=[C]CC1C=CCC=[C]C1(3098)',
    structure = SMILES('O=[C]CC1C=CCC=[C]C1'),
    E0 = (349.819,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2850,1437.5,1250,1305,750,350,1855,455,950,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.478119,0.0241397,0.000189453,-3.00692e-07,1.28701e-10,42250.2,36.5479], Tmin=(100,'K'), Tmax=(908.201,'K')), NASAPolynomial(coeffs=[39.7522,-0.00932495,1.43066e-05,-2.98617e-09,1.92773e-13,29362.9,-180.821], Tmin=(908.201,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(349.819,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(1,4-Cycloheptadiene) + radical(CCCJ=O) + radical(Cds_S)"""),
)

species(
    label = 'O=[C]CC1[C]=CCC=CC1(3100)',
    structure = SMILES('O=[C]CC1[C]=CCC=CC1'),
    E0 = (349.819,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2850,1437.5,1250,1305,750,350,1855,455,950,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.478119,0.0241397,0.000189453,-3.00692e-07,1.28701e-10,42250.2,36.5479], Tmin=(100,'K'), Tmax=(908.201,'K')), NASAPolynomial(coeffs=[39.7522,-0.00932495,1.43066e-05,-2.98617e-09,1.92773e-13,29362.9,-180.821], Tmin=(908.201,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(349.819,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(1,4-Cycloheptadiene) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = 'O=CCC1[CH][CH]C=CC=C1(3052)',
    structure = SMILES('O=CCC1[CH]C=C[CH]C=C1'),
    E0 = (194.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07852,-0.000769254,0.000274894,-3.96117e-07,1.64105e-10,23601.5,31.6641], Tmin=(100,'K'), Tmax=(910.155,'K')), NASAPolynomial(coeffs=[42.5115,-0.0137838,1.76911e-05,-3.61604e-09,2.3052e-13,9056.42,-202.797], Tmin=(910.155,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(194.854,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(1,4-Cycloheptadiene) + radical(Allyl_S) + radical(C=CCJC=C)"""),
)

species(
    label = 'O=CCC1[C]=C[CH]C=CC1(3118)',
    structure = SMILES('O=CCC1[C]=C[CH]C=CC1'),
    E0 = (291.583,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.852959,0.00898906,0.000241377,-3.58629e-07,1.50266e-10,35239.1,34.054], Tmin=(100,'K'), Tmax=(908.722,'K')), NASAPolynomial(coeffs=[41.2312,-0.0117171,1.63505e-05,-3.38112e-09,2.17214e-13,21416.9,-192.557], Tmin=(908.722,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(291.583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(1,4-Cycloheptadiene) + radical(Cds_S) + radical(C=CCJC=C)"""),
)

species(
    label = 'O=[C]C[C]1C=CCC=CC1(3097)',
    structure = SMILES('O=[C]C[C]1C=CCC=CC1'),
    E0 = (243.958,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.781541,0.0142152,0.000219546,-3.32706e-07,1.40422e-10,29510.2,33.9322], Tmin=(100,'K'), Tmax=(907.475,'K')), NASAPolynomial(coeffs=[39.8018,-0.00999606,1.52888e-05,-3.1953e-09,2.06608e-13,16343.2,-184.057], Tmin=(907.475,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.958,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(1,4-Cycloheptadiene) + radical(CCCJ=O) + radical(Allyl_T)"""),
)

species(
    label = 'O=[C]CC1[CH]C=CCC=C1(2246)',
    structure = SMILES('O=[C]CC1C=C[CH]CC=C1'),
    E0 = (251.272,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.46551,0.0260414,0.000181888,-2.89927e-07,1.23965e-10,30396.7,32.5997], Tmin=(100,'K'), Tmax=(909.75,'K')), NASAPolynomial(coeffs=[38.5536,-0.00634623,1.25698e-05,-2.64223e-09,1.6908e-13,17876.7,-178.273], Tmin=(909.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(251.272,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(1,4-Cycloheptadiene) + radical(Allyl_S) + radical(CCCJ=O)"""),
)

species(
    label = 'O=CCC1C=C[CH]C=[C]C1(3119)',
    structure = SMILES('O=CCC1C=C[CH]C=[C]C1'),
    E0 = (291.583,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.852678,0.00899259,0.000241364,-3.5861e-07,1.50257e-10,35239.1,34.055], Tmin=(100,'K'), Tmax=(908.732,'K')), NASAPolynomial(coeffs=[41.232,-0.0117184,1.63513e-05,-3.38131e-09,2.1723e-13,21416.6,-192.561], Tmin=(908.732,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(291.583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(1,4-Cycloheptadiene) + radical(C=CCJC=C) + radical(Cds_S)"""),
)

species(
    label = 'O=CCC1C=[C][CH]C=CC1(3120)',
    structure = SMILES('O=CCC1C=[C]C=C[CH]C1'),
    E0 = (274.627,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.716781,0.0951449,-7.69721e-05,3.24444e-08,-5.54937e-12,33206.8,17.8887], Tmin=(100,'K'), Tmax=(1383.4,'K')), NASAPolynomial(coeffs=[18.918,0.0383722,-1.54142e-05,2.77936e-09,-1.88463e-13,27774.2,-83.21], Tmin=(1383.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.627,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3-Cycloheptadiene) + radical(Allyl_S) + radical(C=CJC=C)"""),
)

species(
    label = 'O=[C][CH]C1C=CCC=CC1(3099)',
    structure = SMILES('O=[C][CH]C1C=CCC=CC1'),
    E0 = (279.507,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,1855,455,950,3025,407.5,1350,352.5,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.665341,0.0132882,0.000231355,-3.52058e-07,1.49181e-10,33793.3,35.1638], Tmin=(100,'K'), Tmax=(906.019,'K')), NASAPolynomial(coeffs=[42.7033,-0.0151399,1.82169e-05,-3.76476e-09,2.45128e-13,19725.2,-199.095], Tmin=(906.019,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.507,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(1,4-Cycloheptadiene) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = 'O=CCC1C=C[CH][C]=CC1(3121)',
    structure = SMILES('O=CCC1[CH]C=C[C]=CC1'),
    E0 = (273.557,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.666728,0.0948523,-7.65378e-05,3.22258e-08,-5.52018e-12,33075.6,17.6733], Tmin=(100,'K'), Tmax=(1378.05,'K')), NASAPolynomial(coeffs=[18.535,0.0391159,-1.58685e-05,2.87532e-09,-1.95502e-13,27783.5,-81.1214], Tmin=(1378.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(273.557,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3-Cycloheptadiene) + radical(Allyl_S) + radical(C=CJC=C)"""),
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
    label = 'O=[C]CC1C=CC2[CH]C2C1(3122)',
    structure = SMILES('O=[C]CC1C=CC2[CH]C2C1'),
    E0 = (353.499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.512168,0.0629273,-4.68967e-06,-2.94213e-08,1.2934e-11,42653.3,31.6023], Tmin=(100,'K'), Tmax=(1089.59,'K')), NASAPolynomial(coeffs=[13.709,0.0436376,-1.82742e-05,3.45025e-09,-2.43372e-13,38046.7,-41.1395], Tmin=(1089.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(353.499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + polycyclic(s2_3_6_ene_1) + radical(cyclopropane) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C]CC1CC=CC2[CH]C21(3123)',
    structure = SMILES('O=[C]CC1CC=CC2[CH]C21'),
    E0 = (352.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.511469,0.0632008,-6.08704e-06,-2.74687e-08,1.21221e-11,42524.5,31.5712], Tmin=(100,'K'), Tmax=(1099.63,'K')), NASAPolynomial(coeffs=[13.5961,0.0439487,-1.84899e-05,3.49167e-09,-2.46e-13,37933.1,-40.5894], Tmin=(1099.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.43,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + polycyclic(s2_3_6_ene_1) + radical(CCCJ=O) + radical(cyclopropane)"""),
)

species(
    label = 'O=C1CC2C=C[CH][CH]C1C2(3124)',
    structure = SMILES('O=C1CC2[CH][CH]C=CC1C2'),
    E0 = (143.228,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.179081,0.0908091,-7.11651e-05,3.02088e-08,-5.58391e-12,17358,-2.07354], Tmin=(100,'K'), Tmax=(1206.55,'K')), NASAPolynomial(coeffs=[11.2706,0.0540381,-2.54511e-05,4.94998e-09,-3.50247e-13,14681.5,-57.6663], Tmin=(1206.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(143.228,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclopentane) + ring(Cycloheptene) + ring(Cyclopentane) + radical(Cs_S) + radical(Allyl_S)"""),
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
    label = 'O=C=CC1C=CCC=CC1(3109)',
    structure = SMILES('O=C=CC1C=CCC=CC1'),
    E0 = (75.0592,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.419526,0.0210405,0.000208888,-3.27334e-07,1.39667e-10,9210.5,32.6308], Tmin=(100,'K'), Tmax=(908.507,'K')), NASAPolynomial(coeffs=[42.524,-0.0134266,1.66315e-05,-3.41838e-09,2.20425e-13,-4667.95,-200.735], Tmin=(908.507,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(75.0592,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCd(CCO)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + ring(1,4-Cycloheptadiene)"""),
)

species(
    label = 'O=CCC1C=C=CC=CC1(3126)',
    structure = SMILES('O=CCC1C=C=CC=CC1'),
    E0 = (157.418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.201584,0.0678699,-6.7069e-06,-3.41249e-08,1.60594e-11,19083.3,27.0957], Tmin=(100,'K'), Tmax=(1046.56,'K')), NASAPolynomial(coeffs=[16.0675,0.0410554,-1.67562e-05,3.15995e-09,-2.24474e-13,13909.9,-59.0204], Tmin=(1046.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.418,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene)"""),
)

species(
    label = 'O=CCC1C=CC=C=CC1(3127)',
    structure = SMILES('O=CCC1C=CC=C=CC1'),
    E0 = (157.418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.201584,0.0678699,-6.7069e-06,-3.41249e-08,1.60594e-11,19083.3,27.0957], Tmin=(100,'K'), Tmax=(1046.56,'K')), NASAPolynomial(coeffs=[16.0675,0.0410554,-1.67562e-05,3.15995e-09,-2.24474e-13,13909.9,-59.0204], Tmin=(1046.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.418,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene)"""),
)

species(
    label = 'C=C([O])C1[CH]C=CC=CC1(2239)',
    structure = SMILES('C=C([O])C1C=C[CH]C=CC1'),
    E0 = (187.021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.554994,0.0115722,0.000248281,-3.78125e-07,1.60819e-10,22677.8,32.2663], Tmin=(100,'K'), Tmax=(903.379,'K')), NASAPolynomial(coeffs=[45.8044,-0.0203488,2.1608e-05,-4.45426e-09,2.92929e-13,7629.38,-219.478], Tmin=(903.379,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(187.021,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,4-Cycloheptadiene) + radical(C=C(C)OJ) + radical(C=CCJC=C)"""),
)

species(
    label = 'O=C1CC2C=CC1C=CC2(2250)',
    structure = SMILES('O=C1CC2C=CC1C=CC2'),
    E0 = (-76.4957,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.511618,0.10077,-9.94395e-05,5.4281e-08,-1.22191e-11,-9039.08,0.0667863], Tmin=(100,'K'), Tmax=(1060.33,'K')), NASAPolynomial(coeffs=[14.629,0.043654,-1.86403e-05,3.47997e-09,-2.4156e-13,-12249.9,-73.8649], Tmin=(1060.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-76.4957,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclohexane) + ring(1,4-Cycloheptadiene) + ring(Cyclohexene)"""),
)

species(
    label = '[CH2]C1C=C[CH]C=CC1(636)',
    structure = SMILES('[CH2]C1C=C[CH]C=CC1'),
    E0 = (365.403,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85714,-0.0209364,0.000314066,-4.38338e-07,1.81372e-10,44089.1,27.4842], Tmin=(100,'K'), Tmax=(900.881,'K')), NASAPolynomial(coeffs=[42.8094,-0.0249472,2.46653e-05,-5.07268e-09,3.35443e-13,29494.6,-205.861], Tmin=(900.881,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(365.403,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cycloheptadiene) + radical(Isobutyl) + radical(C=CCJC=C)"""),
)

species(
    label = 'CO(28)',
    structure = SMILES('[C-]#[O+]'),
    E0 = (-119.219,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2084.51],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0101,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(762.44,'J/mol'), sigma=(3.69,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.5971,-0.00102424,2.83336e-06,-1.75825e-09,3.42587e-13,-14343.2,3.45822], Tmin=(100,'K'), Tmax=(1669.93,'K')), NASAPolynomial(coeffs=[2.92796,0.00181931,-8.35308e-07,1.51269e-10,-9.88872e-15,-14292.7,6.51157], Tmin=(1669.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-119.219,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""CO""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'C=CC1[CH]C1C=CC[C]=O(3128)',
    structure = SMILES('C=CC1[CH]C1C=CC[C]=O'),
    E0 = (437.158,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1855,455,950,2750,2950,3150,900,1000,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.209266,0.0640707,7.62422e-06,-5.23051e-08,2.28116e-11,52731.1,38.0457], Tmin=(100,'K'), Tmax=(1040.18,'K')), NASAPolynomial(coeffs=[18.7728,0.0364595,-1.56849e-05,3.09247e-09,-2.26693e-13,46501,-63.6275], Tmin=(1040.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(437.158,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(CCCJ=O) + radical(cyclopropane)"""),
)

species(
    label = '[C]=O(353)',
    structure = SMILES('[C]=O'),
    E0 = (439.086,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3054.48],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0101,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.08916,0.00200416,-1.61661e-05,2.55058e-08,-1.16424e-11,52802.7,4.52505], Tmin=(100,'K'), Tmax=(856.11,'K')), NASAPolynomial(coeffs=[0.961625,0.00569045,-3.48044e-06,7.19202e-10,-5.08041e-14,53738.7,21.4663], Tmin=(856.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(439.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CdCdJ2_triplet)"""),
)

species(
    label = 'O=[C]CC1=CC=CC=CC1(3129)',
    structure = SMILES('O=[C]CC1=CC=CC=CC1'),
    E0 = (187.147,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (133.167,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.684062,0.0830321,-2.95454e-05,-2.84837e-08,1.80634e-11,22694.8,29.8705], Tmin=(100,'K'), Tmax=(1010.33,'K')), NASAPolynomial(coeffs=[24.7192,0.0263936,-1.06859e-05,2.11342e-09,-1.58032e-13,15319.2,-104.044], Tmin=(1010.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(187.147,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3,5-Cycloheptatriene) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C]C[C]1CC=CC=CC1(3130)',
    structure = SMILES('O=[C]C[C]1CC=CC=CC1'),
    E0 = (278.832,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.184279,0.0968461,-9.31987e-05,5.13062e-08,-1.19938e-11,33682.4,19.4168], Tmin=(100,'K'), Tmax=(1006.14,'K')), NASAPolynomial(coeffs=[11.8018,0.0491942,-2.21566e-05,4.23348e-09,-2.9738e-13,31270.5,-38.4822], Tmin=(1006.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(278.832,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3-Cycloheptadiene) + radical(Tertalkyl) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C]CC1C[C]=CC=CC1(3131)',
    structure = SMILES('O=[C]CC1C[C]=CC=CC1'),
    E0 = (331.251,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.681245,0.0986676,-8.80647e-05,4.13144e-08,-7.86783e-12,40012.4,20.5244], Tmin=(100,'K'), Tmax=(1250.43,'K')), NASAPolynomial(coeffs=[18.0474,0.0387564,-1.61958e-05,2.99751e-09,-2.0707e-13,35328.6,-74.0161], Tmin=(1250.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(331.251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3-Cycloheptadiene) + radical(CCCJ=O) + radical(Cds_S)"""),
)

species(
    label = 'O=[C][CH]C1CC=CC=CC1(3132)',
    structure = SMILES('O=[C][CH]C1CC=CC=CC1'),
    E0 = (260.939,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,1855,455,950,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.768561,0.0910969,-5.79661e-05,5.97572e-09,5.35525e-12,31567.5,20.1236], Tmin=(100,'K'), Tmax=(1016.32,'K')), NASAPolynomial(coeffs=[20.7772,0.0332851,-1.24719e-05,2.26107e-09,-1.58103e-13,25794.2,-91.0277], Tmin=(1016.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.939,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3-Cycloheptadiene) + radical(CCCJ=O) + radical(CCJCHO)"""),
)

species(
    label = 'O=[C]CC1CC=[C]C=CC1(3133)',
    structure = SMILES('O=[C]CC1CC=[C]C=CC1'),
    E0 = (292.405,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.744471,0.100832,-9.32084e-05,4.61776e-08,-9.30409e-12,35342,19.8745], Tmin=(100,'K'), Tmax=(1186.94,'K')), NASAPolynomial(coeffs=[17.3187,0.0399591,-1.62812e-05,2.9704e-09,-2.03621e-13,31054,-70.3657], Tmin=(1186.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(292.405,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3-Cycloheptadiene) + radical(C=CJC=C) + radical(CCCJ=O)"""),
)

species(
    label = 'O=CCC1[CH]C=[C]C=CC1(3134)',
    structure = SMILES('O=CCC1[CH]C=[C]C=CC1'),
    E0 = (273.557,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.666727,0.0948523,-7.65378e-05,3.22258e-08,-5.52018e-12,33075.6,17.6733], Tmin=(100,'K'), Tmax=(1378.05,'K')), NASAPolynomial(coeffs=[18.535,0.0391159,-1.58685e-05,2.87532e-09,-1.95502e-13,27783.5,-81.1214], Tmin=(1378.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(273.557,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3-Cycloheptadiene) + radical(C=CJC=C) + radical(Allyl_S)"""),
)

species(
    label = 'O=[C]CC1C[CH]C2C=CC21(3135)',
    structure = SMILES('O=[C]CC1C[CH]C2C=CC21'),
    E0 = (368.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.677141,0.058181,8.01385e-06,-4.17217e-08,1.69567e-11,44492,30.784], Tmin=(100,'K'), Tmax=(1075.7,'K')), NASAPolynomial(coeffs=[13.4603,0.0439503,-1.85819e-05,3.54196e-09,-2.5176e-13,39814.9,-40.7766], Tmin=(1075.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(368.828,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + Estimated bicyclic component: polycyclic(s2_4_5_ane) - ring(Cyclopentane) - ring(Cyclobutane) + ring(Cyclopentane) + ring(Cyclobutene) + radical(CCCJ=O) + radical(Cs_S)"""),
)

species(
    label = 'O=[C]CC1CC2[CH]C=CC21(3136)',
    structure = SMILES('O=[C]CC1CC2[CH]C=CC21'),
    E0 = (272.677,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.686884,0.0539011,3.05778e-05,-7.07241e-08,2.83846e-11,32931.3,27.9195], Tmin=(100,'K'), Tmax=(1021.48,'K')), NASAPolynomial(coeffs=[15.3264,0.0414093,-1.6916e-05,3.24119e-09,-2.34118e-13,27601.5,-54.4681], Tmin=(1021.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(272.677,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + polycyclic(s2_4_5_ene_1) + radical(cyclopentene-allyl) + radical(CCCJ=O)"""),
)

species(
    label = 'O=C1CC2[CH][CH]C1C=CC2(3137)',
    structure = SMILES('O=C1CC2[CH][CH]C1C=CC2'),
    E0 = (202.017,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.102498,0.0916567,-7.43453e-05,3.28651e-08,-6.27969e-12,24432.2,1.80633], Tmin=(100,'K'), Tmax=(1182.45,'K')), NASAPolynomial(coeffs=[11.7213,0.0523525,-2.44858e-05,4.75431e-09,-3.36357e-13,21684.5,-56.195], Tmin=(1182.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(202.017,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclohexane) + ring(Cycloheptene) + ring(Cyclohexane) + radical(Cs_S) + radical(CCJCC=O)"""),
)

species(
    label = '[CH2]C=CC=CC=CC[C]=O(3138)',
    structure = SMILES('[CH2]C=CC=CC=CC[C]=O'),
    E0 = (252.396,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3001,3007,3013,3019,3025,975,980,985,990,995,1000,1300,1315,1330,1345,1360,1375,400,420,440,460,480,500,1630,1640,1650,1660,1670,1680,1855,455,950,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.498941,0.080292,-2.40839e-05,-2.98873e-08,1.76411e-11,30534.6,37.2685], Tmin=(100,'K'), Tmax=(1012.65,'K')), NASAPolynomial(coeffs=[22.1638,0.0316558,-1.25979e-05,2.41778e-09,-1.76351e-13,23848.5,-82.7009], Tmin=(1012.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(252.396,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(C=CC=CCJ)"""),
)

species(
    label = 'O=CCC1=CC=CC=CC1(3139)',
    structure = SMILES('O=CCC1=CC=CC=CC1'),
    E0 = (27.1859,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.625914,0.0802202,-1.49333e-05,-4.23548e-08,2.22009e-11,3455.12,28.8444], Tmin=(100,'K'), Tmax=(1017.96,'K')), NASAPolynomial(coeffs=[24.1304,0.0308809,-1.287e-05,2.5565e-09,-1.90491e-13,-4068.86,-103.231], Tmin=(1017.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(27.1859,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3,5-Cycloheptatriene)"""),
)

species(
    label = 'O=C=CC1CC=CC=CC1(3140)',
    structure = SMILES('O=C=CC1CC=CC=CC1'),
    E0 = (58.4969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.824035,0.0920542,-6.09543e-05,1.05787e-08,3.06872e-12,7221.3,18.7953], Tmin=(100,'K'), Tmax=(1058.92,'K')), NASAPolynomial(coeffs=[21.0187,0.0340235,-1.34267e-05,2.48722e-09,-1.74972e-13,1222.97,-94.3137], Tmin=(1058.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(58.4969,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene)"""),
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
    label = 'O=[C]C[CH]C1C=CC=CC1(3070)',
    structure = SMILES('O=[C]C[CH]C1C=CC=CC1'),
    E0 = (286.123,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1855,455,950,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.207034,0.066462,-1.32103e-06,-4.3768e-08,2.07384e-11,34564.1,34.7209], Tmin=(100,'K'), Tmax=(1007.35,'K')), NASAPolynomial(coeffs=[17.4072,0.0366314,-1.41822e-05,2.65172e-09,-1.89669e-13,29147,-58.0737], Tmin=(1007.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(286.123,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(469.768,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3-Cyclohexadiene) + radical(CCJCC=O) + radical(CCCJ=O)"""),
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
    label = 'O=[C]CC1C=CC=CC=C1(3041)',
    structure = SMILES('O=[C]CC1C=CC=CC=C1'),
    E0 = (185.315,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2750,2850,1437.5,1250,1305,750,350,1855,455,950,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (133.167,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.01898,0.101714,-9.72871e-05,4.77547e-08,-9.28807e-12,22476.2,28.9865], Tmin=(100,'K'), Tmax=(1247.38,'K')), NASAPolynomial(coeffs=[21.1151,0.030736,-1.19346e-05,2.13771e-09,-1.45499e-13,16954.3,-82.6902], Tmin=(1247.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(185.315,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3,5-Cycloheptatriene) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C]CC1C=CC=[C]CC1(3047)',
    structure = SMILES('O=[C]CC1C=CC=[C]CC1'),
    E0 = (332.321,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.724695,0.0988952,-8.83232e-05,4.13615e-08,-7.84243e-12,40143.2,20.7151], Tmin=(100,'K'), Tmax=(1258.12,'K')), NASAPolynomial(coeffs=[18.4231,0.0380176,-1.57414e-05,2.90102e-09,-1.99956e-13,35325.2,-76.0586], Tmin=(1258.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(332.321,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3-Cycloheptadiene) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C]C[C]1C=CC=CCC1(3043)',
    structure = SMILES('O=[C]CC1=C[CH]C=CCC1'),
    E0 = (213.716,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2850,1437.5,1250,1305,750,350,1855,455,950,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.901371,0.000438657,0.000277427,-4.01269e-07,1.657e-10,25878.8,33.8973], Tmin=(100,'K'), Tmax=(916.444,'K')), NASAPolynomial(coeffs=[44.9195,-0.0162834,1.77025e-05,-3.48507e-09,2.15069e-13,10445,-214.811], Tmin=(916.444,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(213.716,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(1,4-Cycloheptadiene) + radical(C=CCJC=C) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C]CC1C=C[C]=CCC1(3051)',
    structure = SMILES('O=[C]CC1C=C[C]=CCC1'),
    E0 = (293.475,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.78392,0.101018,-9.33457e-05,4.60963e-08,-9.23423e-12,35472.6,20.0504], Tmin=(100,'K'), Tmax=(1195.58,'K')), NASAPolynomial(coeffs=[17.6816,0.0392383,-1.58355e-05,2.87572e-09,-1.96642e-13,31057.2,-72.3337], Tmin=(1195.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293.475,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3-Cycloheptadiene) + radical(CCCJ=O) + radical(C=CJC=C)"""),
)

species(
    label = 'O=[C][CH]C1C=CC=CCC1(3045)',
    structure = SMILES('O=[C][CH]C1C=CC=CCC1'),
    E0 = (262.009,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,1855,455,950,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.76438,0.0907713,-5.63336e-05,3.64004e-09,6.36481e-12,31696.2,20.143], Tmin=(100,'K'), Tmax=(1008.62,'K')), NASAPolynomial(coeffs=[20.9452,0.0328888,-1.22106e-05,2.20944e-09,-1.54664e-13,25881.7,-91.8938], Tmin=(1008.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(262.009,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3-Cycloheptadiene) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C]CC1[C]=CC=CCC1(3046)',
    structure = SMILES('O=[C]CC1[C]=CC=CCC1'),
    E0 = (332.321,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2850,1437.5,1250,1305,750,350,1855,455,950,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.724695,0.0988952,-8.83232e-05,4.13615e-08,-7.84243e-12,40143.2,20.7151], Tmin=(100,'K'), Tmax=(1258.12,'K')), NASAPolynomial(coeffs=[18.4231,0.0380176,-1.57414e-05,2.90102e-09,-1.99956e-13,35325.2,-76.0586], Tmin=(1258.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(332.321,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3-Cycloheptadiene) + radical(CCCJ=O) + radical(Cds_S)"""),
)

species(
    label = 'O=[C]CC1C=[C]C=CCC1(3050)',
    structure = SMILES('O=[C]CC1C=[C]C=CCC1'),
    E0 = (293.475,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2850,1437.5,1250,1305,750,350,1855,455,950,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.78392,0.101018,-9.33457e-05,4.60963e-08,-9.23423e-12,35472.6,20.0504], Tmin=(100,'K'), Tmax=(1195.58,'K')), NASAPolynomial(coeffs=[17.6816,0.0392383,-1.58355e-05,2.87572e-09,-1.96642e-13,31057.2,-72.3337], Tmin=(1195.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293.475,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3-Cycloheptadiene) + radical(C=CJC=C) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C]CC1[CH]C2C=CC2C1(3141)',
    structure = SMILES('O=[C]CC1[CH]C2C=CC2C1'),
    E0 = (368.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.677141,0.058181,8.01385e-06,-4.17217e-08,1.69567e-11,44492,30.784], Tmin=(100,'K'), Tmax=(1075.7,'K')), NASAPolynomial(coeffs=[13.4603,0.0439503,-1.85819e-05,3.54196e-09,-2.5176e-13,39814.9,-40.7766], Tmin=(1075.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(368.828,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + Estimated bicyclic component: polycyclic(s2_4_5_ane) - ring(Cyclopentane) - ring(Cyclobutane) + ring(Cyclopentane) + ring(Cyclobutene) + radical(CCCJ=O) + radical(Cs_S)"""),
)

species(
    label = 'O=CCC1C=CC=CC=C1(3068)',
    structure = SMILES('O=CCC1C=CC=CC=C1'),
    E0 = (25.354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.07018,0.100108,-8.64808e-05,3.82098e-08,-6.7315e-12,3241.38,28.358], Tmin=(100,'K'), Tmax=(1364.03,'K')), NASAPolynomial(coeffs=[21.6808,0.0333911,-1.31141e-05,2.35219e-09,-1.59529e-13,-2965.25,-88.4652], Tmin=(1364.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(25.354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3,5-Cycloheptatriene)"""),
)

species(
    label = 'O=C=CC1C=CC=CCC1(3067)',
    structure = SMILES('O=C=CC1C=CC=CCC1'),
    E0 = (57.5614,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.29793,0.101661,-8.8503e-05,3.8992e-08,-6.78531e-12,7126.22,18.6585], Tmin=(100,'K'), Tmax=(1390.83,'K')), NASAPolynomial(coeffs=[23.5928,0.0300749,-1.12965e-05,1.98413e-09,-1.33094e-13,202.546,-109.636], Tmin=(1390.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(57.5614,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCd(CCO)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene)"""),
)

species(
    label = 'O=[C]CC1=CC=CC[CH]C1(3142)',
    structure = SMILES('O=[C]CC1=CC=CC[CH]C1'),
    E0 = (288.952,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.536172,0.0896619,-6.55424e-05,2.33129e-08,-3.31636e-12,34924.2,21.8249], Tmin=(100,'K'), Tmax=(1640.95,'K')), NASAPolynomial(coeffs=[22.122,0.0344299,-1.50546e-05,2.80129e-09,-1.91392e-13,27488,-98.71], Tmin=(1640.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(288.952,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3-Cycloheptadiene) + radical(RCCJCC) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C1C=CC=CC1C[C]=O(3071)',
    structure = SMILES('[CH2]C1C=CC=CC1C[C]=O'),
    E0 = (283.332,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1855,455,950,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0921783,0.0699942,-8.95123e-06,-3.87872e-08,2.01368e-11,34231.8,33.4785], Tmin=(100,'K'), Tmax=(980.026,'K')), NASAPolynomial(coeffs=[17.5951,0.0355214,-1.2767e-05,2.2965e-09,-1.61581e-13,29026,-59.6668], Tmin=(980.026,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(283.332,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(469.768,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3-Cyclohexadiene) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = 'O=C1CC2C=CC=CC1C2(3143)',
    structure = SMILES('O=C1CC2C=CC=CC1C2'),
    E0 = (-95.5158,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.85059,0.17071,-0.000363423,3.93675e-07,-1.53891e-10,-11317.7,-15.4845], Tmin=(100,'K'), Tmax=(856.811,'K')), NASAPolynomial(coeffs=[-8.39068,0.0972592,-5.27943e-05,1.03393e-08,-7.12965e-13,-6380.17,37.3304], Tmin=(856.811,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-95.5158,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclopentane) + ring(1,3-Cycloheptadiene) + ring(Cyclopentane)"""),
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
    E0 = (213.702,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (234.622,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (298.783,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (399.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (545.966,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (545.966,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (331.383,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (372.972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (506.719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (506.719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (327.638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (499.188,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (499.188,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (316.286,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (392.761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (359.436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (443.109,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (328.161,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (318.925,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (341.505,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (335.865,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (499.135,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (459.023,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (459.023,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (256.784,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (263.194,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (252.927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (231.484,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (231.484,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (459.303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (221.233,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (271.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (584.207,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (804.489,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (403.876,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (327.521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (404.77,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (531.447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (385.483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (439.682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (499.962,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (440.751,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (393.107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (239.643,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (279.592,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (291.949,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (302.671,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (334.974,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (446.058,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (221.61,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (397.107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (532.517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (340.491,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (440.752,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (367.864,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (376.63,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (554.477,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (440.751,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (231.484,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (252.927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (405.507,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (443.267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (221.024,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    products = ['C7H8(11)', 'CH2CO(43)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    products = ['O=C1CC2C=C[CH]C1[CH]C2(3111)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(414000,'s^-1'), n=1.14, Ea=(20.92,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7plus;doublebond_intra_HNd_pri;radadd_intra] for rate rule [R8;doublebond_intra_HNd_pri;radadd_intra_CO]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    products = ['O=C1CC2[CH]C1[CH]C=CC2(3112)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.56068e+09,'s^-1'), n=0.543284, Ea=(85.0813,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7plus;doublebond_intra;radadd_intra] for rate rule [R9;doublebond_intra;radadd_intra_CO]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(19)', 'O=C=CC1C=C[CH]C=CC1(3113)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.82e-16,'cm^3/(molecule*s)'), n=1.61, Ea=(10.992,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Ck;HJ] for rate rule [Cds-CsH_Ck;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(19)', 'O=[C]CC1C=CC=C=CC1(3114)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.149e+09,'cm^3/(mol*s)'), n=1.595, Ea=(16.7946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-OneDeH;HJ] for rate rule [Ca_Cds-CdH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(19)', 'O=[C]CC1C=C=CC=CC1(3115)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.149e+09,'cm^3/(mol*s)'), n=1.595, Ea=(16.7946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-OneDeH;HJ] for rate rule [Ca_Cds-CdH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]1C=C[CH]CC=C1(125)', 'CH2CO(43)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.0038925,'m^3/(mol*s)'), n=2.41, Ea=(59.0153,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cdd;CsJ-CdCsH] for rate rule [Cds-HH_Ck;CsJ-CdCsH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    products = ['O=C[CH]C1C=C[CH]C=CC1(3116)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.74568e+08,'s^-1'), n=1.384, Ea=(159.27,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['O=[C]CC1C=CC[C]=CC1(3101)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.448e+10,'s^-1'), n=0.82, Ea=(156.9,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 154 used for R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['O=[C]CC1C=[C]CC=CC1(3102)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.448e+10,'s^-1'), n=0.82, Ea=(156.9,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 154 used for R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    products = ['O=CC[C]1C=C[CH]C=CC1(3117)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(46592.1,'s^-1'), n=2.43015, Ea=(113.936,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;Cs_H_out_noH] for rate rule [R3H_SS_Cs;CO_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['O=[C]CC1C=CCC=[C]C1(3098)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.182e+10,'s^-1'), n=0.86, Ea=(149.369,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['O=[C]CC1[C]=CCC=CC1(3100)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.182e+10,'s^-1'), n=0.86, Ea=(149.369,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    products = ['O=CCC1[CH][CH]C=CC=C1(3052)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(657459,'s^-1'), n=1.87237, Ea=(102.584,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R4H_SSS;CO_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O=CCC1[C]=C[CH]C=CC1(3118)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['O=[C]C[C]1C=CCC=CC1(3097)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.30755e+06,'s^-1'), n=1.805, Ea=(115.478,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;C_rad_out_single;Cs_H_out_H/Cd] for rate rule [R4H_SDS;C_rad_out_Cs2;Cs_H_out_H/Cd]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['O=[C]CC1[CH]C=CCC=C1(2246)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(285077,'s^-1'), n=1.88833, Ea=(191.836,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H;C_rad_out_H/NonDeC;Cs_H_out_H/Cd] + [R4H_SDS;C_rad_out_single;Cs_H_out_H/Cd] + [R4H_SDS;C_rad_out_H/NonDeC;Cs_H_out_1H] for rate rule [R4H_SDS;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['O=CCC1C=C[CH]C=[C]C1(3119)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(116983,'s^-1'), n=1.68042, Ea=(36.5782,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;Y_rad_out;XH_out] for rate rule [R5H_CCC_O;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 3.31662479036
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['O=CCC1C=[C][CH]C=CC1(3120)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(380071,'s^-1'), n=1.62386, Ea=(44.2978,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;XH_out] for rate rule [R5H_DSSS;Cd_rad_out;CO_H_out]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['O=[C][CH]C1C=CCC=CC1(3099)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(60205.5,'s^-1'), n=1.86417, Ea=(61.9987,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R5H_SSMS;Y_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['O=CCC1C=C[CH][C]=CC1(3121)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(13085.7,'s^-1'), n=1.78967, Ea=(62.3074,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_RSSSR;Y_rad_out;XH_out] for rate rule [R6H_DSSSS;Cd_rad_out;CO_H_out]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]1C=C[CH]CC=C1(125)', '[CH2][C]=O(360)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(5.42928e+07,'m^3/(mol*s)'), n=0.107721, Ea=(5.76381,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [C_rad/H/CdCs;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    products = ['O=[C]CC1C=CC2[CH]C2C1(3122)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.00063e+13,'s^-1'), n=-0.283562, Ea=(245.321,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_beta;doublebond_intra_pri_HNd_Cs;radadd_intra_csHDe] for rate rule [Rn0c7_beta_short;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCd]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    products = ['O=[C]CC1CC=CC2[CH]C21(3123)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.00063e+13,'s^-1'), n=-0.283562, Ea=(245.321,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_beta;doublebond_intra_pri_HNd_Cs;radadd_intra_csHDe] for rate rule [Rn0c7_beta_short;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCd]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    products = ['O=C1CC2C=C[CH][CH]C1C2(3124)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.95208e+10,'s^-1'), n=0.228823, Ea=(43.0818,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn2c7_beta_long;doublebond_intra_pri_HNd_Cs;radadd_intra_CO]
Euclidian distance = 3.16227766017
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    products = ['O=C1CC2CC=C[CH][CH]C12(3125)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.99832e+10,'s^-1'), n=0.37247, Ea=(49.4918,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri;radadd_intra] for rate rule [Rn2c7_alpha_long;doublebond_intra_pri;radadd_intra_CO]
Euclidian distance = 3.16227766017
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    products = ['O=C=CC1C=CCC=CC1(3109)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.70659e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    products = ['O=CCC1C=C=CC=CC1(3126)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.07e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_De] for rate rule [R5radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    products = ['O=CCC1C=CC=C=CC1(3127)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.07e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad_De] for rate rule [R6radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    products = ['C=C([O])C1[CH]C=CC=CC1(2239)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction31',
    reactants = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    products = ['O=C1CC2C=CC1C=CC2(2250)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Cpri_rad_out_single] for rate rule [R6_SSSDS;Y_rad_out;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C1C=C[CH]C=CC1(636)', 'CO(28)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1230.59,'m^3/(mol*s)'), n=1.0523, Ea=(25.6182,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [COm;C_rad/H2/Cs]
Euclidian distance = 0
family: R_Addition_COm"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=CC1[CH]C1C=CC[C]=O(3128)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(6.21184e+10,'s^-1'), n=0.288169, Ea=(147.049,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_5_unsaturated_hexane] for rate rule [1_5_hexadiene]
Euclidian distance = 1.0
family: 6_membered_central_C-C_shift"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C1C=C[CH]C=CC1(636)', '[C]=O(353)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(2.23625e+06,'m^3/(mol*s)'), n=0.36814, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction35',
    reactants = ['H(19)', 'O=[C]CC1=CC=CC=CC1(3129)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(7.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(4.93712,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2570 used for Cds-CsCs_Cds-CdH;HJ
Exact match found for rate rule [Cds-CsCs_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C7H8(11)', '[CH2][C]=O(360)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.11106,'m^3/(mol*s)'), n=1.92886, Ea=(1.83261,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cds-CsH_Cds-OneDeH;CJ] + [Cds-CsH_Cds-(Cd-Cd-Cd)H;YJ] for rate rule [Cds-CsH_Cds-(Cd-Cd-Cd)H;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction37',
    reactants = ['O=[C]C[C]1CC=CC=CC1(3130)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(5.88e+08,'s^-1'), n=1.27, Ea=(125.938,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_Cs2;Cs_H_out_H/Cd] for rate rule [R2H_S;C_rad_out_Cs2;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['O=[C]CC1C[C]=CC=CC1(3131)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.9054e+11,'s^-1'), n=0.853, Ea=(200.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['O=[C][CH]C1CC=CC=CC1(3132)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(2.59423e+07,'s^-1'), n=1.52333, Ea=(124.544,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['O=[C]CC1CC=[C]C=CC1(3133)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_single;Cs_H_out_H/NonDeC] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/(NonDeC/Cs/Cs)]
Euclidian distance = 2.82842712475
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['O=CCC1[CH]C=[C]C=CC1(3134)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.61756e+09,'s^-1'), n=1.29078, Ea=(226.405,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleDe_Cd;XH_out] for rate rule [R6HJ_2;Cd_rad_out_singleDe_Cd;CO_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    products = ['O=[C]CC1C[CH]C2C=CC21(3135)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(5.65845e+18,'s^-1'), n=-1.57151, Ea=(227.049,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 131 used for Rn0c7_gamma_long_SDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCd
Exact match found for rate rule [Rn0c7_gamma_long_SDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCd]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction43',
    reactants = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    products = ['O=[C]CC1CC2[CH]C=CC21(3136)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(6.33004e+14,'s^-1'), n=-0.792922, Ea=(179.405,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 112 used for Rn0c7_gamma_long_SSS_D;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs
Exact match found for rate rule [Rn0c7_gamma_long_SSS_D;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction44',
    reactants = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    products = ['O=C1CC2[CH][CH]C1C=CC2(3137)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(3.125e+09,'s^-1'), n=0.76, Ea=(25.9408,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_cyclic;doublebond_intra_pri_HCd;radadd_intra] for rate rule [Rn2c7_gamma_short;doublebond_intra_pri_HCd;radadd_intra_CO]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]C=CC=CC=CC[C]=O(3138)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(110000,'s^-1'), n=1.18, Ea=(27.196,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R7_linear;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H] for rate rule [R7_SDSD_D;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction46',
    reactants = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    products = ['O=CCC1=CC=CC=CC1(3139)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction47',
    reactants = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    products = ['O=C=CC1CC=CC=CC1(3140)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction48',
    reactants = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(4.25e+09,'s^-1'), n=0.991, Ea=(45.9529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_3_pentadiene;CH_end;CdHC_2] for rate rule [1_3_pentadiene;CH(CJ)_1;CdHC_2]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction49',
    reactants = ['O=[C]C[CH]C1C=CC=CC1(3070)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction50',
    reactants = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    products = ['O=C1CC2CC=CC=CC12(2251)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction51',
    reactants = ['H(19)', 'O=[C]CC1C=CC=CC=C1(3041)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(187.105,'m^3/(mol*s)'), n=1.3603, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 20 used for Cds-CsH_Cds-(Cd-Cd-Cd)H;HJ
Exact match found for rate rule [Cds-CsH_Cds-(Cd-Cd-Cd)H;HJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from -1.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction52',
    reactants = ['O=[C]CC1C=CC=[C]CC1(3047)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(1.9054e+11,'s^-1'), n=0.853, Ea=(200.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction53',
    reactants = ['O=[C]C[C]1C=CC=CCC1(3043)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(1.28e+07,'s^-1'), n=1.56, Ea=(126.775,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_noH;Cs_H_out_H/Cd] for rate rule [R3H_SS_Cs;C_rad_out_OneDe/Cs;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction54',
    reactants = ['O=[C]CC1C=C[C]=CCC1(3051)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_single;Cs_H_out_H/NonDeC] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction55',
    reactants = ['O=[C][CH]C1C=CC=CCC1(3045)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(3.56e+06,'s^-1'), n=1.75, Ea=(105.855,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;Y_rad_out;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction56',
    reactants = ['O=[C]CC1[C]=CC=CCC1(3046)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out_H/(Cd-Cd-Cd)] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction57',
    reactants = ['O=[C]CC1C=[C]C=CCC1(3050)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(2.10713e+09,'s^-1'), n=0.59575, Ea=(261.003,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SDS;Y_rad_out;Cs_H_out_H/(NonDeC/Cs)] + [R4H_SDS;Cd_rad_out_Cd;Cs_H_out] for rate rule [R4H_SDS;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction58',
    reactants = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    products = ['O=[C]CC1[CH]C2C=CC2C1(3141)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(5.65845e+18,'s^-1'), n=-1.57151, Ea=(227.049,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 131 used for Rn0c7_gamma_long_SDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCd
Exact match found for rate rule [Rn0c7_gamma_long_SDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCd]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction59',
    reactants = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    products = ['O=CCC1C=CC=CC=C1(3068)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(5.18e+08,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_De] for rate rule [R4radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction60',
    reactants = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    products = ['O=C=CC1C=CC=CCC1(3067)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction61',
    reactants = ['O=[C]CC1=CC=CC[CH]C1(3142)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(3.94565e+09,'s^-1'), n=0.909333, Ea=(116.555,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [1_3_unsaturated_pentane_backbone;CH(CJ)_1;unsaturated_end] + [1_3_pentadiene;CH_end;unsaturated_end] for rate rule [1_3_pentadiene;CH(CJ)_1;Cd(C)C_2]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction62',
    reactants = ['[CH2]C1C=CC=CC1C[C]=O(3071)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction63',
    reactants = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    products = ['O=C1CC2C=CC=CC1C2(3143)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(1.24579e+11,'s^-1'), n=0.1555, Ea=(7.322,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R5_SSSS;C_rad_out_single;Ypri_rad_out] for rate rule [R5_SSSS;C_rad_out_H/OneDe;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

network(
    label = '498',
    isomers = [
        'O=[C]CC1[CH]C=CC=CC1(2245)',
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
    label = '498',
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

