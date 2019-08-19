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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.570022,0.0938456,-7.73464e-05,3.29559e-08,-5.69593e-12,34931.1,22.2447], Tmin=(100,'K'), Tmax=(1367.25,'K')), NASAPolynomial(coeffs=[18.6404,0.037644,-1.56881e-05,2.89157e-09,-1.9871e-13,29678,-76.4436], Tmin=(1367.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(289.021,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3-Cycloheptadiene) + radical(CCCJ=O) + radical(Cs_S)"""),
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
    label = 'O=[C]CC1C=CC=CCC=1(3040)',
    structure = SMILES('O=[C]CC1C=CC=CCC=1'),
    E0 = (188.611,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2750,2850,1437.5,1250,1305,750,350,1855,455,950,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (133.167,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.684062,0.0830321,-2.95454e-05,-2.84837e-08,1.80634e-11,22870.9,29.7145], Tmin=(100,'K'), Tmax=(1010.33,'K')), NASAPolynomial(coeffs=[24.7192,0.0263936,-1.06859e-05,2.11342e-09,-1.58032e-13,15495.3,-104.2], Tmin=(1010.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.611,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3,5-Cycloheptatriene) + radical(CCCJ=O)"""),
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
    label = 'O=C=CC1[CH]CC=CC=C1(3042)',
    structure = SMILES('O=C=CC1[CH]CC=CC=C1'),
    E0 = (252.104,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,2120,512.5,787.5,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (133.167,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.703783,0.0958474,-8.45192e-05,3.8122e-08,-6.88945e-12,30496.6,19.0714], Tmin=(100,'K'), Tmax=(1324.31,'K')), NASAPolynomial(coeffs=[19.9435,0.033484,-1.38828e-05,2.5634e-09,-1.76824e-13,25027.9,-86.3394], Tmin=(1324.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(252.104,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCd(CCO)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-(Cdd-O2d)CsH) + ring(1,3-Cycloheptadiene) + radical(Cs_S)"""),
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
    label = 'O=[C]CC1[CH]C=CC=CC1(2245)',
    structure = SMILES('O=[C]CC1C=C[CH]C=CC1'),
    E0 = (213.702,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.819605,0.00725946,0.000251568,-3.74444e-07,1.57485e-10,25875.7,34.5197], Tmin=(100,'K'), Tmax=(906.275,'K')), NASAPolynomial(coeffs=[43.2142,-0.0160425,1.90032e-05,-3.91762e-09,2.54642e-13,11464.2,-202.952], Tmin=(906.275,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(213.702,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(1,4-Cycloheptadiene) + radical(CCCJ=O) + radical(C=CCJC=C)"""),
)

species(
    label = 'O=C[CH]C1[CH]CC=CC=C1(3044)',
    structure = SMILES('[O]C=CC1[CH]CC=CC=C1'),
    E0 = (271.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.990733,0.094221,-5.92769e-05,3.40592e-09,6.8676e-12,32879,20.2589], Tmin=(100,'K'), Tmax=(1018.35,'K')), NASAPolynomial(coeffs=[22.9374,0.0311488,-1.19112e-05,2.20947e-09,-1.57309e-13,26402.5,-103.486], Tmin=(1018.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(271.764,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(1,3-Cycloheptadiene) + radical(C=COJ) + radical(Cs_S)"""),
)

species(
    label = 'O=[C][CH]C1C=CC=CCC1(3045)',
    structure = SMILES('O=[C][CH]C1C=CC=CCC1'),
    E0 = (262.009,'kJ/mol'),
    spinMultiplicity = 1,
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
    label = 'O=[C]CC1C=CC=[C]CC1(3047)',
    structure = SMILES('O=[C]CC1C=CC=[C]CC1'),
    E0 = (332.321,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.724695,0.0988952,-8.83232e-05,4.13615e-08,-7.84243e-12,40143.2,20.7151], Tmin=(100,'K'), Tmax=(1258.12,'K')), NASAPolynomial(coeffs=[18.4231,0.0380176,-1.57414e-05,2.90102e-09,-1.99956e-13,35325.2,-76.0586], Tmin=(1258.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(332.321,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3-Cycloheptadiene) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = 'O=CC[C]1[CH]CC=CC=C1(3048)',
    structure = SMILES('O=CCC1[CH]CC=C[CH]C=1'),
    E0 = (194.868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16241,-0.00762165,0.0003009,-4.23198e-07,1.72466e-10,23604.6,31.0344], Tmin=(100,'K'), Tmax=(919.774,'K')), NASAPolynomial(coeffs=[44.231,-0.0140479,1.64035e-05,-3.18653e-09,1.91195e-13,8031.06,-214.736], Tmin=(919.774,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(194.868,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(1,4-Cycloheptadiene) + radical(Allyl_S) + radical(C=CCJC=C)"""),
)

species(
    label = 'O=CCC1[C]=CC=CC[CH]1(3049)',
    structure = SMILES('O=CCC1[C]=CC=CC[CH]1'),
    E0 = (366.903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.194408,0.0916056,-7.39112e-05,3.13401e-08,-5.52231e-12,44279.4,20.5462], Tmin=(100,'K'), Tmax=(1313.49,'K')), NASAPolynomial(coeffs=[15.4817,0.0438669,-1.93938e-05,3.66972e-09,-2.55745e-13,40161.3,-59.3568], Tmin=(1313.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(366.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3-Cycloheptadiene) + radical(Cds_S) + radical(Cs_S)"""),
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
    label = 'O=CCC1[CH]CC=C[C]=C1(3053)',
    structure = SMILES('O=CCC1[CH]CC=C[C]=C1'),
    E0 = (328.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.261223,0.0937816,-7.89956e-05,3.60467e-08,-6.87755e-12,39609.3,19.9118], Tmin=(100,'K'), Tmax=(1221.38,'K')), NASAPolynomial(coeffs=[14.4713,0.0455328,-1.97402e-05,3.70328e-09,-2.57267e-13,36010.5,-54.1104], Tmin=(1221.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(328.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3-Cycloheptadiene) + radical(Cs_S) + radical(C=CJC=C)"""),
)

species(
    label = 'O=CCC1[CH]C[C]=CC=C1(3054)',
    structure = SMILES('O=CCC1[CH]C[C]=CC=C1'),
    E0 = (366.903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.194408,0.0916056,-7.39112e-05,3.13401e-08,-5.52231e-12,44279.4,20.5462], Tmin=(100,'K'), Tmax=(1313.49,'K')), NASAPolynomial(coeffs=[15.4817,0.0438669,-1.93938e-05,3.66972e-09,-2.55745e-13,40161.3,-59.3568], Tmin=(1313.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(366.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3-Cycloheptadiene) + radical(Cs_S) + radical(Cds_S)"""),
)

species(
    label = 'O=CCC1[CH]CC=[C]C=C1(3055)',
    structure = SMILES('O=CCC1[CH]CC=[C]C=C1'),
    E0 = (328.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.261223,0.0937816,-7.89956e-05,3.60467e-08,-6.87755e-12,39609.3,19.9118], Tmin=(100,'K'), Tmax=(1221.38,'K')), NASAPolynomial(coeffs=[14.4713,0.0455328,-1.97402e-05,3.70328e-09,-2.57267e-13,36010.5,-54.1104], Tmin=(1221.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(328.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3-Cycloheptadiene) + radical(C=CJC=C) + radical(Cs_S)"""),
)

species(
    label = 'O=[C]CC1[CH]C2C=CCC21(3056)',
    structure = SMILES('O=[C]CC1[CH]C2C=CCC21'),
    E0 = (333.21,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.625367,0.0599528,2.98976e-06,-3.61356e-08,1.49041e-11,40209.5,31.5998], Tmin=(100,'K'), Tmax=(1089.46,'K')), NASAPolynomial(coeffs=[13.3172,0.0446053,-1.8907e-05,3.59311e-09,-2.54388e-13,35589.4,-39.2302], Tmin=(1089.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(333.21,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + polycyclic(s2_4_5_ene_1) + radical(cyclobutane) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C]CC1C=CC2[CH]CC21(3057)',
    structure = SMILES('O=[C]CC1C=CC2[CH]CC21'),
    E0 = (334.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.625208,0.0596878,4.36457e-06,-3.80699e-08,1.57128e-11,40338.4,31.6341], Tmin=(100,'K'), Tmax=(1080.93,'K')), NASAPolynomial(coeffs=[13.445,0.0442706,-1.86784e-05,3.54876e-09,-2.51524e-13,35696.1,-39.865], Tmin=(1080.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(334.28,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + polycyclic(s2_4_5_ene_1) + radical(cyclobutane) + radical(CCCJ=O)"""),
)

species(
    label = 'O=C1CC2[CH]CC=CC1[CH]2(3058)',
    structure = SMILES('O=C1CC2[CH]CC=CC1[CH]2'),
    E0 = (202.017,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.102498,0.0916567,-7.43453e-05,3.28651e-08,-6.27969e-12,24432.2,1.80633], Tmin=(100,'K'), Tmax=(1182.45,'K')), NASAPolynomial(coeffs=[11.7213,0.0523525,-2.44858e-05,4.75431e-09,-3.36357e-13,21684.5,-56.195], Tmin=(1182.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(202.017,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclopentane) + ring(Cycloheptene) + ring(Cyclopentane) + radical(CCJCC=O) + radical(Cs_S)"""),
)

species(
    label = 'O=[C]CC1C=C[CH]C2CC21(3059)',
    structure = SMILES('O=[C]CC1[CH]C=CC2CC12'),
    E0 = (265.119,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.580061,0.0570729,2.17323e-05,-6.23022e-08,2.56821e-11,32025.4,27.8651], Tmin=(100,'K'), Tmax=(1022.18,'K')), NASAPolynomial(coeffs=[15.5257,0.0408809,-1.657e-05,3.15607e-09,-2.27055e-13,26760.4,-55.3749], Tmin=(1022.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(265.119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + polycyclic(s2_3_6_ene_1) + radical(cyclohexene-allyl) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C]CC1C2[CH]C=CCC21(3060)',
    structure = SMILES('O=[C]CC1C2C=C[CH]CC21'),
    E0 = (265.119,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.580061,0.0570729,2.17323e-05,-6.23022e-08,2.56821e-11,32025.4,27.8651], Tmin=(100,'K'), Tmax=(1022.18,'K')), NASAPolynomial(coeffs=[15.5257,0.0408809,-1.657e-05,3.15607e-09,-2.27055e-13,26760.4,-55.3749], Tmin=(1022.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(265.119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + polycyclic(s2_3_6_ene_1) + radical(CCCJ=O) + radical(cyclohexene-allyl)"""),
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
    label = 'O=C1CC2[CH]CC1[CH]C=C2(3062)',
    structure = SMILES('O=C1CC2[CH]C=CC1C[CH]2'),
    E0 = (143.228,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.179081,0.0908091,-7.11651e-05,3.02088e-08,-5.58391e-12,17358,-2.07354], Tmin=(100,'K'), Tmax=(1206.55,'K')), NASAPolynomial(coeffs=[11.2706,0.0540381,-2.54511e-05,4.94998e-09,-3.50247e-13,14681.5,-57.6663], Tmin=(1206.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(143.228,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclohexane) + ring(Cycloheptene) + ring(Cyclohexane) + radical(Cs_S) + radical(Allyl_S)"""),
)

species(
    label = 'O=C1CC2[CH]CC=C[CH]C12(3063)',
    structure = SMILES('O=C1CC2[CH]C[CH]C=CC12'),
    E0 = (143.228,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.179079,0.0908091,-7.11652e-05,3.02088e-08,-5.58393e-12,17358,-2.07353], Tmin=(100,'K'), Tmax=(1206.55,'K')), NASAPolynomial(coeffs=[11.2706,0.0540381,-2.54511e-05,4.94998e-09,-3.50248e-13,14681.5,-57.6663], Tmin=(1206.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(143.228,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclobutane) + ring(Cycloheptene) + ring(Cyclobutane) + radical(Allyl_S) + radical(Cs_S)"""),
)

species(
    label = '[CH]=CC=CC(C=C)C[C]=O(3064)',
    structure = SMILES('[CH]=CC=CC(C=C)C[C]=O'),
    E0 = (407.429,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3120,650,792.5,1650,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.637361,0.104819,-0.000111278,6.50621e-08,-1.55883e-11,49167.1,37.2009], Tmin=(100,'K'), Tmax=(1002.81,'K')), NASAPolynomial(coeffs=[14.8809,0.0429219,-1.86939e-05,3.51365e-09,-2.44611e-13,46054.6,-37.7094], Tmin=(1002.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(407.429,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]=CC=CCC=CC[C]=O(3065)',
    structure = SMILES('[CH]=CC=CCC=CC[C]=O'),
    E0 = (412.29,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,3120,650,792.5,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1855,455,950,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.703319,0.088721,-5.76854e-05,1.06046e-08,1.94435e-12,49768.8,39.8796], Tmin=(100,'K'), Tmax=(1135.85,'K')), NASAPolynomial(coeffs=[21.4664,0.0333055,-1.44252e-05,2.77558e-09,-1.97812e-13,43270.9,-76.3344], Tmin=(1135.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(412.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CCCJ=O)"""),
)

species(
    label = 'O=CCC1C=CC=CCC=1(3066)',
    structure = SMILES('O=CCC1C=CC=CCC=1'),
    E0 = (28.6503,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.625914,0.0802202,-1.49333e-05,-4.23548e-08,2.22009e-11,3631.24,28.6884], Tmin=(100,'K'), Tmax=(1017.96,'K')), NASAPolynomial(coeffs=[24.1304,0.0308809,-1.287e-05,2.5565e-09,-1.90491e-13,-3892.73,-103.387], Tmin=(1017.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.6503,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3,5-Cycloheptatriene)"""),
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
    label = 'O=[C]CC1[CH]CCC=CC=1(3069)',
    structure = SMILES('O=[C]CC1[CH]C=CCCC=1'),
    E0 = (213.716,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.901366,0.000438717,0.000277427,-4.01268e-07,1.657e-10,25878.8,33.8973], Tmin=(100,'K'), Tmax=(916.444,'K')), NASAPolynomial(coeffs=[44.9195,-0.0162834,1.77026e-05,-3.48508e-09,2.1507e-13,10445,-214.811], Tmin=(916.444,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(213.716,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(1,4-Cycloheptadiene) + radical(CCCJ=O) + radical(C=CCJC=C)"""),
)

species(
    label = 'O=[C]C[CH]C1C=CC=CC1(3070)',
    structure = SMILES('O=[C]C[CH]C1C=CC=CC1'),
    E0 = (286.123,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.207034,0.066462,-1.32103e-06,-4.3768e-08,2.07384e-11,34564.1,34.7209], Tmin=(100,'K'), Tmax=(1007.35,'K')), NASAPolynomial(coeffs=[17.4072,0.0366314,-1.41822e-05,2.65172e-09,-1.89669e-13,29147,-58.0737], Tmin=(1007.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(286.123,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(469.768,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3-Cyclohexadiene) + radical(CCJCC=O) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C1C=CC=CC1C[C]=O(3071)',
    structure = SMILES('[CH2]C1C=CC=CC1C[C]=O'),
    E0 = (283.332,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0921783,0.0699942,-8.95123e-06,-3.87872e-08,2.01368e-11,34231.8,33.4785], Tmin=(100,'K'), Tmax=(980.026,'K')), NASAPolynomial(coeffs=[17.5951,0.0355214,-1.2767e-05,2.2965e-09,-1.61581e-13,29026,-59.6668], Tmin=(980.026,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(283.332,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(469.768,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + ring(1,3-Cyclohexadiene) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = 'C=C([O])C1[CH]CC=CC=C1(2238)',
    structure = SMILES('C=C([O])C1[CH]CC=CC=C1'),
    E0 = (262.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.04339,0.10061,-8.92457e-05,4.06621e-08,-7.38355e-12,31742.4,20.7418], Tmin=(100,'K'), Tmax=(1326.26,'K')), NASAPolynomial(coeffs=[21.3094,0.0331933,-1.29969e-05,2.33409e-09,-1.58673e-13,25813.3,-93.409], Tmin=(1326.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(262.341,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(1,3-Cycloheptadiene) + radical(Cs_S) + radical(C=C(C)OJ)"""),
)

species(
    label = 'O=C1CC2C=CC=CCC12(2253)',
    structure = SMILES('O=C1CC2C=CC=CCC12'),
    E0 = (-94.381,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.79144,0.167054,-0.00034929,3.76319e-07,-1.47037e-10,-11181.1,-15.0288], Tmin=(100,'K'), Tmax=(854.275,'K')), NASAPolynomial(coeffs=[-7.28853,0.0951242,-5.14965e-05,1.00923e-08,-6.96947e-13,-6678.07,31.4847], Tmin=(854.275,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-94.381,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclobutane) + ring(1,3-Cycloheptadiene) + ring(Cyclobutane)"""),
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
    label = '[CH2]C1[CH]CC=CC=C1(582)',
    structure = SMILES('[CH2]C1[CH]CC=CC=C1'),
    E0 = (440.723,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.340624,0.0672658,-2.12794e-05,-2.10625e-08,1.31393e-11,53150.1,15.659], Tmin=(100,'K'), Tmax=(993.827,'K')), NASAPolynomial(coeffs=[16.4928,0.0315499,-1.15873e-05,2.09542e-09,-1.4708e-13,48492.9,-69.4445], Tmin=(993.827,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(440.723,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cycloheptadiene) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = 'O=[C]CC1[CH]CC2C=CC12(3072)',
    structure = SMILES('O=[C]CC1[CH]CC2C=CC12'),
    E0 = (368.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.175,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.677141,0.058181,8.01385e-06,-4.17217e-08,1.69567e-11,44492,30.784], Tmin=(100,'K'), Tmax=(1075.7,'K')), NASAPolynomial(coeffs=[13.4603,0.0439503,-1.85819e-05,3.54196e-09,-2.5176e-13,39814.9,-40.7766], Tmin=(1075.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(368.828,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + Estimated bicyclic component: polycyclic(s2_4_5_ane) - ring(Cyclopentane) - ring(Cyclobutane) + ring(Cyclopentane) + ring(Cyclobutene) + radical(CCCJ=O) + radical(Cs_S)"""),
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
    E0 = (289.021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (411.532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (418.264,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (474.888,'kJ/mol'),
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
    E0 = (318.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (406.18,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (388.789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (445.895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (433.062,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (473.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (466.841,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (405.157,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (468.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (337.783,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (337.783,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (349.921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (361.096,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (518.562,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (447.609,'kJ/mol'),
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
    E0 = (419.879,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (419.879,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (369.773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (290.277,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (290.277,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (385.372,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (316.217,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (338.513,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (448.474,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (439.486,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (352.422,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (352.422,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (313.995,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (334.974,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (448.957,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (534.622,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (534.622,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (297.306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (347.122,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (411.434,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (879.808,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    products = ['C7H8(11)', 'CH2CO(43)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(19)', 'O=[C]CC1C=CC=CCC=1(3040)'],
    products = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.11e+08,'cm^3/(mol*s)'), n=1.64, Ea=(11.1294,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2594 used for Cds-CdCs_Cds-CsH;HJ
Exact match found for rate rule [Cds-CdCs_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(19)', 'O=[C]CC1C=CC=CC=C1(3041)'],
    products = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(0.0545849,'m^3/(mol*s)'), n=2.81111, Ea=(21.1569,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 26 used for Cds-CdH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CdH_Cds-CsH;HJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(19)', 'O=C=CC1[CH]CC=CC=C1(3042)'],
    products = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.82e-16,'cm^3/(molecule*s)'), n=1.61, Ea=(10.992,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Ck;HJ] for rate rule [Cds-CsH_Ck;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C7H8(11)', '[CH2][C]=O(360)'],
    products = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.0147295,'m^3/(mol*s)'), n=2.37207, Ea=(16.7338,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CdH_Cds-CsH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]1C=C[CH]CC=C1(125)', 'CH2CO(43)'],
    products = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.568606,'m^3/(mol*s)'), n=1.93802, Ea=(45.6341,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Ck;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O=[C]C[C]1C=CC=CCC1(3043)'],
    products = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.26e+10,'s^-1'), n=0.77, Ea=(192.464,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 164 used for R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    products = ['O=[C]CC1[CH]C=CC=CC1(2245)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.2802e+09,'s^-1'), n=0.983917, Ea=(99.768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;C_rad_out_1H;Cs_H_out_H/(Cd-Cd-Cd)] + [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Cd] for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    products = ['O=C[CH]C1[CH]CC=CC=C1(3044)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(791180,'s^-1'), n=2.19286, Ea=(156.873,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    products = ['O=[C][CH]C1C=CC=CCC1(3045)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.17661e+06,'s^-1'), n=1.79367, Ea=(144.041,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['O=[C]CC1[C]=CC=CCC1(3046)'],
    products = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.32e+09,'s^-1'), n=0.99, Ea=(141.419,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 206 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    products = ['O=[C]CC1C=CC=[C]CC1(3047)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.37e+07,'s^-1'), n=1.41, Ea=(177.82,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 207 used for R3H_SS_Cs;C_rad_out_H/NonDeC;Cd_H_out_doubleC
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    products = ['O=CC[C]1[CH]CC=CC=C1(3048)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(285601,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;XH_out] for rate rule [R3H_SS_Cs;CO_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['O=CCC1[C]=CC=CC[CH]1(3049)'],
    products = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O=[C]CC1C=[C]C=CCC1(3050)'],
    products = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['O=[C]CC1C=C[C]=CCC1(3051)'],
    products = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 2.82842712475
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    products = ['O=CCC1[CH][CH]C=CC=C1(3052)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(61085.3,'s^-1'), n=1.8624, Ea=(60.8997,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R5HJ_3;CO_rad_out;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['O=CCC1[CH]CC=C[C]=C1(3053)'],
    products = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_single;XH_out] for rate rule [R5H_DSSS;Cd_rad_out_singleDe_Cd;CO_H_out]
Euclidian distance = 3.74165738677
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['O=CCC1[CH]C[C]=CC=C1(3054)'],
    products = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.81182e+08,'s^-1'), n=1.25566, Ea=(151.659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;XH_out] for rate rule [R6HJ_2;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['O=CCC1[CH]CC=[C]C=C1(3055)'],
    products = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.62236e+06,'s^-1'), n=1.59783, Ea=(119.552,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_Cd;XH_out] + [R6H_SMSSR;Y_rad_out;XH_out] for rate rule [R6H_SMSSR;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 3.16227766017
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]1C=C[CH]CC=C1(125)', '[CH2][C]=O(360)'],
    products = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.47003e+07,'m^3/(mol*s)'), n=0.0284742, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    products = ['O=[C]CC1[CH]C2C=CCC21(3056)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(130.857,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 116 used for Rn0c7_gamma_short_SSDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCs
Exact match found for rate rule [Rn0c7_gamma_short_SSDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    products = ['O=[C]CC1C=CC2[CH]CC21(3057)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(130.857,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 116 used for Rn0c7_gamma_short_SSDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCs
Exact match found for rate rule [Rn0c7_gamma_short_SSDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    products = ['O=C1CC2[CH]CC=CC1[CH]2(3058)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(9.9e+10,'s^-1'), n=0.06, Ea=(80.7512,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn2cx_beta;doublebond_intra_pri_HCd;radadd_intra] for rate rule [Rn2c7_beta_short;doublebond_intra_pri_HCd;radadd_intra_CO]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    products = ['O=[C]CC1C=C[CH]C2CC21(3059)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.1e+12,'s^-1'), n=0.14, Ea=(1.2552,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    products = ['O=[C]CC1C2[CH]C=CCC21(3060)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.1e+12,'s^-1'), n=0.14, Ea=(1.2552,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    products = ['O=C1CC2[CH]C[CH]C1C=C2(3061)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.54269e+17,'s^-1'), n=-0.500262, Ea=(96.3509,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_cyclic;doublebond_intra_pri_HCd;radadd_intra] for rate rule [Rn2c7_gamma_long;doublebond_intra_pri_HCd;radadd_intra_CO]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    products = ['O=C1CC2[CH]CC1[CH]C=C2(3062)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(110000,'s^-1'), n=1.18, Ea=(27.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn2c7_gamma_long;doublebond_intra_pri_HNd_Cs;radadd_intra_CO]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    products = ['O=C1CC2[CH]CC=C[CH]C12(3063)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.99832e+10,'s^-1'), n=0.37247, Ea=(49.4918,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri;radadd_intra] for rate rule [Rn2c7_alpha_long;doublebond_intra_pri;radadd_intra_CO]
Euclidian distance = 3.16227766017
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]=CC=CC(C=C)C[C]=O(3064)'],
    products = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(4.625e+11,'s^-1'), n=0.16, Ea=(41.045,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;doublebond_intra_pri;radadd_intra_cdsingleH] for rate rule [R7_linear;doublebond_intra_pri_2H;radadd_intra_cdsingleH]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=CC=CCC=CC[C]=O(3065)'],
    products = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(110000,'s^-1'), n=1.18, Ea=(27.196,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R7_linear;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [R7_linear;doublebond_intra_pri_HNd_Cs;radadd_intra_cdsingleH]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction32',
    reactants = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    products = ['O=CCC1C=CC=CCC=1(3066)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    products = ['O=C=CC1C=CC=CCC1(3067)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    products = ['O=CCC1C=CC=CC=C1(3068)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction35',
    reactants = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    products = ['O=[C]CC1[CH]CCC=CC=1(3069)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2.125e+09,'s^-1'), n=0.991, Ea=(45.9529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_3_pentadiene;CH_end;CdHC_2] for rate rule [1_3_pentadiene;CH(CJ)_1;CdHC_2]
Euclidian distance = 1.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction36',
    reactants = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    products = ['O=[C]C[CH]C1C=CC=CC1(3070)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction37',
    reactants = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    products = ['[CH2]C1C=CC=CC1C[C]=O(3071)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HH)CJ;CsJ;C] for rate rule [cCs(-HH)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction38',
    reactants = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    products = ['C=C([O])C1[CH]CC=CC=C1(2238)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction39',
    reactants = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    products = ['O=C1CC2C=CC=CCC12(2253)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['CO(28)', '[CH2]C1[CH]CC=CC=C1(582)'],
    products = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1230.59,'m^3/(mol*s)'), n=1.0523, Ea=(25.6182,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [COm;C_rad/H2/Cs]
Euclidian distance = 0
family: R_Addition_COm"""),
)

reaction(
    label = 'reaction41',
    reactants = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    products = ['O=[C]CC1[CH]CC2C=CC12(3072)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1,3-butadiene_backbone;C=C_1;C=C_2] for rate rule [1,3-butadiene_backbone;CdH(C)_1;CdH(C)_2]
Euclidian distance = 1.41421356237
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[C]=O(353)', '[CH2]C1[CH]CC=CC=C1(582)'],
    products = ['O=[C]CC1[CH]CC=CC=C1(2244)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(2.23625e+06,'m^3/(mol*s)'), n=0.36814, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

network(
    label = '497',
    isomers = [
        'O=[C]CC1[CH]CC=CC=C1(2244)',
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
    label = '497',
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

