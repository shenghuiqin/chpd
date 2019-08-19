species(
    label = 'C=CC=CC(CC=O)O[O](1059)',
    structure = SMILES('C=CC=CC(CC=O)O[O]'),
    E0 = (-20.6265,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,492.5,1135,1000,2782.5,750,1395,475,1775,1000,2750,2850,1437.5,1250,1305,750,350,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4499.23,'J/mol'), sigma=(7.17214,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=702.77 K, Pc=27.67 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.343321,0.0952679,-8.99694e-05,4.51117e-08,-9.2483e-12,-2324.2,37.3712], Tmin=(100,'K'), Tmax=(1160.44,'K')), NASAPolynomial(coeffs=[16.0169,0.0388748,-1.70748e-05,3.23403e-09,-2.26358e-13,-6121.19,-43.9917], Tmin=(1160.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-20.6265,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(ROOJ)"""),
)

species(
    label = 'O2(2)',
    structure = SMILES('[O][O]'),
    E0 = (-8.62178,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1483.7],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(887.157,'J/mol'), sigma=(3.467,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53764,-0.00122828,5.36759e-06,-4.93128e-09,1.45955e-12,-1037.99,4.6718], Tmin=(100,'K'), Tmax=(1087.71,'K')), NASAPolynomial(coeffs=[3.16427,0.00169454,-8.00335e-07,1.5903e-10,-1.14891e-14,-1048.45,6.08303], Tmin=(1087.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.62178,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = '[CH2]C=CC=CCC=O(956)',
    structure = SMILES('[CH2]C=CC=CCC=O'),
    E0 = (40.6633,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,286.105,645.958],'cm^-1')),
        HinderedRotor(inertia=(0.0787087,'amu*angstrom^2'), symmetry=1, barrier=(22.8412,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.994264,'amu*angstrom^2'), symmetry=1, barrier=(22.8601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0789567,'amu*angstrom^2'), symmetry=1, barrier=(22.8332,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.994741,'amu*angstrom^2'), symmetry=1, barrier=(22.8711,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3882.24,'J/mol'), sigma=(6.33784,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=606.40 K, Pc=34.6 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.868578,0.0528089,9.16168e-06,-4.58481e-08,1.96156e-11,5017.29,28.9089], Tmin=(100,'K'), Tmax=(1039.72,'K')), NASAPolynomial(coeffs=[15.7717,0.0320753,-1.37307e-05,2.68843e-09,-1.96067e-13,-60.0845,-53.0844], Tmin=(1039.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(40.6633,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=CC=C[CH]O[O](670)',
    structure = SMILES('[CH2]C=CC=CO[O]'),
    E0 = (243.658,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,492.5,1135,1000,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.784965,'amu*angstrom^2'), symmetry=1, barrier=(18.0479,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.785044,'amu*angstrom^2'), symmetry=1, barrier=(18.0497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.93487,'amu*angstrom^2'), symmetry=1, barrier=(67.4784,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19261,0.0517996,-2.08373e-05,-1.61991e-08,1.20077e-11,29415.5,26.0962], Tmin=(100,'K'), Tmax=(937.741,'K')), NASAPolynomial(coeffs=[15.1311,0.0171574,-5.11526e-06,8.4127e-10,-5.79576e-14,25710.4,-46.0702], Tmin=(937.741,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.658,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(ROOJ) + radical(C=CC=CCJ)"""),
)

species(
    label = 'CH2CHO(51)',
    structure = SMILES('[CH2]C=O'),
    E0 = (1.22925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,526.75,532.597,975.94,1639.13,1641.45],'cm^-1')),
        HinderedRotor(inertia=(0.00114821,'amu*angstrom^2'), symmetry=1, barrier=(2.1986,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66874,0.0096233,1.60617e-05,-2.87682e-08,1.2503e-11,219.438,12.5694], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.91637,0.0088465,-3.14955e-06,5.05413e-10,-3.01305e-14,-1047.8,-6.1065], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(1.22925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CH2CHO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'CH2CHCHCH(170)',
    structure = SMILES('[CH]=CC=C'),
    E0 = (346.45,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(1.31937,'amu*angstrom^2'), symmetry=1, barrier=(30.3349,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (53.0825,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.64255,0.0163337,3.86225e-05,-6.71377e-08,2.83603e-11,41729.6,13.282], Tmin=(100,'K'), Tmax=(937.724,'K')), NASAPolynomial(coeffs=[12.9705,0.00669127,-1.0007e-06,1.67602e-10,-1.71436e-14,38279.7,-43.9476], Tmin=(937.724,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(346.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""CH2CHCHCH""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[O]O[CH]CC=O(1096)',
    structure = SMILES('[O]O[CH]CC=O'),
    E0 = (59.2375,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,3025,407.5,1350,352.5,492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,180],'cm^-1')),
        HinderedRotor(inertia=(0.241917,'amu*angstrom^2'), symmetry=1, barrier=(5.56215,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.241711,'amu*angstrom^2'), symmetry=1, barrier=(5.55741,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.243716,'amu*angstrom^2'), symmetry=1, barrier=(5.60351,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0621,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38092,0.0716649,-0.000143172,1.43316e-07,-5.23474e-11,7205.57,21.5986], Tmin=(100,'K'), Tmax=(883.04,'K')), NASAPolynomial(coeffs=[2.98182,0.0303115,-1.49979e-05,2.81594e-09,-1.88403e-13,8252.39,21.6026], Tmin=(883.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(59.2375,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(CCsJOOH)"""),
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
    label = 'C=CC=C[C](CC=O)O[O](1097)',
    structure = SMILES('[CH2]C=CC=C(CC=O)O[O]'),
    E0 = (102.604,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,492.5,1135,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.364475,'amu*angstrom^2'), symmetry=1, barrier=(8.37999,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0147916,'amu*angstrom^2'), symmetry=1, barrier=(19.2269,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0147916,'amu*angstrom^2'), symmetry=1, barrier=(19.2272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.836252,'amu*angstrom^2'), symmetry=1, barrier=(19.2271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.836256,'amu*angstrom^2'), symmetry=1, barrier=(19.2272,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.085428,0.0822534,-6.29043e-05,2.36171e-08,-3.56029e-12,12493.4,37.1999], Tmin=(100,'K'), Tmax=(1551.9,'K')), NASAPolynomial(coeffs=[19.5045,0.0317604,-1.40999e-05,2.65156e-09,-1.82891e-13,6413.14,-65.9195], Tmin=(1551.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(102.604,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + radical(C=CC=CCJ) + radical(ROOJ)"""),
)

species(
    label = 'HCO(30)',
    structure = SMILES('[CH]=O'),
    E0 = (32.4782,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1131.19,1955.83,1955.83],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.23755,-0.00332075,1.4003e-05,-1.3424e-08,4.37416e-12,3872.41,3.30835], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.92002,0.00252279,-6.71004e-07,1.05616e-10,-7.43798e-15,3653.43,3.58077], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(32.4782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HCO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH2]C(C=CC=C)O[O](1098)',
    structure = SMILES('[CH2]C(C=CC=C)O[O]'),
    E0 = (299.916,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,492.5,1135,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,227.499,227.549],'cm^-1')),
        HinderedRotor(inertia=(0.407757,'amu*angstrom^2'), symmetry=1, barrier=(14.9797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.407992,'amu*angstrom^2'), symmetry=1, barrier=(14.9799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.407714,'amu*angstrom^2'), symmetry=1, barrier=(14.9797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.407677,'amu*angstrom^2'), symmetry=1, barrier=(14.9791,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.077792,0.0765788,-6.29129e-05,2.20909e-08,-1.69906e-12,36221.1,33.5752], Tmin=(100,'K'), Tmax=(1050.03,'K')), NASAPolynomial(coeffs=[18.0054,0.0235206,-8.88194e-06,1.60453e-09,-1.11389e-13,31616.3,-57.7897], Tmin=(1050.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(299.916,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(ROOJ) + radical(CJCOOH)"""),
)

species(
    label = 'C=CC=CC([CH]C=O)O[O](1099)',
    structure = SMILES('C=CC=CC(C=C[O])O[O]'),
    E0 = (118.868,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,492.5,1135,1000,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.89987,'amu*angstrom^2'), symmetry=1, barrier=(20.6898,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.90131,'amu*angstrom^2'), symmetry=1, barrier=(20.7229,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.899561,'amu*angstrom^2'), symmetry=1, barrier=(20.6827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.900123,'amu*angstrom^2'), symmetry=1, barrier=(20.6956,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.85555,0.0923745,-6.78579e-05,1.08125e-08,5.33863e-12,14484.2,38.7541], Tmin=(100,'K'), Tmax=(990.317,'K')), NASAPolynomial(coeffs=[23.9179,0.022616,-8.09777e-06,1.48261e-09,-1.06541e-13,8091.44,-88.0253], Tmin=(990.317,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(118.868,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(ROOJ)"""),
)

species(
    label = 'C=CC=[C]C(CC=O)O[O](1100)',
    structure = SMILES('C=CC=[C]C(CC=O)O[O]'),
    E0 = (217.215,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,492.5,1135,1000,1685,370,2782.5,750,1395,475,1775,1000,2750,2850,1437.5,1250,1305,750,350,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.270749,0.0984534,-0.000109134,6.63616e-08,-1.65693e-11,26275,37.5972], Tmin=(100,'K'), Tmax=(962.213,'K')), NASAPolynomial(coeffs=[13.7167,0.0403052,-1.84841e-05,3.55369e-09,-2.50348e-13,23583.3,-29.3448], Tmin=(962.213,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(217.215,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(ROOJ)"""),
)

species(
    label = 'C2H3(44)',
    structure = SMILES('[CH]=C'),
    E0 = (286.361,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,677.08,1086.68,3788.01],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (27.0452,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36378,0.000265766,2.79621e-05,-3.72987e-08,1.5159e-11,34475,7.9151], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.15027,0.00754021,-2.62998e-06,4.15974e-10,-2.45408e-14,33856.6,1.72812], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(286.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""C2H3""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH]=CC(CC=O)O[O](1101)',
    structure = SMILES('[CH]=CC(CC=O)O[O]'),
    E0 = (174.698,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,492.5,1135,1000,3120,650,792.5,1650,2782.5,750,1395,475,1775,1000,2750,2850,1437.5,1250,1305,750,350,287.915],'cm^-1')),
        HinderedRotor(inertia=(0.188866,'amu*angstrom^2'), symmetry=1, barrier=(11.11,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164504,'amu*angstrom^2'), symmetry=1, barrier=(9.67744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164522,'amu*angstrom^2'), symmetry=1, barrier=(9.67743,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164511,'amu*angstrom^2'), symmetry=1, barrier=(9.67746,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.856155,0.0743939,-8.83775e-05,6.03952e-08,-1.72342e-11,21119.9,30.8468], Tmin=(100,'K'), Tmax=(841.53,'K')), NASAPolynomial(coeffs=[9.38819,0.0338387,-1.60886e-05,3.12687e-09,-2.20935e-13,19683.9,-8.84323], Tmin=(841.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(174.698,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(ROOJ)"""),
)

species(
    label = 'C=C[C]=CC(CC=O)O[O](1102)',
    structure = SMILES('C=C[C]=CC(CC=O)O[O]'),
    E0 = (178.369,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,492.5,1135,1000,1685,370,2782.5,750,1395,475,1775,1000,2750,2850,1437.5,1250,1305,750,350,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.427055,0.101725,-0.000118155,7.61891e-08,-2.00746e-11,21608.5,37.2799], Tmin=(100,'K'), Tmax=(917.406,'K')), NASAPolynomial(coeffs=[13.6249,0.0404561,-1.79756e-05,3.38846e-09,-2.35577e-13,19030.3,-29.3013], Tmin=(917.406,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(178.369,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(ROOJ)"""),
)

species(
    label = 'C=[C]C=CC(CC=O)O[O](1103)',
    structure = SMILES('C=C=C[CH]C(CC=O)O[O]'),
    E0 = (161.653,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,540,610,2055,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,492.5,1135,1000,2782.5,750,1395,475,1775,1000,2750,2850,1437.5,1250,1305,750,350,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.547054,0.107461,-0.000133015,9.35398e-08,-2.72277e-11,19599.5,34.8798], Tmin=(100,'K'), Tmax=(828.769,'K')), NASAPolynomial(coeffs=[12.1536,0.0461596,-2.20594e-05,4.28313e-09,-3.02065e-13,17494.4,-24.0078], Tmin=(828.769,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(161.653,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJCO) + radical(ROOJ)"""),
)

species(
    label = 'C=CC=CC(C[C]=O)O[O](1104)',
    structure = SMILES('C=CC=CC(C[C]=O)O[O]'),
    E0 = (139.334,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.437617,0.0984375,-0.000105567,5.99948e-08,-1.37435e-11,16917.2,38.5323], Tmin=(100,'K'), Tmax=(1054.59,'K')), NASAPolynomial(coeffs=[16.3141,0.034899,-1.51926e-05,2.86358e-09,-2.00006e-13,13384,-43.1752], Tmin=(1054.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(139.334,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]=CC=CC(CC=O)O[O](1105)',
    structure = SMILES('[CH]=CC=CC(CC=O)O[O]'),
    E0 = (226.47,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1380,1390,370,380,2900,435,492.5,1135,1000,3120,650,792.5,1650,2782.5,750,1395,475,1775,1000,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.587562,'amu*angstrom^2'), symmetry=1, barrier=(13.5092,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.587529,'amu*angstrom^2'), symmetry=1, barrier=(13.5085,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.587617,'amu*angstrom^2'), symmetry=1, barrier=(13.5105,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.587401,'amu*angstrom^2'), symmetry=1, barrier=(13.5055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.587416,'amu*angstrom^2'), symmetry=1, barrier=(13.5058,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.364222,0.0980875,-0.000103869,5.86342e-08,-1.34482e-11,27393.5,37.8557], Tmin=(100,'K'), Tmax=(1048.17,'K')), NASAPolynomial(coeffs=[15.5969,0.0371766,-1.67007e-05,3.19219e-09,-2.24573e-13,24047.5,-39.8983], Tmin=(1048.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(226.47,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(ROOJ)"""),
)

species(
    label = 'CC=CC=C(CC=O)O[O](1108)',
    structure = SMILES('CC=CC=C(CC=O)O[O]'),
    E0 = (-15.4515,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2750,2850,1437.5,1250,1305,750,350,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,350,440,435,1725,291.427,291.427],'cm^-1')),
        HinderedRotor(inertia=(0.182369,'amu*angstrom^2'), symmetry=1, barrier=(10.9907,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182357,'amu*angstrom^2'), symmetry=1, barrier=(10.9907,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182354,'amu*angstrom^2'), symmetry=1, barrier=(10.9907,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182354,'amu*angstrom^2'), symmetry=1, barrier=(10.9907,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182361,'amu*angstrom^2'), symmetry=1, barrier=(10.9907,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.136792,0.085488,-6.74415e-05,2.6695e-08,-4.35052e-12,-1720.37,34.8272], Tmin=(100,'K'), Tmax=(1407.64,'K')), NASAPolynomial(coeffs=[16.2087,0.0398174,-1.87743e-05,3.64591e-09,-2.56951e-13,-6245.07,-48.2058], Tmin=(1407.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-15.4515,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + radical(ROOJ)"""),
)

species(
    label = '[O]OC1CC2OCC=CC12(1109)',
    structure = SMILES('[O]OC1CC2OCC=CC12'),
    E0 = (34.5564,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.500102,0.0638391,-1.12414e-05,-2.16879e-08,9.97511e-12,4293.08,26.4064], Tmin=(100,'K'), Tmax=(1129.31,'K')), NASAPolynomial(coeffs=[14.1082,0.0424112,-1.83396e-05,3.49441e-09,-2.46701e-13,-487.639,-48.4583], Tmin=(1129.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(34.5564,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_6_ene_1) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(CC=O)C1C=CC1(1113)',
    structure = SMILES('[O]OC(CC=O)C1C=CC1'),
    E0 = (33.3948,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.18931,0.0835052,-6.24454e-05,2.37492e-08,-3.66292e-12,4174.19,36.0545], Tmin=(100,'K'), Tmax=(1515.88,'K')), NASAPolynomial(coeffs=[18.2196,0.03493,-1.43801e-05,2.61112e-09,-1.76883e-13,-1407.08,-60.4165], Tmin=(1515.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(33.3948,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(440.667,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclobutene) + radical(ROOJ)"""),
)

species(
    label = 'C=C[CH]C1OOC1CC=O(1086)',
    structure = SMILES('C=C[CH]C1OOC1CC=O'),
    E0 = (2.00137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0472272,0.0899857,-7.41351e-05,3.24544e-08,-5.99959e-12,380.217,28.8126], Tmin=(100,'K'), Tmax=(1240.26,'K')), NASAPolynomial(coeffs=[13.555,0.0464217,-2.14483e-05,4.13437e-09,-2.91157e-13,-2970.45,-39.2633], Tmin=(1240.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.00137,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(440.667,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(12dioxetane) + radical(C=CCJCO)"""),
)

species(
    label = 'C=CC=CC1CC([O])OO1(1087)',
    structure = SMILES('C=CC=CC1CC([O])OO1'),
    E0 = (20.0681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0382947,0.0676252,7.32452e-06,-6.69867e-08,3.37842e-11,2574.25,33.9678], Tmin=(100,'K'), Tmax=(929.079,'K')), NASAPolynomial(coeffs=[21.1006,0.0267726,-7.16578e-06,1.13606e-09,-7.94021e-14,-3489.99,-77.6696], Tmin=(929.079,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.0681,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(12dioxolane) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C1C=CC(CC=O)OO1(1088)',
    structure = SMILES('[CH2]C1C=CC(CC=O)OO1'),
    E0 = (-8.79098,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0242605,0.0667856,7.80344e-06,-6.18301e-08,2.87956e-11,-893.164,33.9787], Tmin=(100,'K'), Tmax=(995.609,'K')), NASAPolynomial(coeffs=[22.1493,0.0285524,-1.12083e-05,2.20188e-09,-1.65015e-13,-7828.72,-85.5554], Tmin=(995.609,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.79098,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(36dihydro12dioxin) + radical(CJCOOH)"""),
)

species(
    label = 'C=CC=C[C](CC=O)OO(1089)',
    structure = SMILES('[CH2]C=CC=C(CC=O)OO'),
    E0 = (-49.4004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.627973,0.0888336,-6.73738e-05,2.45753e-08,-3.54958e-12,-5764.21,38.4092], Tmin=(100,'K'), Tmax=(1636.52,'K')), NASAPolynomial(coeffs=[23.642,0.0295127,-1.30016e-05,2.42578e-09,-1.6595e-13,-13707.9,-90.6339], Tmin=(1636.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-49.4004,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=CC=CC([CH]C=O)OO(1090)',
    structure = SMILES('C=CC=CC(C=C[O])OO'),
    E0 = (-33.1365,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.14274,0.0962261,-6.39791e-05,2.43351e-09,8.77247e-12,-3785.27,39.0267], Tmin=(100,'K'), Tmax=(995.267,'K')), NASAPolynomial(coeffs=[25.8047,0.023767,-8.79413e-06,1.65355e-09,-1.20869e-13,-10924.5,-99.7698], Tmin=(995.267,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-33.1365,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ)"""),
)

species(
    label = 'C=CC=[C]C(CC=O)OO(1091)',
    structure = SMILES('C=CC=[C]C(CC=O)OO'),
    E0 = (65.2106,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,1685,370,2782.5,750,1395,475,1775,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.493787,0.10164,-0.000103311,5.58893e-08,-1.23812e-11,8002.57,37.6331], Tmin=(100,'K'), Tmax=(1077.13,'K')), NASAPolynomial(coeffs=[15.8233,0.0410454,-1.89273e-05,3.6619e-09,-2.59285e-13,4487.44,-42.2999], Tmin=(1077.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(65.2106,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C=CC=CC(C[C]=O)OO(1092)',
    structure = SMILES('C=CC=CC(C[C]=O)OO'),
    E0 = (-12.6706,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.748915,0.102629,-0.000103079,5.35957e-08,-1.11887e-11,-1351.35,38.887], Tmin=(100,'K'), Tmax=(1153.44,'K')), NASAPolynomial(coeffs=[18.693,0.0352065,-1.53993e-05,2.91821e-09,-2.04645e-13,-5836.37,-57.6844], Tmin=(1153.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-12.6706,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O)"""),
)

species(
    label = 'C=C[C]=CC(CC=O)OO(1093)',
    structure = SMILES('C=C[C]=CC(CC=O)OO'),
    E0 = (26.3643,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,1685,370,2782.5,750,1395,475,1775,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.619167,0.104507,-0.000110752,6.3488e-08,-1.4877e-11,3334.87,37.2079], Tmin=(100,'K'), Tmax=(1024.29,'K')), NASAPolynomial(coeffs=[15.4822,0.0416288,-1.86717e-05,3.55707e-09,-2.49564e-13,36.3679,-40.8583], Tmin=(1024.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(26.3643,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = 'C=[C]C=CC(CC=O)OO(1094)',
    structure = SMILES('C=C=C[CH]C(CC=O)OO'),
    E0 = (9.64833,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,540,610,2055,3010,987.5,1337.5,450,1655,3615,1310,387.5,850,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.630782,0.108884,-0.000120492,7.36661e-08,-1.8723e-11,1321.28,34.4246], Tmin=(100,'K'), Tmax=(939.363,'K')), NASAPolynomial(coeffs=[13.6958,0.0478776,-2.30739e-05,4.52753e-09,-3.22375e-13,-1370.26,-33.7965], Tmin=(939.363,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(9.64833,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJCO)"""),
)

species(
    label = '[CH]=CC=CC(CC=O)OO(1095)',
    structure = SMILES('[CH]=CC=CC(CC=O)OO'),
    E0 = (74.465,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3120,650,792.5,1650,2782.5,750,1395,475,1775,1000,200,800],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.663758,0.102159,-0.000101036,5.18736e-08,-1.07677e-11,9124.37,38.1668], Tmin=(100,'K'), Tmax=(1153.7,'K')), NASAPolynomial(coeffs=[18.007,0.0374246,-1.68702e-05,3.23753e-09,-2.28409e-13,4816.33,-54.578], Tmin=(1153.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.465,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C=CC1[CH]C(CC=O)OO1(1079)',
    structure = SMILES('C=CC1[CH]C(CC=O)OO1'),
    E0 = (0.234255,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0278075,0.0737231,-1.84093e-05,-3.4976e-08,2.11808e-11,186.821,36.03], Tmin=(100,'K'), Tmax=(933.557,'K')), NASAPolynomial(coeffs=[18.7669,0.0306442,-9.36577e-06,1.53678e-09,-1.04356e-13,-4954.31,-62.0917], Tmin=(933.557,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(0.234255,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(12dioxolane) + radical(CCJCOOH)"""),
)

species(
    label = 'C=CC=CC1C[CH]OOO1(1106)',
    structure = SMILES('C=CC=CC1C[CH]OOO1'),
    E0 = (210.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.466987,0.077868,-7.35774e-06,-6.32026e-08,3.57894e-11,25537.4,30.641], Tmin=(100,'K'), Tmax=(907.347,'K')), NASAPolynomial(coeffs=[24.8783,0.0208187,-3.44815e-06,3.47604e-10,-2.18628e-14,18687,-101.575], Tmin=(907.347,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(210.834,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(123trioxane) + radical(CCsJOO)"""),
)

species(
    label = 'O=CCC1C=C[CH]COO1(1107)',
    structure = SMILES('O=CCC1C=C[CH]COO1'),
    E0 = (-76.0546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0508553,0.0353158,0.000167589,-2.80426e-07,1.2076e-10,-8951.36,33.0022], Tmin=(100,'K'), Tmax=(918.781,'K')), NASAPolynomial(coeffs=[42.1806,-0.00903149,1.22239e-05,-2.42602e-09,1.47584e-13,-22600.1,-199.208], Tmin=(918.781,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-76.0546,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cycloheptane) + radical(C=CCJCO)"""),
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
    label = 'C=CC=CC(C)O[O](1110)',
    structure = SMILES('C=CC=CC(C)O[O]'),
    E0 = (85.9531,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,492.5,1135,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.682316,'amu*angstrom^2'), symmetry=1, barrier=(15.6878,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.681985,'amu*angstrom^2'), symmetry=1, barrier=(15.6802,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.681878,'amu*angstrom^2'), symmetry=1, barrier=(15.6777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.682218,'amu*angstrom^2'), symmetry=1, barrier=(15.6855,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.157551,0.0736837,-4.91122e-05,8.63996e-09,2.57288e-12,10485.5,31.2941], Tmin=(100,'K'), Tmax=(1041.24,'K')), NASAPolynomial(coeffs=[17.2903,0.0273772,-1.05094e-05,1.91929e-09,-1.34129e-13,5860,-57.1324], Tmin=(1041.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(85.9531,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(ROOJ)"""),
)

species(
    label = 'HO2(26)',
    structure = SMILES('[O]O'),
    E0 = (2.67648,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1112.81,1388.53,3298.45],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (33.0067,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02956,-0.00263985,1.5223e-05,-1.71671e-08,6.26738e-12,322.677,4.84428], Tmin=(100,'K'), Tmax=(923.913,'K')), NASAPolynomial(coeffs=[4.15133,0.00191146,-4.11274e-07,6.34957e-11,-4.86385e-15,83.4208,3.09341], Tmin=(923.913,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.67648,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HO2""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'C=CC=CC=CC=O(1013)',
    structure = SMILES('C=CC=CC=CC=O'),
    E0 = (22.2051,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.976937,'amu*angstrom^2'), symmetry=1, barrier=(22.4617,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.976304,'amu*angstrom^2'), symmetry=1, barrier=(22.4472,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.978369,'amu*angstrom^2'), symmetry=1, barrier=(22.4946,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.279345,0.0696767,-4.23342e-05,2.69834e-09,4.17388e-12,2815.02,27.0664], Tmin=(100,'K'), Tmax=(1070.38,'K')), NASAPolynomial(coeffs=[18.1419,0.0249069,-1.04006e-05,1.99581e-09,-1.43366e-13,-2268.19,-66.2073], Tmin=(1070.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(22.2051,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'C=CC=C=CCC=O(1001)',
    structure = SMILES('C=CC=C=CCC=O'),
    E0 = (99.2119,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2750,2850,1437.5,1250,1305,750,350,540,610,2055,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.977632,'amu*angstrom^2'), symmetry=1, barrier=(22.4777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.978595,'amu*angstrom^2'), symmetry=1, barrier=(22.4998,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.976925,'amu*angstrom^2'), symmetry=1, barrier=(22.4614,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.61164,0.0614653,-2.31198e-05,-1.21261e-08,7.8506e-12,12065.4,27.3318], Tmin=(100,'K'), Tmax=(1105.39,'K')), NASAPolynomial(coeffs=[16.7958,0.0283936,-1.28352e-05,2.53455e-09,-1.83671e-13,6930.01,-59.4137], Tmin=(1105.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(99.2119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=CC=CC(O[O])OC=C(1111)',
    structure = SMILES('C=CC=CC(O[O])OC=C'),
    E0 = (-7.58914,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,492.5,1135,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.33793,0.102327,-8.09969e-05,1.74477e-08,4.494e-12,-707.08,37.2085], Tmin=(100,'K'), Tmax=(977.729,'K')), NASAPolynomial(coeffs=[26.3114,0.0217788,-7.38857e-06,1.32677e-09,-9.52644e-14,-7670.48,-103.521], Tmin=(977.729,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-7.58914,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ)"""),
)

species(
    label = 'C=CC=CC(C=CO)O[O](1112)',
    structure = SMILES('C=CC=CC(C=CO)O[O]'),
    E0 = (-22.5945,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,492.5,1135,1000,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,3615,1277.5,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.24912,0.0991465,-7.06446e-05,3.93328e-09,1.02934e-11,-2513.69,38.756], Tmin=(100,'K'), Tmax=(953.47,'K')), NASAPolynomial(coeffs=[27.0215,0.0191784,-5.61618e-06,9.61008e-10,-6.96718e-14,-9660.8,-105.495], Tmin=(953.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-22.5945,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ)"""),
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
    label = 'C=CC=CC([O])CC=O(1114)',
    structure = SMILES('C=CC=CC([O])CC=O'),
    E0 = (-13.7754,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,351.534,351.534,351.534,351.534],'cm^-1')),
        HinderedRotor(inertia=(0.16245,'amu*angstrom^2'), symmetry=1, barrier=(14.2456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16245,'amu*angstrom^2'), symmetry=1, barrier=(14.2456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16245,'amu*angstrom^2'), symmetry=1, barrier=(14.2456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16245,'amu*angstrom^2'), symmetry=1, barrier=(14.2456,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.31885,0.0842774,-6.96587e-05,2.91691e-08,-4.87333e-12,-1492.46,35.579], Tmin=(100,'K'), Tmax=(1431.49,'K')), NASAPolynomial(coeffs=[19.5558,0.0287421,-1.14657e-05,2.06786e-09,-1.40302e-13,-7182.56,-67.4342], Tmin=(1431.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-13.7754,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CC(C)OJ)"""),
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
    E0 = (37.8054,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (244.888,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (405.688,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (314.915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (332.394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (333.622,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (429.007,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (461.058,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (394.358,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (377.642,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (351.126,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (438.262,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (148.263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (309.027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (101.787,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (101.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (20.0681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-4.16942,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (44.351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (81.5188,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (166.388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (29.1155,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (59.4045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (129.201,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (225.986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (68.0946,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (210.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (58.3005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (324.466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (118.743,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (158.031,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (306.211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (121.272,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (229.229,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O2(2)', '[CH2]C=CC=CCC=O(956)'],
    products = ['C=CC=CC(CC=O)O[O](1059)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5.42928e+07,'m^3/(mol*s)'), n=0.107721, Ea=(5.76381,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [C_rad/H/CdCs;Y_rad] for rate rule [C_rad/H/CdCs;O2_birad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=CC=C[CH]O[O](670)', 'CH2CHO(51)'],
    products = ['C=CC=CC(CC=O)O[O](1059)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.3146e+08,'m^3/(mol*s)'), n=-0.359722, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_sec_rad;C_pri_rad] for rate rule [C_rad/H/OneDeO;C_rad/H2/CO]
Euclidian distance = 2.2360679775
family: R_Recombination
Ea raised from -0.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['CH2CHCHCH(170)', '[O]O[CH]CC=O(1096)'],
    products = ['C=CC=CC(CC=O)O[O](1059)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.18e+12,'cm^3/(mol*s)'), n=-0.085, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;C_rad/H/CsO] for rate rule [Cd_pri_rad;C_rad/H/CsO]
Euclidian distance = 2.0
family: R_Recombination
Ea raised from -2.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(19)', 'C=CC=C[C](CC=O)O[O](1097)'],
    products = ['C=CC=CC(CC=O)O[O](1059)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.92e+13,'cm^3/(mol*s)'), n=0.18, Ea=(0.518816,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""Estimated using template [C_rad/OneDe;H_rad] for rate rule [C_rad/OneDeO;H_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['HCO(30)', '[CH2]C(C=CC=C)O[O](1098)'],
    products = ['C=CC=CC(CC=O)O[O](1059)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.81e+13,'cm^3/(mol*s)','*|/',3), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 78 used for C_rad/H2/Cs;CO_pri_rad
Exact match found for rate rule [C_rad/H2/Cs;CO_pri_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(19)', 'C=CC=CC([CH]C=O)O[O](1099)'],
    products = ['C=CC=CC(CC=O)O[O](1059)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.59671e+07,'m^3/(mol*s)'), n=0.0113737, Ea=(2.96199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [C_rad/H/OneDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(19)', 'C=CC=[C]C(CC=O)O[O](1100)'],
    products = ['C=CC=CC(CC=O)O[O](1059)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C2H3(44)', '[CH]=CC(CC=O)O[O](1101)'],
    products = ['C=CC=CC(CC=O)O[O](1059)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.23e+13,'cm^3/(mol*s)','+|-',1.2e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 89 used for Cd_pri_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;Cd_pri_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(19)', 'C=C[C]=CC(CC=O)O[O](1102)'],
    products = ['C=CC=CC(CC=O)O[O](1059)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [Cd_rad/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(19)', 'C=[C]C=CC(CC=O)O[O](1103)'],
    products = ['C=CC=CC(CC=O)O[O](1059)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [Cd_rad/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(19)', 'C=CC=CC(C[C]=O)O[O](1104)'],
    products = ['C=CC=CC(CC=O)O[O](1059)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.53107e+07,'m^3/(mol*s)'), n=-0.133333, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_rad/NonDe;Y_rad] for rate rule [CO_rad/NonDe;H_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(19)', '[CH]=CC=CC(CC=O)O[O](1105)'],
    products = ['C=CC=CC(CC=O)O[O](1059)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.21e+14,'cm^3/(mol*s)','+|-',4.82e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 60 used for H_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CC=CC=C(CC=O)O[O](1108)'],
    products = ['C=CC=CC(CC=O)O[O](1059)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.02873e+09,'s^-1'), n=1.23767, Ea=(163.714,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_3_pentadiene;CH_end;unsaturated_end] for rate rule [1_3_pentadiene;CH3_1;unsaturated_end]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O]OC1CC2OCC=CC12(1109)'],
    products = ['C=CC=CC(CC=O)O[O](1059)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.24e+10,'s^-1'), n=1.27, Ea=(274.47,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [cyclohexene] for rate rule [cyclohexene_2inring]
Euclidian distance = 1.0
family: Intra_Retro_Diels_alder_bicyclic"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=CC=CC(CC=O)O[O](1059)'],
    products = ['[O]OC(CC=O)C1C=CC1(1113)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1,3-butadiene_backbone;C=C_1;C=C_2] for rate rule [1,3-butadiene_backbone;CdH(C)_1;CdH2_2]
Euclidian distance = 1.41421356237
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=CC=CC(CC=O)O[O](1059)'],
    products = ['C=C[CH]C1OOC1CC=O(1086)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra;radadd_intra_O] for rate rule [R5_SS_D;doublebond_intra_HCd_pri;radadd_intra_O]
Euclidian distance = 3.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=CC=CC(CC=O)O[O](1059)'],
    products = ['C=CC=CC1CC([O])OO1(1087)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(132522,'s^-1'), n=1.48406, Ea=(40.6946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSS;multiplebond_intra;radadd_intra] for rate rule [R6_SSS_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Exocyclic
Ea raised from 36.9 to 40.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=CC=CC(CC=O)O[O](1059)'],
    products = ['[CH2]C1C=CC(CC=O)OO1(1088)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(17945.7,'s^-1'), n=1.45333, Ea=(16.4571,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_SSSR;doublebond_intra_2H_pri;radadd_intra] for rate rule [R7_SSSM_D;doublebond_intra_2H_pri;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=CC=CC(CC=O)O[O](1059)'],
    products = ['C=CC=C[C](CC=O)OO(1089)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.3012e-19,'s^-1'), n=9.03667, Ea=(64.9775,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Y_rad_out;Cs_H_out_OneDe] for rate rule [R3H_SS_O;O_rad_out;Cs_H_out_(CdCdCd)]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=CC=CC(CC=O)O[O](1059)'],
    products = ['C=CC=CC([CH]C=O)OO(1090)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(363473,'s^-1'), n=1.92229, Ea=(102.145,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R4H_SSS_OCs;O_rad_out;Cs_H_out_H/CO]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=CC=[C]C(CC=O)OO(1091)'],
    products = ['C=CC=CC(CC=O)O[O](1059)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=CC=CC(C[C]=O)OO(1092)'],
    products = ['C=CC=CC(CC=O)O[O](1059)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(67170.6,'s^-1'), n=1.77845, Ea=(41.7861,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;Y_rad_out;XH_out] for rate rule [R5H_SSSS;CO_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=C[C]=CC(CC=O)OO(1093)'],
    products = ['C=CC=CC(CC=O)O[O](1059)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_single;XH_out] for rate rule [R5H_DSSS;Cd_rad_out_singleDe_Cd;O_H_out]
Euclidian distance = 3.74165738677
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=[C]C=CC(CC=O)OO(1094)'],
    products = ['C=CC=CC(CC=O)O[O](1059)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.62236e+06,'s^-1'), n=1.59783, Ea=(119.552,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_Cd;XH_out] + [R6H_SMSSR;Y_rad_out;XH_out] for rate rule [R6H_SMSSR;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 3.16227766017
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=CC=CC(CC=O)OO(1095)'],
    products = ['C=CC=CC(CC=O)O[O](1059)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R7H;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=CC=CC(CC=O)O[O](1059)'],
    products = ['C=CC1[CH]C(CC=O)OO1(1079)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.81589e+09,'s^-1'), n=0.403324, Ea=(88.7211,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5_linear;doublebond_intra_pri_HCd;radadd_intra] + [R5_SS_D;doublebond_intra_pri;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri_HCd;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=CC=CC(CC=O)O[O](1059)'],
    products = ['C=CC=CC1C[CH]OOO1(1106)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(463580,'s^-1'), n=1.14062, Ea=(231.461,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSS;multiplebond_intra;radadd_intra] for rate rule [R6_SSS_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic
Ea raised from 229.4 to 231.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=CC=CC(CC=O)O[O](1059)'],
    products = ['O=CCC1C=C[CH]COO1(1107)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.76476e+07,'s^-1'), n=0.815689, Ea=(78.927,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;doublebond_intra_pri_2H;radadd_intra] for rate rule [R7_linear;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CO(28)', 'C=CC=CC(C)O[O](1110)'],
    products = ['C=CC=CC(CC=O)O[O](1059)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(274200,'cm^3/(mol*s)'), n=2.53, Ea=(357.732,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 4 used for CO;C_pri/NonDeC
Exact match found for rate rule [CO;C_pri/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=CC=CC(CC=O)O[O](1059)'],
    products = ['HO2(26)', 'C=CC=CC=CC=O(1013)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.91635e+11,'s^-1'), n=0.573333, Ea=(139.369,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2OO] for rate rule [R2OO_HDe_HDe]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=CC=CC(CC=O)O[O](1059)'],
    products = ['HO2(26)', 'C=CC=C=CCC=O(1001)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(3.63e+09,'s^-1'), n=1.11, Ea=(178.657,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R2OO_0H]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=CC=CC(O[O])OC=C(1111)'],
    products = ['C=CC=CC(CC=O)O[O](1059)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 7 used for R_ROR;R1_doublebond_CH2;R2_doublebond_H;R_O_C
Exact match found for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=CC=CC(C=CO)O[O](1112)'],
    products = ['C=CC=CC(CC=O)O[O](1059)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_H;R_O_H]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction34',
    reactants = ['O(20)', 'C=CC=CC([O])CC=O(1114)'],
    products = ['C=CC=CC(CC=O)O[O](1059)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(54738.4,'m^3/(mol*s)'), n=0.884925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.9 to 0 kJ/mol."""),
)

network(
    label = '178',
    isomers = [
        'C=CC=CC(CC=O)O[O](1059)',
    ],
    reactants = [
        ('O2(2)', '[CH2]C=CC=CCC=O(956)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '178',
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

