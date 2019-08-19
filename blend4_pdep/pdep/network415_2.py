species(
    label = 'C4H8O(183)(182)',
    structure = SMILES('CCCC=O'),
    E0 = (-227.103,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.209507,'amu*angstrom^2'), symmetry=1, barrier=(4.81697,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.210933,'amu*angstrom^2'), symmetry=1, barrier=(4.84976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.209775,'amu*angstrom^2'), symmetry=1, barrier=(4.82314,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3282.85,'J/mol'), sigma=(5.66921,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=512.77 K, Pc=40.88 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.11642,0.0389099,-1.83641e-05,3.46285e-09,-1.58355e-13,-27244.3,18.4213], Tmin=(100,'K'), Tmax=(1843.53,'K')), NASAPolynomial(coeffs=[12.0049,0.0224534,-9.04166e-06,1.56257e-09,-1.00129e-13,-31739.7,-37.6376], Tmin=(1843.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-227.103,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH)"""),
)

species(
    label = 'C4H8O(193)(192)',
    structure = SMILES('CCC=CO'),
    E0 = (-225.314,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,1437.5,1250,1305,750,350,180],'cm^-1')),
        HinderedRotor(inertia=(0.844777,'amu*angstrom^2'), symmetry=1, barrier=(19.4231,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.84389,'amu*angstrom^2'), symmetry=1, barrier=(19.4027,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.843515,'amu*angstrom^2'), symmetry=1, barrier=(19.3941,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3571.11,'J/mol'), sigma=(6.00322,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=557.80 K, Pc=37.45 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60692,0.036112,2.60996e-05,-6.96777e-08,3.24567e-11,-26997.6,18.527], Tmin=(100,'K'), Tmax=(923.192,'K')), NASAPolynomial(coeffs=[17.4571,0.0101509,-1.12125e-06,9.70601e-11,-1.0085e-14,-31744.4,-66.5329], Tmin=(923.192,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-225.314,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH)"""),
)

species(
    label = 'C2H5(27)(28)',
    structure = SMILES('C[CH2]'),
    E0 = (107.874,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1190.6,1642.82,1642.96,3622.23,3622.39],'cm^-1')),
        HinderedRotor(inertia=(0.866817,'amu*angstrom^2'), symmetry=1, barrier=(19.9298,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.0611,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2097.75,'J/mol'), sigma=(4.302,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.24186,-0.00356905,4.82667e-05,-5.85401e-08,2.25805e-11,12969,4.44704], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.32196,0.0123931,-4.39681e-06,7.0352e-10,-4.18435e-14,12175.9,0.171104], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(107.874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""C2H5""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'CH2CHO(38)(39)',
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
    label = 'CH3(15)(16)',
    structure = SMILES('[CH3]'),
    E0 = (136.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([604.263,1333.71,1492.19,2836.77,2836.77,3806.92],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (15.0345,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.65718,0.0021266,5.45839e-06,-6.6181e-09,2.46571e-12,16422.7,1.67354], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.97812,0.00579785,-1.97558e-06,3.07298e-10,-1.79174e-14,16509.5,4.72248], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(136.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH3""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'S(337)(336)',
    structure = SMILES('[CH2]CC=O'),
    E0 = (11.2619,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.221237,'amu*angstrom^2'), symmetry=1, barrier=(5.08666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.221178,'amu*angstrom^2'), symmetry=1, barrier=(5.08532,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0712,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3133.67,'J/mol'), sigma=(5.35118,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=489.47 K, Pc=46.4 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.76345,0.0293577,-2.47892e-05,1.49239e-08,-4.38497e-12,1397.08,14.4322], Tmin=(100,'K'), Tmax=(767.858,'K')), NASAPolynomial(coeffs=[4.24224,0.0216537,-9.73869e-06,1.85596e-09,-1.3004e-13,1169.99,7.6886], Tmin=(767.858,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(11.2619,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""CH2CH2CHO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C3H7(69)(68)',
    structure = SMILES('[CH2]CC'),
    E0 = (87.0621,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.0928812,'amu*angstrom^2'), symmetry=1, barrier=(2.13552,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.092914,'amu*angstrom^2'), symmetry=1, barrier=(2.13628,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0877,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2218.31,'J/mol'), sigma=(4.982,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.02815,0.0147023,2.4051e-05,-3.66738e-08,1.38611e-11,10512.1,12.4699], Tmin=(100,'K'), Tmax=(984.464,'K')), NASAPolynomial(coeffs=[6.16543,0.0184495,-6.79029e-06,1.23049e-09,-8.63866e-14,9095.06,-6.67607], Tmin=(984.464,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(87.0621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""npropyl""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'HCO(14)(15)',
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
    label = 'CH2(S)(21)(22)',
    structure = SMILES('[CH2]'),
    E0 = (418.921,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1358.21,2621.43,3089.55],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.19331,-0.00233105,8.15676e-06,-6.62986e-09,1.93233e-12,50366.2,-0.746734], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.13502,0.00289594,-8.16668e-07,1.13573e-10,-6.36263e-15,50504.1,4.06031], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(418.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'C3H6O(201)(200)',
    structure = SMILES('CCC=O'),
    E0 = (-204.33,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.207559,'amu*angstrom^2'), symmetry=1, barrier=(4.77219,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208362,'amu*angstrom^2'), symmetry=1, barrier=(4.79065,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0791,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3133.67,'J/mol'), sigma=(5.35118,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=489.47 K, Pc=46.4 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.90578,0.0240644,-7.06356e-06,-9.81837e-10,5.55825e-13,-24535.9,13.5806], Tmin=(100,'K'), Tmax=(1712.49,'K')), NASAPolynomial(coeffs=[7.69109,0.0189242,-7.84934e-06,1.38273e-09,-8.99057e-14,-27060.1,-14.6647], Tmin=(1712.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-204.33,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""propanal""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CHCHOH(49)(49)',
    structure = SMILES('[CH]=CO'),
    E0 = (118.635,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3120,650,792.5,1650,443.815,444.347,445.091,446.033],'cm^-1')),
        HinderedRotor(inertia=(0.000109588,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1868.27,'J/mol'), sigma=(4.162,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.0819,0.0104649,2.52045e-05,-4.72208e-08,2.08838e-11,14310.7,10.1558], Tmin=(100,'K'), Tmax=(919.505,'K')), NASAPolynomial(coeffs=[11.3957,-0.000441188,1.78863e-06,-3.67188e-10,2.20229e-14,11713.9,-35.0624], Tmin=(919.505,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(118.635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CHCHOH""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'H(3)(3)',
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
    label = 'C4H7OJ_13(6559)',
    structure = SMILES('C[CH]CC=O'),
    E0 = (-29.0891,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,3025,407.5,1350,352.5,2608.8],'cm^-1')),
        HinderedRotor(inertia=(0.636121,'amu*angstrom^2'), symmetry=1, barrier=(14.6257,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147754,'amu*angstrom^2'), symmetry=1, barrier=(3.39716,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147923,'amu*angstrom^2'), symmetry=1, barrier=(3.40103,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2079,0.036465,-1.90363e-05,4.62366e-09,-4.40937e-13,-3431.46,20.959], Tmin=(100,'K'), Tmax=(2385.62,'K')), NASAPolynomial(coeffs=[14.496,0.0158613,-6.08135e-06,1.00335e-09,-6.15459e-14,-9294.39,-49.008], Tmin=(2385.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-29.0891,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), label="""C4H7OJ_13""", comment="""Thermo library: CBS_QB3_1dHR"""),
)

species(
    label = 'CC[CH]C=O(6560)',
    structure = SMILES('CCC=C[O]'),
    E0 = (-83.8518,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.831929,'amu*angstrom^2'), symmetry=1, barrier=(19.1277,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.832282,'amu*angstrom^2'), symmetry=1, barrier=(19.1358,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.02265,0.029098,2.96177e-05,-6.35297e-08,2.76942e-11,-10000.7,18.4445], Tmin=(100,'K'), Tmax=(944.296,'K')), NASAPolynomial(coeffs=[14.1518,0.013926,-3.79527e-06,6.6373e-10,-5.06679e-14,-13905.7,-47.9237], Tmin=(944.296,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-83.8518,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ)"""),
)

species(
    label = '[CH2]CCC=O(3565)',
    structure = SMILES('[CH2]CCC=O'),
    E0 = (-21.857,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.213547,'amu*angstrom^2'), symmetry=1, barrier=(4.90987,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.213568,'amu*angstrom^2'), symmetry=1, barrier=(4.91036,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00471643,'amu*angstrom^2'), symmetry=1, barrier=(4.91084,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.49762,0.0360988,-1.83384e-05,4.20693e-09,-3.76081e-13,-2578.43,18.7989], Tmin=(100,'K'), Tmax=(2478.97,'K')), NASAPolynomial(coeffs=[13.6768,0.0180607,-7.42392e-06,1.27176e-09,-8.00786e-14,-8121.12,-45.2836], Tmin=(2478.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-21.857,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJ)"""),
)

species(
    label = 'CCC[C]=O(6561)',
    structure = SMILES('CCC[C]=O'),
    E0 = (-67.1428,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.225289,'amu*angstrom^2'), symmetry=1, barrier=(5.17985,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.225185,'amu*angstrom^2'), symmetry=1, barrier=(5.17744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.225204,'amu*angstrom^2'), symmetry=1, barrier=(5.17788,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.32545,0.0389501,-2.48868e-05,8.84133e-09,-1.41386e-12,-8017.15,18.4614], Tmin=(100,'K'), Tmax=(1320.47,'K')), NASAPolynomial(coeffs=[6.17416,0.0272916,-1.16432e-05,2.15506e-09,-1.47972e-13,-9033.58,-1.17631], Tmin=(1320.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-67.1428,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCCJ=O)"""),
)

species(
    label = 'CC[CH]C[O](6562)',
    structure = SMILES('CC[CH]C[O]'),
    E0 = (128.691,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,369.546,2216.53,2216.6],'cm^-1')),
        HinderedRotor(inertia=(0.087782,'amu*angstrom^2'), symmetry=1, barrier=(8.50737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0877831,'amu*angstrom^2'), symmetry=1, barrier=(8.50714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0078347,'amu*angstrom^2'), symmetry=1, barrier=(27.3165,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.45303,0.0351607,-1.26562e-05,2.71113e-10,4.32354e-13,15532.6,21.449], Tmin=(100,'K'), Tmax=(1816.99,'K')), NASAPolynomial(coeffs=[11.3564,0.023744,-9.98716e-06,1.77056e-09,-1.15006e-13,10946.2,-30.5387], Tmin=(1816.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(128.691,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo library: CBS_QB3_1dHR + radical(CCJCO) + radical(CCOJ)"""),
)

species(
    label = 'C[CH]CC[O](4512)',
    structure = SMILES('C[CH]CC[O]'),
    E0 = (123.235,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,180,1822.89,1823.25],'cm^-1')),
        HinderedRotor(inertia=(0.227973,'amu*angstrom^2'), symmetry=1, barrier=(5.24154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00222312,'amu*angstrom^2'), symmetry=1, barrier=(5.24358,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.228013,'amu*angstrom^2'), symmetry=1, barrier=(5.24247,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.44243,0.0298368,-6.59903e-06,-1.96417e-09,6.28703e-13,14825.4,18.071], Tmin=(100,'K'), Tmax=(2180.06,'K')), NASAPolynomial(coeffs=[16.562,0.0192354,-8.57322e-06,1.47387e-09,-9.0585e-14,5904.07,-62.7899], Tmin=(2180.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(123.235,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo library: CBS_QB3_1dHR + radical(RCCJC) + radical(CCOJ)"""),
)

species(
    label = 'C[CH]C[CH]O(6555)',
    structure = SMILES('C[CH]C[CH]O'),
    E0 = (77.8274,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3050,390,425,1340,1360,335,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,1822.31,1822.32],'cm^-1')),
        HinderedRotor(inertia=(0.231878,'amu*angstrom^2'), symmetry=1, barrier=(8.39281,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.232194,'amu*angstrom^2'), symmetry=1, barrier=(8.39241,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.232073,'amu*angstrom^2'), symmetry=1, barrier=(8.39204,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00330987,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80373,0.0526741,-6.35307e-05,5.3501e-08,-1.90983e-11,9435.35,22.6847], Tmin=(100,'K'), Tmax=(808.334,'K')), NASAPolynomial(coeffs=[4.40672,0.0330248,-1.4508e-05,2.71096e-09,-1.86245e-13,9235.66,12.0484], Tmin=(808.334,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(77.8274,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo library: CBS_QB3_1dHR + radical(CCsJOH) + radical(RCCJC)"""),
)

species(
    label = '[CH2]CCC[O](6078)',
    structure = SMILES('[CH2]CCC[O]'),
    E0 = (134.035,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3100,440,815,1455,1000,180,1794.8,1794.98],'cm^-1')),
        HinderedRotor(inertia=(0.140914,'amu*angstrom^2'), symmetry=1, barrier=(3.23988,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.140926,'amu*angstrom^2'), symmetry=1, barrier=(3.24018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.140832,'amu*angstrom^2'), symmetry=1, barrier=(3.23801,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.80734,0.0355782,-1.43614e-05,1.69993e-09,5.94106e-14,16154.8,19.429], Tmin=(100,'K'), Tmax=(2227.99,'K')), NASAPolynomial(coeffs=[17.8076,0.0169,-7.34231e-06,1.26213e-09,-7.80086e-14,7422.45,-69.5521], Tmin=(2227.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(134.035,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo library: CBS_QB3_1dHR + radical(RCCJ) + radical(CCOJ)"""),
)

species(
    label = '[CH2]CC[CH]O(3571)',
    structure = SMILES('[CH2]CC[CH]O'),
    E0 = (88.6275,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,334.259,2978.8],'cm^-1')),
        HinderedRotor(inertia=(0.00151418,'amu*angstrom^2'), symmetry=1, barrier=(0.119868,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0875886,'amu*angstrom^2'), symmetry=1, barrier=(6.89694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0879318,'amu*angstrom^2'), symmetry=1, barrier=(6.90675,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.261794,'amu*angstrom^2'), symmetry=1, barrier=(20.5826,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3491.32,'J/mol'), sigma=(6.14151,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=545.34 K, Pc=34.2 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8068,0.0511205,-4.65065e-05,2.54309e-08,-6.04008e-12,10736,21.7335], Tmin=(100,'K'), Tmax=(980.793,'K')), NASAPolynomial(coeffs=[7.24026,0.0289605,-1.26145e-05,2.3931e-09,-1.67688e-13,9670.19,-4.3743], Tmin=(980.793,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(88.6275,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo library: CBS_QB3_1dHR + radical(RCCJ) + radical(CCsJOH)"""),
)

species(
    label = 'CO(10)(11)',
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
    label = 'C3H8(6563)',
    structure = SMILES('CCC'),
    E0 = (-119.987,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.226927,'amu*angstrom^2'), symmetry=1, barrier=(5.2175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.227034,'amu*angstrom^2'), symmetry=1, barrier=(5.21996,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0956,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.06336,0.0129215,3.47033e-05,-4.70946e-08,1.71384e-11,-14390.6,10.7836], Tmin=(100,'K'), Tmax=(989.994,'K')), NASAPolynomial(coeffs=[5.6044,0.0219529,-8.22081e-06,1.50104e-09,-1.05635e-13,-15839.4,-6.22623], Tmin=(989.994,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-119.987,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), label="""C3H8""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=COCC(3569)',
    structure = SMILES('C=COCC'),
    E0 = (-185.507,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.907668,'amu*angstrom^2'), symmetry=1, barrier=(20.8691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.910452,'amu*angstrom^2'), symmetry=1, barrier=(20.9331,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.909928,'amu*angstrom^2'), symmetry=1, barrier=(20.921,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54653,0.0376688,2.2728e-05,-6.5099e-08,3.02379e-11,-22208,17.4901], Tmin=(100,'K'), Tmax=(933.676,'K')), NASAPolynomial(coeffs=[17.2763,0.0117625,-2.29581e-06,3.53792e-10,-2.90291e-14,-26953.4,-67.0005], Tmin=(933.676,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-185.507,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'OH(5)(5)',
    structure = SMILES('[OH]'),
    E0 = (28.372,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3287.46],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (17.0073,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.4858,0.00133397,-4.70043e-06,5.64379e-09,-2.06318e-12,3411.96,1.99788], Tmin=(100,'K'), Tmax=(1005.25,'K')), NASAPolynomial(coeffs=[2.88225,0.00103869,-2.35652e-07,1.40229e-11,6.34581e-16,3669.56,5.59053], Tmin=(1005.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.372,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'buten1yl1(1222)',
    structure = SMILES('[CH]=CCC'),
    E0 = (231.568,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.43749,'amu*angstrom^2'), symmetry=1, barrier=(10.0588,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.437485,'amu*angstrom^2'), symmetry=1, barrier=(10.0586,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.52157,0.0259936,6.08981e-06,-2.23098e-08,9.36555e-12,27910.1,15.2845], Tmin=(100,'K'), Tmax=(1012.54,'K')), NASAPolynomial(coeffs=[7.74082,0.0203541,-7.74605e-06,1.41035e-09,-9.84411e-14,26085.3,-13.7521], Tmin=(1012.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(231.568,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), label="""buten1yl1""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'hydroxyl1allyl(6564)',
    structure = SMILES('[CH2]C=CO'),
    E0 = (-25.0312,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.978343,'amu*angstrom^2'), symmetry=1, barrier=(22.494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.975892,'amu*angstrom^2'), symmetry=1, barrier=(22.4377,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0712,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5477,0.0242083,5.90889e-06,-2.6575e-08,1.23966e-11,-2951.28,13.5324], Tmin=(100,'K'), Tmax=(959.506,'K')), NASAPolynomial(coeffs=[10.1178,0.011534,-3.79865e-06,6.81299e-10,-4.93305e-14,-5273.27,-27.2058], Tmin=(959.506,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-25.0312,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""hydroxyl1allyl""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C4H7O_1_2(5368)',
    structure = SMILES('C[CH]C=CO'),
    E0 = (-53.7053,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,226.636],'cm^-1')),
        HinderedRotor(inertia=(0.285044,'amu*angstrom^2'), symmetry=1, barrier=(10.3913,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.285138,'amu*angstrom^2'), symmetry=1, barrier=(10.3912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.708449,'amu*angstrom^2'), symmetry=1, barrier=(25.8169,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96321,0.0448423,-3.37197e-05,1.39447e-08,-2.46162e-12,-6386.12,17.3708], Tmin=(100,'K'), Tmax=(1288.78,'K')), NASAPolynomial(coeffs=[8.38579,0.0249084,-1.05189e-05,1.94326e-09,-1.33564e-13,-8041.57,-15.2438], Tmin=(1288.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-53.7053,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), label="""C4H7O_1_2""", comment="""Thermo library: CBS_QB3_1dHR"""),
)

species(
    label = 'C4H7OJ_9(6552)',
    structure = SMILES('[CH2]CC=CO'),
    E0 = (20.2278,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,1437.5,1250,1305,750,350,180],'cm^-1')),
        HinderedRotor(inertia=(0.536333,'amu*angstrom^2'), symmetry=1, barrier=(12.3313,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.536142,'amu*angstrom^2'), symmetry=1, barrier=(12.327,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.536463,'amu*angstrom^2'), symmetry=1, barrier=(12.3343,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59785,0.0501872,-4.48597e-05,2.16972e-08,-4.27051e-12,2521.53,18.5107], Tmin=(100,'K'), Tmax=(1214.93,'K')), NASAPolynomial(coeffs=[10.6951,0.020236,-7.88095e-06,1.40601e-09,-9.51517e-14,311.029,-27.1492], Tmin=(1214.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.2278,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), label="""C4H7OJ_9""", comment="""Thermo library: CBS_QB3_1dHR"""),
)

species(
    label = 'CC[C]=CO(6565)',
    structure = SMILES('CC[C]=CO'),
    E0 = (12.5274,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,180],'cm^-1')),
        HinderedRotor(inertia=(0.804828,'amu*angstrom^2'), symmetry=1, barrier=(18.5046,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.805933,'amu*angstrom^2'), symmetry=1, barrier=(18.53,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.805259,'amu*angstrom^2'), symmetry=1, barrier=(18.5145,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58049,0.0406363,1.55919e-06,-4.06894e-08,2.15864e-11,1605.48,19.0955], Tmin=(100,'K'), Tmax=(920.473,'K')), NASAPolynomial(coeffs=[16.2429,0.00969272,-1.42287e-06,1.51705e-10,-1.18792e-14,-2482.2,-57.9691], Tmin=(920.473,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(12.5274,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S)"""),
)

species(
    label = 'CCC=[C]O(6566)',
    structure = SMILES('CCC=[C]O'),
    E0 = (14.4298,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,182.581],'cm^-1')),
        HinderedRotor(inertia=(0.593889,'amu*angstrom^2'), symmetry=1, barrier=(14.0608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.603947,'amu*angstrom^2'), symmetry=1, barrier=(14.0598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.59215,'amu*angstrom^2'), symmetry=1, barrier=(14.058,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81784,0.0378633,-7.9602e-07,-3.08042e-08,1.61556e-11,1823.36,20.8982], Tmin=(100,'K'), Tmax=(935.822,'K')), NASAPolynomial(coeffs=[13.3179,0.0142749,-3.96708e-06,6.48519e-10,-4.57877e-14,-1448.57,-39.8018], Tmin=(935.822,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(14.4298,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C[CH]CO(6554)',
    structure = SMILES('[CH2]C[CH]CO'),
    E0 = (108.232,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,368.593,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0031139,'amu*angstrom^2'), symmetry=1, barrier=(7.40873,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0768567,'amu*angstrom^2'), symmetry=1, barrier=(7.40938,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.256172,'amu*angstrom^2'), symmetry=1, barrier=(24.6931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0103795,'amu*angstrom^2'), symmetry=1, barrier=(24.6933,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75546,0.0430739,-2.65447e-05,8.29314e-09,-1.0487e-12,13103.4,24.8756], Tmin=(100,'K'), Tmax=(1816.1,'K')), NASAPolynomial(coeffs=[12.2608,0.0199356,-7.43372e-06,1.27774e-09,-8.29726e-14,9287.6,-32.0751], Tmin=(1816.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(108.232,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo library: CBS_QB3_1dHR + radical(CCJCO) + radical(RCCJ)"""),
)

species(
    label = 'propen1ol(2833)',
    structure = SMILES('CC=CO'),
    E0 = (-166.942,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680],'cm^-1')),
        HinderedRotor(inertia=(0.636476,'amu*angstrom^2'), symmetry=1, barrier=(14.6338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.637477,'amu*angstrom^2'), symmetry=1, barrier=(14.6568,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0791,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.43518,0.0321429,-2.11365e-05,7.43426e-09,-1.09824e-12,-20020.5,13.0756], Tmin=(100,'K'), Tmax=(1529.28,'K')), NASAPolynomial(coeffs=[8.00322,0.0175792,-6.85177e-06,1.2071e-09,-8.02562e-14,-21723.5,-16.1523], Tmin=(1529.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-166.942,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""propen1ol""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CC[C]CO(6567)',
    structure = SMILES('CC[C]CO'),
    E0 = (124.532,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,3674.16],'cm^-1')),
        HinderedRotor(inertia=(0.545249,'amu*angstrom^2'), symmetry=1, barrier=(12.5364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.54117,'amu*angstrom^2'), symmetry=1, barrier=(58.4264,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.53982,'amu*angstrom^2'), symmetry=1, barrier=(58.3956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.54093,'amu*angstrom^2'), symmetry=1, barrier=(58.421,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.12974,0.0473344,-2.25595e-05,-2.56043e-08,3.17524e-11,15039,28.4089], Tmin=(100,'K'), Tmax=(516.954,'K')), NASAPolynomial(coeffs=[5.03255,0.0349037,-1.5594e-05,2.9452e-09,-2.04871e-13,14604.9,15.0235], Tmin=(516.954,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(124.532,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsOsHH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'CCC[C-]=[OH+](6568)',
    structure = SMILES('CCC[C-]=[OH+]'),
    E0 = (205.569,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,970.714,1208.95,2880,3200],'cm^-1')),
        HinderedRotor(inertia=(0.0482186,'amu*angstrom^2'), symmetry=1, barrier=(1.10864,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0482186,'amu*angstrom^2'), symmetry=1, barrier=(1.10864,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0482186,'amu*angstrom^2'), symmetry=1, barrier=(1.10864,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28378,0.0368221,-1.71893e-05,3.21904e-09,-2.13942e-13,24727.3,14.3205], Tmin=(100,'K'), Tmax=(2971.11,'K')), NASAPolynomial(coeffs=[37.0771,-0.00303094,8.2151e-08,-1.71452e-11,4.5754e-15,2155.98,-189.703], Tmin=(2971.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.569,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(CsJ2_singlet-CsH)"""),
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
    label = 'Ar(8)',
    structure = SMILES('[Ar]'),
    E0 = (-6.19426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (39.348,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1134.93,'J/mol'), sigma=(3.33,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,-745,4.3663], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,-745,4.3663], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-6.19426,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ar""", comment="""Thermo library: BurkeH2O2"""),
)

transitionState(
    label = 'TS1',
    E0 = (109.103,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (147.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (182.703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (119.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (130.902,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (189.935,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (144.649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (150.992,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (186.635,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (86.1954,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (159.009,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (113.601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (214.591,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (118.526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (128.293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-81.4477,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (259.992,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (127.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (111.156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (227.606,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (163.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (232.02,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (224.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (226.222,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (151.552,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (100.973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (152.028,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (133.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (251.99,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (211.972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (242.607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction21',
    reactants = ['C2H5(27)(28)', 'CH2CHO(38)(39)'],
    products = ['C4H8O(183)(182)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.33778e+08,'m^3/(mol*s)'), n=-0.3495, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_pri_rad;C_rad/H2/Cs] for rate rule [C_rad/H2/CO;C_rad/H2/Cs]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['CH3(15)(16)', 'S(337)(336)'],
    products = ['C4H8O(183)(182)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.23e+15,'cm^3/(mol*s)'), n=-0.562, Ea=(0.085772,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 10 used for C_methyl;C_rad/H2/Cs
Exact match found for rate rule [C_rad/H2/Cs;C_methyl]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)(3)', 'C4H7OJ_13(6559)'],
    products = ['C4H8O(183)(182)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C3H7(69)(68)', 'HCO(14)(15)'],
    products = ['C4H8O(183)(182)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.81e+13,'cm^3/(mol*s)','*|/',3), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 78 used for C_rad/H2/Cs;CO_pri_rad
Exact match found for rate rule [CO_pri_rad;C_rad/H2/Cs]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)(3)', 'CC[CH]C=O(6560)'],
    products = ['C4H8O(183)(182)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.59671e+07,'m^3/(mol*s)'), n=0.0113737, Ea=(2.96199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [C_rad/H/OneDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(3)(3)', '[CH2]CCC=O(3565)'],
    products = ['C4H8O(183)(182)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(3)(3)', 'CCC[C]=O(6561)'],
    products = ['C4H8O(183)(182)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.53107e+07,'m^3/(mol*s)'), n=-0.133333, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;CO_rad/NonDe] for rate rule [H_rad;CO_rad/NonDe]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CC[CH]C[O](6562)'],
    products = ['C4H8O(183)(182)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.48818e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C[CH]CC[O](4512)'],
    products = ['C4H8O(183)(182)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C[CH]C[CH]O(6555)'],
    products = ['C4H8O(183)(182)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]CCC[O](6078)'],
    products = ['C4H8O(183)(182)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]CC[CH]O(3571)'],
    products = ['C4H8O(183)(182)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CH2(S)(21)(22)', 'C3H6O(201)(200)'],
    products = ['C4H8O(183)(182)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction14',
    reactants = ['CO(10)(11)', 'C3H8(6563)'],
    products = ['C4H8O(183)(182)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(548400,'cm^3/(mol*s)'), n=2.53, Ea=(357.732,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 4 used for CO;C_pri/NonDeC
Exact match found for rate rule [CO;C_pri/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 6.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=COCC(3569)'],
    products = ['C4H8O(183)(182)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 7 used for R_ROR;R1_doublebond_CH2;R2_doublebond_H;R_O_C
Exact match found for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C4H8O(193)(192)'],
    products = ['C4H8O(183)(182)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_H;R_O_H]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction36',
    reactants = ['OH(5)(5)', 'buten1yl1(1222)'],
    products = ['C4H8O(193)(192)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.44289e+07,'m^3/(mol*s)'), n=0.0225, Ea=(0.0523,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;Cd_pri_rad] + [O_pri_rad;Y_rad] for rate rule [O_pri_rad;Cd_pri_rad]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['H(3)(3)', 'CC[CH]C=O(6560)'],
    products = ['C4H8O(193)(192)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7e+11,'cm^3/(mol*s)'), n=0.493, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2500,'K'), comment="""From training reaction 25 used for H_rad;O_rad/OneDe
Exact match found for rate rule [H_rad;O_rad/OneDe]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['CH3(15)(16)', 'hydroxyl1allyl(6564)'],
    products = ['C4H8O(193)(192)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.02e+14,'cm^3/(mol*s)'), n=-0.32, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 127 used for C_rad/H2/Cd;C_methyl
Exact match found for rate rule [C_rad/H2/Cd;C_methyl]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction39',
    reactants = ['C2H5(27)(28)', 'CHCHOH(49)(49)'],
    products = ['C4H8O(193)(192)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(6.41031e+09,'m^3/(mol*s)'), n=-0.9855, Ea=(1.09644,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_rad;C_pri_rad] for rate rule [Cd_pri_rad;C_rad/H2/Cs]
Euclidian distance = 1.41421356237
family: R_Recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['H(3)(3)', 'C4H7O_1_2(5368)'],
    products = ['C4H8O(193)(192)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.71464e+07,'m^3/(mol*s)'), n=0.107721, Ea=(5.76381,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 36 used for C_rad/H/CdCs;H_rad
Exact match found for rate rule [C_rad/H/CdCs;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)(3)', 'C4H7OJ_9(6552)'],
    products = ['C4H8O(193)(192)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction42',
    reactants = ['H(3)(3)', 'CC[C]=CO(6565)'],
    products = ['C4H8O(193)(192)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction43',
    reactants = ['H(3)(3)', 'CCC=[C]O(6566)'],
    products = ['C4H8O(193)(192)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction44',
    reactants = ['CC[CH]C[O](6562)'],
    products = ['C4H8O(193)(192)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction45',
    reactants = ['C[CH]C[CH]O(6555)'],
    products = ['C4H8O(193)(192)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(5.10299e+10,'s^-1'), n=0.2847, Ea=(23.1459,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad_NDe] + [R2radExo;Y_rad;XH_Rrad_NDe] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]CC[CH]O(3571)'],
    products = ['C4H8O(193)(192)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2]C[CH]CO(6554)'],
    products = ['C4H8O(193)(192)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['CH2(S)(21)(22)', 'propen1ol(2833)'],
    products = ['C4H8O(193)(192)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.87e+13,'cm^3/(mol*s)','*|/',0.25), n=-0.146, Ea=(0.0118826,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 6 used for carbene;C_pri/Cd
Exact match found for rate rule [carbene;C_pri/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene"""),
)

reaction(
    label = 'reaction48',
    reactants = ['CC[C]CO(6567)'],
    products = ['C4H8O(193)(192)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(3.23663e+16,'s^-1'), n=-0.885455, Ea=(87.4392,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [CsJ2-C;CsJ2(CsC);CH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction49',
    reactants = ['CCC[C-]=[OH+](6568)'],
    products = ['C4H8O(193)(192)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.90176e+13,'s^-1'), n=-0.332469, Ea=(37.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [singletcarbene_CH;singletcarbene;CH2(C)] + [CsJ2-C;singletcarbene;CH] for rate rule [CsJ2-C;singletcarbene;CH2(C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

network(
    label = '415',
    isomers = [
        'C4H8O(183)(182)',
        'C4H8O(193)(192)',
    ],
    reactants = [
        ('C2H5(27)(28)', 'CH2CHO(38)(39)'),
        ('CH3(15)(16)', 'S(337)(336)'),
        ('C3H7(69)(68)', 'HCO(14)(15)'),
        ('CH2(S)(21)(22)', 'C3H6O(201)(200)'),
        ('C2H5(27)(28)', 'CHCHOH(49)(49)'),
    ],
    bathGas = {
        'Ne': 0.333333,
        'N2': 0.333333,
        'Ar(8)': 0.333333,
    },
)

pressureDependence(
    label = '415',
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

