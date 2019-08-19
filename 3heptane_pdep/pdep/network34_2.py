species(
    label = 'npropylperoxy(106)',
    structure = SMILES('CCCO[O]'),
    E0 = (-58.7349,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.153611,'amu*angstrom^2'), symmetry=1, barrier=(3.53182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153632,'amu*angstrom^2'), symmetry=1, barrier=(3.53231,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153656,'amu*angstrom^2'), symmetry=1, barrier=(3.53285,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (75.0865,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3494.2,'J/mol'), sigma=(6.02528,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=545.79 K, Pc=36.25 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.15077,0.0341093,-6.28775e-06,-1.16905e-08,5.70826e-12,-6992.04,19.9651], Tmin=(100,'K'), Tmax=(1072.55,'K')), NASAPolynomial(coeffs=[8.78442,0.0227415,-9.09063e-06,1.67579e-09,-1.16717e-13,-9184.16,-16.0886], Tmin=(1072.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-58.7349,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), label="""npropylperoxy""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'QOOH_1(211)',
    structure = SMILES('[CH2]CCOO'),
    E0 = (8.94627,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,755.791],'cm^-1')),
        HinderedRotor(inertia=(0.704177,'amu*angstrom^2'), symmetry=1, barrier=(16.1904,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0126973,'amu*angstrom^2'), symmetry=1, barrier=(5.14686,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223859,'amu*angstrom^2'), symmetry=1, barrier=(5.14697,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223857,'amu*angstrom^2'), symmetry=1, barrier=(5.1469,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (75.0865,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73075,0.0457841,-3.51689e-05,1.42448e-08,-2.36092e-12,1161.01,22.4746], Tmin=(100,'K'), Tmax=(1415.62,'K')), NASAPolynomial(coeffs=[10.8293,0.0200752,-7.92784e-06,1.41612e-09,-9.53847e-14,-1415.05,-24.5834], Tmin=(1415.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(8.94627,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), label="""QOOH_1""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53764,-0.00122827,5.36759e-06,-4.93128e-09,1.45955e-12,-1037.99,4.6718], Tmin=(100,'K'), Tmax=(1087.71,'K')), NASAPolynomial(coeffs=[3.16427,0.00169454,-8.00335e-07,1.5903e-10,-1.14891e-14,-1048.45,6.08303], Tmin=(1087.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.62178,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'npropyl(70)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.02816,0.0147023,2.40511e-05,-3.6674e-08,1.38612e-11,10512.1,12.4699], Tmin=(100,'K'), Tmax=(984.463,'K')), NASAPolynomial(coeffs=[6.16542,0.0184495,-6.7903e-06,1.23049e-09,-8.63868e-14,9095.06,-6.676], Tmin=(984.463,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(87.0621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""npropyl""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CH2(S)(24)',
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
    label = 'CH3CH2OO(54)',
    structure = SMILES('CCO[O]'),
    E0 = (-35.8481,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
        HinderedRotor(inertia=(0.196253,'amu*angstrom^2'), symmetry=1, barrier=(4.51224,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197956,'amu*angstrom^2'), symmetry=1, barrier=(4.5514,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (61.0599,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3284.22,'J/mol'), sigma=(4.037,'angstroms'), dipoleMoment=(1.3,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.88894,0.0185233,1.04542e-05,-2.35199e-08,9.42773e-12,-4266.36,15.0609], Tmin=(100,'K'), Tmax=(1009.3,'K')), NASAPolynomial(coeffs=[7.15149,0.0156281,-6.0467e-06,1.12072e-09,-7.93637e-14,-5839.77,-9.07495], Tmin=(1009.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-35.8481,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""CH3CH2OO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C2H5(32)',
    structure = SMILES('C[CH2]'),
    E0 = (107.874,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1190.59,1642.33,1643.48,3621.68,3622.96],'cm^-1')),
        HinderedRotor(inertia=(0.866827,'amu*angstrom^2'), symmetry=1, barrier=(19.9301,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.0611,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2097.75,'J/mol'), sigma=(4.302,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.24186,-0.00356905,4.82667e-05,-5.85401e-08,2.25805e-11,12969,4.44704], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.32196,0.0123931,-4.39681e-06,7.0352e-10,-4.18435e-14,12175.9,0.171104], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(107.874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""C2H5""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH2]O[O](203)',
    structure = SMILES('[CH2]O[O]'),
    E0 = (206.89,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,291.401,3698.4],'cm^-1')),
        HinderedRotor(inertia=(0.37416,'amu*angstrom^2'), symmetry=1, barrier=(22.5475,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.92061,0.0189601,-2.08419e-05,1.10972e-08,-2.15337e-12,24926,11.4734], Tmin=(100,'K'), Tmax=(1506.17,'K')), NASAPolynomial(coeffs=[8.0388,0.00135553,6.86082e-07,-2.00091e-10,1.53442e-14,23839.3,-13.8045], Tmin=(1506.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(206.89,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CsJOOH) + radical(ROOJ)"""),
)

species(
    label = '[CH2]CO[O](102)',
    structure = SMILES('[CH2]CO[O]'),
    E0 = (188.389,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,492.5,1135,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.346471,'amu*angstrom^2'), symmetry=1, barrier=(7.96605,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000699441,'amu*angstrom^2'), symmetry=1, barrier=(7.94149,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.81263,0.0259367,-2.20523e-05,1.10177e-08,-2.31288e-12,22700.9,17.2501], Tmin=(100,'K'), Tmax=(1124.15,'K')), NASAPolynomial(coeffs=[6.4098,0.0131371,-4.9735e-06,8.89322e-10,-6.04388e-14,21892.2,-0.525196], Tmin=(1124.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.389,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(CJCOOH)"""),
)

species(
    label = 'CH3(18)',
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
    label = 'H(3)',
    structure = SMILES('[H]'),
    E0 = (211.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (1.00794,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1205.6,'J/mol'), sigma=(2.05,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,-2.38914e-13,3.12709e-16,-1.33367e-19,1.7499e-23,25472.7,-0.459566], Tmin=(100,'K'), Tmax=(4383.16,'K')), NASAPolynomial(coeffs=[2.50003,-3.04997e-08,1.01101e-11,-1.48797e-15,8.20356e-20,25472.7,-0.459785], Tmin=(4383.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.792,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'C[CH]CO[O](212)',
    structure = SMILES('C[CH]CO[O]'),
    E0 = (151.534,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.24901,'amu*angstrom^2'), symmetry=1, barrier=(5.72523,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.248616,'amu*angstrom^2'), symmetry=1, barrier=(5.71618,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.248551,'amu*angstrom^2'), symmetry=1, barrier=(5.71468,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0785,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2481,0.0410859,-4.15988e-05,2.81248e-08,-8.40506e-12,18286.1,21.4524], Tmin=(100,'K'), Tmax=(790.433,'K')), NASAPolynomial(coeffs=[5.4752,0.0247549,-1.06073e-05,1.98572e-09,-1.37648e-13,17776,6.6425], Tmin=(790.433,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(151.534,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJCOOH) + radical(ROOJ)"""),
)

species(
    label = 'CC[CH]O[O](182)',
    structure = SMILES('CC[CH]O[O]'),
    E0 = (139.702,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200.873],'cm^-1')),
        HinderedRotor(inertia=(0.242257,'amu*angstrom^2'), symmetry=1, barrier=(7.02304,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.229193,'amu*angstrom^2'), symmetry=1, barrier=(7.02158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000883538,'amu*angstrom^2'), symmetry=1, barrier=(6.99943,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0785,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.20784,0.0420296,-4.33791e-05,3.01325e-08,-9.24763e-12,16864.4,19.4406], Tmin=(100,'K'), Tmax=(771.373,'K')), NASAPolynomial(coeffs=[5.43298,0.0253065,-1.08618e-05,2.03087e-09,-1.40557e-13,16366.8,4.71819], Tmin=(771.373,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(139.702,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(CCsJOOH)"""),
)

species(
    label = '[CH2]CCO[O](201)',
    structure = SMILES('[CH2]CCO[O]'),
    E0 = (156.366,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,492.5,1135,1000,3000,3100,440,815,1455,1000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.220045,'amu*angstrom^2'), symmetry=1, barrier=(5.05926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220031,'amu*angstrom^2'), symmetry=1, barrier=(5.05896,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220024,'amu*angstrom^2'), symmetry=1, barrier=(5.05879,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0785,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.30676,0.0392486,-3.42472e-05,1.87088e-08,-4.51413e-12,18865.7,20.7232], Tmin=(100,'K'), Tmax=(960.537,'K')), NASAPolynomial(coeffs=[6.0396,0.023704,-9.97258e-06,1.86101e-09,-1.29187e-13,18148.6,2.86463], Tmin=(960.537,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(156.366,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(RCCJ)"""),
)

species(
    label = 'CC[CH]OO(209)',
    structure = SMILES('CC[CH]OO'),
    E0 = (-12.2999,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3327.51],'cm^-1')),
        HinderedRotor(inertia=(0.00392121,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.279614,'amu*angstrom^2'), symmetry=1, barrier=(8.68178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.614567,'amu*angstrom^2'), symmetry=1, barrier=(18.7467,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.38401,'amu*angstrom^2'), symmetry=1, barrier=(42.2894,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (75.0865,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0636,0.0442419,-3.39821e-05,1.47137e-08,-2.76751e-12,-1411.09,19.1984], Tmin=(100,'K'), Tmax=(1201.28,'K')), NASAPolynomial(coeffs=[7.51866,0.0260777,-1.1301e-05,2.12648e-09,-1.4797e-13,-2721.7,-8.11948], Tmin=(1201.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-12.2999,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCsJOOH)"""),
)

species(
    label = 'QOOH_2(210)',
    structure = SMILES('C[CH]COO'),
    E0 = (0.506212,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (75.0865,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03865,0.0408811,-2.67535e-05,9.01764e-09,-1.26357e-12,133.04,21.7174], Tmin=(100,'K'), Tmax=(1602.58,'K')), NASAPolynomial(coeffs=[9.65915,0.0218605,-8.95037e-06,1.61158e-09,-1.08232e-13,-2309.45,-18.641], Tmin=(1602.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(0.506212,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), label="""QOOH_2""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'HO2(10)',
    structure = SMILES('[O]O'),
    E0 = (2.67648,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1112.8,1388.53,3298.45],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (33.0067,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02957,-0.00263999,1.52235e-05,-1.71679e-08,6.26771e-12,322.677,4.84424], Tmin=(100,'K'), Tmax=(923.901,'K')), NASAPolynomial(coeffs=[4.1513,0.00191152,-4.11308e-07,6.35038e-11,-4.86452e-15,83.4341,3.09359], Tmin=(923.901,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.67648,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HO2""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'C3H6(59)',
    structure = SMILES('C=CC'),
    E0 = (5.9763,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.497558,'amu*angstrom^2'), symmetry=1, barrier=(11.4398,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31912,0.00817957,3.34737e-05,-4.36194e-08,1.58214e-11,749.325,9.54025], Tmin=(100,'K'), Tmax=(983.754,'K')), NASAPolynomial(coeffs=[5.36755,0.0170743,-6.35108e-06,1.1662e-09,-8.27621e-14,-487.137,-4.54465], Tmin=(983.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.9763,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""C3H6""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'npropoxy(213)',
    structure = SMILES('CCC[O]'),
    E0 = (-50.3873,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,553.578,4000],'cm^-1')),
        HinderedRotor(inertia=(0.953815,'amu*angstrom^2'), symmetry=1, barrier=(21.9301,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.537925,'amu*angstrom^2'), symmetry=1, barrier=(12.368,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (59.0871,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.49,0.0272742,1.6647e-06,-1.66154e-08,7.01483e-12,-6000.76,15.8164], Tmin=(100,'K'), Tmax=(1050.44,'K')), NASAPolynomial(coeffs=[7.46142,0.0213308,-8.39406e-06,1.53869e-09,-1.07024e-13,-7761.72,-11.8231], Tmin=(1050.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-50.3873,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), label="""npropoxy""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'O(4)',
    structure = SMILES('[O]'),
    E0 = (243.005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (15.9994,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,-2.38914e-13,3.12709e-16,-1.33367e-19,1.7499e-23,29226.7,5.11107], Tmin=(100,'K'), Tmax=(4383.16,'K')), NASAPolynomial(coeffs=[2.50003,-3.04997e-08,1.01101e-11,-1.48797e-15,8.20356e-20,29226.7,5.11085], Tmin=(4383.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.005,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'CH2CHCH2OOH(252)',
    structure = SMILES('C=CCOO'),
    E0 = (-69.8212,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(0.489916,'amu*angstrom^2'), symmetry=1, barrier=(11.2641,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.489805,'amu*angstrom^2'), symmetry=1, barrier=(11.2616,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.48987,'amu*angstrom^2'), symmetry=1, barrier=(11.2631,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (74.0785,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.84686,0.0510074,-6.52616e-05,5.25166e-08,-1.77878e-11,-8323.49,18.0257], Tmin=(100,'K'), Tmax=(775.707,'K')), NASAPolynomial(coeffs=[6.14842,0.0261422,-1.19895e-05,2.27271e-09,-1.57396e-13,-8910.1,-1.11392], Tmin=(775.707,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-69.8212,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), label="""CH2CHCH2OOH""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C2H4(29)',
    structure = SMILES('C=C'),
    E0 = (42.0619,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2334.71,'J/mol'), sigma=(3.971,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.9592,-0.00757051,5.7099e-05,-6.91588e-08,2.69884e-11,5089.78,4.0973], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.99183,0.0104834,-3.71721e-06,5.94628e-10,-3.5363e-14,4268.66,-0.269082], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(42.0619,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""C2H4""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH2]OO(239)',
    structure = SMILES('[CH2]OO'),
    E0 = (54.8878,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.878705,'amu*angstrom^2'), symmetry=1, barrier=(20.2032,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (47.0333,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.85643,0.0201434,-7.5271e-06,-9.6693e-09,6.63648e-12,6647.28,10.9504], Tmin=(100,'K'), Tmax=(928.338,'K')), NASAPolynomial(coeffs=[9.91701,0.00263656,-1.08425e-07,-1.03822e-11,-4.89454e-16,4779.82,-25.5853], Tmin=(928.338,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(54.8878,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CsJOOH)"""),
)

species(
    label = 'OH(5)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.48579,0.001334,-4.70054e-06,5.64393e-09,-2.06324e-12,3411.96,1.99789], Tmin=(100,'K'), Tmax=(1005.24,'K')), NASAPolynomial(coeffs=[2.88226,0.00103867,-2.35641e-07,1.40204e-11,6.3479e-16,3669.56,5.59047], Tmin=(1005.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.372,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = '[CH2]CC[O](253)',
    structure = SMILES('[CH2]CC[O]'),
    E0 = (158.075,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,1593.31,1593.43],'cm^-1')),
        HinderedRotor(inertia=(0.00365688,'amu*angstrom^2'), symmetry=1, barrier=(6.58664,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.286463,'amu*angstrom^2'), symmetry=1, barrier=(6.58636,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0791,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.12517,0.020232,6.23149e-07,-6.19731e-09,1.72058e-12,19042.6,16.7651], Tmin=(100,'K'), Tmax=(1581.51,'K')), NASAPolynomial(coeffs=[6.2457,0.0215752,-9.41043e-06,1.72475e-09,-1.15421e-13,16900.5,-3.37162], Tmin=(1581.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(158.075,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJ) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C[CH2](83)',
    structure = SMILES('[CH2]C[CH2]'),
    E0 = (290.499,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.0167608,'amu*angstrom^2'), symmetry=1, barrier=(29.9602,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000703342,'amu*angstrom^2'), symmetry=1, barrier=(7.98579,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.12863,0.0145638,1.31003e-05,-2.16419e-08,7.8594e-12,34974.2,14.113], Tmin=(100,'K'), Tmax=(1043.7,'K')), NASAPolynomial(coeffs=[5.43609,0.0169627,-6.50459e-06,1.20126e-09,-8.4342e-14,33880.3,-0.0512246], Tmin=(1043.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(290.499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = 'C2H4(T)(79)',
    structure = SMILES('[CH2][CH2]'),
    E0 = (318.146,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,180,1436.16,1437.36,2688.4,2689.63],'cm^-1')),
        HinderedRotor(inertia=(0.0257474,'amu*angstrom^2'), symmetry=1, barrier=(17.2422,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.40737,0.0100311,6.40983e-06,-1.41299e-08,5.92706e-12,38288.2,6.11699], Tmin=(100,'K'), Tmax=(954.25,'K')), NASAPolynomial(coeffs=[5.52245,0.0085618,-2.90747e-06,5.02363e-10,-3.4458e-14,37547.8,-5.75253], Tmin=(954.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.146,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""C2H4(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2][CH]COO(254)',
    structure = SMILES('[CH2][CH]COO'),
    E0 = (204.775,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,3615,1310,387.5,850,1000,3025,407.5,1350,352.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.12454,'amu*angstrom^2'), symmetry=1, barrier=(2.86342,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124333,'amu*angstrom^2'), symmetry=1, barrier=(2.85866,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124615,'amu*angstrom^2'), symmetry=1, barrier=(2.86514,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00479815,'amu*angstrom^2'), symmetry=1, barrier=(2.85803,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0785,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.12915,0.0440999,-4.25342e-05,2.44512e-08,-6.06471e-12,24693.7,22.9094], Tmin=(100,'K'), Tmax=(944.857,'K')), NASAPolynomial(coeffs=[6.84545,0.024134,-1.08377e-05,2.08721e-09,-1.4749e-13,23802.4,0.42344], Tmin=(944.857,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.775,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJ) + radical(CCJCOOH)"""),
)

species(
    label = '[CH2]C[CH]OO(255)',
    structure = SMILES('[CH2]C[CH]OO'),
    E0 = (192.943,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,3615,1310,387.5,850,1000,3025,407.5,1350,352.5,2973.75],'cm^-1')),
        HinderedRotor(inertia=(0.11029,'amu*angstrom^2'), symmetry=1, barrier=(4.35625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.95001,'amu*angstrom^2'), symmetry=1, barrier=(38.0557,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.41367,'amu*angstrom^2'), symmetry=1, barrier=(16.3795,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0060699,'amu*angstrom^2'), symmetry=1, barrier=(38.1003,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0785,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.09556,0.044945,-4.38712e-05,2.57352e-08,-6.52733e-12,23271.7,20.875], Tmin=(100,'K'), Tmax=(925.023,'K')), NASAPolynomial(coeffs=[6.76408,0.0247575,-1.1136e-05,2.14314e-09,-1.51319e-13,22407.9,-1.28417], Tmin=(925.023,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.943,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCsJOOH) + radical(RCCJ)"""),
)

species(
    label = 'oxetane(256)',
    structure = SMILES('C1COC1'),
    E0 = (-91.7532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0791,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.4757,-0.00369722,8.70656e-05,-1.0561e-07,3.87164e-11,-11002.7,11.3085], Tmin=(100,'K'), Tmax=(960.535,'K')), NASAPolynomial(coeffs=[8.01388,0.0165131,-5.56907e-06,1.07257e-09,-8.24755e-14,-13678.7,-19.7942], Tmin=(960.535,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-91.7532,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), label="""oxetane""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[O]CCCO(257)',
    structure = SMILES('[O]CCCO'),
    E0 = (-214.474,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (75.0865,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05137,0.0482537,-6.19208e-05,5.72111e-08,-2.18955e-11,-25730.3,20.6204], Tmin=(100,'K'), Tmax=(815.087,'K')), NASAPolynomial(coeffs=[2.69252,0.0341488,-1.5797e-05,2.9916e-09,-2.06359e-13,-25470.8,19.8912], Tmin=(815.087,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-214.474,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + radical(CCOJ)"""),
)

species(
    label = 'CH2CH2OOH(55)',
    structure = SMILES('[CH2]COO'),
    E0 = (34.3631,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.257945,'amu*angstrom^2'), symmetry=1, barrier=(5.93066,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.258128,'amu*angstrom^2'), symmetry=1, barrier=(5.93487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0141801,'amu*angstrom^2'), symmetry=1, barrier=(16.6288,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (61.0599,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.44846,0.0313149,-2.37472e-05,9.50342e-09,-1.55883e-12,4191.03,17.4641], Tmin=(100,'K'), Tmax=(1426.94,'K')), NASAPolynomial(coeffs=[8.63688,0.0139676,-5.51186e-06,9.83892e-10,-6.62178e-14,2424.92,-14.5917], Tmin=(1426.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(34.3631,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), label="""CH2CH2OOH""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CH2(T)(17)',
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
    label = '[CH]CCOO(247)',
    structure = SMILES('[CH]CCOO'),
    E0 = (247.329,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1310,387.5,850,1000,320.321,320.324,1352.47,1352.47],'cm^-1')),
        HinderedRotor(inertia=(0.133672,'amu*angstrom^2'), symmetry=1, barrier=(9.73435,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00164293,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133695,'amu*angstrom^2'), symmetry=1, barrier=(9.73433,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00164264,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0785,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92275,0.0455979,-4.03983e-05,1.93855e-08,-3.83586e-12,29821.7,20.5537], Tmin=(100,'K'), Tmax=(1193.07,'K')), NASAPolynomial(coeffs=[9.4651,0.0203106,-8.60526e-06,1.62007e-09,-1.13197e-13,28022,-17.1652], Tmin=(1193.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(247.329,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJ2_triplet)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.61263,-0.00100893,2.49898e-06,-1.43375e-09,2.58635e-13,-1051.1,2.6527], Tmin=(100,'K'), Tmax=(1817.04,'K')), NASAPolynomial(coeffs=[2.97591,0.0016414,-7.19719e-07,1.25377e-10,-7.91522e-15,-1025.85,5.53754], Tmin=(1817.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.69489,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""N2""", comment="""Thermo library: BurkeH2O2"""),
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
    E0 = (93.5627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (314.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (324.663,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (363.326,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (351.494,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (368.158,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (383.073,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (107.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (67.1621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (29.2128,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (64.2751,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (192.618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (149.544,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (126.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (160.825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (130.374,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (186.447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (293.175,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (368.158,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (373.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (416.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (404.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (84.8022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (125.261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (415.443,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (459.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['npropylperoxy(106)'],
    products = ['O2(2)', 'npropyl(70)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.52e+23,'s^-1'), n=-2.71, Ea=(36.4,'kcal/mol'), T0=(1,'K'), comment="""Matched reaction 2 C3H7O2 <=> O2 + C3H7 in R_Recombination/training
This reaction matched rate rule [C_rad/H2/Cs;O2_birad]
family: R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C2H5(32)', '[CH2]O[O](203)'],
    products = ['npropylperoxy(106)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.33778e+08,'m^3/(mol*s)'), n=-0.3495, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H2/Cs;C_pri_rad] for rate rule [C_rad/H2/Cs;C_rad/H2/O]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]CO[O](102)', 'CH3(18)'],
    products = ['npropylperoxy(106)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.23e+15,'cm^3/(mol*s)'), n=-0.562, Ea=(0.085772,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 10 used for C_methyl;C_rad/H2/Cs
Exact match found for rate rule [C_methyl;C_rad/H2/Cs]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', 'C[CH]CO[O](212)'],
    products = ['npropylperoxy(106)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', 'CC[CH]O[O](182)'],
    products = ['npropylperoxy(106)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.18e+12,'cm^3/(mol*s)'), n=-0.085, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H/CsO;Y_rad] for rate rule [C_rad/H/CsO;H_rad]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -2.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH2]CCO[O](201)'],
    products = ['npropylperoxy(106)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2(S)(24)', 'CH3CH2OO(54)'],
    products = ['npropylperoxy(106)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction2',
    reactants = ['npropylperoxy(106)'],
    products = ['CC[CH]OO(209)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.44e+09,'s^-1'), n=1.17, Ea=(165.937,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Matched reaction 242 C3H7O2 <=> C3H7O2-2 in intra_H_migration/training
This reaction matched rate rule [R3H_SS_O;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
family: intra_H_migration"""),
)

reaction(
    label = 'reaction3',
    reactants = ['npropylperoxy(106)'],
    products = ['QOOH_2(210)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4e+08,'s^-1'), n=1.1, Ea=(125.897,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Matched reaction 249 C3H7O2-7 <=> C3H7O2-8 in intra_H_migration/training
This reaction matched rate rule [R4H_SSS_OCs;O_rad_out;Cs_H_out_H/NonDeC]
family: intra_H_migration"""),
)

reaction(
    label = 'reaction4',
    reactants = ['npropylperoxy(106)'],
    products = ['QOOH_1(211)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.07e+06,'s^-1'), n=1.55, Ea=(87.9477,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Matched reaction 255 C3H7O2-9 <=> C3H7O2-10 in intra_H_migration/training
This reaction matched rate rule [R5H_SSSS_OCC;O_rad_out;Cs_H_out_2H]
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['npropylperoxy(106)'],
    products = ['HO2(10)', 'C3H6(59)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(9.58e+07,'s^-1'), n=1.46, Ea=(123.01,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Matched reaction 3 C3H7O2 <=> C3H6 + HO2 in HO2_Elimination_from_PeroxyRadical/training
This reaction matched rate rule [R2OO_HNd_2H]
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction12',
    reactants = ['npropoxy(213)', 'O(4)'],
    products = ['npropylperoxy(106)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(54738.4,'m^3/(mol*s)'), n=0.884925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(3)', 'CH2CHCH2OOH(252)'],
    products = ['QOOH_1(211)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.73e+08,'cm^3/(mol*s)'), n=1.583, Ea=(7.57304,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2832 used for Cds-Cs\O2s/H_Cds-HH;HJ
Exact match found for rate rule [Cds-Cs\O2s/H_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C2H4(29)', '[CH2]OO(239)'],
    products = ['QOOH_1(211)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.82e+10,'cm^3/(mol*s)','*|/',5), n=0, Ea=(29.1206,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 2750 used for Cds-HH_Cds-HH;CsJ-OsHH
Exact match found for rate rule [Cds-HH_Cds-HH;CsJ-OsHH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction15',
    reactants = ['QOOH_1(211)'],
    products = ['QOOH_2(210)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(718000,'s^-1'), n=2.05, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(500,'K'), Tmax=(2000,'K'), comment="""From training reaction 147 used for R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CC[CH]OO(209)'],
    products = ['QOOH_1(211)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(6e-15,'s^-1'), n=8.23, Ea=(142.674,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 343 used for R3H_SS_Cs;C_rad_out_H/NonDeO;Cs_H_out_2H
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeO;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['OH(5)', '[CH2]CC[O](253)'],
    products = ['QOOH_1(211)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 103 used for O_pri_rad;O_rad/NonDe
Exact match found for rate rule [O_pri_rad;O_rad/NonDe]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C[CH2](83)', 'HO2(10)'],
    products = ['QOOH_1(211)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(5.1779e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [O_rad;C_rad/H2/Cs] + [O_rad/NonDe;Cs_rad] for rate rule [O_rad/NonDe;C_rad/H2/Cs]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(3)', '[CH2]CCO[O](201)'],
    products = ['QOOH_1(211)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(5.21063e+06,'m^3/(mol*s)'), n=0.156446, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;O_rad/NonDe] + [H_rad;O_sec_rad] for rate rule [H_rad;O_rad/NonDe]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['C2H4(T)(79)', '[CH2]OO(239)'],
    products = ['QOOH_1(211)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.08115e+07,'m^3/(mol*s)'), n=0.104005, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_pri_rad;Y_rad] for rate rule [C_rad/H2/O;Y_rad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -0.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(3)', '[CH2][CH]COO(254)'],
    products = ['QOOH_1(211)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.34601e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', '[CH2]C[CH]OO(255)'],
    products = ['QOOH_1(211)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.18e+12,'cm^3/(mol*s)'), n=-0.085, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H/CsO;Y_rad] for rate rule [C_rad/H/CsO;H_rad]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -2.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['QOOH_1(211)'],
    products = ['OH(5)', 'oxetane(256)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.47e+11,'s^-1','*|/',1.74), n=0, Ea=(75.8559,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 7 used for R3OO_SS;C_pri_rad_intra;OOH
Exact match found for rate rule [R3OO_SS;C_pri_rad_intra;OOH]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['QOOH_1(211)'],
    products = ['[O]CCCO(257)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.47e+10,'s^-1'), n=0, Ea=(116.315,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Matched reaction 3 C3H7O2-3 <=> C3H7O2-4 in intra_OH_migration/training
This reaction matched rate rule [R3OOH_SS;C_rad_out_2H]
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CH2CH2OOH(55)', 'CH2(T)(17)'],
    products = ['QOOH_1(211)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.23625e+06,'m^3/(mol*s)'), n=0.36814, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(3)', '[CH]CCOO(247)'],
    products = ['QOOH_1(211)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '34',
    isomers = [
        'npropylperoxy(106)',
        'QOOH_1(211)',
    ],
    reactants = [
        ('O2(2)', 'npropyl(70)'),
        ('CH2(S)(24)', 'CH3CH2OO(54)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '34',
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

