species(
    label = 'CH3C(O)OO(60)',
    structure = SMILES('CC(=O)O[O]'),
    E0 = (-177.162,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,957.025,958.022],'cm^-1')),
        HinderedRotor(inertia=(0.113086,'amu*angstrom^2'), symmetry=1, barrier=(2.60007,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113755,'amu*angstrom^2'), symmetry=1, barrier=(2.61545,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (75.0434,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3620.95,'J/mol'), sigma=(4.86,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.12838,0.0294919,-2.00221e-05,5.32494e-09,3.03452e-14,-21203.4,18.0452], Tmin=(298,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.76604,0.0188598,-9.95181e-06,2.52994e-09,-2.50541e-13,-22126.8,-0.48408], Tmin=(1000,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-177.162,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""CH3C(O)OO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'HO2(10)',
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
    label = 'CH2CO(28)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13241,0.0181319,-1.74093e-05,9.35336e-09,-2.01725e-12,-7148.09,13.3808], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.75871,0.00635124,-2.25955e-06,3.62322e-10,-2.15856e-14,-8085.33,-4.9649], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-60.8183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CH2CO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'CH3C(O)O(68)',
    structure = SMILES('CC([O])=O'),
    E0 = (-205.363,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,180,645.948,645.955,645.966,2494.63],'cm^-1')),
        HinderedRotor(inertia=(0.00694592,'amu*angstrom^2'), symmetry=1, barrier=(2.05657,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (59.044,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3761.4,'J/mol'), sigma=(5.64299,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=587.52 K, Pc=47.5 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30016,0.0237699,-1.23328e-05,5.87355e-11,1.44737e-12,-24433.8,20.8502], Tmin=(298,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.15457,0.0168675,-8.75191e-06,2.18661e-09,-2.13394e-13,-25230.4,5.95053], Tmin=(1000,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-205.363,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CH3C(O)O""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,29226.7,5.11107], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,29226.7,5.11107], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.005,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'CC1([O])OO1(121)',
    structure = SMILES('CC1([O])OO1'),
    E0 = (-17.2638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (75.0434,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.85839,0.0202092,6.29314e-06,-2.44175e-08,1.209e-11,-2030.6,15.3442], Tmin=(100,'K'), Tmax=(880.153,'K')), NASAPolynomial(coeffs=[7.75001,0.0117955,-2.91569e-06,3.93819e-10,-2.36847e-14,-3426.86,-10.6709], Tmin=(880.153,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-17.2638,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsOsOsOs) + group(Cs-CsHHH) + ring(dioxirane) + radical(CCOJ)"""),
)

species(
    label = 'hoo_vinoxy(122)',
    structure = SMILES('[CH2]C(=O)OO'),
    E0 = (-164.177,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (75.0434,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.15507,0.0316121,-7.88167e-06,-1.61009e-08,8.92097e-12,-19671.6,15.8772], Tmin=(100,'K'), Tmax=(1013.67,'K')), NASAPolynomial(coeffs=[12.9803,0.00893269,-3.97261e-06,8.29147e-10,-6.39522e-14,-22895.8,-41.5734], Tmin=(1013.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-164.177,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), label="""hoo_vinoxy""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CO3t2(123)',
    structure = SMILES('[O]O[C]=O'),
    E0 = (106.943,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,262.32,262.473],'cm^-1')),
        HinderedRotor(inertia=(0.00244921,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.0089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.24172,0.0164884,-1.87036e-05,1.04487e-08,-2.29623e-12,12889.7,11.1068], Tmin=(100,'K'), Tmax=(1108.89,'K')), NASAPolynomial(coeffs=[6.67163,0.00411578,-1.96698e-06,3.86453e-10,-2.76577e-14,12129,-5.79498], Tmin=(1108.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(106.943,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(78.9875,'J/(mol*K)'), label="""CO3t2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'CH3(17)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.65718,0.0021266,5.45839e-06,-6.6181e-09,2.46571e-12,16422.7,1.67354], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.97812,0.00579785,-1.97558e-06,3.07298e-10,-1.79174e-14,16509.5,4.72248], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(136.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH3""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,25472.7,-0.459566], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,25472.7,-0.459566], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.792,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = '[CH2]C(=O)O[O](124)',
    structure = SMILES('[CH2]C(=O)O[O]'),
    E0 = (38.1404,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3000,3100,440,815,1455,1000,368.627,370.335,1670.47,1670.97],'cm^-1')),
        HinderedRotor(inertia=(0.0855604,'amu*angstrom^2'), symmetry=1, barrier=(8.42829,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00428116,'amu*angstrom^2'), symmetry=1, barrier=(8.44955,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0355,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.744,0.0267546,-2.56506e-05,1.35178e-08,-2.87591e-12,4633.19,16.5352], Tmin=(100,'K'), Tmax=(1134.57,'K')), NASAPolynomial(coeffs=[7.43341,0.0102217,-3.79238e-06,6.73934e-10,-4.57706e-14,3569.11,-6.6805], Tmin=(1134.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(38.1404,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(C(=O)OOJ) + radical(CJCO)"""),
)

species(
    label = 'C[C]1OOO1(125)',
    structure = SMILES('C[C]1OOO1'),
    E0 = (183.868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (75.0434,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.04173,0.00219293,7.91671e-05,-1.07874e-07,4.17847e-11,22165.8,16.4434], Tmin=(100,'K'), Tmax=(959.262,'K')), NASAPolynomial(coeffs=[13.7121,0.00574212,-1.50793e-06,4.03361e-10,-4.12563e-14,17908.3,-46.1125], Tmin=(959.262,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(183.868,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Cs_P)"""),
)

species(
    label = 'C=C(O)O[O](126)',
    structure = SMILES('C=C(O)O[O]'),
    E0 = (-47.2805,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,492.5,1135,1000,3615,1277.5,1000],'cm^-1')),
        HinderedRotor(inertia=(0.495194,'amu*angstrom^2'), symmetry=1, barrier=(11.3855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.493002,'amu*angstrom^2'), symmetry=1, barrier=(11.3351,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (75.0434,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98587,0.0461222,-6.82074e-05,5.20894e-08,-1.54506e-11,-5615.74,17.0372], Tmin=(100,'K'), Tmax=(895.779,'K')), NASAPolynomial(coeffs=[9.12923,0.0111284,-4.42523e-06,7.62346e-10,-4.91255e-14,-6771.3,-15.9459], Tmin=(895.779,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-47.2805,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(ROOJ)"""),
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
    E0 = (-34.4874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-17.2638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-80.093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (243.13,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (249.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (183.868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (156.898,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (37.642,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction6',
    reactants = ['CH3C(O)OO(60)'],
    products = ['HO2(10)', 'CH2CO(28)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(2.6e+09,'s^-1','*|/',2.51189), n=1.2, Ea=(34.1,'kcal/mol'), T0=(1,'K'), comment="""Matched reaction 0 C2H3O3 <=> C2H2O + HO2 in HO2_Elimination_from_PeroxyRadical/training
This reaction matched rate rule [R2OO_2H]
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction1',
    reactants = ['CH3C(O)OO(60)'],
    products = ['CC1([O])OO1(121)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(7.785e+11,'s^-1'), n=0.342, Ea=(159.898,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonylbond_intra_Nd;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 159.3 to 159.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH3C(O)OO(60)'],
    products = ['hoo_vinoxy(122)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.3e+09,'s^-1','*|/',2.51189), n=0.75, Ea=(23.2,'kcal/mol'), T0=(1,'K'), comment="""Matched reaction 0 C2H3O3 <=> C2H3O3-2 in intra_H_migration/training
This reaction matched rate rule [R4H_SSS;O_rad_out;Cs_H_out_2H]
family: intra_H_migration"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CO3t2(123)', 'CH3(17)'],
    products = ['CH3C(O)OO(60)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.2e+13,'cm^3/(mol*s)','+|-',8.4e+12), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 72 used for C_methyl;CO_rad/NonDe
Exact match found for rate rule [C_methyl;CO_rad/NonDe]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[CH2]C(=O)O[O](124)'],
    products = ['CH3C(O)OO(60)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(8.83648e+06,'m^3/(mol*s)'), n=0.3308, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_pri_rad;H_rad] for rate rule [C_rad/H2/CO;H_rad]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -1.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH3C(O)OO(60)'],
    products = ['C[C]1OOO1(125)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(361.03,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 358.4 to 361.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=C(O)O[O](126)'],
    products = ['CH3C(O)OO(60)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(37989.5,'s^-1'), n=2.515, Ea=(204.179,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH3C(O)O(68)', 'O(4)'],
    products = ['CH3C(O)OO(60)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2711.41,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [O_sec_rad;O_birad] for rate rule [O_rad/OneDe;O_birad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

network(
    label = '219',
    isomers = [
        'CH3C(O)OO(60)',
    ],
    reactants = [
        ('HO2(10)', 'CH2CO(28)'),
        ('CH3C(O)O(68)', 'O(4)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '219',
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

