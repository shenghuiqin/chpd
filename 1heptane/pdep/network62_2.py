species(
    label = 'CCCCO[O](164)',
    structure = SMILES('CCCCO[O]'),
    E0 = (-80.6218,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,232.129,232.181],'cm^-1')),
        HinderedRotor(inertia=(0.157672,'amu*angstrom^2'), symmetry=1, barrier=(6.03016,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157561,'amu*angstrom^2'), symmetry=1, barrier=(6.02839,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157484,'amu*angstrom^2'), symmetry=1, barrier=(6.02954,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157681,'amu*angstrom^2'), symmetry=1, barrier=(6.02962,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (89.1131,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3638.29,'J/mol'), sigma=(6.34008,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=568.29 K, Pc=32.39 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56797,0.0567642,-4.41695e-05,2.06267e-08,-4.35472e-12,-9611.76,22.5652], Tmin=(100,'K'), Tmax=(1060.64,'K')), NASAPolynomial(coeffs=[6.8104,0.0369938,-1.62098e-05,3.05291e-09,-2.12531e-13,-10723.8,-3.03526], Tmin=(1060.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-80.6218,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(ROOJ)"""),
)

species(
    label = 'C[CH]CCOO(267)',
    structure = SMILES('C[CH]CCOO'),
    E0 = (-38.1802,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,277.55,2321.54],'cm^-1')),
        HinderedRotor(inertia=(0.174609,'amu*angstrom^2'), symmetry=1, barrier=(9.54475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174588,'amu*angstrom^2'), symmetry=1, barrier=(9.54476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174621,'amu*angstrom^2'), symmetry=1, barrier=(9.54475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.882932,'amu*angstrom^2'), symmetry=1, barrier=(48.2632,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0126195,'amu*angstrom^2'), symmetry=1, barrier=(48.2635,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (89.1131,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95867,0.05283,-1.61638e-05,-5.25247e-08,5.56416e-11,-4526.53,23.266], Tmin=(100,'K'), Tmax=(485.214,'K')), NASAPolynomial(coeffs=[4.83957,0.0414031,-1.89328e-05,3.61992e-09,-2.53629e-13,-4951.15,9.95596], Tmin=(485.214,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-38.1802,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJC)"""),
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
    label = 'butyl_1(82)',
    structure = SMILES('[CH2]CCC'),
    E0 = (63.0573,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.0977402,'amu*angstrom^2'), symmetry=1, barrier=(2.24724,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0976865,'amu*angstrom^2'), symmetry=1, barrier=(2.246,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0977534,'amu*angstrom^2'), symmetry=1, barrier=(2.24754,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.1143,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.176,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.25388,0.0316763,2.89994e-06,-1.98049e-08,8.20503e-12,7652.64,17.2725], Tmin=(100,'K'), Tmax=(1050.57,'K')), NASAPolynomial(coeffs=[7.59591,0.0260842,-1.01719e-05,1.85189e-09,-1.28169e-13,5716.37,-12.6366], Tmin=(1050.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(63.0573,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), label="""butyl_1""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'npropylperoxy(166)',
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
    label = 'CH2(S)(23)',
    structure = SMILES('[CH2]'),
    E0 = (419.862,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1369.36,2789.41,2993.36],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.19195,-0.00230793,8.0509e-06,-6.60123e-09,1.95638e-12,50484.3,-0.754589], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.28556,0.00460255,-1.97412e-06,4.09548e-10,-3.34695e-14,50922.4,8.67684], Tmin=(1000,'K'), Tmax=(3000,'K'))], Tmin=(200,'K'), Tmax=(3000,'K'), E0=(419.862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[CH2]CO[O](98)',
    structure = SMILES('[CH2]CO[O]'),
    E0 = (188.234,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,492.5,1135,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.0676416,'amu*angstrom^2'), symmetry=1, barrier=(11.1524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000766328,'amu*angstrom^2'), symmetry=1, barrier=(8.70093,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.83949,0.0210173,-4.14404e-06,-8.74986e-09,4.63377e-12,22685.1,18.0793], Tmin=(100,'K'), Tmax=(1017.31,'K')), NASAPolynomial(coeffs=[7.90625,0.0110099,-4.00757e-06,7.40873e-10,-5.2804e-14,21141.1,-8.97353], Tmin=(1017.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.234,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CJCOOH) + radical(ROOJ)"""),
)

species(
    label = 'C2H5(29)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.24186,-0.00356905,4.82667e-05,-5.85401e-08,2.25805e-11,12969,4.44704], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.32196,0.0123931,-4.39681e-06,7.0352e-10,-4.18435e-14,12175.9,0.171104], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(107.874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""C2H5""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'npropyl(83)',
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
    label = '[CH2]O[O](92)',
    structure = SMILES('[CH2]O[O]'),
    E0 = (206.805,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2219.2,2221.31],'cm^-1')),
        HinderedRotor(inertia=(0.156598,'amu*angstrom^2'), symmetry=1, barrier=(14.2064,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.45042,0.0115293,-8.64543e-06,3.98163e-09,-7.73155e-13,24893.2,10.7919], Tmin=(100,'K'), Tmax=(1212.82,'K')), NASAPolynomial(coeffs=[5.07483,0.00617178,-2.01911e-06,3.39169e-10,-2.23136e-14,24499.1,2.64172], Tmin=(1212.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(206.805,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CsJOOH) + radical(ROOJ)"""),
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
    label = 'CC[CH]CO[O](269)',
    structure = SMILES('CC[CH]CO[O]'),
    E0 = (119.793,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,1890.9],'cm^-1')),
        HinderedRotor(inertia=(0.400568,'amu*angstrom^2'), symmetry=1, barrier=(9.20985,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.400643,'amu*angstrom^2'), symmetry=1, barrier=(9.21156,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.400217,'amu*angstrom^2'), symmetry=1, barrier=(9.20177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0570395,'amu*angstrom^2'), symmetry=1, barrier=(1.31145,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.1051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78459,0.0552902,-3.99949e-05,-3.3506e-09,1.8958e-11,14481.2,24.161], Tmin=(100,'K'), Tmax=(545.219,'K')), NASAPolynomial(coeffs=[5.82119,0.036187,-1.63571e-05,3.10671e-09,-2.16676e-13,13884.8,5.70241], Tmin=(545.219,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(119.793,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCJCOOH) + radical(ROOJ)"""),
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
    label = '[CH2]CCO[O](270)',
    structure = SMILES('[CH2]CCO[O]'),
    E0 = (156.366,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,492.5,1135,1000,3000,3100,440,815,1455,1000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.220275,'amu*angstrom^2'), symmetry=1, barrier=(5.06457,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220108,'amu*angstrom^2'), symmetry=1, barrier=(5.06072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.219679,'amu*angstrom^2'), symmetry=1, barrier=(5.05086,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0785,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.30676,0.0392487,-3.42473e-05,1.87089e-08,-4.51416e-12,18865.7,20.7232], Tmin=(100,'K'), Tmax=(960.527,'K')), NASAPolynomial(coeffs=[6.03959,0.023704,-9.97259e-06,1.86101e-09,-1.29188e-13,18148.6,2.86467], Tmin=(960.527,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(156.366,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJ) + radical(ROOJ)"""),
)

species(
    label = 'C[CH]CCO[O](271)',
    structure = SMILES('C[CH]CCO[O]'),
    E0 = (113.824,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,2117.29],'cm^-1')),
        HinderedRotor(inertia=(0.211012,'amu*angstrom^2'), symmetry=1, barrier=(4.85157,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.210796,'amu*angstrom^2'), symmetry=1, barrier=(4.84662,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.210571,'amu*angstrom^2'), symmetry=1, barrier=(4.84145,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.211001,'amu*angstrom^2'), symmetry=1, barrier=(4.85132,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.1051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5302,0.0598843,-7.45167e-05,6.46865e-08,-2.34812e-11,13773.5,25.4286], Tmin=(100,'K'), Tmax=(824.048,'K')), NASAPolynomial(coeffs=[3.73794,0.0388313,-1.73791e-05,3.2396e-09,-2.21491e-13,13760.6,17.3342], Tmin=(824.048,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(113.824,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(ROOJ)"""),
)

species(
    label = 'CCC[CH]O[O](257)',
    structure = SMILES('CCC[CH]O[O]'),
    E0 = (107.961,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.119013,'amu*angstrom^2'), symmetry=1, barrier=(2.73635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120007,'amu*angstrom^2'), symmetry=1, barrier=(2.75919,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121685,'amu*angstrom^2'), symmetry=1, barrier=(2.79779,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.481289,'amu*angstrom^2'), symmetry=1, barrier=(11.0658,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.1051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80219,0.0553016,-3.67985e-05,-1.19632e-08,2.59165e-11,13057,21.9547], Tmin=(100,'K'), Tmax=(531.462,'K')), NASAPolynomial(coeffs=[5.78891,0.0367206,-1.66008e-05,3.14923e-09,-2.19362e-13,12471.9,3.72319], Tmin=(531.462,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(107.961,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CCsJOOH)"""),
)

species(
    label = '[CH2]CCCO[O](272)',
    structure = SMILES('[CH2]CCCO[O]'),
    E0 = (124.625,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,492.5,1135,1000,3000,3100,440,815,1455,1000,211.819,1654.03],'cm^-1')),
        HinderedRotor(inertia=(0.262171,'amu*angstrom^2'), symmetry=1, barrier=(8.34708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00375289,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.261522,'amu*angstrom^2'), symmetry=1, barrier=(8.3433,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.262083,'amu*angstrom^2'), symmetry=1, barrier=(8.34519,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.1051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51474,0.0584856,-5.77493e-05,3.6598e-08,-1.02894e-11,15075,24.5481], Tmin=(100,'K'), Tmax=(833.732,'K')), NASAPolynomial(coeffs=[6.38667,0.0351122,-1.56986e-05,2.97451e-09,-2.07471e-13,14262.6,1.92965], Tmin=(833.732,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(124.625,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(ROOJ)"""),
)

species(
    label = 'CCC[CH]OO(265)',
    structure = SMILES('CCC[CH]OO'),
    E0 = (-44.044,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,315.337,4000],'cm^-1')),
        HinderedRotor(inertia=(0.558143,'amu*angstrom^2'), symmetry=1, barrier=(39.3843,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.503974,'amu*angstrom^2'), symmetry=1, barrier=(35.562,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139417,'amu*angstrom^2'), symmetry=1, barrier=(9.83764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139416,'amu*angstrom^2'), symmetry=1, barrier=(9.83764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.69531,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (89.1131,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35904,0.0624115,-5.35472e-05,2.72538e-08,-6.17117e-12,-5205.8,22.7109], Tmin=(100,'K'), Tmax=(1006.83,'K')), NASAPolynomial(coeffs=[7.45805,0.0381808,-1.74476e-05,3.3506e-09,-2.35878e-13,-6433.93,-6.7548], Tmin=(1006.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-44.044,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOOH)"""),
)

species(
    label = 'CC[CH]COO(266)',
    structure = SMILES('CC[CH]COO'),
    E0 = (-32.2116,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,255.57,3395],'cm^-1')),
        HinderedRotor(inertia=(0.262813,'amu*angstrom^2'), symmetry=1, barrier=(12.1744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00258186,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.26274,'amu*angstrom^2'), symmetry=1, barrier=(12.1742,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.826409,'amu*angstrom^2'), symmetry=1, barrier=(38.2899,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.826502,'amu*angstrom^2'), symmetry=1, barrier=(38.2899,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (89.1131,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39054,0.0615984,-5.23512e-05,2.61865e-08,-5.81263e-12,-3783.71,24.7523], Tmin=(100,'K'), Tmax=(1025.72,'K')), NASAPolynomial(coeffs=[7.56862,0.0375059,-1.71189e-05,3.28735e-09,-2.31433e-13,-5051.12,-5.21033], Tmin=(1025.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-32.2116,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCJCOOH)"""),
)

species(
    label = '[CH2]CCCOO(268)',
    structure = SMILES('[CH2]CCCOO'),
    E0 = (-27.3802,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,289.998,2078.78],'cm^-1')),
        HinderedRotor(inertia=(0.00200452,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.717646,'amu*angstrom^2'), symmetry=1, barrier=(42.8277,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161089,'amu*angstrom^2'), symmetry=1, barrier=(9.61352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161089,'amu*angstrom^2'), symmetry=1, barrier=(9.61355,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.405615,'amu*angstrom^2'), symmetry=1, barrier=(24.2066,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (89.1131,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38392,0.0606317,-4.84653e-05,2.17894e-08,-4.27174e-12,-3201.55,24.2497], Tmin=(100,'K'), Tmax=(1154.34,'K')), NASAPolynomial(coeffs=[8.58867,0.0356656,-1.60228e-05,3.05258e-09,-2.13776e-13,-4864.87,-11.543], Tmin=(1154.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-27.3802,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJ)"""),
)

species(
    label = 'butene1(127)',
    structure = SMILES('C=CCC'),
    E0 = (-16.4325,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,252.555],'cm^-1')),
        HinderedRotor(inertia=(0.178654,'amu*angstrom^2'), symmetry=1, barrier=(7.72883,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185589,'amu*angstrom^2'), symmetry=1, barrier=(7.72103,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.58773,0.0232778,1.93412e-05,-3.55496e-08,1.36906e-11,-1918.73,14.5751], Tmin=(100,'K'), Tmax=(1007.28,'K')), NASAPolynomial(coeffs=[7.20517,0.0236362,-9.0315e-06,1.65393e-09,-1.16019e-13,-3797.34,-12.4426], Tmin=(1007.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-16.4325,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), label="""butene1""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = 'C4H9O-5(273)',
    structure = SMILES('CCCC[O]'),
    E0 = (-76.2318,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,180,990.146,4000],'cm^-1')),
        HinderedRotor(inertia=(0.339916,'amu*angstrom^2'), symmetry=1, barrier=(7.81534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.930862,'amu*angstrom^2'), symmetry=1, barrier=(21.4023,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.647149,'amu*angstrom^2'), symmetry=1, barrier=(14.8792,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (73.1137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61175,0.0460187,-2.42198e-05,4.56076e-09,1.72087e-13,-9077.28,19.886], Tmin=(100,'K'), Tmax=(1312.71,'K')), NASAPolynomial(coeffs=[10.2963,0.0268501,-1.06514e-05,1.90298e-09,-1.27907e-13,-11985.8,-26.7686], Tmin=(1312.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-76.2318,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), label="""C4H9O-5""", comment="""Thermo library: CBS_QB3_1dHR"""),
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
    label = 'CC=CCOO(305)',
    structure = SMILES('CC=CCOO'),
    E0 = (-111.533,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,247.075],'cm^-1')),
        HinderedRotor(inertia=(0.18063,'amu*angstrom^2'), symmetry=1, barrier=(7.82442,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180636,'amu*angstrom^2'), symmetry=1, barrier=(7.82445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.316325,'amu*angstrom^2'), symmetry=1, barrier=(13.7021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.84087,'amu*angstrom^2'), symmetry=1, barrier=(36.4237,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.1051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80979,0.0510569,-3.39648e-05,1.18365e-08,-1.78624e-12,-13338.4,21.3409], Tmin=(100,'K'), Tmax=(1425.32,'K')), NASAPolynomial(coeffs=[8.40348,0.0325527,-1.44913e-05,2.72822e-09,-1.88684e-13,-15218,-12.8067], Tmin=(1425.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-111.533,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH)"""),
)

species(
    label = 'C=CCCOO(306)',
    structure = SMILES('C=CCCOO'),
    E0 = (-104.844,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,3010,987.5,1337.5,450,1655,323.996,323.996],'cm^-1')),
        HinderedRotor(inertia=(0.111179,'amu*angstrom^2'), symmetry=1, barrier=(8.28185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.331823,'amu*angstrom^2'), symmetry=1, barrier=(24.718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111178,'amu*angstrom^2'), symmetry=1, barrier=(8.28185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.419395,'amu*angstrom^2'), symmetry=1, barrier=(31.2413,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.1051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37317,0.05546,-4.0154e-05,1.50002e-08,-2.32081e-12,-12513.7,22.9715], Tmin=(100,'K'), Tmax=(1476.88,'K')), NASAPolynomial(coeffs=[11.7098,0.0274642,-1.17202e-05,2.16524e-09,-1.48175e-13,-15566.9,-30.9276], Tmin=(1476.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-104.844,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'CH2OOH(35)',
    structure = SMILES('[CH2]OO'),
    E0 = (52.1952,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(2.16183,'amu*angstrom^2'), symmetry=1, barrier=(49.7048,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (47.0333,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3467.13,'J/mol'), sigma=(3.69,'angstroms'), dipoleMoment=(1.7,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[5.83127,-0.00351771,4.54551e-05,-5.66903e-08,2.21633e-11,6061.87,-0.579143], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[6.98746,0.00900484,-3.24367e-06,5.24325e-10,-3.13587e-14,5012.58,-10.2619], Tmin=(1000,'K'), Tmax=(2500,'K'))], Tmin=(200,'K'), Tmax=(2500,'K'), E0=(52.1952,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), label="""CH2OOH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C3H6(72)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31912,0.00817959,3.34736e-05,-4.36194e-08,1.58213e-11,749.325,9.54025], Tmin=(100,'K'), Tmax=(983.754,'K')), NASAPolynomial(coeffs=[5.36755,0.0170743,-6.35108e-06,1.1662e-09,-8.2762e-14,-487.138,-4.54468], Tmin=(983.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.9763,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""C3H6""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C[CH]CC[O](307)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.44243,0.0298368,-6.59903e-06,-1.96417e-09,6.28703e-13,14825.4,18.071], Tmin=(100,'K'), Tmax=(2180.06,'K')), NASAPolynomial(coeffs=[16.562,0.0192354,-8.57322e-06,1.47387e-09,-9.0585e-14,5904.07,-62.7899], Tmin=(2180.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(123.235,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo library: CBS_QB3_1dHR + radical(CCOJ) + radical(RCCJC)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.4858,0.00133397,-4.70043e-06,5.64379e-09,-2.06318e-12,3411.96,1.99788], Tmin=(100,'K'), Tmax=(1005.25,'K')), NASAPolynomial(coeffs=[2.88225,0.00103869,-2.35652e-07,1.40229e-11,6.34581e-16,3669.56,5.59053], Tmin=(1005.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.372,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = '[CH2]C[CH]C(129)',
    structure = SMILES('[CH2]C[CH]C'),
    E0 = (255.669,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,2031.24],'cm^-1')),
        HinderedRotor(inertia=(0.244974,'amu*angstrom^2'), symmetry=1, barrier=(5.63244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00192352,'amu*angstrom^2'), symmetry=1, barrier=(5.63177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244928,'amu*angstrom^2'), symmetry=1, barrier=(5.63137,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.98997,0.0287412,-9.51469e-06,4.19232e-10,1.90526e-13,30780.1,16.8971], Tmin=(100,'K'), Tmax=(2154.56,'K')), NASAPolynomial(coeffs=[12.4231,0.0182241,-7.06316e-06,1.16769e-09,-7.11818e-14,25091.4,-39.6212], Tmin=(2154.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(255.669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'C3H6(T)(143)',
    structure = SMILES('[CH2][CH]C'),
    E0 = (284.865,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.238389,'amu*angstrom^2'), symmetry=1, barrier=(5.48103,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00909639,'amu*angstrom^2'), symmetry=1, barrier=(22.1005,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.93778,0.0190991,4.26842e-06,-1.44873e-08,5.74941e-12,34303.2,12.9695], Tmin=(100,'K'), Tmax=(1046.81,'K')), NASAPolynomial(coeffs=[5.93909,0.0171892,-6.69152e-06,1.21546e-09,-8.39795e-14,33151.2,-4.14888], Tmin=(1046.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(284.865,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""C3H6(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C[CH][CH]COO(308)',
    structure = SMILES('C[CH][CH]COO'),
    E0 = (162.235,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,3070.01],'cm^-1')),
        HinderedRotor(inertia=(1.82788,'amu*angstrom^2'), symmetry=1, barrier=(42.0267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00627932,'amu*angstrom^2'), symmetry=1, barrier=(42.0177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161571,'amu*angstrom^2'), symmetry=1, barrier=(3.71484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161691,'amu*angstrom^2'), symmetry=1, barrier=(3.7176,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.1051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29794,0.0654316,-8.55343e-05,7.44939e-08,-2.70526e-11,19603.8,27.8082], Tmin=(100,'K'), Tmax=(804.143,'K')), NASAPolynomial(coeffs=[4.60467,0.0391518,-1.81748e-05,3.44674e-09,-2.38096e-13,19389.9,14.5524], Tmin=(804.143,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(162.235,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCJCOOH) + radical(RCCJC)"""),
)

species(
    label = 'C[CH]C[CH]OO(309)',
    structure = SMILES('C[CH]C[CH]OO'),
    E0 = (150.402,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,215.104,1603.04],'cm^-1')),
        HinderedRotor(inertia=(0.168488,'amu*angstrom^2'), symmetry=1, barrier=(5.53204,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00303367,'amu*angstrom^2'), symmetry=1, barrier=(5.53202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.7442,'amu*angstrom^2'), symmetry=1, barrier=(57.2689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168485,'amu*angstrom^2'), symmetry=1, barrier=(5.53203,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.74446,'amu*angstrom^2'), symmetry=1, barrier=(57.2689,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.1051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26598,0.0662667,-8.68783e-05,7.58369e-08,-2.75572e-11,18181.7,25.7673], Tmin=(100,'K'), Tmax=(806.273,'K')), NASAPolynomial(coeffs=[4.55228,0.0397216,-1.84402e-05,3.49454e-09,-2.41229e-13,17984.7,12.6848], Tmin=(806.273,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(150.402,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(CCsJOOH)"""),
)

species(
    label = '[CH2][CH]CCOO(310)',
    structure = SMILES('[CH2][CH]CCOO'),
    E0 = (167.066,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,1842.41,1843.45],'cm^-1')),
        HinderedRotor(inertia=(0.204753,'amu*angstrom^2'), symmetry=1, barrier=(9.22363,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208519,'amu*angstrom^2'), symmetry=1, barrier=(9.26148,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00383821,'amu*angstrom^2'), symmetry=1, barrier=(9.25006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12134,'amu*angstrom^2'), symmetry=1, barrier=(49.1916,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07868,'amu*angstrom^2'), symmetry=1, barrier=(49.2764,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.1051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33193,0.063861,-7.90196e-05,6.60931e-08,-2.3602e-11,20184.5,27.1691], Tmin=(100,'K'), Tmax=(784.709,'K')), NASAPolynomial(coeffs=[5.06473,0.038298,-1.76625e-05,3.35254e-09,-2.32269e-13,19799.9,11.3475], Tmin=(784.709,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(167.066,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(RCCJC)"""),
)

species(
    label = 'CC1CCO1(311)',
    structure = SMILES('CC1CCO1'),
    E0 = (-132.743,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.60713,0.0163702,6.56459e-05,-9.99103e-08,4.19656e-11,-15901.6,14.375], Tmin=(100,'K'), Tmax=(884.962,'K')), NASAPolynomial(coeffs=[10.4297,0.0189449,-3.01349e-06,2.48205e-10,-1.18947e-14,-18771.5,-30.8008], Tmin=(884.962,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-132.743,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane)"""),
)

species(
    label = '[CH2]C(C)COO(312)',
    structure = SMILES('[CH2]C(C)COO'),
    E0 = (-36.5857,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,199.379],'cm^-1')),
        HinderedRotor(inertia=(0.457221,'amu*angstrom^2'), symmetry=1, barrier=(12.9746,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00423632,'amu*angstrom^2'), symmetry=1, barrier=(0.119751,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.458367,'amu*angstrom^2'), symmetry=1, barrier=(12.976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.64437,'amu*angstrom^2'), symmetry=1, barrier=(74.8234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0267402,'amu*angstrom^2'), symmetry=1, barrier=(74.8241,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (89.1131,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12376,0.0620874,-5.11548e-05,2.37099e-08,-4.6268e-12,-4295.76,24.8951], Tmin=(100,'K'), Tmax=(1197.69,'K')), NASAPolynomial(coeffs=[10.3239,0.031361,-1.26727e-05,2.28973e-09,-1.55656e-13,-6499.55,-21.1501], Tmin=(1197.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-36.5857,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl)"""),
)

species(
    label = 'CC(O)CC[O](313)',
    structure = SMILES('CC(O)CC[O]'),
    E0 = (-250.897,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (89.1131,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07665,0.0502249,-9.66436e-06,-6.36821e-08,6.41349e-11,-30114.7,22.5744], Tmin=(100,'K'), Tmax=(472.896,'K')), NASAPolynomial(coeffs=[4.62333,0.0408781,-1.86969e-05,3.58055e-09,-2.51192e-13,-30492,10.7536], Tmin=(472.896,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-250.897,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ)"""),
)

species(
    label = 'CH2CH2OOH(45)',
    structure = SMILES('[CH2]COO'),
    E0 = (32.3724,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.00650046,'amu*angstrom^2'), symmetry=1, barrier=(11.0723,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2783,'amu*angstrom^2'), symmetry=1, barrier=(29.3906,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.22246,'amu*angstrom^2'), symmetry=1, barrier=(74.0907,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (61.0599,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3284.22,'J/mol'), sigma=(4.037,'angstroms'), dipoleMoment=(1.3,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[5.07456,0.013239,2.53759e-05,-4.3809e-08,1.89062e-11,3710.63,2.7229], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[9.4138,0.0133542,-4.68068e-06,7.43231e-10,-4.39824e-14,1984.6,-22.4517], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(32.3724,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), label="""CH2CH2OOH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'CHCH3(T)(95)',
    structure = SMILES('[CH]C'),
    E0 = (343.893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,592.414,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00438699,'amu*angstrom^2'), symmetry=1, barrier=(26.7685,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.82363,-0.000909515,3.2138e-05,-3.7348e-08,1.3309e-11,41371.4,7.10948], Tmin=(100,'K'), Tmax=(960.812,'K')), NASAPolynomial(coeffs=[4.30487,0.00943069,-3.27559e-06,5.95121e-10,-4.27307e-14,40709.1,1.84202], Tmin=(960.812,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CHCH3(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]CCOO(314)',
    structure = SMILES('[CH]CCOO'),
    E0 = (247.329,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1310,387.5,850,1000,377.738,377.935,378.029,1475.62],'cm^-1')),
        HinderedRotor(inertia=(0.00117972,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00117919,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0639907,'amu*angstrom^2'), symmetry=1, barrier=(6.48327,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0638991,'amu*angstrom^2'), symmetry=1, barrier=(6.48421,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0785,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92275,0.045598,-4.03983e-05,1.93856e-08,-3.83587e-12,29821.7,20.5537], Tmin=(100,'K'), Tmax=(1193.06,'K')), NASAPolynomial(coeffs=[9.46509,0.0203106,-8.60528e-06,1.62007e-09,-1.13197e-13,28022,-17.1651], Tmin=(1193.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(247.329,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJ2_triplet)"""),
)

species(
    label = 'C[C]CCOO(315)',
    structure = SMILES('C[C]CCOO'),
    E0 = (215.589,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.1051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2489,0.063427,-5.89516e-05,3.08357e-08,-6.84935e-12,26026,22.858], Tmin=(100,'K'), Tmax=(1053.35,'K')), NASAPolynomial(coeffs=[9.43763,0.0323308,-1.46694e-05,2.80924e-09,-1.97541e-13,24300.9,-17.0734], Tmin=(1053.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(215.589,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
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
    E0 = (54.4355,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (296.108,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (293.867,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (331.585,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (292.639,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (325.617,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (319.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (336.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (361.127,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (85.3152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (51.8868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-4.59852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (3.26738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (42.3882,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (166.773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (105.992,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (109.458,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (95.2836,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (126.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (120.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (128.546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (151.607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (258.345,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (325.617,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (337.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (374.027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (362.194,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (378.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (37.6757,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (123.349,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (73.1138,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (376.265,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (383.517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (427.381,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O2(2)', 'butyl_1(82)'],
    products = ['CCCCO[O](164)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(8.21765e+07,'m^3/(mol*s)'), n=-0.529675, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for O2_birad;C_rad/H2/Cs
Exact match found for rate rule [C_rad/H2/Cs;O2_birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -1.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]CO[O](98)', 'C2H5(29)'],
    products = ['CCCCO[O](164)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.73e+14,'cm^3/(mol*s)'), n=-0.699, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 11 used for C_rad/H2/Cs;C_rad/H2/Cs
Exact match found for rate rule [C_rad/H2/Cs;C_rad/H2/Cs]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['npropyl(83)', '[CH2]O[O](92)'],
    products = ['CCCCO[O](164)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.33778e+08,'m^3/(mol*s)'), n=-0.3495, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H2/Cs;C_pri_rad] for rate rule [C_rad/H2/Cs;C_rad/H2/O]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', 'CC[CH]CO[O](269)'],
    products = ['CCCCO[O](164)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CH3(17)', '[CH2]CCO[O](270)'],
    products = ['CCCCO[O](164)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.23e+15,'cm^3/(mol*s)'), n=-0.562, Ea=(0.085772,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 10 used for C_methyl;C_rad/H2/Cs
Exact match found for rate rule [C_methyl;C_rad/H2/Cs]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', 'C[CH]CCO[O](271)'],
    products = ['CCCCO[O](164)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', 'CCC[CH]O[O](257)'],
    products = ['CCCCO[O](164)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.18e+12,'cm^3/(mol*s)'), n=-0.085, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H/CsO;Y_rad] for rate rule [C_rad/H/CsO;H_rad]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -2.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)', '[CH2]CCCO[O](272)'],
    products = ['CCCCO[O](164)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction9',
    reactants = ['npropylperoxy(166)', 'CH2(S)(23)'],
    products = ['CCCCO[O](164)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction2',
    reactants = ['CCCCO[O](164)'],
    products = ['CCC[CH]OO(265)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.44e+09,'s^-1'), n=1.17, Ea=(165.937,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 242 used for R3H_SS_O;O_rad_out;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R3H_SS_O;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CC[CH]COO(266)'],
    products = ['CCCCO[O](164)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2960,'s^-1'), n=2.11, Ea=(84.0984,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 323 used for R4H_SSS;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CCCCO[O](164)'],
    products = ['C[CH]CCOO(267)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.788e+07,'s^-1'), n=1.26, Ea=(76.0233,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Matched reaction 260 C4H9O2-13 <=> C4H9O2-14 in intra_H_migration/training
This reaction matched rate rule [R5H_SSSS_OCC;O_rad_out;Cs_H_out_H/NonDeC]
family: intra_H_migration"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CCCCO[O](164)'],
    products = ['[CH2]CCCOO(268)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.107e+06,'s^-1'), n=1.52, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Matched reaction 263 C4H9O2-15 <=> C4H9O2-16 in intra_H_migration/training
This reaction matched rate rule [R6H_SSSSS_OO;O_rad_out;Cs_H_out_2H]
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['CCCCO[O](164)'],
    products = ['butene1(127)', 'HO2(10)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(9.58e+07,'s^-1'), n=1.46, Ea=(123.01,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 3 used for R2OO_HNd_2H
Exact match found for rate rule [R2OO_HNd_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C4H9O-5(273)', 'O(4)'],
    products = ['CCCCO[O](164)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', 'CC=CCOO(305)'],
    products = ['C[CH]CCOO(267)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(3)', 'C=CCCOO(306)'],
    products = ['C[CH]CCOO(267)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.36e+08,'cm^3/(mol*s)'), n=1.56, Ea=(2.5104,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 9 used for Cds-HH_Cds-CsH;HJ
Exact match found for rate rule [Cds-HH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CH2OOH(35)', 'C3H6(72)'],
    products = ['C[CH]CCOO(267)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.381e+06,'cm^3/(mol*s)'), n=1.76, Ea=(37.1121,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2781 used for Cds-HH_Cds-Cs\H3/H;CsJ-OsHH
Exact match found for rate rule [Cds-HH_Cds-Cs\H3/H;CsJ-OsHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CC[CH]COO(266)'],
    products = ['C[CH]CCOO(267)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(6.76e+09,'s^-1'), n=0.88, Ea=(158.992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 357 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]CCCOO(268)'],
    products = ['C[CH]CCOO(267)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(6.48e+07,'s^-1'), n=1.57, Ea=(147.695,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 106 used for R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CCC[CH]OO(265)'],
    products = ['C[CH]CCOO(267)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.85927e+09,'s^-1'), n=0.875, Ea=(172.59,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_1H;Cs_H_out_H/NonDeC] for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeO;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C[CH]CC[O](307)', 'OH(5)'],
    products = ['C[CH]CCOO(267)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 103 used for O_pri_rad;O_rad/NonDe
Exact match found for rate rule [O_pri_rad;O_rad/NonDe]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C[CH]C(129)', 'HO2(10)'],
    products = ['C[CH]CCOO(267)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.58895e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [O_rad;C_rad/H2/Cs] + [O_rad/NonDe;Cs_rad] for rate rule [O_rad/NonDe;C_rad/H2/Cs]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', 'C[CH]CCO[O](271)'],
    products = ['C[CH]CCOO(267)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(5.21063e+06,'m^3/(mol*s)'), n=0.156446, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;O_rad/NonDe] + [H_rad;O_sec_rad] for rate rule [H_rad;O_rad/NonDe]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['C3H6(T)(143)', 'CH2OOH(35)'],
    products = ['C[CH]CCOO(267)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.68901e+06,'m^3/(mol*s)'), n=0.124131, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_pri_rad;Y_rad] for rate rule [C_rad/H2/O;Y_rad]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(3)', 'C[CH][CH]COO(308)'],
    products = ['C[CH]CCOO(267)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.34601e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(3)', 'C[CH]C[CH]OO(309)'],
    products = ['C[CH]CCOO(267)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.18e+12,'cm^3/(mol*s)'), n=-0.085, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H/CsO;Y_rad] for rate rule [C_rad/H/CsO;H_rad]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -2.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['H(3)', '[CH2][CH]CCOO(310)'],
    products = ['C[CH]CCOO(267)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(4.34601e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction29',
    reactants = ['C[CH]CCOO(267)'],
    products = ['OH(5)', 'CC1CCO1(311)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.04e+11,'s^-1','*|/',1.74), n=0, Ea=(75.8559,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 9 used for R3OO_SS;C_rad/H/NonDeC_intra;OOH
Exact match found for rate rule [R3OO_SS;C_rad/H/NonDeC_intra;OOH]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(C)COO(312)'],
    products = ['C[CH]CCOO(267)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C[CH]CCOO(267)'],
    products = ['CC(O)CC[O](313)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.88e+10,'s^-1'), n=0, Ea=(111.294,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Matched reaction 4 C4H9O2 <=> C4H9O2-2 in intra_OH_migration/training
This reaction matched rate rule [R3OOH_SS;C_rad_out_H/NonDeC]
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['CH2CH2OOH(45)', 'CHCH3(T)(95)'],
    products = ['C[CH]CCOO(267)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction33',
    reactants = ['CH3(17)', '[CH]CCOO(314)'],
    products = ['C[CH]CCOO(267)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction34',
    reactants = ['H(3)', 'C[C]CCOO(315)'],
    products = ['C[CH]CCOO(267)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '62',
    isomers = [
        'CCCCO[O](164)',
        'C[CH]CCOO(267)',
    ],
    reactants = [
        ('O2(2)', 'butyl_1(82)'),
        ('npropylperoxy(166)', 'CH2(S)(23)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '62',
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

