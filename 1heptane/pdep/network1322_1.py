species(
    label = '[CH2]OC([CH2])(O)C(=O)CC(5498)',
    structure = SMILES('[CH2]OC([CH2])(O)C(=O)CC'),
    E0 = (-247.331,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,3615,1277.5,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,470.296,666.997,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156542,'amu*angstrom^2'), symmetry=1, barrier=(3.59922,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156542,'amu*angstrom^2'), symmetry=1, barrier=(3.59922,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156542,'amu*angstrom^2'), symmetry=1, barrier=(3.59922,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156542,'amu*angstrom^2'), symmetry=1, barrier=(3.59922,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156542,'amu*angstrom^2'), symmetry=1, barrier=(3.59922,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156542,'amu*angstrom^2'), symmetry=1, barrier=(3.59922,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156542,'amu*angstrom^2'), symmetry=1, barrier=(3.59922,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.59674,0.116279,-0.000132833,7.55921e-08,-1.6688e-11,-29539.6,37.3463], Tmin=(100,'K'), Tmax=(1116.22,'K')), NASAPolynomial(coeffs=[24.1455,0.0240295,-8.86446e-06,1.55052e-09,-1.0467e-13,-35286.3,-89.6749], Tmin=(1116.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-247.331,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-OdCsCs) + radical(CJC(C)OC) + radical(CsJOCH3)"""),
)

species(
    label = 'CH2O(15)',
    structure = SMILES('C=O'),
    E0 = (-119.055,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.79372,-0.00990833,3.7322e-05,-3.79285e-08,1.31773e-11,-14379.2,0.602798], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.16953,0.00619321,-2.25056e-06,3.65976e-10,-2.20149e-14,-14548.7,6.04208], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-119.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH2O""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C=C(O)C(=O)CC(4626)',
    structure = SMILES('C=C(O)C(=O)CC'),
    E0 = (-358.468,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,332.629,332.633],'cm^-1')),
        HinderedRotor(inertia=(0.152178,'amu*angstrom^2'), symmetry=1, barrier=(11.9524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152152,'amu*angstrom^2'), symmetry=1, barrier=(11.9478,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152157,'amu*angstrom^2'), symmetry=1, barrier=(11.9534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152209,'amu*angstrom^2'), symmetry=1, barrier=(11.954,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4071.48,'J/mol'), sigma=(6.49965,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=635.96 K, Pc=33.65 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.785843,0.0702816,-6.54766e-05,3.23543e-08,-6.53762e-12,-42997.7,22.9463], Tmin=(100,'K'), Tmax=(1175.98,'K')), NASAPolynomial(coeffs=[12.9689,0.0288414,-1.26178e-05,2.38822e-09,-1.6711e-13,-45863,-37.8049], Tmin=(1175.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-358.468,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2]C(O)=C([O])CC(4557)',
    structure = SMILES('[CH2]C(O)=C([O])CC'),
    E0 = (-184.943,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,402.631,402.631],'cm^-1')),
        HinderedRotor(inertia=(0.100035,'amu*angstrom^2'), symmetry=1, barrier=(11.5078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100034,'amu*angstrom^2'), symmetry=1, barrier=(11.5078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100035,'amu*angstrom^2'), symmetry=1, barrier=(11.5078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100035,'amu*angstrom^2'), symmetry=1, barrier=(11.5078,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4303.36,'J/mol'), sigma=(7.05412,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=672.18 K, Pc=27.82 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.516209,0.0797578,-7.96062e-05,3.90546e-08,-7.19184e-12,-22064.1,30.4076], Tmin=(100,'K'), Tmax=(1528.49,'K')), NASAPolynomial(coeffs=[20.9225,0.0119894,-1.65417e-06,6.23996e-11,2.30461e-15,-27255.3,-77.6604], Tmin=(1528.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-184.943,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(O)CJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]OC1(O)CC1([O])CC(5883)',
    structure = SMILES('[CH2]OC1(O)CC1([O])CC'),
    E0 = (-143.047,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.08479,0.0881781,-2.75525e-05,-5.01522e-08,3.19808e-11,-16999.7,33.4602], Tmin=(100,'K'), Tmax=(933.002,'K')), NASAPolynomial(coeffs=[30.8273,0.0115526,-1.12765e-06,1.09981e-10,-1.42967e-14,-25574.2,-132.323], Tmin=(933.002,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-143.047,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + ring(Cyclopropane) + radical(CsJOCH3) + radical(CC(C)2OJ)"""),
)

species(
    label = '[CH2]C1(O)OCC1([O])CC(5884)',
    structure = SMILES('[CH2]C1(O)OCC1([O])CC'),
    E0 = (-143.498,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.12862,0.115884,-0.000119024,5.90593e-08,-1.06886e-11,-16966.9,38.5401], Tmin=(100,'K'), Tmax=(1644.55,'K')), NASAPolynomial(coeffs=[27.3545,0.0124138,2.1003e-06,-8.85386e-10,7.24245e-14,-23027.3,-111.63], Tmin=(1644.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-143.498,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(440.667,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CJC(C)OC) + radical(CC(C)2OJ)"""),
)

species(
    label = 'C2H5CO(71)',
    structure = SMILES('CC[C]=O'),
    E0 = (-44.1874,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(0.192008,'amu*angstrom^2'), symmetry=1, barrier=(4.41465,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191563,'amu*angstrom^2'), symmetry=1, barrier=(4.40441,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0712,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.91313,0.0225012,-8.59328e-06,4.80319e-10,2.48743e-13,-5274.24,14.655], Tmin=(100,'K'), Tmax=(1673.18,'K')), NASAPolynomial(coeffs=[7.46306,0.0158512,-6.42141e-06,1.12499e-09,-7.32055e-14,-7388.54,-11.406], Tmin=(1673.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-44.1874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""CH3CH2CO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]OC(=C)O(3587)',
    structure = SMILES('[CH2]OC(=C)O'),
    E0 = (-113.088,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3615,1277.5,1000,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,385.624,385.986],'cm^-1')),
        HinderedRotor(inertia=(0.149807,'amu*angstrom^2'), symmetry=1, barrier=(15.8087,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149921,'amu*angstrom^2'), symmetry=1, barrier=(15.8066,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149522,'amu*angstrom^2'), symmetry=1, barrier=(15.8058,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (73.0706,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.930845,0.0518764,-5.23343e-05,2.50764e-08,-4.45103e-12,-13477.5,20.1549], Tmin=(100,'K'), Tmax=(1612.24,'K')), NASAPolynomial(coeffs=[16.0622,0.00502412,-8.08492e-08,-1.12707e-10,1.03568e-14,-17146.4,-56.3195], Tmin=(1612.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-113.088,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-OsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(C=COCJ)"""),
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
    label = '[CH2]OC(=C)C(=O)CC(5885)',
    structure = SMILES('[CH2]OC(=C)C(=O)CC'),
    E0 = (-147.776,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.585343,0.0775375,-6.83954e-05,3.24397e-08,-6.44822e-12,-17652.5,27.788], Tmin=(100,'K'), Tmax=(1171.43,'K')), NASAPolynomial(coeffs=[12.2864,0.0375831,-1.7235e-05,3.32439e-09,-2.3466e-13,-20393.9,-30.5145], Tmin=(1171.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-147.776,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsHH) + radical(C=COCJ)"""),
)

species(
    label = '[CH2]OC(C)([O])C(=O)CC(5886)',
    structure = SMILES('[CH2]OC(C)([O])C(=O)CC'),
    E0 = (-214.108,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3100,440,815,1455,1000,375,552.5,462.5,1710,180,180,180,180,1086.53,1129.13,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.163454,'amu*angstrom^2'), symmetry=1, barrier=(3.75813,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163454,'amu*angstrom^2'), symmetry=1, barrier=(3.75813,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163454,'amu*angstrom^2'), symmetry=1, barrier=(3.75813,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163454,'amu*angstrom^2'), symmetry=1, barrier=(3.75813,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163454,'amu*angstrom^2'), symmetry=1, barrier=(3.75813,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163454,'amu*angstrom^2'), symmetry=1, barrier=(3.75813,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.27049,0.103665,-0.00010355,5.1973e-08,-1.01704e-11,-25551.1,37.2071], Tmin=(100,'K'), Tmax=(1253.48,'K')), NASAPolynomial(coeffs=[23.5953,0.0243155,-8.59487e-06,1.47086e-09,-9.80052e-14,-31784.9,-88.3738], Tmin=(1253.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-214.108,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-OdCsCs) + radical(CsJOCH3) + radical(C=OCOJ)"""),
)

species(
    label = '[CH2]C([O])(OC)C(=O)CC(5887)',
    structure = SMILES('[CH2]C([O])(OC)C(=O)CC'),
    E0 = (-191.469,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3100,440,815,1455,1000,375,552.5,462.5,1710,180,180,180,180,1054.29,1166.28,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.165852,'amu*angstrom^2'), symmetry=1, barrier=(3.81326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165852,'amu*angstrom^2'), symmetry=1, barrier=(3.81326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165852,'amu*angstrom^2'), symmetry=1, barrier=(3.81326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165852,'amu*angstrom^2'), symmetry=1, barrier=(3.81326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165852,'amu*angstrom^2'), symmetry=1, barrier=(3.81326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165852,'amu*angstrom^2'), symmetry=1, barrier=(3.81326,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.08178,0.109839,-0.000124279,7.27241e-08,-1.68069e-11,-22843.6,35.492], Tmin=(100,'K'), Tmax=(1057.83,'K')), NASAPolynomial(coeffs=[19.7681,0.0309981,-1.24814e-05,2.26632e-09,-1.55239e-13,-27254.7,-66.2685], Tmin=(1057.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-191.469,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-OdCsCs) + radical(CJC(C)OC) + radical(C=OCOJ)"""),
)

species(
    label = '[CH2]OC(C)(O)C(=O)[CH]C(5888)',
    structure = SMILES('[CH2]OC(C)(O)C([O])=CC'),
    E0 = (-306.457,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.1431,0.111811,-0.000116045,5.80873e-08,-1.09555e-11,-36617.1,39.1468], Tmin=(100,'K'), Tmax=(1452.31,'K')), NASAPolynomial(coeffs=[29.0454,0.0138526,-2.41526e-06,2.0993e-10,-8.38331e-15,-44404.5,-118.58], Tmin=(1452.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-306.457,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CsJOCH3) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C(O)(OC)C(=O)[CH]C(5889)',
    structure = SMILES('[CH2]C(O)(OC)C([O])=CC'),
    E0 = (-283.817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.57128,0.113619,-0.000122129,6.05896e-08,-1.008e-11,-33926.7,36.0461], Tmin=(100,'K'), Tmax=(953.156,'K')), NASAPolynomial(coeffs=[25.2076,0.020665,-6.41565e-06,1.03819e-09,-6.86689e-14,-39914,-96.4908], Tmin=(953.156,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-283.817,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]CC(=O)C(C)(O)O[CH2](5890)',
    structure = SMILES('[CH2]CC(=O)C(C)(O)O[CH2]'),
    E0 = (-246.252,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,3615,1277.5,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,411.775,720.695,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15519,'amu*angstrom^2'), symmetry=1, barrier=(3.56813,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15519,'amu*angstrom^2'), symmetry=1, barrier=(3.56813,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15519,'amu*angstrom^2'), symmetry=1, barrier=(3.56813,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15519,'amu*angstrom^2'), symmetry=1, barrier=(3.56813,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15519,'amu*angstrom^2'), symmetry=1, barrier=(3.56813,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15519,'amu*angstrom^2'), symmetry=1, barrier=(3.56813,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15519,'amu*angstrom^2'), symmetry=1, barrier=(3.56813,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.64862,0.111631,-0.000119131,6.2744e-08,-1.27598e-11,-29403.2,38.4157], Tmin=(100,'K'), Tmax=(1213.58,'K')), NASAPolynomial(coeffs=[25.858,0.0209686,-7.07069e-06,1.18481e-09,-7.84671e-14,-36079.5,-99.6125], Tmin=(1213.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-246.252,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-OdCsCs) + radical(CJCC=O) + radical(CsJOCH3)"""),
)

species(
    label = '[CH2]CC(=O)C([CH2])(O)OC(5891)',
    structure = SMILES('[CH2]CC(=O)C([CH2])(O)OC'),
    E0 = (-223.613,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,3615,1277.5,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,508.983,631.523,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157433,'amu*angstrom^2'), symmetry=1, barrier=(3.61969,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157433,'amu*angstrom^2'), symmetry=1, barrier=(3.61969,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157433,'amu*angstrom^2'), symmetry=1, barrier=(3.61969,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157433,'amu*angstrom^2'), symmetry=1, barrier=(3.61969,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157433,'amu*angstrom^2'), symmetry=1, barrier=(3.61969,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157433,'amu*angstrom^2'), symmetry=1, barrier=(3.61969,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157433,'amu*angstrom^2'), symmetry=1, barrier=(3.61969,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.4958,0.118208,-0.000141169,8.50422e-08,-1.99897e-11,-26694.1,36.8307], Tmin=(100,'K'), Tmax=(1046.41,'K')), NASAPolynomial(coeffs=[22.3118,0.0271995,-1.07072e-05,1.92298e-09,-1.31059e-13,-31676.5,-79.1065], Tmin=(1046.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-223.613,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-OdCsCs) + radical(CJC(C)OC) + radical(CJCC=O)"""),
)

species(
    label = '[CH2][O](167)',
    structure = SMILES('[CH2][O]'),
    E0 = (192.903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88409,-0.00363885,3.28543e-05,-4.13611e-08,1.59631e-11,23210.8,7.47983], Tmin=(100,'K'), Tmax=(933.06,'K')), NASAPolynomial(coeffs=[6.69335,0.000289989,8.61416e-07,-1.56351e-10,7.33778e-15,21991.3,-9.6043], Tmin=(933.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(H3COJ) + radical(CsJOH)"""),
)

species(
    label = '[CH2]OC1(O)CO[C]1CC(5279)',
    structure = SMILES('[CH2]OC1(O)CO[C]1CC'),
    E0 = (-151.84,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,3150,900,1100,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4301.97,'J/mol'), sigma=(7.41337,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=671.96 K, Pc=23.96 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.16168,0.106278,-0.000105889,5.27136e-08,-9.76685e-12,-18015.1,37.9977], Tmin=(100,'K'), Tmax=(1564.27,'K')), NASAPolynomial(coeffs=[24.5133,0.0176282,-1.28228e-06,-2.20804e-10,2.80208e-14,-23859.9,-94.6355], Tmin=(1564.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-151.84,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + ring(Oxetane) + radical(C2CsJOCs) + radical(CsJOCH3)"""),
)

species(
    label = '[CH2]C1(O)OCO[C]1CC(5892)',
    structure = SMILES('[CH2]C1(O)OCO[C]1CC'),
    E0 = (-243.579,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.98387,0.126837,-0.000143089,7.34146e-08,-1.31364e-11,-28910.7,43.4122], Tmin=(100,'K'), Tmax=(1755.97,'K')), NASAPolynomial(coeffs=[27.3684,-0.000664578,1.17886e-05,-2.83582e-09,2.0386e-13,-31977.4,-107.264], Tmin=(1755.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-243.579,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(440.667,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-OsOsHH) + ring(1,3-Dioxolane) + radical(C2CsJOCs) + radical(CJC(C)OC)"""),
)

species(
    label = 'CH3CH2OH(54)',
    structure = SMILES('CCO'),
    E0 = (-248.759,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
        HinderedRotor(inertia=(0.00248522,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168417,'amu*angstrom^2'), symmetry=1, barrier=(8.09857,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (46.0684,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3014.83,'J/mol'), sigma=(4.53,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.8587,-0.00374017,6.95554e-05,-8.86548e-08,3.51688e-11,-29996.1,4.80185], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[6.56244,0.0152042,-5.38968e-06,8.6225e-10,-5.12898e-14,-31525.6,-9.47302], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-248.759,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""CH3CH2OH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[CH2]OC([CH2])=C=O(5893)',
    structure = SMILES('[CH2]OC(=C)[C]=O'),
    E0 = (83.4173,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,1855,455,950,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.792728,'amu*angstrom^2'), symmetry=1, barrier=(18.2264,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.792013,'amu*angstrom^2'), symmetry=1, barrier=(18.2099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.792413,'amu*angstrom^2'), symmetry=1, barrier=(18.2191,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22613,0.059634,-6.86022e-05,3.80885e-08,-8.23875e-12,10133.9,20.7083], Tmin=(100,'K'), Tmax=(1131.19,'K')), NASAPolynomial(coeffs=[14.5722,0.0124407,-6.0219e-06,1.20671e-09,-8.76019e-14,7114.51,-45.3241], Tmin=(1131.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(83.4173,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-OsHHH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJ=O) + radical(C=COCJ)"""),
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
    label = '[CH2]OC([CH2])(O)C(C)=O(5894)',
    structure = SMILES('[CH2]OC([CH2])(O)C(C)=O'),
    E0 = (-221.914,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,3615,1277.5,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,180,1600,1607.8,2880.97,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15024,'amu*angstrom^2'), symmetry=1, barrier=(3.45431,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15024,'amu*angstrom^2'), symmetry=1, barrier=(3.45431,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15024,'amu*angstrom^2'), symmetry=1, barrier=(3.45431,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15024,'amu*angstrom^2'), symmetry=1, barrier=(3.45431,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15024,'amu*angstrom^2'), symmetry=1, barrier=(3.45431,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15024,'amu*angstrom^2'), symmetry=1, barrier=(3.45431,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.19317,0.0984671,-0.000108685,5.7434e-08,-1.1478e-11,-26489.7,34.5695], Tmin=(100,'K'), Tmax=(1331.37,'K')), NASAPolynomial(coeffs=[25.7077,0.010793,-2.18576e-06,2.39712e-10,-1.22828e-14,-33045.3,-100.63], Tmin=(1331.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-221.914,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cs-OsHHH) + group(Cds-OdCsCs) + radical(CJC(C)OC) + radical(CsJOCH3)"""),
)

species(
    label = '[CH2]O[C](O)CC(=O)CC(5895)',
    structure = SMILES('[CH2]O[C](O)CC(=O)CC'),
    E0 = (-229.806,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.627294,0.112633,-0.000169589,1.4995e-07,-5.31496e-11,-27483.2,37.5988], Tmin=(100,'K'), Tmax=(814.131,'K')), NASAPolynomial(coeffs=[8.07234,0.0520761,-2.51937e-05,4.83335e-09,-3.34635e-13,-28309.3,1.04273], Tmin=(814.131,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-229.806,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-OdCsCs) + radical(Cs_P) + radical(CsJOCH3)"""),
)

species(
    label = 'CCC(=O)C1(O)CCO1(5870)',
    structure = SMILES('CCC(=O)C1(O)CCO1'),
    E0 = (-506.206,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.268727,0.0733928,-1.00838e-07,-7.19724e-08,4.00025e-11,-60709.5,30.3585], Tmin=(100,'K'), Tmax=(893.398,'K')), NASAPolynomial(coeffs=[24.9503,0.0172215,-1.05808e-06,-1.67833e-10,1.61729e-14,-67480.1,-101.139], Tmin=(893.398,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-506.206,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(440.667,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Oxetane)"""),
)

species(
    label = '[CH2]OO[C](CC)C(=C)O(5868)',
    structure = SMILES('[CH2]OOC(CC)=C([CH2])O'),
    E0 = (-4.73204,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,325,375,415,465,420,450,1700,1750,350,500,795,815,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.8904,0.102277,-0.000102834,5.33044e-08,-1.09785e-11,-388.562,36.2375], Tmin=(100,'K'), Tmax=(1177.95,'K')), NASAPolynomial(coeffs=[19.9379,0.03155,-1.27702e-05,2.33214e-09,-1.60436e-13,-5295.49,-67.6579], Tmin=(1177.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-4.73204,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(CsJOOC) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]OC([CH2])=C(CC)OO(5896)',
    structure = SMILES('[CH2]OC([CH2])=C(CC)OO'),
    E0 = (12.4386,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,3615,1310,387.5,850,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.690028,0.0983639,-9.68064e-05,4.93582e-08,-1.0042e-11,1668.94,37.6869], Tmin=(100,'K'), Tmax=(1188.56,'K')), NASAPolynomial(coeffs=[19.0199,0.0320321,-1.30941e-05,2.40395e-09,-1.65783e-13,-3016.37,-60.8069], Tmin=(1188.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(12.4386,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=COCJ) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]OC([CH2])(O)C(=C)OC(5897)',
    structure = SMILES('[CH2]OC([CH2])(O)C(=C)OC'),
    E0 = (-179.042,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,3615,1277.5,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,180,1011.28,1196.76,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.160384,'amu*angstrom^2'), symmetry=1, barrier=(3.68755,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160384,'amu*angstrom^2'), symmetry=1, barrier=(3.68755,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160384,'amu*angstrom^2'), symmetry=1, barrier=(3.68755,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160384,'amu*angstrom^2'), symmetry=1, barrier=(3.68755,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160384,'amu*angstrom^2'), symmetry=1, barrier=(3.68755,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160384,'amu*angstrom^2'), symmetry=1, barrier=(3.68755,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160384,'amu*angstrom^2'), symmetry=1, barrier=(3.68755,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.65824,0.123943,-0.000138892,7.34154e-08,-1.45272e-11,-21274.7,39.0686], Tmin=(100,'K'), Tmax=(1379.77,'K')), NASAPolynomial(coeffs=[32.3899,0.00954422,-6.15972e-07,-1.15494e-10,1.34116e-14,-29728.6,-136.888], Tmin=(1379.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-179.042,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CsJOCH3) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]OC([CH2])(O)C(O)=CC(5898)',
    structure = SMILES('[CH2]OC([CH2])(O)C(O)=CC'),
    E0 = (-233.751,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,180,1600,1649.82,2843.12,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152517,'amu*angstrom^2'), symmetry=1, barrier=(3.50666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152517,'amu*angstrom^2'), symmetry=1, barrier=(3.50666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152517,'amu*angstrom^2'), symmetry=1, barrier=(3.50666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152517,'amu*angstrom^2'), symmetry=1, barrier=(3.50666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152517,'amu*angstrom^2'), symmetry=1, barrier=(3.50666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152517,'amu*angstrom^2'), symmetry=1, barrier=(3.50666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152517,'amu*angstrom^2'), symmetry=1, barrier=(3.50666,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.70503,0.126885,-0.000146988,8.02368e-08,-1.63792e-11,-27854.7,39.4027], Tmin=(100,'K'), Tmax=(1336.11,'K')), NASAPolynomial(coeffs=[32.4133,0.00891226,-1.33163e-07,-2.29055e-10,2.22688e-14,-36093.3,-135.911], Tmin=(1336.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-233.751,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CsJOCH3) + radical(CJC(C)OC)"""),
)

species(
    label = 'CH2(19)',
    structure = SMILES('[CH2]'),
    E0 = (381.563,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1032.72,2936.3,3459],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.8328,0.000224446,4.68033e-06,-6.04743e-09,2.59009e-12,45920.8,1.40666], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.16229,0.00281798,-7.56235e-07,5.05446e-11,5.65236e-15,46099.1,4.77656], Tmin=(1000,'K'), Tmax=(3000,'K'))], Tmin=(200,'K'), Tmax=(3000,'K'), E0=(381.563,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[CH2]C([O])(O)C(=O)CC(4996)',
    structure = SMILES('[CH2]C([O])(O)C(=O)CC'),
    E0 = (-206.342,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,180,1600,1704.87,2798.05,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153863,'amu*angstrom^2'), symmetry=1, barrier=(3.53761,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153863,'amu*angstrom^2'), symmetry=1, barrier=(3.53761,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153863,'amu*angstrom^2'), symmetry=1, barrier=(3.53761,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153863,'amu*angstrom^2'), symmetry=1, barrier=(3.53761,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153863,'amu*angstrom^2'), symmetry=1, barrier=(3.53761,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.520056,0.0941117,-0.000106227,6.02089e-08,-1.32794e-11,-24649.9,32.6512], Tmin=(100,'K'), Tmax=(1115.68,'K')), NASAPolynomial(coeffs=[19.9145,0.0208496,-7.7298e-06,1.35337e-09,-9.12887e-14,-29209.7,-68.1708], Tmin=(1115.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-206.342,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJC(O)2C) + radical(C=OCOJ)"""),
)

species(
    label = '[CH2]O[C](O)C(=O)CC(5899)',
    structure = SMILES('[CH2]OC(O)=C([O])CC'),
    E0 = (-248.513,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.803558,0.0919494,-9.87095e-05,5.20579e-08,-1.04705e-11,-29704.6,32.6445], Tmin=(100,'K'), Tmax=(1311.28,'K')), NASAPolynomial(coeffs=[22.9136,0.0141377,-3.44923e-06,4.49243e-10,-2.53658e-14,-35454.8,-86.4135], Tmin=(1311.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-248.513,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsCs) + radical(C=C(C)OJ) + radical(C=COCJ)"""),
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
    E0 = (-247.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-143.047,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-143.498,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-247.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-63.3449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-107.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-75.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-94.4001,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-160.722,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-200.889,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-193.534,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-173.823,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (7.96041,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-60.9059,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-209.463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (163.939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (197.948,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-45.6307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-239.047,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (138.382,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (155.552,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (134.758,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-105.904,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (175.221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (133.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]OC([CH2])(O)C(=O)CC(5498)'],
    products = ['CH2O(15)', 'C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]OC([CH2])(O)C(=O)CC(5498)'],
    products = ['[CH2]OC1(O)CC1([O])CC(5883)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.34238e+09,'s^-1'), n=0.889391, Ea=(104.284,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonylbond_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 101.0 to 104.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]OC([CH2])(O)C(=O)CC(5498)'],
    products = ['[CH2]C1(O)OCC1([O])CC(5884)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.46159e+06,'s^-1'), n=1.55572, Ea=(103.833,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5_SS_CO;carbonylbond_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 101.4 to 103.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2O(15)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['[CH2]OC([CH2])(O)C(=O)CC(5498)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2330,'cm^3/(mol*s)'), n=3.17, Ea=(56.6665,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Od_CO-HH;YJ] for rate rule [Od_CO-HH;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 54.5 to 56.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction5',
    reactants = ['C2H5CO(71)', '[CH2]OC(=C)O(3587)'],
    products = ['[CH2]OC([CH2])(O)C(=O)CC(5498)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=0, Ea=(93.9308,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_R;CO_rad] for rate rule [Cds-OsOs_Cds;CO_rad/NonDe]
Euclidian distance = 2.82842712475
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['OH(5)', '[CH2]OC(=C)C(=O)CC(5885)'],
    products = ['[CH2]OC([CH2])(O)C(=O)CC(5498)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.812049,'m^3/(mol*s)'), n=2.01336, Ea=(12.3541,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;YJ] for rate rule [Cds-COOs_Cds;OJ_pri]
Euclidian distance = 2.82842712475
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]OC(C)([O])C(=O)CC(5886)'],
    products = ['[CH2]OC([CH2])(O)C(=O)CC(5498)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C([O])(OC)C(=O)CC(5887)'],
    products = ['[CH2]OC([CH2])(O)C(=O)CC(5498)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.3e+09,'s^-1','*|/',2.51189), n=0.75, Ea=(97.0688,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 0 used for R4H_SSS;O_rad_out;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]OC([CH2])(O)C(=O)CC(5498)'],
    products = ['[CH2]OC(C)(O)C(=O)[CH]C(5888)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(6.44e+09,'s^-1'), n=0.13, Ea=(86.6088,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 131 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]OC([CH2])(O)C(=O)CC(5498)'],
    products = ['[CH2]C(O)(OC)C(=O)[CH]C(5889)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(262000,'s^-1'), n=1.62, Ea=(46.4424,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5H_SSSS;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]CC(=O)C(C)(O)O[CH2](5890)'],
    products = ['[CH2]OC([CH2])(O)C(=O)CC(5498)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(68850,'s^-1'), n=1.68, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_2H;Cs_H_out_2H] for rate rule [R5H_CCC(O2d)C;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]CC(=O)C([CH2])(O)OC(5891)'],
    products = ['[CH2]OC([CH2])(O)C(=O)CC(5498)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3690,'s^-1'), n=1.79, Ea=(49.7896,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 113 used for R6H_SSSSS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R6H_SSSSS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][O](167)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['[CH2]OC([CH2])(O)C(=O)CC(5498)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]OC([CH2])(O)C(=O)CC(5498)'],
    products = ['[CH2]OC1(O)CO[C]1CC(5279)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.01137e+09,'s^-1'), n=0.572544, Ea=(186.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4_S;multiplebond_intra;radadd_intra_cs2H] + [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]OC([CH2])(O)C(=O)CC(5498)'],
    products = ['[CH2]C1(O)OCO[C]1CC(5892)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.64784e+07,'s^-1'), n=0.990488, Ea=(37.8683,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5_SS_CO;carbonyl_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CH3CH2OH(54)', '[CH2]OC([CH2])=C=O(5893)'],
    products = ['[CH2]OC([CH2])(O)C(=O)CC(5498)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.79e-11,'m^3/(mol*s)'), n=3.97, Ea=(329.281,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [doublebond;Cs_OH] for rate rule [Cdd_Cd;C_pri_OH]
Euclidian distance = 1.41421356237
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CH2(S)(23)', '[CH2]OC([CH2])(O)C(C)=O(5894)'],
    products = ['[CH2]OC([CH2])(O)C(=O)CC(5498)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;C_pri/De]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]O[C](O)CC(=O)CC(5895)'],
    products = ['[CH2]OC([CH2])(O)C(=O)CC(5498)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.77139e+09,'s^-1'), n=0.981, Ea=(184.175,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ;CO] + [cCs(-HH)CJ;CsJ;C] for rate rule [cCs(-HH)CJ;CsJ;CO]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]OC([CH2])(O)C(=O)CC(5498)'],
    products = ['CCC(=O)C1(O)CCO1(5870)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]OO[C](CC)C(=C)O(5868)'],
    products = ['[CH2]OC([CH2])(O)C(=O)CC(5498)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O] for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O_R]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]OC([CH2])=C(CC)OO(5896)'],
    products = ['[CH2]OC([CH2])(O)C(=O)CC(5498)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O] for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O_R]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]OC([CH2])(O)C(=C)OC(5897)'],
    products = ['[CH2]OC([CH2])(O)C(=O)CC(5498)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C] for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_Cs;R_O_C]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]OC([CH2])(O)C(O)=CC(5898)'],
    products = ['[CH2]OC([CH2])(O)C(=O)CC(5498)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1713.16,'s^-1'), n=2.9, Ea=(127.847,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_Cs;R_O_H] for rate rule [R_ROR;R1_doublebond_CHCH3;R2_doublebond_Cs;R_O_H]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CH2(19)', '[CH2]C([O])(O)C(=O)CC(4996)'],
    products = ['[CH2]OC([CH2])(O)C(=O)CC(5498)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['CH2(19)', '[CH2]O[C](O)C(=O)CC(5899)'],
    products = ['[CH2]OC([CH2])(O)C(=O)CC(5498)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/ODMustO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '1322',
    isomers = [
        '[CH2]OC([CH2])(O)C(=O)CC(5498)',
    ],
    reactants = [
        ('CH2O(15)', 'C=C(O)C(=O)CC(4626)'),
        ('CH2O(15)', '[CH2]C(O)=C([O])CC(4557)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '1322',
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

