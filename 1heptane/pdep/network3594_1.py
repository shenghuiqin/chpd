species(
    label = '[CH2]C(C=C)C([O])(CC)C(=C)O(19804)',
    structure = SMILES('[CH2]C(C=C)C([O])(CC)C(=C)O'),
    E0 = (6.16714,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,180,891.493,1315.34,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15743,'amu*angstrom^2'), symmetry=1, barrier=(3.61962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15743,'amu*angstrom^2'), symmetry=1, barrier=(3.61962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15743,'amu*angstrom^2'), symmetry=1, barrier=(3.61962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15743,'amu*angstrom^2'), symmetry=1, barrier=(3.61962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15743,'amu*angstrom^2'), symmetry=1, barrier=(3.61962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15743,'amu*angstrom^2'), symmetry=1, barrier=(3.61962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15743,'amu*angstrom^2'), symmetry=1, barrier=(3.61962,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.10588,0.119884,-9.2049e-05,2.58179e-08,1.27997e-12,974.314,48.0474], Tmin=(100,'K'), Tmax=(991.066,'K')), NASAPolynomial(coeffs=[24.8373,0.0397338,-1.4018e-05,2.44045e-09,-1.66611e-13,-5770.46,-88.7809], Tmin=(991.066,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.16714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'butadiene13(2459)',
    structure = SMILES('C=CC=C'),
    E0 = (96.4553,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(1.30712,'amu*angstrom^2'), symmetry=1, barrier=(30.0532,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.18,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.80599,0.0102584,6.1726e-05,-9.01643e-08,3.59117e-11,11658.5,12.0621], Tmin=(100,'K'), Tmax=(946.047,'K')), NASAPolynomial(coeffs=[12.4694,0.0100554,-2.41207e-06,4.57077e-10,-3.93161e-14,8010.78,-43.6375], Tmin=(946.047,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(96.4553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), label="""butadiene13""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2]C(C=C)C1(CC)OC1([CH2])O(20144)',
    structure = SMILES('[CH2]C(C=C)C1(CC)OC1([CH2])O'),
    E0 = (52.3695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.60858,0.147322,-0.000148209,7.37469e-08,-1.3608e-11,6644.66,49.9713], Tmin=(100,'K'), Tmax=(1575.19,'K')), NASAPolynomial(coeffs=[33.1712,0.0219983,-8.83515e-07,-4.49265e-10,4.74732e-14,-1611.66,-137.888], Tmin=(1575.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(52.3695,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CJC(O)2C) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1CC1C([O])(CC)C(=C)O(20145)',
    structure = SMILES('[CH2]C1CC1C([O])(CC)C(=C)O'),
    E0 = (30.8702,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.6794,0.103105,-3.47141e-05,-4.00955e-08,2.6328e-11,3937.12,44.5889], Tmin=(100,'K'), Tmax=(950.546,'K')), NASAPolynomial(coeffs=[27.074,0.0350834,-1.09698e-05,1.88245e-09,-1.3288e-13,-3922.43,-105.26], Tmin=(950.546,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(30.8702,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(Isobutyl) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH2]C1(O)CC(C=C)C1([O])CC(19852)',
    structure = SMILES('[CH2]C1(O)CC(C=C)C1([O])CC'),
    E0 = (74.8061,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.09621,0.117031,-7.99055e-05,1.41415e-08,4.52408e-12,9231.45,41.2692], Tmin=(100,'K'), Tmax=(1023.66,'K')), NASAPolynomial(coeffs=[25.6225,0.0408599,-1.53864e-05,2.79477e-09,-1.95453e-13,1872.53,-101.331], Tmin=(1023.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.8061,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(CC(C)2OJ) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH2]C1OC(CC)(C(=C)O)C1[CH2](19966)',
    structure = SMILES('[CH2]C1OC(CC)(C(=C)O)C1[CH2]'),
    E0 = (30.1116,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.86498,0.135599,-0.000129861,6.34361e-08,-1.16173e-11,3936.88,48.0807], Tmin=(100,'K'), Tmax=(1578.1,'K')), NASAPolynomial(coeffs=[28.8917,0.0276058,-3.48301e-06,2.36145e-11,1.64686e-14,-3293.17,-115.047], Tmin=(1578.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(30.1116,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Oxetane) + radical(Isobutyl) + radical(CJC(C)OC)"""),
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
    label = 'C=CC(=C)C([O])(CC)C(=C)O(20146)',
    structure = SMILES('C=CC(=C)C([O])(CC)C(=C)O'),
    E0 = (-99.382,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,2980,3010,3040,3070,3100,1330,1380,1430,900,975,1050,1000,1025,1050,1600,1650,1700,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,180,180,180,180,954.919,1251.21,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15946,'amu*angstrom^2'), symmetry=1, barrier=(3.66629,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15946,'amu*angstrom^2'), symmetry=1, barrier=(3.66629,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15946,'amu*angstrom^2'), symmetry=1, barrier=(3.66629,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15946,'amu*angstrom^2'), symmetry=1, barrier=(3.66629,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15946,'amu*angstrom^2'), symmetry=1, barrier=(3.66629,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15946,'amu*angstrom^2'), symmetry=1, barrier=(3.66629,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (153.198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.08441,0.116425,-8.34353e-05,1.71188e-08,3.67275e-12,-11718.7,42.6119], Tmin=(100,'K'), Tmax=(1029.75,'K')), NASAPolynomial(coeffs=[26.7911,0.0365114,-1.40083e-05,2.58603e-09,-1.8297e-13,-19375.6,-105.845], Tmin=(1029.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-99.382,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(557.07,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH2][CH]C=C(2458)',
    structure = SMILES('[CH2]C=C[CH2]'),
    E0 = (274.714,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.210311,'amu*angstrom^2'), symmetry=1, barrier=(25.2351,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.779031,'amu*angstrom^2'), symmetry=1, barrier=(93.4717,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.56318,0.0223429,1.87067e-05,-3.93099e-08,1.63982e-11,33100.5,13.4097], Tmin=(100,'K'), Tmax=(974.264,'K')), NASAPolynomial(coeffs=[9.82995,0.0151966,-5.22272e-06,9.67656e-10,-7.0786e-14,30607.7,-26.985], Tmin=(974.264,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Allyl_P) + radical(Allyl_P)"""),
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
    label = '[CH2]C(C=C)C(=O)C(=C)O(20147)',
    structure = SMILES('[CH2]C(C=C)C(=O)C(=C)O'),
    E0 = (-74.6305,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,350,440,435,1725,375,552.5,462.5,1710,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.491839,0.103006,-0.000120387,7.54606e-08,-1.91377e-11,-8817.74,32.0542], Tmin=(100,'K'), Tmax=(954.552,'K')), NASAPolynomial(coeffs=[15.2047,0.0372296,-1.7023e-05,3.26894e-09,-2.30118e-13,-11814.3,-42.942], Tmin=(954.552,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-74.6305,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CJC(C)C=O)"""),
)

species(
    label = 'CH2COH(99)',
    structure = SMILES('C=[C]O'),
    E0 = (103.269,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3615,1277.5,1000,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.989114,'amu*angstrom^2'), symmetry=1, barrier=(22.7417,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.1624,0.0134245,5.56346e-06,-1.95511e-08,9.36369e-12,12455.2,10.1544], Tmin=(100,'K'), Tmax=(925.618,'K')), NASAPolynomial(coeffs=[8.19875,0.00453462,-8.93448e-07,1.26083e-10,-9.46513e-15,10971.3,-16.733], Tmin=(925.618,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(103.269,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CH2COH""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C(C=C)C(=O)CC(20148)',
    structure = SMILES('[CH2]C(C=C)C(=O)CC'),
    E0 = (0.226361,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.162,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.920075,0.0806318,-4.27935e-05,-6.78125e-08,8.80461e-11,124.932,28.0904], Tmin=(100,'K'), Tmax=(483.302,'K')), NASAPolynomial(coeffs=[7.21443,0.051536,-2.38707e-05,4.54754e-09,-3.15957e-13,-752.087,-0.47831], Tmin=(483.302,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(0.226361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJC(C)C=O)"""),
)

species(
    label = 'C2H3(30)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36378,0.000265766,2.79621e-05,-3.72987e-08,1.5159e-11,34475,7.9151], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.15027,0.00754021,-2.62998e-06,4.15974e-10,-2.45408e-14,33856.6,1.72812], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(286.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""C2H3""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C=CC([O])(CC)C(=C)O(6516)',
    structure = SMILES('C=CC([O])(CC)C(=C)O'),
    E0 = (-149.589,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,422.355,712.271,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155253,'amu*angstrom^2'), symmetry=1, barrier=(3.56957,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155253,'amu*angstrom^2'), symmetry=1, barrier=(3.56957,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155253,'amu*angstrom^2'), symmetry=1, barrier=(3.56957,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155253,'amu*angstrom^2'), symmetry=1, barrier=(3.56957,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155253,'amu*angstrom^2'), symmetry=1, barrier=(3.56957,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.161,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.623671,0.0878525,-5.65832e-05,8.57472e-09,3.28341e-12,-17813.1,35.9415], Tmin=(100,'K'), Tmax=(1069.66,'K')), NASAPolynomial(coeffs=[20.4152,0.0330793,-1.32915e-05,2.48319e-09,-1.75318e-13,-23681.3,-73.3671], Tmin=(1069.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-149.589,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C[C](C)C([O])(CC)C(=C)O(20149)',
    structure = SMILES('[CH2]C=C(C)C([O])(CC)C(=C)O'),
    E0 = (-73.1703,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,262.04,866.396,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15197,'amu*angstrom^2'), symmetry=1, barrier=(3.49409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15197,'amu*angstrom^2'), symmetry=1, barrier=(3.49409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15197,'amu*angstrom^2'), symmetry=1, barrier=(3.49409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15197,'amu*angstrom^2'), symmetry=1, barrier=(3.49409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15197,'amu*angstrom^2'), symmetry=1, barrier=(3.49409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15197,'amu*angstrom^2'), symmetry=1, barrier=(3.49409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15197,'amu*angstrom^2'), symmetry=1, barrier=(3.49409,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.24133,0.122242,-0.00010271,4.37536e-08,-7.43058e-12,-8563.24,44.6026], Tmin=(100,'K'), Tmax=(1410.73,'K')), NASAPolynomial(coeffs=[26.3926,0.0410524,-1.63822e-05,2.95775e-09,-2.00987e-13,-16642.1,-103.392], Tmin=(1410.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-73.1703,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH2][C](C=C)C(O)(CC)C(=C)O(20150)',
    structure = SMILES('[CH2]C=C([CH2])C(O)(CC)C(=C)O'),
    E0 = (-150.774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.46663,0.124453,-9.31103e-05,2.1403e-08,3.42576e-12,-17885.5,44.4623], Tmin=(100,'K'), Tmax=(1011.41,'K')), NASAPolynomial(coeffs=[28.4254,0.0365467,-1.356e-05,2.46658e-09,-1.73752e-13,-25887.1,-113.589], Tmin=(1011.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-150.774,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)',
    structure = SMILES('[CH2]C(C=C)C(O)([CH]C)C(=C)O'),
    E0 = (-23.0335,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,180,180,180,344.984,784.383,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153972,'amu*angstrom^2'), symmetry=1, barrier=(3.54011,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153972,'amu*angstrom^2'), symmetry=1, barrier=(3.54011,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153972,'amu*angstrom^2'), symmetry=1, barrier=(3.54011,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153972,'amu*angstrom^2'), symmetry=1, barrier=(3.54011,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153972,'amu*angstrom^2'), symmetry=1, barrier=(3.54011,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153972,'amu*angstrom^2'), symmetry=1, barrier=(3.54011,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153972,'amu*angstrom^2'), symmetry=1, barrier=(3.54011,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153972,'amu*angstrom^2'), symmetry=1, barrier=(3.54011,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.29486,0.132516,-0.000125856,6.01436e-08,-1.10024e-11,-2483.82,53.6144], Tmin=(100,'K'), Tmax=(1483.82,'K')), NASAPolynomial(coeffs=[31.0249,0.0273455,-6.7481e-06,8.82588e-10,-4.95903e-14,-11275.8,-120.808], Tmin=(1483.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-23.0335,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(CCJCO)"""),
)

species(
    label = 'C=[C]C(C)C([O])(CC)C(=C)O(20152)',
    structure = SMILES('C=[C]C(C)C([O])(CC)C(=C)O'),
    E0 = (38.9266,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,180,1004.78,1208.81,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.160351,'amu*angstrom^2'), symmetry=1, barrier=(3.68678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160351,'amu*angstrom^2'), symmetry=1, barrier=(3.68678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160351,'amu*angstrom^2'), symmetry=1, barrier=(3.68678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160351,'amu*angstrom^2'), symmetry=1, barrier=(3.68678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160351,'amu*angstrom^2'), symmetry=1, barrier=(3.68678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160351,'amu*angstrom^2'), symmetry=1, barrier=(3.68678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160351,'amu*angstrom^2'), symmetry=1, barrier=(3.68678,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.39537,0.125449,-0.000112778,5.22387e-08,-9.59442e-12,4924.78,47.788], Tmin=(100,'K'), Tmax=(1319.05,'K')), NASAPolynomial(coeffs=[26.1448,0.0389017,-1.43578e-05,2.49583e-09,-1.66652e-13,-2604.41,-97.805], Tmin=(1319.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(38.9266,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=CC(C)C([O])([CH]C)C(=C)O(20153)',
    structure = SMILES('C=CC(C)C([O])([CH]C)C(=C)O'),
    E0 = (0.986866,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.9345,0.11304,-7.10988e-05,4.06188e-09,8.6232e-12,347.804,48.1475], Tmin=(100,'K'), Tmax=(996.605,'K')), NASAPolynomial(coeffs=[25.6259,0.0387884,-1.40757e-05,2.53058e-09,-1.77273e-13,-6951.51,-93.7823], Tmin=(996.605,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(0.986866,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCJCO) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH2]C([C]=C)C(O)(CC)C(=C)O(20154)',
    structure = SMILES('[CH2]C([C]=C)C(O)(CC)C(=C)O'),
    E0 = (14.9063,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,411.055,724.254,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155331,'amu*angstrom^2'), symmetry=1, barrier=(3.57137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155331,'amu*angstrom^2'), symmetry=1, barrier=(3.57137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155331,'amu*angstrom^2'), symmetry=1, barrier=(3.57137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155331,'amu*angstrom^2'), symmetry=1, barrier=(3.57137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155331,'amu*angstrom^2'), symmetry=1, barrier=(3.57137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155331,'amu*angstrom^2'), symmetry=1, barrier=(3.57137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155331,'amu*angstrom^2'), symmetry=1, barrier=(3.57137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155331,'amu*angstrom^2'), symmetry=1, barrier=(3.57137,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.46676,0.130155,-0.000117921,4.69382e-08,-4.32561e-12,2036.52,48.6023], Tmin=(100,'K'), Tmax=(957.092,'K')), NASAPolynomial(coeffs=[26.3948,0.0360113,-1.18734e-05,1.9759e-09,-1.31513e-13,-4700.86,-95.7079], Tmin=(957.092,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(14.9063,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH2]CC(O)(C(=C)O)C([CH2])C=C(20155)',
    structure = SMILES('[CH2]CC(O)(C(=C)O)C([CH2])C=C'),
    E0 = (-17.6892,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,180,309.411,820.877,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153355,'amu*angstrom^2'), symmetry=1, barrier=(3.52593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153355,'amu*angstrom^2'), symmetry=1, barrier=(3.52593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153355,'amu*angstrom^2'), symmetry=1, barrier=(3.52593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153355,'amu*angstrom^2'), symmetry=1, barrier=(3.52593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153355,'amu*angstrom^2'), symmetry=1, barrier=(3.52593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153355,'amu*angstrom^2'), symmetry=1, barrier=(3.52593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153355,'amu*angstrom^2'), symmetry=1, barrier=(3.52593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153355,'amu*angstrom^2'), symmetry=1, barrier=(3.52593,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.11048,0.134576,-0.00013193,6.5739e-08,-1.26611e-11,-1852.96,52.2328], Tmin=(100,'K'), Tmax=(1360.54,'K')), NASAPolynomial(coeffs=[29.9995,0.029864,-8.3614e-06,1.20958e-09,-7.23992e-14,-10180.5,-115.191], Tmin=(1360.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-17.6892,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(RCCJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C=C)C(O)(CC)C(=C)[O](20156)',
    structure = SMILES('[CH2]C(C=C)C(O)(CC)C(=C)[O]'),
    E0 = (-85.1307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.5603,0.127299,-0.000120168,5.93212e-08,-1.15225e-11,-9988.29,49.9702], Tmin=(100,'K'), Tmax=(1276.34,'K')), NASAPolynomial(coeffs=[26.1857,0.0360116,-1.14763e-05,1.81319e-09,-1.14218e-13,-17228.6,-95.3445], Tmin=(1276.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-85.1307,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C(O)C(O)(CC)C([CH2])C=C(20157)',
    structure = SMILES('[CH]=C(O)C(O)(CC)C([CH2])C=C'),
    E0 = (24.1606,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,180,1556.5,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.150267,'amu*angstrom^2'), symmetry=1, barrier=(3.45494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150267,'amu*angstrom^2'), symmetry=1, barrier=(3.45494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150267,'amu*angstrom^2'), symmetry=1, barrier=(3.45494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150267,'amu*angstrom^2'), symmetry=1, barrier=(3.45494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150267,'amu*angstrom^2'), symmetry=1, barrier=(3.45494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150267,'amu*angstrom^2'), symmetry=1, barrier=(3.45494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150267,'amu*angstrom^2'), symmetry=1, barrier=(3.45494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150267,'amu*angstrom^2'), symmetry=1, barrier=(3.45494,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.18286,0.136944,-0.000136814,6.93583e-08,-1.35791e-11,3182.36,51.1069], Tmin=(100,'K'), Tmax=(1335.15,'K')), NASAPolynomial(coeffs=[30.2989,0.0295609,-8.2254e-06,1.18296e-09,-7.05231e-14,-5127.73,-117.74], Tmin=(1335.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(24.1606,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC(C)C([O])(CC)C(=C)O(20158)',
    structure = SMILES('[CH]=CC(C)C([O])(CC)C(=C)O'),
    E0 = (48.181,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,454.304,681.32,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156165,'amu*angstrom^2'), symmetry=1, barrier=(3.59055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156165,'amu*angstrom^2'), symmetry=1, barrier=(3.59055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156165,'amu*angstrom^2'), symmetry=1, barrier=(3.59055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156165,'amu*angstrom^2'), symmetry=1, barrier=(3.59055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156165,'amu*angstrom^2'), symmetry=1, barrier=(3.59055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156165,'amu*angstrom^2'), symmetry=1, barrier=(3.59055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156165,'amu*angstrom^2'), symmetry=1, barrier=(3.59055,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.25135,0.122328,-9.8184e-05,3.30406e-08,-1.91601e-12,6033.04,47.1925], Tmin=(100,'K'), Tmax=(1037.77,'K')), NASAPolynomial(coeffs=[25.6991,0.0395897,-1.47213e-05,2.63245e-09,-1.81647e-13,-1114.11,-95.1734], Tmin=(1037.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(48.181,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH2]CC([O])(C(=C)O)C(C)C=C(20159)',
    structure = SMILES('[CH2]CC([O])(C(=C)O)C(C)C=C'),
    E0 = (6.33117,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,180,961.091,1246.59,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.159149,'amu*angstrom^2'), symmetry=1, barrier=(3.65914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159149,'amu*angstrom^2'), symmetry=1, barrier=(3.65914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159149,'amu*angstrom^2'), symmetry=1, barrier=(3.65914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159149,'amu*angstrom^2'), symmetry=1, barrier=(3.65914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159149,'amu*angstrom^2'), symmetry=1, barrier=(3.65914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159149,'amu*angstrom^2'), symmetry=1, barrier=(3.65914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159149,'amu*angstrom^2'), symmetry=1, barrier=(3.65914,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.13233,0.119442,-9.16307e-05,2.74485e-08,-2.34878e-13,995.628,48.1488], Tmin=(100,'K'), Tmax=(1035.85,'K')), NASAPolynomial(coeffs=[25.2102,0.0402001,-1.50282e-05,2.69832e-09,-1.86709e-13,-6082.2,-91.5478], Tmin=(1035.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.33117,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(RCCJ) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=CC(C)C([O])(CC)C(=C)[O](20160)',
    structure = SMILES('C=CC(C)C([O])(CC)C(=C)[O]'),
    E0 = (-61.1104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.0554,0.117492,-9.72816e-05,4.18401e-08,-7.20398e-12,-7118.8,47.6006], Tmin=(100,'K'), Tmax=(1393.41,'K')), NASAPolynomial(coeffs=[24.4317,0.041457,-1.54313e-05,2.67984e-09,-1.78063e-13,-14500.3,-88.9716], Tmin=(1393.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-61.1104,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C(O)C([O])(CC)C(C)C=C(20161)',
    structure = SMILES('[CH]=C(O)C([O])(CC)C(C)C=C'),
    E0 = (48.181,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,454.305,681.32,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156165,'amu*angstrom^2'), symmetry=1, barrier=(3.59055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156165,'amu*angstrom^2'), symmetry=1, barrier=(3.59055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156165,'amu*angstrom^2'), symmetry=1, barrier=(3.59055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156165,'amu*angstrom^2'), symmetry=1, barrier=(3.59055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156165,'amu*angstrom^2'), symmetry=1, barrier=(3.59055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156165,'amu*angstrom^2'), symmetry=1, barrier=(3.59055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156165,'amu*angstrom^2'), symmetry=1, barrier=(3.59055,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.25134,0.122328,-9.81835e-05,3.30399e-08,-1.91574e-12,6033.04,47.1924], Tmin=(100,'K'), Tmax=(1037.76,'K')), NASAPolynomial(coeffs=[25.6991,0.0395898,-1.47213e-05,2.63246e-09,-1.81648e-13,-1114.08,-95.173], Tmin=(1037.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(48.181,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC([CH2])C(O)(CC)C(=C)O(20162)',
    structure = SMILES('[CH]=CC([CH2])C(O)(CC)C(=C)O'),
    E0 = (24.1606,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.18286,0.136944,-0.000136814,6.93583e-08,-1.35791e-11,3182.36,51.1069], Tmin=(100,'K'), Tmax=(1335.15,'K')), NASAPolynomial(coeffs=[30.2988,0.0295609,-8.22542e-06,1.18296e-09,-7.05235e-14,-5127.72,-117.739], Tmin=(1335.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(24.1606,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C=C)C1(CC)OC[C]1O(20163)',
    structure = SMILES('[CH2]C(C=C)C1(CC)OC[C]1O'),
    E0 = (59.4226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.43543,0.124593,-0.000116992,5.99293e-08,-1.20733e-11,7393.02,45.2473], Tmin=(100,'K'), Tmax=(1317.22,'K')), NASAPolynomial(coeffs=[22.8282,0.039655,-1.09065e-05,1.50002e-09,-8.46103e-14,1450.62,-80.8886], Tmin=(1317.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(59.4226,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Oxetane) + radical(Isobutyl) + radical(C2CsJOH)"""),
)

species(
    label = 'C=C(O)C([O])(CC)C1[CH]CC1(20164)',
    structure = SMILES('C=C(O)C([O])(CC)C1[CH]CC1'),
    E0 = (18.2484,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.32529,0.0962944,-2.6978e-05,-3.61947e-08,2.12244e-11,2404.89,44.8718], Tmin=(100,'K'), Tmax=(998.603,'K')), NASAPolynomial(coeffs=[23.6314,0.0426654,-1.60253e-05,2.96047e-09,-2.11117e-13,-4889.86,-87.063], Tmin=(998.603,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(18.2484,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(cyclobutane) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=CC1CC[C](O)C1([O])CC(20064)',
    structure = SMILES('C=CC1CC[C](O)C1([O])CC'),
    E0 = (-24.3652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.32994,0.0991218,-4.05522e-05,-1.90808e-08,1.46697e-11,-2722.67,41.2222], Tmin=(100,'K'), Tmax=(1007.81,'K')), NASAPolynomial(coeffs=[22.0771,0.0446803,-1.67672e-05,3.05281e-09,-2.14288e-13,-9393.83,-81.5754], Tmin=(1007.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-24.3652,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(590.328,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(CC(C)2OJ) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C1[CH]COC1(CC)C(=C)O(20014)',
    structure = SMILES('[CH2]C1[CH]COC1(CC)C(=C)O'),
    E0 = (-46.3159,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.16086,0.0867709,2.07476e-05,-1.09837e-07,5.651e-11,-5359.47,39.9526], Tmin=(100,'K'), Tmax=(893.38,'K')), NASAPolynomial(coeffs=[28.4299,0.0276986,-3.33716e-06,1.22071e-10,-1.62242e-15,-13576.4,-115.867], Tmin=(893.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-46.3159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Tetrahydrofuran) + radical(CCJCO) + radical(Isobutyl)"""),
)

species(
    label = 'C=CC(=C)C(O)(CC)C(=C)O(20165)',
    structure = SMILES('C=CC(=C)C(O)(CC)C(=C)O'),
    E0 = (-328.485,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.41239,0.122132,-8.48259e-05,9.65445e-09,8.60192e-12,-39259.9,42.5728], Tmin=(100,'K'), Tmax=(983.692,'K')), NASAPolynomial(coeffs=[29.3514,0.0335909,-1.17542e-05,2.11114e-09,-1.49685e-13,-47474.4,-120.137], Tmin=(983.692,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-328.485,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
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
    label = '[CH2]C(C=C)C(C)([O])C(=C)O(20166)',
    structure = SMILES('[CH2]C(C=C)C(C)([O])C(=C)O'),
    E0 = (29.9474,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,180,180,180,297.22,832.645,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152159,'amu*angstrom^2'), symmetry=1, barrier=(3.49844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152159,'amu*angstrom^2'), symmetry=1, barrier=(3.49844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152159,'amu*angstrom^2'), symmetry=1, barrier=(3.49844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152159,'amu*angstrom^2'), symmetry=1, barrier=(3.49844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152159,'amu*angstrom^2'), symmetry=1, barrier=(3.49844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152159,'amu*angstrom^2'), symmetry=1, barrier=(3.49844,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.34658,0.103752,-7.3395e-05,1.13767e-08,6.20671e-12,3806.8,43.0779], Tmin=(100,'K'), Tmax=(961.627,'K')), NASAPolynomial(coeffs=[23.3845,0.032229,-1.07285e-05,1.83233e-09,-1.25292e-13,-2399.07,-82.8041], Tmin=(961.627,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(29.9474,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(Isobutyl)"""),
)

species(
    label = 'C=C[CH]CC([O])(CC)C(=C)O(19803)',
    structure = SMILES('[CH2]C=CCC([O])(CC)C(=C)O'),
    E0 = (-51.6895,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,180,180,180,180,902.467,1303.41,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157691,'amu*angstrom^2'), symmetry=1, barrier=(3.62563,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157691,'amu*angstrom^2'), symmetry=1, barrier=(3.62563,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157691,'amu*angstrom^2'), symmetry=1, barrier=(3.62563,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157691,'amu*angstrom^2'), symmetry=1, barrier=(3.62563,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157691,'amu*angstrom^2'), symmetry=1, barrier=(3.62563,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157691,'amu*angstrom^2'), symmetry=1, barrier=(3.62563,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157691,'amu*angstrom^2'), symmetry=1, barrier=(3.62563,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.15688,0.119489,-8.97057e-05,2.5025e-08,6.07937e-13,-5981.34,46.5382], Tmin=(100,'K'), Tmax=(1037.02,'K')), NASAPolynomial(coeffs=[25.4516,0.0405121,-1.52672e-05,2.75569e-09,-1.91319e-13,-13186.9,-94.7938], Tmin=(1037.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-51.6895,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(Allyl_P)"""),
)

species(
    label = 'C=CC[CH]C([O])(CC)C(=C)O(20167)',
    structure = SMILES('C=CC[CH]C([O])(CC)C(=C)O'),
    E0 = (8.95876,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
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
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.02309,0.11672,-8.56518e-05,2.23318e-08,1.29885e-12,1308,48.7121], Tmin=(100,'K'), Tmax=(1035.48,'K')), NASAPolynomial(coeffs=[24.8559,0.0405006,-1.52385e-05,2.75023e-09,-1.90965e-13,-5738.86,-89.0487], Tmin=(1035.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(8.95876,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCJCO)"""),
)

species(
    label = 'C=CC1COC1(CC)C(=C)O(19816)',
    structure = SMILES('C=CC1COC1(CC)C(=C)O'),
    E0 = (-243.986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.42024,0.0941325,-6.507e-08,-8.65448e-08,4.71704e-11,-29126.1,39.3122], Tmin=(100,'K'), Tmax=(902.609,'K')), NASAPolynomial(coeffs=[28.5352,0.0292808,-5.12867e-06,5.36639e-10,-3.24213e-14,-37299.5,-117.457], Tmin=(902.609,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-243.986,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Oxetane)"""),
)

species(
    label = '[CH2]C(C=C)C([C]=C)(CC)OO(20168)',
    structure = SMILES('[CH2]C(C=C)C([C]=C)(CC)OO'),
    E0 = (300.97,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
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
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.9761,0.126348,-0.000119143,6.02787e-08,-1.23509e-11,36417.8,47.5485], Tmin=(100,'K'), Tmax=(1172.32,'K')), NASAPolynomial(coeffs=[20.9065,0.0482727,-1.92447e-05,3.46947e-09,-2.36272e-13,31052.6,-66.4849], Tmin=(1172.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(300.97,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(C=C)C([O])(CC)C(C)=O(20169)',
    structure = SMILES('[CH2]C(C=C)C([O])(CC)C(C)=O'),
    E0 = (9.7711,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.36057,0.120896,-0.000120211,6.88339e-08,-1.64643e-11,1365.74,45.6833], Tmin=(100,'K'), Tmax=(996.696,'K')), NASAPolynomial(coeffs=[14.5899,0.0568838,-2.3876e-05,4.39897e-09,-3.02527e-13,-1813.88,-31.2159], Tmin=(996.696,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(9.7711,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(CC(C)(C=O)OJ)"""),
)

species(
    label = 'C=C[CH]CC[C](O)C(=O)CC(19807)',
    structure = SMILES('[CH2]C=CCCC(O)=C([O])CC'),
    E0 = (-170.81,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    molecularWeight = (154.206,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4928.5,'J/mol'), sigma=(8.15945,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=769.82 K, Pc=20.59 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.16487,0.119094,-8.57441e-05,1.68379e-08,4.777e-12,-20307.1,45.5882], Tmin=(100,'K'), Tmax=(993.397,'K')), NASAPolynomial(coeffs=[26.4502,0.037523,-1.33857e-05,2.37746e-09,-1.65324e-13,-27652.7,-100.631], Tmin=(993.397,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-170.81,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(Allyl_P)"""),
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
    label = '[CH2]C(C=C)[C](CC)C(=C)O(20170)',
    structure = SMILES('[CH2]C(O)=C(CC)C([CH2])C=C'),
    E0 = (106.894,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1600],'cm^-1')),
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
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.00515,0.122615,-0.000120041,6.25618e-08,-1.29703e-11,13080.7,41.9668], Tmin=(100,'K'), Tmax=(1174.46,'K')), NASAPolynomial(coeffs=[22.3781,0.0395698,-1.39759e-05,2.35478e-09,-1.54265e-13,7353.28,-79.5891], Tmin=(1174.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(106.894,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(C=C(O)CJ)"""),
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
    label = 'C=C[CH]C([O])(CC)C(=C)O(20171)',
    structure = SMILES('C=C[CH]C([O])(CC)C(=C)O'),
    E0 = (-97.0196,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.15212,0.0910576,-2.3837e-05,-4.49703e-08,2.65515e-11,-11463.2,41.6421], Tmin=(100,'K'), Tmax=(966.812,'K')), NASAPolynomial(coeffs=[26.5265,0.0292851,-9.82682e-06,1.79446e-09,-1.31617e-13,-19280.2,-103.705], Tmin=(966.812,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-97.0196,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(C=CCJC(O)C=C)"""),
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
    E0 = (6.16714,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (53.8503,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (43.4047,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (74.8061,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (128.21,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (112.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (6.16714,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (61.2764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (151.653,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (6.16714,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (150.495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (130.172,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (81.4478,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (82.5285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (190.806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (92.7759,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (116.084,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (65.9908,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (107.975,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (68.4692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (92.4895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (59.0496,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (48.8439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (81.2212,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (117.612,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (89.7708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (131.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (130.432,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (121.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (62.2037,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (69.5673,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (449.809,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (166.102,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (254.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (14.4515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (408.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (210.346,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (153.217,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (349.899,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (284.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    products = ['butadiene13(2459)', 'C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    products = ['[CH2]C(C=C)C1(CC)OC1([CH2])O(20144)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.52e+08,'s^-1'), n=0.89, Ea=(47.6831,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H_secNd;radadd_intra] for rate rule [R4_S_D;doublebond_intra_2H_secNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    products = ['[CH2]C1CC1C([O])(CC)C(=C)O(20145)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.32e+08,'s^-1'), n=0.97, Ea=(37.2376,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 335 used for R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H
Exact match found for rate rule [R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    products = ['[CH2]C1(O)CC(C=C)C1([O])CC(19852)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.48e+06,'s^-1'), n=1.25, Ea=(68.639,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 349 used for R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 68.3 to 68.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    products = ['[CH2]C1OC(CC)(C(=C)O)C1[CH2](19966)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 33 used for R5_SS_D;doublebond_intra_2H_pri;radadd_intra_O
Exact match found for rate rule [R5_SS_D;doublebond_intra_2H_pri;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', 'C=CC(=C)C([O])(CC)C(=C)O(20146)'],
    products = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.88727,'m^3/(mol*s)'), n=1.9687, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 18 used for Cds-CdCs_Cds-HH;HJ
Exact match found for rate rule [Cds-CdCs_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -8.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2][CH]C=C(2458)', 'C=C(O)C(=O)CC(4626)'],
    products = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.94e+07,'cm^3/(mol*s)'), n=1.88, Ea=(89.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO-CdCs_O;YJ] for rate rule [CO-CdCs_O;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 84.1 to 89.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['C2H5(29)', '[CH2]C(C=C)C(=O)C(=C)O(20147)'],
    products = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.94e+10,'cm^3/(mol*s)'), n=0, Ea=(28.0328,'kJ/mol'), T0=(1,'K'), Tmin=(333,'K'), Tmax=(363,'K'), comment="""Estimated using template [CO_O;CsJ-CsHH] for rate rule [CO-CdCs_O;CsJ-CsHH]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CH2COH(99)', '[CH2]C(C=C)C(=O)CC(20148)'],
    products = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.16e+10,'cm^3/(mol*s)'), n=0, Ea=(48.1578,'kJ/mol'), T0=(1,'K'), Tmin=(413,'K'), Tmax=(563,'K'), comment="""Estimated using template [CO-CsCs_O;CJ] for rate rule [CO-CsCs_O;CdsJ-O2s]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['butadiene13(2459)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(0.0136943,'m^3/(mol*s)'), n=2.49, Ea=(94.6549,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CdH_Cds-HH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 90.5 to 94.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction11',
    reactants = ['C2H3(30)', 'C=CC([O])(CC)C(=C)O(6516)'],
    products = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6870,'cm^3/(mol*s)'), n=2.41, Ea=(13.7235,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 432 used for Cds-CsH_Cds-HH;CdsJ-H
Exact match found for rate rule [Cds-CsH_Cds-HH;CdsJ-H]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=C[C](C)C([O])(CC)C(=C)O(20149)'],
    products = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.307e+09,'s^-1'), n=1.31, Ea=(203.342,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 163 used for R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    products = ['[CH2][C](C=C)C(O)(CC)C(=C)O(20150)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(C=C)C(O)([CH]C)C(=C)O(20151)'],
    products = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(5.71,'s^-1'), n=3.021, Ea=(105.562,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 319 used for R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    products = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 204 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    products = ['C=CC(C)C([O])([CH]C)C(=C)O(20153)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(6.44e+09,'s^-1'), n=0.13, Ea=(86.6088,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 131 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C([C]=C)C(O)(CC)C(=C)O(20154)'],
    products = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]CC(O)(C(=C)O)C([CH2])C=C(20155)'],
    products = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    products = ['[CH2]C(C=C)C(O)(CC)C(=C)[O](20156)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1070.11,'s^-1'), n=2.50856, Ea=(101.808,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;O_H_out] + [R4H_SSS;O_rad_out;XH_out] + [R4H_SS(Cd)S;Y_rad_out;XH_out] for rate rule [R4H_SS(Cd)S;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C(O)C(O)(CC)C([CH2])C=C(20157)'],
    products = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=CC(C)C([O])(CC)C(=C)O(20158)'],
    products = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]CC([O])(C(=C)O)C(C)C=C(20159)'],
    products = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(68850,'s^-1'), n=1.68, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 111 used for R5H_CCC;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R5H_CCC;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    products = ['C=CC(C)C([O])(CC)C(=C)[O](20160)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2300,'s^-1'), n=1.98, Ea=(42.6768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC(Cd);C_rad_out_2H;XH_out] for rate rule [R5H_CCC(Cd);C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=C(O)C([O])(CC)C(C)C=C(20161)'],
    products = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5H;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    products = ['[CH]=CC([CH2])C(O)(CC)C(=C)O(20162)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.468e+06,'s^-1','*|/',3), n=1.554, Ea=(111.445,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 290 used for R5H_SSSD;O_rad_out;Cd_H_out_singleH
Exact match found for rate rule [R5H_SSSD;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][CH]C=C(2458)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.49215e+07,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    products = ['[CH2]C(C=C)C1(CC)OC[C]1O(20163)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.21748e+08,'s^-1'), n=0.95, Ea=(124.892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R4_S_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    products = ['C=C(O)C([O])(CC)C1[CH]CC1(20164)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.01e+08,'s^-1'), n=1.02, Ea=(124.265,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 18 used for R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    products = ['C=CC1CC[C](O)C1([O])CC(20064)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.71e+11,'s^-1'), n=0.2, Ea=(115.06,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 4 used for R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    products = ['[CH2]C1[CH]COC1(CC)C(=C)O(20014)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.07178e+07,'s^-1'), n=1.01592, Ea=(56.0365,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri_2H;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    products = ['C=CC(=C)C(O)(CC)C(=C)O(20165)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['CH2(S)(23)', '[CH2]C(C=C)C(C)([O])C(=C)O(20166)'],
    products = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    products = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=CC[CH]C([O])(CC)C(=C)O(20167)'],
    products = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HH)CJ;CsJ;C] for rate rule [cCs(-HH)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    products = ['C=CC1COC1(CC)C(=C)O(19816)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R4_SSS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C(C=C)C([C]=C)(CC)OO(20168)'],
    products = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(3.01978e+11,'s^-1'), n=0, Ea=(107.111,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2OOH_S;Y_rad_out] for rate rule [R2OOH_S;Cd_rad_out_double]
Euclidian distance = 2.0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    products = ['[CH2]C(C=C)C([O])(CC)C(C)=O(20169)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(204.179,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    products = ['C=C[CH]CC[C](O)C(=O)CC(19807)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(6.21184e+10,'s^-1'), n=0.288169, Ea=(147.049,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_5_unsaturated_hexane] for rate rule [1_5_hexadiene]
Euclidian distance = 1.0
family: 6_membered_central_C-C_shift"""),
)

reaction(
    label = 'reaction39',
    reactants = ['O(4)', '[CH2]C(C=C)[C](CC)C(=C)O(20170)'],
    products = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/Cs2;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction40',
    reactants = ['CH2(19)', 'C=C[CH]C([O])(CC)C(=C)O(20171)'],
    products = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeC;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '3594',
    isomers = [
        '[CH2]C(C=C)C([O])(CC)C(=C)O(19804)',
    ],
    reactants = [
        ('butadiene13(2459)', 'C=C(O)C(=O)CC(4626)'),
        ('butadiene13(2459)', '[CH2]C(O)=C([O])CC(4557)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '3594',
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

