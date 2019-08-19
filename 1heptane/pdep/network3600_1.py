species(
    label = '[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)',
    structure = SMILES('[CH2]C(C=C)C([CH2])(O)C(=O)CC'),
    E0 = (-22.2369,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,428.544,718.295,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155916,'amu*angstrom^2'), symmetry=1, barrier=(3.58481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155916,'amu*angstrom^2'), symmetry=1, barrier=(3.58481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155916,'amu*angstrom^2'), symmetry=1, barrier=(3.58481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155916,'amu*angstrom^2'), symmetry=1, barrier=(3.58481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155916,'amu*angstrom^2'), symmetry=1, barrier=(3.58481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155916,'amu*angstrom^2'), symmetry=1, barrier=(3.58481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155916,'amu*angstrom^2'), symmetry=1, barrier=(3.58481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155916,'amu*angstrom^2'), symmetry=1, barrier=(3.58481,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.45095,0.14396,-0.000169003,1.07491e-07,-2.73862e-11,-2443.92,45.8821], Tmin=(100,'K'), Tmax=(957.061,'K')), NASAPolynomial(coeffs=[20.3257,0.0487695,-1.98159e-05,3.57396e-09,-2.42303e-13,-6803.79,-63.003], Tmin=(957.061,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-22.2369,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(CJC(C)(C=O)O)"""),
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
    label = '[CH2]C(C=C)C1(O)CC1([O])CC(19909)',
    structure = SMILES('[CH2]C(C=C)C1(O)CC1([O])CC'),
    E0 = (72.2917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.03031,0.116898,-8.12505e-05,1.33507e-08,5.99862e-12,8925.81,43.6529], Tmin=(100,'K'), Tmax=(978.722,'K')), NASAPolynomial(coeffs=[25.2274,0.0391234,-1.35889e-05,2.36688e-09,-1.62662e-13,1979.74,-95.4911], Tmin=(978.722,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(72.2917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(Isobutyl) + radical(CC(C)2OJ)"""),
)

species(
    label = '[CH2]C1CC1C([CH2])(O)C(=O)CC(19910)',
    structure = SMILES('[CH2]C1CC1C([CH2])(O)C(=O)CC'),
    E0 = (2.46618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.5047,0.132986,-0.000132528,6.92275e-08,-1.43225e-11,539.38,44.1358], Tmin=(100,'K'), Tmax=(1178.37,'K')), NASAPolynomial(coeffs=[24.7906,0.0403318,-1.45838e-05,2.50007e-09,-1.65777e-13,-5893.42,-92.0284], Tmin=(1178.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.46618,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclopropane) + radical(Isobutyl) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = '[CH2]C1(O)C(C=C)CC1([O])CC(19911)',
    structure = SMILES('[CH2]C1(O)C(C=C)CC1([O])CC'),
    E0 = (74.8061,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.09621,0.117031,-7.99055e-05,1.41415e-08,4.52408e-12,9231.45,41.2692], Tmin=(100,'K'), Tmax=(1023.66,'K')), NASAPolynomial(coeffs=[25.6225,0.0408599,-1.53864e-05,2.79477e-09,-1.95453e-13,1872.53,-101.331], Tmin=(1023.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.8061,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(CJC(C)2O) + radical(CC(C)2OJ)"""),
)

species(
    label = '[CH2]C1CC(O)(C(=O)CC)C1[CH2](19893)',
    structure = SMILES('[CH2]C1CC(O)(C(=O)CC)C1[CH2]'),
    E0 = (-8.48219,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.08959,0.118427,-0.000101587,4.65253e-08,-8.53567e-12,-787.903,44.1392], Tmin=(100,'K'), Tmax=(1316.72,'K')), NASAPolynomial(coeffs=[23.1068,0.0418839,-1.4389e-05,2.37636e-09,-1.53291e-13,-7423.22,-84.3515], Tmin=(1316.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.48219,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutane) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    label = '[CH2]C(O)(C(=C)C=C)C(=O)CC(19912)',
    structure = SMILES('[CH2]C(O)(C(=C)C=C)C(=O)CC'),
    E0 = (-121.103,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3615,1277.5,1000,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,515.753,631.951,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157751,'amu*angstrom^2'), symmetry=1, barrier=(3.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157751,'amu*angstrom^2'), symmetry=1, barrier=(3.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157751,'amu*angstrom^2'), symmetry=1, barrier=(3.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157751,'amu*angstrom^2'), symmetry=1, barrier=(3.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157751,'amu*angstrom^2'), symmetry=1, barrier=(3.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157751,'amu*angstrom^2'), symmetry=1, barrier=(3.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157751,'amu*angstrom^2'), symmetry=1, barrier=(3.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (153.198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.32811,0.144294,-0.000177487,1.18125e-07,-3.15824e-11,-14341.8,40.6211], Tmin=(100,'K'), Tmax=(911.422,'K')), NASAPolynomial(coeffs=[19.2712,0.0494957,-2.14629e-05,3.9941e-09,-2.75365e-13,-18278.9,-61.579], Tmin=(911.422,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-121.103,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-OdCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(O)(C=O)CJ)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3133.67,'J/mol'), sigma=(5.35118,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=489.47 K, Pc=46.4 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.91313,0.0225012,-8.59328e-06,4.80319e-10,2.48743e-13,-5274.24,14.655], Tmin=(100,'K'), Tmax=(1673.18,'K')), NASAPolynomial(coeffs=[7.46306,0.0158512,-6.42141e-06,1.12499e-09,-7.32055e-14,-7388.54,-11.406], Tmin=(1673.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-44.1874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""CH3CH2CO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C(C=C)C(=C)O(19913)',
    structure = SMILES('[CH2]C(C=C)C(=C)O'),
    E0 = (45.7038,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,308.967,562.42],'cm^-1')),
        HinderedRotor(inertia=(0.627096,'amu*angstrom^2'), symmetry=1, barrier=(14.4182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.627087,'amu*angstrom^2'), symmetry=1, barrier=(14.418,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.627088,'amu*angstrom^2'), symmetry=1, barrier=(14.418,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0642351,'amu*angstrom^2'), symmetry=1, barrier=(14.4179,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0827877,0.0775295,-7.31474e-05,3.58698e-08,-6.86001e-12,5654.3,28.9069], Tmin=(100,'K'), Tmax=(1359.35,'K')), NASAPolynomial(coeffs=[18.099,0.0203041,-5.89137e-06,8.69965e-10,-5.24693e-14,1055.32,-63.1258], Tmin=(1359.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(45.7038,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl)"""),
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
    label = '[CH2]C(O)(C=C)C(=O)CC(6377)',
    structure = SMILES('[CH2]C(O)(C=C)C(=O)CC'),
    E0 = (-171.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1706.04,2810.13,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154561,'amu*angstrom^2'), symmetry=1, barrier=(3.55366,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154561,'amu*angstrom^2'), symmetry=1, barrier=(3.55366,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154561,'amu*angstrom^2'), symmetry=1, barrier=(3.55366,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154561,'amu*angstrom^2'), symmetry=1, barrier=(3.55366,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154561,'amu*angstrom^2'), symmetry=1, barrier=(3.55366,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154561,'amu*angstrom^2'), symmetry=1, barrier=(3.55366,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.161,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.749509,0.113763,-0.000140295,8.89869e-08,-1.85042e-11,-20441.3,33.559], Tmin=(100,'K'), Tmax=(647.718,'K')), NASAPolynomial(coeffs=[12.7354,0.0463508,-2.09198e-05,3.93392e-09,-2.71362e-13,-22521,-28.2104], Tmin=(647.718,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-171.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(O)(C=O)CJ)"""),
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
    label = '[CH2]C(C=C)C(=C)C(=O)CC(19914)',
    structure = SMILES('[CH2]C(C=C)C(=C)C(=O)CC'),
    E0 = (58.8847,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.199,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.250164,0.0964778,-5.37062e-05,-5.48548e-08,7.57733e-11,7203.19,34.7878], Tmin=(100,'K'), Tmax=(500.064,'K')), NASAPolynomial(coeffs=[6.95259,0.0657948,-3.04486e-05,5.83342e-09,-4.08224e-13,6246.17,4.23082], Tmin=(500.064,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(58.8847,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(O)([C](C)C=C)C(=O)CC(19915)',
    structure = SMILES('[CH2]C=C(C)C([CH2])(O)C(=O)CC'),
    E0 = (-94.891,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3615,1277.5,1000,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1589.73,1600,2920.46,3200],'cm^-1')),
        HinderedRotor(inertia=(0.151573,'amu*angstrom^2'), symmetry=1, barrier=(3.48495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151573,'amu*angstrom^2'), symmetry=1, barrier=(3.48495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151573,'amu*angstrom^2'), symmetry=1, barrier=(3.48495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151573,'amu*angstrom^2'), symmetry=1, barrier=(3.48495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151573,'amu*angstrom^2'), symmetry=1, barrier=(3.48495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151573,'amu*angstrom^2'), symmetry=1, barrier=(3.48495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151573,'amu*angstrom^2'), symmetry=1, barrier=(3.48495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151573,'amu*angstrom^2'), symmetry=1, barrier=(3.48495,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.39668,0.149333,-0.000195449,1.45429e-07,-4.41529e-11,-11190.4,42.2789], Tmin=(100,'K'), Tmax=(801.412,'K')), NASAPolynomial(coeffs=[15.6661,0.0591689,-2.66711e-05,5.01446e-09,-3.45918e-13,-14085.2,-40.863], Tmin=(801.412,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-94.891,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(C=CC(O)(C=O)CJ)"""),
)

species(
    label = '[CH2]C(C=C)C(C)([O])C(=O)CC(19916)',
    structure = SMILES('[CH2]C(C=C)C(C)([O])C(=O)CC'),
    E0 = (8.13447,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1056.18,1174.62,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.164861,'amu*angstrom^2'), symmetry=1, barrier=(3.79047,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164861,'amu*angstrom^2'), symmetry=1, barrier=(3.79047,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164861,'amu*angstrom^2'), symmetry=1, barrier=(3.79047,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164861,'amu*angstrom^2'), symmetry=1, barrier=(3.79047,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164861,'amu*angstrom^2'), symmetry=1, barrier=(3.79047,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164861,'amu*angstrom^2'), symmetry=1, barrier=(3.79047,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164861,'amu*angstrom^2'), symmetry=1, barrier=(3.79047,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.59373,0.129422,-0.000149225,1.02593e-07,-2.94342e-11,1174.01,45.6101], Tmin=(100,'K'), Tmax=(839.84,'K')), NASAPolynomial(coeffs=[12.8871,0.0604544,-2.60507e-05,4.82062e-09,-3.30664e-13,-1258.4,-21.7246], Tmin=(839.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(8.13447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(CC(C)(C=O)OJ)"""),
)

species(
    label = '[CH2][C](C=C)C(C)(O)C(=O)CC(19917)',
    structure = SMILES('[CH2]C=C([CH2])C(C)(O)C(=O)CC'),
    E0 = (-156.615,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.95759,0.134641,-0.000143882,8.46632e-08,-2.03974e-11,-18624.9,40.9302], Tmin=(100,'K'), Tmax=(998.066,'K')), NASAPolynomial(coeffs=[17.9749,0.0547577,-2.38264e-05,4.47216e-09,-3.11021e-13,-22603.7,-55.1943], Tmin=(998.066,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-156.615,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)',
    structure = SMILES('[CH2]C(O)(C(=O)CC)C(C)[C]=C'),
    E0 = (10.5226,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,564.351,583.131,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15842,'amu*angstrom^2'), symmetry=1, barrier=(3.64238,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15842,'amu*angstrom^2'), symmetry=1, barrier=(3.64238,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15842,'amu*angstrom^2'), symmetry=1, barrier=(3.64238,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15842,'amu*angstrom^2'), symmetry=1, barrier=(3.64238,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15842,'amu*angstrom^2'), symmetry=1, barrier=(3.64238,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15842,'amu*angstrom^2'), symmetry=1, barrier=(3.64238,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15842,'amu*angstrom^2'), symmetry=1, barrier=(3.64238,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15842,'amu*angstrom^2'), symmetry=1, barrier=(3.64238,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.35289,0.145021,-0.000174398,1.14841e-07,-3.05544e-11,1489.77,44.2293], Tmin=(100,'K'), Tmax=(913.62,'K')), NASAPolynomial(coeffs=[18.7293,0.052718,-2.28506e-05,4.25556e-09,-2.93651e-13,-2362.39,-55.5751], Tmin=(913.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(10.5226,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = '[CH2]C([O])(C(=O)CC)C(C)C=C(19919)',
    structure = SMILES('[CH2]C([O])(C(=O)CC)C(C)C=C'),
    E0 = (14.6635,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.93985,0.137359,-0.000161429,1.08063e-07,-2.97611e-11,1971.6,44.6915], Tmin=(100,'K'), Tmax=(877.285,'K')), NASAPolynomial(coeffs=[15.4915,0.0578799,-2.55334e-05,4.79285e-09,-3.3201e-13,-1086.83,-37.1222], Tmin=(877.285,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(14.6635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJC(C)(C=O)O) + radical(CC(C)(C=O)OJ)"""),
)

species(
    label = '[CH2]C([C]=C)C(C)(O)C(=O)CC(19920)',
    structure = SMILES('[CH2]C([C]=C)C(C)(O)C(=O)CC'),
    E0 = (3.99358,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,493.199,658.879,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157197,'amu*angstrom^2'), symmetry=1, barrier=(3.61427,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157197,'amu*angstrom^2'), symmetry=1, barrier=(3.61427,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157197,'amu*angstrom^2'), symmetry=1, barrier=(3.61427,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157197,'amu*angstrom^2'), symmetry=1, barrier=(3.61427,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157197,'amu*angstrom^2'), symmetry=1, barrier=(3.61427,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157197,'amu*angstrom^2'), symmetry=1, barrier=(3.61427,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157197,'amu*angstrom^2'), symmetry=1, barrier=(3.61427,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157197,'amu*angstrom^2'), symmetry=1, barrier=(3.61427,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.98537,0.136807,-0.000161109,1.07782e-07,-2.94613e-11,691.276,45.0727], Tmin=(100,'K'), Tmax=(886.928,'K')), NASAPolynomial(coeffs=[16.067,0.0553958,-2.34294e-05,4.29821e-09,-2.93561e-13,-2511.09,-39.8542], Tmin=(886.928,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(3.99358,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(C=C)C(C)(O)C(=O)[CH]C(19921)',
    structure = SMILES('[CH2]C(C=C)C(C)(O)C([O])=CC'),
    E0 = (-97.3761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.07503,0.122458,-0.000107632,4.53975e-08,-6.07364e-12,-11483.1,47.7843], Tmin=(100,'K'), Tmax=(998.088,'K')), NASAPolynomial(coeffs=[23.0445,0.0410772,-1.43176e-05,2.43235e-09,-1.61999e-13,-17458.2,-78.168], Tmin=(998.088,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-97.3761,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=CC(C)C([CH2])(O)C(=O)CC(19922)',
    structure = SMILES('[CH]=CC(C)C([CH2])(O)C(=O)CC'),
    E0 = (19.777,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1685.23,2821.25,3200],'cm^-1')),
        HinderedRotor(inertia=(0.1542,'amu*angstrom^2'), symmetry=1, barrier=(3.54537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1542,'amu*angstrom^2'), symmetry=1, barrier=(3.54537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1542,'amu*angstrom^2'), symmetry=1, barrier=(3.54537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1542,'amu*angstrom^2'), symmetry=1, barrier=(3.54537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1542,'amu*angstrom^2'), symmetry=1, barrier=(3.54537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1542,'amu*angstrom^2'), symmetry=1, barrier=(3.54537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1542,'amu*angstrom^2'), symmetry=1, barrier=(3.54537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1542,'amu*angstrom^2'), symmetry=1, barrier=(3.54537,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.41069,0.144192,-0.000167332,1.04566e-07,-2.62702e-11,2606.83,44.363], Tmin=(100,'K'), Tmax=(967.548,'K')), NASAPolynomial(coeffs=[20.3717,0.0500043,-2.1311e-05,3.95243e-09,-2.72768e-13,-1801.73,-64.797], Tmin=(967.548,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(19.777,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJC(C)(C=O)O) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C(O)(C(=O)[CH]C)C(C)C=C(19923)',
    structure = SMILES('[CH2]C(O)(C([O])=CC)C(C)C=C'),
    E0 = (-89.0149,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.22063,0.126743,-0.000120399,6.01559e-08,-1.19798e-11,-10473.6,47.4032], Tmin=(100,'K'), Tmax=(1217.39,'K')), NASAPolynomial(coeffs=[23.694,0.0415939,-1.54817e-05,2.70081e-09,-1.80829e-13,-16783.2,-82.7174], Tmin=(1217.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-89.0149,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH]=CC([CH2])C(C)(O)C(=O)CC(19924)',
    structure = SMILES('[CH]=CC([CH2])C(C)(O)C(=O)CC'),
    E0 = (13.2479,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1630.36,2879.01,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152693,'amu*angstrom^2'), symmetry=1, barrier=(3.51071,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152693,'amu*angstrom^2'), symmetry=1, barrier=(3.51071,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152693,'amu*angstrom^2'), symmetry=1, barrier=(3.51071,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152693,'amu*angstrom^2'), symmetry=1, barrier=(3.51071,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152693,'amu*angstrom^2'), symmetry=1, barrier=(3.51071,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152693,'amu*angstrom^2'), symmetry=1, barrier=(3.51071,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152693,'amu*angstrom^2'), symmetry=1, barrier=(3.51071,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152693,'amu*angstrom^2'), symmetry=1, barrier=(3.51071,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.02463,0.135742,-0.000153138,9.62226e-08,-2.45814e-11,1807.56,45.141], Tmin=(100,'K'), Tmax=(948.392,'K')), NASAPolynomial(coeffs=[17.6391,0.0528055,-2.19626e-05,4.01257e-09,-2.74147e-13,-1922.19,-48.6831], Tmin=(948.392,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(13.2479,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]CC(=O)C(C)(O)C([CH2])C=C(19925)',
    structure = SMILES('[CH2]CC(=O)C(C)(O)C([CH2])C=C'),
    E0 = (-22.2591,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.12168,0.137638,-0.000157089,9.91912e-08,-2.53649e-11,-2459.26,45.7191], Tmin=(100,'K'), Tmax=(949.129,'K')), NASAPolynomial(coeffs=[18.2472,0.0517963,-2.14257e-05,3.90254e-09,-2.66127e-13,-6325.83,-51.4857], Tmin=(949.129,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-22.2591,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJCC=O) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]CC(=O)C([CH2])(O)C(C)C=C(19926)',
    structure = SMILES('[CH2]CC(=O)C([CH2])(O)C(C)C=C'),
    E0 = (-15.7301,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,521.013,624.13,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157581,'amu*angstrom^2'), symmetry=1, barrier=(3.62309,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157581,'amu*angstrom^2'), symmetry=1, barrier=(3.62309,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157581,'amu*angstrom^2'), symmetry=1, barrier=(3.62309,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157581,'amu*angstrom^2'), symmetry=1, barrier=(3.62309,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157581,'amu*angstrom^2'), symmetry=1, barrier=(3.62309,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157581,'amu*angstrom^2'), symmetry=1, barrier=(3.62309,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157581,'amu*angstrom^2'), symmetry=1, barrier=(3.62309,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157581,'amu*angstrom^2'), symmetry=1, barrier=(3.62309,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.50845,0.146096,-0.000171306,1.07562e-07,-2.7064e-11,-1659.96,44.9437], Tmin=(100,'K'), Tmax=(967.569,'K')), NASAPolynomial(coeffs=[20.9784,0.0489978,-2.07758e-05,3.84283e-09,-2.64786e-13,-6204.92,-67.5922], Tmin=(967.569,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-15.7301,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJC(C)(C=O)O) + radical(CJCC=O)"""),
)

species(
    label = '[CH2]C(C=C)C1(O)CO[C]1CC(19927)',
    structure = SMILES('[CH2]C(C=C)C1(O)CO[C]1CC'),
    E0 = (63.4986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.56872,0.116783,-9.53157e-05,3.15792e-08,1.13121e-12,7843.21,42.6735], Tmin=(100,'K'), Tmax=(854.918,'K')), NASAPolynomial(coeffs=[18.7837,0.0461603,-1.45699e-05,2.27434e-09,-1.42616e-13,3464.21,-57.5833], Tmin=(854.918,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(63.4986,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Oxetane) + radical(C2CsJOCs) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(O)(C(=O)CC)C1[CH]CC1(19928)',
    structure = SMILES('[CH2]C(O)(C(=O)CC)C1[CH]CC1'),
    E0 = (-10.1556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.73786,0.121304,-0.000107681,5.07276e-08,-9.73663e-12,-1010.86,42.9377], Tmin=(100,'K'), Tmax=(1239.09,'K')), NASAPolynomial(coeffs=[20.6046,0.0491777,-2.03667e-05,3.74953e-09,-2.58203e-13,-6547.68,-69.6414], Tmin=(1239.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-10.1556,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutane) + radical(cyclobutane) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = '[CH2]C1(O)[C](CC)OCC1C=C(19929)',
    structure = SMILES('[CH2]C1(O)[C](CC)OCC1C=C'),
    E0 = (-8.02995,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.12614,0.120862,-0.000109601,5.47959e-08,-1.097e-11,-733.476,39.826], Tmin=(100,'K'), Tmax=(1240.26,'K')), NASAPolynomial(coeffs=[21.469,0.0431048,-1.35524e-05,2.08918e-09,-1.28445e-13,-6458.68,-78.5731], Tmin=(1240.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.02995,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Tetrahydrofuran) + radical(CJC(C)2O) + radical(C2CsJOCs)"""),
)

species(
    label = '[CH2]C1[CH]CCC1(O)C(=O)CC(19930)',
    structure = SMILES('[CH2]C1[CH]CCC1(O)C(=O)CC'),
    E0 = (-91.7475,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.33295,0.103972,-6.84534e-05,1.7646e-08,2.22253e-13,-10831.4,43.359], Tmin=(100,'K'), Tmax=(1089.48,'K')), NASAPolynomial(coeffs=[19.1796,0.0480599,-1.81815e-05,3.22714e-09,-2.19276e-13,-16452.3,-62.644], Tmin=(1089.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-91.7475,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclopentane) + radical(Cs_S) + radical(Isobutyl)"""),
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
    label = '[CH2]C(=C=O)C([CH2])C=C(19931)',
    structure = SMILES('[CH2]C(C=C)C(=C)[C]=O'),
    E0 = (293.06,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,323.028,323.028],'cm^-1')),
        HinderedRotor(inertia=(0.118552,'amu*angstrom^2'), symmetry=1, barrier=(8.77842,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118552,'amu*angstrom^2'), symmetry=1, barrier=(8.77842,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118552,'amu*angstrom^2'), symmetry=1, barrier=(8.77842,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.34795,'amu*angstrom^2'), symmetry=1, barrier=(25.7647,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.620073,0.0803907,-9.20355e-05,6.38003e-08,-1.89123e-11,35363.3,29.5289], Tmin=(100,'K'), Tmax=(804.459,'K')), NASAPolynomial(coeffs=[8.4272,0.0415704,-1.96487e-05,3.81064e-09,-2.6893e-13,34107.2,-6.43703], Tmin=(804.459,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293.06,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + radical(Isobutyl) + radical(C=C(C)CJ=O)"""),
)

species(
    label = 'C=CC(=C)C(C)(O)C(=O)CC(19932)',
    structure = SMILES('C=CC(=C)C(C)(O)C(=O)CC'),
    E0 = (-334.326,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.04752,0.13408,-0.000141993,8.14602e-08,-1.89412e-11,-39993.2,39.5534], Tmin=(100,'K'), Tmax=(1038.06,'K')), NASAPolynomial(coeffs=[19.6283,0.0505544,-2.12965e-05,3.94485e-09,-2.72646e-13,-44493.3,-65.8295], Tmin=(1038.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-334.326,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-OdCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
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
    label = '[CH2]C(C=C)C([CH2])(O)C(C)=O(19933)',
    structure = SMILES('[CH2]C(C=C)C([CH2])(O)C(C)=O'),
    E0 = (3.18003,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,375,552.5,462.5,1710,180,180,180,180,1588.91,1600,2910.45,3200],'cm^-1')),
        HinderedRotor(inertia=(0.150685,'amu*angstrom^2'), symmetry=1, barrier=(3.46454,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150685,'amu*angstrom^2'), symmetry=1, barrier=(3.46454,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150685,'amu*angstrom^2'), symmetry=1, barrier=(3.46454,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150685,'amu*angstrom^2'), symmetry=1, barrier=(3.46454,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150685,'amu*angstrom^2'), symmetry=1, barrier=(3.46454,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150685,'amu*angstrom^2'), symmetry=1, barrier=(3.46454,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150685,'amu*angstrom^2'), symmetry=1, barrier=(3.46454,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.8,0.123405,-0.000136017,7.88222e-08,-1.80555e-11,594.95,42.2057], Tmin=(100,'K'), Tmax=(1068.31,'K')), NASAPolynomial(coeffs=[21.5152,0.0361071,-1.34418e-05,2.3303e-09,-1.55175e-13,-4386.58,-71.817], Tmin=(1068.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(3.18003,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJC(C)(C=O)O) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C=C)C[C](O)C(=O)CC(19808)',
    structure = SMILES('[CH2]C(C=C)CC(O)=C([O])CC'),
    E0 = (-112.953,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(4873.63,'J/mol'), sigma=(8.12413,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=761.25 K, Pc=20.62 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.12093,0.119539,-8.80957e-05,1.73632e-08,5.69966e-12,-13351.1,47.1247], Tmin=(100,'K'), Tmax=(957.103,'K')), NASAPolynomial(coeffs=[26.0227,0.036442,-1.1968e-05,2.02345e-09,-1.37465e-13,-20319.6,-95.6787], Tmin=(957.103,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-112.953,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C(C=C)[C](O)CC(=O)CC(19934)',
    structure = SMILES('[CH2]C(C=C)[C](O)CC(=O)CC'),
    E0 = (-49.6817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.94634,0.140223,-0.000187749,1.51011e-07,-4.98139e-11,-5770.22,46.4316], Tmin=(100,'K'), Tmax=(814.389,'K')), NASAPolynomial(coeffs=[11.491,0.0640132,-2.85761e-05,5.31646e-09,-3.63061e-13,-7620.29,-13.5581], Tmin=(814.389,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-49.6817,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)',
    structure = SMILES('[CH2]C=CCC([CH2])(O)C(=O)CC'),
    E0 = (-80.0934,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,180,180,180,431.394,712.751,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155923,'amu*angstrom^2'), symmetry=1, barrier=(3.58497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155923,'amu*angstrom^2'), symmetry=1, barrier=(3.58497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155923,'amu*angstrom^2'), symmetry=1, barrier=(3.58497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155923,'amu*angstrom^2'), symmetry=1, barrier=(3.58497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155923,'amu*angstrom^2'), symmetry=1, barrier=(3.58497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155923,'amu*angstrom^2'), symmetry=1, barrier=(3.58497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155923,'amu*angstrom^2'), symmetry=1, barrier=(3.58497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155923,'amu*angstrom^2'), symmetry=1, barrier=(3.58497,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4713.9,'J/mol'), sigma=(7.86321,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=736.30 K, Pc=22 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.28638,0.141016,-0.000157743,9.51904e-08,-2.31948e-11,-9408.87,43.6005], Tmin=(100,'K'), Tmax=(994.152,'K')), NASAPolynomial(coeffs=[20.0925,0.0509721,-2.18793e-05,4.08028e-09,-2.82779e-13,-13858.4,-64.2329], Tmin=(994.152,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-80.0934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = '[CH2]C(O)([CH]CC=C)C(=O)CC(19854)',
    structure = SMILES('[CH2]C(O)([CH]CC=C)C(=O)CC'),
    E0 = (-19.4452,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,180,180,180,180,1106.03,1115.1,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.162321,'amu*angstrom^2'), symmetry=1, barrier=(3.73209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162321,'amu*angstrom^2'), symmetry=1, barrier=(3.73209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162321,'amu*angstrom^2'), symmetry=1, barrier=(3.73209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162321,'amu*angstrom^2'), symmetry=1, barrier=(3.73209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162321,'amu*angstrom^2'), symmetry=1, barrier=(3.73209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162321,'amu*angstrom^2'), symmetry=1, barrier=(3.73209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162321,'amu*angstrom^2'), symmetry=1, barrier=(3.73209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162321,'amu*angstrom^2'), symmetry=1, barrier=(3.73209,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.1599,0.138332,-0.000153979,9.28619e-08,-2.26542e-11,-2119.22,45.8008], Tmin=(100,'K'), Tmax=(992.211,'K')), NASAPolynomial(coeffs=[19.5245,0.0509148,-2.18247e-05,4.06881e-09,-2.81933e-13,-6422.36,-58.6445], Tmin=(992.211,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-19.4452,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCJCO) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = 'C=CC1CCC1(O)C(=O)CC(19814)',
    structure = SMILES('C=CC1CCC1(O)C(=O)CC'),
    E0 = (-280.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.56216,0.106414,-6.50864e-05,9.05437e-09,3.93009e-12,-33553,40.1613], Tmin=(100,'K'), Tmax=(1071.22,'K')), NASAPolynomial(coeffs=[21.9687,0.0456254,-1.78811e-05,3.27311e-09,-2.27682e-13,-40147.9,-82.2316], Tmin=(1071.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-280.753,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutane)"""),
)

species(
    label = '[CH2]C(=C(CC)OO)C([CH2])C=C(19935)',
    structure = SMILES('[CH2]C(=C(CC)OO)C([CH2])C=C'),
    E0 = (223.972,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1390,370,380,2900,435,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.58492,0.126708,-0.000127799,7.27131e-08,-1.72136e-11,27135.5,44.1094], Tmin=(100,'K'), Tmax=(1006.87,'K')), NASAPolynomial(coeffs=[15.7705,0.0577598,-2.50821e-05,4.70196e-09,-3.26707e-13,23640.6,-39.7393], Tmin=(1006.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(223.972,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(C=C)O[C](CC)C(=C)O(19806)',
    structure = SMILES('[CH2]C(O)=C(CC)OC([CH2])C=C'),
    E0 = (-63.8589,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.75635,0.144522,-0.000145669,7.20805e-08,-1.36345e-11,-7379,49.2849], Tmin=(100,'K'), Tmax=(1397.02,'K')), NASAPolynomial(coeffs=[35.0738,0.0241172,-6.48372e-06,9.33827e-10,-5.68284e-14,-17328.1,-147.809], Tmin=(1397.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-63.8589,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(O)CJ) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]C(C=C)C([CH2])(O)C(=C)OC(19936)',
    structure = SMILES('[CH2]C(C=C)C([CH2])(O)C(=C)OC'),
    E0 = (32.9723,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,180,180,895.621,1311.24,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157346,'amu*angstrom^2'), symmetry=1, barrier=(3.6177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157346,'amu*angstrom^2'), symmetry=1, barrier=(3.6177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157346,'amu*angstrom^2'), symmetry=1, barrier=(3.6177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157346,'amu*angstrom^2'), symmetry=1, barrier=(3.6177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157346,'amu*angstrom^2'), symmetry=1, barrier=(3.6177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157346,'amu*angstrom^2'), symmetry=1, barrier=(3.6177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157346,'amu*angstrom^2'), symmetry=1, barrier=(3.6177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157346,'amu*angstrom^2'), symmetry=1, barrier=(3.6177,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.4142,0.128305,-0.00011384,4.52019e-08,-4.62975e-12,4207.92,48.2053], Tmin=(100,'K'), Tmax=(984.366,'K')), NASAPolynomial(coeffs=[25.975,0.0375103,-1.29183e-05,2.20474e-09,-1.48466e-13,-2571.34,-94.3544], Tmin=(984.366,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(32.9723,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)(O)CJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)',
    structure = SMILES('[CH2]C(C=C)C([CH2])(O)C(O)=CC'),
    E0 = (-21.7374,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,180,180,1558.78,1600,2932.02,3200],'cm^-1')),
        HinderedRotor(inertia=(0.150348,'amu*angstrom^2'), symmetry=1, barrier=(3.4568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150348,'amu*angstrom^2'), symmetry=1, barrier=(3.4568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150348,'amu*angstrom^2'), symmetry=1, barrier=(3.4568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150348,'amu*angstrom^2'), symmetry=1, barrier=(3.4568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150348,'amu*angstrom^2'), symmetry=1, barrier=(3.4568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150348,'amu*angstrom^2'), symmetry=1, barrier=(3.4568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150348,'amu*angstrom^2'), symmetry=1, barrier=(3.4568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150348,'amu*angstrom^2'), symmetry=1, barrier=(3.4568,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.51634,0.131825,-0.00012362,5.37432e-08,-7.02178e-12,-2369.51,48.7432], Tmin=(100,'K'), Tmax=(967.443,'K')), NASAPolynomial(coeffs=[26.1957,0.0365554,-1.22553e-05,2.0498e-09,-1.36254e-13,-9022.06,-94.4956], Tmin=(967.443,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-21.7374,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)(O)CJ) + radical(Isobutyl)"""),
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
    label = '[CH2]C(C=C)[C](O)C(=O)CC(19938)',
    structure = SMILES('[CH2]C(C=C)C(O)=C([O])CC'),
    E0 = (-89.7211,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.67764,0.115875,-0.000113209,5.44546e-08,-9.26479e-12,-10578.6,40.9027], Tmin=(100,'K'), Tmax=(981.401,'K')), NASAPolynomial(coeffs=[22.5822,0.0330759,-1.12327e-05,1.87611e-09,-1.2375e-13,-16114.6,-79.6265], Tmin=(981.401,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-89.7211,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C(O)([CH]C=C)C(=O)CC(19939)',
    structure = SMILES('[CH2]C(O)([CH]C=C)C(=O)CC'),
    E0 = (-78.6503,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,245.003,892.17,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.151835,'amu*angstrom^2'), symmetry=1, barrier=(3.49099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151835,'amu*angstrom^2'), symmetry=1, barrier=(3.49099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151835,'amu*angstrom^2'), symmetry=1, barrier=(3.49099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151835,'amu*angstrom^2'), symmetry=1, barrier=(3.49099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151835,'amu*angstrom^2'), symmetry=1, barrier=(3.49099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151835,'amu*angstrom^2'), symmetry=1, barrier=(3.49099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151835,'amu*angstrom^2'), symmetry=1, barrier=(3.49099,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.89882,0.132384,-0.000155266,9.63781e-08,-2.39279e-11,-9249.18,37.3212], Tmin=(100,'K'), Tmax=(980.428,'K')), NASAPolynomial(coeffs=[19.9861,0.0430959,-1.86607e-05,3.48956e-09,-2.42108e-13,-13540.5,-67.8283], Tmin=(980.428,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-78.6503,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CCJCO) + radical(CJC(C)(C=O)O)"""),
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
    E0 = (-22.2369,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (72.2917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (15.0007,'kJ/mol'),
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
    E0 = (36.7575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (90.6895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-22.2369,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (95.4472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-22.2369,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (128.775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (101.661,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (108.451,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (146.812,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (95.7023,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (162.402,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (61.4431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (48.3021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (64.3719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (64.0855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (24.2055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (46.2881,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (30.4815,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (34.0595,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (89.7708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (164.188,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (102.028,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (15.6314,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (92.4047,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (373.582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (41.1633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (423.042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (135.082,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (135.082,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (137.698,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (226.156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (-13.9525,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (367.085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (79.2547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (291.962,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (121.376,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (291.842,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (302.913,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    products = ['butadiene13(2459)', 'C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    products = ['[CH2]C(C=C)C1(O)CC1([O])CC(19909)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.34238e+09,'s^-1'), n=0.889391, Ea=(94.5285,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonylbond_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 90.6 to 94.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    products = ['[CH2]C1CC1C([CH2])(O)C(=O)CC(19910)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.32e+08,'s^-1'), n=0.97, Ea=(37.2376,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 335 used for R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H
Exact match found for rate rule [R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    products = ['[CH2]C1(O)C(C=C)CC1([O])CC(19911)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.46159e+06,'s^-1'), n=1.55572, Ea=(97.043,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5_SS_CO;carbonylbond_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 93.1 to 97.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    products = ['[CH2]C1CC(O)(C(=O)CC)C1[CH2](19893)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.07e+06,'s^-1'), n=1.46, Ea=(58.9944,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 343 used for R5_SS_D;doublebond_intra_2H_pri;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[CH2]C(O)(C(=C)C=C)C(=O)CC(19912)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
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
    products = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00784329,'m^3/(mol*s)'), n=2.41519, Ea=(61.517,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;CJ] for rate rule [Cds-COOs_Cds;CJ]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 59.3 to 61.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['C2H5CO(71)', '[CH2]C(C=C)C(=C)O(19913)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=0, Ea=(93.9308,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_R;CO_rad] for rate rule [Cds-OsCs_Cds;CO_rad/NonDe]
Euclidian distance = 2.82842712475
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['butadiene13(2459)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.0136943,'m^3/(mol*s)'), n=2.49, Ea=(66.2509,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CdH_Cds-HH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 65.7 to 66.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['C2H3(30)', '[CH2]C(O)(C=C)C(=O)CC(6377)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6870,'cm^3/(mol*s)'), n=2.41, Ea=(13.7235,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;CdsJ-H] for rate rule [Cds-Cs\O2s/H_Cds-HH;CdsJ-H]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['OH(5)', '[CH2]C(C=C)C(=C)C(=O)CC(19914)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(0.52832,'m^3/(mol*s)'), n=2.02802, Ea=(14.4047,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDeCs_Cds;YJ] for rate rule [Cds-COCs_Cds;OJ_pri]
Euclidian distance = 2.2360679775
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(O)([C](C)C=C)C(=O)CC(19915)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.307e+09,'s^-1'), n=1.31, Ea=(203.342,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 163 used for R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(C=C)C(C)([O])C(=O)CC(19916)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    products = ['[CH2][C](C=C)C(C)(O)C(=O)CC(19917)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 204 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    products = ['[CH2]C([O])(C(=O)CC)C(C)C=C(19919)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C([C]=C)C(C)(O)C(=O)CC(19920)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    products = ['[CH2]C(C=C)C(C)(O)C(=O)[CH]C(19921)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(6.44e+09,'s^-1'), n=0.13, Ea=(86.6088,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 131 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=CC(C)C([CH2])(O)C(=O)CC(19922)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    products = ['[CH2]C(O)(C(=O)[CH]C)C(C)C=C(19923)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(262000,'s^-1'), n=1.62, Ea=(46.4424,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_2H;Cs_H_out_H/NonDeC] for rate rule [R5H_CC(O2d)CC;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=CC([CH2])C(C)(O)C(=O)CC(19924)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5H;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    products = ['[CH2]CC(=O)C(C)(O)C([CH2])C=C(19925)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(68850,'s^-1'), n=1.68, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_2H;Cs_H_out_2H] for rate rule [R5H_CCC(O2d)C;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]CC(=O)C([CH2])(O)C(C)C=C(19926)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3690,'s^-1'), n=1.79, Ea=(49.7896,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 113 used for R6H_SSSSS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R6H_SSSSS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2][CH]C=C(2458)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.49215e+07,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    products = ['[CH2]C(C=C)C1(O)CO[C]1CC(19927)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(7.01137e+09,'s^-1'), n=0.572544, Ea=(186.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4_S;multiplebond_intra;radadd_intra_cs2H] + [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    products = ['[CH2]C(O)(C(=O)CC)C1[CH]CC1(19928)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.01e+08,'s^-1'), n=1.02, Ea=(124.265,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 18 used for R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    products = ['[CH2]C1(O)[C](CC)OCC1C=C(19929)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.64784e+07,'s^-1'), n=0.990488, Ea=(37.8683,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5_SS_CO;carbonyl_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    products = ['[CH2]C1[CH]CCC1(O)C(=O)CC(19930)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.12e+09,'s^-1'), n=0.63, Ea=(114.642,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 3 used for R5_SS_D;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['CH3CH2OH(54)', '[CH2]C(=C=O)C([CH2])C=C(19931)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.79e-11,'m^3/(mol*s)'), n=3.97, Ea=(329.281,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [doublebond;Cs_OH] for rate rule [Cdd_Cd;C_pri_OH]
Euclidian distance = 1.41421356237
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    products = ['C=CC(=C)C(C)(O)C(=O)CC(19932)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['CH2(S)(23)', '[CH2]C(C=C)C([CH2])(O)C(C)=O(19933)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;C_pri/De]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction55',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    products = ['[CH2]C(C=C)C[C](O)C(=O)CC(19808)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    products = ['[CH2]C(C=C)[C](O)CC(=O)CC(19934)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ-HH;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;CO]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    products = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C(O)([CH]CC=C)C(=O)CC(19854)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HH)CJ;CsJ;C] for rate rule [cCs(-HH)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    products = ['C=CC1CCC1(O)C(=O)CC(19814)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C(=C(CC)OO)C([CH2])C=C(19935)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O] for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O_R]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C(C=C)O[C](CC)C(=C)O(19806)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O] for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O_C]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(=C)OC(19936)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(37989.5,'s^-1'), n=2.515, Ea=(258.99,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C] + [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O] for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_C]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(O)=CC(19937)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O_H] for rate rule [R_ROR;R1_doublebond_CHCH3;R2_doublebond_CsC;R_O_H]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction41',
    reactants = ['CH2(19)', '[CH2]C(C=C)[C](O)C(=O)CC(19938)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/ODMustO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction42',
    reactants = ['CH2(19)', '[CH2]C(O)([CH]C=C)C(=O)CC(19939)'],
    products = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeC;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '3600',
    isomers = [
        '[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)',
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
    label = '3600',
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

