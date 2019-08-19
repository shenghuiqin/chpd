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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.28638,0.141016,-0.000157743,9.51904e-08,-2.31948e-11,-9408.87,43.6005], Tmin=(100,'K'), Tmax=(994.152,'K')), NASAPolynomial(coeffs=[20.0925,0.0509721,-2.18793e-05,4.08028e-09,-2.82779e-13,-13858.4,-64.2329], Tmin=(994.152,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-80.0934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(CJC(C)(C=O)O)"""),
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
    label = '[CH2]C=CCC1(O)CC1([O])CC(19821)',
    structure = SMILES('[CH2]C=CCC1(O)CC1([O])CC'),
    E0 = (14.4351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.07014,0.116384,-7.85573e-05,1.22154e-08,5.41784e-12,1969.67,42.1028], Tmin=(100,'K'), Tmax=(1014.98,'K')), NASAPolynomial(coeffs=[25.7249,0.0400928,-1.49453e-05,2.70691e-09,-1.89394e-13,-5385.2,-100.842], Tmin=(1014.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(14.4351,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopropane) + radical(CC(C)2OJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]C1CC(O)(C1)C(=O)CC(19822)',
    structure = SMILES('[CH2][CH]C1CC(O)(C1)C(=O)CC'),
    E0 = (-9.81685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.34548,0.110493,-8.55651e-05,3.45734e-08,-5.73137e-12,-982.762,43.2166], Tmin=(100,'K'), Tmax=(1407.13,'K')), NASAPolynomial(coeffs=[20.1162,0.0494843,-2.05295e-05,3.76082e-09,-2.56977e-13,-7022.6,-67.6537], Tmin=(1407.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-9.81685,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutane) + radical(RCCJ) + radical(Cs_S)"""),
)

species(
    label = '[CH2]C1(O)CC=CCC1([O])CC(19823)',
    structure = SMILES('[CH2]C1(O)CC=CCC1([O])CC'),
    E0 = (-32.5857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.91028,0.112845,-7.13855e-05,8.97001e-09,5.3574e-12,-3691.54,37.9232], Tmin=(100,'K'), Tmax=(1040.64,'K')), NASAPolynomial(coeffs=[24.3451,0.0434348,-1.67545e-05,3.06809e-09,-2.1474e-13,-10862.2,-97.9871], Tmin=(1040.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-32.5857,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(590.328,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclohexene) + radical(CC(C)2OJ) + radical(CJC(C)2O)"""),
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
    label = '[CH2]C(O)(CC=C=C)C(=O)CC(19824)',
    structure = SMILES('[CH2]C(O)(CC=C=C)C(=O)CC'),
    E0 = (-54.9883,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,180,180,180,543.609,604.543,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.158349,'amu*angstrom^2'), symmetry=1, barrier=(3.64075,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158349,'amu*angstrom^2'), symmetry=1, barrier=(3.64075,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158349,'amu*angstrom^2'), symmetry=1, barrier=(3.64075,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158349,'amu*angstrom^2'), symmetry=1, barrier=(3.64075,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158349,'amu*angstrom^2'), symmetry=1, barrier=(3.64075,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158349,'amu*angstrom^2'), symmetry=1, barrier=(3.64075,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158349,'amu*angstrom^2'), symmetry=1, barrier=(3.64075,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (153.198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.25421,0.143177,-0.00017515,1.16428e-07,-3.11948e-11,-6393.22,42.2627], Tmin=(100,'K'), Tmax=(908.076,'K')), NASAPolynomial(coeffs=[18.7603,0.0506089,-2.22416e-05,4.16844e-09,-2.88621e-13,-10209.7,-57.0934], Tmin=(908.076,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-54.9883,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CJC(C)(C=O)O)"""),
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
    label = '[CH2]C=CCC(=C)O(19825)',
    structure = SMILES('[CH2]C=CCC(=C)O'),
    E0 = (-10.3352,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,343.785,343.786],'cm^-1')),
        HinderedRotor(inertia=(0.208803,'amu*angstrom^2'), symmetry=1, barrier=(17.512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208802,'amu*angstrom^2'), symmetry=1, barrier=(17.512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208803,'amu*angstrom^2'), symmetry=1, barrier=(17.512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208803,'amu*angstrom^2'), symmetry=1, barrier=(17.512,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.574876,0.0600409,-1.13075e-05,-3.61043e-08,2.03683e-11,-1105.84,27.261], Tmin=(100,'K'), Tmax=(954.154,'K')), NASAPolynomial(coeffs=[19.1322,0.0194544,-5.9984e-06,1.05733e-09,-7.71983e-14,-6340.94,-70.2723], Tmin=(954.154,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-10.3352,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P)"""),
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
    label = '[CH2]C=CCC(=C)C(=O)CC(19826)',
    structure = SMILES('[CH2]C=CCC(=C)C(=O)CC'),
    E0 = (2.84571,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0917281,0.0924965,-6.39243e-05,2.31059e-08,-3.61961e-12,476.641,35.8533], Tmin=(100,'K'), Tmax=(1379.96,'K')), NASAPolynomial(coeffs=[11.9308,0.0581796,-2.66224e-05,5.08527e-09,-3.54917e-13,-2790.86,-25.0762], Tmin=(1379.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.84571,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(O)(CC=[C]C)C(=O)CC(19827)',
    structure = SMILES('[CH2]C(O)(CC=[C]C)C(=O)CC'),
    E0 = (6.24915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.4157,0.149235,-0.000195039,1.4321e-07,-4.26962e-11,975.177,44.2594], Tmin=(100,'K'), Tmax=(816.956,'K')), NASAPolynomial(coeffs=[16.461,0.0568149,-2.53555e-05,4.74827e-09,-3.26924e-13,-2109.24,-42.9941], Tmin=(816.956,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.24915,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CJC(C)(C=O)O) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=CCC(C)([O])C(=O)CC(19828)',
    structure = SMILES('[CH2]C=CCC(C)([O])C(=O)CC'),
    E0 = (-49.7221,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3100,440,815,1455,1000,375,552.5,462.5,1710,180,180,180,180,1081.7,1146.81,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.164144,'amu*angstrom^2'), symmetry=1, barrier=(3.774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164144,'amu*angstrom^2'), symmetry=1, barrier=(3.774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164144,'amu*angstrom^2'), symmetry=1, barrier=(3.774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164144,'amu*angstrom^2'), symmetry=1, barrier=(3.774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164144,'amu*angstrom^2'), symmetry=1, barrier=(3.774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164144,'amu*angstrom^2'), symmetry=1, barrier=(3.774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164144,'amu*angstrom^2'), symmetry=1, barrier=(3.774,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.38745,0.12591,-0.000135631,8.67877e-08,-2.35419e-11,-5792.63,43.184], Tmin=(100,'K'), Tmax=(878.437,'K')), NASAPolynomial(coeffs=[12.4057,0.0631013,-2.83799e-05,5.39142e-09,-3.76597e-13,-8215.9,-21.5722], Tmin=(878.437,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-49.7221,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(CC(C)(C=O)OJ)"""),
)

species(
    label = '[CH2]C=C[CH]C(C)(O)C(=O)CC(19829)',
    structure = SMILES('[CH2]C=C[CH]C(C)(O)C(=O)CC'),
    E0 = (-174.788,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.03184,0.13365,-0.000134842,7.32507e-08,-1.62272e-11,-20805.7,40.2646], Tmin=(100,'K'), Tmax=(1082.16,'K')), NASAPolynomial(coeffs=[19.811,0.0529125,-2.29306e-05,4.30816e-09,-3.00199e-13,-25533.2,-66.8393], Tmin=(1082.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-174.788,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(C=CCJCO)"""),
)

species(
    label = '[CH2]C(O)(C[C]=CC)C(=O)CC(19830)',
    structure = SMILES('[CH2]C(O)(C[C]=CC)C(=O)CC'),
    E0 = (6.24915,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,180,180,180,513.982,636.666,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.159811,'amu*angstrom^2'), symmetry=1, barrier=(3.67437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159811,'amu*angstrom^2'), symmetry=1, barrier=(3.67437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159811,'amu*angstrom^2'), symmetry=1, barrier=(3.67437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159811,'amu*angstrom^2'), symmetry=1, barrier=(3.67437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159811,'amu*angstrom^2'), symmetry=1, barrier=(3.67437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159811,'amu*angstrom^2'), symmetry=1, barrier=(3.67437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159811,'amu*angstrom^2'), symmetry=1, barrier=(3.67437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159811,'amu*angstrom^2'), symmetry=1, barrier=(3.67437,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.4157,0.149235,-0.000195039,1.4321e-07,-4.26962e-11,975.177,44.2594], Tmin=(100,'K'), Tmax=(816.956,'K')), NASAPolynomial(coeffs=[16.461,0.0568149,-2.53555e-05,4.74827e-09,-3.26924e-13,-2109.24,-42.9941], Tmin=(816.956,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.24915,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CJC(C)(C=O)O) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=[C]CC(C)(O)C(=O)CC(19831)',
    structure = SMILES('[CH2]C=[C]CC(C)(O)C(=O)CC'),
    E0 = (-53.863,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,180,180,180,488.561,660.355,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157059,'amu*angstrom^2'), symmetry=1, barrier=(3.6111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157059,'amu*angstrom^2'), symmetry=1, barrier=(3.6111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157059,'amu*angstrom^2'), symmetry=1, barrier=(3.6111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157059,'amu*angstrom^2'), symmetry=1, barrier=(3.6111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157059,'amu*angstrom^2'), symmetry=1, barrier=(3.6111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157059,'amu*angstrom^2'), symmetry=1, barrier=(3.6111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157059,'amu*angstrom^2'), symmetry=1, barrier=(3.6111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157059,'amu*angstrom^2'), symmetry=1, barrier=(3.6111,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.7944,0.133513,-0.000148452,9.34466e-08,-2.43116e-11,-6274.75,42.6991], Tmin=(100,'K'), Tmax=(923.946,'K')), NASAPolynomial(coeffs=[15.6761,0.0578777,-2.56587e-05,4.84455e-09,-3.37411e-13,-9503.05,-40.2036], Tmin=(923.946,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-53.863,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=CCC(C)(O)C(=O)[CH]C(19832)',
    structure = SMILES('[CH2]C=CCC(C)(O)C([O])=CC'),
    E0 = (-155.233,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.31627,0.124203,-0.000112243,5.27936e-08,-9.85996e-12,-18430.4,46.964], Tmin=(100,'K'), Tmax=(1298.29,'K')), NASAPolynomial(coeffs=[25.2432,0.0392936,-1.41424e-05,2.42005e-09,-1.60089e-13,-25586.5,-93.1892], Tmin=(1298.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-155.233,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(O)([CH]C=CC)C(=O)CC(19833)',
    structure = SMILES('[CH2]C(O)([CH]C=CC)C(=O)CC'),
    E0 = (-114.676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.44027,0.14664,-0.000170873,1.07951e-07,-2.75899e-11,-13564.6,41.0766], Tmin=(100,'K'), Tmax=(948.498,'K')), NASAPolynomial(coeffs=[19.6413,0.0535197,-2.36105e-05,4.44758e-09,-3.09483e-13,-17753.5,-64.2869], Tmin=(948.498,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-114.676,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C=CCJCO) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = '[CH2][C]=CCC(C)(O)C(=O)CC(19834)',
    structure = SMILES('[CH2][C]=CCC(C)(O)C(=O)CC'),
    E0 = (-53.863,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,180,180,180,488.561,660.355,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157059,'amu*angstrom^2'), symmetry=1, barrier=(3.6111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157059,'amu*angstrom^2'), symmetry=1, barrier=(3.6111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157059,'amu*angstrom^2'), symmetry=1, barrier=(3.6111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157059,'amu*angstrom^2'), symmetry=1, barrier=(3.6111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157059,'amu*angstrom^2'), symmetry=1, barrier=(3.6111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157059,'amu*angstrom^2'), symmetry=1, barrier=(3.6111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157059,'amu*angstrom^2'), symmetry=1, barrier=(3.6111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157059,'amu*angstrom^2'), symmetry=1, barrier=(3.6111,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.7944,0.133513,-0.000148452,9.34466e-08,-2.43116e-11,-6274.75,42.6991], Tmin=(100,'K'), Tmax=(923.946,'K')), NASAPolynomial(coeffs=[15.6761,0.0578777,-2.56587e-05,4.84455e-09,-3.37411e-13,-9503.05,-40.2036], Tmin=(923.946,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-53.863,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=CCC(C)(O)C(=O)C[CH2](19835)',
    structure = SMILES('[CH2]C=CCC(C)(O)C(=O)C[CH2]'),
    E0 = (-80.1157,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.95113,0.134626,-0.0001456,8.6606e-08,-2.10556e-11,-9424.46,43.416], Tmin=(100,'K'), Tmax=(990.394,'K')), NASAPolynomial(coeffs=[18.0104,0.0540036,-2.34911e-05,4.40919e-09,-3.06622e-13,-13378.3,-52.6936], Tmin=(990.394,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-80.1157,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CJCC=O) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([O])(CC=CC)C(=O)CC(19836)',
    structure = SMILES('[CH2]C([O])(CC=CC)C(=O)CC'),
    E0 = (10.39,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3100,440,815,1455,1000,375,552.5,462.5,1710,180,180,180,180,965.576,1262.15,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.166937,'amu*angstrom^2'), symmetry=1, barrier=(3.8382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166937,'amu*angstrom^2'), symmetry=1, barrier=(3.8382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166937,'amu*angstrom^2'), symmetry=1, barrier=(3.8382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166937,'amu*angstrom^2'), symmetry=1, barrier=(3.8382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166937,'amu*angstrom^2'), symmetry=1, barrier=(3.8382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166937,'amu*angstrom^2'), symmetry=1, barrier=(3.8382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166937,'amu*angstrom^2'), symmetry=1, barrier=(3.8382,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.64891,0.136361,-0.000157514,9.15868e-08,-1.42864e-11,1441.84,43.5103], Tmin=(100,'K'), Tmax=(626.327,'K')), NASAPolynomial(coeffs=[13.1203,0.062199,-2.81876e-05,5.32467e-09,-3.68776e-13,-803.683,-23.9894], Tmin=(626.327,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(10.39,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CJC(C)(C=O)O) + radical(CC(C)(C=O)OJ)"""),
)

species(
    label = '[CH2]C(O)(CC=CC)C(=O)[CH]C(19837)',
    structure = SMILES('[CH2]C(O)(CC=CC)C([O])=CC'),
    E0 = (-93.2884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.86455,0.125971,-0.000123382,6.52485e-08,-1.39761e-11,-11006.2,45.934], Tmin=(100,'K'), Tmax=(1123.85,'K')), NASAPolynomial(coeffs=[20.049,0.0479755,-1.92803e-05,3.49494e-09,-2.38838e-13,-15931.7,-62.3444], Tmin=(1123.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-93.2884,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]CC(=O)C([CH2])(O)CC=CC(19838)',
    structure = SMILES('[CH2]CC(=O)C([CH2])(O)CC=CC'),
    E0 = (-20.0035,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,180,180,180,553.142,596.05,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.159041,'amu*angstrom^2'), symmetry=1, barrier=(3.65666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159041,'amu*angstrom^2'), symmetry=1, barrier=(3.65666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159041,'amu*angstrom^2'), symmetry=1, barrier=(3.65666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159041,'amu*angstrom^2'), symmetry=1, barrier=(3.65666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159041,'amu*angstrom^2'), symmetry=1, barrier=(3.65666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159041,'amu*angstrom^2'), symmetry=1, barrier=(3.65666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159041,'amu*angstrom^2'), symmetry=1, barrier=(3.65666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159041,'amu*angstrom^2'), symmetry=1, barrier=(3.65666,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.4762,0.149084,-0.000187158,1.28954e-07,-3.5865e-11,-2178.54,44.6396], Tmin=(100,'K'), Tmax=(875.607,'K')), NASAPolynomial(coeffs=[18.4246,0.0536016,-2.35821e-05,4.40834e-09,-3.04197e-13,-5838.62,-53.4177], Tmin=(875.607,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-20.0035,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CJCC=O) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = '[CH]=C[CH2](2461)',
    structure = SMILES('[CH]C=C'),
    E0 = (376.654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,229.711,230.18,230.787],'cm^-1')),
        HinderedRotor(inertia=(1.33306,'amu*angstrom^2'), symmetry=1, barrier=(50.5153,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31912,0.00817959,3.34736e-05,-4.36194e-08,1.58213e-11,45331.5,10.6389], Tmin=(100,'K'), Tmax=(983.754,'K')), NASAPolynomial(coeffs=[5.36755,0.0170743,-6.35108e-06,1.1662e-09,-8.2762e-14,44095,-3.44606], Tmin=(983.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(376.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C([CH2])(O)C(=O)CC(5988)',
    structure = SMILES('[CH2]C([CH2])(O)C(=O)CC'),
    E0 = (-64.1774,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,3615,1277.5,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.24882,0.119041,-0.00016428,1.17739e-07,-3.31244e-11,-7533.04,31.7431], Tmin=(100,'K'), Tmax=(875.483,'K')), NASAPolynomial(coeffs=[18.1268,0.0305211,-1.26247e-05,2.26293e-09,-1.51404e-13,-10925.9,-59.1577], Tmin=(875.483,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-64.1774,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJC(C)(C=O)O) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = '[CH2]C(O)(CC1[CH]C1)C(=O)CC(19839)',
    structure = SMILES('[CH2]C(O)(CC1[CH]C1)C(=O)CC'),
    E0 = (32.3382,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.87413,0.128491,-0.00012743,6.80747e-08,-1.47907e-11,4101.57,43.4176], Tmin=(100,'K'), Tmax=(1104.65,'K')), NASAPolynomial(coeffs=[19.7444,0.0502084,-2.1129e-05,3.92066e-09,-2.71453e-13,-674.574,-63.0309], Tmin=(1104.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(32.3382,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclopropane) + radical(cyclopropane) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = '[CH2]C=CCC1(O)CO[C]1CC(19840)',
    structure = SMILES('[CH2]C=CCC1(O)CO[C]1CC'),
    E0 = (5.64199,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.69004,0.117406,-9.75285e-05,3.83829e-08,-3.65542e-12,890.493,41.4051], Tmin=(100,'K'), Tmax=(927.406,'K')), NASAPolynomial(coeffs=[19.0822,0.0474328,-1.60863e-05,2.64964e-09,-1.72113e-13,-3806.1,-61.7924], Tmin=(927.406,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.64199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Oxetane) + radical(Allyl_P) + radical(C2CsJOCs)"""),
)

species(
    label = '[CH2]C1[CH]CC(O)(C1)C(=O)CC(19841)',
    structure = SMILES('[CH2]C1[CH]CC(O)(C1)C(=O)CC'),
    E0 = (-91.7475,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.33293,0.103972,-6.84527e-05,1.76451e-08,2.2261e-13,-10831.4,43.359], Tmin=(100,'K'), Tmax=(1089.48,'K')), NASAPolynomial(coeffs=[19.1795,0.0480601,-1.81816e-05,3.22717e-09,-2.19278e-13,-16452.3,-62.6433], Tmin=(1089.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-91.7475,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclopentane) + radical(Isobutyl) + radical(Cs_S)"""),
)

species(
    label = '[CH2]C1(O)CC=CCO[C]1CC(19842)',
    structure = SMILES('[CH2]C1(O)CC=CCO[C]1CC'),
    E0 = (-0.256597,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.60066,0.0799701,7.34021e-05,-1.88792e-07,8.95783e-11,211.22,41.529], Tmin=(100,'K'), Tmax=(909.216,'K')), NASAPolynomial(coeffs=[41.1215,0.00625942,6.53653e-06,-1.57033e-09,1.01506e-13,-12279.5,-186.482], Tmin=(909.216,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-0.256597,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(590.328,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptane) + radical(C2CsJOCs) + radical(CJC(C)2O)"""),
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
    label = '[CH2]C=CCC([CH2])=C=O(19843)',
    structure = SMILES('[CH2]C=CCC(=C)[C]=O'),
    E0 = (237.021,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,263.167,263.167],'cm^-1')),
        HinderedRotor(inertia=(0.00243408,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.400323,'amu*angstrom^2'), symmetry=1, barrier=(19.6744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.400314,'amu*angstrom^2'), symmetry=1, barrier=(19.6744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.40032,'amu*angstrom^2'), symmetry=1, barrier=(19.6744,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.794058,0.0687665,-5.10686e-05,1.86868e-08,-2.78958e-12,28623.3,29.6028], Tmin=(100,'K'), Tmax=(1532.45,'K')), NASAPolynomial(coeffs=[15.1157,0.0313841,-1.44778e-05,2.76854e-09,-1.92714e-13,24233.9,-45.6043], Tmin=(1532.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(237.021,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_P) + radical(C=C(C)CJ=O)"""),
)

species(
    label = 'C=C=CCC(C)(O)C(=O)CC(19844)',
    structure = SMILES('C=C=CCC(C)(O)C(=O)CC'),
    E0 = (-266.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.6996,0.130809,-0.000140507,8.50921e-08,-2.13398e-11,-31863.9,40.6569], Tmin=(100,'K'), Tmax=(955.597,'K')), NASAPolynomial(coeffs=[15.8145,0.0574965,-2.5427e-05,4.80582e-09,-3.3526e-13,-35211.1,-43.0428], Tmin=(955.597,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-266.6,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
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
    label = '[CH2]C=CCC([CH2])(O)C(C)=O(19845)',
    structure = SMILES('[CH2]C=CCC([CH2])(O)C(C)=O'),
    E0 = (-54.6766,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1589.57,1600,2907.16,3200],'cm^-1')),
        HinderedRotor(inertia=(0.150668,'amu*angstrom^2'), symmetry=1, barrier=(3.46416,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150668,'amu*angstrom^2'), symmetry=1, barrier=(3.46416,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150668,'amu*angstrom^2'), symmetry=1, barrier=(3.46416,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150668,'amu*angstrom^2'), symmetry=1, barrier=(3.46416,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150668,'amu*angstrom^2'), symmetry=1, barrier=(3.46416,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150668,'amu*angstrom^2'), symmetry=1, barrier=(3.46416,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150668,'amu*angstrom^2'), symmetry=1, barrier=(3.46416,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.68399,0.121054,-0.000126899,6.93502e-08,-1.5074e-11,-6367.96,40.0966], Tmin=(100,'K'), Tmax=(1118.42,'K')), NASAPolynomial(coeffs=[21.6041,0.0377638,-1.5191e-05,2.7625e-09,-1.8951e-13,-11577.1,-74.8611], Tmin=(1118.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-54.6766,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CJC(C)(C=O)O) + radical(Allyl_P)"""),
)

species(
    label = 'C=C[CH]CC[C](O)C(=O)CC(19807)',
    structure = SMILES('[CH2]C=CCCC(O)=C([O])CC'),
    E0 = (-170.81,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.16487,0.119094,-8.57441e-05,1.68379e-08,4.777e-12,-20307.1,45.5882], Tmin=(100,'K'), Tmax=(993.397,'K')), NASAPolynomial(coeffs=[26.4502,0.037523,-1.33857e-05,2.37746e-09,-1.65324e-13,-27652.7,-100.631], Tmin=(993.397,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-170.81,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=CC[C](O)CC(=O)CC(19846)',
    structure = SMILES('[CH2]C=CC[C](O)CC(=O)CC'),
    E0 = (-107.538,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.82698,0.137803,-0.000178363,1.41421e-07,-4.70304e-11,-12733.1,44.3134], Tmin=(100,'K'), Tmax=(759.502,'K')), NASAPolynomial(coeffs=[11.0119,0.0666648,-3.09121e-05,5.88963e-09,-4.09241e-13,-14581.8,-13.4262], Tmin=(759.502,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-107.538,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C2CsJOH) + radical(Allyl_P)"""),
)

species(
    label = 'CCC(=O)C1(O)CC=CCC1(19812)',
    structure = SMILES('CCC(=O)C1(O)CC=CCC1'),
    E0 = (-388.145,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.3716,0.10216,-5.62589e-05,3.38622e-09,5.01774e-12,-46476.2,36.7997], Tmin=(100,'K'), Tmax=(1090.36,'K')), NASAPolynomial(coeffs=[20.7566,0.0481,-1.91957e-05,3.53454e-09,-2.46027e-13,-52913.7,-79.2627], Tmin=(1090.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-388.145,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(590.328,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclohexene)"""),
)

species(
    label = '[CH2]C=CCC([CH2])=C(CC)OO(19847)',
    structure = SMILES('[CH2]C=CCC([CH2])=C(CC)OO'),
    E0 = (167.933,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1310,387.5,850,1000,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.70688,0.118479,-9.8393e-05,4.24133e-08,-7.44172e-12,20408.7,45.2527], Tmin=(100,'K'), Tmax=(1344.6,'K')), NASAPolynomial(coeffs=[21.7636,0.0486577,-2.05022e-05,3.7943e-09,-2.61354e-13,14097,-74.9285], Tmin=(1344.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(167.933,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = 'C=C[CH]CO[C](CC)C(=C)O(19805)',
    structure = SMILES('[CH2]C(O)=C(CC)OC[CH]C=C'),
    E0 = (-144.473,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,200,800,1000,1200,1400,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.83233,0.131458,-0.000100543,1.99598e-08,6.3655e-12,-17113.2,42.8377], Tmin=(100,'K'), Tmax=(973.609,'K')), NASAPolynomial(coeffs=[31.3436,0.0321506,-1.08674e-05,1.9155e-09,-1.34994e-13,-25716.1,-131.131], Tmin=(973.609,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-144.473,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CCJCO) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]C=CCC([CH2])(O)C(=C)OC(19848)',
    structure = SMILES('[CH2]C=CCC([CH2])(O)C(=C)OC'),
    E0 = (-24.8843,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,180,180,893.259,1313.5,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157309,'amu*angstrom^2'), symmetry=1, barrier=(3.61685,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157309,'amu*angstrom^2'), symmetry=1, barrier=(3.61685,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157309,'amu*angstrom^2'), symmetry=1, barrier=(3.61685,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157309,'amu*angstrom^2'), symmetry=1, barrier=(3.61685,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157309,'amu*angstrom^2'), symmetry=1, barrier=(3.61685,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157309,'amu*angstrom^2'), symmetry=1, barrier=(3.61685,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157309,'amu*angstrom^2'), symmetry=1, barrier=(3.61685,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157309,'amu*angstrom^2'), symmetry=1, barrier=(3.61685,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.7906,0.131601,-0.000123679,5.91127e-08,-1.10843e-11,-2733.44,47.8727], Tmin=(100,'K'), Tmax=(1302.52,'K')), NASAPolynomial(coeffs=[28.7179,0.034839,-1.22463e-05,2.07775e-09,-1.37205e-13,-10941.5,-112.465], Tmin=(1302.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-24.8843,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]C=CCC([CH2])(O)C(O)=CC(19849)',
    structure = SMILES('[CH2]C=CCC([CH2])(O)C(O)=CC'),
    E0 = (-79.594,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,180,180,1557.28,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15032,'amu*angstrom^2'), symmetry=1, barrier=(3.45615,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15032,'amu*angstrom^2'), symmetry=1, barrier=(3.45615,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15032,'amu*angstrom^2'), symmetry=1, barrier=(3.45615,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15032,'amu*angstrom^2'), symmetry=1, barrier=(3.45615,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15032,'amu*angstrom^2'), symmetry=1, barrier=(3.45615,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15032,'amu*angstrom^2'), symmetry=1, barrier=(3.45615,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15032,'amu*angstrom^2'), symmetry=1, barrier=(3.45615,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15032,'amu*angstrom^2'), symmetry=1, barrier=(3.45615,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.82935,0.134454,-0.000131501,6.56377e-08,-1.28369e-11,-9313.74,48.1778], Tmin=(100,'K'), Tmax=(1251.75,'K')), NASAPolynomial(coeffs=[28.3564,0.0347983,-1.20804e-05,2.03508e-09,-1.33991e-13,-17121,-109.277], Tmin=(1251.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-79.594,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C=CC(C)(O)CJ) + radical(Allyl_P)"""),
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
    label = '[CH2]C=CC[C](O)C(=O)CC(19850)',
    structure = SMILES('[CH2]C=CCC(O)=C([O])CC'),
    E0 = (-145.76,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.4118,0.102927,-6.68693e-05,2.0323e-09,9.85956e-12,-17321.6,40.6678], Tmin=(100,'K'), Tmax=(967.6,'K')), NASAPolynomial(coeffs=[25.1013,0.0297671,-9.94981e-06,1.74008e-09,-1.21965e-13,-24158.4,-95.1848], Tmin=(967.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-145.76,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=CCC([CH2])(O)C(=O)CC(19851)',
    structure = SMILES('[CH]=CCC([CH2])(O)C(=O)CC'),
    E0 = (51.5291,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,180,180,180,180,1600,1708.48,2801.23,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154745,'amu*angstrom^2'), symmetry=1, barrier=(3.55789,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154745,'amu*angstrom^2'), symmetry=1, barrier=(3.55789,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154745,'amu*angstrom^2'), symmetry=1, barrier=(3.55789,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154745,'amu*angstrom^2'), symmetry=1, barrier=(3.55789,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154745,'amu*angstrom^2'), symmetry=1, barrier=(3.55789,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154745,'amu*angstrom^2'), symmetry=1, barrier=(3.55789,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154745,'amu*angstrom^2'), symmetry=1, barrier=(3.55789,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.80066,0.132469,-0.000165868,1.11999e-07,-3.03143e-11,6402.14,40.1749], Tmin=(100,'K'), Tmax=(901.365,'K')), NASAPolynomial(coeffs=[18.0509,0.0443788,-1.92815e-05,3.58712e-09,-2.47079e-13,2823.24,-53.5368], Tmin=(901.365,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(51.5291,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CJC(C)(C=O)O)"""),
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
    label = '[CH2]C(O)(C=CC=C)C(=O)CC(19853)',
    structure = SMILES('[CH2]C(O)(C=CC=C)C(=O)CC'),
    E0 = (-119.538,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,545.449,603.295,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.158399,'amu*angstrom^2'), symmetry=1, barrier=(3.6419,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158399,'amu*angstrom^2'), symmetry=1, barrier=(3.6419,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158399,'amu*angstrom^2'), symmetry=1, barrier=(3.6419,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158399,'amu*angstrom^2'), symmetry=1, barrier=(3.6419,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158399,'amu*angstrom^2'), symmetry=1, barrier=(3.6419,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158399,'amu*angstrom^2'), symmetry=1, barrier=(3.6419,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158399,'amu*angstrom^2'), symmetry=1, barrier=(3.6419,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (153.198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.12687,0.139811,-0.000167269,1.09146e-07,-2.87668e-11,-14160.8,41.1048], Tmin=(100,'K'), Tmax=(921.878,'K')), NASAPolynomial(coeffs=[18.4446,0.0505522,-2.20356e-05,4.11918e-09,-2.85009e-13,-17953.6,-56.4674], Tmin=(921.878,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-119.538,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC(O)(C=O)CJ)"""),
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
    label = '[CH2]C(O)(CC[C]=C)C(=O)CC(19855)',
    structure = SMILES('[CH2]C(O)(CC[C]=C)C(=O)CC'),
    E0 = (18.4945,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1277.5,1000,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,375,552.5,462.5,1710,180,180,180,180,1105.04,1119.6,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.162485,'amu*angstrom^2'), symmetry=1, barrier=(3.73586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162485,'amu*angstrom^2'), symmetry=1, barrier=(3.73586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162485,'amu*angstrom^2'), symmetry=1, barrier=(3.73586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162485,'amu*angstrom^2'), symmetry=1, barrier=(3.73586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162485,'amu*angstrom^2'), symmetry=1, barrier=(3.73586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162485,'amu*angstrom^2'), symmetry=1, barrier=(3.73586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162485,'amu*angstrom^2'), symmetry=1, barrier=(3.73586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162485,'amu*angstrom^2'), symmetry=1, barrier=(3.73586,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.41654,0.148419,-0.000188063,1.32209e-07,-3.76539e-11,2448.9,44.7039], Tmin=(100,'K'), Tmax=(854.496,'K')), NASAPolynomial(coeffs=[17.494,0.0552205,-2.44695e-05,4.58303e-09,-3.16297e-13,-953.988,-48.2233], Tmin=(854.496,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(18.4945,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJC(C)(C=O)O) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CCCC([CH2])(O)C(=O)CC(19856)',
    structure = SMILES('[CH]=CCCC([CH2])(O)C(=O)CC'),
    E0 = (27.7489,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.42032,0.146902,-0.000178359,1.18169e-07,-3.16058e-11,3563.69,44.6471], Tmin=(100,'K'), Tmax=(909.069,'K')), NASAPolynomial(coeffs=[18.9658,0.0528071,-2.31072e-05,4.32253e-09,-2.98995e-13,-324.83,-56.4908], Tmin=(909.069,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(27.7489,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJC(C)(C=O)O) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C([O])(CCC=C)C(=O)CC(19857)',
    structure = SMILES('[CH2]C([O])(CCC=C)C(=O)CC'),
    E0 = (22.6354,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,375,552.5,462.5,1710,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3100,440,815,1455,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.05319,0.141402,-0.000177651,1.29221e-07,-3.87135e-11,2932.83,45.3406], Tmin=(100,'K'), Tmax=(808.952,'K')), NASAPolynomial(coeffs=[14.3754,0.060169,-2.70247e-05,5.08941e-09,-3.52042e-13,274.81,-30.435], Tmin=(808.952,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(22.6354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJC(C)(C=O)O) + radical(CC(C)(C=O)OJ)"""),
)

species(
    label = '[CH2]C(O)(CCC=C)C(=O)[CH]C(19858)',
    structure = SMILES('[CH2]C(O)(CCC=C)C([O])=CC'),
    E0 = (-81.043,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.02731,0.127129,-0.000123618,6.40707e-08,-1.33548e-11,-9525.59,46.955], Tmin=(100,'K'), Tmax=(1158.28,'K')), NASAPolynomial(coeffs=[21.6168,0.0454758,-1.78738e-05,3.20709e-09,-2.18061e-13,-15002.8,-70.588], Tmin=(1158.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-81.043,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)(O)CJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]CC(=O)C([CH2])(O)CCC=C(19859)',
    structure = SMILES('[CH2]CC(=O)C([CH2])(O)CCC=C'),
    E0 = (-7.7582,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,375,552.5,462.5,1710,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,180,180,1084.65,1137.86,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.161703,'amu*angstrom^2'), symmetry=1, barrier=(3.71787,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161703,'amu*angstrom^2'), symmetry=1, barrier=(3.71787,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161703,'amu*angstrom^2'), symmetry=1, barrier=(3.71787,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161703,'amu*angstrom^2'), symmetry=1, barrier=(3.71787,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161703,'amu*angstrom^2'), symmetry=1, barrier=(3.71787,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161703,'amu*angstrom^2'), symmetry=1, barrier=(3.71787,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161703,'amu*angstrom^2'), symmetry=1, barrier=(3.71787,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161703,'amu*angstrom^2'), symmetry=1, barrier=(3.71787,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.51513,0.148773,-0.000182224,1.21031e-07,-3.23443e-11,-703.223,45.2172], Tmin=(100,'K'), Tmax=(911.136,'K')), NASAPolynomial(coeffs=[19.5745,0.0517961,-2.25689e-05,4.21212e-09,-2.9094e-13,-4728.51,-59.2964], Tmin=(911.136,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-7.7582,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJCC=O) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = '[CH]=C[CH]CC(C)(O)C(=O)CC(19860)',
    structure = SMILES('[CH]C=CCC(C)(O)C(=O)CC'),
    E0 = (-72.5195,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,375,552.5,462.5,1710,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.69998,0.130639,-0.000129525,7.40815e-08,-1.78927e-11,-8521.24,42.7127], Tmin=(100,'K'), Tmax=(980.861,'K')), NASAPolynomial(coeffs=[14.5192,0.064497,-2.83771e-05,5.33494e-09,-3.70907e-13,-11703,-35.2224], Tmin=(980.861,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-72.5195,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C1(O)CC(C=C)O[C]1CC(19861)',
    structure = SMILES('[CH2]C1(O)CC(C=C)O[C]1CC'),
    E0 = (-13.0373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.32811,0.120967,-0.000106918,5.09723e-08,-9.64104e-12,-1324.77,41.3665], Tmin=(100,'K'), Tmax=(1352.84,'K')), NASAPolynomial(coeffs=[23.5844,0.0399085,-1.21182e-05,1.82893e-09,-1.11073e-13,-7929.42,-89.9751], Tmin=(1352.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-13.0373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Tetrahydrofuran) + radical(C2CsJOCs) + radical(CJC(C)2O)"""),
)

species(
    label = 'CCC(=O)C1(O)C[CH][CH]CC1(19862)',
    structure = SMILES('CCC(=O)C1(O)C[CH][CH]CC1'),
    E0 = (-144.513,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.388379,0.0919152,-5.00896e-05,1.24127e-08,-1.2246e-12,-17265.3,33.9411], Tmin=(100,'K'), Tmax=(2139.99,'K')), NASAPolynomial(coeffs=[19.8915,0.055461,-2.45377e-05,4.45268e-09,-2.9469e-13,-25612.7,-74.9883], Tmin=(2139.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-144.513,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(590.328,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclohexane) + radical(cyclohexane) + radical(cyclohexane)"""),
)

species(
    label = 'C=CC=CC(C)(O)C(=O)CC(19863)',
    structure = SMILES('C=CC=CC(C)(O)C(=O)CC'),
    E0 = (-332.761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.8643,0.129833,-0.000132696,7.37866e-08,-1.67204e-11,-39811.3,40.1001], Tmin=(100,'K'), Tmax=(1060.64,'K')), NASAPolynomial(coeffs=[18.9412,0.0513681,-2.17268e-05,4.03585e-09,-2.79435e-13,-44224.7,-61.499], Tmin=(1060.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-332.761,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

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
    label = 'C=CC1CC(O)(C1)C(=O)CC(19813)',
    structure = SMILES('C=CC1CC(O)(C1)C(=O)CC'),
    E0 = (-280.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.56216,0.106414,-6.50864e-05,9.05437e-09,3.93009e-12,-33553,40.1613], Tmin=(100,'K'), Tmax=(1071.22,'K')), NASAPolynomial(coeffs=[21.9687,0.0456254,-1.78811e-05,3.27311e-09,-2.27682e-13,-40147.9,-82.2316], Tmin=(1071.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-280.753,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutane)"""),
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
    E0 = (-80.0934,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (14.4351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-9.81685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (61.3258,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (168.603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-30.0784,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (39.4082,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (45.6224,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (127.638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (88.9559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (39.1506,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (168.17,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-9.55446,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (6.51536,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (338.307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-13.9637,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-27.375,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (97.5218,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-24.8646,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (120.997,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (89.7708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (312.476,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (145.843,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (106.332,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-21.0293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (3.64619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (317.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-55.1202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (365.185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (77.225,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (77.225,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-71.8928,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (311.046,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (-1.35945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (234.105,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (63.5196,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (235.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (433.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (74.8061,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (93.8445,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (-80.0934,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (105.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (218.691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (188.101,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (125.219,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (-51.0564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (19.8562,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (59.7573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (76.3317,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (7.38423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (-1.84637,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (137.698,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (-71.8091,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    products = ['butadiene13(2459)', 'C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    products = ['[CH2]C=CCC1(O)CC1([O])CC(19821)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.34238e+09,'s^-1'), n=0.889391, Ea=(94.5285,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonylbond_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 90.6 to 94.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    products = ['[CH2][CH]C1CC(O)(C1)C(=O)CC(19822)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.20551e+07,'s^-1'), n=1.225, Ea=(70.2766,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 65.5 to 70.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    products = ['[CH2]C1(O)CC=CCC1([O])CC(19823)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(8.21347e+10,'s^-1'), n=0.43, Ea=(141.419,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_SMSR;multiplebond_intra;radadd_intra_cs2H] for rate rule [R7_SMSS_CO;carbonylbond_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.82842712475
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', '[CH2]C(O)(CC=C=C)C(=O)CC(19824)'],
    products = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.42e+08,'cm^3/(mol*s)'), n=1.64, Ea=(11.7989,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2713 used for Ca_Cds-HH;HJ
Exact match found for rate rule [Ca_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2][CH]C=C(2458)', 'C=C(O)C(=O)CC(4626)'],
    products = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.0191526,'m^3/(mol*s)'), n=2.41, Ea=(53.6755,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;CsJ-CdHH] for rate rule [Cds-COOs_Cds;CsJ-CdHH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C2H5CO(71)', '[CH2]C=CCC(=C)O(19825)'],
    products = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=0, Ea=(93.9308,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_R;CO_rad] for rate rule [Cds-OsCs_Cds;CO_rad/NonDe]
Euclidian distance = 2.82842712475
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['OH(5)', '[CH2]C=CCC(=C)C(=O)CC(19826)'],
    products = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.52832,'m^3/(mol*s)'), n=2.02802, Ea=(14.4047,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDeCs_Cds;YJ] for rate rule [Cds-COCs_Cds;OJ_pri]
Euclidian distance = 2.2360679775
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    products = ['[CH2]C(O)(CC=[C]C)C(=O)CC(19827)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C=CCC(C)([O])C(=O)CC(19828)'],
    products = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    products = ['[CH2]C=C[CH]C(C)(O)C(=O)CC(19829)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(25000,'s^-1'), n=2.28, Ea=(119.244,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 85 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/Cd
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    products = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C=[C]CC(C)(O)C(=O)CC(19831)'],
    products = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    products = ['[CH2]C=CCC(C)(O)C(=O)[CH]C(19832)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(6.44e+09,'s^-1'), n=0.13, Ea=(86.6088,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 131 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    products = ['[CH2]C(O)([CH]C=CC)C(=O)CC(19833)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2e+10,'s^-1'), n=0, Ea=(418.4,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;C_rad_out_single;Cs_H_out_H/NonDeC] for rate rule [R4H_SDS;C_rad_out_2H;Cs_H_out_H/(NonDeC/O)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][C]=CCC(C)(O)C(=O)CC(19834)'],
    products = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(263079,'s^-1'), n=1.73643, Ea=(39.8993,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;Cs_H_out_2H] for rate rule [R5H_DSSS;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    products = ['[CH2]C=CCC(C)(O)C(=O)C[CH2](19835)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(68850,'s^-1'), n=1.68, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_2H;Cs_H_out_2H] for rate rule [R5H_CCC(O2d)C;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C([O])(CC=CC)C(=O)CC(19836)'],
    products = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.32924e+07,'s^-1'), n=0.9625, Ea=(87.1318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6H;O_rad_out;Cs_H_out_2H] + [R6H_RSSMS;Y_rad_out;Cs_H_out_2H] for rate rule [R6H_RSSMS;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    products = ['[CH2]C(O)(CC=CC)C(=O)[CH]C(19837)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1062,'s^-1'), n=1.81, Ea=(55.2288,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 116 used for R7H;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R7H;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]CC(=O)C([CH2])(O)CC=CC(19838)'],
    products = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(6.81e+06,'s^-1'), n=1.5, Ea=(141.001,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R8H;C_rad_out_2H;Cs_H_out] for rate rule [R8H;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2][CH]C=C(2458)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.13324e+07,'m^3/(mol*s)'), n=0.074875, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_rad/H2/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -0.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C[CH2](2461)', '[CH2]C([CH2])(O)C(=O)CC(5988)'],
    products = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.67111e+08,'m^3/(mol*s)'), n=-0.424942, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [C_pri_rad;Cd_rad] + [C_rad/H2/Cs;Y_rad] for rate rule [C_rad/H2/Cs;Cd_rad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -0.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    products = ['[CH2]C(O)(CC1[CH]C1)C(=O)CC(19839)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra_cs2H] for rate rule [R3_D;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    products = ['[CH2]C=CCC1(O)CO[C]1CC(19840)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(7.01137e+09,'s^-1'), n=0.572544, Ea=(186.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4_S;multiplebond_intra;radadd_intra_cs2H] + [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    products = ['[CH2]C1[CH]CC(O)(C1)C(=O)CC(19841)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.47614e+07,'s^-1'), n=0.947292, Ea=(59.0641,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_CsCs_RR_D;doublebond_intra_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    products = ['[CH2]C1(O)CC=CCO[C]1CC(19842)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.57149e+06,'s^-1'), n=1.0129, Ea=(83.7396,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;multiplebond_intra;radadd_intra_cs2H] for rate rule [R7_linear;carbonyl_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CH3CH2OH(54)', '[CH2]C=CCC([CH2])=C=O(19843)'],
    products = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.79e-11,'m^3/(mol*s)'), n=3.97, Ea=(329.281,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [doublebond;Cs_OH] for rate rule [Cdd_Cd;C_pri_OH]
Euclidian distance = 1.41421356237
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    products = ['C=C=CCC(C)(O)C(=O)CC(19844)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['CH2(S)(23)', '[CH2]C=CCC([CH2])(O)C(C)=O(19845)'],
    products = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;C_pri/De]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    products = ['C=C[CH]CC[C](O)C(=O)CC(19807)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    products = ['[CH2]C=CC[C](O)CC(=O)CC(19846)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ-HH;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;CO]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    products = ['CCC(=O)C1(O)CC=CCC1(19812)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.53377e+11,'s^-1'), n=0.0685, Ea=(8.20064,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6;C_rad_out_2H;Cpri_rad_out_2H] + [R6_SSSDS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R6_SSSDS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C=CCC([CH2])=C(CC)OO(19847)'],
    products = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O] for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O_R]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=C[CH]CO[C](CC)C(=C)O(19805)'],
    products = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O] for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O_C]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C=CCC([CH2])(O)C(=C)OC(19848)'],
    products = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(37989.5,'s^-1'), n=2.515, Ea=(258.99,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C] + [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O] for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_C]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C=CCC([CH2])(O)C(O)=CC(19849)'],
    products = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O_H] for rate rule [R_ROR;R1_doublebond_CHCH3;R2_doublebond_CsC;R_O_H]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction37',
    reactants = ['CH2(19)', '[CH2]C=CC[C](O)C(=O)CC(19850)'],
    products = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/ODMustO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction38',
    reactants = ['CH2(19)', '[CH]=CCC([CH2])(O)C(=O)CC(19851)'],
    products = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    products = ['[CH2]C1(O)CC(C=C)C1([O])CC(19852)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(902977,'s^-1'), n=1.63829, Ea=(154.9,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs] for rate rule [R5_SS_CO;carbonylbond_intra_Nd;radadd_intra_csHCd]
Euclidian distance = 3.0
family: Intra_R_Add_Exocyclic
Ea raised from 151.1 to 154.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction40',
    reactants = ['H(3)', '[CH2]C(O)(C=CC=C)C(=O)CC(19853)'],
    products = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.35e+08,'cm^3/(mol*s)'), n=1.64, Ea=(1.58992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2557 used for Cds-CsH_Cds-CdH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction41',
    reactants = ['butadiene13(2459)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(0.0534234,'m^3/(mol*s)'), n=2.459, Ea=(8.39435,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CdH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 7.7 to 8.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2]C(O)([CH]CC=C)C(=O)CC(19854)'],
    products = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.682e+10,'s^-1'), n=0.35, Ea=(125.102,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 160 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2]C(O)(CC[C]=C)C(=O)CC(19855)'],
    products = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.9054e+11,'s^-1'), n=0.853, Ea=(200.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    products = ['[CH]=CCCC([CH2])(O)C(=O)CC(19856)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(8.32e+10,'s^-1'), n=0.77, Ea=(268.194,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 195 used for R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleH
Exact match found for rate rule [R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]C([O])(CCC=C)C(=O)CC(19857)'],
    products = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(657459,'s^-1'), n=1.87237, Ea=(102.584,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R4H_SSS;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    products = ['[CH2]C(O)(CCC=C)C(=O)[CH]C(19858)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(0.1016,'s^-1'), n=3.24, Ea=(29.037,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_1H;Cs_H_out_1H] for rate rule [R5H_CC(O2d)CC;C_rad_out_H/Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2]CC(=O)C([CH2])(O)CCC=C(19859)'],
    products = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1160,'s^-1'), n=1.94, Ea=(27.6144,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 88 used for R6H_SSSSS;C_rad_out_2H;Cs_H_out_H/Cd
Exact match found for rate rule [R6H_SSSSS;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH]=C[CH]CC(C)(O)C(=O)CC(19860)'],
    products = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(22.7193,'s^-1'), n=3.21897, Ea=(132.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R6HJ_2;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    products = ['[CH2]C1(O)CC(C=C)O[C]1CC(19861)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(6.8435e+15,'s^-1'), n=-1.17677, Ea=(156.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_csHCd] for rate rule [R5_SS_CO;carbonyl_intra_Nd;radadd_intra_csHCd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    products = ['CCC(=O)C1(O)C[CH][CH]CC1(19862)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(9.63396e+08,'s^-1'), n=0.483333, Ea=(87.4777,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    products = ['C=CC=CC(C)(O)C(=O)CC(19863)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH2]C(C=C)C([CH2])(O)C(=O)CC(19810)'],
    products = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)'],
    products = ['C=CC1CC(O)(C1)C(=O)CC(19813)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

network(
    label = '3599',
    isomers = [
        '[CH2]C(O)(C[CH]C=C)C(=O)CC(19809)',
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
    label = '3599',
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

