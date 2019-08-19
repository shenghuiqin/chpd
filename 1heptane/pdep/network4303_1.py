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
    label = 'CH3CHCCH2(18175)',
    structure = SMILES('C=C=CC'),
    E0 = (145.615,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.759584,'amu*angstrom^2'), symmetry=1, barrier=(17.4643,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2996.71,'J/mol'), sigma=(5.18551,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=468.08 K, Pc=48.77 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.74635,0.0218189,8.22353e-06,-2.14768e-08,8.55624e-12,17563.6,12.7381], Tmin=(100,'K'), Tmax=(1025.6,'K')), NASAPolynomial(coeffs=[6.82078,0.0192338,-7.45622e-06,1.36536e-09,-9.53195e-14,16028,-10.4333], Tmin=(1025.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.615,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), label="""CH3CHCCH2""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = 'C=[C]C(C)C1(O)CC1([O])CC(25674)',
    structure = SMILES('C=[C]C(C)C1(O)CC1([O])CC'),
    E0 = (105.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.17221,0.120842,-9.69101e-05,3.41204e-08,-2.88127e-12,12869.8,42.8566], Tmin=(100,'K'), Tmax=(1059.13,'K')), NASAPolynomial(coeffs=[24.9399,0.0408474,-1.53408e-05,2.74537e-09,-1.88862e-13,5870.37,-95.4325], Tmin=(1059.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(105.051,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(Cds_S) + radical(CC(C)2OJ)"""),
)

species(
    label = '[CH2]C1(O)C(C)C(=C)C1([O])CC(25675)',
    structure = SMILES('[CH2]C1(O)C(C)C(=C)C1([O])CC'),
    E0 = (50.079,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[8.58962,0.0480643,6.21627e-05,-8.80534e-08,2.11534e-11,5726.93,-4.72753], Tmin=(100,'K'), Tmax=(1797.08,'K')), NASAPolynomial(coeffs=[101.578,0.022656,-6.81826e-05,1.6523e-08,-1.21599e-12,-57013.7,-589.425], Tmin=(1797.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(50.079,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(C=CC(C)2OJ) + radical(CJC(C)2O)"""),
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
    label = '[CH2]C(O)(C(=O)CC)C(C)=C=C(25676)',
    structure = SMILES('[CH2]C(O)(C(=O)CC)C(C)=C=C'),
    E0 = (-69.7858,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3615,1277.5,1000,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1685.35,2830.27,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154112,'amu*angstrom^2'), symmetry=1, barrier=(3.54333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154112,'amu*angstrom^2'), symmetry=1, barrier=(3.54333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154112,'amu*angstrom^2'), symmetry=1, barrier=(3.54333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154112,'amu*angstrom^2'), symmetry=1, barrier=(3.54333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154112,'amu*angstrom^2'), symmetry=1, barrier=(3.54333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154112,'amu*angstrom^2'), symmetry=1, barrier=(3.54333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154112,'amu*angstrom^2'), symmetry=1, barrier=(3.54333,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (153.198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.51151,0.15348,-0.000221029,1.79095e-07,-5.83116e-11,-8168.61,41.4521], Tmin=(100,'K'), Tmax=(816.861,'K')), NASAPolynomial(coeffs=[14.9812,0.0576363,-2.63292e-05,4.93098e-09,-3.37185e-13,-10686.6,-37.3218], Tmin=(816.861,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-69.7858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CC(O)(C=O)CJ)"""),
)

species(
    label = 'C#CC(C)C([CH2])(O)C(=O)CC(25677)',
    structure = SMILES('C#CC(C)C([CH2])(O)C(=O)CC'),
    E0 = (-62.3691,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,2175,525,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,750,770,3400,2100,180,180,180,234.619,1510.47,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.149432,'amu*angstrom^2'), symmetry=1, barrier=(3.43575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149432,'amu*angstrom^2'), symmetry=1, barrier=(3.43575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149432,'amu*angstrom^2'), symmetry=1, barrier=(3.43575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149432,'amu*angstrom^2'), symmetry=1, barrier=(3.43575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149432,'amu*angstrom^2'), symmetry=1, barrier=(3.43575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149432,'amu*angstrom^2'), symmetry=1, barrier=(3.43575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149432,'amu*angstrom^2'), symmetry=1, barrier=(3.43575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149432,'amu*angstrom^2'), symmetry=1, barrier=(3.43575,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (153.198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.51724,0.146351,-0.000183028,1.23135e-07,-3.29602e-11,-7269.26,42.1453], Tmin=(100,'K'), Tmax=(915.057,'K')), NASAPolynomial(coeffs=[20.4319,0.0460287,-1.85675e-05,3.31185e-09,-2.22168e-13,-11469,-66.533], Tmin=(915.057,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-62.3691,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CtCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Ct-CtCs) + group(Ct-CtH) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = 'C=[C][CH]C(18176)',
    structure = SMILES('[CH2][C]=CC'),
    E0 = (361.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.352622,'amu*angstrom^2'), symmetry=1, barrier=(8.10748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.828631,'amu*angstrom^2'), symmetry=1, barrier=(19.0519,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.42015,0.030446,-1.69076e-05,4.64684e-09,-5.12013e-13,43485.7,14.8304], Tmin=(100,'K'), Tmax=(2065.83,'K')), NASAPolynomial(coeffs=[10.7464,0.014324,-5.20136e-06,8.69079e-10,-5.48385e-14,40045.6,-31.3799], Tmin=(2065.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Cds_S) + radical(Allyl_P)"""),
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
    label = 'C=[C]C(C)C(=C)O(25093)',
    structure = SMILES('C=[C]C(C)C(=C)O'),
    E0 = (78.4632,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,350,440,435,1725,191.069,191.255],'cm^-1')),
        HinderedRotor(inertia=(0.537604,'amu*angstrom^2'), symmetry=1, barrier=(13.9898,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.529071,'amu*angstrom^2'), symmetry=1, barrier=(13.9882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.540397,'amu*angstrom^2'), symmetry=1, barrier=(13.9907,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.538676,'amu*angstrom^2'), symmetry=1, barrier=(13.9898,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.199685,0.076551,-7.19768e-05,3.53677e-08,-6.90749e-12,9579.69,26.5827], Tmin=(100,'K'), Tmax=(1242.4,'K')), NASAPolynomial(coeffs=[16.4084,0.0243651,-8.96989e-06,1.55807e-09,-1.04111e-13,5552.19,-55.1329], Tmin=(1242.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(78.4632,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S)"""),
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
    label = '[CH2]C(O)(C=C=C)C(=O)CC(25678)',
    structure = SMILES('[CH2]C(O)(C=C=C)C(=O)CC'),
    E0 = (-30.7308,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1672.99,2841.8,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153306,'amu*angstrom^2'), symmetry=1, barrier=(3.52481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153306,'amu*angstrom^2'), symmetry=1, barrier=(3.52481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153306,'amu*angstrom^2'), symmetry=1, barrier=(3.52481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153306,'amu*angstrom^2'), symmetry=1, barrier=(3.52481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153306,'amu*angstrom^2'), symmetry=1, barrier=(3.52481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153306,'amu*angstrom^2'), symmetry=1, barrier=(3.52481,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (139.172,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.69952,0.133837,-0.000191635,1.53303e-07,-4.93124e-11,-3498.89,37.7815], Tmin=(100,'K'), Tmax=(809.728,'K')), NASAPolynomial(coeffs=[14.2511,0.0487121,-2.22166e-05,4.1621e-09,-2.84867e-13,-5874.49,-34.5231], Tmin=(809.728,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-30.7308,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CC(O)(C=O)CJ)"""),
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
    label = 'C=[C]C(C)C(=C)C(=O)CC(25679)',
    structure = SMILES('C=[C]C(C)C(=C)C(=O)CC'),
    E0 = (91.6442,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.765449,0.114611,-0.000144681,1.23153e-07,-4.50517e-11,11184.5,36.9185], Tmin=(100,'K'), Tmax=(757.153,'K')), NASAPolynomial(coeffs=[5.95775,0.0686101,-3.2781e-05,6.34028e-09,-4.4451e-13,10466.8,8.33755], Tmin=(757.153,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(91.6442,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S)"""),
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
    label = 'C=[C]C(C)C(C)([O])C(=O)CC(25680)',
    structure = SMILES('C=[C]C(C)C(C)([O])C(=O)CC'),
    E0 = (40.8939,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,945.039,1282.25,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.167414,'amu*angstrom^2'), symmetry=1, barrier=(3.84917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167414,'amu*angstrom^2'), symmetry=1, barrier=(3.84917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167414,'amu*angstrom^2'), symmetry=1, barrier=(3.84917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167414,'amu*angstrom^2'), symmetry=1, barrier=(3.84917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167414,'amu*angstrom^2'), symmetry=1, barrier=(3.84917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167414,'amu*angstrom^2'), symmetry=1, barrier=(3.84917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167414,'amu*angstrom^2'), symmetry=1, barrier=(3.84917,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.32214,0.127753,-0.000140924,8.37016e-08,-1.59209e-11,5100.24,43.3723], Tmin=(100,'K'), Tmax=(638.867,'K')), NASAPolynomial(coeffs=[11.3006,0.0644279,-2.91188e-05,5.51356e-09,-3.83179e-13,3166.86,-14.3778], Tmin=(638.867,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(40.8939,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)(C=O)OJ) + radical(Cds_S)"""),
)

species(
    label = 'C=[C][C](C)C(C)(O)C(=O)CC(25681)',
    structure = SMILES('[CH2][C]=C(C)C(C)(O)C(=O)CC'),
    E0 = (-70.2726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.13552,0.143455,-0.000183382,1.35774e-07,-4.13522e-11,-8238.74,41.7621], Tmin=(100,'K'), Tmax=(796.671,'K')), NASAPolynomial(coeffs=[14.3761,0.0605478,-2.72737e-05,5.1336e-09,-3.54642e-13,-10869.5,-34.1431], Tmin=(796.671,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-70.2726,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + radical(Cds_S) + radical(Allyl_P)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(4659.56,'J/mol'), sigma=(7.82797,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=727.81 K, Pc=22.04 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.45095,0.14396,-0.000169003,1.07491e-07,-2.73862e-11,-2443.92,45.8821], Tmin=(100,'K'), Tmax=(957.061,'K')), NASAPolynomial(coeffs=[20.3257,0.0487695,-1.98159e-05,3.57396e-09,-2.42303e-13,-6803.79,-63.003], Tmin=(957.061,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-22.2369,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = '[CH2]C([O])(C(=O)CC)C(C)C=C(19919)',
    structure = SMILES('[CH2]C([O])(C(=O)CC)C(C)C=C'),
    E0 = (14.6635,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1021.97,1203.88,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.165481,'amu*angstrom^2'), symmetry=1, barrier=(3.80474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165481,'amu*angstrom^2'), symmetry=1, barrier=(3.80474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165481,'amu*angstrom^2'), symmetry=1, barrier=(3.80474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165481,'amu*angstrom^2'), symmetry=1, barrier=(3.80474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165481,'amu*angstrom^2'), symmetry=1, barrier=(3.80474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165481,'amu*angstrom^2'), symmetry=1, barrier=(3.80474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165481,'amu*angstrom^2'), symmetry=1, barrier=(3.80474,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
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
    label = 'C=[C]C(C)C(C)(O)C(=O)[CH]C(25682)',
    structure = SMILES('C=[C]C(C)C(C)(O)C([O])=CC'),
    E0 = (-64.6166,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.07946,0.124783,-0.000117615,5.87426e-08,-1.17447e-11,-7545.17,46.4945], Tmin=(100,'K'), Tmax=(1209.6,'K')), NASAPolynomial(coeffs=[22.6574,0.0429826,-1.61783e-05,2.83716e-09,-1.90418e-13,-13529.6,-77.554], Tmin=(1209.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-64.6166,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_S)"""),
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
    label = '[CH]=[C]C(C)C(C)(O)C(=O)CC(25683)',
    structure = SMILES('[CH]=[C]C(C)C(C)(O)C(=O)CC'),
    E0 = (46.0074,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,3615,1277.5,1000,1685,370,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1740.94,2770.93,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155853,'amu*angstrom^2'), symmetry=1, barrier=(3.58337,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155853,'amu*angstrom^2'), symmetry=1, barrier=(3.58337,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155853,'amu*angstrom^2'), symmetry=1, barrier=(3.58337,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155853,'amu*angstrom^2'), symmetry=1, barrier=(3.58337,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155853,'amu*angstrom^2'), symmetry=1, barrier=(3.58337,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155853,'amu*angstrom^2'), symmetry=1, barrier=(3.58337,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155853,'amu*angstrom^2'), symmetry=1, barrier=(3.58337,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155853,'amu*angstrom^2'), symmetry=1, barrier=(3.58337,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.94384,0.137004,-0.000159231,1.04477e-07,-2.814e-11,5742,43.5502], Tmin=(100,'K'), Tmax=(896.883,'K')), NASAPolynomial(coeffs=[16.0525,0.056741,-2.49915e-05,4.69308e-09,-3.25424e-13,2513.95,-41.3126], Tmin=(896.883,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(46.0074,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]CC(=O)C(C)(O)C(C)[C]=C(25684)',
    structure = SMILES('[CH2]CC(=O)C(C)(O)C(C)[C]=C'),
    E0 = (10.5003,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.03753,0.138862,-0.000163051,1.07281e-07,-2.88538e-11,1475.04,44.1163], Tmin=(100,'K'), Tmax=(899.723,'K')), NASAPolynomial(coeffs=[16.6613,0.0557296,-2.44528e-05,4.58255e-09,-3.17358e-13,-1889.67,-44.1185], Tmin=(899.723,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(10.5003,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CJCC=O)"""),
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
    label = 'C=[C]C(C)C1(O)CO[C]1CC(23254)',
    structure = SMILES('C=[C]C(C)C1(O)CO[C]1CC'),
    E0 = (96.258,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.70489,0.120838,-0.000112328,5.57707e-08,-1.00814e-11,11786.8,41.846], Tmin=(100,'K'), Tmax=(909.153,'K')), NASAPolynomial(coeffs=[17.8105,0.0490015,-1.69458e-05,2.79675e-09,-1.80536e-13,7658.7,-53.6334], Tmin=(909.153,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(96.258,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Oxetane) + radical(C2CsJOCs) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C1(O)[C](CC)OC(=C)C1C(25685)',
    structure = SMILES('[CH2]C1(O)[C](CC)OC(=C)C1C'),
    E0 = (-43.7114,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.9623,0.109705,-5.13643e-05,-2.34482e-08,2.00959e-11,-5023.33,40.2858], Tmin=(100,'K'), Tmax=(968.593,'K')), NASAPolynomial(coeffs=[28.3091,0.0343412,-1.1541e-05,2.06224e-09,-1.4758e-13,-13216.4,-116.812], Tmin=(968.593,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-43.7114,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(C2CsJOC(O)) + radical(CJC(C)2O)"""),
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
    label = '[CH2]C(=C=O)C(C)[C]=C(25686)',
    structure = SMILES('C=[C]C(C)C(=C)[C]=O'),
    E0 = (325.82,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,1855,455,950,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,350,440,435,1725,250.113,250.113],'cm^-1')),
        HinderedRotor(inertia=(0.27201,'amu*angstrom^2'), symmetry=1, barrier=(12.0749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.27201,'amu*angstrom^2'), symmetry=1, barrier=(12.0749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.27201,'amu*angstrom^2'), symmetry=1, barrier=(12.0749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.27201,'amu*angstrom^2'), symmetry=1, barrier=(12.0749,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00968,0.0763001,-6.78039e-05,5.26222e-09,2.69283e-11,39284.4,26.9193], Tmin=(100,'K'), Tmax=(519.127,'K')), NASAPolynomial(coeffs=[7.06001,0.0451472,-2.24783e-05,4.44558e-09,-3.1653e-13,38447.8,-0.310646], Tmin=(519.127,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(325.82,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)CJ=O) + radical(Cds_S)"""),
)

species(
    label = 'C=C=C(C)C(C)(O)C(=O)CC(25687)',
    structure = SMILES('C=C=C(C)C(C)(O)C(=O)CC'),
    E0 = (-283.009,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.00058,0.140195,-0.000173092,1.23777e-07,-3.65429e-11,-33829.5,39.5812], Tmin=(100,'K'), Tmax=(819.564,'K')), NASAPolynomial(coeffs=[14.3276,0.0605064,-2.72475e-05,5.14511e-09,-3.56767e-13,-36506,-35.9443], Tmin=(819.564,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-283.009,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
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
    label = '[CH2]C(O)(C[C]=C)C(=O)CC(25688)',
    structure = SMILES('[CH2]C(O)(C[C]=C)C(=O)CC'),
    E0 = (42.2748,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1277.5,1000,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,180,180,180,555.264,594.351,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.159434,'amu*angstrom^2'), symmetry=1, barrier=(3.6657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159434,'amu*angstrom^2'), symmetry=1, barrier=(3.6657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159434,'amu*angstrom^2'), symmetry=1, barrier=(3.6657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159434,'amu*angstrom^2'), symmetry=1, barrier=(3.6657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159434,'amu*angstrom^2'), symmetry=1, barrier=(3.6657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159434,'amu*angstrom^2'), symmetry=1, barrier=(3.6657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159434,'amu*angstrom^2'), symmetry=1, barrier=(3.6657,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.80138,0.134055,-0.000175888,1.26558e-07,-3.66325e-11,5287.53,40.2469], Tmin=(100,'K'), Tmax=(843.587,'K')), NASAPolynomial(coeffs=[16.6162,0.0467241,-2.06022e-05,3.83739e-09,-2.63506e-13,2180.2,-45.4743], Tmin=(843.587,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(42.2748,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJC(C)(C=O)O) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(O)(C(C)=O)C(C)[C]=C(25689)',
    structure = SMILES('[CH2]C(O)(C(C)=O)C(C)[C]=C'),
    E0 = (35.9395,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,375,552.5,462.5,1710,180,180,180,180,1600,1702.54,2802.55,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154456,'amu*angstrom^2'), symmetry=1, barrier=(3.55125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154456,'amu*angstrom^2'), symmetry=1, barrier=(3.55125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154456,'amu*angstrom^2'), symmetry=1, barrier=(3.55125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154456,'amu*angstrom^2'), symmetry=1, barrier=(3.55125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154456,'amu*angstrom^2'), symmetry=1, barrier=(3.55125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154456,'amu*angstrom^2'), symmetry=1, barrier=(3.55125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154456,'amu*angstrom^2'), symmetry=1, barrier=(3.55125,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.61264,0.123429,-0.000137858,8.16512e-08,-1.93219e-11,4524.76,40.2315], Tmin=(100,'K'), Tmax=(1028.18,'K')), NASAPolynomial(coeffs=[19.7474,0.0403307,-1.66278e-05,3.04631e-09,-2.09292e-13,132.356,-63.4119], Tmin=(1028.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(35.9395,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = 'C=[C]C(C)CC(O)=C([O])CC(19870)',
    structure = SMILES('C=[C]C(C)CC(O)=C([O])CC'),
    E0 = (-80.1938,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.24914,0.123359,-0.000103514,3.8135e-08,-3.33349e-12,-9407.81,46.277], Tmin=(100,'K'), Tmax=(1022.17,'K')), NASAPolynomial(coeffs=[25.4916,0.0385609,-1.39397e-05,2.45249e-09,-1.67774e-13,-16320.1,-94.2354], Tmin=(1022.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-80.1938,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=[C]C(C)[C](O)CC(=O)CC(25690)',
    structure = SMILES('C=[C]C(C)[C](O)CC(=O)CC'),
    E0 = (-16.9223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.97851,0.14297,-0.000199806,1.68228e-07,-5.78118e-11,-1831.02,45.2365], Tmin=(100,'K'), Tmax=(808.725,'K')), NASAPolynomial(coeffs=[10.2073,0.0674027,-3.12769e-05,5.91714e-09,-4.07569e-13,-3301.84,-7.87378], Tmin=(808.725,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-16.9223,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)',
    structure = SMILES('[CH2]C(=CC)C([CH2])(O)C(=O)CC'),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(4732.62,'J/mol'), sigma=(7.88591,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=739.22 K, Pc=21.9 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.39668,0.149333,-0.000195449,1.45429e-07,-4.41529e-11,-11190.4,42.2789], Tmin=(100,'K'), Tmax=(801.412,'K')), NASAPolynomial(coeffs=[15.6661,0.0591689,-2.66711e-05,5.01446e-09,-3.45918e-13,-14085.2,-40.863], Tmin=(801.412,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-94.891,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(C=CC(O)(C=O)CJ)"""),
)

species(
    label = '[CH2]C(O)([CH]C(=C)C)C(=O)CC(25691)',
    structure = SMILES('[CH2]C(O)([CH]C(=C)C)C(=O)CC'),
    E0 = (-117.705,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.63936,0.151091,-0.00018095,1.16709e-07,-3.02966e-11,-13921.9,40.7416], Tmin=(100,'K'), Tmax=(936.619,'K')), NASAPolynomial(coeffs=[20.4457,0.0525033,-2.30618e-05,4.32836e-09,-3.00336e-13,-18246.3,-69.1188], Tmin=(936.619,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-117.705,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsHH) + radical(CJC(C)(C=O)O) + radical(C=CCJCO)"""),
)

species(
    label = 'C=C1CC(O)(C(=O)CC)C1C(25692)',
    structure = SMILES('C=C1CC(O)(C(=O)CC)C1C'),
    E0 = (-298.275,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[9.13131,0.0411926,7.21588e-05,-9.3506e-08,2.21474e-11,-36194.3,-8.44765], Tmin=(100,'K'), Tmax=(1796.01,'K')), NASAPolynomial(coeffs=[99.7485,0.0253111,-6.98687e-05,1.68567e-08,-1.2385e-12,-98732.7,-582.172], Tmin=(1796.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-298.275,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane)"""),
)

species(
    label = '[CH2]C(=C(CC)OO)C(C)[C]=C(25693)',
    structure = SMILES('[CH2]C(=C(CC)OO)C(C)[C]=C'),
    E0 = (256.731,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.47772,0.127623,-0.00013252,7.89896e-08,-1.98506e-11,31068.9,42.4265], Tmin=(100,'K'), Tmax=(944.659,'K')), NASAPolynomial(coeffs=[14.027,0.0619714,-2.82742e-05,5.42172e-09,-3.81285e-13,28139.5,-31.4921], Tmin=(944.659,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(256.731,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(O)=C(CC)OC(C)[C]=C(19945)',
    structure = SMILES('[CH2]C(O)=C(CC)OC(C)[C]=C'),
    E0 = (-36.5272,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.60648,0.126789,-9.45921e-05,1.62281e-08,7.31419e-12,-4138.73,45.8654], Tmin=(100,'K'), Tmax=(971.629,'K')), NASAPolynomial(coeffs=[30.4892,0.0314319,-1.05068e-05,1.84798e-09,-1.30356e-13,-12500.2,-122.782], Tmin=(971.629,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-36.5272,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(O)CJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(O)(C(=C)OC)C(C)[C]=C(25694)',
    structure = SMILES('[CH2]C(O)(C(=C)OC)C(C)[C]=C'),
    E0 = (65.7318,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,180,180,180,180,1013.1,1200,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.159983,'amu*angstrom^2'), symmetry=1, barrier=(3.67832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159983,'amu*angstrom^2'), symmetry=1, barrier=(3.67832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159983,'amu*angstrom^2'), symmetry=1, barrier=(3.67832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159983,'amu*angstrom^2'), symmetry=1, barrier=(3.67832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159983,'amu*angstrom^2'), symmetry=1, barrier=(3.67832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159983,'amu*angstrom^2'), symmetry=1, barrier=(3.67832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159983,'amu*angstrom^2'), symmetry=1, barrier=(3.67832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159983,'amu*angstrom^2'), symmetry=1, barrier=(3.67832,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.5416,0.132059,-0.000128713,6.47207e-08,-1.28567e-11,8151.21,47.3578], Tmin=(100,'K'), Tmax=(1225.96,'K')), NASAPolynomial(coeffs=[26.1245,0.0385302,-1.42793e-05,2.49342e-09,-1.6737e-13,1122.42,-96.7799], Tmin=(1225.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(65.7318,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]C(O)(C(O)=CC)C(C)[C]=C(25695)',
    structure = SMILES('[CH2]C(O)(C(O)=CC)C(C)[C]=C'),
    E0 = (11.022,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,180,180,180,180,1600,1661.6,2837.62,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153399,'amu*angstrom^2'), symmetry=1, barrier=(3.52694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153399,'amu*angstrom^2'), symmetry=1, barrier=(3.52694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153399,'amu*angstrom^2'), symmetry=1, barrier=(3.52694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153399,'amu*angstrom^2'), symmetry=1, barrier=(3.52694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153399,'amu*angstrom^2'), symmetry=1, barrier=(3.52694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153399,'amu*angstrom^2'), symmetry=1, barrier=(3.52694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153399,'amu*angstrom^2'), symmetry=1, barrier=(3.52694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153399,'amu*angstrom^2'), symmetry=1, barrier=(3.52694,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.61621,0.135288,-0.000137649,7.24423e-08,-1.50272e-11,1572.57,47.7949], Tmin=(100,'K'), Tmax=(1177.66,'K')), NASAPolynomial(coeffs=[25.979,0.0381634,-1.39416e-05,2.41297e-09,-1.61192e-13,-5162.59,-94.8369], Tmin=(1177.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(11.022,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CC(C)(O)CJ)"""),
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
    label = 'C=[C]C(C)[C](O)C(=O)CC(25696)',
    structure = SMILES('C=[C]C(C)C(O)=C([O])CC'),
    E0 = (-56.9616,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.62689,0.117562,-0.000121053,6.52102e-08,-1.3929e-11,-6642.98,39.4145], Tmin=(100,'K'), Tmax=(1140.7,'K')), NASAPolynomial(coeffs=[21.6829,0.0358231,-1.35671e-05,2.39086e-09,-1.61166e-13,-11960.9,-76.1105], Tmin=(1140.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-56.9616,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=C(C)OJ)"""),
)

species(
    label = 'H2CC(41)',
    structure = SMILES('[C]=C'),
    E0 = (401.202,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2480.69,'J/mol'), sigma=(4.48499,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=387.48 K, Pc=62.39 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28155,0.00697643,-2.38528e-06,-1.21078e-09,9.82042e-13,48319.2,5.92036], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.27807,0.00475623,-1.63007e-06,2.54623e-10,-1.4886e-14,48014,0.639979], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(401.202,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""H2CC""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[CH2]C(O)([CH]C)C(=O)CC(6379)',
    structure = SMILES('[CH2]C(O)([CH]C)C(=O)CC'),
    E0 = (-99.667,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1754.16,2755.36,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156202,'amu*angstrom^2'), symmetry=1, barrier=(3.5914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156202,'amu*angstrom^2'), symmetry=1, barrier=(3.5914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156202,'amu*angstrom^2'), symmetry=1, barrier=(3.5914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156202,'amu*angstrom^2'), symmetry=1, barrier=(3.5914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156202,'amu*angstrom^2'), symmetry=1, barrier=(3.5914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156202,'amu*angstrom^2'), symmetry=1, barrier=(3.5914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156202,'amu*angstrom^2'), symmetry=1, barrier=(3.5914,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (128.169,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.16691,0.116456,-0.000135482,8.5667e-08,-2.18244e-11,-11803.5,37.2714], Tmin=(100,'K'), Tmax=(953.997,'K')), NASAPolynomial(coeffs=[16.727,0.0414301,-1.75175e-05,3.23263e-09,-2.22295e-13,-15217.7,-48.2139], Tmin=(953.997,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-99.667,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CCJCO) + radical(CJC(C)(C=O)O)"""),
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
    E0 = (10.5226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (105.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (90.9999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (153.847,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (162.573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (18.1955,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (128.207,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (10.5226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (137.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (134.421,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (193.104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (125.214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (179.572,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (128.462,'kJ/mol'),
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
    E0 = (160.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (93.3658,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (97.1314,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (46.6096,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (79.0476,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (63.241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (76.5949,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (176.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (196.948,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (61.4628,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (406.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (88.7697,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (462.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (455.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (167.841,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (167.841,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (104.997,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (155.143,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (18.8069,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (399.845,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (106.586,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (324.721,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (154.136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (324.601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (301.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    products = ['C=C(O)C(=O)CC(4626)', 'CH3CHCCH2(18175)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    products = ['C=[C]C(C)C1(O)CC1([O])CC(25674)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.34238e+09,'s^-1'), n=0.889391, Ea=(94.5285,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonylbond_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 90.6 to 94.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    products = ['[CH2]C1(O)C(C)C(=C)C1([O])CC(25675)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.98674e+07,'s^-1'), n=1.31443, Ea=(80.4773,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra] for rate rule [R5_SS_CO;carbonylbond_intra_Nd;radadd_intra_cddouble]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[CH2]C(O)(C(=O)CC)C(C)=C=C(25676)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.51e+07,'cm^3/(mol*s)'), n=1.64, Ea=(11.8407,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2579 used for Cds-CsCs_Ca;HJ
Exact match found for rate rule [Cds-CsCs_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C#CC(C)C([CH2])(O)C(=O)CC(25677)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.255e+11,'cm^3/(mol*s)'), n=1.005, Ea=(13.1503,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 138 used for Ct-H_Ct-Cs;HJ
Exact match found for rate rule [Ct-H_Ct-Cs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=C(O)C(=O)CC(4626)', 'C=[C][CH]C(18176)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.00392164,'m^3/(mol*s)'), n=2.41519, Ea=(15.6067,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;CJ] for rate rule [Cds-COOs_Cds;CJ]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C2H5CO(71)', 'C=[C]C(C)C(=C)O(25093)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=0, Ea=(93.9308,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_R;CO_rad] for rate rule [Cds-OsCs_Cds;CO_rad/NonDe]
Euclidian distance = 2.82842712475
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C(O)=C([O])CC(4557)', 'CH3CHCCH2(18175)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.00472174,'m^3/(mol*s)'), n=2.41, Ea=(49.8505,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Ca;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 48.4 to 49.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['CH3(17)', '[CH2]C(O)(C=C=C)C(=O)CC(25678)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(10800,'cm^3/(mol*s)'), n=2.41, Ea=(32.1331,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 597 used for Cds-CsH_Ca;CsJ-HHH
Exact match found for rate rule [Cds-CsH_Ca;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['OH(5)', 'C=[C]C(C)C(=C)C(=O)CC(25679)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(0.52832,'m^3/(mol*s)'), n=2.02802, Ea=(14.4047,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDeCs_Cds;YJ] for rate rule [Cds-COCs_Cds;OJ_pri]
Euclidian distance = 2.2360679775
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    products = ['[CH2]C(O)([C](C)C=C)C(=O)CC(19915)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_noH] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_Cs2]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=CC(C)C([CH2])(O)C(=O)CC(19922)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=[C]C(C)C(C)([O])C(=O)CC(25680)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    products = ['C=[C][C](C)C(C)(O)C(=O)CC(25681)'],
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
    reactants = ['[CH2]C([O])(C(=O)CC)C(C)C=C(19919)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(274,'s^-1'), n=3.09, Ea=(145.603,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    products = ['[CH2]C([C]=C)C(C)(O)C(=O)CC(19920)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(114000,'s^-1'), n=1.74, Ea=(82.8432,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 109 used for R4H_SSS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    products = ['C=[C]C(C)C(C)(O)C(=O)[CH]C(25682)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(6.44e+09,'s^-1'), n=0.13, Ea=(86.6088,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 131 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    products = ['[CH2]C(O)(C(=O)[CH]C)C(C)C=C(19923)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.44776e+08,'s^-1'), n=0.81, Ea=(36.087,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5H_CC(O2d)CC;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 3.16227766017
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=[C]C(C)C(C)(O)C(=O)CC(25683)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5HJ_1;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    products = ['[CH2]CC(=O)C(C)(O)C(C)[C]=C(25684)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(68850,'s^-1'), n=1.68, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_2H;Cs_H_out_2H] for rate rule [R5H_CCC(O2d)C;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    products = ['[CH2]CC(=O)C([CH2])(O)C(C)C=C(19926)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(60905.8,'s^-1'), n=1.5925, Ea=(66.0723,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_SSSSS;Y_rad_out;Cs_H_out_2H] for rate rule [R6H_SSSSS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C(O)=C([O])CC(4557)', 'C=[C][CH]C(18176)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    products = ['C=[C]C(C)C1(O)CO[C]1CC(23254)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(7.01137e+09,'s^-1'), n=0.572544, Ea=(186.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4_S;multiplebond_intra;radadd_intra_cs2H] + [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    products = ['[CH2]C1(O)[C](CC)OC(=C)C1C(25685)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.67521e+10,'s^-1'), n=0.355, Ea=(50.9402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cddouble] for rate rule [R5_SS_CO;carbonyl_intra_Nd;radadd_intra_cddouble]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CH3CH2OH(54)', '[CH2]C(=C=O)C(C)[C]=C(25686)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.79e-11,'m^3/(mol*s)'), n=3.97, Ea=(329.281,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [doublebond;Cs_OH] for rate rule [Cdd_Cd;C_pri_OH]
Euclidian distance = 1.41421356237
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    products = ['C=C=C(C)C(C)(O)C(=O)CC(25687)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CH2(S)(23)', '[CH2]C(O)(C[C]=C)C(=O)CC(25688)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction29',
    reactants = ['CH2(S)(23)', '[CH2]C(O)(C(C)=O)C(C)[C]=C(25689)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;C_pri/De]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    products = ['C=[C]C(C)CC(O)=C([O])CC(19870)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    products = ['C=[C]C(C)[C](O)CC(=O)CC(25690)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ-HH;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;CO]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    products = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 5 used for cCs(-HC)CJ;CdsJ;C
Exact match found for rate rule [cCs(-HC)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    products = ['[CH2]C(O)([CH]C(=C)C)C(=O)CC(25691)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(6.95888e+10,'s^-1'), n=0.7315, Ea=(144.62,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCs(-HC)CJ;CJ;CH3] + [cCs(-HC)CJ;CdsJ;C] for rate rule [cCs(-HC)CJ;CdsJ;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    products = ['C=C1CC(O)(C(=O)CC)C1C(25692)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C(=C(CC)OO)C(C)[C]=C(25693)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O] for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O_R]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C(O)=C(CC)OC(C)[C]=C(19945)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O] for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O_C]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C(O)(C(=C)OC)C(C)[C]=C(25694)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(37989.5,'s^-1'), n=2.515, Ea=(258.99,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C] + [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O] for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_C]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C(O)(C(O)=CC)C(C)[C]=C(25695)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O_H] for rate rule [R_ROR;R1_doublebond_CHCH3;R2_doublebond_CsC;R_O_H]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction39',
    reactants = ['CH2(19)', 'C=[C]C(C)[C](O)C(=O)CC(25696)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/ODMustO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction40',
    reactants = ['H2CC(41)', '[CH2]C(O)([CH]C)C(=O)CC(6379)'],
    products = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4303',
    isomers = [
        '[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)',
    ],
    reactants = [
        ('C=C(O)C(=O)CC(4626)', 'CH3CHCCH2(18175)'),
        ('[CH2]C(O)=C([O])CC(4557)', 'CH3CHCCH2(18175)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4303',
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

