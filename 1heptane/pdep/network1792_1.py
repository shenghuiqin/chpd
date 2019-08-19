species(
    label = '[CH2]C(O)(C(=C)[O])C(=O)CC(8846)',
    structure = SMILES('[CH2]C(O)(C(=C)[O])C(=O)CC'),
    E0 = (-248.063,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,463.975,687.011,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.162293,'amu*angstrom^2'), symmetry=1, barrier=(3.73143,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162293,'amu*angstrom^2'), symmetry=1, barrier=(3.73143,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162293,'amu*angstrom^2'), symmetry=1, barrier=(3.73143,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162293,'amu*angstrom^2'), symmetry=1, barrier=(3.73143,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162293,'amu*angstrom^2'), symmetry=1, barrier=(3.73143,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162293,'amu*angstrom^2'), symmetry=1, barrier=(3.73143,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.78334,0.134036,-0.000194366,1.51022e-07,-4.64553e-11,-29633.3,38.7807], Tmin=(100,'K'), Tmax=(837.148,'K')), NASAPolynomial(coeffs=[16.8055,0.0402935,-1.75764e-05,3.2102e-09,-2.15971e-13,-32573.2,-46.5652], Tmin=(837.148,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-248.063,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CC(O)(C=O)CJ)"""),
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
    label = 'C=C([O])C1(O)CC1([O])CC(9508)',
    structure = SMILES('C=C([O])C1(O)CC1([O])CC'),
    E0 = (-161.405,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.11852,0.0964975,-6.08587e-05,-4.88005e-09,1.3262e-11,-19213.5,37.3388], Tmin=(100,'K'), Tmax=(942.858,'K')), NASAPolynomial(coeffs=[25.506,0.0229628,-6.58168e-06,1.08273e-09,-7.59267e-14,-25986.2,-98.834], Tmin=(942.858,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-161.405,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(CC(C)2OJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C1([O])CC1(O)C(=O)CC(9509)',
    structure = SMILES('[CH2]C1([O])CC1(O)C(=O)CC'),
    E0 = (-98.1059,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.37471,0.118528,-0.000136589,8.2202e-08,-1.96537e-11,-11606,34.7862], Tmin=(100,'K'), Tmax=(1020.21,'K')), NASAPolynomial(coeffs=[19.7073,0.0358699,-1.50563e-05,2.78416e-09,-1.92341e-13,-15907.6,-67.3441], Tmin=(1020.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-98.1059,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclopropane) + radical(CC(C)2OJ) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH2]C1(O)C(=C)OC1([O])CC(9487)',
    structure = SMILES('[CH2]C1(O)C(=C)OC1([O])CC'),
    E0 = (-103.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.24596,0.0929469,-3.06201e-05,-5.16118e-08,3.44474e-11,-12274.5,32.9656], Tmin=(100,'K'), Tmax=(907.619,'K')), NASAPolynomial(coeffs=[30.4149,0.0140412,-4.13096e-07,-2.01565e-10,1.44416e-14,-20518.9,-130.468], Tmin=(907.619,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-103.801,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(C=CC(C)(O)CJ) + radical(CC(C)(O)OJ)"""),
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
    label = 'C=C([O])C(=C)O(4628)',
    structure = SMILES('C=C([O])C(=C)O'),
    E0 = (-196.06,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.916777,'amu*angstrom^2'), symmetry=1, barrier=(21.0785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.920477,'amu*angstrom^2'), symmetry=1, barrier=(21.1636,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.284274,0.0664439,-7.46756e-05,3.95333e-08,-7.70645e-12,-23433.6,21.3173], Tmin=(100,'K'), Tmax=(1467.22,'K')), NASAPolynomial(coeffs=[18.8952,0.00393949,1.25465e-06,-4.33245e-10,3.4763e-14,-27628.4,-71.2884], Tmin=(1467.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-196.06,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=[C][O](173)',
    structure = SMILES('[CH2][C]=O'),
    E0 = (160.185,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,539.612,539.669],'cm^-1')),
        HinderedRotor(inertia=(0.000578908,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39563,0.0101365,2.30741e-06,-8.97566e-09,3.68242e-12,19290.3,10.0703], Tmin=(100,'K'), Tmax=(1068.9,'K')), NASAPolynomial(coeffs=[6.35055,0.00638951,-2.69368e-06,5.4221e-10,-4.02476e-14,18240.9,-6.33602], Tmin=(1068.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CsCJ=O) + radical(CJC=O)"""),
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
    label = 'C=C([O])C(=C)C(=O)CC(9510)',
    structure = SMILES('C=C([O])C(=C)C(=O)CC'),
    E0 = (-166.692,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,375,552.5,462.5,1710,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,351.27,351.27,351.271,351.271],'cm^-1')),
        HinderedRotor(inertia=(0.13069,'amu*angstrom^2'), symmetry=1, barrier=(11.4433,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13069,'amu*angstrom^2'), symmetry=1, barrier=(11.4433,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13069,'amu*angstrom^2'), symmetry=1, barrier=(11.4433,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13069,'amu*angstrom^2'), symmetry=1, barrier=(11.4433,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.234966,0.093935,-9.61713e-05,5.27534e-08,-1.17633e-11,-19896.5,29.4879], Tmin=(100,'K'), Tmax=(1077.1,'K')), NASAPolynomial(coeffs=[15.3185,0.036174,-1.57312e-05,2.96503e-09,-2.07065e-13,-23247,-46.7038], Tmin=(1077.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-166.692,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])C(C)([O])C(=O)CC(9511)',
    structure = SMILES('C=C([O])C(C)([O])C(=O)CC'),
    E0 = (-217,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,375,552.5,462.5,1710,180,180,180,180,859.178,1365.91,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.174744,'amu*angstrom^2'), symmetry=1, barrier=(4.0177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174744,'amu*angstrom^2'), symmetry=1, barrier=(4.0177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174744,'amu*angstrom^2'), symmetry=1, barrier=(4.0177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174744,'amu*angstrom^2'), symmetry=1, barrier=(4.0177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174744,'amu*angstrom^2'), symmetry=1, barrier=(4.0177,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.29625,0.126682,-0.000190379,1.60373e-07,-5.34463e-11,-25918.1,37.967], Tmin=(100,'K'), Tmax=(847.128,'K')), NASAPolynomial(coeffs=[11.5512,0.0488627,-2.22085e-05,4.12151e-09,-2.79097e-13,-27479.3,-18.25], Tmin=(847.128,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-217,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CC(C)(C=O)OJ)"""),
)

species(
    label = '[CH]=C(O)C([CH2])(O)C(=O)CC(9512)',
    structure = SMILES('[CH]=C(O)C([CH2])(O)C(=O)CC'),
    E0 = (-138.772,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.17612,0.140909,-0.000200956,1.47366e-07,-4.23166e-11,-16472.7,39.0966], Tmin=(100,'K'), Tmax=(858.659,'K')), NASAPolynomial(coeffs=[20.7291,0.0342039,-1.45483e-05,2.63494e-09,-1.76974e-13,-20406.1,-67.917], Tmin=(858.659,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-138.772,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(O)(C=O)CJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C([O])(C(=C)O)C(=O)CC(9513)',
    structure = SMILES('[CH2]C([O])(C(=C)O)C(=O)CC'),
    E0 = (-141.582,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,464.529,688.273,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.16237,'amu*angstrom^2'), symmetry=1, barrier=(3.7332,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16237,'amu*angstrom^2'), symmetry=1, barrier=(3.7332,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16237,'amu*angstrom^2'), symmetry=1, barrier=(3.7332,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16237,'amu*angstrom^2'), symmetry=1, barrier=(3.7332,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16237,'amu*angstrom^2'), symmetry=1, barrier=(3.7332,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16237,'amu*angstrom^2'), symmetry=1, barrier=(3.7332,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.94863,0.141019,-0.000219198,1.81307e-07,-5.86705e-11,-16823.8,38.8618], Tmin=(100,'K'), Tmax=(855.543,'K')), NASAPolynomial(coeffs=[15.3946,0.0436182,-1.98234e-05,3.65868e-09,-2.46192e-13,-19194.3,-38.6139], Tmin=(855.543,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-141.582,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)(C=O)OJ) + radical(C=CC(O)(C=O)CJ)"""),
)

species(
    label = 'C=C([O])C(C)(O)C(=O)[CH]C(9514)',
    structure = SMILES('C=C([O])C(C)(O)C([O])=CC'),
    E0 = (-329.886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.25894,0.105577,-0.000103756,5.21823e-08,-1.03288e-11,-39478.2,39.6261], Tmin=(100,'K'), Tmax=(1232.77,'K')), NASAPolynomial(coeffs=[22.3046,0.0291191,-1.07225e-05,1.87062e-09,-1.25654e-13,-45287.9,-78.985], Tmin=(1232.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-329.886,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C([O])C(C)(O)C(=O)CC(9515)',
    structure = SMILES('[CH]=C([O])C(C)(O)C(=O)CC'),
    E0 = (-214.191,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1851.99,2662.34,3200],'cm^-1')),
        HinderedRotor(inertia=(0.16001,'amu*angstrom^2'), symmetry=1, barrier=(3.67894,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16001,'amu*angstrom^2'), symmetry=1, barrier=(3.67894,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16001,'amu*angstrom^2'), symmetry=1, barrier=(3.67894,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16001,'amu*angstrom^2'), symmetry=1, barrier=(3.67894,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16001,'amu*angstrom^2'), symmetry=1, barrier=(3.67894,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16001,'amu*angstrom^2'), symmetry=1, barrier=(3.67894,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.49595,0.126209,-0.000170696,1.24282e-07,-3.60329e-11,-25568.2,38.1043], Tmin=(100,'K'), Tmax=(845.625,'K')), NASAPolynomial(coeffs=[16.8213,0.0395647,-1.70033e-05,3.11477e-09,-2.11321e-13,-28666.1,-47.1946], Tmin=(845.625,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-214.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C(O)(C(=C)O)C(=O)[CH]C(9516)',
    structure = SMILES('[CH2]C(O)(C(=C)O)C([O])=CC'),
    E0 = (-254.247,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.7888,0.11601,-0.000123571,6.56439e-08,-1.35266e-11,-30360.8,40.9012], Tmin=(100,'K'), Tmax=(1194.97,'K')), NASAPolynomial(coeffs=[25.5939,0.0243486,-8.5114e-06,1.45204e-09,-9.68676e-14,-36905,-96.0817], Tmin=(1194.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-254.247,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)(O)CJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]CC(=O)C(C)(O)C(=C)[O](9517)',
    structure = SMILES('[CH2]CC(=O)C(C)(O)C(=C)[O]'),
    E0 = (-249.698,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.58629,0.128024,-0.000174348,1.26837e-07,-3.66267e-11,-29835.3,38.6586], Tmin=(100,'K'), Tmax=(849.983,'K')), NASAPolynomial(coeffs=[17.4219,0.038568,-1.64734e-05,3.00636e-09,-2.03433e-13,-33066.4,-49.9546], Tmin=(849.983,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-249.698,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CJCC=O) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]CC(=O)C([CH2])(O)C(=C)O(9518)',
    structure = SMILES('[CH2]CC(=O)C([CH2])(O)C(=C)O'),
    E0 = (-174.279,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,375,552.5,462.5,1710,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,180,180,1600,1686.66,2820.57,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153989,'amu*angstrom^2'), symmetry=1, barrier=(3.54051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153989,'amu*angstrom^2'), symmetry=1, barrier=(3.54051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153989,'amu*angstrom^2'), symmetry=1, barrier=(3.54051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153989,'amu*angstrom^2'), symmetry=1, barrier=(3.54051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153989,'amu*angstrom^2'), symmetry=1, barrier=(3.54051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153989,'amu*angstrom^2'), symmetry=1, barrier=(3.54051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153989,'amu*angstrom^2'), symmetry=1, barrier=(3.54051,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.26795,0.142742,-0.000204677,1.5002e-07,-4.29569e-11,-20739.7,39.6562], Tmin=(100,'K'), Tmax=(861.854,'K')), NASAPolynomial(coeffs=[21.3312,0.0332048,-1.40171e-05,2.52621e-09,-1.6906e-13,-24807.1,-70.6853], Tmin=(861.854,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-174.279,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(O)(C=O)CJ) + radical(CJCC=O)"""),
)

species(
    label = '[CH2]C(O)([C]1CO1)C(=O)CC(9519)',
    structure = SMILES('[CH2]C(O)([C]1CO1)C(=O)CC'),
    E0 = (-102.371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.92724,0.134325,-0.000190858,1.42361e-07,-4.08085e-11,-12102.7,38.0024], Tmin=(100,'K'), Tmax=(973.302,'K')), NASAPolynomial(coeffs=[18.1087,0.0350369,-1.17251e-05,1.7752e-09,-1.03187e-13,-15200.2,-53.9937], Tmin=(973.302,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-102.371,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Ethylene_oxide) + radical(C2CsJO) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = 'C=C([O])C1(O)CO[C]1CC(9520)',
    structure = SMILES('C=C([O])C1(O)CO[C]1CC'),
    E0 = (-170.198,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.41549,0.105646,-0.000108904,5.9921e-08,-1.26966e-11,-20263.5,39.0602], Tmin=(100,'K'), Tmax=(1298.85,'K')), NASAPolynomial(coeffs=[20.2409,0.0278698,-6.28561e-06,6.80502e-10,-2.971e-14,-24954.5,-67.4839], Tmin=(1298.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-170.198,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Oxetane) + radical(C=C(C)OJ) + radical(C2CsJOCs)"""),
)

species(
    label = 'CC[C]([O])C1(O)CCC1=O(8973)',
    structure = SMILES('CC[C]([O])C1(O)CCC1=O'),
    E0 = (-132.761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.909053,0.0989668,-8.81924e-05,4.07324e-08,-7.52729e-12,-15783.3,36.5104], Tmin=(100,'K'), Tmax=(1301.39,'K')), NASAPolynomial(coeffs=[20.2128,0.0340443,-1.33602e-05,2.39698e-09,-1.62809e-13,-21280.7,-70.9541], Tmin=(1301.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-132.761,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(CC(C)OJ) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C1(O)[C](CC)OOC1=C(9470)',
    structure = SMILES('[CH2]C1(O)[C](CC)OOC1=C'),
    E0 = (97.395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.506628,0.0860661,-5.32974e-05,8.4671e-09,2.57204e-12,11887.2,37.7654], Tmin=(100,'K'), Tmax=(1091.39,'K')), NASAPolynomial(coeffs=[18.7209,0.0368039,-1.47392e-05,2.71853e-09,-1.895e-13,6427.22,-62.4642], Tmin=(1091.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(97.395,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(C=CC(C)(O)CJ) + radical(C2CsJOOC)"""),
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
    label = '[CH2]C(=C=O)C(=C)[O](9521)',
    structure = SMILES('C=C([O])C(=C)[C]=O'),
    E0 = (67.4834,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,1855,455,950,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.952057,'amu*angstrom^2'), symmetry=1, barrier=(21.8897,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.954388,'amu*angstrom^2'), symmetry=1, barrier=(21.9433,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.645614,0.0683587,-7.79362e-05,4.26589e-08,-8.97899e-12,8241.81,22.5795], Tmin=(100,'K'), Tmax=(1172.93,'K')), NASAPolynomial(coeffs=[17.5333,0.0107662,-4.28294e-06,7.95362e-10,-5.59939e-14,4280.25,-61.5873], Tmin=(1172.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(67.4834,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)CJ=O) + radical(C=C(C)OJ)"""),
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
    label = '[CH2]C(O)(C(=C)[O])C(C)=O(9522)',
    structure = SMILES('[CH2]C(O)(C(=C)[O])C(C)=O'),
    E0 = (-222.646,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,375,552.5,462.5,1710,180,180,180,180,1600,1810.43,2702.18,3200],'cm^-1')),
        HinderedRotor(inertia=(0.158877,'amu*angstrom^2'), symmetry=1, barrier=(3.65289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158877,'amu*angstrom^2'), symmetry=1, barrier=(3.65289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158877,'amu*angstrom^2'), symmetry=1, barrier=(3.65289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158877,'amu*angstrom^2'), symmetry=1, barrier=(3.65289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158877,'amu*angstrom^2'), symmetry=1, barrier=(3.65289,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.866419,0.110156,-0.000148855,1.04788e-07,-2.90144e-11,-26605.8,34.1626], Tmin=(100,'K'), Tmax=(888.646,'K')), NASAPolynomial(coeffs=[17.1345,0.0291282,-1.20799e-05,2.17635e-09,-1.46402e-13,-29805,-50.556], Tmin=(888.646,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-222.646,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(O)(C=O)CJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])[C](O)CC(=O)CC(9523)',
    structure = SMILES('[CH2]C([O])=C(O)CC(=O)CC'),
    E0 = (-341.848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.20529,0.108585,-0.000107554,5.38034e-08,-1.06678e-11,-40922.4,37.3843], Tmin=(100,'K'), Tmax=(1220.82,'K')), NASAPolynomial(coeffs=[22.0897,0.0322607,-1.37789e-05,2.59564e-09,-1.81758e-13,-46610.4,-79.6491], Tmin=(1220.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-341.848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-OdCsCs) + group(Cds-CdsCsOs) + radical(C=C(O)CJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])C[C](O)C(=O)CC(8845)',
    structure = SMILES('C=C([O])CC(O)=C([O])CC'),
    E0 = (-337.987,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.16894,0.113535,-0.000115587,5.78171e-08,-1.09746e-11,-40409.5,42.8655], Tmin=(100,'K'), Tmax=(1430.46,'K')), NASAPolynomial(coeffs=[28.2173,0.0177505,-3.80482e-06,4.35562e-10,-2.23884e-14,-47996.2,-110.741], Tmin=(1430.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-337.987,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C1OCC1(O)C(=O)CC(9371)',
    structure = SMILES('C=C1OCC1(O)C(=O)CC'),
    E0 = (-423.249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.768,0.110759,-0.000113043,5.8893e-08,-1.18249e-11,-50683.6,31.946], Tmin=(100,'K'), Tmax=(1316.9,'K')), NASAPolynomial(coeffs=[24.55,0.0239072,-6.24141e-06,8.39899e-10,-4.74329e-14,-57015.9,-99.9924], Tmin=(1316.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-423.249,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane)"""),
)

species(
    label = '[CH2][C](C(=C)OO)C(=O)CC(9524)',
    structure = SMILES('[CH2]C(C(=C)OO)=C([O])CC'),
    E0 = (-16.9931,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,325,350,375,415,440,465,420,435,450,1700,1725,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.33022,0.113534,-0.00012077,6.64787e-08,-1.45126e-11,-1848.65,37.5733], Tmin=(100,'K'), Tmax=(1115.17,'K')), NASAPolynomial(coeffs=[20.8476,0.0339852,-1.37708e-05,2.51295e-09,-1.7277e-13,-6795.07,-71.8393], Tmin=(1115.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-16.9931,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C(C(=C)[O])=C(CC)OO(9525)',
    structure = SMILES('[CH2]C(C(=C)[O])=C(CC)OO'),
    E0 = (-16.9931,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,325,350,375,415,440,465,420,435,450,1700,1725,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.33022,0.113534,-0.00012077,6.64787e-08,-1.45126e-11,-1848.65,37.5733], Tmin=(100,'K'), Tmax=(1115.17,'K')), NASAPolynomial(coeffs=[20.8476,0.0339852,-1.37708e-05,2.51295e-09,-1.7277e-13,-6795.07,-71.8393], Tmin=(1115.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-16.9931,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = 'C=C([O])O[C](CC)C(=C)O(8844)',
    structure = SMILES('[CH2]C(=O)OC(CC)=C([CH2])O'),
    E0 = (-294.584,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,3615,1277.5,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
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
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.49281,0.112045,-0.000114954,5.96899e-08,-1.21638e-11,-35225.1,39.9511], Tmin=(100,'K'), Tmax=(1199.16,'K')), NASAPolynomial(coeffs=[23.3791,0.0290818,-1.11781e-05,1.99683e-09,-1.36078e-13,-41190.3,-84.5591], Tmin=(1199.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-294.584,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsOs) + radical(C=C(O)CJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C(O)(C(=C)[O])C(=C)OC(9526)',
    structure = SMILES('[CH2]C(O)(C(=C)[O])C(=C)OC'),
    E0 = (-199.537,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,180,1057.92,1155.05,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.162564,'amu*angstrom^2'), symmetry=1, barrier=(3.73766,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162564,'amu*angstrom^2'), symmetry=1, barrier=(3.73766,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162564,'amu*angstrom^2'), symmetry=1, barrier=(3.73766,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162564,'amu*angstrom^2'), symmetry=1, barrier=(3.73766,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162564,'amu*angstrom^2'), symmetry=1, barrier=(3.73766,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162564,'amu*angstrom^2'), symmetry=1, barrier=(3.73766,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.72606,0.112904,-0.000114995,5.83038e-08,-1.14885e-11,-23781.6,40.508], Tmin=(100,'K'), Tmax=(1246.65,'K')), NASAPolynomial(coeffs=[25.7736,0.0246674,-8.82548e-06,1.52764e-09,-1.02686e-13,-30638.1,-98.2245], Tmin=(1246.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-199.537,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)(O)CJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C(O)(C(=C)[O])C(O)=CC(9527)',
    structure = SMILES('[CH2]C(O)(C(=C)[O])C(O)=CC'),
    E0 = (-254.247,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,180,180,180,180,1600,1703.58,2796.42,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154109,'amu*angstrom^2'), symmetry=1, barrier=(3.54327,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154109,'amu*angstrom^2'), symmetry=1, barrier=(3.54327,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154109,'amu*angstrom^2'), symmetry=1, barrier=(3.54327,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154109,'amu*angstrom^2'), symmetry=1, barrier=(3.54327,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154109,'amu*angstrom^2'), symmetry=1, barrier=(3.54327,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154109,'amu*angstrom^2'), symmetry=1, barrier=(3.54327,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.7888,0.11601,-0.000123571,6.56439e-08,-1.35266e-11,-30360.8,40.9012], Tmin=(100,'K'), Tmax=(1194.97,'K')), NASAPolynomial(coeffs=[25.5939,0.0243486,-8.5114e-06,1.45204e-09,-9.68676e-14,-36905,-96.0817], Tmin=(1194.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-254.247,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)(O)CJ) + radical(C=C(C)OJ)"""),
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
    label = 'C=C([O])[C](O)C(=O)CC(9528)',
    structure = SMILES('C=C([O])C(O)=C([O])CC'),
    E0 = (-331.485,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,325,350,375,415,440,465,420,435,450,1700,1725,1750,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,208.541,208.669,208.681,208.749],'cm^-1')),
        HinderedRotor(inertia=(0.429466,'amu*angstrom^2'), symmetry=1, barrier=(13.2563,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.42854,'amu*angstrom^2'), symmetry=1, barrier=(13.2568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.428109,'amu*angstrom^2'), symmetry=1, barrier=(13.2564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.428697,'amu*angstrom^2'), symmetry=1, barrier=(13.2565,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.52798,0.107321,-0.00012339,6.89722e-08,-1.45623e-11,-39657.1,34.0946], Tmin=(100,'K'), Tmax=(1276.75,'K')), NASAPolynomial(coeffs=[25.681,0.0130805,-2.10173e-06,1.21849e-10,-1.74264e-16,-45871.7,-100.949], Tmin=(1276.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-331.485,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=C(C)OJ)"""),
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
    label = '[CH2]C(O)([C]=C)C(=O)CC(9529)',
    structure = SMILES('[CH2]C(O)([C]=C)C(=O)CC'),
    E0 = (66.5323,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1277.5,1000,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1716.33,2802.01,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154991,'amu*angstrom^2'), symmetry=1, barrier=(3.56355,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154991,'amu*angstrom^2'), symmetry=1, barrier=(3.56355,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154991,'amu*angstrom^2'), symmetry=1, barrier=(3.56355,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154991,'amu*angstrom^2'), symmetry=1, barrier=(3.56355,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154991,'amu*angstrom^2'), symmetry=1, barrier=(3.56355,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154991,'amu*angstrom^2'), symmetry=1, barrier=(3.56355,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.18847,0.124008,-0.000189765,1.59755e-07,-5.29925e-11,8179.34,35.558], Tmin=(100,'K'), Tmax=(847.227,'K')), NASAPolynomial(coeffs=[12.1305,0.0447765,-2.0543e-05,3.82204e-09,-2.58945e-13,6509.25,-23.0275], Tmin=(847.227,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(66.5323,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CC(O)(C=O)CJ)"""),
)

species(
    label = '[CH2]C1(O)C(=O)CC1([O])CC(9319)',
    structure = SMILES('[CH2]C1(O)C(=O)CC1([O])CC'),
    E0 = (-117.206,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.39288,0.106077,-9.2262e-05,3.37367e-08,-2.48605e-12,-13891.4,35.7626], Tmin=(100,'K'), Tmax=(1007.44,'K')), NASAPolynomial(coeffs=[23.9752,0.0281196,-1.00844e-05,1.78647e-09,-1.23664e-13,-20158.1,-92.5454], Tmin=(1007.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-117.206,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(CC(C)2OJ) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = '[CH2]C([O])(C(C)=O)C(=O)CC(9530)',
    structure = SMILES('[CH2]C([O])(C(C)=O)C(=O)CC'),
    E0 = (-178.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.03904,0.119644,-0.000169893,1.40246e-07,-4.69977e-11,-21276.4,40.1854], Tmin=(100,'K'), Tmax=(804.968,'K')), NASAPolynomial(coeffs=[10.9838,0.0504824,-2.34639e-05,4.43938e-09,-3.05737e-13,-22906.9,-13.3142], Tmin=(804.968,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-178.341,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-OdCsCs) + radical(CC(C)(C=O)OJ) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = '[CH2]C(O)(C(C)=O)C(=O)[CH]C(9531)',
    structure = SMILES('[CH2]C(O)(C(C)=O)C(=O)[CH]C'),
    E0 = (-258.195,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.09962,0.112787,-0.000122399,6.98567e-08,-1.60085e-11,-30870.4,38.3789], Tmin=(100,'K'), Tmax=(1056.27,'K')), NASAPolynomial(coeffs=[18.5539,0.0383611,-1.67065e-05,3.14882e-09,-2.1996e-13,-35022.3,-57.5137], Tmin=(1056.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-258.195,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-OdCsCs) + radical(CJC(C)(C=O)O) + radical(CCJC=O)"""),
)

species(
    label = '[CH2]CC(=O)C([CH2])(O)C(C)=O(9532)',
    structure = SMILES('[CH2]CC(=O)C([CH2])(O)C(C)=O'),
    E0 = (-208.734,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,365,385,505,600,445,480,1700,1720,3615,1277.5,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,180,1600,1657.12,2852.63,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153081,'amu*angstrom^2'), symmetry=1, barrier=(3.51963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153081,'amu*angstrom^2'), symmetry=1, barrier=(3.51963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153081,'amu*angstrom^2'), symmetry=1, barrier=(3.51963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153081,'amu*angstrom^2'), symmetry=1, barrier=(3.51963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153081,'amu*angstrom^2'), symmetry=1, barrier=(3.51963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153081,'amu*angstrom^2'), symmetry=1, barrier=(3.51963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153081,'amu*angstrom^2'), symmetry=1, barrier=(3.51963,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.41421,0.125808,-0.00016931,1.2389e-07,-3.64112e-11,-24916.1,39.7622], Tmin=(100,'K'), Tmax=(831.312,'K')), NASAPolynomial(coeffs=[15.8552,0.0427138,-1.93771e-05,3.65295e-09,-2.52407e-13,-27787.3,-40.3622], Tmin=(831.312,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-208.734,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-OdCsCs) + radical(CJC(C)(C=O)O) + radical(CJCC=O)"""),
)

species(
    label = 'C=C1OCC1(O)[C]([O])CC(9340)',
    structure = SMILES('C=C1OCC1(O)[C]([O])CC'),
    E0 = (-89.9066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.40662,0.113693,-0.000113356,5.59993e-08,-1.03995e-11,-10559.3,37.7935], Tmin=(100,'K'), Tmax=(1517.31,'K')), NASAPolynomial(coeffs=[27.4023,0.0184176,-2.66716e-06,1.15336e-10,2.76637e-15,-17683.8,-112.114], Tmin=(1517.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-89.9066,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(C2CsJOH) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C1(O)[C](CC)OCC1=O(9019)',
    structure = SMILES('[CH2]C1(O)[C](CC)OCC1=O'),
    E0 = (-167.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.25054,0.106508,-0.000104653,5.30037e-08,-1.06057e-11,-20003.8,35.9847], Tmin=(100,'K'), Tmax=(1217.1,'K')), NASAPolynomial(coeffs=[21.7218,0.0310103,-1.16077e-05,2.03866e-09,-1.37321e-13,-25595.8,-79.3573], Tmin=(1217.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-167.956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclopentane) + radical(C2CsJOCs) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = '[CH2]C(O)(C[C]=O)C(=O)CC(8854)',
    structure = SMILES('[CH2]C(O)(C[C]=O)C(=O)CC'),
    E0 = (-222.408,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1855,455,950,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,180,180,180,180,1600,1777.09,2739.24,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157068,'amu*angstrom^2'), symmetry=1, barrier=(3.6113,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157068,'amu*angstrom^2'), symmetry=1, barrier=(3.6113,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157068,'amu*angstrom^2'), symmetry=1, barrier=(3.6113,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157068,'amu*angstrom^2'), symmetry=1, barrier=(3.6113,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157068,'amu*angstrom^2'), symmetry=1, barrier=(3.6113,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157068,'amu*angstrom^2'), symmetry=1, barrier=(3.6113,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157068,'amu*angstrom^2'), symmetry=1, barrier=(3.6113,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4709.65,'J/mol'), sigma=(7.59808,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=735.64 K, Pc=24.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.88299,0.13861,-0.000212102,1.72889e-07,-5.52672e-11,-26546.4,40.186], Tmin=(100,'K'), Tmax=(855.013,'K')), NASAPolynomial(coeffs=[15.7712,0.0422671,-1.89583e-05,3.48228e-09,-2.33891e-13,-29062.7,-39.2807], Tmin=(855.013,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-222.408,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(CJC(C)(C=O)O) + radical(CCCJ=O)"""),
)

species(
    label = 'CCC(=O)C1(O)CCC1=O(8855)',
    structure = SMILES('CCC(=O)C1(O)CCC1=O'),
    E0 = (-489.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.497262,0.0859157,-6.20011e-05,2.25653e-08,-3.30352e-12,-58721,36.0152], Tmin=(100,'K'), Tmax=(1605.09,'K')), NASAPolynomial(coeffs=[20.3845,0.0338765,-1.33687e-05,2.36597e-09,-1.5736e-13,-65424.4,-74.608], Tmin=(1605.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-489.67,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-OdCsCs) + ring(Cyclobutanone)"""),
)

species(
    label = '[CH2]C(O)=C([CH2])OC(=O)CC(9533)',
    structure = SMILES('[CH2]C(O)=C([CH2])OC(=O)CC'),
    E0 = (-350.028,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,3615,1277.5,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
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
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.35728,0.115934,-0.00012895,7.46602e-08,-1.71411e-11,-41904,37.2483], Tmin=(100,'K'), Tmax=(1062.5,'K')), NASAPolynomial(coeffs=[20.2645,0.0345342,-1.40329e-05,2.55521e-09,-1.75194e-13,-46498.7,-68.3753], Tmin=(1062.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-350.028,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsOs) + radical(C=C(O)CJ) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]C(O)([C]=O)C(=O)CC(9534)',
    structure = SMILES('[CH2]C(O)([C]=O)C(=O)CC'),
    E0 = (-209.121,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1855,455,950,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.478282,0.108348,-0.000143834,9.19645e-08,-1.74808e-11,-24999.3,34.299], Tmin=(100,'K'), Tmax=(629.153,'K')), NASAPolynomial(coeffs=[13.4902,0.0374454,-1.74798e-05,3.31168e-09,-2.28405e-13,-27111.4,-29.4347], Tmin=(629.153,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-209.121,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(CJC(C)(C=O)O) + radical(CC(C)(O)CJ=O)"""),
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
    E0 = (-248.063,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-161.405,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-98.1059,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-103.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-224.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-182.676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-202.423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-138.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-78.3221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (53.3569,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-39.774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-161.454,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-169.882,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-222.332,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-195.345,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-129.51,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-24.758,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-16.8471,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-61.6381,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-122.299,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (97.395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (148.005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (197.216,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-90.7449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-90.7449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-239.779,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (90.1174,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (126.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-151.471,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (59.452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-111.134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (50.0783,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (309.537,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (-117.206,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (-164.383,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (-201.621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (-158.945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (-89.9066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (-167.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (23.193,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (-239.779,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (-36.2282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (172.442,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    products = ['CH2CO(28)', 'C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    products = ['C=C([O])C1(O)CC1([O])CC(9508)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.34238e+09,'s^-1'), n=0.889391, Ea=(86.6587,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonylbond_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 81.9 to 86.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    products = ['[CH2]C1([O])CC1(O)C(=O)CC(9509)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.24059e+09,'s^-1'), n=0.736667, Ea=(149.957,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 148.4 to 150.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    products = ['[CH2]C1(O)C(=C)OC1([O])CC(9487)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(144.262,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_O] for rate rule [R5_SS_CO;carbonylbond_intra_Nd;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 139.5 to 144.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction5',
    reactants = ['C2H5CO(71)', 'C=C([O])C(=C)O(4628)'],
    products = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00392164,'m^3/(mol*s)'), n=2.41519, Ea=(15.6067,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;CJ] for rate rule [Cds-CdOs_Cds;CO_rad/NonDe]
Euclidian distance = 3.60555127546
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=[C][O](173)', 'C=C(O)C(=O)CC(4626)'],
    products = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.00392164,'m^3/(mol*s)'), n=2.41519, Ea=(15.6067,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;CJ] for rate rule [Cds-COOs_Cds;CJ]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2CO(28)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3525.68,'m^3/(mol*s)'), n=1.07323, Ea=(43.3386,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cdd_Od;CJ] + [Ck_O;YJ] for rate rule [Ck_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['OH(5)', 'C=C([O])C(=C)C(=O)CC(9510)'],
    products = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(931.236,'m^3/(mol*s)'), n=1.015, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;OJ_pri] for rate rule [Cds-CdCO_Cds;OJ_pri]
Euclidian distance = 2.0
family: R_Addition_MultipleBond
Ea raised from -7.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=C([O])C(C)([O])C(=O)CC(9511)'],
    products = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=C(O)C([CH2])(O)C(=O)CC(9512)'],
    products = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C([O])(C(=C)O)C(=O)CC(9513)'],
    products = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1070.11,'s^-1'), n=2.50856, Ea=(101.808,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;O_H_out] + [R4H_SSS;O_rad_out;XH_out] + [R4H_SS(Cd)S;Y_rad_out;XH_out] for rate rule [R4H_SS(Cd)S;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    products = ['C=C([O])C(C)(O)C(=O)[CH]C(9514)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(6.44e+09,'s^-1'), n=0.13, Ea=(86.6088,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 131 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=C([O])C(C)(O)C(=O)CC(9515)'],
    products = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    products = ['[CH2]C(O)(C(=C)O)C(=O)[CH]C(9516)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(8e+10,'s^-1'), n=0, Ea=(25.7316,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(1000,'K'), comment="""Estimated using template [R5H_CCC;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5H_C(Cd)CC;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    products = ['[CH2]CC(=O)C(C)(O)C(=C)[O](9517)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(68850,'s^-1'), n=1.68, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_2H;Cs_H_out_2H] for rate rule [R5H_CCC(O2d)C;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]CC(=O)C([CH2])(O)C(=C)O(9518)'],
    products = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(104,'s^-1'), n=2.1, Ea=(44.7688,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_SSSS(Cd)S;C_rad_out_2H;XH_out] for rate rule [R6H_SSSS(Cd)S;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=[C][O](173)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    products = ['[CH2]C(O)([C]1CO1)C(=O)CC(9519)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    products = ['C=C([O])C1(O)CO[C]1CC(9520)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.01137e+09,'s^-1'), n=0.572544, Ea=(186.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4_S;multiplebond_intra;radadd_intra_cs2H] + [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    products = ['CC[C]([O])C1(O)CCC1=O(8973)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.51071e+08,'s^-1'), n=0.996667, Ea=(125.764,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    products = ['[CH2]C1(O)[C](CC)OOC1=C(9470)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.81184e+09,'s^-1'), n=0.551229, Ea=(345.458,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra] for rate rule [R5_SS_CO;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic
Ea raised from 338.9 to 345.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction22',
    reactants = ['CH3CH2OH(54)', '[CH2]C(=C=O)C(=C)[O](9521)'],
    products = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.79e-11,'m^3/(mol*s)'), n=3.97, Ea=(329.281,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [doublebond;Cs_OH] for rate rule [Cdd_Cd;C_pri_OH]
Euclidian distance = 1.41421356237
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CH2(S)(23)', '[CH2]C(O)(C(=C)[O])C(C)=O(9522)'],
    products = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;C_pri/De]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    products = ['C=C([O])[C](O)CC(=O)CC(9523)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ-HH;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;CO]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    products = ['C=C([O])C[C](O)C(=O)CC(8845)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    products = ['C=C1OCC1(O)C(=O)CC(9371)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][C](C(=C)OO)C(=O)CC(9524)'],
    products = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.01978e+11,'s^-1'), n=0, Ea=(107.111,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OOH_S;Y_rad_out]
Euclidian distance = 0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C(C(=C)[O])=C(CC)OO(9525)'],
    products = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O] for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O_R]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=C([O])O[C](CC)C(=C)O(8844)'],
    products = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O] for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O_C]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(O)(C(=C)[O])C(=C)OC(9526)'],
    products = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(37989.5,'s^-1'), n=2.515, Ea=(258.99,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C] + [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O] for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_C]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(O)(C(=C)[O])C(O)=CC(9527)'],
    products = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O_H] for rate rule [R_ROR;R1_doublebond_CHCH3;R2_doublebond_CsC;R_O_H]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction32',
    reactants = ['CH2(19)', 'C=C([O])[C](O)C(=O)CC(9528)'],
    products = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/TDMustO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction33',
    reactants = ['O(4)', '[CH2]C(O)([C]=C)C(=O)CC(9529)'],
    products = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    products = ['[CH2]C1(O)C(=O)CC1([O])CC(9319)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.46159e+06,'s^-1'), n=1.55572, Ea=(130.857,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5_SS_CO;carbonylbond_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 127.3 to 130.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    products = ['[CH2]C([O])(C(C)=O)C(=O)CC(9530)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    products = ['[CH2]C(O)(C(C)=O)C(=O)[CH]C(9531)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(262000,'s^-1'), n=1.62, Ea=(46.4424,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_2H;Cs_H_out_H/NonDeC] for rate rule [R5H_CC(O2d)CC;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]CC(=O)C([CH2])(O)C(C)=O(9532)'],
    products = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(3690,'s^-1'), n=1.79, Ea=(49.7896,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 113 used for R6H_SSSSS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R6H_SSSSS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    products = ['C=C1OCC1(O)[C]([O])CC(9340)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(3.27074e+08,'s^-1'), n=0.924088, Ea=(158.157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_cs2H]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic
Ea raised from 154.1 to 158.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    products = ['[CH2]C1(O)[C](CC)OCC1=O(9019)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(2.64784e+07,'s^-1'), n=0.990488, Ea=(80.1071,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5_SS_CO;carbonyl_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic
Ea raised from 76.3 to 80.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(O)(C[C]=O)C(=O)CC(8854)'],
    products = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    products = ['CCC(=O)C1(O)CCC1=O(8855)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2]C(O)=C([CH2])OC(=O)CC(9533)'],
    products = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction43',
    reactants = ['CH2(19)', '[CH2]C(O)([C]=O)C(=O)CC(9534)'],
    products = ['[CH2]C(O)(C(=C)[O])C(=O)CC(8846)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '1792',
    isomers = [
        '[CH2]C(O)(C(=C)[O])C(=O)CC(8846)',
    ],
    reactants = [
        ('CH2CO(28)', 'C=C(O)C(=O)CC(4626)'),
        ('CH2CO(28)', '[CH2]C(O)=C([O])CC(4557)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '1792',
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

