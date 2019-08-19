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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.39668,0.149333,-0.000195449,1.45429e-07,-4.41529e-11,-11190.4,42.2789], Tmin=(100,'K'), Tmax=(801.412,'K')), NASAPolynomial(coeffs=[15.6661,0.0591689,-2.66711e-05,5.01446e-09,-3.45918e-13,-14085.2,-40.863], Tmin=(801.412,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-94.891,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(C=CC(O)(C=O)CJ)"""),
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
    label = '[CH2]C1([CH]C)CC1(O)C(=O)CC(25638)',
    structure = SMILES('[CH2]C1([CH]C)CC1(O)C(=O)CC'),
    E0 = (-17.7057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.65243,0.119336,-0.000103757,4.72642e-08,-8.76721e-12,-1922.03,42.2294], Tmin=(100,'K'), Tmax=(1277.4,'K')), NASAPolynomial(coeffs=[20.9883,0.048441,-2.05085e-05,3.81756e-09,-2.64343e-13,-7706.34,-72.5422], Tmin=(1277.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-17.7057,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsCs) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclopropane) + radical(Cs_S) + radical(Neopentyl)"""),
)

species(
    label = '[CH2]C(=CC)C1(O)CC1([O])CC(25639)',
    structure = SMILES('[CH2]C(=CC)C1(O)CC1([O])CC'),
    E0 = (-8.23227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.03024,0.115477,-7.5537e-05,8.05486e-09,7.32279e-12,-757.858,41.8961], Tmin=(100,'K'), Tmax=(997.449,'K')), NASAPolynomial(coeffs=[25.694,0.0395355,-1.43288e-05,2.56519e-09,-1.7896e-13,-8041.58,-100.574], Tmin=(997.449,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.23227,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclopropane) + radical(Allyl_P) + radical(CC(C)2OJ)"""),
)

species(
    label = '[CH2]C1(O)C(=CC)CC1([O])CC(25640)',
    structure = SMILES('[CH2]C1(O)C(=CC)CC1([O])CC'),
    E0 = (47.0461,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[8.85565,0.0465354,6.31884e-05,-8.80389e-08,2.10357e-11,5348.73,-4.61778], Tmin=(100,'K'), Tmax=(1804.37,'K')), NASAPolynomial(coeffs=[102.358,0.0217628,-6.79402e-05,1.64672e-08,-1.21087e-12,-58103.9,-593.229], Tmin=(1804.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(47.0461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(methylenecyclobutane) + radical(C=CC(C)(O)CJ) + radical(CC(C)2OJ)"""),
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
    label = '[CH2]C(=CC)C(=C)O(25079)',
    structure = SMILES('[CH2]C(=CC)C(=C)O'),
    E0 = (-42.0887,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.902368,'amu*angstrom^2'), symmetry=1, barrier=(20.7472,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.902031,'amu*angstrom^2'), symmetry=1, barrier=(20.7395,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.902103,'amu*angstrom^2'), symmetry=1, barrier=(20.7411,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.902314,'amu*angstrom^2'), symmetry=1, barrier=(20.746,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.121751,0.073074,-4.45307e-05,-4.84265e-09,1.02418e-11,-4911.33,23.8451], Tmin=(100,'K'), Tmax=(945.744,'K')), NASAPolynomial(coeffs=[19.9712,0.0188285,-5.6114e-06,9.36059e-10,-6.54532e-14,-9994.37,-77.8334], Tmin=(945.744,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-42.0887,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P)"""),
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
    label = '[CH2]C(=CC)C(=C)C(=O)CC(25641)',
    structure = SMILES('[CH2]C(=CC)C(=C)C(=O)CC'),
    E0 = (-12.7211,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.983317,0.10755,-9.06778e-05,4.04319e-08,-7.43817e-12,-1349.11,34.1106], Tmin=(100,'K'), Tmax=(1273.47,'K')), NASAPolynomial(coeffs=[17.8597,0.0483632,-2.09633e-05,3.93612e-09,-2.73535e-13,-6148.32,-61.3513], Tmin=(1273.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-12.7211,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=CC)C(C)([O])C(=O)CC(25642)',
    structure = SMILES('[CH2]C(=CC)C(C)([O])C(=O)CC'),
    E0 = (-63.8278,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,545.241,613.389,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.160027,'amu*angstrom^2'), symmetry=1, barrier=(3.67934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160027,'amu*angstrom^2'), symmetry=1, barrier=(3.67934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160027,'amu*angstrom^2'), symmetry=1, barrier=(3.67934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160027,'amu*angstrom^2'), symmetry=1, barrier=(3.67934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160027,'amu*angstrom^2'), symmetry=1, barrier=(3.67934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160027,'amu*angstrom^2'), symmetry=1, barrier=(3.67934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160027,'amu*angstrom^2'), symmetry=1, barrier=(3.67934,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.02572,0.143511,-0.00019767,1.64247e-07,-5.5921e-11,-7470.26,41.8719], Tmin=(100,'K'), Tmax=(805.198,'K')), NASAPolynomial(coeffs=[10.6481,0.0673082,-3.10432e-05,5.8622e-09,-4.03636e-13,-9081.96,-13.8602], Tmin=(805.198,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-63.8278,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + radical(C=CC(C)(C=O)OJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(O)(C(=O)CC)C(C)=[C]C(25643)',
    structure = SMILES('[CH2]C(O)(C(=O)CC)C(C)=[C]C'),
    E0 = (-8.54836,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3100,440,815,1455,1000,1685,370,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1764.52,2755.57,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156646,'amu*angstrom^2'), symmetry=1, barrier=(3.60161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156646,'amu*angstrom^2'), symmetry=1, barrier=(3.60161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156646,'amu*angstrom^2'), symmetry=1, barrier=(3.60161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156646,'amu*angstrom^2'), symmetry=1, barrier=(3.60161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156646,'amu*angstrom^2'), symmetry=1, barrier=(3.60161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156646,'amu*angstrom^2'), symmetry=1, barrier=(3.60161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156646,'amu*angstrom^2'), symmetry=1, barrier=(3.60161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156646,'amu*angstrom^2'), symmetry=1, barrier=(3.60161,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.6219,0.158934,-0.000238771,2.02955e-07,-6.84528e-11,-802.435,43.2655], Tmin=(100,'K'), Tmax=(839.825,'K')), NASAPolynomial(coeffs=[12.7434,0.063721,-2.93656e-05,5.4911e-09,-3.73761e-13,-2606.38,-23.5557], Tmin=(839.825,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.54836,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + radical(C=CC(O)(C=O)CJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([O])(C(=O)CC)C(C)=CC(25644)',
    structure = SMILES('[CH2]C([O])(C(=O)CC)C(C)=CC'),
    E0 = (-2.10361,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,408.797,749.007,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.16284,'amu*angstrom^2'), symmetry=1, barrier=(3.74402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16284,'amu*angstrom^2'), symmetry=1, barrier=(3.74402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16284,'amu*angstrom^2'), symmetry=1, barrier=(3.74402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16284,'amu*angstrom^2'), symmetry=1, barrier=(3.74402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16284,'amu*angstrom^2'), symmetry=1, barrier=(3.74402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16284,'amu*angstrom^2'), symmetry=1, barrier=(3.74402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16284,'amu*angstrom^2'), symmetry=1, barrier=(3.74402,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.38181,0.157271,-0.0002461,2.20829e-07,-7.76816e-11,-39.4743,42.9191], Tmin=(100,'K'), Tmax=(846.576,'K')), NASAPolynomial(coeffs=[8.7428,0.0709772,-3.34349e-05,6.29299e-09,-4.28992e-13,-714.3,-1.75883], Tmin=(846.576,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-2.10361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + radical(C=CC(C)(C=O)OJ) + radical(C=CC(O)(C=O)CJ)"""),
)

species(
    label = '[CH2]C(=[C]C)C(C)(O)C(=O)CC(25645)',
    structure = SMILES('[CH2]C(=[C]C)C(C)(O)C(=O)CC'),
    E0 = (-70.2726,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3100,440,815,1455,1000,1685,370,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1661.29,2852.14,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153633,'amu*angstrom^2'), symmetry=1, barrier=(3.53233,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153633,'amu*angstrom^2'), symmetry=1, barrier=(3.53233,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153633,'amu*angstrom^2'), symmetry=1, barrier=(3.53233,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153633,'amu*angstrom^2'), symmetry=1, barrier=(3.53233,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153633,'amu*angstrom^2'), symmetry=1, barrier=(3.53233,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153633,'amu*angstrom^2'), symmetry=1, barrier=(3.53233,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153633,'amu*angstrom^2'), symmetry=1, barrier=(3.53233,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153633,'amu*angstrom^2'), symmetry=1, barrier=(3.53233,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.13552,0.143455,-0.000183382,1.35774e-07,-4.13522e-11,-8238.74,41.7621], Tmin=(100,'K'), Tmax=(796.671,'K')), NASAPolynomial(coeffs=[14.3761,0.0605478,-2.72737e-05,5.1336e-09,-3.54642e-13,-10869.5,-34.1431], Tmin=(796.671,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-70.2726,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=CC)C(C)(O)C(=O)[CH]C(25646)',
    structure = SMILES('[CH2]C(=CC)C(C)(O)C([O])=CC'),
    E0 = (-176.714,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.87037,0.120999,-0.000105849,4.84601e-08,-8.96665e-12,-21035.7,43.1064], Tmin=(100,'K'), Tmax=(1289.64,'K')), NASAPolynomial(coeffs=[22.2875,0.0460703,-1.86987e-05,3.40886e-09,-2.33411e-13,-27266.8,-79.5865], Tmin=(1289.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-176.714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(Allyl_P)"""),
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
    label = '[CH2]C(O)(C(=O)[CH]C)C(C)=CC(25647)',
    structure = SMILES('[CH2]C(O)(C(C)=CC)C([O])=CC'),
    E0 = (-114.769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.5149,0.123746,-0.000119783,6.38247e-08,-1.40816e-11,-13607,42.4334], Tmin=(100,'K'), Tmax=(1076.97,'K')), NASAPolynomial(coeffs=[17.004,0.0549648,-2.39846e-05,4.52315e-09,-3.15706e-13,-17595.8,-48.2827], Tmin=(1076.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-114.769,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(C=CC(C)(O)CJ)"""),
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
    label = '[CH2]CC(=O)C(C)(O)C([CH2])=CC(25648)',
    structure = SMILES('[CH2]CC(=O)C(C)(O)C([CH2])=CC'),
    E0 = (-96.5253,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.17203,0.142992,-0.000174258,1.19655e-07,-3.35951e-11,-11393.5,42.058], Tmin=(100,'K'), Tmax=(863.219,'K')), NASAPolynomial(coeffs=[16.2966,0.057412,-2.55466e-05,4.80492e-09,-3.32867e-13,-14582,-44.3259], Tmin=(863.219,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-96.5253,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + radical(CJCC=O) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]CC(=O)C([CH2])(O)C(C)=CC(25649)',
    structure = SMILES('[CH2]CC(=O)C([CH2])(O)C(C)=CC'),
    E0 = (-34.8011,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3615,1277.5,1000,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1729.74,2788.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155639,'amu*angstrom^2'), symmetry=1, barrier=(3.57844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155639,'amu*angstrom^2'), symmetry=1, barrier=(3.57844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155639,'amu*angstrom^2'), symmetry=1, barrier=(3.57844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155639,'amu*angstrom^2'), symmetry=1, barrier=(3.57844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155639,'amu*angstrom^2'), symmetry=1, barrier=(3.57844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155639,'amu*angstrom^2'), symmetry=1, barrier=(3.57844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155639,'amu*angstrom^2'), symmetry=1, barrier=(3.57844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155639,'amu*angstrom^2'), symmetry=1, barrier=(3.57844,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.72811,0.159344,-0.000232971,1.91633e-07,-6.30129e-11,-3954.2,43.8083], Tmin=(100,'K'), Tmax=(829.242,'K')), NASAPolynomial(coeffs=[14.7193,0.0604919,-2.75854e-05,5.15004e-09,-3.50972e-13,-6342.72,-34.0531], Tmin=(829.242,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-34.8011,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + radical(C=CC(O)(C=O)CJ) + radical(CJCC=O)"""),
)

species(
    label = '[CH2]C(O)([C]1CC1C)C(=O)CC(25650)',
    structure = SMILES('[CH2]C(O)([C]1CC1C)C(=O)CC'),
    E0 = (-50.0422,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.40588,0.126357,-0.000114046,5.30792e-08,-9.80559e-12,-5775.89,42.8953], Tmin=(100,'K'), Tmax=(1310.56,'K')), NASAPolynomial(coeffs=[26.0099,0.0396276,-1.47798e-05,2.58338e-09,-1.73047e-13,-13224,-101.879], Tmin=(1310.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-50.0422,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclopropane) + radical(CJC(C)(C=O)O) + radical(CCJ(C)CO)"""),
)

species(
    label = '[CH2][C]1C(C)CC1(O)C(=O)CC(25651)',
    structure = SMILES('[CH2][C]1C(C)CC1(O)C(=O)CC'),
    E0 = (-60.9906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.51233,0.106389,-6.53439e-05,9.08775e-09,4.28696e-12,-7124.2,41.1676], Tmin=(100,'K'), Tmax=(1040,'K')), NASAPolynomial(coeffs=[20.9721,0.0465962,-1.75936e-05,3.15089e-09,-2.16698e-13,-13244.1,-75.1268], Tmin=(1040,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-60.9906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutane) + radical(CCJ(C)CO) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(=CC)C1(O)CO[C]1CC(25652)',
    structure = SMILES('[CH2]C(=CC)C1(O)CO[C]1CC'),
    E0 = (-17.0253,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.09223,0.12175,-0.000113014,5.83515e-08,-1.20802e-11,-1817.95,42.7812], Tmin=(100,'K'), Tmax=(1199.08,'K')), NASAPolynomial(coeffs=[20.7278,0.0440547,-1.38556e-05,2.1288e-09,-1.30372e-13,-7177.66,-70.9845], Tmin=(1199.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-17.0253,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Oxetane) + radical(Allyl_P) + radical(C2CsJOCs)"""),
)

species(
    label = '[CH2]C1(O)[C](CC)OCC1=CC(25653)',
    structure = SMILES('[CH2]C1(O)[C](CC)OCC1=CC'),
    E0 = (-25.2522,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.35216,0.108725,-8.40464e-05,3.47122e-08,-5.89292e-12,-2837.2,42.809], Tmin=(100,'K'), Tmax=(1381.84,'K')), NASAPolynomial(coeffs=[19.3,0.0489434,-1.91531e-05,3.40456e-09,-2.28816e-13,-8544.81,-63.5051], Tmin=(1381.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-25.2522,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclopentane) + radical(C=CC(C)(O)CJ) + radical(C2CsJOCs)"""),
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
    label = '[CH2]C(=C=O)C([CH2])=CC(25654)',
    structure = SMILES('[CH2]C(=CC)C(=C)[C]=O'),
    E0 = (221.454,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.945205,'amu*angstrom^2'), symmetry=1, barrier=(21.7321,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.944509,'amu*angstrom^2'), symmetry=1, barrier=(21.7161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.942945,'amu*angstrom^2'), symmetry=1, barrier=(21.6802,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.943418,'amu*angstrom^2'), symmetry=1, barrier=(21.691,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.242108,0.0834175,-7.66663e-05,3.48309e-08,-6.22063e-12,26795.8,27.7167], Tmin=(100,'K'), Tmax=(1356.62,'K')), NASAPolynomial(coeffs=[20.423,0.0224858,-9.29424e-06,1.72283e-09,-1.19375e-13,21188.9,-78.2831], Tmin=(1356.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(221.454,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_P) + radical(C=C(C)CJ=O)"""),
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
    label = '[CH2]C(=CC)C([CH2])(O)C(C)=O(25655)',
    structure = SMILES('[CH2]C(=CC)C([CH2])(O)C(C)=O'),
    E0 = (-69.4741,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,375,552.5,462.5,1710,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.45019,0.125124,-0.000148875,9.78991e-08,-2.61713e-11,-8164.08,37.5535], Tmin=(100,'K'), Tmax=(906.851,'K')), NASAPolynomial(coeffs=[16.0655,0.0478667,-2.10879e-05,3.95872e-09,-2.74443e-13,-11341,-45.2372], Tmin=(906.851,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-69.4741,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-O2d)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + radical(C=CC(O)(C=O)CJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=C)C([CH2])(O)C(=O)CC(25656)',
    structure = SMILES('[CH2]C(=C)C([CH2])(O)C(=O)CC'),
    E0 = (-58.8654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,186.472,1569.07,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.150282,'amu*angstrom^2'), symmetry=1, barrier=(3.45529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150282,'amu*angstrom^2'), symmetry=1, barrier=(3.45529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150282,'amu*angstrom^2'), symmetry=1, barrier=(3.45529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150282,'amu*angstrom^2'), symmetry=1, barrier=(3.45529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150282,'amu*angstrom^2'), symmetry=1, barrier=(3.45529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150282,'amu*angstrom^2'), symmetry=1, barrier=(3.45529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150282,'amu*angstrom^2'), symmetry=1, barrier=(3.45529,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.77289,0.134027,-0.000175787,1.27994e-07,-3.76922e-11,-6878.42,37.5401], Tmin=(100,'K'), Tmax=(827.924,'K')), NASAPolynomial(coeffs=[15.8026,0.0491121,-2.19383e-05,4.10862e-09,-2.82929e-13,-9788.6,-43.9323], Tmin=(827.924,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-58.8654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=CC(O)(C=O)CJ)"""),
)

species(
    label = 'C=C([CH]C)C[C](O)C(=O)CC(25634)',
    structure = SMILES('[CH2]C(=CC)CC(O)=C([O])CC'),
    E0 = (-184.815,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.29682,0.123414,-9.9085e-05,3.1192e-08,-3.96361e-13,-21988.2,44.6011], Tmin=(100,'K'), Tmax=(1005.85,'K')), NASAPolynomial(coeffs=[26.2908,0.0379147,-1.36164e-05,2.40395e-09,-1.65592e-13,-29165,-100.572], Tmin=(1005.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-184.815,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=CC)[C](O)CC(=O)CC(25657)',
    structure = SMILES('[CH2]C([CH]C)=C(O)CC(=O)CC'),
    E0 = (-193.099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.04431,0.119692,-9.61201e-05,3.83082e-08,-6.10108e-12,-22996.2,41.0074], Tmin=(100,'K'), Tmax=(1487.86,'K')), NASAPolynomial(coeffs=[26.7125,0.0423815,-1.8179e-05,3.38515e-09,-2.33089e-13,-31553.4,-109.154], Tmin=(1487.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-193.099,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-OdCsCs) + radical(Allyl_S) + radical(Allyl_P)"""),
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
    label = 'CC=C1CCC1(O)C(=O)CC(25658)',
    structure = SMILES('CC=C1CCC1(O)C(=O)CC'),
    E0 = (-303.683,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[10.0022,0.0379217,7.41309e-05,-9.30223e-08,2.16944e-11,-36890.3,-11.0922], Tmin=(100,'K'), Tmax=(1821.76,'K')), NASAPolynomial(coeffs=[103.764,0.0220986,-6.93215e-05,1.6737e-08,-1.22607e-12,-102589,-606.233], Tmin=(1821.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-303.683,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + ring(methylenecyclobutane)"""),
)

species(
    label = '[CH2]C(=CC)C([CH2])=C(CC)OO(25659)',
    structure = SMILES('[CH2]C(=CC)C([CH2])=C(CC)OO'),
    E0 = (136.978,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1310,387.5,850,1000,325,350,375,415,440,465,420,435,450,1700,1725,1750,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.12621,0.127625,-0.000116606,5.55103e-08,-1.06444e-11,16701.1,42.3733], Tmin=(100,'K'), Tmax=(1249.37,'K')), NASAPolynomial(coeffs=[23.2151,0.0464913,-1.9195e-05,3.53108e-09,-2.43248e-13,10369,-85.5255], Tmin=(1249.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(136.978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = 'C=C([CH]C)O[C](CC)C(=C)O(25633)',
    structure = SMILES('[CH2]C(=CC)OC(CC)=C([CH2])O'),
    E0 = (-105.191,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,325,350,375,415,440,465,420,435,450,1700,1725,1750,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.2723,0.130167,-0.0001281,6.62641e-08,-1.36647e-11,-12419.1,45.5225], Tmin=(100,'K'), Tmax=(1176.41,'K')), NASAPolynomial(coeffs=[23.5268,0.0424455,-1.62481e-05,2.87814e-09,-1.9446e-13,-18489.2,-83.1347], Tmin=(1176.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-105.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(O)CJ) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]C(=CC)C([CH2])(O)C(=C)OC(25660)',
    structure = SMILES('[CH2]C(=CC)C([CH2])(O)C(=C)OC'),
    E0 = (-46.3652,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,325,375,415,465,420,450,1700,1750,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,180,282.337,845.691,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15284,'amu*angstrom^2'), symmetry=1, barrier=(3.5141,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15284,'amu*angstrom^2'), symmetry=1, barrier=(3.5141,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15284,'amu*angstrom^2'), symmetry=1, barrier=(3.5141,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15284,'amu*angstrom^2'), symmetry=1, barrier=(3.5141,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15284,'amu*angstrom^2'), symmetry=1, barrier=(3.5141,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15284,'amu*angstrom^2'), symmetry=1, barrier=(3.5141,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15284,'amu*angstrom^2'), symmetry=1, barrier=(3.5141,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15284,'amu*angstrom^2'), symmetry=1, barrier=(3.5141,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.34285,0.12838,-0.000117242,5.47393e-08,-1.01792e-11,-5338.85,44.0082], Tmin=(100,'K'), Tmax=(1295.89,'K')), NASAPolynomial(coeffs=[25.7733,0.0415959,-1.67909e-05,3.06372e-09,-2.10289e-13,-12626.1,-98.9244], Tmin=(1295.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-46.3652,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]C(=CC)C([CH2])(O)C(O)=CC(25661)',
    structure = SMILES('[CH2]C(=CC)C([CH2])(O)C(O)=CC'),
    E0 = (-101.075,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2995,3025,975,1000,1300,1375,400,500,1630,1680,3580,3650,1210,1345,900,1100,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.39742,0.131385,-0.000125465,6.16457e-08,-1.20503e-11,-11918.4,44.3726], Tmin=(100,'K'), Tmax=(1237.34,'K')), NASAPolynomial(coeffs=[25.3593,0.0416547,-1.66864e-05,3.03627e-09,-2.08381e-13,-18787.2,-95.4484], Tmin=(1237.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-101.075,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C=CC(C)(O)CJ) + radical(Allyl_P)"""),
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
    label = '[CH2]C(=CC)[C](O)C(=O)CC(25662)',
    structure = SMILES('[CH2]C(=CC)C(O)=C([O])CC'),
    E0 = (-177.514,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,325,350,375,415,440,465,420,435,450,1700,1725,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.43497,0.122563,-0.000122598,6.16003e-08,-1.19501e-11,-21102.1,39.3045], Tmin=(100,'K'), Tmax=(1332.6,'K')), NASAPolynomial(coeffs=[28.3154,0.0252184,-7.34859e-06,1.10407e-09,-6.8045e-14,-28850,-116.197], Tmin=(1332.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-177.514,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C(O)([C]=CC)C(=O)CC(25663)',
    structure = SMILES('[CH2]C(O)([C]=CC)C(=O)CC'),
    E0 = (30.5067,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1753.41,2766.31,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156364,'amu*angstrom^2'), symmetry=1, barrier=(3.59512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156364,'amu*angstrom^2'), symmetry=1, barrier=(3.59512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156364,'amu*angstrom^2'), symmetry=1, barrier=(3.59512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156364,'amu*angstrom^2'), symmetry=1, barrier=(3.59512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156364,'amu*angstrom^2'), symmetry=1, barrier=(3.59512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156364,'amu*angstrom^2'), symmetry=1, barrier=(3.59512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156364,'amu*angstrom^2'), symmetry=1, barrier=(3.59512,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.81583,0.139369,-0.000209685,1.77619e-07,-5.96777e-11,3867.54,39.6156], Tmin=(100,'K'), Tmax=(839.639,'K')), NASAPolynomial(coeffs=[12.0277,0.0547706,-2.52372e-05,4.71837e-09,-3.21117e-13,2200.16,-20.8375], Tmin=(839.639,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(30.5067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(C=CC(O)(C=O)CJ)"""),
)

species(
    label = '[CH2]C1(O)C(=C)C(C)C1([O])CC(25664)',
    structure = SMILES('[CH2]C1(O)C(=C)C(C)C1([O])CC'),
    E0 = (51.3196,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[8.30878,0.0492268,6.07581e-05,-8.75673e-08,2.1155e-11,5890.71,-3.82484], Tmin=(100,'K'), Tmax=(1789.8,'K')), NASAPolynomial(coeffs=[100.402,0.0231125,-6.7963e-05,1.64774e-08,-1.21394e-12,-55857.9,-582.138], Tmin=(1789.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(51.3196,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(C=CC(C)(O)CJ) + radical(CC(C)2OJ)"""),
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
    label = '[CH2]CC(=C)C([CH2])(O)C(=O)CC(25665)',
    structure = SMILES('[CH2]CC(=C)C([CH2])(O)C(=O)CC'),
    E0 = (-27.7637,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,180,180,180,533.45,619.076,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157954,'amu*angstrom^2'), symmetry=1, barrier=(3.63168,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157954,'amu*angstrom^2'), symmetry=1, barrier=(3.63168,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157954,'amu*angstrom^2'), symmetry=1, barrier=(3.63168,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157954,'amu*angstrom^2'), symmetry=1, barrier=(3.63168,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157954,'amu*angstrom^2'), symmetry=1, barrier=(3.63168,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157954,'amu*angstrom^2'), symmetry=1, barrier=(3.63168,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157954,'amu*angstrom^2'), symmetry=1, barrier=(3.63168,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157954,'amu*angstrom^2'), symmetry=1, barrier=(3.63168,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.38538,0.149172,-0.000197596,1.48988e-07,-4.58029e-11,-3117.36,44.8486], Tmin=(100,'K'), Tmax=(791.94,'K')), NASAPolynomial(coeffs=[15.4949,0.0588598,-2.65345e-05,4.98449e-09,-3.43431e-13,-5949.36,-37.2425], Tmin=(791.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-27.7637,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsHH) + radical(RCCJ) + radical(C=CC(O)(C=O)CJ)"""),
)

species(
    label = '[CH]=C(CC)C([CH2])(O)C(=O)CC(25666)',
    structure = SMILES('[CH]=C(CC)C([CH2])(O)C(=O)CC'),
    E0 = (14.0861,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,180,180,180,180,1600,1704.86,2809.64,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154879,'amu*angstrom^2'), symmetry=1, barrier=(3.56097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154879,'amu*angstrom^2'), symmetry=1, barrier=(3.56097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154879,'amu*angstrom^2'), symmetry=1, barrier=(3.56097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154879,'amu*angstrom^2'), symmetry=1, barrier=(3.56097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154879,'amu*angstrom^2'), symmetry=1, barrier=(3.56097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154879,'amu*angstrom^2'), symmetry=1, barrier=(3.56097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154879,'amu*angstrom^2'), symmetry=1, barrier=(3.56097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154879,'amu*angstrom^2'), symmetry=1, barrier=(3.56097,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.26154,0.14861,-0.000188705,1.28022e-07,-3.21711e-11,1909.69,43.0545], Tmin=(100,'K'), Tmax=(675.814,'K')), NASAPolynomial(coeffs=[15.8005,0.0586071,-2.6454e-05,4.97573e-09,-3.43338e-13,-917.617,-39.863], Tmin=(675.814,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(14.0861,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsHH) + radical(C=CC(O)(C=O)CJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C([O])(C(=C)CC)C(=O)CC(25667)',
    structure = SMILES('[CH2]C([O])(C(=C)CC)C(=O)CC'),
    E0 = (11.2765,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,180,180,180,180,1036.38,1196.16,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.165424,'amu*angstrom^2'), symmetry=1, barrier=(3.80343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165424,'amu*angstrom^2'), symmetry=1, barrier=(3.80343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165424,'amu*angstrom^2'), symmetry=1, barrier=(3.80343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165424,'amu*angstrom^2'), symmetry=1, barrier=(3.80343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165424,'amu*angstrom^2'), symmetry=1, barrier=(3.80343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165424,'amu*angstrom^2'), symmetry=1, barrier=(3.80343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165424,'amu*angstrom^2'), symmetry=1, barrier=(3.80343,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.35497,0.153194,-0.000226505,1.94684e-07,-6.6873e-11,1572.18,43.9309], Tmin=(100,'K'), Tmax=(831.187,'K')), NASAPolynomial(coeffs=[10.8843,0.0672392,-3.12472e-05,5.88008e-09,-4.023e-13,139.643,-12.8711], Tmin=(831.187,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(11.2765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsHH) + radical(C=CC(O)(C=O)CJ) + radical(C=CC(C)(C=O)OJ)"""),
)

species(
    label = '[CH]=C([CH]C)C(C)(O)C(=O)CC(25668)',
    structure = SMILES('[CH]C(=CC)C(C)(O)C(=O)CC'),
    E0 = (-88.9291,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.99027,0.139851,-0.000161281,1.11391e-07,-3.2382e-11,-10487.3,41.602], Tmin=(100,'K'), Tmax=(825.993,'K')), NASAPolynomial(coeffs=[12.898,0.0677515,-3.03457e-05,5.71046e-09,-3.955e-13,-12946.8,-27.3789], Tmin=(825.993,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-88.9291,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C(O)(C(=C)CC)C(=O)[CH]C(25669)',
    structure = SMILES('[CH2]C(O)(C(=C)CC)C([O])=CC'),
    E0 = (-101.389,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.78302,0.123277,-0.000113317,5.52613e-08,-1.0953e-11,-11982.9,44.4934], Tmin=(100,'K'), Tmax=(1204.22,'K')), NASAPolynomial(coeffs=[20.7063,0.0485764,-2.02693e-05,3.74948e-09,-2.59094e-13,-17399.3,-68.1839], Tmin=(1204.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-101.389,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)(O)CJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]CC(=O)C([CH2])(O)C(=C)CC(25670)',
    structure = SMILES('[CH2]CC(=O)C([CH2])(O)C(=C)CC'),
    E0 = (-21.4209,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,180,180,180,545.09,606.555,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.158152,'amu*angstrom^2'), symmetry=1, barrier=(3.63623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158152,'amu*angstrom^2'), symmetry=1, barrier=(3.63623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158152,'amu*angstrom^2'), symmetry=1, barrier=(3.63623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158152,'amu*angstrom^2'), symmetry=1, barrier=(3.63623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158152,'amu*angstrom^2'), symmetry=1, barrier=(3.63623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158152,'amu*angstrom^2'), symmetry=1, barrier=(3.63623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158152,'amu*angstrom^2'), symmetry=1, barrier=(3.63623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158152,'amu*angstrom^2'), symmetry=1, barrier=(3.63623,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.60979,0.154067,-0.000208553,1.58217e-07,-4.85836e-11,-2346.42,44.4995], Tmin=(100,'K'), Tmax=(795.095,'K')), NASAPolynomial(coeffs=[16.653,0.0571306,-2.5625e-05,4.79261e-09,-3.28995e-13,-5408.7,-44.0102], Tmin=(795.095,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-21.4209,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsHH) + radical(CJCC=O) + radical(C=CC(O)(C=O)CJ)"""),
)

species(
    label = 'C[CH][C]1CCC1(O)C(=O)CC(25671)',
    structure = SMILES('C[CH][C]1CCC1(O)C(=O)CC'),
    E0 = (-62.4893,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.49089,0.107287,-7.70551e-05,2.81335e-08,-4.16132e-12,-7307.07,42.8611], Tmin=(100,'K'), Tmax=(1577.99,'K')), NASAPolynomial(coeffs=[22.9709,0.0452801,-1.81132e-05,3.23207e-09,-2.16234e-13,-15027.2,-86.3115], Tmin=(1577.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-62.4893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutane) + radical(CCJ(C)CO) + radical(Cs_S)"""),
)

species(
    label = '[CH2]C1(O)[C](CC)OC(C)C1=C(25672)',
    structure = SMILES('[CH2]C1(O)[C](CC)OC(C)C1=C'),
    E0 = (-31.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.53714,0.105519,-6.01953e-05,1.07333e-09,7.68127e-12,-3580.23,42.702], Tmin=(100,'K'), Tmax=(1022.88,'K')), NASAPolynomial(coeffs=[22.3166,0.0440019,-1.65637e-05,2.99508e-09,-2.08371e-13,-10121.8,-81.0408], Tmin=(1022.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-31.543,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(C=CC(C)(O)CJ) + radical(C2CsJOCs)"""),
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
    label = 'C=C1C(C)CC1(O)C(=O)CC(25673)',
    structure = SMILES('C=C1C(C)CC1(O)C(=O)CC'),
    E0 = (-299.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[9.44757,0.040685,7.15202e-05,-9.23916e-08,2.17689e-11,-36347.9,-10.2701], Tmin=(100,'K'), Tmax=(1808.2,'K')), NASAPolynomial(coeffs=[101.784,0.0234542,-6.93375e-05,1.67441e-08,-1.22883e-12,-100316,-594.981], Tmin=(1808.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-299.41,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane)"""),
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
    E0 = (-94.891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-17.7057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-8.23227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (47.9493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (18.1955,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-70.6693,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-7.30068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (15.6508,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (74.8502,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (153.373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (87.6432,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-25.964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-8.28216,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (54.8962,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-37.7794,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-39.2438,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-42.1726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (9.96774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (176.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (136.325,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (30.8731,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (91.5342,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-25.2522,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (301.976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (350.388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (360.997,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (62.4274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (62.4274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (176.287,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-86.6066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (280.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (37.9224,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (212.624,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (42.0387,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (204.049,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (412.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (51.9241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (90.6895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (86.0411,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (159.271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (65.1664,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (-44.6205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (-65.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (21.2559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (30.8731,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (-31.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (-69.9177,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (104.997,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (-86.6066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (410.425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    products = ['C=C(O)C(=O)CC(4626)', 'CH3CHCCH2(18175)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    products = ['[CH2]C1([CH]C)CC1(O)C(=O)CC(25638)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.68e+09,'s^-1'), n=0.84, Ea=(77.1852,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_HNd;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 73.1 to 77.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    products = ['[CH2]C(=CC)C1(O)CC1([O])CC(25639)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.34238e+09,'s^-1'), n=0.889391, Ea=(86.6587,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonylbond_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 81.9 to 86.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    products = ['[CH2]C1(O)C(=CC)CC1([O])CC(25640)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.46159e+06,'s^-1'), n=1.55572, Ea=(142.84,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5_SS_CO;carbonylbond_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=C(O)C(=O)CC(4626)', 'C=[C][CH]C(18176)'],
    products = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00392164,'m^3/(mol*s)'), n=2.41519, Ea=(15.6067,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;CJ] for rate rule [Cds-COOs_Cds;CJ]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C2H5CO(71)', '[CH2]C(=CC)C(=C)O(25079)'],
    products = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.00392164,'m^3/(mol*s)'), n=2.41519, Ea=(15.6067,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;CJ] for rate rule [Cds-CdOs_Cds;CO_rad/NonDe]
Euclidian distance = 3.60555127546
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C(O)=C([O])CC(4557)', 'CH3CHCCH2(18175)'],
    products = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00086947,'m^3/(mol*s)'), n=2.67356, Ea=(32.0272,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ca_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['OH(5)', '[CH2]C(=CC)C(=C)C(=O)CC(25641)'],
    products = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(931.236,'m^3/(mol*s)'), n=1.015, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;OJ_pri] for rate rule [Cds-CdCO_Cds;OJ_pri]
Euclidian distance = 2.0
family: R_Addition_MultipleBond
Ea raised from -7.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(=CC)C(C)([O])C(=O)CC(25642)'],
    products = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)=[C]C(25643)'],
    products = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C([O])(C(=O)CC)C(C)=CC(25644)'],
    products = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(8.6325e+06,'s^-1'), n=1.395, Ea=(89.7468,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;O_rad_out;Cs_H_out_2H] + [R4H_SS(Cd)S;Y_rad_out;Cs_H_out_2H] for rate rule [R4H_SS(Cd)S;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(=[C]C)C(C)(O)C(=O)CC(25645)'],
    products = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    products = ['[CH2]C(=CC)C(C)(O)C(=O)[CH]C(25646)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(6.44e+09,'s^-1'), n=0.13, Ea=(86.6088,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 131 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(O)([C](C)C=C)C(=O)CC(19915)'],
    products = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(800000,'s^-1'), n=1.81, Ea=(149.787,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 101 used for R4H_SDS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SDS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    products = ['[CH2]C(O)(C(=O)[CH]C)C(C)=CC(25647)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(49679.4,'s^-1'), n=1.775, Ea=(57.1116,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5H_CCC;C_rad_out_2H;Cs_H_out_H/NonDeC] + [R5H_C(Cd)CC;C_rad_out_2H;Cs_H_out_1H] for rate rule [R5H_C(Cd)CC;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    products = ['[CH2][C](C=C)C(C)(O)C(=O)CC(19917)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(121000,'s^-1'), n=1.9, Ea=(55.6472,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 92 used for R5H_SSMS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R5H_SSMS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    products = ['[CH2]CC(=O)C(C)(O)C([CH2])=CC(25648)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(68850,'s^-1'), n=1.68, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_2H;Cs_H_out_2H] for rate rule [R5H_CCC(O2d)C;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]CC(=O)C([CH2])(O)C(C)=CC(25649)'],
    products = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(312,'s^-1'), n=2.1, Ea=(44.7688,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 91 used for R6H_SSSS(Cd)S;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R6H_SSSS(Cd)S;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(O)=C([O])CC(4557)', 'C=[C][CH]C(18176)'],
    products = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    products = ['[CH2]C(O)([C]1CC1C)C(=O)CC(25650)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_secNd_HNd;radadd_intra_cs2H]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    products = ['[CH2][C]1C(C)CC1(O)C(=O)CC(25651)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.51071e+08,'s^-1'), n=0.996667, Ea=(125.764,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    products = ['[CH2]C(=CC)C1(O)CO[C]1CC(25652)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7.01137e+09,'s^-1'), n=0.572544, Ea=(186.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4_S;multiplebond_intra;radadd_intra_cs2H] + [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    products = ['[CH2]C1(O)[C](CC)OCC1=CC(25653)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.64784e+07,'s^-1'), n=0.990488, Ea=(69.6388,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5_SS_CO;carbonyl_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic
Ea raised from 63.6 to 69.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction24',
    reactants = ['CH3CH2OH(54)', '[CH2]C(=C=O)C([CH2])=CC(25654)'],
    products = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.79e-11,'m^3/(mol*s)'), n=3.97, Ea=(329.281,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [doublebond;Cs_OH] for rate rule [Cdd_Cd;C_pri_OH]
Euclidian distance = 1.41421356237
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CH2(S)(23)', '[CH2]C(=CC)C([CH2])(O)C(C)=O(25655)'],
    products = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;C_pri/De]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['CH2(S)(23)', '[CH2]C(=C)C([CH2])(O)C(=O)CC(25656)'],
    products = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(7.94e+13,'cm^3/(mol*s)','*|/',0.25), n=-0.324, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 4 used for carbene;Cd_pri
Exact match found for rate rule [carbene;Cd_pri]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: 1,2_Insertion_carbene
Ea raised from -3.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    products = ['C=C([CH]C)C[C](O)C(=O)CC(25634)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    products = ['[CH2]C(=CC)[C](O)CC(=O)CC(25657)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ-HH;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;CO]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C(O)(C[C]=CC)C(=O)CC(19830)'],
    products = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    products = ['CC=C1CCC1(O)C(=O)CC(25658)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(=CC)C([CH2])=C(CC)OO(25659)'],
    products = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O] for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O_R]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=C([CH]C)O[C](CC)C(=C)O(25633)'],
    products = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O] for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O_C]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(=CC)C([CH2])(O)C(=C)OC(25660)'],
    products = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(37989.5,'s^-1'), n=2.515, Ea=(258.99,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C] + [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O] for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_C]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C(=CC)C([CH2])(O)C(O)=CC(25661)'],
    products = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O_H] for rate rule [R_ROR;R1_doublebond_CHCH3;R2_doublebond_CsC;R_O_H]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction35',
    reactants = ['CH2(19)', '[CH2]C(=CC)[C](O)C(=O)CC(25662)'],
    products = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/TDMustO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction36',
    reactants = ['CH2(19)', '[CH2]C(O)([C]=CC)C(=O)CC(25663)'],
    products = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    products = ['[CH2]C1(O)C(=C)C(C)C1([O])CC(25664)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2.13771e+06,'s^-1'), n=1.58803, Ea=(146.815,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_csHNd] for rate rule [R5_SS_CO;carbonylbond_intra_Nd;radadd_intra_csHNd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction38',
    reactants = ['H(3)', '[CH2]C(O)(C(=C)C=C)C(=O)CC(19912)'],
    products = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(2.31e+08,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2544 used for Cds-HH_Cds-CdH;HJ
Exact match found for rate rule [Cds-HH_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]CC(=C)C([CH2])(O)C(=O)CC(25665)'],
    products = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.72e+06,'s^-1'), n=1.99, Ea=(113.805,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 84 used for R2H_S;C_rad_out_2H;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]=C(CC)C([CH2])(O)C(=O)CC(25666)'],
    products = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 194 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]C([O])(C(=C)CC)C(=O)CC(25667)'],
    products = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(210000,'s^-1'), n=1.76, Ea=(53.8899,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4H_SSS;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R4H_SS(Cd)S;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH]=C([CH]C)C(C)(O)C(=O)CC(25668)'],
    products = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    products = ['[CH2]C(O)(C(=C)CC)C(=O)[CH]C(25669)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(0.1016,'s^-1'), n=3.24, Ea=(29.037,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_H/NonDeC;Cs_H_out_1H] for rate rule [R5H_C(Cd)CC;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2]CC(=O)C([CH2])(O)C(=C)CC(25670)'],
    products = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(25800,'s^-1'), n=1.67, Ea=(42.6768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_SSSSS;C_rad_out_2H;Cs_H_out_H/NonDeC] for rate rule [R6H_SSSS(Cd)S;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    products = ['C[CH][C]1CCC1(O)C(=O)CC(25671)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.51071e+08,'s^-1'), n=0.996667, Ea=(125.764,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    products = ['[CH2]C1(O)[C](CC)OC(C)C1=C(25672)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(4.64e+06,'s^-1'), n=1.15, Ea=(63.348,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_csHCs] for rate rule [R5_SS_CO;carbonyl_intra_Nd;radadd_intra_csHCs]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic
Ea raised from 58.2 to 63.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    products = ['C=CC(=C)C(C)(O)C(=O)CC(19932)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[C]=C(19918)'],
    products = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 5 used for cCs(-HC)CJ;CdsJ;C
Exact match found for rate rule [cCs(-HC)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    products = ['C=C1C(C)CC1(O)C(=O)CC(25673)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction50',
    reactants = ['CHCH3(T)(95)', '[CH2]C(O)([C]=C)C(=O)CC(9529)'],
    products = ['[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4302',
    isomers = [
        '[CH2]C(O)(C(=C)[CH]C)C(=O)CC(25635)',
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
    label = '4302',
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

