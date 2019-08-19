species(
    label = 'C=C(O)C([O])(CC)C(=C)[CH]C(25632)',
    structure = SMILES('[CH2]C(=CC)C([O])(CC)C(=C)O'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.24133,0.122242,-0.00010271,4.37536e-08,-7.43059e-12,-8563.24,44.6026], Tmin=(100,'K'), Tmax=(1410.71,'K')), NASAPolynomial(coeffs=[26.3925,0.0410524,-1.63822e-05,2.95775e-09,-2.00988e-13,-16642.1,-103.392], Tmin=(1410.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-73.1703,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(Allyl_P)"""),
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
    label = '[CH2]C1([CH]C)OC1(CC)C(=C)O(26036)',
    structure = SMILES('[CH2]C1([CH]C)OC1(CC)C(=C)O'),
    E0 = (21.521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.00938,0.145105,-0.000150046,7.7804e-08,-1.51251e-11,2903.4,47.6688], Tmin=(100,'K'), Tmax=(1462.36,'K')), NASAPolynomial(coeffs=[31.763,0.0243418,-2.66914e-06,-9.80743e-11,2.46242e-14,-5108.85,-130.13], Tmin=(1462.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(21.521,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCJCO) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]C(=CC)C1(CC)OC1([CH2])O(26080)',
    structure = SMILES('[CH2]C(=CC)C1(CC)OC1([CH2])O'),
    E0 = (-28.1545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.68875,0.14668,-0.000144503,7.01281e-08,-1.26374e-11,-3035.14,48.5144], Tmin=(100,'K'), Tmax=(1609.83,'K')), NASAPolynomial(coeffs=[33.802,0.0220808,-1.41949e-06,-3.00743e-10,3.53711e-14,-11675.3,-143.853], Tmin=(1609.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-28.1545,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Ethylene_oxide) + radical(Allyl_P) + radical(CJC(O)2C)"""),
)

species(
    label = '[CH2]C1(O)CC(=CC)C1([O])CC(25698)',
    structure = SMILES('[CH2]C1(O)CC(=CC)C1([O])CC'),
    E0 = (45.8055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[9.13945,0.0453452,6.46628e-05,-8.85869e-08,2.10516e-11,5184.8,-5.53165], Tmin=(100,'K'), Tmax=(1811.35,'K')), NASAPolynomial(coeffs=[103.555,0.0212895,-6.81555e-05,1.65126e-08,-1.21294e-12,-59276.7,-600.644], Tmin=(1811.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(45.8055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(methylenecyclobutane) + radical(CJC(C)2O) + radical(C=CC(C)2OJ)"""),
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
    label = '[CH2]C(=CC)C(=O)C(=C)O(26081)',
    structure = SMILES('[CH2]C(=CC)C(=O)C(=C)O'),
    E0 = (-152.73,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,375,552.5,462.5,1710,292.721,292.721],'cm^-1')),
        HinderedRotor(inertia=(0.266989,'amu*angstrom^2'), symmetry=1, barrier=(16.2345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.266992,'amu*angstrom^2'), symmetry=1, barrier=(16.2345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00196737,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.266993,'amu*angstrom^2'), symmetry=1, barrier=(16.2345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.266992,'amu*angstrom^2'), symmetry=1, barrier=(16.2345,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.343692,0.091552,-8.65445e-05,4.22651e-08,-8.28827e-12,-18209.5,30.2088], Tmin=(100,'K'), Tmax=(1224.31,'K')), NASAPolynomial(coeffs=[17.8811,0.0320089,-1.35934e-05,2.54137e-09,-1.76819e-13,-22672.1,-61.4037], Tmin=(1224.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-152.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-Cds)(Cds-Cds)) + group(Cds-CdsHH) + radical(C=C(C=O)CJ)"""),
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
    label = '[CH2]C(=CC)C(=O)CC(24970)',
    structure = SMILES('[CH2]C(=CC)C(=O)CC'),
    E0 = (-74.928,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,227.53,227.849],'cm^-1')),
        HinderedRotor(inertia=(0.13575,'amu*angstrom^2'), symmetry=1, barrier=(4.98767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00325207,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00325561,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135484,'amu*angstrom^2'), symmetry=1, barrier=(4.98597,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.32788,'amu*angstrom^2'), symmetry=1, barrier=(12.0546,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.162,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.90034,0.0754599,-5.94478e-05,2.89016e-08,-6.66614e-12,-8906.28,26.5503], Tmin=(100,'K'), Tmax=(948.49,'K')), NASAPolynomial(coeffs=[6.01178,0.053904,-2.53583e-05,4.94128e-09,-3.50812e-13,-9875.92,2.16079], Tmin=(948.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-74.928,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + radical(C=C(C=O)CJ)"""),
)

species(
    label = '[CH2]C(=CC)C(O)([CH]C)C(=C)O(26082)',
    structure = SMILES('[CH2]C(=CC)C(O)([CH]C)C(=C)O'),
    E0 = (-102.371,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.30415,0.12211,-9.43383e-05,2.6951e-08,5.45538e-13,-12071,46.0955], Tmin=(100,'K'), Tmax=(1032.03,'K')), NASAPolynomial(coeffs=[26.9401,0.0378991,-1.42898e-05,2.59715e-09,-1.81658e-13,-19658.8,-103.43], Tmin=(1032.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-102.371,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(CCJCO)"""),
)

species(
    label = 'C=C(O)C([O])(CC)C(C)=[C]C(26083)',
    structure = SMILES('C=C(O)C([O])(CC)C(C)=[C]C'),
    E0 = (13.1723,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1277.5,1000,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,2950,3100,1380,975,1025,1650,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,180,180,180,520.958,621.408,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157729,'amu*angstrom^2'), symmetry=1, barrier=(3.62651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157729,'amu*angstrom^2'), symmetry=1, barrier=(3.62651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157729,'amu*angstrom^2'), symmetry=1, barrier=(3.62651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157729,'amu*angstrom^2'), symmetry=1, barrier=(3.62651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157729,'amu*angstrom^2'), symmetry=1, barrier=(3.62651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157729,'amu*angstrom^2'), symmetry=1, barrier=(3.62651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157729,'amu*angstrom^2'), symmetry=1, barrier=(3.62651,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.52872,0.120778,-0.000107237,5.06347e-08,-9.84649e-12,1783.87,42.2246], Tmin=(100,'K'), Tmax=(1212.33,'K')), NASAPolynomial(coeffs=[19.0822,0.0527729,-2.30949e-05,4.3641e-09,-3.04758e-13,-3213.52,-61.1793], Tmin=(1212.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(13.1723,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(Cds_S)"""),
)

species(
    label = 'C=C(O)C([O])([CH]C)C(C)=CC(26084)',
    structure = SMILES('C=C(O)C([O])([CH]C)C(C)=CC'),
    E0 = (-24.7675,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,325,375,415,465,420,450,1700,1750,180,180,180,467.361,670.276,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156484,'amu*angstrom^2'), symmetry=1, barrier=(3.59788,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156484,'amu*angstrom^2'), symmetry=1, barrier=(3.59788,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156484,'amu*angstrom^2'), symmetry=1, barrier=(3.59788,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156484,'amu*angstrom^2'), symmetry=1, barrier=(3.59788,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156484,'amu*angstrom^2'), symmetry=1, barrier=(3.59788,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156484,'amu*angstrom^2'), symmetry=1, barrier=(3.59788,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156484,'amu*angstrom^2'), symmetry=1, barrier=(3.59788,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.89479,0.117816,-9.70415e-05,4.08702e-08,-6.92039e-12,-2756.97,45.5693], Tmin=(100,'K'), Tmax=(1404.22,'K')), NASAPolynomial(coeffs=[24.1521,0.0436201,-1.77851e-05,3.2426e-09,-2.21398e-13,-10072.1,-88.9344], Tmin=(1404.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-24.7675,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCJCO) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH2]CC(O)(C(=C)O)C([CH2])=CC(26085)',
    structure = SMILES('[CH2]CC(O)(C(=C)O)C([CH2])=CC'),
    E0 = (-97.0266,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.65787,0.130216,-0.00012016,5.61741e-08,-1.03439e-11,-11416.3,46.6645], Tmin=(100,'K'), Tmax=(1320.08,'K')), NASAPolynomial(coeffs=[28.2584,0.0365352,-1.37107e-05,2.41471e-09,-1.6276e-13,-19578.6,-111.074], Tmin=(1320.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-97.0266,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=[C]C)C(O)(CC)C(=C)O(26086)',
    structure = SMILES('[CH2]C(=[C]C)C(O)(CC)C(=C)O'),
    E0 = (-64.4312,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3580,3650,1210,1345,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,180,1555.39,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.150222,'amu*angstrom^2'), symmetry=1, barrier=(3.45389,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150222,'amu*angstrom^2'), symmetry=1, barrier=(3.45389,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150222,'amu*angstrom^2'), symmetry=1, barrier=(3.45389,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150222,'amu*angstrom^2'), symmetry=1, barrier=(3.45389,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150222,'amu*angstrom^2'), symmetry=1, barrier=(3.45389,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150222,'amu*angstrom^2'), symmetry=1, barrier=(3.45389,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150222,'amu*angstrom^2'), symmetry=1, barrier=(3.45389,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150222,'amu*angstrom^2'), symmetry=1, barrier=(3.45389,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.45192,0.13093,-0.000123957,6.01856e-08,-1.1597e-11,-7507.8,44.6059], Tmin=(100,'K'), Tmax=(1256.85,'K')), NASAPolynomial(coeffs=[25.9976,0.0403878,-1.58975e-05,2.86796e-09,-1.95914e-13,-14659.1,-99.1504], Tmin=(1256.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-64.4312,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=CC)C(O)(CC)C(=C)[O](26087)',
    structure = SMILES('[CH2]C(=CC)C(O)(CC)C(=C)[O]'),
    E0 = (-164.468,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.08048,0.122653,-0.000107553,4.88617e-08,-8.90069e-12,-19552.9,44.3018], Tmin=(100,'K'), Tmax=(1316.7,'K')), NASAPolynomial(coeffs=[24.0399,0.0433018,-1.71555e-05,3.09181e-09,-2.10405e-13,-26431.4,-88.9004], Tmin=(1316.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-164.468,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C(O)C(O)(CC)C([CH2])=CC(26088)',
    structure = SMILES('[CH]=C(O)C(O)(CC)C([CH2])=CC'),
    E0 = (-55.1768,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.72576,0.132532,-0.000124879,5.96073e-08,-1.1196e-11,-6381.17,45.5224], Tmin=(100,'K'), Tmax=(1295.88,'K')), NASAPolynomial(coeffs=[28.3858,0.0364998,-1.37194e-05,2.42071e-09,-1.63495e-13,-14444.5,-112.636], Tmin=(1295.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-55.1768,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_P)"""),
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
    label = '[CH2]CC([O])(C(=C)O)C(C)=CC(26089)',
    structure = SMILES('[CH2]CC([O])(C(=C)O)C(C)=CC'),
    E0 = (-19.4232,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,445.023,693.583,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156052,'amu*angstrom^2'), symmetry=1, barrier=(3.58794,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156052,'amu*angstrom^2'), symmetry=1, barrier=(3.58794,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156052,'amu*angstrom^2'), symmetry=1, barrier=(3.58794,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156052,'amu*angstrom^2'), symmetry=1, barrier=(3.58794,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156052,'amu*angstrom^2'), symmetry=1, barrier=(3.58794,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156052,'amu*angstrom^2'), symmetry=1, barrier=(3.58794,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156052,'amu*angstrom^2'), symmetry=1, barrier=(3.58794,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.6889,0.119619,-0.000102258,4.54893e-08,-8.23717e-12,-2126.91,44.1116], Tmin=(100,'K'), Tmax=(1305.59,'K')), NASAPolynomial(coeffs=[21.4652,0.048681,-2.0757e-05,3.87296e-09,-2.68324e-13,-8172.87,-73.7677], Tmin=(1305.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-19.4232,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(RCCJ)"""),
)

species(
    label = 'C=C([O])C([O])(CC)C(C)=CC(26090)',
    structure = SMILES('C=C([O])C([O])(CC)C(C)=CC'),
    E0 = (-86.8647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.11774,0.112114,-8.97963e-05,3.83062e-08,-6.8315e-12,-10263.2,41.7725], Tmin=(100,'K'), Tmax=(1294.97,'K')), NASAPolynomial(coeffs=[17.2034,0.0555234,-2.42465e-05,4.56077e-09,-3.16867e-13,-15008.3,-51.3522], Tmin=(1294.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-86.8647,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH]=C(O)C([O])(CC)C(C)=CC(26091)',
    structure = SMILES('[CH]=C(O)C([O])(CC)C(C)=CC'),
    E0 = (22.4266,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,180,180,180,180,1600,1642.84,2856.38,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152433,'amu*angstrom^2'), symmetry=1, barrier=(3.50474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152433,'amu*angstrom^2'), symmetry=1, barrier=(3.50474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152433,'amu*angstrom^2'), symmetry=1, barrier=(3.50474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152433,'amu*angstrom^2'), symmetry=1, barrier=(3.50474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152433,'amu*angstrom^2'), symmetry=1, barrier=(3.50474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152433,'amu*angstrom^2'), symmetry=1, barrier=(3.50474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152433,'amu*angstrom^2'), symmetry=1, barrier=(3.50474,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.76993,0.122061,-0.000107303,4.92255e-08,-9.18146e-12,2908.89,43.019], Tmin=(100,'K'), Tmax=(1271.78,'K')), NASAPolynomial(coeffs=[21.5426,0.0487387,-2.08229e-05,3.89304e-09,-2.70262e-13,-3020.81,-75.0553], Tmin=(1271.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(22.4266,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CC(C)2OJ)"""),
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
    label = 'C=C(O)C([O])(CC)[C]1CC1C(26092)',
    structure = SMILES('C=C(O)C([O])(CC)[C]1CC1C'),
    E0 = (-21.6383,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.43913,0.0949353,-1.14848e-05,-6.14156e-08,3.26025e-11,-2384.42,42.833], Tmin=(100,'K'), Tmax=(966.649,'K')), NASAPolynomial(coeffs=[26.8789,0.0366458,-1.24176e-05,2.25209e-09,-1.63324e-13,-10610.5,-107.056], Tmin=(966.649,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-21.6383,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(C=CC(C)2OJ) + radical(CCJ(C)CO)"""),
)

species(
    label = '[CH2][C](O)C1(CC)OC(C)C1=C(25966)',
    structure = SMILES('[CH2][C](O)C1(CC)OC(C)C1=C'),
    E0 = (44.2933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.98741,0.114198,-7.40209e-05,8.95309e-09,6.12078e-12,5558.07,43.3627], Tmin=(100,'K'), Tmax=(1024.05,'K')), NASAPolynomial(coeffs=[25.4462,0.0406157,-1.54196e-05,2.81979e-09,-1.98129e-13,-1821.06,-98.236], Tmin=(1024.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(44.2933,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(C2CsJOH) + radical(CJCO)"""),
)

species(
    label = '[CH2]C(=CC)C1(CC)OC[C]1O(26093)',
    structure = SMILES('[CH2]C(=CC)C1(CC)OC[C]1O'),
    E0 = (-21.1013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.46129,0.123397,-0.000111702,5.46782e-08,-1.05579e-11,-2289.38,43.589], Tmin=(100,'K'), Tmax=(1357.16,'K')), NASAPolynomial(coeffs=[23.6623,0.039498,-1.13407e-05,1.62996e-09,-9.54743e-14,-8744.27,-88.0773], Tmin=(1357.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-21.1013,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Oxetane) + radical(Allyl_P) + radical(C2CsJOH)"""),
)

species(
    label = 'CC=C1CC[C](O)C1([O])CC(25762)',
    structure = SMILES('CC=C1CC[C](O)C1([O])CC'),
    E0 = (-55.0371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.27805,0.100472,-5.48324e-05,5.15634e-09,3.57206e-12,-6416.37,41.4735], Tmin=(100,'K'), Tmax=(1136.67,'K')), NASAPolynomial(coeffs=[19.9725,0.0500142,-2.03445e-05,3.75522e-09,-2.6047e-13,-12818.7,-70.6826], Tmin=(1136.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-55.0371,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(590.328,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(methylenecyclopentane) + radical(C=CC(C)2OJ) + radical(C2CsJOH)"""),
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
    label = '[CH2]C(=CC)C(C)([O])C(=C)O(26094)',
    structure = SMILES('[CH2]C(=CC)C(C)([O])C(=C)O'),
    E0 = (-49.3901,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.68094,0.108385,-9.17425e-05,3.90015e-08,-6.55559e-12,-5721.88,40.3522], Tmin=(100,'K'), Tmax=(1434.84,'K')), NASAPolynomial(coeffs=[25.4614,0.0327168,-1.26362e-05,2.24565e-09,-1.5128e-13,-13510.7,-100.393], Tmin=(1434.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-49.3901,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH2]C(=C)C([O])(CC)C(=C)O(26095)',
    structure = SMILES('[CH2]C(=C)C([O])(CC)C(=C)O'),
    E0 = (-37.1447,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.48386,0.105539,-7.91202e-05,2.28198e-08,-2.53845e-13,-4257.18,39.3735], Tmin=(100,'K'), Tmax=(1080.51,'K')), NASAPolynomial(coeffs=[23.7151,0.0354866,-1.41241e-05,2.61728e-09,-1.83728e-13,-11059,-90.4238], Tmin=(1080.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-37.1447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C(O)C([O])(CC)C[C]=CC(20179)',
    structure = SMILES('C=C(O)C([O])(CC)C[C]=CC'),
    E0 = (34.6531,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,180,1085.24,1133.56,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.162437,'amu*angstrom^2'), symmetry=1, barrier=(3.73475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162437,'amu*angstrom^2'), symmetry=1, barrier=(3.73475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162437,'amu*angstrom^2'), symmetry=1, barrier=(3.73475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162437,'amu*angstrom^2'), symmetry=1, barrier=(3.73475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162437,'amu*angstrom^2'), symmetry=1, barrier=(3.73475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162437,'amu*angstrom^2'), symmetry=1, barrier=(3.73475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162437,'amu*angstrom^2'), symmetry=1, barrier=(3.73475,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.94699,0.123714,-0.000112923,5.42935e-08,-1.0529e-11,4387.81,45.9784], Tmin=(100,'K'), Tmax=(1236.57,'K')), NASAPolynomial(coeffs=[22.1353,0.0458134,-1.8426e-05,3.3475e-09,-2.29054e-13,-1568.05,-75.3184], Tmin=(1236.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(34.6531,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(Cds_S)"""),
)

species(
    label = 'C=C(O)C1(CC)OCC1=CC(26017)',
    structure = SMILES('C=C(O)C1(CC)OCC1=CC'),
    E0 = (-258.145,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.34032,0.094616,-1.84462e-05,-4.63342e-08,2.47758e-11,-30835.2,39.1578], Tmin=(100,'K'), Tmax=(1003.2,'K')), NASAPolynomial(coeffs=[24.7633,0.0417917,-1.61024e-05,3.03811e-09,-2.19731e-13,-38652,-99.7157], Tmin=(1003.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-258.145,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutane)"""),
)

species(
    label = '[CH2]C(=CC)C([C]=C)(CC)OO(26096)',
    structure = SMILES('[CH2]C(=CC)C([C]=C)(CC)OO'),
    E0 = (221.632,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.47449,0.121457,-0.000105734,4.89021e-08,-9.3907e-12,26852.3,41.8015], Tmin=(100,'K'), Tmax=(1217.48,'K')), NASAPolynomial(coeffs=[18.4142,0.0561122,-2.5225e-05,4.81659e-09,-3.37981e-13,22009.5,-58.0634], Tmin=(1217.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(221.632,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=CC)C([O])(CC)C(C)=O(26097)',
    structure = SMILES('[CH2]C(=CC)C([O])(CC)C(C)=O'),
    E0 = (-62.1912,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.64708,0.132942,-0.000159895,1.16686e-07,-3.59053e-11,-7284.53,41.4441], Tmin=(100,'K'), Tmax=(781.828,'K')), NASAPolynomial(coeffs=[11.6148,0.0650826,-2.96849e-05,5.6407e-09,-3.92566e-13,-9357.97,-19.2711], Tmin=(781.828,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-62.1912,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-O2d)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + radical(C=CC(C)(C=O)OJ) + radical(Allyl_P)"""),
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
    label = '[CH2]C(=CC)[C](CC)C(=C)O(26098)',
    structure = SMILES('[CH2]C([CH]C)=C(CC)C(=C)O'),
    E0 = (-18.7655,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,325,350,375,415,440,465,420,435,450,1700,1725,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.58798,0.105749,-5.82926e-05,-7.33959e-09,1.27203e-11,-2040.44,36.7334], Tmin=(100,'K'), Tmax=(973.801,'K')), NASAPolynomial(coeffs=[24.2374,0.0380844,-1.32389e-05,2.32747e-09,-1.61737e-13,-8891.67,-96.526], Tmin=(973.801,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-18.7655,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(Allyl_S) + radical(C=CC=CCJ)"""),
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
    label = 'C=C(O)C([O])([C]=CC)CC(26099)',
    structure = SMILES('C=C(O)C([O])([C]=CC)CC'),
    E0 = (52.2273,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,489.64,651.799,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157218,'amu*angstrom^2'), symmetry=1, barrier=(3.61476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157218,'amu*angstrom^2'), symmetry=1, barrier=(3.61476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157218,'amu*angstrom^2'), symmetry=1, barrier=(3.61476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157218,'amu*angstrom^2'), symmetry=1, barrier=(3.61476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157218,'amu*angstrom^2'), symmetry=1, barrier=(3.61476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157218,'amu*angstrom^2'), symmetry=1, barrier=(3.61476,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.986928,0.104324,-8.89923e-05,3.9313e-08,-7.05326e-12,6465.27,39.5227], Tmin=(100,'K'), Tmax=(1318.71,'K')), NASAPolynomial(coeffs=[19.648,0.0417333,-1.77981e-05,3.32157e-09,-2.30119e-13,1022.91,-65.7378], Tmin=(1318.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(52.2273,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(Cds_S)"""),
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
    label = '[CH2]CC(=C)C([O])(CC)C(=C)O(26100)',
    structure = SMILES('[CH2]CC(=C)C([O])(CC)C(=C)O'),
    E0 = (-6.04305,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,180,180,180,180,921.865,1284.89,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.158177,'amu*angstrom^2'), symmetry=1, barrier=(3.6368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158177,'amu*angstrom^2'), symmetry=1, barrier=(3.6368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158177,'amu*angstrom^2'), symmetry=1, barrier=(3.6368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158177,'amu*angstrom^2'), symmetry=1, barrier=(3.6368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158177,'amu*angstrom^2'), symmetry=1, barrier=(3.6368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158177,'amu*angstrom^2'), symmetry=1, barrier=(3.6368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158177,'amu*angstrom^2'), symmetry=1, barrier=(3.6368,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.16636,0.121322,-0.000102145,4.36724e-08,-7.45113e-12,-493.011,46.9442], Tmin=(100,'K'), Tmax=(1403.4,'K')), NASAPolynomial(coeffs=[25.9498,0.0411853,-1.64928e-05,2.98488e-09,-2.03159e-13,-8384.71,-98.2286], Tmin=(1403.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-6.04305,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(RCCJ) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH]=C(CC)C([O])(CC)C(=C)O(26101)',
    structure = SMILES('[CH]=C(CC)C([O])(CC)C(=C)O'),
    E0 = (35.8067,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,414.908,720.221,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155333,'amu*angstrom^2'), symmetry=1, barrier=(3.5714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155333,'amu*angstrom^2'), symmetry=1, barrier=(3.5714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155333,'amu*angstrom^2'), symmetry=1, barrier=(3.5714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155333,'amu*angstrom^2'), symmetry=1, barrier=(3.5714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155333,'amu*angstrom^2'), symmetry=1, barrier=(3.5714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155333,'amu*angstrom^2'), symmetry=1, barrier=(3.5714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155333,'amu*angstrom^2'), symmetry=1, barrier=(3.5714,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.22339,0.123525,-0.000106531,4.67569e-08,-8.18584e-12,4541.63,45.7622], Tmin=(100,'K'), Tmax=(1370.26,'K')), NASAPolynomial(coeffs=[25.9063,0.0414089,-1.66388e-05,3.02127e-09,-2.06292e-13,-3167.26,-98.8078], Tmin=(1370.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(35.8067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(Cds_P)"""),
)

species(
    label = 'C=C(O)C([O])([CH]C)C(=C)CC(26102)',
    structure = SMILES('C=C(O)C([O])([CH]C)C(=C)CC'),
    E0 = (-11.3874,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,180,180,180,180,968.043,1236.82,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.159239,'amu*angstrom^2'), symmetry=1, barrier=(3.66122,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159239,'amu*angstrom^2'), symmetry=1, barrier=(3.66122,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159239,'amu*angstrom^2'), symmetry=1, barrier=(3.66122,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159239,'amu*angstrom^2'), symmetry=1, barrier=(3.66122,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159239,'amu*angstrom^2'), symmetry=1, barrier=(3.66122,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159239,'amu*angstrom^2'), symmetry=1, barrier=(3.66122,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159239,'amu*angstrom^2'), symmetry=1, barrier=(3.66122,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.99454,0.11533,-8.3569e-05,2.36237e-08,-3.66959e-13,-1139.76,47.03], Tmin=(100,'K'), Tmax=(1101.82,'K')), NASAPolynomial(coeffs=[24.7767,0.0422378,-1.68685e-05,3.11547e-09,-2.17533e-13,-8501.89,-91.3594], Tmin=(1101.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-11.3874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCJCO)"""),
)

species(
    label = '[CH]=C([CH]C)C(O)(CC)C(=C)O(26103)',
    structure = SMILES('[CH]C(=CC)C(O)(CC)C(=C)O'),
    E0 = (-83.0877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.58088,0.130573,-0.000113278,5.07479e-08,-9.09799e-12,-9744.45,45.4285], Tmin=(100,'K'), Tmax=(1339.7,'K')), NASAPolynomial(coeffs=[26.049,0.0450918,-1.75679e-05,3.12039e-09,-2.1027e-13,-17415.6,-101.067], Tmin=(1339.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-83.0877,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]CC([O])(C(=C)O)C(=C)CC(26104)',
    structure = SMILES('[CH2]CC([O])(C(=C)O)C(=C)CC'),
    E0 = (-6.04305,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,180,180,180,180,921.865,1284.89,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.158177,'amu*angstrom^2'), symmetry=1, barrier=(3.6368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158177,'amu*angstrom^2'), symmetry=1, barrier=(3.6368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158177,'amu*angstrom^2'), symmetry=1, barrier=(3.6368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158177,'amu*angstrom^2'), symmetry=1, barrier=(3.6368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158177,'amu*angstrom^2'), symmetry=1, barrier=(3.6368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158177,'amu*angstrom^2'), symmetry=1, barrier=(3.6368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158177,'amu*angstrom^2'), symmetry=1, barrier=(3.6368,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.16636,0.121322,-0.000102145,4.36724e-08,-7.45113e-12,-493.011,46.9442], Tmin=(100,'K'), Tmax=(1403.4,'K')), NASAPolynomial(coeffs=[25.9498,0.0411853,-1.64928e-05,2.98488e-09,-2.03159e-13,-8384.71,-98.2286], Tmin=(1403.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-6.04305,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(RCCJ)"""),
)

species(
    label = 'C=C([O])C([O])(CC)C(=C)CC(26105)',
    structure = SMILES('C=C([O])C([O])(CC)C(=C)CC'),
    E0 = (-73.4846,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.57833,0.113666,-8.93303e-05,3.61995e-08,-5.96858e-12,-8630.18,44.5406], Tmin=(100,'K'), Tmax=(1423.58,'K')), NASAPolynomial(coeffs=[21.9184,0.047645,-1.9765e-05,3.62188e-09,-2.47519e-13,-15320.1,-77.116], Tmin=(1423.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-73.4846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C(O)C([O])(CC)C(=C)CC(26106)',
    structure = SMILES('[CH]=C(O)C([O])(CC)C(=C)CC'),
    E0 = (35.8067,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,414.908,720.221,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155333,'amu*angstrom^2'), symmetry=1, barrier=(3.5714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155333,'amu*angstrom^2'), symmetry=1, barrier=(3.5714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155333,'amu*angstrom^2'), symmetry=1, barrier=(3.5714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155333,'amu*angstrom^2'), symmetry=1, barrier=(3.5714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155333,'amu*angstrom^2'), symmetry=1, barrier=(3.5714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155333,'amu*angstrom^2'), symmetry=1, barrier=(3.5714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155333,'amu*angstrom^2'), symmetry=1, barrier=(3.5714,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.22339,0.123525,-0.000106531,4.67569e-08,-8.18584e-12,4541.63,45.7622], Tmin=(100,'K'), Tmax=(1370.26,'K')), NASAPolynomial(coeffs=[25.9063,0.0414089,-1.66388e-05,3.02127e-09,-2.06292e-13,-3167.26,-98.8078], Tmin=(1370.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(35.8067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C](O)C1(CC)OCC1=CC(25993)',
    structure = SMILES('[CH2][C](O)C1(CC)OCC1=CC'),
    E0 = (50.5841,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.81122,0.117499,-9.81763e-05,4.29648e-08,-7.60754e-12,6301.53,43.502], Tmin=(100,'K'), Tmax=(1343.86,'K')), NASAPolynomial(coeffs=[22.2822,0.045786,-1.81326e-05,3.2571e-09,-2.20801e-13,-174.216,-79.8558], Tmin=(1343.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(50.5841,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclobutane) + radical(C2CsJOH) + radical(CJCO)"""),
)

species(
    label = 'C=C1C(C)C[C](O)C1([O])CC(25729)',
    structure = SMILES('C=C1C(C)C[C](O)C1([O])CC'),
    E0 = (-50.7636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.21013,0.0963624,-3.54222e-05,-2.01503e-08,1.36644e-11,-5902.14,40.0287], Tmin=(100,'K'), Tmax=(1040.4,'K')), NASAPolynomial(coeffs=[21.0406,0.0478696,-1.89304e-05,3.51433e-09,-2.47727e-13,-12537.5,-77.8375], Tmin=(1040.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-50.7636,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(590.328,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclopentane) + radical(C2CsJOH) + radical(C=CC(C)2OJ)"""),
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
    label = 'C=C(O)C1(CC)OC(C)C1=C(25992)',
    structure = SMILES('C=C(O)C1(CC)OC(C)C1=C'),
    E0 = (-264.436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.43785,0.090325,9.56518e-06,-8.60808e-08,4.13231e-11,-31582.2,38.7392], Tmin=(100,'K'), Tmax=(975.516,'K')), NASAPolynomial(coeffs=[29.2208,0.0345655,-1.226e-05,2.34335e-09,-1.76271e-13,-40892.3,-125.471], Tmin=(975.516,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-264.436,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutane)"""),
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
    label = 'C=[C]C([O])(CC)C(=C)O(9633)',
    structure = SMILES('C=[C]C([O])(CC)C(=C)O'),
    E0 = (88.2529,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,420.019,717.559,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155273,'amu*angstrom^2'), symmetry=1, barrier=(3.57004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155273,'amu*angstrom^2'), symmetry=1, barrier=(3.57004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155273,'amu*angstrom^2'), symmetry=1, barrier=(3.57004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155273,'amu*angstrom^2'), symmetry=1, barrier=(3.57004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155273,'amu*angstrom^2'), symmetry=1, barrier=(3.57004,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.628126,0.092044,-7.95451e-05,3.48078e-08,-6.07079e-12,10789,36.4352], Tmin=(100,'K'), Tmax=(1374.89,'K')), NASAPolynomial(coeffs=[20.5059,0.0305585,-1.24648e-05,2.28149e-09,-1.56468e-13,4977.61,-72.2528], Tmin=(1374.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(88.2529,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CC(C)2OJ)"""),
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
    E0 = (-73.1703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (23.2849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-26.9022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (51.1647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-16.8232,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (34.8055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (86.565,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-7.30068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (3.19106,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (175.093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (58.2326,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-13.3466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-20.1227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (28.6373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-10.8683,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (76.6169,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (23.2536,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-30.4935,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (55.4668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (16.0535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (176.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (158.046,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (53.0652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (51.7221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (41.8897,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (370.472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (382.717,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (204.691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-64.886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (328.743,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (131.009,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (224.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (433.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (55.1395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (112.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (107.762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (180.992,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (35.3479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (346.237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (40.3994,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (-44.1333,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (74.0295,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (53.0652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (-1.20554,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (-48.1971,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (133.401,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (-64.886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (432.145,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    products = ['C=C(O)C(=O)CC(4626)', 'CH3CHCCH2(18175)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    products = ['[CH2]C1([CH]C)OC1(CC)C(=C)O(26036)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.68e+09,'s^-1'), n=0.84, Ea=(96.4553,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_HNd;radadd_intra] for rate rule [R4_S_D;doublebond_intra_HNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    products = ['[CH2]C(=CC)C1(CC)OC1([CH2])O(26080)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.52e+08,'s^-1'), n=0.89, Ea=(46.2682,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H_secNd;radadd_intra] for rate rule [R4_S_D;doublebond_intra_2H_secNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    products = ['[CH2]C1(O)CC(=CC)C1([O])CC(25698)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.48e+06,'s^-1'), n=1.25, Ea=(124.335,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 349 used for R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C2H5(29)', '[CH2]C(=CC)C(=O)C(=C)O(26081)'],
    products = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7.94e+10,'cm^3/(mol*s)'), n=0, Ea=(28.0328,'kJ/mol'), T0=(1,'K'), Tmin=(333,'K'), Tmax=(363,'K'), comment="""Estimated using template [CO_O;CsJ-CsHH] for rate rule [CO-DeDe_O;CsJ-CsHH]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=C(O)C(=O)CC(4626)', 'C=[C][CH]C(18176)'],
    products = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.97e+07,'cm^3/(mol*s)'), n=1.88, Ea=(32.2168,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO-CdCs_O;YJ] for rate rule [CO-CdCs_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2COH(99)', '[CH2]C(=CC)C(=O)CC(24970)'],
    products = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2339.95,'m^3/(mol*s)'), n=0.573452, Ea=(58.2237,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;CJ] for rate rule [CO-CdCs_O;CdsJ-O2s]
Euclidian distance = 2.82842712475
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C(O)=C([O])CC(4557)', 'CH3CHCCH2(18175)'],
    products = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.00086947,'m^3/(mol*s)'), n=2.67356, Ea=(32.0272,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ca_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(=CC)C(O)([CH]C)C(=C)O(26082)'],
    products = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.71,'s^-1'), n=3.021, Ea=(105.562,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 319 used for R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=C(O)C([O])(CC)C(C)=[C]C(26083)'],
    products = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=C(O)C([O])([CH]C)C(C)=CC(26084)'],
    products = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1213.05,'s^-1'), n=2.57583, Ea=(83.0001,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;C_rad_out_1H;Cs_H_out_2H] + [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out] + [R4H_SS(Cd)S;C_rad_out_1H;Cs_H_out] for rate rule [R4H_SS(Cd)S;C_rad_out_H/NonDeC;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]CC(O)(C(=C)O)C([CH2])=CC(26085)'],
    products = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(=[C]C)C(O)(CC)C(=C)O(26086)'],
    products = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;XH_out] for rate rule [R4H_DSS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    products = ['[CH2]C(=CC)C(O)(CC)C(=C)[O](26087)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1070.11,'s^-1'), n=2.50856, Ea=(101.808,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;O_H_out] + [R4H_SSS;O_rad_out;XH_out] + [R4H_SS(Cd)S;Y_rad_out;XH_out] for rate rule [R4H_SS(Cd)S;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=C(O)C(O)(CC)C([CH2])=CC(26088)'],
    products = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=C[C](C)C([O])(CC)C(=C)O(20149)'],
    products = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(800000,'s^-1'), n=1.81, Ea=(149.787,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 101 used for R4H_SDS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SDS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]CC([O])(C(=C)O)C(C)=CC(26089)'],
    products = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(6900,'s^-1'), n=1.98, Ea=(42.6768,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 90 used for R5H_CCC(Cd);C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R5H_CCC(Cd);C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    products = ['C=C([O])C([O])(CC)C(C)=CC(26090)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2300,'s^-1'), n=1.98, Ea=(42.6768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC(Cd);C_rad_out_2H;XH_out] for rate rule [R5H_CCC(Cd);C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C(O)C([O])(CC)C(C)=CC(26091)'],
    products = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5H;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    products = ['[CH2][C](C=C)C(O)(CC)C(=C)O(20150)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.45388e+06,'s^-1'), n=1.705, Ea=(89.2238,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSMS;Y_rad_out;Cs_H_out_2H] for rate rule [R5H_SSMS;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(O)=C([O])CC(4557)', 'C=[C][CH]C(18176)'],
    products = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    products = ['C=C(O)C([O])(CC)[C]1CC1C(26092)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_secNd_HNd;radadd_intra_cs2H]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    products = ['[CH2][C](O)C1(CC)OC(C)C1=C(25966)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    products = ['[CH2]C(=CC)C1(CC)OC[C]1O(26093)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.21748e+08,'s^-1'), n=0.95, Ea=(124.892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R4_S_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    products = ['CC=C1CC[C](O)C1([O])CC(25762)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.71e+11,'s^-1'), n=0.2, Ea=(115.06,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 4 used for R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CH2(S)(23)', '[CH2]C(=CC)C(C)([O])C(=C)O(26094)'],
    products = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['CH2(S)(23)', '[CH2]C(=C)C([O])(CC)C(=C)O(26095)'],
    products = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(7.94e+13,'cm^3/(mol*s)','*|/',0.25), n=-0.324, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 4 used for carbene;Cd_pri
Exact match found for rate rule [carbene;Cd_pri]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: 1,2_Insertion_carbene
Ea raised from -3.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    products = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    products = ['C=C(O)C1(CC)OCC1=CC(26017)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R4_SSS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(=CC)C([C]=C)(CC)OO(26096)'],
    products = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(3.01978e+11,'s^-1'), n=0, Ea=(107.111,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2OOH_S;Y_rad_out] for rate rule [R2OOH_S;Cd_rad_out_double]
Euclidian distance = 2.0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    products = ['[CH2]C(=CC)C([O])(CC)C(C)=O(26097)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(204.179,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction32',
    reactants = ['O(4)', '[CH2]C(=CC)[C](CC)C(=C)O(26098)'],
    products = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/Cs;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction33',
    reactants = ['CH2(19)', 'C=C(O)C([O])([C]=CC)CC(26099)'],
    products = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    products = ['[CH2]C1(O)C(C)C(=C)C1([O])CC(25675)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(5.53628e+06,'s^-1'), n=1.28, Ea=(128.31,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5_SS_D;doublebond_intra_2H;radadd_intra_csHNd] + [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_cs] for rate rule [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_csHNd]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction35',
    reactants = ['H(3)', 'C=CC(=C)C([O])(CC)C(=C)O(20146)'],
    products = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2.31e+08,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2544 used for Cds-HH_Cds-CdH;HJ
Exact match found for rate rule [Cds-HH_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]CC(=C)C([O])(CC)C(=C)O(26100)'],
    products = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.72e+06,'s^-1'), n=1.99, Ea=(113.805,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 84 used for R2H_S;C_rad_out_2H;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]=C(CC)C([O])(CC)C(=C)O(26101)'],
    products = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 194 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C=C(O)C([O])([CH]C)C(=C)CC(26102)'],
    products = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(0.00228,'s^-1'), n=3.95, Ea=(46.7353,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC] for rate rule [R4H_SS(Cd)S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    products = ['[CH]=C([CH]C)C(O)(CC)C(=C)O(26103)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(3.4207e+12,'s^-1'), n=1.11009, Ea=(419.408,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSD;Y_rad_out;Cd_H_out_singleH] for rate rule [R4H_SSD;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]CC([O])(C(=C)O)C(=C)CC(26104)'],
    products = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(262000,'s^-1'), n=1.62, Ea=(46.4424,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_2H;Cs_H_out_H/NonDeC] for rate rule [R5H_CCC(Cd);C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    products = ['C=C([O])C([O])(CC)C(=C)CC(26105)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(0.0508,'s^-1'), n=3.24, Ea=(29.037,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_H/NonDeC;XH_out] for rate rule [R5H_CCC(Cd);C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH]=C(O)C([O])(CC)C(=C)CC(26106)'],
    products = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(315594,'s^-1'), n=1.73223, Ea=(38.2227,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5H;Cd_rad_out_singleH;Cs_H_out] + [R5H_RSSR;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 3.60555127546
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    products = ['[CH2][C](O)C1(CC)OCC1=CC(25993)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction44',
    reactants = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    products = ['C=C1C(C)C[C](O)C1([O])CC(25729)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(9.51884e+07,'s^-1'), n=0.875, Ea=(71.9648,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5_SS_D;doublebond_intra;radadd_intra_csHCs] + [R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_cs] for rate rule [R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction45',
    reactants = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    products = ['C=CC(=C)C(O)(CC)C(=C)O(20165)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction46',
    reactants = ['C=[C]C(C)C([O])(CC)C(=C)O(20152)'],
    products = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 5 used for cCs(-HC)CJ;CdsJ;C
Exact match found for rate rule [cCs(-HC)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction47',
    reactants = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    products = ['C=C(O)C1(CC)OC(C)C1=C(25992)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;O_rad;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction48',
    reactants = ['CHCH3(T)(95)', 'C=[C]C([O])(CC)C(=C)O(9633)'],
    products = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4293',
    isomers = [
        'C=C(O)C([O])(CC)C(=C)[CH]C(25632)',
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
    label = '4293',
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

