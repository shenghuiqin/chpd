species(
    label = 'S(228)(227)',
    structure = SMILES('[O]OCC=CC=CO'),
    E0 = (-80.5231,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3615,1277.5,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.869004,'amu*angstrom^2'), symmetry=1, barrier=(19.9801,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.868435,'amu*angstrom^2'), symmetry=1, barrier=(19.967,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.868703,'amu*angstrom^2'), symmetry=1, barrier=(19.9732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.866613,'amu*angstrom^2'), symmetry=1, barrier=(19.9251,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4481.69,'J/mol'), sigma=(7.05185,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=700.03 K, Pc=29 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.807315,0.0859974,-8.81813e-05,4.37494e-08,-8.14715e-12,-9494.72,32.1986], Tmin=(100,'K'), Tmax=(1494.57,'K')), NASAPolynomial(coeffs=[23.057,0.0108932,-1.52852e-06,6.74665e-11,9.47456e-16,-15373.3,-88.3245], Tmin=(1494.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-80.5231,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(ROOJ)"""),
)

species(
    label = 'O2(2)(2)',
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
    label = 'C5H7O(224)(223)',
    structure = SMILES('[CH2]C=CC=CO'),
    E0 = (-32.8413,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.31979,'amu*angstrom^2'), symmetry=1, barrier=(30.3447,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32385,'amu*angstrom^2'), symmetry=1, barrier=(30.4378,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31983,'amu*angstrom^2'), symmetry=1, barrier=(30.3455,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3794.25,'J/mol'), sigma=(6.17562,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=592.65 K, Pc=36.55 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25554,0.0396093,3.48643e-05,-9.07575e-08,4.31029e-11,-3831.66,21.0706], Tmin=(100,'K'), Tmax=(908.805,'K')), NASAPolynomial(coeffs=[21.7542,0.00426328,2.6289e-06,-6.6847e-10,4.32897e-14,-9823.72,-88.3316], Tmin=(908.805,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-32.8413,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CC=CCJ)"""),
)

species(
    label = 'OC=C[CH]C1COO1(7291)',
    structure = SMILES('OC=C[CH]C1COO1'),
    E0 = (-63.7885,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0254379,0.0731306,-4.04386e-05,-1.10217e-08,1.24571e-11,-7515.97,23.2658], Tmin=(100,'K'), Tmax=(962.046,'K')), NASAPolynomial(coeffs=[21.5304,0.0172267,-5.52217e-06,9.88194e-10,-7.23809e-14,-13204.5,-87.7108], Tmin=(962.046,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-63.7885,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(12dioxetane) + radical(C=CCJCO)"""),
)

species(
    label = 'O[CH]C1C=CCOO1(7292)',
    structure = SMILES('O[CH]C1C=CCOO1'),
    E0 = (-53.3871,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.856462,0.0530259,1.01326e-06,-4.39645e-08,2.15783e-11,-6293.4,25.0653], Tmin=(100,'K'), Tmax=(980.297,'K')), NASAPolynomial(coeffs=[18.178,0.0204444,-7.42684e-06,1.41952e-09,-1.05981e-13,-11520,-67.4932], Tmin=(980.297,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-53.3871,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(36dihydro12dioxin) + radical(CCsJOH)"""),
)

species(
    label = 'OC=CC=C[CH]OO(7293)',
    structure = SMILES('OC=CC=C[CH]OO'),
    E0 = (-115.233,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.17086,0.0704856,-1.5825e-05,-5.14879e-08,3.08282e-11,-13689.5,29.745], Tmin=(100,'K'), Tmax=(928.828,'K')), NASAPolynomial(coeffs=[27.5839,0.0055353,9.28867e-07,-2.52719e-10,1.09836e-14,-21199.5,-114.78], Tmin=(928.828,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-115.233,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CCJO)"""),
)

species(
    label = 'OC=CC=[C]COO(7294)',
    structure = SMILES('OC=CC=[C]COO'),
    E0 = (5.31397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.84922,0.0912549,-9.82833e-05,5.10634e-08,-1.00522e-11,826.865,32.0586], Tmin=(100,'K'), Tmax=(1357.02,'K')), NASAPolynomial(coeffs=[23.9217,0.0114118,-2.481e-06,2.90899e-10,-1.54915e-14,-5267.43,-92.693], Tmin=(1357.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.31397,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(Cds_S)"""),
)

species(
    label = 'OC=C[C]=CCOO(7295)',
    structure = SMILES('OC=C[C]=CCOO'),
    E0 = (-33.5324,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,1685,370,2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1277.5,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.15833,'amu*angstrom^2'), symmetry=1, barrier=(26.6323,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15793,'amu*angstrom^2'), symmetry=1, barrier=(26.6231,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15817,'amu*angstrom^2'), symmetry=1, barrier=(26.6287,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15815,'amu*angstrom^2'), symmetry=1, barrier=(26.6281,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15795,'amu*angstrom^2'), symmetry=1, barrier=(26.6236,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.879288,0.0930879,-0.000102506,5.49961e-08,-1.11792e-11,-3845.15,31.2851], Tmin=(100,'K'), Tmax=(1322.51,'K')), NASAPolynomial(coeffs=[23.2721,0.0124673,-2.47585e-06,2.41421e-10,-1.01238e-14,-9570.92,-89.4783], Tmin=(1322.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-33.5324,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CJC=C)"""),
)

species(
    label = 'OC=[C]C=CCOO(7296)',
    structure = SMILES('OC=C=C[CH]COO'),
    E0 = (-56.1415,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.379,'amu*angstrom^2'), symmetry=1, barrier=(31.706,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.37897,'amu*angstrom^2'), symmetry=1, barrier=(31.7053,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.37944,'amu*angstrom^2'), symmetry=1, barrier=(31.7161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.37887,'amu*angstrom^2'), symmetry=1, barrier=(31.703,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.379,'amu*angstrom^2'), symmetry=1, barrier=(31.7059,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.06027,0.097078,-0.000105318,5.49186e-08,-1.09287e-11,-6557.96,29.628], Tmin=(100,'K'), Tmax=(1286.98,'K')), NASAPolynomial(coeffs=[25.1264,0.0128194,-3.76861e-06,5.83171e-10,-3.73835e-14,-13060.7,-102.391], Tmin=(1286.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-56.1415,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=CCJCO)"""),
)

species(
    label = 'O[C]=CC=CCOO(7297)',
    structure = SMILES('O[C]=CC=CCOO'),
    E0 = (7.21638,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,1685,370,2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1277.5,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.868176,'amu*angstrom^2'), symmetry=1, barrier=(19.9611,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.869011,'amu*angstrom^2'), symmetry=1, barrier=(19.9803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.868379,'amu*angstrom^2'), symmetry=1, barrier=(19.9657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.868724,'amu*angstrom^2'), symmetry=1, barrier=(19.9737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.869261,'amu*angstrom^2'), symmetry=1, barrier=(19.986,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.154029,0.0831439,-8.22473e-05,3.7466e-08,-5.60813e-12,1024.75,32.2147], Tmin=(100,'K'), Tmax=(999.731,'K')), NASAPolynomial(coeffs=[19.6017,0.0183321,-6.35799e-06,1.09966e-09,-7.50855e-14,-3636.58,-66.6473], Tmin=(999.731,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(7.21638,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CJO)"""),
)

species(
    label = '[O]C=CC=CCOO(7298)',
    structure = SMILES('O=CC=C[CH]COO'),
    E0 = (-111.416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.74109,0.0792805,-7.96442e-05,4.3692e-08,-1.0239e-11,-13289.2,25.5404], Tmin=(100,'K'), Tmax=(993.626,'K')), NASAPolynomial(coeffs=[10.425,0.0402968,-2.07944e-05,4.20764e-09,-3.04743e-13,-15213.6,-21.1172], Tmin=(993.626,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-111.416,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJCO)"""),
)

species(
    label = 'OH(5)(5)',
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
    label = '[CH]=CC=CCO[O](2610)',
    structure = SMILES('[CH]=CC=CCO[O]'),
    E0 = (375.366,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,180],'cm^-1')),
        HinderedRotor(inertia=(0.660561,'amu*angstrom^2'), symmetry=1, barrier=(15.1876,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.663518,'amu*angstrom^2'), symmetry=1, barrier=(15.2556,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.663598,'amu*angstrom^2'), symmetry=1, barrier=(15.2574,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10283,0.0626934,-6.20181e-05,3.29906e-08,-7.11315e-12,45251.3,25.8806], Tmin=(100,'K'), Tmax=(1115.62,'K')), NASAPolynomial(coeffs=[11.9293,0.0238755,-9.82536e-06,1.80133e-09,-1.23872e-13,42835.7,-27.5351], Tmin=(1115.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(375.366,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_P) + radical(ROOJ)"""),
)

species(
    label = 'H(3)(3)',
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
    label = '[O]OC[CH]C=CC=O(7129)',
    structure = SMILES('[O]OC[CH]C=CC=O'),
    E0 = (40.5891,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2782.5,750,1395,475,1775,1000,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.42599,'amu*angstrom^2'), symmetry=1, barrier=(32.7864,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.646963,'amu*angstrom^2'), symmetry=1, barrier=(14.875,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.646122,'amu*angstrom^2'), symmetry=1, barrier=(14.8556,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12569,'amu*angstrom^2'), symmetry=1, barrier=(25.8818,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.766893,0.0784694,-9.39986e-05,6.56001e-08,-1.95106e-11,4991.65,26.2082], Tmin=(100,'K'), Tmax=(800.872,'K')), NASAPolynomial(coeffs=[8.67507,0.0389727,-2.00251e-05,4.02441e-09,-2.89717e-13,3724.93,-10.1883], Tmin=(800.872,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(40.5891,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(ROOJ) + radical(C=CCJCO)"""),
)

species(
    label = '[CH2]O[O](1428)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.45042,0.0115293,-8.64543e-06,3.98163e-09,-7.73155e-13,24893.2,10.7919], Tmin=(100,'K'), Tmax=(1212.82,'K')), NASAPolynomial(coeffs=[5.07483,0.00617178,-2.01911e-06,3.39169e-10,-2.23136e-14,24499.1,2.64172], Tmin=(1212.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(206.805,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(ROOJ) + radical(CsJOOH)"""),
)

species(
    label = '[CH]=CC=CO(6201)',
    structure = SMILES('[CH]=CC=CO'),
    E0 = (132.225,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(1.18593,'amu*angstrom^2'), symmetry=1, barrier=(27.2668,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19042,'amu*angstrom^2'), symmetry=1, barrier=(27.3702,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56381,0.0362679,1.62217e-05,-6.68286e-08,3.43645e-11,16007,17.0038], Tmin=(100,'K'), Tmax=(899.999,'K')), NASAPolynomial(coeffs=[21.2117,-0.00423318,5.68525e-06,-1.2177e-09,8.19777e-14,10574,-86.2509], Tmin=(899.999,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(132.225,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[O]O[CH]C=CC=CO(7299)',
    structure = SMILES('[O]O[CH]C=CC=CO'),
    E0 = (36.772,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3025,407.5,1350,352.5,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3615,1277.5,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.09155,'amu*angstrom^2'), symmetry=1, barrier=(25.097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09215,'amu*angstrom^2'), symmetry=1, barrier=(25.1106,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09167,'amu*angstrom^2'), symmetry=1, barrier=(25.0996,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09144,'amu*angstrom^2'), symmetry=1, barrier=(25.0943,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.118038,0.0665948,-1.94631e-05,-4.36018e-08,2.76995e-11,4579.9,29.4674], Tmin=(100,'K'), Tmax=(918.581,'K')), NASAPolynomial(coeffs=[25.7713,0.00426273,1.69346e-06,-4.39458e-10,2.6602e-14,-2216.16,-103.455], Tmin=(918.581,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(36.772,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CCJO) + radical(ROOJ)"""),
)

species(
    label = '[O]OC[C]=CC=CO(7300)',
    structure = SMILES('[O]OC[C]=CC=CO'),
    E0 = (157.319,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1685,370,2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1277.5,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.834603,'amu*angstrom^2'), symmetry=1, barrier=(19.1892,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.835128,'amu*angstrom^2'), symmetry=1, barrier=(19.2012,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.835003,'amu*angstrom^2'), symmetry=1, barrier=(19.1984,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.83508,'amu*angstrom^2'), symmetry=1, barrier=(19.2001,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.427957,0.0858996,-9.72588e-05,5.35784e-08,-1.11936e-11,19090.3,31.2997], Tmin=(100,'K'), Tmax=(1273.03,'K')), NASAPolynomial(coeffs=[21.5185,0.0111086,-2.26042e-06,2.30005e-10,-1.01421e-14,13975.3,-78.0207], Tmin=(1273.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.319,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = 'CHCHOH(49)(49)',
    structure = SMILES('[CH]=CO'),
    E0 = (120.933,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3120,650,792.5,1650,458.887,459.577,459.893,460.324],'cm^-1')),
        HinderedRotor(inertia=(0.000141718,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1868.27,'J/mol'), sigma=(4.162,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.99181,0.0473233,-6.66059e-05,4.68997e-08,-1.30686e-11,15200.1,31.4259], Tmin=(298,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[8.78246,0.00524797,-1.71857e-06,2.59722e-10,-1.48227e-14,12883.6,-21.0851], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(298,'K'), Tmax=(6000,'K'), E0=(120.933,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CHCHOH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[CH]=CCO[O](2607)',
    structure = SMILES('[CH]=CCO[O]'),
    E0 = (329.273,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350],'cm^-1')),
        HinderedRotor(inertia=(0.376082,'amu*angstrom^2'), symmetry=1, barrier=(8.64688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.375142,'amu*angstrom^2'), symmetry=1, barrier=(8.62525,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.0627,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96452,0.0519072,-9.03046e-05,8.42143e-08,-2.96214e-11,39668.9,18.753], Tmin=(100,'K'), Tmax=(878.278,'K')), NASAPolynomial(coeffs=[5.14607,0.0210528,-9.65991e-06,1.78528e-09,-1.19369e-13,39741.2,7.40998], Tmin=(878.278,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.273,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = '[O]OCC=[C]C=CO(7301)',
    structure = SMILES('[O]OCC=[C]C=CO'),
    E0 = (118.472,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1685,370,2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1277.5,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.984691,'amu*angstrom^2'), symmetry=1, barrier=(22.64,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.984959,'amu*angstrom^2'), symmetry=1, barrier=(22.6462,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.984612,'amu*angstrom^2'), symmetry=1, barrier=(22.6382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.984786,'amu*angstrom^2'), symmetry=1, barrier=(22.6422,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0764494,0.0830944,-8.45924e-05,3.45939e-08,-2.06633e-12,14401.9,29.1666], Tmin=(100,'K'), Tmax=(888.726,'K')), NASAPolynomial(coeffs=[20.1199,0.0134833,-3.03363e-06,3.67284e-10,-2.04381e-14,9971.41,-70.6165], Tmin=(888.726,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(118.472,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CJC=C) + radical(ROOJ)"""),
)

species(
    label = '[O]OCC=C[C]=CO(7302)',
    structure = SMILES('[O]OC[CH]C=C=CO'),
    E0 = (95.8632,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,1437.5,1250,1305,750,350,492.5,1135,1000,3025,407.5,1350,352.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.25031,'amu*angstrom^2'), symmetry=1, barrier=(28.747,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24909,'amu*angstrom^2'), symmetry=1, barrier=(28.719,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.25081,'amu*angstrom^2'), symmetry=1, barrier=(28.7586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24947,'amu*angstrom^2'), symmetry=1, barrier=(28.7278,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.415767,0.0890943,-9.51579e-05,4.57323e-08,-7.16296e-12,11695.9,28.0685], Tmin=(100,'K'), Tmax=(967.236,'K')), NASAPolynomial(coeffs=[21.6029,0.014366,-4.59286e-06,7.65394e-10,-5.19731e-14,6672.56,-81.3741], Tmin=(967.236,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(95.8632,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=CCJCO) + radical(ROOJ)"""),
)

species(
    label = '[O]OCC=CC=[C]O(7303)',
    structure = SMILES('[O]OCC=CC=[C]O'),
    E0 = (159.221,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1685,370,2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1277.5,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.695235,'amu*angstrom^2'), symmetry=1, barrier=(15.9848,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.69578,'amu*angstrom^2'), symmetry=1, barrier=(15.9974,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.69556,'amu*angstrom^2'), symmetry=1, barrier=(15.9923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.695448,'amu*angstrom^2'), symmetry=1, barrier=(15.9897,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.14939,0.0790964,-8.54222e-05,4.49097e-08,-8.6347e-12,19293.5,31.8842], Tmin=(100,'K'), Tmax=(974.833,'K')), NASAPolynomial(coeffs=[17.654,0.0172851,-5.72181e-06,9.42965e-10,-6.19412e-14,15404.8,-54.5601], Tmin=(974.833,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(159.221,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CJO) + radical(ROOJ)"""),
)

species(
    label = 'OC=CC1[CH]COO1(7304)',
    structure = SMILES('OC=CC1[CH]COO1'),
    E0 = (-65.5556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.632245,0.048697,4.47556e-05,-1.18543e-07,5.78327e-11,-7739.21,28.0429], Tmin=(100,'K'), Tmax=(888.409,'K')), NASAPolynomial(coeffs=[26.4697,0.00192401,6.28422e-06,-1.54405e-09,1.09004e-13,-15075.1,-109], Tmin=(888.409,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-65.5556,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(12dioxolane) + radical(CCJCOOH)"""),
)

species(
    label = 'OC1[CH]C=CCOO1(7305)',
    structure = SMILES('OC1C=C[CH]COO1'),
    E0 = (-134.081,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.89044,0.0160801,0.000182109,-2.88263e-07,1.23453e-10,-15965.4,25.0312], Tmin=(100,'K'), Tmax=(911.314,'K')), NASAPolynomial(coeffs=[40.3276,-0.0204779,1.75379e-05,-3.46134e-09,2.21254e-13,-28823.2,-192.676], Tmin=(911.314,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-134.081,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptane) + radical(C=CCJCO)"""),
)

species(
    label = '[O]OC=CC=CCO(7306)',
    structure = SMILES('[O]OC=CC=CCO'),
    E0 = (-27.5333,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,492.5,1135,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.519792,'amu*angstrom^2'), symmetry=1, barrier=(11.9511,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.519768,'amu*angstrom^2'), symmetry=1, barrier=(11.9505,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.519792,'amu*angstrom^2'), symmetry=1, barrier=(11.951,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.519676,'amu*angstrom^2'), symmetry=1, barrier=(11.9484,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.306287,0.075935,-7.7203e-05,4.0931e-08,-8.57613e-12,-3173.97,31.4739], Tmin=(100,'K'), Tmax=(1164.98,'K')), NASAPolynomial(coeffs=[15.9293,0.0222917,-8.13163e-06,1.40369e-09,-9.35363e-14,-6813.99,-46.2834], Tmin=(1164.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-27.5333,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(ROOJ)"""),
)

species(
    label = 'HO2(8)(9)',
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
    label = 'C=C=CC=CO(6189)',
    structure = SMILES('C=C=CC=CO'),
    E0 = (25.7073,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,610,2055,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.24386,'amu*angstrom^2'), symmetry=1, barrier=(28.5988,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24416,'amu*angstrom^2'), symmetry=1, barrier=(28.6057,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02364,0.048058,2.85443e-06,-5.66469e-08,3.08127e-11,3215.32,19.3983], Tmin=(100,'K'), Tmax=(909.916,'K')), NASAPolynomial(coeffs=[22.142,0.00160986,2.95307e-06,-6.91045e-10,4.50216e-14,-2548.21,-91.0442], Tmin=(909.916,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(25.7073,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = '[O]OCC=CCC=O(7307)',
    structure = SMILES('[O]OCC=CCC=O'),
    E0 = (-58.7905,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62988,0.0571025,-3.56387e-05,1.00352e-08,-1.1246e-12,-6991.49,27.2298], Tmin=(100,'K'), Tmax=(1949.54,'K')), NASAPolynomial(coeffs=[15.2004,0.0292587,-1.42152e-05,2.7091e-09,-1.85125e-13,-12282.7,-47.2996], Tmin=(1949.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-58.7905,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(ROOJ)"""),
)

species(
    label = '[O]OCC1C=CC1O(7308)',
    structure = SMILES('[O]OCC1C=CC1O'),
    E0 = (3.5909,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.834598,0.0573055,-1.93554e-05,-1.87589e-08,1.19194e-11,556.744,28.7661], Tmin=(100,'K'), Tmax=(996.077,'K')), NASAPolynomial(coeffs=[16.2704,0.0227076,-8.49848e-06,1.57922e-09,-1.13469e-13,-3876.99,-52.4621], Tmin=(996.077,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(3.5909,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(ROOJ)"""),
)

species(
    label = 'O(4)(4)',
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
    label = '[O]CC=CC=CO(7309)',
    structure = SMILES('[O]CC=CC=CO'),
    E0 = (-78.3277,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2750,2850,1437.5,1250,1305,750,350,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.808808,'amu*angstrom^2'), symmetry=1, barrier=(18.5961,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.808453,'amu*angstrom^2'), symmetry=1, barrier=(18.5879,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.809117,'amu*angstrom^2'), symmetry=1, barrier=(18.6032,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.736564,0.0609652,-3.38523e-05,-8.13604e-09,1.00592e-11,-9293.29,25.9889], Tmin=(100,'K'), Tmax=(946.545,'K')), NASAPolynomial(coeffs=[17.5016,0.0164096,-4.90909e-06,8.23819e-10,-5.78704e-14,-13644.9,-60.1931], Tmin=(946.545,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-78.3277,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(CCOJ)"""),
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
    label = 'Ar(8)',
    structure = SMILES('[Ar]'),
    E0 = (-6.19426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (39.348,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1134.93,'J/mol'), sigma=(3.33,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,-745,4.3663], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,-745,4.3663], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-6.19426,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ar""", comment="""Thermo library: BurkeH2O2"""),
)

transitionState(
    label = 'TS1',
    E0 = (-41.463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (41.5199,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-53.3871,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (78.8869,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (65.0799,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-0.492151,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (63.4109,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (184.458,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (60.4777,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (403.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (252.381,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (340.126,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (251.526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (369.111,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (450.206,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (334.461,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (311.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (371.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (8.19799,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-53.3271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (136.181,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (98.1339,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (63.3436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (41.8899,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (164.677,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O2(2)(2)', 'C5H7O(224)(223)'],
    products = ['S(228)(227)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(98782.6,'m^3/(mol*s)'), n=0.576325, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [C_pri_rad;O2_birad] + [C_rad/H2/Cd;Y_rad] for rate rule [C_rad/H2/Cd;O2_birad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction2',
    reactants = ['S(228)(227)'],
    products = ['OC=C[CH]C1COO1(7291)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra;radadd_intra_O] for rate rule [R5_SS_D;doublebond_intra_HCd_pri;radadd_intra_O]
Euclidian distance = 3.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['S(228)(227)'],
    products = ['O[CH]C1C=CCOO1(7292)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(414000,'s^-1'), n=1.14, Ea=(27.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_SSSR;doublebond_intra_HNd_pri;radadd_intra] for rate rule [R7_SSSM_D;doublebond_intra_HNd_pri;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 23.8 to 27.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['S(228)(227)'],
    products = ['OC=CC=C[CH]OO(7293)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.32e+07,'s^-1'), n=1.69, Ea=(159.41,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_SS_O;O_rad_out;Cs_H_out_H/Cd] for rate rule [R3H_SS_O;O_rad_out;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction5',
    reactants = ['S(228)(227)'],
    products = ['OC=CC=[C]COO(7294)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(274,'s^-1'), n=3.09, Ea=(145.603,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 337 used for R4H_SSS_OCs;O_rad_out;Cd_H_out_doubleC
Exact match found for rate rule [R4H_SSS_OCs;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['OC=C[C]=CCOO(7295)'],
    products = ['S(228)(227)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_single;XH_out] for rate rule [R5H_DSSS;Cd_rad_out_singleDe_Cd;O_H_out]
Euclidian distance = 3.74165738677
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['OC=[C]C=CCOO(7296)'],
    products = ['S(228)(227)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.62236e+06,'s^-1'), n=1.59783, Ea=(119.552,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_Cd;XH_out] + [R6H_SMSSR;Y_rad_out;XH_out] for rate rule [R6H_SMSSR;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 3.16227766017
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O[C]=CC=CCOO(7297)'],
    products = ['S(228)(227)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(8.44313e+09,'s^-1'), n=0.985167, Ea=(177.241,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleNd;XH_out] for rate rule [R7H;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['S(228)(227)'],
    products = ['[O]C=CC=CCOO(7298)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.27e+06,'s^-1'), n=1.5, Ea=(141.001,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R8H;Y_rad_out;XH_out] for rate rule [R8H;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['OH(5)(5)', '[CH]=CC=CCO[O](2610)'],
    products = ['S(228)(227)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.44289e+07,'m^3/(mol*s)'), n=0.0225, Ea=(0.0523,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;Cd_pri_rad] + [O_pri_rad;Y_rad] for rate rule [O_pri_rad;Cd_pri_rad]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)(3)', '[O]OC[CH]C=CC=O(7129)'],
    products = ['S(228)(227)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7e+11,'cm^3/(mol*s)'), n=0.493, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2500,'K'), comment="""From training reaction 25 used for H_rad;O_rad/OneDe
Exact match found for rate rule [H_rad;O_rad/OneDe]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]O[O](1428)', '[CH]=CC=CO(6201)'],
    products = ['S(228)(227)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(6.41031e+09,'m^3/(mol*s)'), n=-0.9855, Ea=(1.09644,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_rad;C_pri_rad] for rate rule [Cd_pri_rad;C_rad/H2/O]
Euclidian distance = 1.41421356237
family: R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(3)(3)', '[O]O[CH]C=CC=CO(7299)'],
    products = ['S(228)(227)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.59671e+07,'m^3/(mol*s)'), n=0.0113737, Ea=(2.96199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H/OneDe;H_rad] for rate rule [C_rad/H/OneDeO;H_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(3)(3)', '[O]OC[C]=CC=CO(7300)'],
    products = ['S(228)(227)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CHCHOH(49)(49)', '[CH]=CCO[O](2607)'],
    products = ['S(228)(227)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.23e+13,'cm^3/(mol*s)','+|-',1.2e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 89 used for Cd_pri_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;Cd_pri_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)(3)', '[O]OCC=[C]C=CO(7301)'],
    products = ['S(228)(227)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [Cd_rad/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(3)(3)', '[O]OCC=C[C]=CO(7302)'],
    products = ['S(228)(227)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [Cd_rad/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(3)(3)', '[O]OCC=CC=[C]O(7303)'],
    products = ['S(228)(227)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['S(228)(227)'],
    products = ['OC=CC1[CH]COO1(7304)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.81589e+09,'s^-1'), n=0.403324, Ea=(88.7211,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5_linear;doublebond_intra_pri_HCd;radadd_intra] + [R5_SS_D;doublebond_intra_pri;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri_HCd;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['S(228)(227)'],
    products = ['OC1[CH]C=CCOO1(7305)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(110000,'s^-1'), n=1.18, Ea=(27.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;doublebond_intra_pri_HNd;radadd_intra] for rate rule [R7_linear;doublebond_intra_pri_HNd_O;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]OC=CC=CCO(7306)'],
    products = ['S(228)(227)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(6.8582e+08,'s^-1'), n=1.23767, Ea=(163.714,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [1_3_pentadiene;CH_end;unsaturated_end]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction22',
    reactants = ['S(228)(227)'],
    products = ['HO2(8)(9)', 'C=C=CC=CO(6189)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.63e+09,'s^-1'), n=1.11, Ea=(178.657,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 15 used for R2OO_0H_2H
Exact match found for rate rule [R2OO_0H_2H]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction23',
    reactants = ['S(228)(227)'],
    products = ['[O]OCC=CCC=O(7307)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_H;R_O_H]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction24',
    reactants = ['S(228)(227)'],
    products = ['[O]OCC1C=CC1O(7308)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1,3-butadiene_backbone;C=C_1;C=C_2] for rate rule [1,3-butadiene_backbone;CdH(C)_1;C=C_2]
Euclidian distance = 1.0
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction25',
    reactants = ['O(4)(4)', '[O]CC=CC=CO(7309)'],
    products = ['S(228)(227)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

network(
    label = '600',
    isomers = [
        'S(228)(227)',
    ],
    reactants = [
        ('O2(2)(2)', 'C5H7O(224)(223)'),
    ],
    bathGas = {
        'Ne': 0.333333,
        'N2': 0.333333,
        'Ar(8)': 0.333333,
    },
)

pressureDependence(
    label = '600',
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

