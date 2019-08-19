species(
    label = '[CH2]CC([O])C(OO)C1CCCC1(2014)',
    structure = SMILES('[CH2]CC([O])C(OO)C1CCCC1'),
    E0 = (-33.665,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,3000,3100,440,815,1455,1000,3615,1310,387.5,850,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.02852,0.112493,-4.90892e-05,-1.37579e-08,1.22666e-11,-3814.68,49.0221], Tmin=(100,'K'), Tmax=(1054.24,'K')), NASAPolynomial(coeffs=[23.6472,0.0541325,-2.16264e-05,4.01879e-09,-2.82733e-13,-11398.9,-86.4987], Tmin=(1054.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-33.665,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(656.843,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(CC(C)OJ) + radical(RCCJ)"""),
)

species(
    label = 'C2H4(45)',
    structure = SMILES('C=C'),
    E0 = (42.0619,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2334.71,'J/mol'), sigma=(3.971,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.9592,-0.00757051,5.7099e-05,-6.91588e-08,2.69884e-11,5089.78,4.0973], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.99183,0.0104834,-3.71721e-06,5.94628e-10,-3.5363e-14,4268.66,-0.269082], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(42.0619,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""C2H4""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'O=CC(OO)C1CCCC1(1993)',
    structure = SMILES('O=CC(OO)C1CCCC1'),
    E0 = (-329.452,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3615,1310,387.5,850,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4681.46,'J/mol'), sigma=(7.80337,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=731.23 K, Pc=22.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0321018,0.063013,3.67252e-05,-9.13792e-08,3.86252e-11,-39459.5,37.6681], Tmin=(100,'K'), Tmax=(989.421,'K')), NASAPolynomial(coeffs=[20.4293,0.0387808,-1.48147e-05,2.82869e-09,-2.07746e-13,-46345.9,-74.9226], Tmin=(989.421,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-329.452,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(515.497,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-OdCsH) + ring(Cyclopentane)"""),
)

species(
    label = 'H(19)',
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
    label = '[CH2]CC(=O)C(OO)C1CCCC1(2040)',
    structure = SMILES('[CH2]CC(=O)C(OO)C1CCCC1'),
    E0 = (-198.15,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (171.214,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.83132,0.10885,-4.93763e-05,-1.04841e-08,1.05749e-11,-23605.5,47.9178], Tmin=(100,'K'), Tmax=(1070.14,'K')), NASAPolynomial(coeffs=[23.1682,0.0520694,-2.11773e-05,3.96349e-09,-2.79441e-13,-31055.4,-84.1934], Tmin=(1070.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-198.15,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(631.9,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclopentane) + radical(CJCC=O)"""),
)

species(
    label = 'C=CC([O])C(OO)C1CCCC1(2041)',
    structure = SMILES('C=CC([O])C(OO)C1CCCC1'),
    E0 = (-111.465,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3615,1310,387.5,850,1000,3010,987.5,1337.5,450,1655,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (171.214,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.66148,0.100251,-1.51868e-05,-5.20248e-08,2.66406e-11,-13181.1,47.2985], Tmin=(100,'K'), Tmax=(1010.15,'K')), NASAPolynomial(coeffs=[25.0387,0.0486806,-1.90283e-05,3.58553e-09,-2.57884e-13,-21338.4,-95.4606], Tmin=(1010.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-111.465,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(636.057,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(CC(C)OJ)"""),
)

species(
    label = 'OO[CH]C1CCCC1(1997)',
    structure = SMILES('OO[CH]C1CCCC1'),
    E0 = (-24.2356,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1310,387.5,850,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.15,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.743437,0.0552412,1.5949e-05,-5.51478e-08,2.35154e-11,-2783.25,28.8778], Tmin=(100,'K'), Tmax=(1007,'K')), NASAPolynomial(coeffs=[14.777,0.037646,-1.46669e-05,2.74074e-09,-1.956e-13,-7543.87,-48.5281], Tmin=(1007,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-24.2356,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cyclopentane) + radical(CCsJOOH)"""),
)

species(
    label = 'CH2CH2CHO(2042)',
    structure = SMILES('[CH2]CC=O'),
    E0 = (11.2619,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.221237,'amu*angstrom^2'), symmetry=1, barrier=(5.08666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.221178,'amu*angstrom^2'), symmetry=1, barrier=(5.08532,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0712,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.76345,0.0293577,-2.47892e-05,1.49239e-08,-4.38497e-12,1397.08,14.4322], Tmin=(100,'K'), Tmax=(767.858,'K')), NASAPolynomial(coeffs=[4.24224,0.0216537,-9.73869e-06,1.85596e-09,-1.3004e-13,1169.99,7.6886], Tmin=(767.858,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(11.2619,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""CH2CH2CHO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C2H4(T)(169)',
    structure = SMILES('[CH2][CH2]'),
    E0 = (318.146,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,180,1436.54,1437.15,2688.96,2689.16],'cm^-1')),
        HinderedRotor(inertia=(0.0257549,'amu*angstrom^2'), symmetry=1, barrier=(17.2441,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.40736,0.0100312,6.40927e-06,-1.41291e-08,5.92671e-12,38288.2,6.11703], Tmin=(100,'K'), Tmax=(954.26,'K')), NASAPolynomial(coeffs=[5.52249,0.00856173,-2.90743e-06,5.02353e-10,-3.44572e-14,37547.8,-5.75276], Tmin=(954.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.146,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""C2H4(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[O][CH]C(OO)C1CCCC1(2032)',
    structure = SMILES('[O][CH]C(OO)C1CCCC1'),
    E0 = (-3.06575,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3615,1310,387.5,850,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.720628,0.0925877,-6.32804e-05,2.17495e-08,-3.03318e-12,-189.696,39.0194], Tmin=(100,'K'), Tmax=(1656.29,'K')), NASAPolynomial(coeffs=[20.6219,0.0410454,-1.66022e-05,2.96149e-09,-1.97353e-13,-7259.65,-74.715], Tmin=(1656.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-3.06575,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(515.497,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cyclopentane) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C[C](O)C(OO)C1CCCC1(2043)',
    structure = SMILES('[CH2]C[C](O)C(OO)C1CCCC1'),
    E0 = (-87.398,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.53011,0.127118,-9.7303e-05,3.82587e-08,-6.06928e-12,-10262.9,50.4562], Tmin=(100,'K'), Tmax=(1490.09,'K')), NASAPolynomial(coeffs=[26.2464,0.0498699,-1.9541e-05,3.46791e-09,-2.32231e-13,-18838.8,-99.851], Tmin=(1490.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-87.398,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(652.686,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(C2CsJOH) + radical(RCCJ)"""),
)

species(
    label = 'C[CH]C([O])C(OO)C1CCCC1(2044)',
    structure = SMILES('C[CH]C([O])C(OO)C1CCCC1'),
    E0 = (-39.0094,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.84436,0.106271,-2.9277e-05,-3.60754e-08,2.06074e-11,-4461.93,49.0687], Tmin=(100,'K'), Tmax=(1023.08,'K')), NASAPolynomial(coeffs=[24.0242,0.0527768,-2.07023e-05,3.85707e-09,-2.73754e-13,-12248.6,-88.5088], Tmin=(1023.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-39.0094,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(656.843,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(CC(C)OJ) + radical(CCJCO)"""),
)

species(
    label = '[CH2]CC(O)[C](OO)C1CCCC1(2045)',
    structure = SMILES('[CH2]CC(O)[C](OO)C1CCCC1'),
    E0 = (-77.1267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.30382,0.121494,-8.65773e-05,3.13616e-08,-4.59364e-12,-9035.11,50.5252], Tmin=(100,'K'), Tmax=(1597.04,'K')), NASAPolynomial(coeffs=[26.065,0.0504409,-1.98414e-05,3.50349e-09,-2.32742e-13,-18096.3,-99.6187], Tmin=(1597.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-77.1267,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(652.686,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(C2CsJOOH) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C(O)C(OO)C1CCCC1(2046)',
    structure = SMILES('[CH2][CH]C(O)C(OO)C1CCCC1'),
    E0 = (-64.1238,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.98003,0.111111,-4.59894e-05,-1.83941e-08,1.44311e-11,-7479.33,51.5344], Tmin=(100,'K'), Tmax=(1034.51,'K')), NASAPolynomial(coeffs=[23.9688,0.0521336,-2.04395e-05,3.78326e-09,-2.66596e-13,-15061.2,-85.2301], Tmin=(1034.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-64.1238,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(652.686,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(CCJCO) + radical(RCCJ)"""),
)

species(
    label = 'CC[C]([O])C(OO)C1CCCC1(2047)',
    structure = SMILES('CC[C]([O])C(OO)C1CCCC1'),
    E0 = (-62.2836,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.33594,0.121728,-7.93357e-05,1.99022e-08,5.54137e-15,-7248.18,47.7716], Tmin=(100,'K'), Tmax=(1157.9,'K')), NASAPolynomial(coeffs=[23.8832,0.0542325,-2.1796e-05,3.98707e-09,-2.74805e-13,-14867.1,-89.2452], Tmin=(1157.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-62.2836,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(656.843,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(C2CsJOH) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]CC(O)C(OO)[C]1CCCC1(2048)',
    structure = SMILES('[CH2]CC(O)C(OO)[C]1CCCC1'),
    E0 = (-73.115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.60096,0.103137,-3.20239e-05,-2.63924e-08,1.56691e-11,-8574.86,52.5246], Tmin=(100,'K'), Tmax=(1053.01,'K')), NASAPolynomial(coeffs=[21.5963,0.0556869,-2.23624e-05,4.16692e-09,-2.93521e-13,-15714.9,-71.2928], Tmin=(1053.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-73.115,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(652.686,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(C2CJCOOH) + radical(RCCJ)"""),
)

species(
    label = 'CCC([O])[C](OO)C1CCCC1(2049)',
    structure = SMILES('CCC([O])[C](OO)C1CCCC1'),
    E0 = (-52.0123,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.03797,0.115379,-6.65759e-05,1.09866e-08,2.11866e-12,-6023.8,47.5748], Tmin=(100,'K'), Tmax=(1168.9,'K')), NASAPolynomial(coeffs=[22.1981,0.0570316,-2.32557e-05,4.276e-09,-2.95102e-13,-13369.6,-80.3182], Tmin=(1168.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-52.0123,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(656.843,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(C2CsJOOH) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]CC(O)C(O[O])C1CCCC1(2050)',
    structure = SMILES('[CH2]CC(O)C(O[O])C1CCCC1'),
    E0 = (-112.021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.93017,0.112936,-5.98877e-05,9.02967e-10,6.72969e-12,-13244.6,49.6195], Tmin=(100,'K'), Tmax=(1070.74,'K')), NASAPolynomial(coeffs=[21.9962,0.0544918,-2.13561e-05,3.89813e-09,-2.70337e-13,-20142,-75.7289], Tmin=(1070.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-112.021,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(656.843,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(RCCJ) + radical(ROOJ)"""),
)

species(
    label = '[CH2]CC(O)C(OO)C1[CH]CCC1(2051)',
    structure = SMILES('[CH2]CC(O)C(OO)C1[CH]CCC1'),
    E0 = (-69.4837,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.07935,0.116064,-6.84984e-05,1.11634e-08,2.51547e-12,-8123.36,51.9609], Tmin=(100,'K'), Tmax=(1144.21,'K')), NASAPolynomial(coeffs=[22.8824,0.0548119,-2.22998e-05,4.11384e-09,-2.85169e-13,-15538.4,-79.2684], Tmin=(1144.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-69.4837,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(652.686,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(RCCJ) + radical(Cs_S)"""),
)

species(
    label = 'CCC([O])C(OO)[C]1CCCC1(2052)',
    structure = SMILES('CCC([O])C(OO)[C]1CCCC1'),
    E0 = (-48.0006,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.46298,0.0982758,-1.52648e-05,-4.40891e-08,2.18312e-11,-5557.57,50.0502], Tmin=(100,'K'), Tmax=(1038.36,'K')), NASAPolynomial(coeffs=[21.6069,0.0564021,-2.2665e-05,4.24984e-09,-3.01416e-13,-12882.1,-74.3172], Tmin=(1038.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-48.0006,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(656.843,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(C2CJCOOH) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]CC(O)C(OO)C1C[CH]CC1(2053)',
    structure = SMILES('[CH2]CC(O)C(OO)C1C[CH]CC1'),
    E0 = (-78.1447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.11446,0.118637,-8.28087e-05,2.95195e-08,-4.27434e-12,-9165.52,52.1125], Tmin=(100,'K'), Tmax=(1606.62,'K')), NASAPolynomial(coeffs=[24.7858,0.0516634,-2.02799e-05,3.57329e-09,-2.36956e-13,-17809.3,-90.4203], Tmin=(1606.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-78.1447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(652.686,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(cyclopentane) + radical(RCCJ)"""),
)

species(
    label = 'CCC([O])C(O[O])C1CCCC1(2054)',
    structure = SMILES('CCC([O])C(O[O])C1CCCC1'),
    E0 = (-86.9067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.77931,0.107929,-4.26555e-05,-1.73508e-08,1.31026e-11,-10227.9,47.0985], Tmin=(100,'K'), Tmax=(1046.57,'K')), NASAPolynomial(coeffs=[21.9093,0.0553666,-2.17482e-05,4.00175e-09,-2.79923e-13,-17266,-78.1997], Tmin=(1046.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-86.9067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(661.001,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(ROOJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'CCC([O])C(OO)C1[CH]CCC1(2055)',
    structure = SMILES('CCC([O])C(OO)C1[CH]CCC1'),
    E0 = (-44.3693,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.90436,0.110815,-5.06189e-05,-7.63794e-09,8.99696e-12,-5107.73,49.3507], Tmin=(100,'K'), Tmax=(1093.8,'K')), NASAPolynomial(coeffs=[22.4079,0.0563003,-2.30276e-05,4.29371e-09,-3.00893e-13,-12483.8,-79.527], Tmin=(1093.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-44.3693,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(656.843,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(CC(C)OJ) + radical(Cs_S)"""),
)

species(
    label = 'CCC([O])C(OO)C1C[CH]CC1(2056)',
    structure = SMILES('CCC([O])C(OO)C1C[CH]CC1'),
    E0 = (-53.0303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.91092,0.113193,-6.48824e-05,1.14777e-08,1.58403e-12,-6151.36,49.3901], Tmin=(100,'K'), Tmax=(1186.67,'K')), NASAPolynomial(coeffs=[21.2999,0.0576512,-2.33647e-05,4.27103e-09,-2.93305e-13,-13258.1,-73.2946], Tmin=(1186.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-53.0303,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(656.843,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(cyclopentane) + radical(CC(C)OJ)"""),
)

species(
    label = 'CCC(=O)C(OO)C1CCCC1(2057)',
    structure = SMILES('CCC(=O)C(OO)C1CCCC1'),
    E0 = (-409.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.66024,0.103444,-2.88919e-05,-3.05584e-08,1.71045e-11,-49058.5,46.6697], Tmin=(100,'K'), Tmax=(1055.95,'K')), NASAPolynomial(coeffs=[22.2423,0.056033,-2.28141e-05,4.28715e-09,-3.03516e-13,-56511.2,-81.3342], Tmin=(1055.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-409.74,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(656.843,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclopentane)"""),
)

species(
    label = 'C=CC(O)C(OO)C1CCCC1(2058)',
    structure = SMILES('C=CC(O)C(OO)C1CCCC1'),
    E0 = (-341.826,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.82682,0.104276,-2.12052e-05,-4.68407e-08,2.4971e-11,-40881.5,48.0837], Tmin=(100,'K'), Tmax=(1011.34,'K')), NASAPolynomial(coeffs=[25.1188,0.0504474,-1.95989e-05,3.66937e-09,-2.62564e-13,-49029.2,-95.5532], Tmin=(1011.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-341.826,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(656.843,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopentane)"""),
)

species(
    label = 'OH(21)',
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
    label = '[CH2]CC1OOC1C1CCCC1(2059)',
    structure = SMILES('[CH2]CC1OOC1C1CCCC1'),
    E0 = (65.1566,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (155.214,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.572559,0.078106,1.94625e-05,-7.32832e-08,3.12528e-11,8020.79,37.0891], Tmin=(100,'K'), Tmax=(1013.23,'K')), NASAPolynomial(coeffs=[18.8106,0.0535168,-2.10151e-05,3.93374e-09,-2.80462e-13,1427.2,-69.8324], Tmin=(1013.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(65.1566,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(12dioxetane) + ring(Cyclopentane) + radical(RCCJ)"""),
)

species(
    label = '[O]C1CCOC1C1CCCC1(2060)',
    structure = SMILES('[O]C1CCOC1C1CCCC1'),
    E0 = (-204.847,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (155.214,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.10113,0.0521664,0.000111277,-1.85266e-07,7.68839e-11,-24466.1,33.3413], Tmin=(100,'K'), Tmax=(929.496,'K')), NASAPolynomial(coeffs=[22.3368,0.0422443,-1.112e-05,1.79264e-09,-1.2849e-13,-32304.7,-92.2374], Tmin=(929.496,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-204.847,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(627.743,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Tetrahydrofuran) + ring(Cyclopentane) + radical(CC(C)OJ)"""),
)

species(
    label = 'OOC(C1CCCC1)C1CCO1(2016)',
    structure = SMILES('OOC(C1CCCC1)C1CCO1'),
    E0 = (-285.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.13937,0.085514,3.68118e-05,-1.12058e-07,5.08983e-11,-34096.3,41.6403], Tmin=(100,'K'), Tmax=(941.509,'K')), NASAPolynomial(coeffs=[23.5728,0.0494032,-1.53936e-05,2.61049e-09,-1.82293e-13,-41802.5,-92.3048], Tmin=(941.509,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-285.24,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(665.158,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Oxetane) + ring(Cyclopentane)"""),
)

species(
    label = '[O]C(CCO)C([O])C1CCCC1(2061)',
    structure = SMILES('[O]C(CCO)C([O])C1CCCC1'),
    E0 = (-239.883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.89707,0.109982,-4.51837e-05,-1.78004e-08,1.40307e-11,-28621.9,48.5103], Tmin=(100,'K'), Tmax=(1032.45,'K')), NASAPolynomial(coeffs=[23.0771,0.0531484,-2.06145e-05,3.78782e-09,-2.65627e-13,-35906.6,-83.0778], Tmin=(1032.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-239.883,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(661.001,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cyclopentane) + radical(CC(C)OJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'O(20)',
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
    label = '[CH2]C[CH]C(OO)C1CCCC1(2062)',
    structure = SMILES('[CH2]C[CH]C(OO)C1CCCC1'),
    E0 = (108.86,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.42706,0.100792,-3.57883e-05,-1.93463e-08,1.27391e-11,13304,47.0032], Tmin=(100,'K'), Tmax=(1065.29,'K')), NASAPolynomial(coeffs=[20.5764,0.0545748,-2.19676e-05,4.08095e-09,-2.86277e-13,6550.45,-70.237], Tmin=(1065.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(108.86,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(631.9,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(CCJCOOH) + radical(RCCJ)"""),
)

species(
    label = 'CH2(T)(33)',
    structure = SMILES('[CH2]'),
    E0 = (381.08,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([971.045,2816.03,3444.23],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.71758,0.00127391,2.17347e-06,-3.48858e-09,1.65209e-12,45872.4,1.75298], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.14632,0.00303671,-9.96474e-07,1.50484e-10,-8.57336e-15,46041.3,4.72342], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(381.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH2]C([O])C(OO)C1CCCC1(2063)',
    structure = SMILES('[CH2]C([O])C(OO)C1CCCC1'),
    E0 = (-3.54207,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3100,440,815,1455,1000,3615,1310,387.5,850,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.46508,0.100943,-4.04415e-05,-1.99494e-08,1.46972e-11,-212.496,43.5997], Tmin=(100,'K'), Tmax=(1025.31,'K')), NASAPolynomial(coeffs=[22.9773,0.0455202,-1.77796e-05,3.30093e-09,-2.33726e-13,-7323.72,-85.1676], Tmin=(1025.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-3.54207,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(CC(C)OJ) + radical(CJCO)"""),
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
    E0 = (-33.665,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (85.0835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (110.368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (18.9504,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (44.7809,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (60.4542,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (132.005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (114.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (49.9941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (41.6156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (84.2741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (67.169,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (63.1455,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (25.1176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-7.93345,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (46.2494,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (46.0402,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (22.8399,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (9.01175,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (21.5638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (315.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (29.7351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (29.7351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (97.0171,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (28.2582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-25.3807,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (34.5342,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (351.864,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (377.538,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]CC([O])C(OO)C1CCCC1(2014)'],
    products = ['C2H4(45)', 'O=CC(OO)C1CCCC1(1993)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(19)', '[CH2]CC(=O)C(OO)C1CCCC1(2040)'],
    products = ['[CH2]CC([O])C(OO)C1CCCC1(2014)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.0366254,'m^3/(mol*s)'), n=1.743, Ea=(71.4418,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsCs_O;YJ] for rate rule [CO-CsCs_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(19)', 'C=CC([O])C(OO)C1CCCC1(2041)'],
    products = ['[CH2]CC([O])C(OO)C1CCCC1(2014)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.01e+08,'cm^3/(mol*s)'), n=1.6, Ea=(10.0416,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 10 used for Cds-CsH_Cds-HH;HJ
Exact match found for rate rule [Cds-CsH_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['OO[CH]C1CCCC1(1997)', 'CH2CH2CHO(2042)'],
    products = ['[CH2]CC([O])C(OO)C1CCCC1(2014)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(8.09579,'m^3/(mol*s)'), n=1.29734, Ea=(31.9242,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;CsJ] for rate rule [CO-CsH_O;CsJ-OsCsH]
Euclidian distance = 2.2360679775
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C2H4(T)(169)', 'O=CC(OO)C1CCCC1(1993)'],
    products = ['[CH2]CC([O])C(OO)C1CCCC1(2014)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.0403742,'m^3/(mol*s)'), n=2.2105, Ea=(56.0866,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsH_O;YJ] for rate rule [CO-CsH_O;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C2H4(45)', '[O][CH]C(OO)C1CCCC1(2032)'],
    products = ['[CH2]CC([O])C(OO)C1CCCC1(2014)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.00310016,'m^3/(mol*s)'), n=2.49247, Ea=(21.458,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-HH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]CC([O])C(OO)C1CCCC1(2014)'],
    products = ['[CH2]C[C](O)C(OO)C1CCCC1(2043)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.56178e+08,'s^-1'), n=1.25272, Ea=(165.67,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_Cs2] for rate rule [R2H_S;O_rad_out;Cs_H_out_Cs2]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]CC([O])C(OO)C1CCCC1(2014)'],
    products = ['C[CH]C([O])C(OO)C1CCCC1(2044)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6.48e+07,'s^-1'), n=1.57, Ea=(147.695,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 106 used for R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]CC([O])C(OO)C1CCCC1(2014)'],
    products = ['[CH2]CC(O)[C](OO)C1CCCC1(2045)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.0250998,'s^-1'), n=4.2935, Ea=(83.6591,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_OOH/Cs]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]CC([O])C(OO)C1CCCC1(2014)'],
    products = ['[CH2][CH]C(O)C(OO)C1CCCC1(2046)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(223829,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]CC([O])C(OO)C1CCCC1(2014)'],
    products = ['CC[C]([O])C(OO)C1CCCC1(2047)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]CC([O])C(OO)C1CCCC1(2014)'],
    products = ['[CH2]CC(O)C(OO)[C]1CCCC1(2048)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1e+11,'s^-1'), n=0, Ea=(100.834,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4H_SSS;O_rad_out;Cs_H_out_Cs2] for rate rule [R4H_SSS;O_rad_out;Cs_H_out_Cs2_cy5]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]CC([O])C(OO)C1CCCC1(2014)'],
    products = ['CCC([O])[C](OO)C1CCCC1(2049)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.43385e+06,'s^-1'), n=1.65912, Ea=(96.8105,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;Cs_H_out_OOH/Cs] + [R4H_SSS;C_rad_out_2H;Cs_H_out] for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_OOH/Cs]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]CC([O])C(OO)C1CCCC1(2014)'],
    products = ['[CH2]CC(O)C(O[O])C1CCCC1(2050)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.37227e+06,'s^-1'), n=1.56745, Ea=(58.7826,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;O_rad_out;XH_out] for rate rule [R5H_SSSS;O_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]CC([O])C(OO)C1CCCC1(2014)'],
    products = ['[CH2]CC(O)C(OO)C1[CH]CCC1(2051)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.6e+11,'s^-1'), n=0, Ea=(25.7316,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(1000,'K'), comment="""Estimated using template [R5H_CCC;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5H_CCC;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]CC([O])C(OO)C1CCCC1(2014)'],
    products = ['CCC([O])C(OO)[C]1CCCC1(2052)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.25e+10,'s^-1'), n=0, Ea=(79.9144,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5H_CCC;Y_rad_out;Cs_H_out_Cs2] for rate rule [R5H_CCC;C_rad_out_2H;Cs_H_out_Cs2_cy5]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]CC([O])C(OO)C1CCCC1(2014)'],
    products = ['[CH2]CC(O)C(OO)C1C[CH]CC1(2053)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(6.24e+09,'s^-1'), n=0, Ea=(79.7052,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6H_SSSSS;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R6H_SSSSS;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]CC([O])C(OO)C1CCCC1(2014)'],
    products = ['CCC([O])C(O[O])C1CCCC1(2054)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(135.127,'s^-1'), n=2.18479, Ea=(56.5049,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_SSSSS;C_rad_out_2H;XH_out] for rate rule [R6H_SSSSS;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]CC([O])C(OO)C1CCCC1(2014)'],
    products = ['CCC([O])C(OO)C1[CH]CCC1(2055)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(51600,'s^-1'), n=1.67, Ea=(42.6768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_SSSSS;C_rad_out_2H;Cs_H_out_H/NonDeC] for rate rule [R6H_SSSSS;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]CC([O])C(OO)C1CCCC1(2014)'],
    products = ['CCC([O])C(OO)C1C[CH]CC1(2056)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2124,'s^-1'), n=1.81, Ea=(55.2288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7H;C_rad_out_2H;Cs_H_out_H/NonDeC] for rate rule [R7H;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C2H4(T)(169)', '[O][CH]C(OO)C1CCCC1(2032)'],
    products = ['[CH2]CC([O])C(OO)C1CCCC1(2014)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.47003e+07,'m^3/(mol*s)'), n=0.0284742, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]CC([O])C(OO)C1CCCC1(2014)'],
    products = ['CCC(=O)C(OO)C1CCCC1(2057)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]CC([O])C(OO)C1CCCC1(2014)'],
    products = ['C=CC(O)C(OO)C1CCCC1(2058)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]CC([O])C(OO)C1CCCC1(2014)'],
    products = ['OH(21)', '[CH2]CC1OOC1C1CCCC1(2059)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.11355e+11,'s^-1'), n=0, Ea=(130.682,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3OO_SS;Y_rad_intra;OOH]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]CC([O])C(OO)C1CCCC1(2014)'],
    products = ['OH(21)', '[O]C1CCOC1C1CCCC1(2060)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(5.13e+10,'s^-1','*|/',1.41), n=0, Ea=(61.9232,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 14 used for R4OO_SSS;C_pri_rad_intra;OOH
Exact match found for rate rule [R4OO_SSS;C_pri_rad_intra;OOH]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]CC([O])C(OO)C1CCCC1(2014)'],
    products = ['OOC(C1CCCC1)C1CCO1(2016)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R4_SSS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]CC([O])C(OO)C1CCCC1(2014)'],
    products = ['[O]C(CCO)C([O])C1CCCC1(2061)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.12e+10,'s^-1'), n=0, Ea=(68.1992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 6 used for R4OOH_SSS;C_rad_out_2H
Exact match found for rate rule [R4OOH_SSS;C_rad_out_2H]
Euclidian distance = 0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['O(20)', '[CH2]C[CH]C(OO)C1CCCC1(2062)'],
    products = ['[CH2]CC([O])C(OO)C1CCCC1(2014)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/NonDeC;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction29',
    reactants = ['CH2(T)(33)', '[CH2]C([O])C(OO)C1CCCC1(2063)'],
    products = ['[CH2]CC([O])C(OO)C1CCCC1(2014)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.23625e+06,'m^3/(mol*s)'), n=0.36814, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

network(
    label = '322',
    isomers = [
        '[CH2]CC([O])C(OO)C1CCCC1(2014)',
    ],
    reactants = [
        ('C2H4(45)', 'O=CC(OO)C1CCCC1(1993)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '322',
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

