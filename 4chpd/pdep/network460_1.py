species(
    label = '[O]CC([O])C(OO)C1CCCC1(2207)',
    structure = SMILES('[O]CC([O])C(OO)C1CCCC1'),
    E0 = (-149.253,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,3615,1310,387.5,850,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.80919,0.113895,-8.02138e-05,2.8415e-08,-4.0754e-12,-17730.7,47.7912], Tmin=(100,'K'), Tmax=(1619.49,'K')), NASAPolynomial(coeffs=[24.5847,0.0487044,-1.9833e-05,3.55906e-09,-2.38397e-13,-26279.6,-92.2687], Tmin=(1619.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-149.253,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(611.114,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cyclopentane) + radical(CCOJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'CH2O(35)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.79372,-0.00990833,3.7322e-05,-3.79285e-08,1.31773e-11,-14379.2,0.602798], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.16953,0.00619321,-2.25056e-06,3.65976e-10,-2.20149e-14,-14548.7,6.04208], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-119.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH2O""", comment="""Thermo library: FFCM1(-)"""),
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
    label = '[O]CC(=O)C(OO)C1CCCC1(3889)',
    structure = SMILES('[O]CC(=O)C(OO)C1CCCC1'),
    E0 = (-285.843,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3615,1310,387.5,850,1000,375,552.5,462.5,1710,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (173.186,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.35467,0.0968795,-2.79713e-05,-3.28143e-08,1.90286e-11,-34168,46.9854], Tmin=(100,'K'), Tmax=(1023.72,'K')), NASAPolynomial(coeffs=[23.3065,0.0452844,-1.79623e-05,3.3818e-09,-2.41908e-13,-41562.8,-84.0248], Tmin=(1023.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-285.843,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsCs) + ring(Cyclopentane) + radical(C=OCOJ)"""),
)

species(
    label = '[O]C(C=O)C(OO)C1CCCC1(3890)',
    structure = SMILES('[O]C(C=O)C(OO)C1CCCC1'),
    E0 = (-281.969,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3615,1310,387.5,850,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (173.186,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.24284,0.0941016,-2.11751e-05,-3.89535e-08,2.1003e-11,-33705.8,47.6285], Tmin=(100,'K'), Tmax=(1020.29,'K')), NASAPolynomial(coeffs=[22.8312,0.0458832,-1.81528e-05,3.41644e-09,-2.44504e-13,-41021.1,-80.7733], Tmin=(1020.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-281.969,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-OdCsH) + ring(Cyclopentane) + radical(C=OCOJ)"""),
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
    label = 'OCH2CHO(3891)',
    structure = SMILES('[O]CC=O'),
    E0 = (-90.1149,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2750,2850,1437.5,1250,1305,750,350,584.193],'cm^-1')),
        HinderedRotor(inertia=(0.069063,'amu*angstrom^2'), symmetry=1, barrier=(16.7258,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (59.044,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.04366,0.0136816,1.69893e-05,-3.2265e-08,1.33002e-11,-10797.3,13.358], Tmin=(100,'K'), Tmax=(973.595,'K')), NASAPolynomial(coeffs=[8.6869,0.00860882,-3.10055e-06,5.99612e-10,-4.54054e-14,-12754.6,-18.1249], Tmin=(973.595,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-90.1149,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""OCH2CHO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2][O](354)',
    structure = SMILES('[CH2][O]'),
    E0 = (192.903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88409,-0.00363885,3.28543e-05,-4.13611e-08,1.59631e-11,23210.8,7.47983], Tmin=(100,'K'), Tmax=(933.06,'K')), NASAPolynomial(coeffs=[6.69335,0.000289989,8.61416e-07,-1.56351e-10,7.33778e-15,21991.3,-9.6043], Tmin=(933.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(H3COJ) + radical(CsJOH)"""),
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
    label = '[O]C[C](O)C(OO)C1CCCC1(3892)',
    structure = SMILES('[O]C[C](O)C(OO)C1CCCC1'),
    E0 = (-202.986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.31263,0.117492,-9.33515e-05,3.98878e-08,-7.18046e-12,-24223.3,45.5955], Tmin=(100,'K'), Tmax=(1276.76,'K')), NASAPolynomial(coeffs=[16.8482,0.0605953,-2.65075e-05,4.9851e-09,-3.46274e-13,-28860.7,-46.4574], Tmin=(1276.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-202.986,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(606.956,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cyclopentane) + radical(C2CsJOH) + radical(CCOJ)"""),
)

species(
    label = '[O]C([CH]O)C(OO)C1CCCC1(3893)',
    structure = SMILES('[O]C([CH]O)C(OO)C1CCCC1'),
    E0 = (-194.661,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.5454,0.126771,-0.000104798,4.42823e-08,-7.46817e-12,-23162.3,49.1127], Tmin=(100,'K'), Tmax=(1421.17,'K')), NASAPolynomial(coeffs=[27.1132,0.0432964,-1.66956e-05,2.95464e-09,-1.98349e-13,-31592.5,-104.398], Tmin=(1421.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-194.661,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(606.956,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cyclopentane) + radical(CCsJOH) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]CC(O)[C](OO)C1CCCC1(3894)',
    structure = SMILES('[O]CC(O)[C](OO)C1CCCC1'),
    E0 = (-192.715,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.974523,0.110788,-7.97777e-05,3.03092e-08,-4.88759e-12,-23001.1,45.2447], Tmin=(100,'K'), Tmax=(1389.27,'K')), NASAPolynomial(coeffs=[15.9165,0.0621553,-2.72684e-05,5.11151e-09,-3.53244e-13,-27694.3,-41.7981], Tmin=(1389.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-192.715,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(606.956,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cyclopentane) + radical(C2CsJOOH) + radical(CCOJ)"""),
)

species(
    label = '[O][CH]C(O)C(OO)C1CCCC1(3895)',
    structure = SMILES('[O][CH]C(O)C(OO)C1CCCC1'),
    E0 = (-199.316,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.63453,0.119898,-9.5698e-05,4.02339e-08,-6.99227e-12,-23766.4,46.877], Tmin=(100,'K'), Tmax=(1338.97,'K')), NASAPolynomial(coeffs=[19.8625,0.0556793,-2.37561e-05,4.41461e-09,-3.04473e-13,-29523.2,-63.1089], Tmin=(1338.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-199.316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(606.956,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cyclopentane) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = '[O][C](CO)C(OO)C1CCCC1(3896)',
    structure = SMILES('[O][C](CO)C(OO)C1CCCC1'),
    E0 = (-198.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.17339,0.123886,-0.00010121,4.27804e-08,-7.30585e-12,-23621.7,47.6425], Tmin=(100,'K'), Tmax=(1389.52,'K')), NASAPolynomial(coeffs=[24.1775,0.04803,-1.93225e-05,3.49249e-09,-2.37245e-13,-30944.8,-88.1537], Tmin=(1389.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-198.331,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(606.956,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cyclopentane) + radical(C2CsJOH) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]CC(O)C(OO)[C]1CCCC1(3897)',
    structure = SMILES('[O]CC(O)C(OO)[C]1CCCC1'),
    E0 = (-188.703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.59351,0.106699,-6.93382e-05,2.22139e-08,-2.86061e-12,-22480.7,52.0794], Tmin=(100,'K'), Tmax=(1794.57,'K')), NASAPolynomial(coeffs=[25.9171,0.0453792,-1.80829e-05,3.17283e-09,-2.07995e-13,-32354.5,-96.7301], Tmin=(1794.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-188.703,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(606.956,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cyclopentane) + radical(C2CJCOOH) + radical(CCOJ)"""),
)

species(
    label = '[O]C(CO)[C](OO)C1CCCC1(3898)',
    structure = SMILES('[O]C(CO)[C](OO)C1CCCC1'),
    E0 = (-188.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.90578,0.117863,-8.9429e-05,3.4894e-08,-5.5319e-12,-22396,47.5565], Tmin=(100,'K'), Tmax=(1478.78,'K')), NASAPolynomial(coeffs=[23.3772,0.0494737,-2.0058e-05,3.61985e-09,-2.44713e-13,-29873.6,-84.3102], Tmin=(1478.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-188.06,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(606.956,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cyclopentane) + radical(CC(C)OJ) + radical(C2CsJOOH)"""),
)

species(
    label = '[O]CC(O)C(O[O])C1CCCC1(3899)',
    structure = SMILES('[O]CC(O)C(O[O])C1CCCC1'),
    E0 = (-227.61,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.35508,0.110484,-7.90181e-05,2.94342e-08,-4.53356e-12,-27176.8,47.088], Tmin=(100,'K'), Tmax=(1490.75,'K')), NASAPolynomial(coeffs=[19.69,0.0540161,-2.22004e-05,4.02542e-09,-2.72518e-13,-33451.5,-62.8451], Tmin=(1490.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-227.61,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(611.114,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cyclopentane) + radical(ROOJ) + radical(CCOJ)"""),
)

species(
    label = '[O]CC(O)C(OO)C1[CH]CCC1(3900)',
    structure = SMILES('[O]CC(O)C(OO)C1[CH]CCC1'),
    E0 = (-185.072,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.22429,0.110505,-7.76007e-05,2.77918e-08,-4.1172e-12,-22068.2,48.4111], Tmin=(100,'K'), Tmax=(1531.62,'K')), NASAPolynomial(coeffs=[19.773,0.0556678,-2.38962e-05,4.41586e-09,-3.01653e-13,-28500.2,-61.8403], Tmin=(1531.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-185.072,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(606.956,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cyclopentane) + radical(Cs_S) + radical(CCOJ)"""),
)

species(
    label = '[O]C(CO)C(OO)[C]1CCCC1(3901)',
    structure = SMILES('[O]C(CO)C(OO)[C]1CCCC1'),
    E0 = (-184.048,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.49844,0.102695,-4.47563e-05,-1.16768e-08,1.05867e-11,-21922.3,50.6371], Tmin=(100,'K'), Tmax=(1064.53,'K')), NASAPolynomial(coeffs=[21.6584,0.0505286,-2.03499e-05,3.78771e-09,-2.66255e-13,-28827,-71.8031], Tmin=(1064.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-184.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(606.956,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cyclopentane) + radical(CC(C)OJ) + radical(C2CJCOOH)"""),
)

species(
    label = '[O]CC(O)C(OO)C1C[CH]CC1(3902)',
    structure = SMILES('[O]CC(O)C(OO)C1C[CH]CC1'),
    E0 = (-193.733,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.79385,0.108008,-7.61887e-05,2.86148e-08,-4.6078e-12,-23131,46.8653], Tmin=(100,'K'), Tmax=(1377.39,'K')), NASAPolynomial(coeffs=[14.4965,0.0636043,-2.78326e-05,5.21021e-09,-3.59809e-13,-27343.2,-31.7975], Tmin=(1377.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-193.733,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(606.956,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cyclopentane) + radical(cyclopentane) + radical(CCOJ)"""),
)

species(
    label = '[O]OC(C([O])CO)C1CCCC1(3903)',
    structure = SMILES('[O]OC(C([O])CO)C1CCCC1'),
    E0 = (-222.954,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.8457,0.112695,-7.32653e-05,1.63636e-08,1.37135e-12,-26591.3,47.7975], Tmin=(100,'K'), Tmax=(1097.73,'K')), NASAPolynomial(coeffs=[22.2079,0.0490909,-1.92085e-05,3.48778e-09,-2.40538e-13,-33320.9,-77.0888], Tmin=(1097.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-222.954,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(611.114,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cyclopentane) + radical(CC(C)OJ) + radical(ROOJ)"""),
)

species(
    label = '[O]C(CO)C(OO)C1[CH]CCC1(3904)',
    structure = SMILES('[O]C(CO)C(OO)C1[CH]CCC1'),
    E0 = (-180.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.22055,0.118167,-8.86675e-05,3.35822e-08,-5.09408e-12,-21459.6,50.9707], Tmin=(100,'K'), Tmax=(1561.54,'K')), NASAPolynomial(coeffs=[26.8745,0.04364,-1.7079e-05,3.01966e-09,-2.0119e-13,-30546.4,-102.364], Tmin=(1561.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-180.417,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(606.956,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cyclopentane) + radical(Cs_S) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]C(CO)C(OO)C1C[CH]CC1(3905)',
    structure = SMILES('[O]C(CO)C(OO)C1C[CH]CC1'),
    E0 = (-189.078,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.71715,0.11501,-8.5665e-05,3.30498e-08,-5.21043e-12,-22526.4,49.1467], Tmin=(100,'K'), Tmax=(1479.04,'K')), NASAPolynomial(coeffs=[22.0193,0.0508151,-2.05595e-05,3.70361e-09,-2.50032e-13,-29547.8,-74.6581], Tmin=(1479.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-189.078,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(606.956,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cyclopentane) + radical(CC(C)OJ) + radical(cyclopentane)"""),
)

species(
    label = 'O=C(CO)C(OO)C1CCCC1(3906)',
    structure = SMILES('O=C(CO)C(OO)C1CCCC1'),
    E0 = (-529.576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.56451,0.0994455,-2.29562e-05,-4.24603e-08,2.32022e-11,-63472.9,46.9577], Tmin=(100,'K'), Tmax=(1013.14,'K')), NASAPolynomial(coeffs=[24.8747,0.0455373,-1.78771e-05,3.37464e-09,-2.42853e-13,-71420.8,-93.7266], Tmin=(1013.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-529.576,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(606.956,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsCs) + ring(Cyclopentane)"""),
)

species(
    label = 'O=CC(O)C(OO)C1CCCC1(3907)',
    structure = SMILES('O=CC(O)C(OO)C1CCCC1'),
    E0 = (-525.703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.45352,0.0966767,-1.61872e-05,-4.85718e-08,2.51684e-11,-63010.6,47.6039], Tmin=(100,'K'), Tmax=(1010.51,'K')), NASAPolynomial(coeffs=[24.4082,0.0461218,-1.80595e-05,3.40741e-09,-2.45297e-13,-70882.9,-90.5248], Tmin=(1010.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-525.703,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(606.956,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-OdCsH) + ring(Cyclopentane)"""),
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
    label = '[O]CC1OOC1C1CCCC1(3908)',
    structure = SMILES('[O]CC1OOC1C1CCCC1'),
    E0 = (-50.4319,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (157.187,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.334496,0.079515,-1.27809e-05,-2.80956e-08,1.29862e-11,-5896.27,35.7763], Tmin=(100,'K'), Tmax=(1104.59,'K')), NASAPolynomial(coeffs=[15.3971,0.0546916,-2.27232e-05,4.25085e-09,-2.97538e-13,-11332.7,-50.561], Tmin=(1104.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-50.4319,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(12dioxetane) + ring(Cyclopentane) + radical(CCOJ)"""),
)

species(
    label = '[O]C1COOC1C1CCCC1(3909)',
    structure = SMILES('[O]C1COOC1C1CCCC1'),
    E0 = (-130.705,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (157.187,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.106648,0.0579335,8.84662e-05,-1.65999e-07,7.16682e-11,-15542.4,34.8296], Tmin=(100,'K'), Tmax=(923.123,'K')), NASAPolynomial(coeffs=[24.7512,0.0335895,-7.44321e-06,1.09726e-09,-7.95679e-14,-23683.9,-102.347], Tmin=(923.123,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-130.705,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(12dioxolane) + ring(Cyclopentane) + radical(CC(C)OJ)"""),
)

species(
    label = 'OOC(C1CCCC1)C1COO1(2211)',
    structure = SMILES('OOC(C1CCCC1)C1COO1'),
    E0 = (-204.632,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.39813,0.0995611,-3.32921e-05,-2.20726e-08,1.36364e-11,-24400.9,40.1843], Tmin=(100,'K'), Tmax=(1067.83,'K')), NASAPolynomial(coeffs=[21.0549,0.0534366,-2.1855e-05,4.09724e-09,-2.88992e-13,-31361.6,-79.7522], Tmin=(1067.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-204.632,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(615.271,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(12dioxetane) + ring(Cyclopentane)"""),
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
    label = '[O]C[CH]C(OO)C1CCCC1(3910)',
    structure = SMILES('[O]C[CH]C(OO)C1CCCC1'),
    E0 = (-6.72892,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3615,1310,387.5,850,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.19528,0.101964,-6.57737e-05,2.09862e-08,-2.7091e-12,-612.315,45.734], Tmin=(100,'K'), Tmax=(1772.3,'K')), NASAPolynomial(coeffs=[23.3502,0.0465661,-1.88877e-05,3.34969e-09,-2.21325e-13,-9312.75,-86.7306], Tmin=(1772.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-6.72892,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cyclopentane) + radical(CCOJ) + radical(CCJCOOH)"""),
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
    E0 = (-149.253,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-2.60929,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-52.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-82.4264,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-80.4618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-63.2381,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (16.4168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-60.2598,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-65.5944,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-73.9728,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-73.9728,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-48.4195,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-41.7245,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-90.4709,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-123.522,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-69.3391,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-69.5483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-79.7108,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-69.5483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-42.3525,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (189.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-85.8533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-85.8533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-18.5714,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-94.4431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-140.969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (236.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (239.463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]CC([O])C(OO)C1CCCC1(2207)'],
    products = ['CH2O(35)', 'O=CC(OO)C1CCCC1(1993)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(19)', '[O]CC(=O)C(OO)C1CCCC1(3889)'],
    products = ['[O]CC([O])C(OO)C1CCCC1(2207)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.0366254,'m^3/(mol*s)'), n=1.743, Ea=(71.4418,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsCs_O;YJ] for rate rule [CO-CsCs_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(19)', '[O]C(C=O)C(OO)C1CCCC1(3890)'],
    products = ['[O]CC([O])C(OO)C1CCCC1(2207)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(9.6e+09,'cm^3/(mol*s)'), n=0.935, Ea=(17.4473,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2782 used for CO-CsH_O;HJ
Exact match found for rate rule [CO-CsH_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['OO[CH]C1CCCC1(1997)', 'OCH2CHO(3891)'],
    products = ['[O]CC([O])C(OO)C1CCCC1(2207)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(8.09579,'m^3/(mol*s)'), n=1.29734, Ea=(31.9242,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;CsJ] for rate rule [CO-CsH_O;CsJ-OsCsH]
Euclidian distance = 2.2360679775
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2][O](354)', 'O=CC(OO)C1CCCC1(1993)'],
    products = ['[O]CC([O])C(OO)C1CCCC1(2207)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.0201871,'m^3/(mol*s)'), n=2.2105, Ea=(56.0866,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsH_O;YJ] for rate rule [CO-CsH_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2O(35)', '[O][CH]C(OO)C1CCCC1(2032)'],
    products = ['[O]CC([O])C(OO)C1CCCC1(2207)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(225.36,'m^3/(mol*s)'), n=0.996465, Ea=(58.8821,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO-HH_O;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]CC([O])C(OO)C1CCCC1(2207)'],
    products = ['[O]C[C](O)C(OO)C1CCCC1(3892)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.56178e+08,'s^-1'), n=1.25272, Ea=(165.67,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_Cs2] for rate rule [R2H_S;O_rad_out;Cs_H_out_Cs2]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]CC([O])C(OO)C1CCCC1(2207)'],
    products = ['[O]C([CH]O)C(OO)C1CCCC1(3893)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(153000,'s^-1'), n=2.26, Ea=(88.9937,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 325 used for R2H_S;O_rad_out;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R2H_S;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]CC([O])C(OO)C1CCCC1(2207)'],
    products = ['[O]CC(O)[C](OO)C1CCCC1(3894)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.0250998,'s^-1'), n=4.2935, Ea=(83.6591,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_OOH/Cs]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]CC([O])C(OO)C1CCCC1(2207)'],
    products = ['[O][CH]C(O)C(OO)C1CCCC1(3895)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(223829,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]CC([O])C(OO)C1CCCC1(2207)'],
    products = ['[O][C](CO)C(OO)C1CCCC1(3896)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]CC([O])C(OO)C1CCCC1(2207)'],
    products = ['[O]CC(O)C(OO)[C]1CCCC1(3897)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1e+11,'s^-1'), n=0, Ea=(100.834,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4H_SSS;O_rad_out;Cs_H_out_Cs2] for rate rule [R4H_SSS;O_rad_out;Cs_H_out_Cs2_cy5]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]CC([O])C(OO)C1CCCC1(2207)'],
    products = ['[O]C(CO)[C](OO)C1CCCC1(3898)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.76e+08,'s^-1'), n=1.2, Ea=(107.529,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;O_rad_out;Cs_H_out_OOH/Cs]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]CC([O])C(OO)C1CCCC1(2207)'],
    products = ['[O]CC(O)C(O[O])C1CCCC1(3899)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.37227e+06,'s^-1'), n=1.56745, Ea=(58.7826,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;O_rad_out;XH_out] for rate rule [R5H_SSSS;O_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O]CC([O])C(OO)C1CCCC1(2207)'],
    products = ['[O]CC(O)C(OO)C1[CH]CCC1(3900)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.6e+11,'s^-1'), n=0, Ea=(25.7316,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(1000,'K'), comment="""Estimated using template [R5H_CCC;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5H_CCC;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]CC([O])C(OO)C1CCCC1(2207)'],
    products = ['[O]C(CO)C(OO)[C]1CCCC1(3901)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.25e+10,'s^-1'), n=0, Ea=(79.9144,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5H_CCC;O_rad_out;Cs_H_out_Cs2] for rate rule [R5H_CCC;O_rad_out;Cs_H_out_Cs2_cy5]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]CC([O])C(OO)C1CCCC1(2207)'],
    products = ['[O]CC(O)C(OO)C1C[CH]CC1(3902)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(6.24e+09,'s^-1'), n=0, Ea=(79.7052,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6H_SSSSS;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R6H_SSSSS;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]CC([O])C(OO)C1CCCC1(2207)'],
    products = ['[O]OC(C([O])CO)C1CCCC1(3903)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(123903,'s^-1'), n=1.46258, Ea=(69.5427,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_SSSSS;O_rad_out;XH_out] for rate rule [R6H_SSSSS;O_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]CC([O])C(OO)C1CCCC1(2207)'],
    products = ['[O]C(CO)C(OO)C1[CH]CCC1(3904)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(6.24e+09,'s^-1'), n=0, Ea=(79.7052,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6H_SSSSS;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R6H_SSSSS;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]CC([O])C(OO)C1CCCC1(2207)'],
    products = ['[O]C(CO)C(OO)C1C[CH]CC1(3905)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.8e+08,'s^-1'), n=0, Ea=(106.901,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R7H;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R7H;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2][O](354)', '[O][CH]C(OO)C1CCCC1(2032)'],
    products = ['[O]CC([O])C(OO)C1CCCC1(2207)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(7.35017e+06,'m^3/(mol*s)'), n=0.0284742, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O]CC([O])C(OO)C1CCCC1(2207)'],
    products = ['O=C(CO)C(OO)C1CCCC1(3906)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O]CC([O])C(OO)C1CCCC1(2207)'],
    products = ['O=CC(O)C(OO)C1CCCC1(3907)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]CC([O])C(OO)C1CCCC1(2207)'],
    products = ['OH(21)', '[O]CC1OOC1C1CCCC1(3908)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.11355e+11,'s^-1'), n=0, Ea=(130.682,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3OO_SS;Y_rad_intra;OOH]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]CC([O])C(OO)C1CCCC1(2207)'],
    products = ['OH(21)', '[O]C1COOC1C1CCCC1(3909)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.63066e+10,'s^-1'), n=0, Ea=(54.8104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4OO_SSS;Y_rad_intra;OOH]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O]CC([O])C(OO)C1CCCC1(2207)'],
    products = ['OOC(C1CCCC1)C1COO1(2211)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['O(20)', '[O]C[CH]C(OO)C1CCCC1(3910)'],
    products = ['[O]CC([O])C(OO)C1CCCC1(2207)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/NonDeC;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['O(20)', '[CH2]C([O])C(OO)C1CCCC1(2063)'],
    products = ['[O]CC([O])C(OO)C1CCCC1(2207)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H2/Cs;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

network(
    label = '460',
    isomers = [
        '[O]CC([O])C(OO)C1CCCC1(2207)',
    ],
    reactants = [
        ('CH2O(35)', 'O=CC(OO)C1CCCC1(1993)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '460',
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

