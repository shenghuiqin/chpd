species(
    label = 'CC(C[CH]CC[O])OO(5265)',
    structure = SMILES('CC(C[CH]CC[O])OO'),
    E0 = (-32.494,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.438269,0.109394,-0.000153679,1.4041e-07,-5.23057e-11,-3759.73,39.2694], Tmin=(100,'K'), Tmax=(809.318,'K')), NASAPolynomial(coeffs=[4.03651,0.0652997,-3.12182e-05,5.97984e-09,-4.14574e-13,-3764.25,23.0747], Tmin=(809.318,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-32.494,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(CCOJ)"""),
)

species(
    label = 'CH2O(15)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.79372,-0.00990833,3.7322e-05,-3.79285e-08,1.31773e-11,-14379.2,0.602798], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.16953,0.00619321,-2.25056e-06,3.65976e-10,-2.20149e-14,-14548.7,6.04208], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-119.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH2O""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C=CCC(C)OO(1790)',
    structure = SMILES('C=CCC(C)OO'),
    E0 = (-141.267,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,245.622,245.624],'cm^-1')),
        HinderedRotor(inertia=(0.420552,'amu*angstrom^2'), symmetry=1, barrier=(18.0056,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00279424,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.201245,'amu*angstrom^2'), symmetry=1, barrier=(8.61544,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.609594,'amu*angstrom^2'), symmetry=1, barrier=(26.0973,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.6096,'amu*angstrom^2'), symmetry=1, barrier=(26.0973,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.132,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3800.23,'J/mol'), sigma=(6.47665,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=593.59 K, Pc=31.74 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.481618,0.070817,-5.22184e-05,1.96663e-08,-3.01737e-12,-16858.7,28.0792], Tmin=(100,'K'), Tmax=(1516.09,'K')), NASAPolynomial(coeffs=[15.5922,0.0309501,-1.27749e-05,2.32203e-09,-1.57362e-13,-21440.5,-51.1086], Tmin=(1516.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-141.267,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH)"""),
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
    label = 'CC(C=CCC[O])OO(5429)',
    structure = SMILES('CC(C=CCC[O])OO'),
    E0 = (-110.617,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (131.15,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.311926,0.0878457,-7.3546e-05,3.48578e-08,-7.27105e-12,-13177,35.6418], Tmin=(100,'K'), Tmax=(1084.34,'K')), NASAPolynomial(coeffs=[9.89393,0.0524986,-2.4649e-05,4.79509e-09,-3.3989e-13,-15255,-11.3617], Tmin=(1084.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-110.617,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ)"""),
)

species(
    label = 'CC(CC=CC[O])OO(5430)',
    structure = SMILES('CC(CC=CC[O])OO'),
    E0 = (-104.724,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (131.15,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.180138,0.0916733,-9.20246e-05,5.95706e-08,-1.74026e-11,-12464.6,35.6033], Tmin=(100,'K'), Tmax=(796.997,'K')), NASAPolynomial(coeffs=[7.05056,0.0571917,-2.71276e-05,5.28585e-09,-3.74692e-13,-13559.7,4.01638], Tmin=(796.997,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-104.724,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ)"""),
)

species(
    label = 'CC(C[CH]CC=O)OO(5431)',
    structure = SMILES('CC(C[CH]CC=O)OO'),
    E0 = (-175.727,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1310,387.5,850,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (131.15,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0662013,0.0925616,-8.55039e-05,4.54261e-08,-1.04373e-11,-20998.5,35.8798], Tmin=(100,'K'), Tmax=(1009.89,'K')), NASAPolynomial(coeffs=[10.5866,0.0508921,-2.36119e-05,4.56894e-09,-3.23095e-13,-23123.4,-14.9788], Tmin=(1009.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-175.727,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCC=O)"""),
)

species(
    label = 'CH3CHOOH(63)',
    structure = SMILES('C[CH]OO'),
    E0 = (10.3272,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
        HinderedRotor(inertia=(0.322978,'amu*angstrom^2'), symmetry=1, barrier=(53.3012,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.32071,'amu*angstrom^2'), symmetry=1, barrier=(53.2874,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.329436,'amu*angstrom^2'), symmetry=1, barrier=(53.2775,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (61.0599,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3284.22,'J/mol'), sigma=(4.037,'angstroms'), dipoleMoment=(1.3,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[6.26854,0.000111351,6.15944e-05,-7.95035e-08,3.12333e-11,962.571,-0.180971], Tmin=(200,'K'), Tmax=(999.99,'K')), NASAPolynomial(coeffs=[9.10837,0.0152964,-5.54864e-06,9.02286e-10,-5.4412e-14,-932.703,-20.3914], Tmin=(999.99,'K'), Tmax=(2500,'K'))], Tmin=(200,'K'), Tmax=(2500,'K'), E0=(10.3272,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), label="""CH3CHOOH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C=CCC[O](1706)',
    structure = SMILES('C=CCC[O]'),
    E0 = (49.3559,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,263.602,264.656,1655.13],'cm^-1')),
        HinderedRotor(inertia=(0.107573,'amu*angstrom^2'), symmetry=1, barrier=(5.36373,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109574,'amu*angstrom^2'), symmetry=1, barrier=(5.3866,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.54294,0.034052,-1.45097e-05,2.03585e-09,1.5555e-14,5986.72,18.8853], Tmin=(100,'K'), Tmax=(2019.44,'K')), NASAPolynomial(coeffs=[12.9947,0.0193782,-8.08802e-06,1.39415e-09,-8.80057e-14,536.132,-41.9277], Tmin=(2019.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(49.3559,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ)"""),
)

species(
    label = '[CH2][O](167)',
    structure = SMILES('[CH2][O]'),
    E0 = (192.903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88409,-0.00363885,3.28543e-05,-4.13611e-08,1.59631e-11,23210.8,7.47983], Tmin=(100,'K'), Tmax=(933.06,'K')), NASAPolynomial(coeffs=[6.69335,0.000289989,8.61416e-07,-1.56351e-10,7.33778e-15,21991.3,-9.6043], Tmin=(933.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(H3COJ) + radical(CsJOH)"""),
)

species(
    label = '[CH2][CH]CC(C)OO(1793)',
    structure = SMILES('[CH2][CH]CC(C)OO'),
    E0 = (130.643,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,200,800],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.132,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.878853,0.0738968,-7.16431e-05,4.40711e-08,-1.20927e-11,15820.5,30.7126], Tmin=(100,'K'), Tmax=(849.381,'K')), NASAPolynomial(coeffs=[7.02985,0.04493,-2.04882e-05,3.92062e-09,-2.7515e-13,14775.6,2.04174], Tmin=(849.381,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(130.643,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'CC([CH]CCC[O])OO(5432)',
    structure = SMILES('CC([CH]CCC[O])OO'),
    E0 = (-26.5373,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.556198,0.091311,-4.11257e-05,-1.00054e-07,1.21735e-10,-3083.79,35.7178], Tmin=(100,'K'), Tmax=(468.201,'K')), NASAPolynomial(coeffs=[7.32107,0.0604993,-2.88597e-05,5.57131e-09,-3.90112e-13,-4013,5.05643], Tmin=(468.201,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-26.5373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CCJCOOH)"""),
)

species(
    label = 'CC(CC[CH]C[O])OO(5406)',
    structure = SMILES('CC(CC[CH]C[O])OO'),
    E0 = (-27.0502,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.15004,0.0991426,-0.000102774,6.71478e-08,-1.94126e-11,-3110.9,38.1975], Tmin=(100,'K'), Tmax=(811.733,'K')), NASAPolynomial(coeffs=[8.11508,0.0584154,-2.75166e-05,5.34103e-09,-3.77652e-13,-4452.75,0.0468175], Tmin=(811.733,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-27.0502,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CCJCO)"""),
)

species(
    label = 'CC(C[CH]C[CH]O)OO(5433)',
    structure = SMILES('CC(C[CH]C[CH]O)OO'),
    E0 = (-77.9012,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1310,387.5,850,1000,200,800,1200,1600],'cm^-1')),
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
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.0233,0.120276,-0.000170143,1.43782e-07,-4.93626e-11,-9197.88,40.0609], Tmin=(100,'K'), Tmax=(814.273,'K')), NASAPolynomial(coeffs=[9.30354,0.0557605,-2.58988e-05,4.8927e-09,-3.36441e-13,-10422.6,-4.8318], Tmin=(814.273,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-77.9012,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(CCsJOH)"""),
)

species(
    label = 'C[C](CCCC[O])OO(5434)',
    structure = SMILES('C[C](CCCC[O])OO'),
    E0 = (-40.0532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.557697,0.111471,-0.000154926,1.39323e-07,-5.15544e-11,-4664.04,37.4283], Tmin=(100,'K'), Tmax=(797.839,'K')), NASAPolynomial(coeffs=[5.00007,0.0646013,-3.10768e-05,5.97935e-09,-4.16054e-13,-4946,15.6612], Tmin=(797.839,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-40.0532,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(C2CsJOOH) + radical(CCOJ)"""),
)

species(
    label = 'CC(CCC[CH][O])OO(5435)',
    structure = SMILES('CC(CCC[CH][O])OO'),
    E0 = (-46.6544,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.165692,0.0999408,-6.91638e-05,-5.3431e-08,8.78535e-11,-5489.04,34.3239], Tmin=(100,'K'), Tmax=(491.747,'K')), NASAPolynomial(coeffs=[8.65173,0.0586613,-2.78888e-05,5.36176e-09,-3.74033e-13,-6659.13,-4.00432], Tmin=(491.747,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-46.6544,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = 'CC(C[CH][CH]CO)OO(5436)',
    structure = SMILES('CC(C[CH][CH]CO)OO'),
    E0 = (-58.2971,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.420355,0.104466,-0.000122737,9.05023e-08,-2.86106e-11,-6859.02,40.8582], Tmin=(100,'K'), Tmax=(757.396,'K')), NASAPolynomial(coeffs=[8.88256,0.0553345,-2.5431e-05,4.85097e-09,-3.38424e-13,-8268.19,-1.43779], Tmin=(757.396,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-58.2971,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCJCO) + radical(RCCJCC)"""),
)

species(
    label = '[CH2]C(CCCC[O])OO(5437)',
    structure = SMILES('[CH2]C(CCCC[O])OO'),
    E0 = (-12.9897,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0945854,0.101862,-9.66306e-05,3.13448e-08,1.30226e-11,-1426.13,37.5773], Tmin=(100,'K'), Tmax=(559.827,'K')), NASAPolynomial(coeffs=[8.4078,0.0586511,-2.7845e-05,5.39447e-09,-3.79916e-13,-2652.94,-0.964001], Tmin=(559.827,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-12.9897,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CJCOOH) + radical(CCOJ)"""),
)

species(
    label = 'CC(CCCC[O])O[O](5438)',
    structure = SMILES('CC(CCCC[O])O[O]'),
    E0 = (-74.9476,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.379535,0.10469,-0.0001319,1.09814e-07,-3.89401e-11,-8864.42,37.2582], Tmin=(100,'K'), Tmax=(774.51,'K')), NASAPolynomial(coeffs=[6.52457,0.0600282,-2.79629e-05,5.33768e-09,-3.7114e-13,-9663.78,7.45765], Tmin=(774.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-74.9476,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(ROOJ)"""),
)

species(
    label = 'CC([CH][CH]CCO)OO(5439)',
    structure = SMILES('CC([CH][CH]CCO)OO'),
    E0 = (-57.7842,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.7454,0.113767,-0.000155626,1.31852e-07,-4.60189e-11,-6787.99,41.8191], Tmin=(100,'K'), Tmax=(801.336,'K')), NASAPolynomial(coeffs=[7.98604,0.0575798,-2.68605e-05,5.1004e-09,-3.52387e-13,-7782.72,4.1535], Tmin=(801.336,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-57.7842,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(CCJCOOH)"""),
)

species(
    label = 'C[C](C[CH]CCO)OO(5440)',
    structure = SMILES('C[C](C[CH]CCO)OO'),
    E0 = (-71.3,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.668879,0.114707,-0.000166472,1.49836e-07,-5.42326e-11,-8418.91,39.531], Tmin=(100,'K'), Tmax=(824.657,'K')), NASAPolynomial(coeffs=[5.53348,0.0619439,-2.92461e-05,5.55141e-09,-3.82097e-13,-8670.72,15.4796], Tmin=(824.657,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-71.3,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(C2CsJOOH)"""),
)

species(
    label = '[CH2]C(C[CH]CCO)OO(5441)',
    structure = SMILES('[CH2]C(C[CH]CCO)OO'),
    E0 = (-44.2366,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.815921,0.114181,-0.000151679,1.2334e-07,-4.16815e-11,-5154.88,41.7651], Tmin=(100,'K'), Tmax=(784.851,'K')), NASAPolynomial(coeffs=[9.41409,0.0550879,-2.54469e-05,4.82436e-09,-3.33641e-13,-6546.45,-3.74559], Tmin=(784.851,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-44.2366,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CJCOOH) + radical(RCCJCC)"""),
)

species(
    label = 'CC(C[CH]CCO)O[O](5442)',
    structure = SMILES('CC(C[CH]CCO)O[O]'),
    E0 = (-106.194,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.516278,0.10826,-0.000144767,1.22276e-07,-4.25575e-11,-12618.2,39.4506], Tmin=(100,'K'), Tmax=(815.405,'K')), NASAPolynomial(coeffs=[7.12869,0.0572423,-2.60546e-05,4.89078e-09,-3.3557e-13,-13415.7,6.8832], Tmin=(815.405,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-106.194,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(ROOJ)"""),
)

species(
    label = 'CC(CC=CCO)OO(5443)',
    structure = SMILES('CC(CC=CCO)OO'),
    E0 = (-330.429,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.184464,0.0947757,-8.23812e-05,3.9584e-08,-8.04873e-12,-39593.1,36.0286], Tmin=(100,'K'), Tmax=(1143.25,'K')), NASAPolynomial(coeffs=[13.0769,0.048377,-2.15038e-05,4.0844e-09,-2.85867e-13,-42625.3,-29.725], Tmin=(1143.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-330.429,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH)"""),
)

species(
    label = 'CC(CCCC=O)OO(5444)',
    structure = SMILES('CC(CCCC=O)OO'),
    E0 = (-375.629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0989805,0.0973327,-9.2158e-05,5.18326e-08,-1.27709e-11,-45036.2,33.9829], Tmin=(100,'K'), Tmax=(942.216,'K')), NASAPolynomial(coeffs=[9.68484,0.0557975,-2.60345e-05,5.0469e-09,-3.57162e-13,-46879.9,-12.6361], Tmin=(942.216,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-375.629,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH)"""),
)

species(
    label = 'CC(C=CCCO)OO(5445)',
    structure = SMILES('CC(C=CCCO)OO'),
    E0 = (-336.322,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.479596,0.095838,-8.03266e-05,3.53529e-08,-6.39548e-12,-40286.9,37.6074], Tmin=(100,'K'), Tmax=(1296.51,'K')), NASAPolynomial(coeffs=[16.9213,0.0421527,-1.82153e-05,3.41528e-09,-2.37095e-13,-44799,-50.8605], Tmin=(1296.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-336.322,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH)"""),
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
    label = '[O]CC[CH]CCOO(5446)',
    structure = SMILES('[O]CC[CH]CCOO'),
    E0 = (3.9292,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,3615,1310,387.5,850,1000,3025,407.5,1350,352.5,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (118.131,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.224688,0.0966163,-0.000150051,1.45744e-07,-5.54005e-11,595.348,34.9893], Tmin=(100,'K'), Tmax=(833.756,'K')), NASAPolynomial(coeffs=[1.74454,0.059258,-2.87476e-05,5.49823e-09,-3.79027e-13,1386.95,34.2003], Tmin=(833.756,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(3.9292,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + radical(CCOJ) + radical(RCCJCC)"""),
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
    label = 'CC1CC(CC[O])O1(5447)',
    structure = SMILES('CC1CC(CC[O])O1'),
    E0 = (-127.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.15,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.747732,0.0632235,-1.7976e-05,-1.70858e-08,1.08388e-11,-15158.2,28.0594], Tmin=(100,'K'), Tmax=(946.388,'K')), NASAPolynomial(coeffs=[10.7063,0.0409767,-1.41679e-05,2.38784e-09,-1.58226e-13,-17931.8,-24.1314], Tmin=(946.388,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-127.068,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C(CC[O])C(C)OO(5448)',
    structure = SMILES('[CH2]C(CC[O])C(C)OO'),
    E0 = (-24.2171,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,200,800,1200,1600],'cm^-1')),
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
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.495066,0.105833,-0.000122744,8.95182e-08,-2.79882e-11,-2757.2,38.4285], Tmin=(100,'K'), Tmax=(765.672,'K')), NASAPolynomial(coeffs=[9.00991,0.0561788,-2.54694e-05,4.8237e-09,-3.35105e-13,-4212.77,-4.88983], Tmin=(765.672,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-24.2171,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C[O])CC(C)OO(5266)',
    structure = SMILES('[CH2]C(C[O])CC(C)OO'),
    E0 = (-27.5643,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,200,800,1200,1600],'cm^-1')),
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
    molecularWeight = (132.158,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4518.77,'J/mol'), sigma=(7.72648,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=705.82 K, Pc=22.23 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.495066,0.105833,-0.000122744,8.95182e-08,-2.79882e-11,-3159.77,38.4285], Tmin=(100,'K'), Tmax=(765.672,'K')), NASAPolynomial(coeffs=[9.00991,0.0561788,-2.54694e-05,4.8237e-09,-3.35105e-13,-4615.34,-4.88983], Tmin=(765.672,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-27.5643,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Isobutyl)"""),
)

species(
    label = 'CC(CC1CCO1)OO(5269)',
    structure = SMILES('CC(CC1CCO1)OO'),
    E0 = (-281.269,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.306844,0.0830982,-3.75657e-05,-1.27588e-08,1.24456e-11,-33663.1,32.4389], Tmin=(100,'K'), Tmax=(938.991,'K')), NASAPolynomial(coeffs=[16.6967,0.0391969,-1.30135e-05,2.16939e-09,-1.44473e-13,-38114.2,-55.2206], Tmin=(938.991,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-281.269,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(486.397,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane)"""),
)

species(
    label = 'CC([O])CC(O)CC[O](5449)',
    structure = SMILES('CC([O])CC(O)CC[O]'),
    E0 = (-240.567,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.305749,0.101946,-0.000112582,7.75945e-08,-2.31378e-11,-28785,37.658], Tmin=(100,'K'), Tmax=(797.027,'K')), NASAPolynomial(coeffs=[8.89263,0.055783,-2.57064e-05,4.92921e-09,-3.45666e-13,-30251.3,-4.63225], Tmin=(797.027,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-240.567,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C[O](96)',
    structure = SMILES('[CH2]C[O]'),
    E0 = (188.892,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1398.33],'cm^-1')),
        HinderedRotor(inertia=(0.00547724,'amu*angstrom^2'), symmetry=1, barrier=(7.58298,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.57171,0.0102136,5.90913e-06,-7.99869e-09,2.07078e-12,22733,11.7517], Tmin=(100,'K'), Tmax=(1490.84,'K')), NASAPolynomial(coeffs=[4.741,0.01502,-6.91914e-06,1.31179e-09,-8.9824e-14,21501.6,2.68291], Tmin=(1490.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.892,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CJCO) + radical(CCOJ)"""),
)

species(
    label = '[CH]CC(C)OO(1806)',
    structure = SMILES('[CH]CC(C)OO'),
    E0 = (202.946,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.1051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10082,0.0630657,-5.47474e-05,2.52728e-08,-4.81693e-12,24513.6,24.135], Tmin=(100,'K'), Tmax=(1232.86,'K')), NASAPolynomial(coeffs=[11.8571,0.028167,-1.22867e-05,2.31217e-09,-1.60955e-13,21861.4,-30.0093], Tmin=(1232.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(202.946,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
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
    label = '[CH2]C[CH]CC(C)OO(5450)',
    structure = SMILES('[CH2]C[CH]CC(C)OO'),
    E0 = (106.875,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1310,387.5,850,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,200,800,1600],'cm^-1')),
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
    molecularWeight = (116.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.232502,0.0895718,-9.20837e-05,6.24051e-08,-1.89798e-11,12983.7,35.1079], Tmin=(100,'K'), Tmax=(771.288,'K')), NASAPolynomial(coeffs=[6.8284,0.0553663,-2.5564e-05,4.91109e-09,-3.44998e-13,11966.2,4.99908], Tmin=(771.288,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(106.875,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJ)"""),
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
    E0 = (-32.494,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (106.907,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (112.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (53.5122,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (89.0123,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (66.2712,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (70.4706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (132.455,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (126.498,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (51.3848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (116.247,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (111.547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (42.7867,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (61.9039,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-5.61177,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (18.3835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (37.0127,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (74.407,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (84.5032,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (323.546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (30.9062,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (30.9062,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-7.52072,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (423.791,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (43.3619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (135.718,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (132.371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-24.2097,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (78.8,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (391.837,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (349.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['CC(C[CH]CC[O])OO(5265)'],
    products = ['CH2O(15)', 'C=CCC(C)OO(1790)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(3)', 'CC(C=CCC[O])OO(5429)'],
    products = ['CC(C[CH]CC[O])OO(5265)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'CC(CC=CC[O])OO(5430)'],
    products = ['CC(C[CH]CC[O])OO(5265)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'CC(C[CH]CC=O)OO(5431)'],
    products = ['CC(C[CH]CC[O])OO(5265)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(9.6e+09,'cm^3/(mol*s)'), n=0.935, Ea=(17.4473,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2782 used for CO-CsH_O;HJ
Exact match found for rate rule [CO-CsH_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH3CHOOH(63)', 'C=CCC[O](1706)'],
    products = ['CC(C[CH]CC[O])OO(5265)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.0140176,'m^3/(mol*s)'), n=2.23645, Ea=(29.3292,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-HH_Cds-CsH;CsJ] for rate rule [Cds-HH_Cds-CsH;CsJ-OsCsH]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2][O](167)', 'C=CCC(C)OO(1790)'],
    products = ['CC(C[CH]CC[O])OO(5265)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.0172287,'m^3/(mol*s)'), n=2.32603, Ea=(14.6351,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CsH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2O(15)', '[CH2][CH]CC(C)OO(1793)'],
    products = ['CC(C[CH]CC[O])OO(5265)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(225.36,'m^3/(mol*s)'), n=0.996465, Ea=(58.8821,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO-HH_O;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CC([CH]CCC[O])OO(5432)'],
    products = ['CC(C[CH]CC[O])OO(5265)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6.76e+09,'s^-1'), n=0.88, Ea=(158.992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC] for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CC(C[CH]CC[O])OO(5265)'],
    products = ['CC(CC[CH]C[O])OO(5406)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(6.76e+09,'s^-1'), n=0.88, Ea=(158.992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 357 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CC(C[CH]C[CH]O)OO(5433)'],
    products = ['CC(C[CH]CC[O])OO(5265)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4500,'s^-1'), n=2.62, Ea=(129.286,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 322 used for R2H_S;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CC(C[CH]CC[O])OO(5265)'],
    products = ['C[C](CCCC[O])OO(5434)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.05815e+09,'s^-1'), n=0.95, Ea=(148.741,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out] for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_OOH/Cs]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CC(C[CH]CC[O])OO(5265)'],
    products = ['CC(CCC[CH][O])OO(5435)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.17661e+06,'s^-1'), n=1.79367, Ea=(144.041,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CC(C[CH]CC[O])OO(5265)'],
    products = ['CC(C[CH][CH]CO)OO(5436)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(223829,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(CCCC[O])OO(5437)'],
    products = ['CC(C[CH]CC[O])OO(5265)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(754000,'s^-1'), n=1.63, Ea=(74.8936,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 110 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CC(C[CH]CC[O])OO(5265)'],
    products = ['CC(CCCC[O])O[O](5438)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(0.0378492,'s^-1'), n=3.26, Ea=(26.8822,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;C_rad_out_H/NonDeC;XH_out] for rate rule [R5H_SSSS;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CC(C[CH]CC[O])OO(5265)'],
    products = ['CC([CH][CH]CCO)OO(5439)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.19599e+09,'s^-1'), n=0.63, Ea=(50.8774,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5HJ_3;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CC(C[CH]CC[O])OO(5265)'],
    products = ['C[C](C[CH]CCO)OO(5440)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2330.97,'s^-1'), n=1.70781, Ea=(69.5067,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;O_rad_out;Cs_H_out_OOH] for rate rule [R6HJ_3;O_rad_out;Cs_H_out_OOH/Cs]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CC(C[CH]CC[O])OO(5265)'],
    products = ['[CH2]C(C[CH]CCO)OO(5441)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(5.85e+08,'s^-1'), n=0, Ea=(106.901,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R7Hall;O_rad_out;Cs_H_out_2H] for rate rule [R7HJ_3;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CC(C[CH]CC[O])OO(5265)'],
    products = ['CC(C[CH]CCO)O[O](5442)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(40268.4,'s^-1'), n=2.03024, Ea=(116.997,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Y_rad_out;O_H_out] + [RnH;O_rad_out;XH_out] + [R8Hall;Y_rad_out;XH_out] for rate rule [R8Hall;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][O](167)', '[CH2][CH]CC(C)OO(1793)'],
    products = ['CC(C[CH]CC[O])OO(5265)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['CC(C[CH]CC[O])OO(5265)'],
    products = ['CC(CC=CCO)OO(5443)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CC(C[CH]CC[O])OO(5265)'],
    products = ['CC(CCCC=O)OO(5444)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CC(C[CH]CC[O])OO(5265)'],
    products = ['CC(C=CCCO)OO(5445)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CH2(S)(23)', '[O]CC[CH]CCOO(5446)'],
    products = ['CC(C[CH]CC[O])OO(5265)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(873476,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;Cs_H] for rate rule [carbene;C/H2/CsO]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['CC(C[CH]CC[O])OO(5265)'],
    products = ['OH(5)', 'CC1CC(CC[O])O1(5447)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.04e+11,'s^-1','*|/',1.74), n=0, Ea=(75.8559,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 9 used for R3OO_SS;C_rad/H/NonDeC_intra;OOH
Exact match found for rate rule [R3OO_SS;C_rad/H/NonDeC_intra;OOH]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(CC[O])C(C)OO(5448)'],
    products = ['CC(C[CH]CC[O])OO(5265)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(C[O])CC(C)OO(5266)'],
    products = ['CC(C[CH]CC[O])OO(5265)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CC(C[CH]CC[O])OO(5265)'],
    products = ['CC(CC1CCO1)OO(5269)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['CC(C[CH]CC[O])OO(5265)'],
    products = ['CC([O])CC(O)CC[O](5449)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.88e+10,'s^-1'), n=0, Ea=(111.294,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 4 used for R3OOH_SS;C_rad_out_H/NonDeC
Exact match found for rate rule [R3OOH_SS;C_rad_out_H/NonDeC]
Euclidian distance = 0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C[O](96)', '[CH]CC(C)OO(1806)'],
    products = ['CC(C[CH]CC[O])OO(5265)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction31',
    reactants = ['O(4)', '[CH2]C[CH]CC(C)OO(5450)'],
    products = ['CC(C[CH]CC[O])OO(5265)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H2/Cs;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

network(
    label = '1251',
    isomers = [
        'CC(C[CH]CC[O])OO(5265)',
    ],
    reactants = [
        ('CH2O(15)', 'C=CCC(C)OO(1790)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '1251',
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

