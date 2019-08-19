species(
    label = 'C=CC(C=CCC=O)O[O](1058)',
    structure = SMILES('C=CC(C=CCC=O)O[O]'),
    E0 = (-0.861799,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,492.5,1135,1000,2782.5,750,1395,475,1775,1000,2750,2850,1437.5,1250,1305,750,350,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4499.23,'J/mol'), sigma=(7.17214,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=702.77 K, Pc=27.67 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.121622,0.0822866,-5.79767e-05,1.92151e-08,-2.55063e-12,37.4582,37.6538], Tmin=(100,'K'), Tmax=(1718.24,'K')), NASAPolynomial(coeffs=[20.4137,0.0350468,-1.67365e-05,3.21395e-09,-2.22475e-13,-6935.79,-71.2276], Tmin=(1718.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-0.861799,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(ROOJ)"""),
)

species(
    label = 'O2(2)',
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
    label = '[CH2]C=CC=CCC=O(956)',
    structure = SMILES('[CH2]C=CC=CCC=O'),
    E0 = (40.6633,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,286.105,645.958],'cm^-1')),
        HinderedRotor(inertia=(0.0787087,'amu*angstrom^2'), symmetry=1, barrier=(22.8412,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.994264,'amu*angstrom^2'), symmetry=1, barrier=(22.8601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0789567,'amu*angstrom^2'), symmetry=1, barrier=(22.8332,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.994741,'amu*angstrom^2'), symmetry=1, barrier=(22.8711,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3882.24,'J/mol'), sigma=(6.33784,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=606.40 K, Pc=34.6 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.868578,0.0528089,9.16168e-06,-4.58481e-08,1.96156e-11,5017.29,28.9089], Tmin=(100,'K'), Tmax=(1039.72,'K')), NASAPolynomial(coeffs=[15.7717,0.0320753,-1.37307e-05,2.68843e-09,-1.96067e-13,-60.0845,-53.0844], Tmin=(1039.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(40.6633,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=C[CH]O[O](661)',
    structure = SMILES('C=C[CH]O[O]'),
    E0 = (199.472,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,492.5,1135,1000,3025,407.5,1350,352.5,316.51],'cm^-1')),
        HinderedRotor(inertia=(0.157895,'amu*angstrom^2'), symmetry=1, barrier=(11.1941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157042,'amu*angstrom^2'), symmetry=1, barrier=(11.2,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.0627,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.21608,0.0376502,-3.45386e-05,1.68249e-08,-3.30513e-12,24056.5,18.0455], Tmin=(100,'K'), Tmax=(1220.77,'K')), NASAPolynomial(coeffs=[9.39745,0.0141196,-5.62587e-06,1.03554e-09,-7.16583e-14,22303.1,-18.0331], Tmin=(1220.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(199.472,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(C=CCJO)"""),
)

species(
    label = '[CH]=CCC=O(986)',
    structure = SMILES('[CH]=CCC=O'),
    E0 = (153.957,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.12336,'amu*angstrom^2'), symmetry=1, barrier=(18.0046,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122606,'amu*angstrom^2'), symmetry=1, barrier=(18.0324,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.51056,0.0244663,1.05068e-05,-2.60457e-08,9.69719e-12,18577.3,17.4244], Tmin=(100,'K'), Tmax=(1119.61,'K')), NASAPolynomial(coeffs=[9.62674,0.019207,-9.46251e-06,1.93136e-09,-1.41855e-13,15720,-23.3551], Tmin=(1119.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(153.957,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C2H3(44)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36378,0.000265766,2.79621e-05,-3.72987e-08,1.5159e-11,34475,7.9151], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.15027,0.00754021,-2.62998e-06,4.15974e-10,-2.45408e-14,33856.6,1.72812], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(286.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""C2H3""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[O]O[CH]C=CCC=O(1070)',
    structure = SMILES('[O]O[CH]C=CCC=O'),
    E0 = (58.5047,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,950.295,951.652],'cm^-1')),
        HinderedRotor(inertia=(0.756491,'amu*angstrom^2'), symmetry=1, barrier=(17.3932,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.756228,'amu*angstrom^2'), symmetry=1, barrier=(17.3872,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.756034,'amu*angstrom^2'), symmetry=1, barrier=(17.3827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.186961,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00488,0.0553461,-2.6325e-05,-2.5933e-09,3.54175e-12,7152.97,30.1124], Tmin=(100,'K'), Tmax=(1201.11,'K')), NASAPolynomial(coeffs=[15.4835,0.0256248,-1.23068e-05,2.44712e-09,-1.75985e-13,2340.7,-47.9458], Tmin=(1201.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(58.5047,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(ROOJ) + radical(C=CCJO)"""),
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
    label = 'C=C[C](C=CCC=O)O[O](1071)',
    structure = SMILES('[CH2]C=C(C=CCC=O)O[O]'),
    E0 = (103.27,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,492.5,1135,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.00507727,'amu*angstrom^2'), symmetry=1, barrier=(4.13621,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179903,'amu*angstrom^2'), symmetry=1, barrier=(4.13631,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.925853,'amu*angstrom^2'), symmetry=1, barrier=(21.2872,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.925842,'amu*angstrom^2'), symmetry=1, barrier=(21.2869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.925819,'amu*angstrom^2'), symmetry=1, barrier=(21.2864,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.279095,0.0817101,-6.36646e-05,2.50406e-08,-4.0635e-12,12554,35.0607], Tmin=(100,'K'), Tmax=(1412.48,'K')), NASAPolynomial(coeffs=[15.4914,0.0386303,-1.79155e-05,3.44778e-09,-2.41716e-13,8256.56,-43.5835], Tmin=(1412.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(103.27,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + radical(ROOJ) + radical(C=CC=CCJ)"""),
)

species(
    label = 'CH2CHO(51)',
    structure = SMILES('[CH2]C=O'),
    E0 = (1.22925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,526.75,532.597,975.94,1639.13,1641.45],'cm^-1')),
        HinderedRotor(inertia=(0.00114821,'amu*angstrom^2'), symmetry=1, barrier=(2.1986,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66874,0.0096233,1.60617e-05,-2.87682e-08,1.2503e-11,219.438,12.5694], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.91637,0.0088465,-3.14955e-06,5.05413e-10,-3.01305e-14,-1047.8,-6.1065], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(1.22925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CH2CHO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH]=CC(C=C)O[O](664)',
    structure = SMILES('[CH]=CC(C=C)O[O]'),
    E0 = (381.522,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3120,650,792.5,1650,492.5,1135,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(0.602836,'amu*angstrom^2'), symmetry=1, barrier=(13.8604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.603153,'amu*angstrom^2'), symmetry=1, barrier=(13.8677,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.603051,'amu*angstrom^2'), symmetry=1, barrier=(13.8653,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31348,0.0588054,-5.20174e-05,2.42616e-08,-4.65637e-12,45983.4,27.4687], Tmin=(100,'K'), Tmax=(1226.44,'K')), NASAPolynomial(coeffs=[11.5055,0.0255643,-1.13616e-05,2.16191e-09,-1.51506e-13,43483.5,-23.7822], Tmin=(1226.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.522,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = 'HCO(30)',
    structure = SMILES('[CH]=O'),
    E0 = (32.4782,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1131.19,1955.83,1955.83],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.23755,-0.00332075,1.4003e-05,-1.3424e-08,4.37416e-12,3872.41,3.30835], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.92002,0.00252279,-6.71004e-07,1.05616e-10,-7.43798e-15,3653.43,3.58077], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(32.4782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HCO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH2]C=CC(C=C)O[O](1072)',
    structure = SMILES('C=C[CH]C(C=C)O[O]'),
    E0 = (184.546,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,492.5,1135,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,481.172,481.815,482.371],'cm^-1')),
        HinderedRotor(inertia=(0.140412,'amu*angstrom^2'), symmetry=1, barrier=(23.133,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.140035,'amu*angstrom^2'), symmetry=1, barrier=(23.1061,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138139,'amu*angstrom^2'), symmetry=1, barrier=(23.1117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13962,'amu*angstrom^2'), symmetry=1, barrier=(23.1724,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.629575,0.0579577,-3.71624e-06,-3.98636e-08,1.98936e-11,22331.5,32.7255], Tmin=(100,'K'), Tmax=(992.999,'K')), NASAPolynomial(coeffs=[18.2648,0.0242857,-9.29672e-06,1.77798e-09,-1.30799e-13,16986.9,-61.5065], Tmin=(992.999,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(184.546,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = 'C=CC(C=C[CH]C=O)O[O](1073)',
    structure = SMILES('C=CC([CH]C=CC=O)O[O]'),
    E0 = (61.0584,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,492.5,1135,1000,2782.5,750,1395,475,1775,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.654416,0.0878171,-6.88985e-05,2.57608e-08,-3.78524e-12,7523.34,40.4484], Tmin=(100,'K'), Tmax=(1623.01,'K')), NASAPolynomial(coeffs=[24.5909,0.0256001,-1.13985e-05,2.1428e-09,-1.47334e-13,-671.558,-93.5723], Tmin=(1623.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(61.0584,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(ROOJ) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = 'C=CC([C]=CCC=O)O[O](1074)',
    structure = SMILES('C=CC([C]=CCC=O)O[O]'),
    E0 = (236.98,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,492.5,1135,1000,1685,370,2782.5,750,1395,475,1775,1000,2750,2850,1437.5,1250,1305,750,350,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.637267,0.0806911,-6.24452e-05,2.40698e-08,-3.90363e-12,28616.8,36.2607], Tmin=(100,'K'), Tmax=(1372.24,'K')), NASAPolynomial(coeffs=[13.3805,0.0435455,-2.18414e-05,4.34361e-09,-3.09846e-13,25119.4,-29.2506], Tmin=(1372.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(236.98,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(ROOJ)"""),
)

species(
    label = 'C=CC(C=[C]CC=O)O[O](1075)',
    structure = SMILES('C=CC(C=[C]CC=O)O[O]'),
    E0 = (236.98,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,492.5,1135,1000,1685,370,2782.5,750,1395,475,1775,1000,2750,2850,1437.5,1250,1305,750,350,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.637267,0.0806911,-6.24452e-05,2.40698e-08,-3.90363e-12,28616.8,36.2607], Tmin=(100,'K'), Tmax=(1372.24,'K')), NASAPolynomial(coeffs=[13.3805,0.0435455,-2.18414e-05,4.34361e-09,-3.09846e-13,25119.4,-29.2506], Tmin=(1372.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(236.98,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C(C=CCC=O)O[O](1076)',
    structure = SMILES('C=[C]C(C=CCC=O)O[O]'),
    E0 = (236.98,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,492.5,1135,1000,1685,370,2782.5,750,1395,475,1775,1000,2750,2850,1437.5,1250,1305,750,350,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.637267,0.0806911,-6.24452e-05,2.40698e-08,-3.90363e-12,28616.8,36.2607], Tmin=(100,'K'), Tmax=(1372.24,'K')), NASAPolynomial(coeffs=[13.3805,0.0435455,-2.18414e-05,4.34361e-09,-3.09846e-13,25119.4,-29.2506], Tmin=(1372.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(236.98,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = 'C=CC(C=CC[C]=O)O[O](1077)',
    structure = SMILES('C=CC(C=CC[C]=O)O[O]'),
    E0 = (159.099,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.255404,0.083,-6.6093e-05,2.59157e-08,-4.15605e-12,19268.7,37.9818], Tmin=(100,'K'), Tmax=(1432.55,'K')), NASAPolynomial(coeffs=[16.7037,0.0370725,-1.80031e-05,3.53606e-09,-2.50489e-13,14556.1,-47.2844], Tmin=(1432.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(159.099,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]=CC(C=CCC=O)O[O](1078)',
    structure = SMILES('[CH]=CC(C=CCC=O)O[O]'),
    E0 = (246.234,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1380,1390,370,380,2900,435,492.5,1135,1000,3120,650,792.5,1650,2782.5,750,1395,475,1775,1000,2750,2850,1437.5,1250,1305,750,350,283.106,283.108],'cm^-1')),
        HinderedRotor(inertia=(0.27465,'amu*angstrom^2'), symmetry=1, barrier=(15.621,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.274646,'amu*angstrom^2'), symmetry=1, barrier=(15.621,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00210327,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.274646,'amu*angstrom^2'), symmetry=1, barrier=(15.621,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.274648,'amu*angstrom^2'), symmetry=1, barrier=(15.621,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.356095,0.0824041,-6.38077e-05,2.40542e-08,-3.71821e-12,29743.5,37.201], Tmin=(100,'K'), Tmax=(1469.73,'K')), NASAPolynomial(coeffs=[16.4614,0.0385723,-1.90736e-05,3.76303e-09,-2.66729e-13,25009.4,-46.6999], Tmin=(1469.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(246.234,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = 'C=CC1OOC1[CH]CC=O(1060)',
    structure = SMILES('C=CC1OOC1[CH]CC=O'),
    E0 = (85.1633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.418037,0.0839118,-6.6651e-05,2.85707e-08,-5.28187e-12,10367.2,33.4977], Tmin=(100,'K'), Tmax=(1217.63,'K')), NASAPolynomial(coeffs=[11.3555,0.0479813,-2.23882e-05,4.33631e-09,-3.06151e-13,7703.64,-21.4231], Tmin=(1217.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(85.1633,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(440.667,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(12dioxetane) + radical(CCJCOOH)"""),
)

species(
    label = '[CH2]C1OOC1C=CCC=O(1061)',
    structure = SMILES('[CH2]C1OOC1C=CCC=O'),
    E0 = (93.7828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.244378,0.0801989,-4.76438e-05,8.62665e-09,7.43503e-13,11443,35.3812], Tmin=(100,'K'), Tmax=(1306.54,'K')), NASAPolynomial(coeffs=[20.3471,0.035078,-1.64156e-05,3.19025e-09,-2.24993e-13,4532.75,-75.32], Tmin=(1306.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(93.7828,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(440.667,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(12dioxetane) + radical(CJCOOH)"""),
)

species(
    label = 'C=CC1C=CCC([O])OO1(1062)',
    structure = SMILES('C=CC1C=CCC([O])OO1'),
    E0 = (39.3275,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.414962,0.0265339,0.000178124,-2.81738e-07,1.18666e-10,4907.62,35.7368], Tmin=(100,'K'), Tmax=(923.757,'K')), NASAPolynomial(coeffs=[38.7877,-0.00462955,9.52064e-06,-1.85853e-09,1.06589e-13,-7941.56,-177.522], Tmin=(923.757,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(39.3275,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cycloheptane) + radical(CCOJ)"""),
)

species(
    label = 'C=C[C](C=CCC=O)OO(1063)',
    structure = SMILES('[CH2]C=C(C=CCC=O)OO'),
    E0 = (-48.7348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.2089,0.0877909,-6.69123e-05,2.49405e-08,-3.75672e-12,-5706.59,36.0625], Tmin=(100,'K'), Tmax=(1537.54,'K')), NASAPolynomial(coeffs=[19.541,0.036412,-1.67893e-05,3.20824e-09,-2.23212e-13,-11780,-67.7161], Tmin=(1537.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-48.7348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=CC([C]=CCC=O)OO(1064)',
    structure = SMILES('C=CC([C]=CCC=O)OO'),
    E0 = (84.9753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.193836,0.0863742,-6.47559e-05,2.31912e-08,-3.38735e-12,10353.8,37.092], Tmin=(100,'K'), Tmax=(1545.52,'K')), NASAPolynomial(coeffs=[17.9083,0.0405267,-2.02585e-05,3.99698e-09,-2.82508e-13,4878.18,-56.0819], Tmin=(1545.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.9753,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C(C=CCC=O)OO(1065)',
    structure = SMILES('C=[C]C(C=CCC=O)OO'),
    E0 = (84.9753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.193836,0.0863742,-6.47559e-05,2.31912e-08,-3.38735e-12,10353.8,37.092], Tmin=(100,'K'), Tmax=(1545.52,'K')), NASAPolynomial(coeffs=[17.9083,0.0405267,-2.02585e-05,3.99698e-09,-2.82508e-13,4878.18,-56.0819], Tmin=(1545.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.9753,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C=CC(C=[C]CC=O)OO(1066)',
    structure = SMILES('C=CC(C=[C]CC=O)OO'),
    E0 = (84.9753,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,1685,370,2782.5,750,1395,475,1775,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.193836,0.0863742,-6.47559e-05,2.31912e-08,-3.38735e-12,10353.8,37.092], Tmin=(100,'K'), Tmax=(1545.52,'K')), NASAPolynomial(coeffs=[17.9083,0.0405267,-2.02585e-05,3.99698e-09,-2.82508e-13,4878.18,-56.0819], Tmin=(1545.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.9753,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC(C=CCC=O)OO(1067)',
    structure = SMILES('[CH]=CC(C=CCC=O)OO'),
    E0 = (94.2297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.144084,0.0886077,-6.73941e-05,2.4283e-08,-3.51305e-12,11483.5,38.248], Tmin=(100,'K'), Tmax=(1591.59,'K')), NASAPolynomial(coeffs=[20.8309,0.0358933,-1.77133e-05,3.47334e-09,-2.44372e-13,4806.8,-72.692], Tmin=(1591.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(94.2297,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C=CC(C=C[CH]C=O)OO(1068)',
    structure = SMILES('C=CC([CH]C=CC=O)OO'),
    E0 = (-90.9463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.802163,0.0902796,-6.12697e-05,1.39311e-08,5.9333e-13,-10752.7,40.2026], Tmin=(100,'K'), Tmax=(1196.68,'K')), NASAPolynomial(coeffs=[23.008,0.0319646,-1.48378e-05,2.91846e-09,-2.09227e-13,-17974.5,-85.3068], Tmin=(1196.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-90.9463,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = 'C=CC(C=CC[C]=O)OO(1069)',
    structure = SMILES('C=CC(C=CC[C]=O)OO'),
    E0 = (7.0941,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.236878,0.0891235,-6.94569e-05,2.59276e-08,-3.88374e-12,1008.37,38.9996], Tmin=(100,'K'), Tmax=(1549.76,'K')), NASAPolynomial(coeffs=[20.8342,0.0347385,-1.68186e-05,3.28409e-09,-2.31036e-13,-5522.71,-71.8876], Tmin=(1549.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(7.0941,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O)"""),
)

species(
    label = 'C=CC1[CH]C(CC=O)OO1(1079)',
    structure = SMILES('C=CC1[CH]C(CC=O)OO1'),
    E0 = (0.234255,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0278075,0.0737231,-1.84093e-05,-3.4976e-08,2.11808e-11,186.821,36.03], Tmin=(100,'K'), Tmax=(933.557,'K')), NASAPolynomial(coeffs=[18.7669,0.0306442,-9.36577e-06,1.53678e-09,-1.04356e-13,-4954.31,-62.0917], Tmin=(933.557,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(0.234255,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(12dioxolane) + radical(CCJCOOH)"""),
)

species(
    label = 'O=CCC=CC1[CH]COO1(1080)',
    structure = SMILES('O=CCC=CC1[CH]COO1'),
    E0 = (7.94899,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.251502,0.0619811,1.78887e-05,-7.0599e-08,3.22464e-11,1109.41,35.8497], Tmin=(100,'K'), Tmax=(967.939,'K')), NASAPolynomial(coeffs=[19.7317,0.0309784,-1.07745e-05,1.97502e-09,-1.43619e-13,-4980.5,-69.4739], Tmin=(967.939,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(7.94899,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(12dioxolane) + radical(CCJCOOH)"""),
)

species(
    label = 'C=CC1C=CC[CH]OOO1(1081)',
    structure = SMILES('C=CC1C=CC[CH]OOO1'),
    E0 = (252.908,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.954972,0.0399867,7.97733e-05,-1.27988e-07,4.97448e-11,30551.4,32.3506], Tmin=(100,'K'), Tmax=(993.441,'K')), NASAPolynomial(coeffs=[18.6495,0.035537,-1.43632e-05,2.86467e-09,-2.16221e-13,23739.6,-69.4883], Tmin=(993.441,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(252.908,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclooctane) + radical(CCsJOO)"""),
)

species(
    label = 'CO(28)',
    structure = SMILES('[C-]#[O+]'),
    E0 = (-119.219,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2084.51],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0101,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(762.44,'J/mol'), sigma=(3.69,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.5971,-0.00102424,2.83336e-06,-1.75825e-09,3.42587e-13,-14343.2,3.45822], Tmin=(100,'K'), Tmax=(1669.93,'K')), NASAPolynomial(coeffs=[2.92796,0.00181931,-8.35308e-07,1.51269e-10,-9.88872e-15,-14292.7,6.51157], Tmin=(1669.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-119.219,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""CO""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'C=CC(C=CC)O[O](1082)',
    structure = SMILES('C=CC(C=CC)O[O]'),
    E0 = (98.4007,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,492.5,1135,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,255.377,255.378],'cm^-1')),
        HinderedRotor(inertia=(0.244869,'amu*angstrom^2'), symmetry=1, barrier=(11.3326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244872,'amu*angstrom^2'), symmetry=1, barrier=(11.3326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.24487,'amu*angstrom^2'), symmetry=1, barrier=(11.3326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244869,'amu*angstrom^2'), symmetry=1, barrier=(11.3326,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.834605,0.069804,-5.23622e-05,2.06657e-08,-3.42921e-12,11948.3,30.5855], Tmin=(100,'K'), Tmax=(1368.86,'K')), NASAPolynomial(coeffs=[12.1333,0.0367878,-1.61829e-05,3.0456e-09,-2.11181e-13,8855,-27.4715], Tmin=(1368.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(98.4007,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ)"""),
)

species(
    label = 'HO2(26)',
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
    label = 'C=CC=C=CCC=O(1001)',
    structure = SMILES('C=CC=C=CCC=O'),
    E0 = (99.2119,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2750,2850,1437.5,1250,1305,750,350,540,610,2055,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.977632,'amu*angstrom^2'), symmetry=1, barrier=(22.4777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.978595,'amu*angstrom^2'), symmetry=1, barrier=(22.4998,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.976925,'amu*angstrom^2'), symmetry=1, barrier=(22.4614,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.61164,0.0614653,-2.31198e-05,-1.21261e-08,7.8506e-12,12065.4,27.3318], Tmin=(100,'K'), Tmax=(1105.39,'K')), NASAPolynomial(coeffs=[16.7958,0.0283936,-1.28352e-05,2.53455e-09,-1.83671e-13,6930.01,-59.4137], Tmin=(1105.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(99.2119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=C=CC=CCC=O(976)',
    structure = SMILES('C=C=CC=CCC=O'),
    E0 = (99.2119,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2750,2850,1437.5,1250,1305,750,350,540,610,2055,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.977632,'amu*angstrom^2'), symmetry=1, barrier=(22.4777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.978595,'amu*angstrom^2'), symmetry=1, barrier=(22.4998,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.976925,'amu*angstrom^2'), symmetry=1, barrier=(22.4614,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.138,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.61164,0.0614653,-2.31198e-05,-1.21261e-08,7.8506e-12,12065.4,27.3318], Tmin=(100,'K'), Tmax=(1105.39,'K')), NASAPolynomial(coeffs=[16.7958,0.0283936,-1.28352e-05,2.53455e-09,-1.83671e-13,6930.01,-59.4137], Tmin=(1105.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(99.2119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=COC=CC(C=C)O[O](1083)',
    structure = SMILES('C=COC=CC(C=C)O[O]'),
    E0 = (62.0907,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,492.5,1135,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.21065,0.0997867,-9.1631e-05,4.20842e-08,-7.56951e-12,7667.85,41.8166], Tmin=(100,'K'), Tmax=(1355.34,'K')), NASAPolynomial(coeffs=[23.9579,0.0255065,-9.42197e-06,1.64675e-09,-1.1055e-13,845.511,-87.2595], Tmin=(1355.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(62.0907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ)"""),
)

species(
    label = 'C=CC(C=CC=CO)O[O](1084)',
    structure = SMILES('C=CC(C=CC=CO)O[O]'),
    E0 = (-22.5945,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,492.5,1135,1000,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,3615,1277.5,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.24912,0.0991465,-7.06446e-05,3.93328e-09,1.02934e-11,-2513.69,38.756], Tmin=(100,'K'), Tmax=(953.47,'K')), NASAPolynomial(coeffs=[27.0215,0.0191784,-5.61618e-06,9.61008e-10,-6.96718e-14,-9660.8,-105.495], Tmin=(953.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-22.5945,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ)"""),
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
    label = 'C=CC([O])C=CCC=O(1085)',
    structure = SMILES('C=CC([O])C=CCC=O'),
    E0 = (5.98929,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,409.879,409.906,409.947,409.989],'cm^-1')),
        HinderedRotor(inertia=(0.00100377,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143932,'amu*angstrom^2'), symmetry=1, barrier=(17.1646,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143909,'amu*angstrom^2'), symmetry=1, barrier=(17.1642,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144156,'amu*angstrom^2'), symmetry=1, barrier=(17.1636,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.234215,0.0701684,-3.36424e-05,-1.59749e-09,3.63715e-12,866.261,35.5564], Tmin=(100,'K'), Tmax=(1231.07,'K')), NASAPolynomial(coeffs=[17.9595,0.0341809,-1.6119e-05,3.16936e-09,-2.26007e-13,-5135.16,-60.2923], Tmin=(1231.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.98929,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CC(C)OJ)"""),
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
    E0 = (32.0416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (356.391,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (347.827,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (315.062,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (383.848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (217.025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (272.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (448.772,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (448.772,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (448.772,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (370.891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (458.026,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (85.8655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (121.181,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (95.0187,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (117.964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (144.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (144.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (118.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (110.583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (76.1985,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (89.5216,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (57.2958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (55.1747,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (252.908,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (336.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (177.795,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (177.795,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (375.891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (121.272,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (248.994,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O2(2)', '[CH2]C=CC=CCC=O(956)'],
    products = ['C=CC(C=CCC=O)O[O](1058)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5.9277e+07,'m^3/(mol*s)'), n=-0.00388889, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H/CdCd;Y_rad] for rate rule [C_rad/H/CdCd;O2_birad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -0.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=C[CH]O[O](661)', '[CH]=CCC=O(986)'],
    products = ['C=CC(C=CCC=O)O[O](1058)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.59671e+07,'m^3/(mol*s)'), n=0.0113737, Ea=(2.96199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H/OneDe;Y_rad] for rate rule [C_rad/H/OneDeO;Cd_pri_rad]
Euclidian distance = 2.2360679775
family: R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C2H3(44)', '[O]O[CH]C=CCC=O(1070)'],
    products = ['C=CC(C=CCC=O)O[O](1058)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.59671e+07,'m^3/(mol*s)'), n=0.0113737, Ea=(2.96199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;C_rad/H/OneDe] for rate rule [Cd_pri_rad;C_rad/H/OneDeO]
Euclidian distance = 2.2360679775
family: R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(19)', 'C=C[C](C=CCC=O)O[O](1071)'],
    products = ['C=CC(C=CCC=O)O[O](1058)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.62e+13,'cm^3/(mol*s)'), n=0.228, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/TwoDe;H_rad] for rate rule [C_rad/TwoDeO;H_rad]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['CH2CHO(51)', '[CH]=CC(C=C)O[O](664)'],
    products = ['C=CC(C=CCC=O)O[O](1058)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(6.41031e+09,'m^3/(mol*s)'), n=-0.9855, Ea=(1.09644,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_rad;C_pri_rad] for rate rule [Cd_pri_rad;C_rad/H2/CO]
Euclidian distance = 1.41421356237
family: R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['HCO(30)', '[CH2]C=CC(C=C)O[O](1072)'],
    products = ['C=CC(C=CCC=O)O[O](1058)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.81e+13,'cm^3/(mol*s)','*|/',3), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [C_pri_rad;CO_pri_rad] for rate rule [C_rad/H2/Cd;CO_pri_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(19)', 'C=CC(C=C[CH]C=O)O[O](1073)'],
    products = ['C=CC(C=CCC=O)O[O](1058)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(78817,'m^3/(mol*s)'), n=0.288419, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [C_rad/H/TwoDe;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(19)', 'C=CC([C]=CCC=O)O[O](1074)'],
    products = ['C=CC(C=CCC=O)O[O](1058)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(19)', 'C=CC(C=[C]CC=O)O[O](1075)'],
    products = ['C=CC(C=CCC=O)O[O](1058)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(19)', 'C=[C]C(C=CCC=O)O[O](1076)'],
    products = ['C=CC(C=CCC=O)O[O](1058)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(19)', 'C=CC(C=CC[C]=O)O[O](1077)'],
    products = ['C=CC(C=CCC=O)O[O](1058)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.53107e+07,'m^3/(mol*s)'), n=-0.133333, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_rad/NonDe;Y_rad] for rate rule [CO_rad/NonDe;H_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(19)', '[CH]=CC(C=CCC=O)O[O](1078)'],
    products = ['C=CC(C=CCC=O)O[O](1058)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.21e+14,'cm^3/(mol*s)','+|-',4.82e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 60 used for H_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=CC(C=CCC=O)O[O](1058)'],
    products = ['C=CC1OOC1[CH]CC=O(1060)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(8.65e+06,'s^-1'), n=1.3, Ea=(86.7273,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_HNd_pri;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_HNd_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=CC(C=CCC=O)O[O](1058)'],
    products = ['[CH2]C1OOC1C=CCC=O(1061)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 33 used for R5_SS_D;doublebond_intra_2H_pri;radadd_intra_O
Exact match found for rate rule [R5_SS_D;doublebond_intra_2H_pri;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=CC(C=CCC=O)O[O](1058)'],
    products = ['C=CC1C=CCC([O])OO1(1062)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.3251e+11,'s^-1'), n=0.157, Ea=(95.8805,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7plus;multiplebond_intra;radadd_intra_O] + [R8;multiplebond_intra;radadd_intra] for rate rule [R8;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=CC(C=CCC=O)O[O](1058)'],
    products = ['C=C[C](C=CCC=O)OO(1063)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.71966e+06,'s^-1'), n=1.65, Ea=(118.826,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Y_rad_out;Cs_H_out_TwoDe] for rate rule [R3H_SS_O;O_rad_out;Cs_H_out_CdCd]
Euclidian distance = 1.73205080757
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=CC(C=CCC=O)O[O](1058)'],
    products = ['C=CC([C]=CCC=O)OO(1064)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(274,'s^-1'), n=3.09, Ea=(145.603,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 337 used for R4H_SSS_OCs;O_rad_out;Cd_H_out_doubleC
Exact match found for rate rule [R4H_SSS_OCs;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=CC(C=CCC=O)O[O](1058)'],
    products = ['C=[C]C(C=CCC=O)OO(1065)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(274,'s^-1'), n=3.09, Ea=(145.603,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 337 used for R4H_SSS_OCs;O_rad_out;Cd_H_out_doubleC
Exact match found for rate rule [R4H_SSS_OCs;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=CC(C=[C]CC=O)OO(1066)'],
    products = ['C=CC(C=CCC=O)O[O](1058)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_single;XH_out] for rate rule [R5H_DSSS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 3.74165738677
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=CC(C=CCC=O)O[O](1058)'],
    products = ['[CH]=CC(C=CCC=O)OO(1067)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.468e+06,'s^-1','*|/',3), n=1.554, Ea=(111.445,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 290 used for R5H_SSSD;O_rad_out;Cd_H_out_singleH
Exact match found for rate rule [R5H_SSSD;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=CC(C=CCC=O)O[O](1058)'],
    products = ['C=CC(C=C[CH]C=O)OO(1068)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(16979.7,'s^-1'), n=1.95635, Ea=(77.0603,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R6H_RSSMS;O_rad_out;Cs_H_out_H/CO]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=CC(C=CC[C]=O)OO(1069)'],
    products = ['C=CC(C=CCC=O)O[O](1058)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.21847e+06,'s^-1'), n=1.22418, Ea=(82.4275,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7H;Y_rad_out;XH_out] for rate rule [R7H;CO_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=CC(C=CCC=O)O[O](1058)'],
    products = ['C=CC1[CH]C(CC=O)OO1(1079)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.94e+07,'s^-1'), n=0.93, Ea=(58.1576,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=CC(C=CCC=O)O[O](1058)'],
    products = ['O=CCC=CC1[CH]COO1(1080)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.07178e+07,'s^-1'), n=1.01592, Ea=(56.0365,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri_2H;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=CC(C=CCC=O)O[O](1058)'],
    products = ['C=CC1C=CC[CH]OOO1(1081)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(5.649e+12,'s^-1'), n=0.287, Ea=(253.77,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;carbonyl_intra_H;radadd_intra] for rate rule [R8_linear;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic
Ea raised from 248.1 to 253.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction26',
    reactants = ['CO(28)', 'C=CC(C=CC)O[O](1082)'],
    products = ['C=CC(C=CCC=O)O[O](1058)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(274200,'cm^3/(mol*s)'), n=2.53, Ea=(357.732,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO;C_pri] for rate rule [CO;C_pri/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=CC(C=CCC=O)O[O](1058)'],
    products = ['HO2(26)', 'C=CC=C=CCC=O(1001)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.63e+09,'s^-1'), n=1.11, Ea=(178.657,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R2OO_0H]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=CC(C=CCC=O)O[O](1058)'],
    products = ['HO2(26)', 'C=C=CC=CCC=O(976)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(3.63e+09,'s^-1'), n=1.11, Ea=(178.657,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R2OO_0H]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=COC=CC(C=C)O[O](1083)'],
    products = ['C=CC(C=CCC=O)O[O](1058)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 7 used for R_ROR;R1_doublebond_CH2;R2_doublebond_H;R_O_C
Exact match found for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=CC(C=CC=CO)O[O](1084)'],
    products = ['C=CC(C=CCC=O)O[O](1058)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_H;R_O_H]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction31',
    reactants = ['O(20)', 'C=CC([O])C=CCC=O(1085)'],
    products = ['C=CC(C=CCC=O)O[O](1058)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(54738.4,'m^3/(mol*s)'), n=0.884925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.9 to 0 kJ/mol."""),
)

network(
    label = '283',
    isomers = [
        'C=CC(C=CCC=O)O[O](1058)',
    ],
    reactants = [
        ('O2(2)', '[CH2]C=CC=CCC=O(956)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '283',
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

