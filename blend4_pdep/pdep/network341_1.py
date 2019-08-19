species(
    label = 'S(221)(220)',
    structure = SMILES('C=C[CH]C(C=O)O[O]'),
    E0 = (47.4434,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2782.5,750,1395,475,1775,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,571.413,804.445],'cm^-1')),
        HinderedRotor(inertia=(0.991122,'amu*angstrom^2'), symmetry=1, barrier=(22.7878,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.991039,'amu*angstrom^2'), symmetry=1, barrier=(22.7859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.99106,'amu*angstrom^2'), symmetry=1, barrier=(22.7864,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4151.01,'J/mol'), sigma=(6.69512,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=648.38 K, Pc=31.39 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.411623,0.069693,-5.96877e-05,2.53074e-08,-4.24461e-12,5842.95,30.7795], Tmin=(100,'K'), Tmax=(1432.96,'K')), NASAPolynomial(coeffs=[17.8423,0.0210371,-8.75595e-06,1.61229e-09,-1.10709e-13,847.415,-59.584], Tmin=(1432.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(47.4434,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJCO)"""),
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
    label = 'C5H6O(217)(216)',
    structure = SMILES('C=CC=CC=O'),
    E0 = (-29.5668,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.9508,'amu*angstrom^2'), symmetry=1, barrier=(21.8608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.95303,'amu*angstrom^2'), symmetry=1, barrier=(21.912,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3519.98,'J/mol'), sigma=(5.63814,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=549.81 K, Pc=44.56 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58677,0.0449815,-2.33922e-05,-1.17435e-10,2.4605e-12,-3462.46,19.7432], Tmin=(100,'K'), Tmax=(1159.43,'K')), NASAPolynomial(coeffs=[12.6828,0.020301,-9.05764e-06,1.75766e-09,-1.25367e-13,-6949.62,-39.3724], Tmin=(1159.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-29.5668,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'C5H6O(219)(218)',
    structure = SMILES('[CH2]C=CC=C[O]'),
    E0 = (108.621,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.59991,'amu*angstrom^2'), symmetry=1, barrier=(36.785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59845,'amu*angstrom^2'), symmetry=1, barrier=(36.7515,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3794.25,'J/mol'), sigma=(6.17562,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=592.65 K, Pc=36.55 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67303,0.0325845,3.83629e-05,-8.44772e-08,3.82244e-11,13165.2,20.9813], Tmin=(100,'K'), Tmax=(922.018,'K')), NASAPolynomial(coeffs=[18.3993,0.00812139,-9.24666e-08,-9.07048e-11,1.79206e-15,8036.21,-69.4428], Tmin=(922.018,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(108.621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CC=CCJ) + radical(C=COJ)"""),
)

species(
    label = 'O2(S)(162)(161)',
    structure = SMILES('O=O'),
    E0 = (85.6848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1487.4],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1857.18,'J/mol'), sigma=(4.34667,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=290.09 K, Pc=51.31 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53732,-0.00121571,5.31618e-06,-4.89443e-09,1.45845e-12,10304.5,4.68368], Tmin=(100,'K'), Tmax=(1074.56,'K')), NASAPolynomial(coeffs=[3.15382,0.00167804,-7.69971e-07,1.51275e-10,-1.08782e-14,10302.3,6.16754], Tmin=(1074.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(85.6848,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C=CC1C([O])C1O[O](5032)',
    structure = SMILES('C=CC1C([O])C1O[O]'),
    E0 = (223.252,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.85392,0.0565529,-1.95645e-05,-2.04283e-08,1.32133e-11,26975.6,27.5817], Tmin=(100,'K'), Tmax=(978.429,'K')), NASAPolynomial(coeffs=[17.1265,0.0191205,-6.77935e-06,1.2502e-09,-9.0806e-14,22398.7,-57.6856], Tmin=(978.429,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(223.252,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(ROOJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=C[CH]C1OOC1[O](5033)',
    structure = SMILES('C=C[CH]C1OOC1[O]'),
    E0 = (170.478,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49447,0.0608169,-4.46989e-05,1.68975e-08,-2.75916e-12,20588.7,21.4554], Tmin=(100,'K'), Tmax=(1333.14,'K')), NASAPolynomial(coeffs=[9.4254,0.0370205,-1.79238e-05,3.50794e-09,-2.48222e-13,18474.1,-19.0871], Tmin=(1333.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(170.478,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(12dioxetane) + radical(CCOJ) + radical(C=CCJCO)"""),
)

species(
    label = '[CH2]C1[CH]C(C=O)OO1(5034)',
    structure = SMILES('[CH2]C1[CH]C(C=O)OO1'),
    E0 = (136.9,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24967,0.0416582,3.43263e-05,-8.73735e-08,4.11877e-11,16582,30.0259], Tmin=(100,'K'), Tmax=(903.082,'K')), NASAPolynomial(coeffs=[19.1759,0.0124351,-4.77599e-07,-1.56135e-10,1.14705e-14,11298.1,-65.9588], Tmin=(903.082,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(136.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(12dioxolane) + radical(CCJCOOH) + radical(CJCOOH)"""),
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
    label = 'C=CC=C(C=O)O[O](5035)',
    structure = SMILES('C=CC=C(C=O)O[O]'),
    E0 = (40.3756,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(0.665654,'amu*angstrom^2'), symmetry=1, barrier=(15.3047,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.664335,'amu*angstrom^2'), symmetry=1, barrier=(15.2744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.66735,'amu*angstrom^2'), symmetry=1, barrier=(15.3437,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.913956,0.0714568,-7.96885e-05,4.67108e-08,-1.11397e-11,4964.26,25.6306], Tmin=(100,'K'), Tmax=(1007.03,'K')), NASAPolynomial(coeffs=[12.1799,0.0267082,-1.30348e-05,2.5857e-09,-1.85537e-13,2695.21,-28.8001], Tmin=(1007.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(40.3756,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + radical(ROOJ)"""),
)

species(
    label = 'C=C=CC(C=O)O[O](5036)',
    structure = SMILES('C=C=CC(C=O)O[O]'),
    E0 = (93.7509,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2782.5,750,1395,475,1775,1000,540,610,2055,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,180],'cm^-1')),
        HinderedRotor(inertia=(0.677665,'amu*angstrom^2'), symmetry=1, barrier=(15.5809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.677638,'amu*angstrom^2'), symmetry=1, barrier=(15.5802,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.677655,'amu*angstrom^2'), symmetry=1, barrier=(15.5806,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17481,0.0632921,-5.9566e-05,2.96031e-08,-6.05569e-12,11376.4,28.7805], Tmin=(100,'K'), Tmax=(1155.35,'K')), NASAPolynomial(coeffs=[11.6614,0.0269855,-1.24286e-05,2.40332e-09,-1.70021e-13,8953.26,-23.3254], Tmin=(1155.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(93.7509,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ)"""),
)

species(
    label = 'HCO(14)(15)',
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
    label = 'C=CC=CO[O](1797)',
    structure = SMILES('C=CC=CO[O]'),
    E0 = (161.628,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.819705,'amu*angstrom^2'), symmetry=1, barrier=(18.8466,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.820984,'amu*angstrom^2'), symmetry=1, barrier=(18.876,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57635,0.0449234,-2.28271e-05,-9.55747e-09,9.09409e-12,19534.3,21.3549], Tmin=(100,'K'), Tmax=(935.266,'K')), NASAPolynomial(coeffs=[14.3024,0.0115713,-3.13734e-06,5.01253e-10,-3.50008e-14,16232.1,-44.1169], Tmin=(935.266,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(161.628,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ)"""),
)

species(
    label = 'C=CC[C](C=O)O[O](5037)',
    structure = SMILES('C=CCC(=C[O])O[O]'),
    E0 = (83.3607,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,1437.5,1250,1305,750,350,306.594,307.699,308.432],'cm^-1')),
        HinderedRotor(inertia=(0.173472,'amu*angstrom^2'), symmetry=1, barrier=(11.7135,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173692,'amu*angstrom^2'), symmetry=1, barrier=(11.7093,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.176341,'amu*angstrom^2'), symmetry=1, barrier=(11.7256,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.507787,0.0682024,-5.68025e-05,1.89581e-08,-7.52708e-13,10159.2,31.4865], Tmin=(100,'K'), Tmax=(1016.23,'K')), NASAPolynomial(coeffs=[16.9693,0.0189072,-6.91856e-06,1.24166e-09,-8.65466e-14,6013.18,-52.1337], Tmin=(1016.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(83.3607,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=COJ)"""),
)

species(
    label = 'C=[C]CC(C=O)O[O](5038)',
    structure = SMILES('C=[C]CC(C=O)O[O]'),
    E0 = (168.368,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2782.5,750,1395,475,1775,1000,2750,2850,1437.5,1250,1305,750,350,1685,370,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,180,436.687],'cm^-1')),
        HinderedRotor(inertia=(0.0409036,'amu*angstrom^2'), symmetry=1, barrier=(12.4803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158547,'amu*angstrom^2'), symmetry=1, barrier=(3.64531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.542805,'amu*angstrom^2'), symmetry=1, barrier=(12.4802,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0409068,'amu*angstrom^2'), symmetry=1, barrier=(12.481,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08552,0.0647403,-5.79421e-05,2.7523e-08,-5.40827e-12,20354.4,31.6256], Tmin=(100,'K'), Tmax=(1195.43,'K')), NASAPolynomial(coeffs=[11.7705,0.0289874,-1.308e-05,2.50422e-09,-1.76083e-13,17799.8,-21.8305], Tmin=(1195.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(168.368,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(ROOJ)"""),
)

species(
    label = 'C=C[CH][C](C=O)OO(5039)',
    structure = SMILES('[CH2]C=CC(=C[O])OO'),
    E0 = (19.2233,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.24957,0.0717151,-5.28315e-05,7.27808e-09,5.15717e-12,2456.77,29.3752], Tmin=(100,'K'), Tmax=(963.047,'K')), NASAPolynomial(coeffs=[19.3787,0.0167924,-5.49285e-06,9.56321e-10,-6.75213e-14,-2365.19,-68.0972], Tmin=(963.047,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(19.2233,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CC=CCJ) + radical(C=COJ)"""),
)

species(
    label = 'C=CCC([C]=O)O[O](5040)',
    structure = SMILES('C=CCC([C]=O)O[O]'),
    E0 = (90.4873,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,1855,455,950,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,426.286,426.289],'cm^-1')),
        HinderedRotor(inertia=(0.0834857,'amu*angstrom^2'), symmetry=1, barrier=(10.7657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0834844,'amu*angstrom^2'), symmetry=1, barrier=(10.7656,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0834856,'amu*angstrom^2'), symmetry=1, barrier=(10.7656,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0834865,'amu*angstrom^2'), symmetry=1, barrier=(10.7656,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.749865,0.0665799,-6.02635e-05,2.80209e-08,-5.21594e-12,11004.2,33.1754], Tmin=(100,'K'), Tmax=(1289.53,'K')), NASAPolynomial(coeffs=[14.8424,0.0228657,-9.41417e-06,1.73235e-09,-1.19377e-13,7369.63,-38.3965], Tmin=(1289.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(90.4873,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]=CCC(C=O)O[O](5041)',
    structure = SMILES('[CH]=CCC(C=O)O[O]'),
    E0 = (177.623,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2782.5,750,1395,475,1775,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.0250002,'amu*angstrom^2'), symmetry=1, barrier=(14.5466,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.632595,'amu*angstrom^2'), symmetry=1, barrier=(14.5446,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169707,'amu*angstrom^2'), symmetry=1, barrier=(3.90189,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.632633,'amu*angstrom^2'), symmetry=1, barrier=(14.5455,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.841943,0.066049,-5.80848e-05,2.61976e-08,-4.77247e-12,21479.5,32.4287], Tmin=(100,'K'), Tmax=(1304.71,'K')), NASAPolynomial(coeffs=[14.2623,0.0249046,-1.07819e-05,2.0273e-09,-1.4112e-13,17977.6,-35.8865], Tmin=(1304.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(177.623,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(ROOJ)"""),
)

species(
    label = 'C=C[CH]C([C]=O)OO(5042)',
    structure = SMILES('C=C[CH]C([C]=O)OO'),
    E0 = (55.3993,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1855,455,950,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,337.334,337.334],'cm^-1')),
        HinderedRotor(inertia=(0.322074,'amu*angstrom^2'), symmetry=1, barrier=(26.0073,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.322067,'amu*angstrom^2'), symmetry=1, barrier=(26.0073,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.322067,'amu*angstrom^2'), symmetry=1, barrier=(26.0072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.322057,'amu*angstrom^2'), symmetry=1, barrier=(26.0072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48142,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0303874,0.0767888,-7.19783e-05,3.28858e-08,-5.86407e-12,6814.71,32.2065], Tmin=(100,'K'), Tmax=(1365.52,'K')), NASAPolynomial(coeffs=[20.0864,0.0180383,-7.44116e-06,1.37749e-09,-9.54639e-14,1337.38,-70.8004], Tmin=(1365.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(55.3993,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O) + radical(C=CCJCO)"""),
)

species(
    label = 'C=[C][CH]C(C=O)OO(5043)',
    structure = SMILES('C=[C][CH]C(C=O)OO'),
    E0 = (133.28,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,3025,407.5,1350,352.5,1685,370,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,281.951,411.685],'cm^-1')),
        HinderedRotor(inertia=(0.218669,'amu*angstrom^2'), symmetry=1, barrier=(26.2992,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218668,'amu*angstrom^2'), symmetry=1, barrier=(26.2992,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218668,'amu*angstrom^2'), symmetry=1, barrier=(26.2992,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218668,'amu*angstrom^2'), symmetry=1, barrier=(26.2992,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.336845,'amu*angstrom^2'), symmetry=1, barrier=(40.512,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.415828,0.0744591,-6.83249e-05,3.10823e-08,-5.63803e-12,16162.5,30.4706], Tmin=(100,'K'), Tmax=(1318.85,'K')), NASAPolynomial(coeffs=[17.0172,0.0241082,-1.10582e-05,2.13462e-09,-1.50749e-13,11783.5,-54.2161], Tmin=(1318.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(133.28,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=CCJCO) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C[CH]C(C=O)OO(5044)',
    structure = SMILES('[CH]C=CC(C=O)OO'),
    E0 = (135.826,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.684123,0.069361,-5.28634e-05,2.0291e-08,-3.19002e-12,16458.1,31.8359], Tmin=(100,'K'), Tmax=(1468.92,'K')), NASAPolynomial(coeffs=[14.7023,0.0311882,-1.38829e-05,2.59979e-09,-1.79094e-13,12339.7,-41.1843], Tmin=(1468.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(135.826,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[O]OC(C=O)C1[CH]C1(5045)',
    structure = SMILES('[O]OC(C=O)C1[CH]C1'),
    E0 = (182.212,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26869,0.0488613,-1.0735e-05,-1.99051e-08,1.04884e-11,22023.1,31.6418], Tmin=(100,'K'), Tmax=(1037.73,'K')), NASAPolynomial(coeffs=[13.859,0.0243375,-9.98734e-06,1.90747e-09,-1.37141e-13,18117.4,-35.7935], Tmin=(1037.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(182.212,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(cyclopropane) + radical(ROOJ)"""),
)

species(
    label = 'C=CC1O[CH]C1O[O](5046)',
    structure = SMILES('C=CC1O[CH]C1O[O]'),
    E0 = (202.418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.425952,0.0733366,-6.89348e-05,3.30114e-08,-5.88591e-12,24525.4,31.5293], Tmin=(100,'K'), Tmax=(1648.76,'K')), NASAPolynomial(coeffs=[17.0703,0.01505,-1.49671e-06,-8.3502e-11,1.57497e-14,20908.9,-55.0995], Tmin=(1648.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(202.418,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Oxetane) + radical(CCsJOCs) + radical(ROOJ)"""),
)

species(
    label = 'C=C[CH]C1[CH]OOO1(5047)',
    structure = SMILES('C=C[CH]C1[CH]OOO1'),
    E0 = (298.821,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.658602,0.0595687,-1.11459e-05,-4.14738e-08,2.50628e-11,36073.3,23.4592], Tmin=(100,'K'), Tmax=(895.527,'K')), NASAPolynomial(coeffs=[19.1104,0.0157618,-2.44251e-06,1.92142e-10,-9.3593e-15,31220.2,-72.1684], Tmin=(895.527,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(298.821,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(123trioxolane) + radical(C=CCJCO) + radical(CCsJOO)"""),
)

species(
    label = 'O=CC1[CH][CH]COO1(5048)',
    structure = SMILES('O=CC1[CH][CH]COO1'),
    E0 = (127.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48639,0.0365197,4.46389e-05,-9.36996e-08,4.21901e-11,15471.7,28.6285], Tmin=(100,'K'), Tmax=(909.215,'K')), NASAPolynomial(coeffs=[17.4585,0.0155384,-2.05798e-06,1.60271e-10,-1.11752e-14,10530.1,-58.1103], Tmin=(909.215,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(127.741,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(12dioxane) + radical(CCJCOOH) + radical(CCJCOOH)"""),
)

species(
    label = 'C=CC=C(C=O)OO(5049)',
    structure = SMILES('C=CC=C(C=O)OO'),
    E0 = (-111.629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.645602,0.0751906,-7.58062e-05,3.87473e-08,-8.00312e-12,-13306.3,25.8277], Tmin=(100,'K'), Tmax=(1155.57,'K')), NASAPolynomial(coeffs=[14.5787,0.0269609,-1.32009e-05,2.62913e-09,-1.89144e-13,-16526.4,-43.4062], Tmin=(1155.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-111.629,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C=CC(C=O)OO(5050)',
    structure = SMILES('C=C=CC(C=O)OO'),
    E0 = (-58.2538,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.793123,0.0682854,-5.97117e-05,2.63277e-08,-4.69664e-12,-6889.2,29.3889], Tmin=(100,'K'), Tmax=(1324,'K')), NASAPolynomial(coeffs=[14.723,0.0262011,-1.20333e-05,2.32048e-09,-1.63549e-13,-10577.8,-41.7243], Tmin=(1324,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-58.2538,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'CO(10)(11)',
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
    label = 'C=C[CH]CO[O](2752)',
    structure = SMILES('C=C[CH]CO[O]'),
    E0 = (164.077,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,222.742,224.411],'cm^-1')),
        HinderedRotor(inertia=(0.865045,'amu*angstrom^2'), symmetry=1, barrier=(29.9089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.45475,'amu*angstrom^2'), symmetry=1, barrier=(16.0987,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.09715,'amu*angstrom^2'), symmetry=1, barrier=(74.4978,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52178,0.0546104,-4.85103e-05,2.32517e-08,-4.62538e-12,19823,20.3981], Tmin=(100,'K'), Tmax=(1182.16,'K')), NASAPolynomial(coeffs=[10.2752,0.024992,-1.09285e-05,2.05781e-09,-1.43343e-13,17753.4,-23.2968], Tmin=(1182.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(164.077,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CCJCO) + radical(ROOJ)"""),
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
    label = 'C=CC1OC1C=O(5051)',
    structure = SMILES('C=CC1OC1C=O'),
    E0 = (-116.399,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80385,0.0294689,5.26188e-05,-1.00687e-07,4.49928e-11,-13902.8,23.9649], Tmin=(100,'K'), Tmax=(900.611,'K')), NASAPolynomial(coeffs=[17.2226,0.0109598,2.16104e-07,-2.86516e-10,2.04043e-14,-18706.6,-60.0588], Tmin=(900.611,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-116.399,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide)"""),
)

species(
    label = 'C=CC([CH]O[O])C=O(5052)',
    structure = SMILES('C=CC([CH]O[O])C=O'),
    E0 = (122.268,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2782.5,750,1395,475,1775,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.204442,'amu*angstrom^2'), symmetry=1, barrier=(4.70053,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204909,'amu*angstrom^2'), symmetry=1, barrier=(4.71125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204621,'amu*angstrom^2'), symmetry=1, barrier=(4.70464,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.573983,0.0839723,-0.000129012,1.15116e-07,-4.06823e-11,14820.5,30.2182], Tmin=(100,'K'), Tmax=(832.968,'K')), NASAPolynomial(coeffs=[6.67123,0.0386397,-1.8469e-05,3.50692e-09,-2.40758e-13,14361.6,5.25977], Tmin=(832.968,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(122.268,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(CCsJOOH)"""),
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
    label = 'C=CC=C[C]=O(4549)',
    structure = SMILES('C=C[CH]C=C=O'),
    E0 = (103.397,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,2120,512.5,787.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,292.831,292.831],'cm^-1')),
        HinderedRotor(inertia=(0.484811,'amu*angstrom^2'), symmetry=1, barrier=(29.5007,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.48481,'amu*angstrom^2'), symmetry=1, barrier=(29.5007,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78546,0.0413554,-1.91973e-05,-5.00058e-09,4.9111e-12,12521.9,18.5543], Tmin=(100,'K'), Tmax=(1025.01,'K')), NASAPolynomial(coeffs=[11.5734,0.0180947,-7.01485e-06,1.29194e-09,-9.09454e-14,9730.71,-32.7362], Tmin=(1025.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(103.397,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-Cd(CCO)HH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + radical(C=CCJC=C=O)"""),
)

species(
    label = 'C=CC1OOC1C=O(5053)',
    structure = SMILES('C=CC1OOC1C=O'),
    E0 = (-65.1023,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31487,0.0452958,4.67689e-06,-3.734e-08,1.66497e-11,-7721.35,25.9413], Tmin=(100,'K'), Tmax=(1027.5,'K')), NASAPolynomial(coeffs=[14.7988,0.0241486,-1.021e-05,2.00813e-09,-1.47618e-13,-12146.9,-47.5283], Tmin=(1027.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-65.1023,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(12dioxetane)"""),
)

species(
    label = 'C=CC=C[CH]OO[O](4882)',
    structure = SMILES('[CH2]C=CC=COO[O]'),
    E0 = (279.635,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.817945,0.0555014,-1.01637e-05,-3.79923e-08,2.21536e-11,33760.3,31.5561], Tmin=(100,'K'), Tmax=(925.081,'K')), NASAPolynomial(coeffs=[19.1676,0.0135546,-2.78526e-06,3.89133e-10,-2.83037e-14,28765.2,-64.1902], Tmin=(925.081,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CC=CCJ) + radical(ROOJ)"""),
)

species(
    label = 'C=C[CH]C(=CO)O[O](5054)',
    structure = SMILES('[CH2]C=CC(=CO)O[O]'),
    E0 = (29.7654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3000,3100,440,815,1455,1000,350,440,435,1725,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1277.5,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.925843,'amu*angstrom^2'), symmetry=1, barrier=(21.2869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.926163,'amu*angstrom^2'), symmetry=1, barrier=(21.2943,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.92584,'amu*angstrom^2'), symmetry=1, barrier=(21.2869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.925829,'amu*angstrom^2'), symmetry=1, barrier=(21.2866,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.556183,0.0828208,-8.78495e-05,4.5187e-08,-8.70714e-12,3758.88,31.6178], Tmin=(100,'K'), Tmax=(1449.19,'K')), NASAPolynomial(coeffs=[21.7621,0.0100903,-1.05035e-06,-4.19022e-11,9.42282e-15,-1541.25,-80.3035], Tmin=(1449.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(29.7654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CC=CCJ) + radical(ROOJ)"""),
)

species(
    label = 'C=C[CH]C([O])C=O(5055)',
    structure = SMILES('C=C[CH]C([O])C=O'),
    E0 = (67.6668,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,180,882.026,884.321],'cm^-1')),
        HinderedRotor(inertia=(0.0406242,'amu*angstrom^2'), symmetry=1, barrier=(22.4166,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.975059,'amu*angstrom^2'), symmetry=1, barrier=(22.4185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.35286,'amu*angstrom^2'), symmetry=1, barrier=(22.4541,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1591,0.0524749,-2.51031e-05,-6.23621e-09,6.20834e-12,8249.38,27.3465], Tmin=(100,'K'), Tmax=(1040.75,'K')), NASAPolynomial(coeffs=[14.5644,0.020963,-8.52556e-06,1.61832e-09,-1.15965e-13,4375.37,-43.0681], Tmin=(1040.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(67.6668,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(C=CCJCO)"""),
)

species(
    label = '[CH]=C[CH2](924)',
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
    label = '[O]O[CH]C=O(4488)',
    structure = SMILES('[O]O[CH]C=O'),
    E0 = (38.8345,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.304071,'amu*angstrom^2'), symmetry=1, barrier=(6.9912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17371,'amu*angstrom^2'), symmetry=1, barrier=(56.7408,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0355,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.48924,0.0359331,-5.02567e-05,4.15061e-08,-1.40856e-11,4722.53,16.7287], Tmin=(100,'K'), Tmax=(782.809,'K')), NASAPolynomial(coeffs=[5.94539,0.0159797,-7.6283e-06,1.46013e-09,-1.01349e-13,4251.69,1.34978], Tmin=(782.809,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(38.8345,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + radical(OCJC=O) + radical(ROOJ)"""),
)

species(
    label = '[CH2][CH]C1OOC1C=O(5056)',
    structure = SMILES('[CH2][CH]C1OOC1C=O'),
    E0 = (213.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02218,0.0588185,-4.17082e-05,1.44395e-08,-2.00646e-12,25744.1,29.7136], Tmin=(100,'K'), Tmax=(1670.62,'K')), NASAPolynomial(coeffs=[15.6574,0.0237774,-1.02462e-05,1.88464e-09,-1.27703e-13,20854.1,-48.4039], Tmin=(1670.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(213.113,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(12dioxetane) + radical(RCCJ) + radical(CCJCOOH)"""),
)

species(
    label = '[O]OC1C=CCC1[O](5057)',
    structure = SMILES('[O]OC1C=CCC1[O]'),
    E0 = (130.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24032,0.0451966,1.09749e-05,-4.90786e-08,2.23223e-11,15819.6,26.1385], Tmin=(100,'K'), Tmax=(985.935,'K')), NASAPolynomial(coeffs=[16.3577,0.0206984,-7.79324e-06,1.5046e-09,-1.12217e-13,11048.3,-55.6592], Tmin=(985.935,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(130.59,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(CC(C)OJ) + radical(ROOJ)"""),
)

species(
    label = 'C[C]=CC(C=O)O[O](5058)',
    structure = SMILES('C[C]=CC(C=O)O[O]'),
    E0 = (154.988,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17019,0.0672584,-7.10642e-05,4.4083e-08,-1.17166e-11,18738.4,30.2313], Tmin=(100,'K'), Tmax=(888.956,'K')), NASAPolynomial(coeffs=[8.34393,0.0349794,-1.6598e-05,3.23683e-09,-2.29643e-13,17463,-3.53346], Tmin=(888.956,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(154.988,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = 'CC=[C]C(C=O)O[O](5059)',
    structure = SMILES('CC=[C]C(C=O)O[O]'),
    E0 = (154.988,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2782.5,750,1395,475,1775,1000,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,294.577],'cm^-1')),
        HinderedRotor(inertia=(0.151114,'amu*angstrom^2'), symmetry=1, barrier=(9.30525,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151114,'amu*angstrom^2'), symmetry=1, barrier=(9.30524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151113,'amu*angstrom^2'), symmetry=1, barrier=(9.30525,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151113,'amu*angstrom^2'), symmetry=1, barrier=(9.30525,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17019,0.0672584,-7.10643e-05,4.40831e-08,-1.17167e-11,18738.4,30.2313], Tmin=(100,'K'), Tmax=(888.95,'K')), NASAPolynomial(coeffs=[8.34392,0.0349794,-1.6598e-05,3.23683e-09,-2.29643e-13,17463,-3.53343], Tmin=(888.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(154.988,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Cds_S) + radical(ROOJ)"""),
)

species(
    label = '[CH2]C=[C]C(C=O)OO(5060)',
    structure = SMILES('[CH2]C=[C]C(C=O)OO'),
    E0 = (154.483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.80357,0.0697735,-6.35533e-05,2.96037e-08,-5.61163e-12,18695.4,31.0524], Tmin=(100,'K'), Tmax=(1247.17,'K')), NASAPolynomial(coeffs=[13.9537,0.0275973,-1.28267e-05,2.48785e-09,-1.76118e-13,15415.3,-35.2938], Tmin=(1247.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(154.483,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'CC=C[C](C=O)O[O](5061)',
    structure = SMILES('CC=CC(=C[O])O[O]'),
    E0 = (53.1722,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.698463,'amu*angstrom^2'), symmetry=1, barrier=(16.059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.694287,'amu*angstrom^2'), symmetry=1, barrier=(15.963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.695516,'amu*angstrom^2'), symmetry=1, barrier=(15.9913,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.202975,0.0773826,-8.20211e-05,4.40911e-08,-9.25946e-12,6537.09,28.7436], Tmin=(100,'K'), Tmax=(1168.78,'K')), NASAPolynomial(coeffs=[17.3359,0.0187469,-6.76811e-06,1.16682e-09,-7.79687e-14,2532.2,-56.5848], Tmin=(1168.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(53.1722,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(ROOJ)"""),
)

species(
    label = 'CC=CC([C]=O)O[O](5062)',
    structure = SMILES('CC=CC([C]=O)O[O]'),
    E0 = (77.1072,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1855,455,950,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,242.046],'cm^-1')),
        HinderedRotor(inertia=(0.234335,'amu*angstrom^2'), symmetry=1, barrier=(9.51132,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.234953,'amu*angstrom^2'), symmetry=1, barrier=(9.53032,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.231811,'amu*angstrom^2'), symmetry=1, barrier=(9.51347,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.229473,'amu*angstrom^2'), symmetry=1, barrier=(9.52989,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06357,0.0665647,-6.5268e-05,3.498e-08,-7.76933e-12,9377.97,30.9481], Tmin=(100,'K'), Tmax=(1068.48,'K')), NASAPolynomial(coeffs=[10.9578,0.0295247,-1.32694e-05,2.53631e-09,-1.7833e-13,7263.59,-17.4412], Tmin=(1068.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(77.1072,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(ROOJ) + radical(CCCJ=O)"""),
)

species(
    label = '[O]OC1[CH]OCC=C1(5063)',
    structure = SMILES('[O]OC1[CH]OCC=C1'),
    E0 = (110.913,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00368,0.0537757,-1.06815e-05,-3.1875e-08,1.8711e-11,13458.8,23.7551], Tmin=(100,'K'), Tmax=(923.949,'K')), NASAPolynomial(coeffs=[16.2626,0.0190321,-5.11659e-06,7.92922e-10,-5.38701e-14,9302.46,-55.8868], Tmin=(923.949,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(110.913,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(36dihydro2hpyran) + radical(CCsJOCs) + radical(ROOJ)"""),
)

species(
    label = 'O=CC1C=CCO1(5064)',
    structure = SMILES('O=CC1C=CCO1'),
    E0 = (-202.911,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.18695,0.0213917,6.18868e-05,-9.45278e-08,3.69486e-11,-24322.7,21.3409], Tmin=(100,'K'), Tmax=(975.366,'K')), NASAPolynomial(coeffs=[13.7287,0.020282,-7.49284e-06,1.48127e-09,-1.13629e-13,-28772.9,-45.3248], Tmin=(975.366,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-202.911,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(25dihydrofuran)"""),
)

species(
    label = 'C=C[C]=CC=O(4714)',
    structure = SMILES('C=C[C]=CC=O'),
    E0 = (169.429,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.18163,'amu*angstrom^2'), symmetry=1, barrier=(27.168,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17982,'amu*angstrom^2'), symmetry=1, barrier=(27.1265,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83685,0.0478006,-4.01103e-05,1.75891e-08,-3.19807e-12,20455.1,18.4329], Tmin=(100,'K'), Tmax=(1277.85,'K')), NASAPolynomial(coeffs=[10.001,0.0222446,-1.01115e-05,1.93838e-09,-1.36142e-13,18368.6,-22.9561], Tmin=(1277.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(169.429,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CJC=C)"""),
)

species(
    label = 'O=CC1C=CCOO1(5065)',
    structure = SMILES('O=CC1C=CCOO1'),
    E0 = (-154.867,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72814,0.0287652,6.00928e-05,-1.0008e-07,4.01538e-11,-18525.3,22.5951], Tmin=(100,'K'), Tmax=(974.812,'K')), NASAPolynomial(coeffs=[16.7477,0.0198441,-7.28708e-06,1.46967e-09,-1.15028e-13,-23958,-62.328], Tmin=(974.812,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-154.867,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(36dihydro12dioxin)"""),
)

species(
    label = '[CH2]C=COC=CO[O](5066)',
    structure = SMILES('C=C[CH]OC=CO[O]'),
    E0 = (102.698,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3025,407.5,1350,352.5,310.899,311.026,311.105,311.263],'cm^-1')),
        HinderedRotor(inertia=(0.28455,'amu*angstrom^2'), symmetry=1, barrier=(19.5528,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284889,'amu*angstrom^2'), symmetry=1, barrier=(19.5511,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284757,'amu*angstrom^2'), symmetry=1, barrier=(19.5531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284938,'amu*angstrom^2'), symmetry=1, barrier=(19.5515,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.320819,0.0700771,-5.08699e-05,5.51535e-09,5.82723e-12,12494,30.5873], Tmin=(100,'K'), Tmax=(959.34,'K')), NASAPolynomial(coeffs=[19.3959,0.0156251,-4.94804e-06,8.56638e-10,-6.08562e-14,7679.88,-66.6629], Tmin=(959.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(102.698,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJ(O)C)"""),
)

species(
    label = 'CH2(17)(18)',
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
    label = '[CH]=CC(C=O)O[O](5067)',
    structure = SMILES('[CH]=CC(C=O)O[O]'),
    E0 = (200.268,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3120,650,792.5,1650,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.664291,'amu*angstrom^2'), symmetry=1, barrier=(15.2734,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.126798,'amu*angstrom^2'), symmetry=1, barrier=(6.75815,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.664893,'amu*angstrom^2'), symmetry=1, barrier=(15.2872,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61891,0.0526661,-5.04425e-05,2.52312e-08,-5.14705e-12,24172.2,26.7287], Tmin=(100,'K'), Tmax=(1165.5,'K')), NASAPolynomial(coeffs=[10.8757,0.0208966,-9.55477e-06,1.84328e-09,-1.30301e-13,22014.4,-19.3475], Tmin=(1165.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(200.268,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_P)"""),
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
    E0 = (47.4434,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (223.252,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (170.478,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (146.576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (261.622,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (321.359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (288.037,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (206.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (326.048,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (212.237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (215.031,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (322.808,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (139.288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (284.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (287.347,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (99.9996,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (273.379,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (202.418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (298.821,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (127.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (125.69,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (72.4166,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (56.118,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (350.165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (177.147,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (245.018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (186.813,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (55.7277,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (593.435,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (173.632,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (310.672,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (415.488,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (213.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (130.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (255.175,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (316.909,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (193.046,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (202.959,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (166.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (146.477,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (109.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (226.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (54.9746,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (416.498,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (581.348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O2(2)(2)', 'C5H6O(217)(216)'],
    products = ['S(221)(220)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.233e+09,'cm^3/(mol*s)'), n=1.533, Ea=(85.632,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-COH_Cds;YJ] for rate rule [Cds-COH_Cds;O2b]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 81.3 to 85.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['S(221)(220)'],
    products = ['C=CC1C([O])C1O[O](5032)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(175.809,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_csHCd] for rate rule [R4_S_CO;carbonylbond_intra_H;radadd_intra_csHCd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 174.2 to 175.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['S(221)(220)'],
    products = ['C=C[CH]C1OOC1[O](5033)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(123.035,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_O] for rate rule [R5_SS_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 123.0 to 123.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['S(221)(220)'],
    products = ['[CH2]C1[CH]C(C=O)OO1(5034)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.82166e+09,'s^-1'), n=0.527281, Ea=(99.1325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_2H_pri;radadd_intra_O] + [R6;doublebond_intra_2H_pri;radadd_intra] for rate rule [R6;doublebond_intra_2H_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)(3)', 'C=CC=C(C=O)O[O](5035)'],
    products = ['S(221)(220)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(87.5179,'m^3/(mol*s)'), n=1.66467, Ea=(9.4546,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;HJ] for rate rule [Cds-COOs_Cds;HJ]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)(3)', 'C=C=CC(C=O)O[O](5036)'],
    products = ['S(221)(220)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(5.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(15.8155,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2714 used for Ca_Cds-CsH;HJ
Exact match found for rate rule [Ca_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['HCO(14)(15)', 'C=CC=CO[O](1797)'],
    products = ['S(221)(220)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(5.2e+11,'cm^3/(mol*s)'), n=0, Ea=(93.9308,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Cd_R;CO_pri_rad] for rate rule [Cds-OsH_Cds;CO_pri_rad]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=CC[C](C=O)O[O](5037)'],
    products = ['S(221)(220)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(11288.6,'s^-1'), n=2.7035, Ea=(123.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_noH;Cs_H_out_H/Cd] for rate rule [R2H_S;C_rad_out_OneDe/O;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=[C]CC(C=O)O[O](5038)'],
    products = ['S(221)(220)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.66329e+10,'s^-1'), n=0.993, Ea=(157.679,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['S(221)(220)'],
    products = ['C=C[CH][C](C=O)OO(5039)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.83109e+08,'s^-1'), n=1.32333, Ea=(164.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_O;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=CCC([C]=O)O[O](5040)'],
    products = ['S(221)(220)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.29711e+07,'s^-1'), n=1.52333, Ea=(124.544,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R3H_SS_Cs;CO_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=CCC(C=O)O[O](5041)'],
    products = ['S(221)(220)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 194 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C[CH]C([C]=O)OO(5042)'],
    products = ['S(221)(220)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;Y_rad_out;O_H_out] for rate rule [R4H_SSS;CO_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=[C][CH]C(C=O)OO(5043)'],
    products = ['S(221)(220)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.81182e+08,'s^-1'), n=1.25566, Ea=(151.659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;XH_out] for rate rule [R5HJ_1;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=C[CH]C(C=O)OO(5044)'],
    products = ['S(221)(220)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R6HJ_2;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['O2(2)(2)', 'C5H6O(219)(218)'],
    products = ['S(221)(220)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(6.34662e+06,'m^3/(mol*s)'), n=-0.0521589, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [C_sec_rad;O2_birad] + [C_rad/H/TwoDe;Y_rad] for rate rule [C_rad/H/TwoDe;O2_birad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -5.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['S(221)(220)'],
    products = ['[O]OC(C=O)C1[CH]C1(5045)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_pri_2H;radadd_intra_csHCs]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['S(221)(220)'],
    products = ['C=CC1O[CH]C1O[O](5046)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(6.89861e+07,'s^-1'), n=1.13751, Ea=(154.975,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs] for rate rule [R4_S_CO;carbonyl_intra_H;radadd_intra_csHCd]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic
Ea raised from 154.0 to 155.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction19',
    reactants = ['S(221)(220)'],
    products = ['C=C[CH]C1[CH]OOO1(5047)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.81184e+09,'s^-1'), n=0.551229, Ea=(251.377,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra] for rate rule [R5_SS_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic
Ea raised from 250.8 to 251.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction20',
    reactants = ['S(221)(220)'],
    products = ['O=CC1[CH][CH]COO1(5048)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(9.91671e+09,'s^-1'), n=0.30082, Ea=(80.2973,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra_pri_2H;radadd_intra] for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 76.4 to 80.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction21',
    reactants = ['S(221)(220)'],
    products = ['C=CC=C(C=O)OO(5049)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['S(221)(220)'],
    products = ['C=C=CC(C=O)OO(5050)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['S(221)(220)'],
    products = ['O2(S)(162)(161)', 'C5H6O(217)(216)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(8.67463,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 8.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction24',
    reactants = ['CO(10)(11)', 'C=C[CH]CO[O](2752)'],
    products = ['S(221)(220)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.118397,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction25',
    reactants = ['S(221)(220)'],
    products = ['O(4)(4)', 'C=CC1OC1C=O(5051)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.27e+09,'s^-1'), n=1.06, Ea=(129.704,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2OO_S;C_sec_rad_intra;OOJ] for rate rule [R2OO_S;C_rad/H/OneDe_intra;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=CC([CH]O[O])C=O(5052)'],
    products = ['S(221)(220)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(8.889e+11,'s^-1'), n=0.232, Ea=(122.75,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ;CO]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction27',
    reactants = ['S(221)(220)'],
    products = ['HO2(8)(9)', 'C=CC=C[C]=O(4549)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(9.58174e+10,'s^-1'), n=0.573333, Ea=(139.369,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction28',
    reactants = ['S(221)(220)'],
    products = ['C=CC1OOC1C=O(5053)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;O_rad;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=CC=C[CH]OO[O](4882)'],
    products = ['S(221)(220)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_R]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=C[CH]C(=CO)O[O](5054)'],
    products = ['S(221)(220)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction31',
    reactants = ['O(4)(4)', 'C=C[CH]C([O])C=O(5055)'],
    products = ['S(221)(220)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(54738.4,'m^3/(mol*s)'), n=0.884925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]=C[CH2](924)', '[O]O[CH]C=O(4488)'],
    products = ['S(221)(220)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.23625e+06,'m^3/(mol*s)'), n=0.36814, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction33',
    reactants = ['S(221)(220)'],
    products = ['[CH2][CH]C1OOC1C=O(5056)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(165.67,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 164.1 to 165.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction34',
    reactants = ['S(221)(220)'],
    products = ['[O]OC1C=CCC1[O](5057)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(3.49749e+08,'s^-1'), n=0.656505, Ea=(83.1471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SMS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R6_SMS_CO;carbonylbond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 80.0 to 83.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction35',
    reactants = ['S(221)(220)'],
    products = ['C[C]=CC(C=O)O[O](5058)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['CC=[C]C(C=O)O[O](5059)'],
    products = ['S(221)(220)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['S(221)(220)'],
    products = ['[CH2]C=[C]C(C=O)OO(5060)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(274,'s^-1'), n=3.09, Ea=(145.603,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 337 used for R4H_SSS_OCs;O_rad_out;Cd_H_out_doubleC
Exact match found for rate rule [R4H_SSS_OCs;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['CC=C[C](C=O)O[O](5061)'],
    products = ['S(221)(220)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(800000,'s^-1'), n=1.81, Ea=(149.787,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;C_rad_out_single;Cs_H_out_2H] for rate rule [R4H_SDS;C_rad_out_OneDe/O;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['CC=CC([C]=O)O[O](5062)'],
    products = ['S(221)(220)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.45388e+06,'s^-1'), n=1.705, Ea=(89.2238,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSMS;Y_rad_out;Cs_H_out_2H] for rate rule [R5H_SSMS;CO_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['S(221)(220)'],
    products = ['[O]OC1[CH]OCC=C1(5063)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(2.39072e+10,'s^-1'), n=0.346137, Ea=(99.0339,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SMS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R6_SMS_CO;carbonyl_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction41',
    reactants = ['S(221)(220)'],
    products = ['O(4)(4)', 'O=CC1C=CCO1(5064)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(5.13e+10,'s^-1'), n=0, Ea=(61.9232,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4OO;C_pri_rad_intra;OO] for rate rule [R4OO_SDS;C_pri_rad_intra;OOJ]
Euclidian distance = 1.41421356237
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction42',
    reactants = ['S(221)(220)'],
    products = ['HO2(8)(9)', 'C=C[C]=CC=O(4714)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(3.63e+09,'s^-1'), n=1.11, Ea=(178.657,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R2OO_0H]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction43',
    reactants = ['S(221)(220)'],
    products = ['O=CC1C=CCOO1(5065)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), Tmin=(550,'K'), Tmax=(650,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R6_SSSDS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2]C=COC=CO[O](5066)'],
    products = ['S(221)(220)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_H;R_O_C]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction45',
    reactants = ['CH2(17)(18)', '[CH]=CC(C=O)O[O](5067)'],
    products = ['S(221)(220)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(2.23625e+06,'m^3/(mol*s)'), n=0.36814, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

network(
    label = '341',
    isomers = [
        'S(221)(220)',
    ],
    reactants = [
        ('O2(2)(2)', 'C5H6O(217)(216)'),
        ('O2(2)(2)', 'C5H6O(219)(218)'),
        ('O2(S)(162)(161)', 'C5H6O(217)(216)'),
    ],
    bathGas = {
        'Ne': 0.333333,
        'N2': 0.333333,
        'Ar(8)': 0.333333,
    },
)

pressureDependence(
    label = '341',
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

