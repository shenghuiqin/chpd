species(
    label = 'S(459)(458)',
    structure = SMILES('C=CCC(=O)C=CC=O'),
    E0 = (-165.214,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,375,552.5,462.5,1710,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,368.992,372.476,372.574],'cm^-1')),
        HinderedRotor(inertia=(0.0221847,'amu*angstrom^2'), symmetry=1, barrier=(2.09348,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.682691,'amu*angstrom^2'), symmetry=1, barrier=(15.6964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16382,'amu*angstrom^2'), symmetry=1, barrier=(15.6983,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161514,'amu*angstrom^2'), symmetry=1, barrier=(15.6953,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4218.9,'J/mol'), sigma=(6.48785,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=658.98 K, Pc=35.05 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39102,0.06549,-3.82654e-05,8.80535e-09,-7.19387e-13,-19784.4,27.477], Tmin=(100,'K'), Tmax=(2264.17,'K')), NASAPolynomial(coeffs=[35.3576,0.0135168,-9.15587e-06,1.8014e-09,-1.19076e-13,-37224.9,-168.697], Tmin=(2264.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-165.214,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'C3H5(89)(88)',
    structure = SMILES('[CH2]C=C'),
    E0 = (156.927,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.332071,'amu*angstrom^2'), symmetry=1, barrier=(25.4373,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2161.77,'J/mol'), sigma=(4.85,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.29606,0.00579326,4.33883e-05,-5.99838e-08,2.33791e-11,18908.2,9.02024], Tmin=(100,'K'), Tmax=(942.201,'K')), NASAPolynomial(coeffs=[8.0689,0.0101832,-2.84768e-06,5.00815e-10,-3.79574e-14,16914.6,-19.5287], Tmin=(942.201,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(156.927,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""allyl""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'S(780)(779)',
    structure = SMILES('O=[C]C=CC=O'),
    E0 = (-44.2053,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.882295,'amu*angstrom^2'), symmetry=1, barrier=(20.2857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.879664,'amu*angstrom^2'), symmetry=1, barrier=(20.2252,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3741.66,'J/mol'), sigma=(5.7541,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=584.44 K, Pc=44.56 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70479,0.0574491,-9.74937e-05,9.58172e-08,-3.8062e-11,-5240.81,19.2682], Tmin=(100,'K'), Tmax=(736.114,'K')), NASAPolynomial(coeffs=[5.25854,0.0278824,-1.6346e-05,3.39827e-09,-2.46552e-13,-5486.14,5.0995], Tmin=(736.114,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-44.2053,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJ=O)"""),
)

species(
    label = 'C2H3(28)(29)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36378,0.000265766,2.79621e-05,-3.72987e-08,1.5159e-11,34475,7.9151], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.15027,0.00754021,-2.62998e-06,4.15974e-10,-2.45408e-14,33856.6,1.72812], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(286.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""C2H3""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[CH2]C(=O)C=CC=O(2803)',
    structure = SMILES('C=C([O])C=CC=O'),
    E0 = (-105.655,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2782.5,750,1395,475,1775,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.804583,'amu*angstrom^2'), symmetry=1, barrier=(18.4989,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.805883,'amu*angstrom^2'), symmetry=1, barrier=(18.5288,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16866,0.0608176,-5.94524e-05,2.9485e-08,-5.85639e-12,-12604.2,21.9805], Tmin=(100,'K'), Tmax=(1208.63,'K')), NASAPolynomial(coeffs=[13.3958,0.0203512,-9.23039e-06,1.78305e-09,-1.26336e-13,-15559.8,-39.3255], Tmin=(1208.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-105.655,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)OJ)"""),
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
    label = 'C=C[CH]C(=O)C=CC=O(2764)',
    structure = SMILES('C=CC=C([O])C=CC=O'),
    E0 = (-53.883,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2782.5,750,1395,475,1775,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.885181,'amu*angstrom^2'), symmetry=1, barrier=(20.3521,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.886434,'amu*angstrom^2'), symmetry=1, barrier=(20.3809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.884396,'amu*angstrom^2'), symmetry=1, barrier=(20.334,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.129,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.374117,0.0881774,-8.71591e-05,4.28208e-08,-8.24952e-12,-6316.41,30.155], Tmin=(100,'K'), Tmax=(1264.29,'K')), NASAPolynomial(coeffs=[20.4015,0.0224468,-9.17397e-06,1.69887e-09,-1.18107e-13,-11569.7,-74.9477], Tmin=(1264.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-53.883,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C4H5O(93)(92)',
    structure = SMILES('C=CC[C]=O'),
    E0 = (66.8219,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,458.926],'cm^-1')),
        HinderedRotor(inertia=(0.0997865,'amu*angstrom^2'), symmetry=1, barrier=(14.9157,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.099798,'amu*angstrom^2'), symmetry=1, barrier=(14.9167,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3285.42,'J/mol'), sigma=(5.46087,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=513.18 K, Pc=45.78 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.51804,0.0238835,1.19491e-05,-2.85418e-08,1.09388e-11,8097.53,17.8098], Tmin=(100,'K'), Tmax=(1083.61,'K')), NASAPolynomial(coeffs=[9.78041,0.0178579,-8.47799e-06,1.72441e-09,-1.27255e-13,5303.46,-23.4402], Tmin=(1083.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(66.8219,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O)"""),
)

species(
    label = 'CHCHCHO(2768)',
    structure = SMILES('[CH]=CC=O'),
    E0 = (171.79,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,688.575,691.322,692.403,693.509,695.331],'cm^-1')),
        HinderedRotor(inertia=(0.0100227,'amu*angstrom^2'), symmetry=1, barrier=(3.38495,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (55.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.95468,0.0145151,1.89382e-05,-3.59857e-08,1.47962e-11,20706.8,12.171], Tmin=(100,'K'), Tmax=(981.049,'K')), NASAPolynomial(coeffs=[9.63601,0.00810119,-3.09994e-06,6.30252e-10,-4.9152e-14,18393.6,-25.043], Tmin=(981.049,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(171.79,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CHCHCHO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=[C]CC(=O)C=CC=O(2804)',
    structure = SMILES('C=[C]CC(=O)C=CC=O'),
    E0 = (72.6276,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,375,552.5,462.5,1710,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,356.019,356.019,356.019],'cm^-1')),
        HinderedRotor(inertia=(0.00133001,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165884,'amu*angstrom^2'), symmetry=1, barrier=(14.9203,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165884,'amu*angstrom^2'), symmetry=1, barrier=(14.9203,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165884,'amu*angstrom^2'), symmetry=1, barrier=(14.9203,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.129,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.25843,0.0608742,-3.59562e-05,8.30304e-09,-7.03262e-13,8774.3,24.7235], Tmin=(100,'K'), Tmax=(2845.92,'K')), NASAPolynomial(coeffs=[43.0428,0.00355249,-5.74445e-06,1.22603e-09,-8.15987e-14,-14440.1,-214.694], Tmin=(2845.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(72.6276,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_S)"""),
)

species(
    label = 'C=CCC(=O)[C]=CC=O(2805)',
    structure = SMILES('C=CCC([O])=C=CC=O'),
    E0 = (29.0867,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.614687,'amu*angstrom^2'), symmetry=1, barrier=(14.1329,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.614772,'amu*angstrom^2'), symmetry=1, barrier=(14.1348,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.614994,'amu*angstrom^2'), symmetry=1, barrier=(14.1399,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.129,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.46311,0.0795325,-7.4101e-05,3.59708e-08,-7.1766e-12,3624.15,31.0178], Tmin=(100,'K'), Tmax=(1181.38,'K')), NASAPolynomial(coeffs=[14.0124,0.0336565,-1.58522e-05,3.10024e-09,-2.20643e-13,422.785,-36.6079], Tmin=(1181.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(29.0867,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds) + radical(C=C(C)OJ)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.23755,-0.00332075,1.4003e-05,-1.3424e-08,4.37416e-12,3872.41,3.30835], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.92002,0.00252279,-6.71004e-07,1.05616e-10,-7.43798e-15,3653.43,3.58077], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(32.4782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HCO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[CH]=CC(=O)CC=C(2806)',
    structure = SMILES('[CH]=CC(=O)CC=C'),
    E0 = (205.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,375,552.5,462.5,1710,3120,650,792.5,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,625.128],'cm^-1')),
        HinderedRotor(inertia=(0.0115181,'amu*angstrom^2'), symmetry=1, barrier=(3.3326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.799554,'amu*angstrom^2'), symmetry=1, barrier=(18.3833,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.800305,'amu*angstrom^2'), symmetry=1, barrier=(18.4006,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51193,0.0511228,-2.78333e-05,5.11544e-09,-6.08057e-14,24792.8,24.3932], Tmin=(100,'K'), Tmax=(1699.06,'K')), NASAPolynomial(coeffs=[17.9658,0.0225202,-1.15282e-05,2.22816e-09,-1.52492e-13,17738.9,-68.0134], Tmin=(1699.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C=CCC(=O)C=[C]C=O(2807)',
    structure = SMILES('C=CCC(=O)C=C=C[O]'),
    E0 = (31.5225,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,375,552.5,462.5,1710,540,610,2055,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,240.585,240.775,240.78,241.18],'cm^-1')),
        HinderedRotor(inertia=(0.553511,'amu*angstrom^2'), symmetry=1, barrier=(22.5684,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.551293,'amu*angstrom^2'), symmetry=1, barrier=(22.5756,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.551594,'amu*angstrom^2'), symmetry=1, barrier=(22.569,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.129,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.032432,0.0761582,-5.77989e-05,2.02814e-08,-2.7778e-12,3946.96,31.6609], Tmin=(100,'K'), Tmax=(1734.13,'K')), NASAPolynomial(coeffs=[23.8698,0.0210248,-1.01094e-05,1.94781e-09,-1.34759e-13,-4342.99,-96.812], Tmin=(1734.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.5225,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ)"""),
)

species(
    label = '[CH]=CCC(=O)C=CC=O(2808)',
    structure = SMILES('[CH]=CCC(=O)C=CC=O'),
    E0 = (81.882,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,375,552.5,462.5,1710,3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,797.808],'cm^-1')),
        HinderedRotor(inertia=(0.0328771,'amu*angstrom^2'), symmetry=1, barrier=(14.8352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.645393,'amu*angstrom^2'), symmetry=1, barrier=(14.8389,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.645133,'amu*angstrom^2'), symmetry=1, barrier=(14.8329,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.417323,'amu*angstrom^2'), symmetry=1, barrier=(14.8396,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.129,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90134,0.063268,-3.89451e-05,9.67149e-09,-9.01491e-13,9905.26,25.9542], Tmin=(100,'K'), Tmax=(2449.31,'K')), NASAPolynomial(coeffs=[28.2305,0.0202695,-1.26122e-05,2.50406e-09,-1.69916e-13,-2992.44,-124.655], Tmin=(2449.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(81.882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_P)"""),
)

species(
    label = 'C=CCC(=O)C=C[C]=O(2809)',
    structure = SMILES('C=CCC(=O)[CH]C=C=O'),
    E0 = (-5.49354,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,936.904,937.649,951.75],'cm^-1')),
        HinderedRotor(inertia=(0.240612,'amu*angstrom^2'), symmetry=1, barrier=(5.53215,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.867992,'amu*angstrom^2'), symmetry=1, barrier=(19.9569,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0598732,'amu*angstrom^2'), symmetry=1, barrier=(19.8942,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.866778,'amu*angstrom^2'), symmetry=1, barrier=(19.9289,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.129,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09899,0.0678411,-4.71278e-05,1.52695e-08,-2.00849e-12,-561.097,30.0519], Tmin=(100,'K'), Tmax=(1681.6,'K')), NASAPolynomial(coeffs=[15.2988,0.0340645,-1.6999e-05,3.32503e-09,-2.32744e-13,-5336.8,-45.8343], Tmin=(1681.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-5.49354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + radical(CCJC(C)=C=O)"""),
)

species(
    label = 'S(458)(457)',
    structure = SMILES('C=C[CH]C([O])C=CC=O'),
    E0 = (67.9095,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3025,407.5,1350,352.5,466.122,466.138,466.658,466.797],'cm^-1')),
        HinderedRotor(inertia=(0.150317,'amu*angstrom^2'), symmetry=1, barrier=(23.1705,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150455,'amu*angstrom^2'), symmetry=1, barrier=(23.1786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150011,'amu*angstrom^2'), symmetry=1, barrier=(23.1625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150418,'amu*angstrom^2'), symmetry=1, barrier=(23.174,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4360.96,'J/mol'), sigma=(6.98043,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=681.17 K, Pc=29.09 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.165628,0.0680478,-2.07106e-05,-2.23407e-08,1.26385e-11,8319.62,35.7619], Tmin=(100,'K'), Tmax=(1066.29,'K')), NASAPolynomial(coeffs=[19.6886,0.0286735,-1.29576e-05,2.59562e-09,-1.91068e-13,2231.15,-68.705], Tmin=(1066.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(67.9095,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJC(O)C=C) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C[CH]C(=O)C=CC=O(2810)',
    structure = SMILES('[CH2]CC=C([O])C=CC=O'),
    E0 = (40.9204,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,350,440,435,1725,299.002,299.002,299.003],'cm^-1')),
        HinderedRotor(inertia=(0.193211,'amu*angstrom^2'), symmetry=1, barrier=(12.2578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193212,'amu*angstrom^2'), symmetry=1, barrier=(12.2578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193213,'amu*angstrom^2'), symmetry=1, barrier=(12.2578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193215,'amu*angstrom^2'), symmetry=1, barrier=(12.2578,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.117418,0.0865188,-8.40596e-05,4.29285e-08,-8.95424e-12,5060.45,32.2622], Tmin=(100,'K'), Tmax=(1140.68,'K')), NASAPolynomial(coeffs=[14.9103,0.0346451,-1.58457e-05,3.06122e-09,-2.16646e-13,1685.66,-41.052], Tmin=(1140.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(40.9204,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)OJ) + radical(RCCJ)"""),
)

species(
    label = 'C=CCC([O])[C]=CC=O(2772)',
    structure = SMILES('C=CCC([O])[C]=CC=O'),
    E0 = (235.608,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,298.163,298.173,298.175,298.176],'cm^-1')),
        HinderedRotor(inertia=(0.00189609,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223185,'amu*angstrom^2'), symmetry=1, barrier=(14.0809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223194,'amu*angstrom^2'), symmetry=1, barrier=(14.0809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223197,'amu*angstrom^2'), symmetry=1, barrier=(14.0809,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.654947,0.077786,-6.43548e-05,2.73372e-08,-4.86994e-12,28453.6,33.6719], Tmin=(100,'K'), Tmax=(1280.46,'K')), NASAPolynomial(coeffs=[13.0233,0.0391487,-1.90929e-05,3.77178e-09,-2.68952e-13,25286.2,-29.0558], Tmin=(1280.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(235.608,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'C=CCC(=O)C=[C]C[O](2811)',
    structure = SMILES('C=CCC(=O)C=[C]C[O]'),
    E0 = (232.659,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1685,370,2950,3100,1380,975,1025,1650,375,552.5,462.5,1710,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.80417,0.0570863,-2.87131e-05,5.08418e-09,-2.52727e-13,27994.9,26.083], Tmin=(100,'K'), Tmax=(2533.12,'K')), NASAPolynomial(coeffs=[45.2833,0.00371789,-5.22889e-06,1.04013e-09,-6.44717e-14,2075.53,-227.018], Tmin=(2533.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(232.659,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + radical(Cds_S) + radical(CCOJ)"""),
)

species(
    label = 'C=CCC(=O)[CH]C[C]=O(2812)',
    structure = SMILES('C=CCC([O])=CC[C]=O'),
    E0 = (36.6685,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1855,455,950,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,385.972,385.982,386.031,386.034],'cm^-1')),
        HinderedRotor(inertia=(0.00113113,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115804,'amu*angstrom^2'), symmetry=1, barrier=(12.2474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115875,'amu*angstrom^2'), symmetry=1, barrier=(12.2472,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115856,'amu*angstrom^2'), symmetry=1, barrier=(12.2474,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0191158,0.0756098,-5.3755e-05,1.49313e-08,-5.1272e-13,4563.48,35.8072], Tmin=(100,'K'), Tmax=(1181.62,'K')), NASAPolynomial(coeffs=[18.8737,0.027157,-1.17626e-05,2.24997e-09,-1.59229e-13,-965.571,-62.8434], Tmin=(1181.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(36.6685,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(CCCJ=O)"""),
)

species(
    label = 'C=[C]CC([O])C=CC=O(2770)',
    structure = SMILES('C=[C]CC([O])C=CC=O'),
    E0 = (235.608,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,298.163,298.173,298.175,298.176],'cm^-1')),
        HinderedRotor(inertia=(0.00189609,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223185,'amu*angstrom^2'), symmetry=1, barrier=(14.0809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223194,'amu*angstrom^2'), symmetry=1, barrier=(14.0809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223197,'amu*angstrom^2'), symmetry=1, barrier=(14.0809,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.654947,0.077786,-6.43548e-05,2.73372e-08,-4.86994e-12,28453.6,33.6719], Tmin=(100,'K'), Tmax=(1280.46,'K')), NASAPolynomial(coeffs=[13.0233,0.0391487,-1.90929e-05,3.77178e-09,-2.68952e-13,25286.2,-29.0558], Tmin=(1280.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(235.608,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'C=C[CH]C(=O)C[CH]C=O(2813)',
    structure = SMILES('C=CC=C([O])CC=C[O]'),
    E0 = (-3.56209,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,350,440,435,1725,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.847619,'amu*angstrom^2'), symmetry=1, barrier=(19.4884,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.846683,'amu*angstrom^2'), symmetry=1, barrier=(19.4669,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.847147,'amu*angstrom^2'), symmetry=1, barrier=(19.4776,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.318776,0.0760891,-2.51903e-05,-3.88268e-08,2.53343e-11,-255.512,33.7881], Tmin=(100,'K'), Tmax=(932.108,'K')), NASAPolynomial(coeffs=[25.3822,0.0139103,-2.55479e-06,3.60804e-10,-2.88539e-14,-7136.82,-99.6089], Tmin=(932.108,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-3.56209,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = 'C=CCC([O])C=[C]C=O(2776)',
    structure = SMILES('C=CCC([O])C=C=C[O]'),
    E0 = (194.503,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,540,610,2055,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,203.415,203.44,203.443,203.473,203.497],'cm^-1')),
        HinderedRotor(inertia=(0.723344,'amu*angstrom^2'), symmetry=1, barrier=(21.2503,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.723702,'amu*angstrom^2'), symmetry=1, barrier=(21.25,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.723536,'amu*angstrom^2'), symmetry=1, barrier=(21.2503,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.280118,0.0787841,-4.26309e-05,-1.01905e-08,1.17895e-11,23561.1,35.6077], Tmin=(100,'K'), Tmax=(984.817,'K')), NASAPolynomial(coeffs=[22.1665,0.0217682,-7.81048e-06,1.45407e-09,-1.06216e-13,17483.6,-80.7498], Tmin=(984.817,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(194.503,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=CCC(=O)[C]=CC[O](2814)',
    structure = SMILES('C=CCC([O])=C=CC[O]'),
    E0 = (189.118,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,350,440,435,1725,238.598,238.598,238.598,238.598,3531.16],'cm^-1')),
        HinderedRotor(inertia=(0.00296118,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.354247,'amu*angstrom^2'), symmetry=1, barrier=(14.311,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.35425,'amu*angstrom^2'), symmetry=1, barrier=(14.311,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.563478,0.0804315,-8.08593e-05,4.80446e-08,-1.22322e-11,22865.3,34.0156], Tmin=(100,'K'), Tmax=(926.459,'K')), NASAPolynomial(coeffs=[9.38506,0.0423438,-1.9192e-05,3.66917e-09,-2.57571e-13,21230.8,-7.86952], Tmin=(926.459,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(189.118,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCOJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=[C]C[C](O)C=CC=O(2815)',
    structure = SMILES('C=[C]CC(O)=CC=C[O]'),
    E0 = (96.4749,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.915269,'amu*angstrom^2'), symmetry=1, barrier=(21.0438,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.913694,'amu*angstrom^2'), symmetry=1, barrier=(21.0076,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.914834,'amu*angstrom^2'), symmetry=1, barrier=(21.0338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.914668,'amu*angstrom^2'), symmetry=1, barrier=(21.03,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.705948,0.0845129,-4.18917e-05,-2.74746e-08,2.27989e-11,11790.3,34.1508], Tmin=(100,'K'), Tmax=(926.87,'K')), NASAPolynomial(coeffs=[27.886,0.0101276,-8.19636e-07,2.82049e-11,-5.58605e-15,4385.01,-112.973], Tmin=(926.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(96.4749,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = 'C=C[CH]C(=O)[CH]CC=O(2816)',
    structure = SMILES('[CH2]C=CC([O])=CCC=O'),
    E0 = (-35.4248,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,350,440,435,1725,351.009,351.133,351.257],'cm^-1')),
        HinderedRotor(inertia=(0.00136737,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.276587,'amu*angstrom^2'), symmetry=1, barrier=(24.1895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.276301,'amu*angstrom^2'), symmetry=1, barrier=(24.1884,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.276243,'amu*angstrom^2'), symmetry=1, barrier=(24.1854,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.120753,0.0723982,-3.93744e-05,-1.02366e-09,5.22804e-12,-4109.96,32.3379], Tmin=(100,'K'), Tmax=(1082.72,'K')), NASAPolynomial(coeffs=[18.0518,0.0295664,-1.2471e-05,2.3829e-09,-1.70053e-13,-9365.11,-61.9312], Tmin=(1082.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-35.4248,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + radical(C=CC=CCJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=CC[C](O)C=[C]C=O(2817)',
    structure = SMILES('C=CCC(O)=C[C]=C[O]'),
    E0 = (57.6286,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.08593,'amu*angstrom^2'), symmetry=1, barrier=(24.9676,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08649,'amu*angstrom^2'), symmetry=1, barrier=(24.9805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08582,'amu*angstrom^2'), symmetry=1, barrier=(24.9652,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08595,'amu*angstrom^2'), symmetry=1, barrier=(24.9681,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.25555,0.103776,-0.000105018,5.03148e-08,-8.88412e-12,7185.47,38.8643], Tmin=(100,'K'), Tmax=(1639.51,'K')), NASAPolynomial(coeffs=[27.8768,0.00934983,5.04373e-07,-3.72783e-10,3.11987e-14,115.418,-112.833], Tmin=(1639.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(57.6286,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C=CCC(=O)[C]=C[CH]O(2818)',
    structure = SMILES('C=CCC([O])=[C]C=CO'),
    E0 = (53.9707,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.973677,'amu*angstrom^2'), symmetry=1, barrier=(22.3867,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.97422,'amu*angstrom^2'), symmetry=1, barrier=(22.3992,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.973241,'amu*angstrom^2'), symmetry=1, barrier=(22.3767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.974311,'amu*angstrom^2'), symmetry=1, barrier=(22.4013,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.10918,0.104768,-0.000110051,5.50396e-08,-1.01512e-11,6736.57,38.4121], Tmin=(100,'K'), Tmax=(1563.61,'K')), NASAPolynomial(coeffs=[27.4056,0.00917618,9.23397e-07,-4.92225e-10,4.12909e-14,-37.714,-109.32], Tmin=(1563.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(53.9707,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C=[C]CC(=O)C[CH]C=O(2819)',
    structure = SMILES('C=[C]CC(=O)CC=C[O]'),
    E0 = (107.27,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1685,370,2950,3100,1380,975,1025,1650,375,552.5,462.5,1710,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0358925,0.0715668,-3.16427e-05,-9.52788e-09,7.37389e-12,13057.3,32.852], Tmin=(100,'K'), Tmax=(1150.75,'K')), NASAPolynomial(coeffs=[20.9044,0.0287361,-1.45378e-05,2.99718e-09,-2.21066e-13,6287.36,-79.3037], Tmin=(1150.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(107.27,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]CCC(=O)[C]=CC=O(2820)',
    structure = SMILES('[CH2]CCC([O])=C=CC=O'),
    E0 = (105.281,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3010,987.5,1337.5,450,1655,350,440,435,1725,3000,3100,440,815,1455,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,278.195,278.205,278.212],'cm^-1')),
        HinderedRotor(inertia=(0.00217812,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200202,'amu*angstrom^2'), symmetry=1, barrier=(10.9965,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200213,'amu*angstrom^2'), symmetry=1, barrier=(10.9966,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200204,'amu*angstrom^2'), symmetry=1, barrier=(10.9966,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.315543,0.0865361,-8.82386e-05,4.94719e-08,-1.16267e-11,12790.5,32.8465], Tmin=(100,'K'), Tmax=(1005.26,'K')), NASAPolynomial(coeffs=[11.9099,0.0404007,-1.93969e-05,3.81698e-09,-2.72554e-13,10459.5,-23.1504], Tmin=(1005.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(105.281,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds) + radical(C=C(C)OJ) + radical(RCCJ)"""),
)

species(
    label = '[CH]=CCC([O])C=CC=O(2773)',
    structure = SMILES('[CH]=CCC([O])C=CC=O'),
    E0 = (244.862,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,352.358,352.358,352.358],'cm^-1')),
        HinderedRotor(inertia=(0.155497,'amu*angstrom^2'), symmetry=1, barrier=(13.6999,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155497,'amu*angstrom^2'), symmetry=1, barrier=(13.6999,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155497,'amu*angstrom^2'), symmetry=1, barrier=(13.6999,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155497,'amu*angstrom^2'), symmetry=1, barrier=(13.6999,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.397267,0.0792587,-6.50371e-05,2.66359e-08,-4.46215e-12,29579.3,34.5254], Tmin=(100,'K'), Tmax=(1382.2,'K')), NASAPolynomial(coeffs=[15.8483,0.0345444,-1.65118e-05,3.23098e-09,-2.28873e-13,25308,-45.0179], Tmin=(1382.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(244.862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(CC(C)OJ) + radical(Cds_P)"""),
)

species(
    label = 'C=CCC([O])C=C[C]=O(2779)',
    structure = SMILES('C=CCC([O])[CH]C=C=O'),
    E0 = (129.565,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0599547,0.0863159,-7.78899e-05,3.63229e-08,-6.8642e-12,15731.6,31.9605], Tmin=(100,'K'), Tmax=(1258.12,'K')), NASAPolynomial(coeffs=[16.6583,0.0331629,-1.4518e-05,2.74277e-09,-1.91525e-13,11524.9,-52.5342], Tmin=(1258.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(129.565,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + radical(C=CCJCO) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=[C]CC(=O)[CH]CC=O(2821)',
    structure = SMILES('C=[C]CC([O])=CCC=O'),
    E0 = (114.55,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,422.097,422.102,422.105,422.106],'cm^-1')),
        HinderedRotor(inertia=(0.00094612,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0985207,'amu*angstrom^2'), symmetry=1, barrier=(12.4564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0985167,'amu*angstrom^2'), symmetry=1, barrier=(12.4563,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0985215,'amu*angstrom^2'), symmetry=1, barrier=(12.4564,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.193356,0.0755718,-5.71635e-05,2.09299e-08,-3.05616e-12,13920.6,34.8419], Tmin=(100,'K'), Tmax=(1605.66,'K')), NASAPolynomial(coeffs=[19.438,0.0276305,-1.23776e-05,2.33518e-09,-1.61025e-13,7740.47,-67.1159], Tmin=(1605.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(114.55,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'C[CH]CC(=O)[C]=CC=O(2822)',
    structure = SMILES('C[CH]CC([O])=C=CC=O'),
    E0 = (94.4811,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,540,610,2055,263.732,263.734,4000],'cm^-1')),
        HinderedRotor(inertia=(0.356164,'amu*angstrom^2'), symmetry=1, barrier=(17.5793,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.356159,'amu*angstrom^2'), symmetry=1, barrier=(17.5793,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160261,'amu*angstrom^2'), symmetry=1, barrier=(7.91021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160261,'amu*angstrom^2'), symmetry=1, barrier=(7.91021,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.331466,0.0876841,-0.000102974,7.33361e-08,-2.23428e-11,11489.3,33.7412], Tmin=(100,'K'), Tmax=(783.502,'K')), NASAPolynomial(coeffs=[8.63814,0.0452769,-2.17878e-05,4.25764e-09,-3.01635e-13,10187.6,-4.30716], Tmin=(783.502,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(94.4811,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds) + radical(C=C(C)OJ) + radical(RCCJC)"""),
)

species(
    label = '[CH]=CC[C](O)C=CC=O(2823)',
    structure = SMILES('[CH]=CCC(O)=CC=C[O]'),
    E0 = (105.729,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3120,650,792.5,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.957669,'amu*angstrom^2'), symmetry=1, barrier=(22.0187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.956663,'amu*angstrom^2'), symmetry=1, barrier=(21.9956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.956609,'amu*angstrom^2'), symmetry=1, barrier=(21.9943,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.955314,'amu*angstrom^2'), symmetry=1, barrier=(21.9646,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.749247,0.083459,-3.38008e-05,-3.93914e-08,2.79072e-11,12906.8,34.236], Tmin=(100,'K'), Tmax=(923.917,'K')), NASAPolynomial(coeffs=[29.3428,0.00774963,5.17611e-07,-2.25536e-10,1.10997e-14,5017.07,-121.164], Tmin=(923.917,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(105.729,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Cds_P)"""),
)

species(
    label = 'C=CC[C](O)C=C[C]=O(2824)',
    structure = SMILES('C=CCC(O)=C[CH][C]=O'),
    E0 = (60.9925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1855,455,950,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.136119,0.0712849,-1.8086e-05,-3.70725e-08,2.10986e-11,7502.12,35.2372], Tmin=(100,'K'), Tmax=(993.259,'K')), NASAPolynomial(coeffs=[23.9975,0.0191824,-7.49175e-06,1.51803e-09,-1.17313e-13,483.858,-92.2262], Tmin=(993.259,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(60.9925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O) + radical(CCJC=O)"""),
)

species(
    label = '[CH]=CCC(=O)C[CH]C=O(2825)',
    structure = SMILES('[CH]=CCC(=O)CC=C[O]'),
    E0 = (116.524,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,408.816,410.261,410.633,410.728],'cm^-1')),
        HinderedRotor(inertia=(0.169283,'amu*angstrom^2'), symmetry=1, barrier=(20.2398,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169662,'amu*angstrom^2'), symmetry=1, barrier=(20.1951,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172409,'amu*angstrom^2'), symmetry=1, barrier=(20.232,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172391,'amu*angstrom^2'), symmetry=1, barrier=(20.2269,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00122202,0.0704993,-2.39432e-05,-2.02593e-08,1.16644e-11,14173.3,32.9006], Tmin=(100,'K'), Tmax=(1103.1,'K')), NASAPolynomial(coeffs=[21.7226,0.0273707,-1.37551e-05,2.86951e-09,-2.14534e-13,7212.96,-83.8518], Tmin=(1103.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(116.524,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]CCC(=O)C=[C]C=O(2826)',
    structure = SMILES('[CH2]CCC(=O)C=C=C[O]'),
    E0 = (102.804,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,277.736,277.737,277.739,4000],'cm^-1')),
        HinderedRotor(inertia=(0.417686,'amu*angstrom^2'), symmetry=1, barrier=(22.864,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18866,'amu*angstrom^2'), symmetry=1, barrier=(10.3271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.417692,'amu*angstrom^2'), symmetry=1, barrier=(22.864,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.417694,'amu*angstrom^2'), symmetry=1, barrier=(22.864,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.316892,0.0862895,-8.66784e-05,4.82361e-08,-1.12922e-11,12492.7,32.1371], Tmin=(100,'K'), Tmax=(1007.45,'K')), NASAPolynomial(coeffs=[11.6672,0.0412226,-1.9576e-05,3.83063e-09,-2.72648e-13,10205.8,-22.7055], Tmin=(1007.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(102.804,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(RCCJ) + radical(C=COJ)"""),
)

species(
    label = 'C=C[CH]C(=O)C=CC[O](2827)',
    structure = SMILES('C=CC=C([O])C=CC[O]'),
    E0 = (106.149,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,350,440,435,1725,265.943,266.083,266.125,266.431,266.652],'cm^-1')),
        HinderedRotor(inertia=(0.241144,'amu*angstrom^2'), symmetry=1, barrier=(12.0906,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.239113,'amu*angstrom^2'), symmetry=1, barrier=(12.0875,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.240596,'amu*angstrom^2'), symmetry=1, barrier=(12.0923,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0363341,0.086497,-8.58514e-05,4.56095e-08,-9.7745e-12,12914.1,32.2856], Tmin=(100,'K'), Tmax=(1125.54,'K')), NASAPolynomial(coeffs=[15.3964,0.0316516,-1.27594e-05,2.31657e-09,-1.58496e-13,9440.02,-43.9934], Tmin=(1125.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(106.149,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(CCOJ)"""),
)

species(
    label = '[CH]=CCC(=O)[CH]CC=O(2828)',
    structure = SMILES('[CH]=CCC([O])=CCC=O'),
    E0 = (123.804,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3120,650,792.5,1650,350,440,435,1725,251.431,251.437,251.438],'cm^-1')),
        HinderedRotor(inertia=(0.00266632,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00266623,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.373607,'amu*angstrom^2'), symmetry=1, barrier=(16.7607,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.373615,'amu*angstrom^2'), symmetry=1, barrier=(16.7606,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.162903,0.0779646,-6.01622e-05,2.23075e-08,-3.25482e-12,15051.5,36.0685], Tmin=(100,'K'), Tmax=(1634.37,'K')), NASAPolynomial(coeffs=[22.1571,0.0233379,-1.00266e-05,1.85695e-09,-1.26625e-13,7755.64,-82.5777], Tmin=(1634.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(123.804,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C[CH]CC(=O)C=[C]C=O(2829)',
    structure = SMILES('C[CH]CC(=O)C=C=C[O]'),
    E0 = (97.4599,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,3025,407.5,1350,352.5,411.592,411.71,411.747,411.807],'cm^-1')),
        HinderedRotor(inertia=(0.116393,'amu*angstrom^2'), symmetry=1, barrier=(14.0018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0208263,'amu*angstrom^2'), symmetry=1, barrier=(2.50607,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.608992,'amu*angstrom^2'), symmetry=1, barrier=(14.0019,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116358,'amu*angstrom^2'), symmetry=1, barrier=(14.0017,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.359742,0.0818709,-7.36661e-05,3.51517e-08,-6.96106e-12,11851.2,32.6797], Tmin=(100,'K'), Tmax=(1183.14,'K')), NASAPolynomial(coeffs=[13.5408,0.037308,-1.7169e-05,3.31728e-09,-2.34413e-13,8732.17,-33.1281], Tmin=(1183.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(97.4599,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(CCJCC=O) + radical(C=COJ)"""),
)

species(
    label = 'C=C[CH]C(=O)C=C[CH]O(2830)',
    structure = SMILES('[CH2]C=CC([O])=CC=CO'),
    E0 = (-57.1574,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,3000,3100,440,815,1455,1000,350,440,435,1725,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.15639,'amu*angstrom^2'), symmetry=1, barrier=(26.5877,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14897,'amu*angstrom^2'), symmetry=1, barrier=(26.417,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15684,'amu*angstrom^2'), symmetry=1, barrier=(26.598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15,'amu*angstrom^2'), symmetry=1, barrier=(26.4407,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.727654,0.0829658,-2.89863e-05,-4.84005e-08,3.29517e-11,-6684.47,31.5694], Tmin=(100,'K'), Tmax=(901.662,'K')), NASAPolynomial(coeffs=[29.5731,0.00626649,2.583e-06,-7.41832e-10,5.1628e-14,-14495.1,-124.489], Tmin=(901.662,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-57.1574,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=[C]CC(=O)C=CC[O](2831)',
    structure = SMILES('C=[C]CC(=O)C=CC[O]'),
    E0 = (232.659,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1685,370,2950,3100,1380,975,1025,1650,375,552.5,462.5,1710,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.80417,0.0570863,-2.87131e-05,5.08418e-09,-2.52727e-13,27994.9,26.083], Tmin=(100,'K'), Tmax=(2533.12,'K')), NASAPolynomial(coeffs=[45.2833,0.00371789,-5.22889e-06,1.04013e-09,-6.44717e-14,2075.53,-227.018], Tmin=(2533.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(232.659,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + radical(Cds_S) + radical(CCOJ)"""),
)

species(
    label = '[CH2]CCC(=O)C=C[C]=O(2832)',
    structure = SMILES('[CH2]CCC(=O)[CH]C=C=O'),
    E0 = (65.7882,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13311,0.0765428,-3.72576e-05,-8.7626e-08,1.09809e-10,8001.7,31.9463], Tmin=(100,'K'), Tmax=(464.229,'K')), NASAPolynomial(coeffs=[7.26298,0.0479467,-2.31228e-05,4.46782e-09,-3.12246e-13,7171.57,4.26608], Tmin=(464.229,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(65.7882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-(Cdd-O2d)CsH) + radical(CCJC(C)=C=O) + radical(RCCJ)"""),
)

species(
    label = 'C=[C]CC(=O)C=C[CH]O(2833)',
    structure = SMILES('C=[C]CC([O])=CC=CO'),
    E0 = (92.8171,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.835352,'amu*angstrom^2'), symmetry=1, barrier=(19.2064,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.835032,'amu*angstrom^2'), symmetry=1, barrier=(19.199,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.835386,'amu*angstrom^2'), symmetry=1, barrier=(19.2072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.835631,'amu*angstrom^2'), symmetry=1, barrier=(19.2128,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.11863,0.103339,-0.000106984,5.23013e-08,-9.42349e-12,11410.5,39.332], Tmin=(100,'K'), Tmax=(1599.51,'K')), NASAPolynomial(coeffs=[27.9497,0.00823509,8.7432e-07,-4.35635e-10,3.55135e-14,4338.52,-111.891], Tmin=(1599.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(92.8171,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'C[CH]CC(=O)C=C[C]=O(2834)',
    structure = SMILES('C[CH]CC(=O)[CH]C=C=O'),
    E0 = (60.4439,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.447834,0.0848502,-9.86054e-05,7.22439e-08,-2.2975e-11,7391.48,34.8912], Tmin=(100,'K'), Tmax=(748.827,'K')), NASAPolynomial(coeffs=[7.57616,0.0467738,-2.23353e-05,4.34385e-09,-3.06743e-13,6323.87,2.56282], Tmin=(748.827,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(60.4439,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-(Cdd-O2d)CsH) + radical(CCJCC=O) + radical(CCJC(C)=C=O)"""),
)

species(
    label = '[CH]=CCC(=O)C=CC[O](2835)',
    structure = SMILES('[CH]=CCC(=O)C=CC[O]'),
    E0 = (241.914,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180,890.762],'cm^-1')),
        HinderedRotor(inertia=(0.00752853,'amu*angstrom^2'), symmetry=1, barrier=(4.23899,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00752878,'amu*angstrom^2'), symmetry=1, barrier=(4.23901,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184314,'amu*angstrom^2'), symmetry=1, barrier=(4.23773,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.611237,'amu*angstrom^2'), symmetry=1, barrier=(14.0536,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.44295,0.0595183,-3.17902e-05,6.51874e-09,-4.65847e-13,29126.1,27.3292], Tmin=(100,'K'), Tmax=(2634.43,'K')), NASAPolynomial(coeffs=[46.0046,0.00241192,-4.4196e-06,8.94287e-10,-5.56536e-14,3038.57,-230.978], Tmin=(2634.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(241.914,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_P)"""),
)

species(
    label = 'C[C]CC(=O)C=CC=O(2840)',
    structure = SMILES('C[C]CC(=O)C=CC=O'),
    E0 = (131.275,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,375,552.5,462.5,1710,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.203458,0.0941451,-0.000135463,1.29452e-07,-5.15175e-11,15915.2,41.8574], Tmin=(100,'K'), Tmax=(753.9,'K')), NASAPolynomial(coeffs=[3.53582,0.0589974,-3.07774e-05,6.14788e-09,-4.37708e-13,15909.1,30.0142], Tmin=(753.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(131.275,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'C=CCC(=O)[C]CC=O(2841)',
    structure = SMILES('C=CCC(=O)[C]CC=O'),
    E0 = (152.239,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,375,552.5,462.5,1710,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23544,0.0694094,-4.3886e-05,1.25525e-08,-1.4442e-12,18400,39.8013], Tmin=(100,'K'), Tmax=(1865.87,'K')), NASAPolynomial(coeffs=[15.6868,0.0384288,-1.89802e-05,3.65374e-09,-2.51876e-13,13007.2,-38.9316], Tmin=(1865.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(152.239,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'C=CCC(=O)C[C]C=O(2842)',
    structure = SMILES('C=CCC(=O)C[C]C=O'),
    E0 = (152.239,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,375,552.5,462.5,1710,3010,987.5,1337.5,450,1655,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23544,0.0694094,-4.3886e-05,1.25525e-08,-1.4442e-12,18400,39.8013], Tmin=(100,'K'), Tmax=(1865.87,'K')), NASAPolynomial(coeffs=[15.6868,0.0384288,-1.89802e-05,3.65375e-09,-2.51876e-13,13007.2,-38.9315], Tmin=(1865.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(152.239,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-CsCsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-OdCsH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = '[CH]CCC(=O)C=CC=O(2843)',
    structure = SMILES('[CH]CCC(=O)C=CC=O'),
    E0 = (151.921,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,375,552.5,462.5,1710,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0332997,0.101017,-0.000164241,1.63084e-07,-6.42386e-11,18401.3,32.4304], Tmin=(100,'K'), Tmax=(788.725,'K')), NASAPolynomial(coeffs=[2.92192,0.0587278,-3.12508e-05,6.24321e-09,-4.42324e-13,18805.3,24.63], Tmin=(788.725,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(151.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-CsCsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'O=CC=C[C]1C[CH]CO1(2845)',
    structure = SMILES('[O]C=CC=C1C[CH]CO1'),
    E0 = (34.5907,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.531048,0.0391628,0.000101223,-1.83046e-07,7.97267e-11,4319.63,30.6623], Tmin=(100,'K'), Tmax=(923.831,'K')), NASAPolynomial(coeffs=[31.0635,0.00283317,4.55003e-06,-9.53486e-10,5.23523e-14,-5412.79,-136.362], Tmin=(923.831,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(34.5907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(Cyclopentane) + radical(C=COJ) + radical(CCJCO)"""),
)

species(
    label = 'O=CC1[CH]C(=O)C[CH]C1(2846)',
    structure = SMILES('[O]C1=CC(C=O)C[CH]C1'),
    E0 = (-17.0393,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2782.5,750,1395,475,1775,1000,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.890382,0.0586368,-2.12412e-05,-6.91648e-09,4.8767e-12,-1929.31,29.2788], Tmin=(100,'K'), Tmax=(1136.7,'K')), NASAPolynomial(coeffs=[12.0726,0.0365662,-1.49183e-05,2.74811e-09,-1.90045e-13,-5587.78,-31.012], Tmin=(1136.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-17.0393,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclohexene) + radical(RCCJCC) + radical(C=C(C)OJ)"""),
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
    label = 'C=CCC(=O)C=C(2836)',
    structure = SMILES('C=CCC(=O)C=C'),
    E0 = (-41.7261,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,375,552.5,462.5,1710,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,663.989,664.51],'cm^-1')),
        HinderedRotor(inertia=(0.164597,'amu*angstrom^2'), symmetry=1, barrier=(3.78441,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0572404,'amu*angstrom^2'), symmetry=1, barrier=(17.9824,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.782279,'amu*angstrom^2'), symmetry=1, barrier=(17.9861,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46718,0.0490731,-1.65731e-05,-5.0927e-09,2.79109e-12,-4922.08,24.1439], Tmin=(100,'K'), Tmax=(1399.71,'K')), NASAPolynomial(coeffs=[13.8719,0.0307146,-1.52144e-05,2.98344e-09,-2.09419e-13,-10068.9,-45.8533], Tmin=(1399.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-41.7261,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=COC(=C)C=CC=O(2837)',
    structure = SMILES('C=COC(=C)C=CC=O'),
    E0 = (-107.003,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2782.5,750,1395,475,1775,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.840964,'amu*angstrom^2'), symmetry=1, barrier=(19.3354,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.84164,'amu*angstrom^2'), symmetry=1, barrier=(19.351,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.841223,'amu*angstrom^2'), symmetry=1, barrier=(19.3414,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.841032,'amu*angstrom^2'), symmetry=1, barrier=(19.337,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0889732,0.0857918,-7.71986e-05,3.4954e-08,-6.35214e-12,-12719,30.8148], Tmin=(100,'K'), Tmax=(1310.93,'K')), NASAPolynomial(coeffs=[18.1281,0.0302061,-1.35954e-05,2.60855e-09,-1.83667e-13,-17495.2,-62.0038], Tmin=(1310.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-107.003,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

species(
    label = 'S(778)(777)',
    structure = SMILES('C=CC=C(O)C=CC=O'),
    E0 = (-191.688,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.95864,'amu*angstrom^2'), symmetry=1, barrier=(22.041,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.961854,'amu*angstrom^2'), symmetry=1, barrier=(22.1149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.958529,'amu*angstrom^2'), symmetry=1, barrier=(22.0385,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.958395,'amu*angstrom^2'), symmetry=1, barrier=(22.0354,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4437.48,'J/mol'), sigma=(6.82012,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=693.12 K, Pc=31.74 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.785439,0.0927288,-8.19227e-05,2.91418e-08,-1.90073e-12,-22871.5,30.1277], Tmin=(100,'K'), Tmax=(1031.1,'K')), NASAPolynomial(coeffs=[23.4878,0.0201354,-7.69741e-06,1.43995e-09,-1.03445e-13,-29023.8,-93.2806], Tmin=(1031.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-191.688,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CCC(O)=C=CC=O(2838)',
    structure = SMILES('C=CCC(O)=C=CC=O'),
    E0 = (-108.718,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.788706,'amu*angstrom^2'), symmetry=1, barrier=(18.1339,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.788832,'amu*angstrom^2'), symmetry=1, barrier=(18.1368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.789678,'amu*angstrom^2'), symmetry=1, barrier=(18.1563,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.788374,'amu*angstrom^2'), symmetry=1, barrier=(18.1263,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.261835,0.0877121,-8.11087e-05,3.73437e-08,-6.83097e-12,-12917.4,32.1187], Tmin=(100,'K'), Tmax=(1313.29,'K')), NASAPolynomial(coeffs=[19.5727,0.0273001,-1.21077e-05,2.31655e-09,-1.63128e-13,-18127.1,-68.977], Tmin=(1313.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-108.718,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=CCC(=O)C=C=CO(2839)',
    structure = SMILES('C=CCC(=O)C=C=CO'),
    E0 = (-109.94,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.927847,'amu*angstrom^2'), symmetry=1, barrier=(21.333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.927603,'amu*angstrom^2'), symmetry=1, barrier=(21.3274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.928123,'amu*angstrom^2'), symmetry=1, barrier=(21.3394,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.926834,'amu*angstrom^2'), symmetry=1, barrier=(21.3097,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.392564,0.0829515,-6.24901e-05,1.82704e-08,-8.33914e-13,-13053.2,31.5134], Tmin=(100,'K'), Tmax=(1178.76,'K')), NASAPolynomial(coeffs=[21.9318,0.0250375,-1.14965e-05,2.27031e-09,-1.63739e-13,-19555.7,-85.1179], Tmin=(1178.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-109.94,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=CCC(C=O)C=C=O(2844)',
    structure = SMILES('C=CCC(C=O)C=C=O'),
    E0 = (-135.665,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,264.65,264.65,264.651],'cm^-1')),
        HinderedRotor(inertia=(0.00204861,'amu*angstrom^2'), symmetry=1, barrier=(23.26,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.467991,'amu*angstrom^2'), symmetry=1, barrier=(23.26,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.467989,'amu*angstrom^2'), symmetry=1, barrier=(23.26,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.209396,'amu*angstrom^2'), symmetry=1, barrier=(10.4074,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.512859,0.0739656,-5.75625e-05,2.3238e-08,-3.8681e-12,-16189.2,33.0393], Tmin=(100,'K'), Tmax=(1391.68,'K')), NASAPolynomial(coeffs=[14.3058,0.0343223,-1.48341e-05,2.7698e-09,-1.9126e-13,-20028.3,-38.0624], Tmin=(1391.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-135.665,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH) + group(Cds-CdsHH)"""),
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
    E0 = (112.722,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (181.802,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (157.909,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (238.612,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (284.42,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (245.075,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (237.848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (247.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (293.674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (206.299,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (90.7709,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (63.7817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (258.469,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (255.521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (64.5214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (324.576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (89.5319,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (283.471,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (278.087,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (114.257,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (3.80022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (118.297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (93.1957,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (125.052,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (144.506,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (269.835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (154.538,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (153.775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (133.706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (123.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (95.5105,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (134.306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (142.029,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (145.374,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (148.777,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (136.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (-17.9324,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (271.884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (90.7614,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (132.042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (85.4171,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (266.887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (218.715,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (189.277,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (189.277,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (188.959,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (184.903,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (-17.0393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (144.361,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (206.797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (-52.0142,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (34.3955,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (33.9267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (11.3849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction194',
    reactants = ['C3H5(89)(88)', 'S(780)(779)'],
    products = ['S(459)(458)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(6.24e+08,'m^3/(mol*s)'), n=-0.5, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_pri_rad;CO_sec_rad] for rate rule [C_rad/H2/Cd;CO_rad/OneDe]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction195',
    reactants = ['C2H3(28)(29)', '[CH2]C(=O)C=CC=O(2803)'],
    products = ['S(459)(458)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.41031e+09,'m^3/(mol*s)'), n=-0.9855, Ea=(1.09644,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_rad;C_pri_rad] for rate rule [Cd_pri_rad;C_rad/H2/CO]
Euclidian distance = 1.41421356237
family: R_Recombination"""),
)

reaction(
    label = 'reaction196',
    reactants = ['H(3)(3)', 'C=C[CH]C(=O)C=CC=O(2764)'],
    products = ['S(459)(458)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(78817,'m^3/(mol*s)'), n=0.288419, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [C_rad/H/TwoDe;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['C4H5O(93)(92)', 'CHCHCHO(2768)'],
    products = ['S(459)(458)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.6647e+07,'m^3/(mol*s)'), n=-0.0666667, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_rad;Cd_pri_rad] + [CO_rad/NonDe;Y_rad] for rate rule [CO_rad/NonDe;Cd_pri_rad]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction198',
    reactants = ['H(3)(3)', 'C=[C]CC(=O)C=CC=O(2804)'],
    products = ['S(459)(458)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction199',
    reactants = ['H(3)(3)', 'C=CCC(=O)[C]=CC=O(2805)'],
    products = ['S(459)(458)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad/OneDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction200',
    reactants = ['HCO(14)(15)', '[CH]=CC(=O)CC=C(2806)'],
    products = ['S(459)(458)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.81e+13,'cm^3/(mol*s)','*|/',3), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 91 used for Cd_pri_rad;CO_pri_rad
Exact match found for rate rule [Cd_pri_rad;CO_pri_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction201',
    reactants = ['H(3)(3)', 'C=CCC(=O)C=[C]C=O(2807)'],
    products = ['S(459)(458)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad/OneDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction202',
    reactants = ['H(3)(3)', '[CH]=CCC(=O)C=CC=O(2808)'],
    products = ['S(459)(458)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.21e+14,'cm^3/(mol*s)','+|-',4.82e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 60 used for H_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction203',
    reactants = ['H(3)(3)', 'C=CCC(=O)C=C[C]=O(2809)'],
    products = ['S(459)(458)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.53107e+07,'m^3/(mol*s)'), n=-0.133333, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_sec_rad;Y_rad] for rate rule [CO_rad/OneDe;H_rad]
Euclidian distance = 1.41421356237
family: R_Recombination"""),
)

reaction(
    label = 'reaction171',
    reactants = ['S(458)(457)'],
    products = ['S(459)(458)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction204',
    reactants = ['[CH2]C[CH]C(=O)C=CC=O(2810)'],
    products = ['S(459)(458)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction205',
    reactants = ['C=CCC([O])[C]=CC=O(2772)'],
    products = ['S(459)(458)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction206',
    reactants = ['C=CCC(=O)C=[C]C[O](2811)'],
    products = ['S(459)(458)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction207',
    reactants = ['C=CCC(=O)[CH]C[C]=O(2812)'],
    products = ['S(459)(458)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.94659e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction208',
    reactants = ['C=[C]CC([O])C=CC=O(2770)'],
    products = ['S(459)(458)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction209',
    reactants = ['C=C[CH]C(=O)C[CH]C=O(2813)'],
    products = ['S(459)(458)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.08e+10,'s^-1'), n=-0.305, Ea=(93.094,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3;Y_rad_De;XH_Rrad_De] for rate rule [R3radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction210',
    reactants = ['C=CCC([O])C=[C]C=O(2776)'],
    products = ['S(459)(458)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction211',
    reactants = ['C=CCC(=O)[C]=CC[O](2814)'],
    products = ['S(459)(458)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction212',
    reactants = ['C=[C]C[C](O)C=CC=O(2815)'],
    products = ['S(459)(458)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.59e+08,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad_De] for rate rule [R4radEndo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction213',
    reactants = ['C=C[CH]C(=O)[CH]CC=O(2816)'],
    products = ['S(459)(458)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction214',
    reactants = ['C=CC[C](O)C=[C]C=O(2817)'],
    products = ['S(459)(458)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(60.668,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction215',
    reactants = ['C=CCC(=O)[C]=C[CH]O(2818)'],
    products = ['S(459)(458)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.48312e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction216',
    reactants = ['C=[C]CC(=O)C[CH]C=O(2819)'],
    products = ['S(459)(458)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(5.18e+08,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad_De] for rate rule [R4radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction217',
    reactants = ['[CH2]CCC(=O)[C]=CC=O(2820)'],
    products = ['S(459)(458)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction218',
    reactants = ['[CH]=CCC([O])C=CC=O(2773)'],
    products = ['S(459)(458)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction219',
    reactants = ['C=CCC([O])C=C[C]=O(2779)'],
    products = ['S(459)(458)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction220',
    reactants = ['C=[C]CC(=O)[CH]CC=O(2821)'],
    products = ['S(459)(458)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(3.70659e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction221',
    reactants = ['C[CH]CC(=O)[C]=CC=O(2822)'],
    products = ['S(459)(458)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(5.55988e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction222',
    reactants = ['[CH]=CC[C](O)C=CC=O(2823)'],
    products = ['S(459)(458)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.07e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_De] for rate rule [R5radEndo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction223',
    reactants = ['C=CC[C](O)C=C[C]=O(2824)'],
    products = ['S(459)(458)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(3.21e+09,'s^-1'), n=0.137, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_NDe] for rate rule [R5radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction224',
    reactants = ['[CH]=CCC(=O)C[CH]C=O(2825)'],
    products = ['S(459)(458)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.14e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_De] for rate rule [R5radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction225',
    reactants = ['[CH2]CCC(=O)C=[C]C=O(2826)'],
    products = ['S(459)(458)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(3.70659e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction226',
    reactants = ['C=C[CH]C(=O)C=CC[O](2827)'],
    products = ['S(459)(458)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(3.70659e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction227',
    reactants = ['[CH]=CCC(=O)[CH]CC=O(2828)'],
    products = ['S(459)(458)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction228',
    reactants = ['C[CH]CC(=O)C=[C]C=O(2829)'],
    products = ['S(459)(458)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(5.55988e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad_De;XH_Rrad] for rate rule [R6radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction229',
    reactants = ['C=C[CH]C(=O)C=C[CH]O(2830)'],
    products = ['S(459)(458)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad_De;XH_Rrad] for rate rule [R6radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction230',
    reactants = ['C=[C]CC(=O)C=CC[O](2831)'],
    products = ['S(459)(458)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(3.70659e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad_De;XH_Rrad] for rate rule [R6radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction231',
    reactants = ['[CH2]CCC(=O)C=C[C]=O(2832)'],
    products = ['S(459)(458)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction232',
    reactants = ['C=[C]CC(=O)C=C[CH]O(2833)'],
    products = ['S(459)(458)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad_De;XH_Rrad] for rate rule [R7radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction233',
    reactants = ['C[CH]CC(=O)C=C[C]=O(2834)'],
    products = ['S(459)(458)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction234',
    reactants = ['[CH]=CCC(=O)C=CC[O](2835)'],
    products = ['S(459)(458)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction240',
    reactants = ['C[C]CC(=O)C=CC=O(2840)'],
    products = ['S(459)(458)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(4.85495e+16,'s^-1'), n=-0.885455, Ea=(87.4392,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [CsJ2-C;CsJ2(CsC);CH] for rate rule [CsJ2-C;CsJ2(CsC);CH3]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction241',
    reactants = ['C=CCC(=O)[C]CC=O(2841)'],
    products = ['S(459)(458)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(2.90176e+13,'s^-1'), n=-0.332469, Ea=(37.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [singletcarbene_CH;CsJ2C;CH2(C)] + [CsJ2-C;CsJ2C;CH] for rate rule [CsJ2-C;CsJ2C;CH2(C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction242',
    reactants = ['C=CCC(=O)C[C]C=O(2842)'],
    products = ['S(459)(458)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(2.90176e+13,'s^-1'), n=-0.332469, Ea=(37.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [singletcarbene_CH;CsJ2C;CH2(C)] + [CsJ2-C;CsJ2C;CH] for rate rule [CsJ2-C;CsJ2C;CH2(C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction243',
    reactants = ['[CH]CCC(=O)C=CC=O(2843)'],
    products = ['S(459)(458)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(2.90176e+13,'s^-1'), n=-0.332469, Ea=(37.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [singletcarbene_CH;singletcarbene;CH2(C)] + [CsJ2-C;singletcarbene;CH] for rate rule [CsJ2-C;CsJ2H;CH2(C)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction245',
    reactants = ['O=CC=C[C]1C[CH]CO1(2845)'],
    products = ['S(459)(458)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(7.69248e+14,'s^-1'), n=-0.917475, Ea=(150.312,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for R5JJ
Exact match found for rate rule [R5JJ]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction246',
    reactants = ['O=CC1[CH]C(=O)C[CH]C1(2846)'],
    products = ['S(459)(458)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RJJ] for rate rule [R6JJ]
Euclidian distance = 1.0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction235',
    reactants = ['CO(10)(11)', 'C=CCC(=O)C=C(2836)'],
    products = ['S(459)(458)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(0.118397,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO;R_H] for rate rule [CO;Cd_pri]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction236',
    reactants = ['C=COC(=C)C=CC=O(2837)'],
    products = ['S(459)(458)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction237',
    reactants = ['S(778)(777)'],
    products = ['S(459)(458)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(1290.48,'s^-1'), n=2.90375, Ea=(139.674,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond;R_O_H] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond;R_O_H]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction238',
    reactants = ['C=CCC(O)=C=CC=O(2838)'],
    products = ['S(459)(458)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction239',
    reactants = ['C=CCC(=O)C=C=CO(2839)'],
    products = ['S(459)(458)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction244',
    reactants = ['C=CCC(C=O)C=C=O(2844)'],
    products = ['S(459)(458)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(6.21184e+10,'s^-1'), n=0.288169, Ea=(147.049,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_5_unsaturated_hexane] for rate rule [1_5_hexadiene]
Euclidian distance = 1.0
family: 6_membered_central_C-C_shift"""),
)

network(
    label = '159',
    isomers = [
        'S(459)(458)',
    ],
    reactants = [
        ('C3H5(89)(88)', 'S(780)(779)'),
    ],
    bathGas = {
        'Ne': 0.333333,
        'N2': 0.333333,
        'Ar(8)': 0.333333,
    },
)

pressureDependence(
    label = '159',
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

