species(
    label = 'CCCC[CH]OC([O])CC(13628)',
    structure = SMILES('CCCC[CH]OC([O])CC'),
    E0 = (-162.564,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,3025,407.5,1350,352.5,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,180,180,180,1600,1675.86,2834.7,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15397,'amu*angstrom^2'), symmetry=1, barrier=(3.54008,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15397,'amu*angstrom^2'), symmetry=1, barrier=(3.54008,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15397,'amu*angstrom^2'), symmetry=1, barrier=(3.54008,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15397,'amu*angstrom^2'), symmetry=1, barrier=(3.54008,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15397,'amu*angstrom^2'), symmetry=1, barrier=(3.54008,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15397,'amu*angstrom^2'), symmetry=1, barrier=(3.54008,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15397,'amu*angstrom^2'), symmetry=1, barrier=(3.54008,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15397,'amu*angstrom^2'), symmetry=1, barrier=(3.54008,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.211,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.07934,0.120801,-0.000119126,7.26186e-08,-1.95438e-11,-19376.8,41.6946], Tmin=(100,'K'), Tmax=(867.569,'K')), NASAPolynomial(coeffs=[9.73791,0.0709261,-3.28929e-05,6.35324e-09,-4.48422e-13,-21253.8,-8.95555], Tmin=(867.569,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-162.564,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(598.642,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCsJOCs) + radical(CCOJ)"""),
)

species(
    label = 'C2H5CHO(70)',
    structure = SMILES('CCC=O'),
    E0 = (-204.33,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.207559,'amu*angstrom^2'), symmetry=1, barrier=(4.77219,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208362,'amu*angstrom^2'), symmetry=1, barrier=(4.79065,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0791,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3133.67,'J/mol'), sigma=(5.35118,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=489.47 K, Pc=46.4 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.90578,0.0240644,-7.06356e-06,-9.81837e-10,5.55825e-13,-24535.9,13.5806], Tmin=(100,'K'), Tmax=(1712.49,'K')), NASAPolynomial(coeffs=[7.69109,0.0189242,-7.84934e-06,1.38273e-09,-8.99057e-14,-27060.1,-14.6647], Tmin=(1712.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-204.33,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""propanal""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CCCCC=O(1733)',
    structure = SMILES('CCCCC=O'),
    E0 = (-250.884,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,198.883,199.064],'cm^-1')),
        HinderedRotor(inertia=(0.00426375,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.518882,'amu*angstrom^2'), symmetry=1, barrier=(14.5862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.334975,'amu*angstrom^2'), symmetry=1, barrier=(9.40848,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00310926,'amu*angstrom^2'), symmetry=1, barrier=(9.40929,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.1323,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3428.95,'J/mol'), sigma=(5.98513,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=535.59 K, Pc=36.29 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63555,0.0520732,-2.77334e-05,6.90636e-09,-6.84311e-13,-30090.2,22.3653], Tmin=(100,'K'), Tmax=(2214.2,'K')), NASAPolynomial(coeffs=[15.0687,0.0278058,-1.12935e-05,1.95647e-09,-1.25427e-13,-36038.9,-53.1197], Tmin=(2214.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-250.884,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH)"""),
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
    label = 'CCCC=COC([O])CC(14052)',
    structure = SMILES('CCCC=COC([O])CC'),
    E0 = (-266.264,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    molecularWeight = (143.203,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.93356,0.113278,-8.77633e-05,3.45118e-08,-5.42951e-12,-31796.3,43.2366], Tmin=(100,'K'), Tmax=(1512.39,'K')), NASAPolynomial(coeffs=[25.4956,0.0407322,-1.58121e-05,2.79543e-09,-1.8675e-13,-40093,-100.441], Tmin=(1512.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-266.264,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCOJ)"""),
)

species(
    label = 'CCCC[CH]OC(=O)CC(14053)',
    structure = SMILES('CCCC[CH]OC(=O)CC'),
    E0 = (-366.318,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
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
    molecularWeight = (143.203,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.646131,0.109042,-9.0324e-05,4.30991e-08,-9.05759e-12,-43896.3,38.0669], Tmin=(100,'K'), Tmax=(1080.09,'K')), NASAPolynomial(coeffs=[11.1356,0.0654098,-2.97285e-05,5.69767e-09,-4.00595e-13,-46441.4,-19.681], Tmin=(1080.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-366.318,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(CCsJOC(O))"""),
)

species(
    label = 'CC[CH][O](563)',
    structure = SMILES('CC[CH][O]'),
    E0 = (133.127,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,298.357,1774.23],'cm^-1')),
        HinderedRotor(inertia=(0.129074,'amu*angstrom^2'), symmetry=1, barrier=(8.14273,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00364816,'amu*angstrom^2'), symmetry=1, barrier=(8.14268,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0791,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.1585,0.0245341,-8.42945e-06,1.83944e-10,2.32791e-13,16036.2,14.3859], Tmin=(100,'K'), Tmax=(2077.96,'K')), NASAPolynomial(coeffs=[11.8474,0.0146996,-6.30487e-06,1.09829e-09,-6.9226e-14,10937.4,-37.4679], Tmin=(2077.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(133.127,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = 'npropyl(83)',
    structure = SMILES('[CH2]CC'),
    E0 = (87.0621,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.0928812,'amu*angstrom^2'), symmetry=1, barrier=(2.13552,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.092914,'amu*angstrom^2'), symmetry=1, barrier=(2.13628,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0877,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2218.31,'J/mol'), sigma=(4.982,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.02815,0.0147023,2.4051e-05,-3.66738e-08,1.38611e-11,10512.1,12.4699], Tmin=(100,'K'), Tmax=(984.464,'K')), NASAPolynomial(coeffs=[6.16543,0.0184495,-6.79029e-06,1.23049e-09,-8.63866e-14,9095.06,-6.67607], Tmin=(984.464,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(87.0621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""npropyl""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=COC([O])CC(11164)',
    structure = SMILES('C=COC([O])CC'),
    E0 = (-183.813,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.124,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.466143,0.0660473,-3.53366e-05,-3.87649e-09,6.69865e-12,-21970.2,27.6339], Tmin=(100,'K'), Tmax=(1025.56,'K')), NASAPolynomial(coeffs=[16.9806,0.0251041,-9.7775e-06,1.82168e-09,-1.29563e-13,-26591.7,-58.4729], Tmin=(1025.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-183.813,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CCOJ)"""),
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
    label = 'CCCC[CH]OC=O(11683)',
    structure = SMILES('CCCC[CH]OC=O'),
    E0 = (-284.28,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2782.5,750,1395,475,1775,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.15,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0706182,0.0805914,-6.01144e-05,2.30215e-08,-3.60922e-12,-34045.1,31.9175], Tmin=(100,'K'), Tmax=(1479.21,'K')), NASAPolynomial(coeffs=[16.4378,0.0363318,-1.52325e-05,2.79353e-09,-1.90486e-13,-38887.2,-53.4526], Tmin=(1479.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-284.28,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdOsH) + radical(CCsJOC(O)H)"""),
)

species(
    label = 'CCC[CH]COC([O])CC(14054)',
    structure = SMILES('CCC[CH]COC([O])CC'),
    E0 = (-143.118,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,3025,407.5,1350,352.5,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,180,180,180,1600,1677.84,2833.07,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154008,'amu*angstrom^2'), symmetry=1, barrier=(3.54094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154008,'amu*angstrom^2'), symmetry=1, barrier=(3.54094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154008,'amu*angstrom^2'), symmetry=1, barrier=(3.54094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154008,'amu*angstrom^2'), symmetry=1, barrier=(3.54094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154008,'amu*angstrom^2'), symmetry=1, barrier=(3.54094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154008,'amu*angstrom^2'), symmetry=1, barrier=(3.54094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154008,'amu*angstrom^2'), symmetry=1, barrier=(3.54094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154008,'amu*angstrom^2'), symmetry=1, barrier=(3.54094,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.211,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.404437,0.106517,-8.44027e-05,4.06815e-08,-9.14432e-12,-17062.6,42.3194], Tmin=(100,'K'), Tmax=(980.605,'K')), NASAPolynomial(coeffs=[7.61465,0.073806,-3.43662e-05,6.66413e-09,-4.71776e-13,-18635.3,3.78892], Tmin=(980.605,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-143.118,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(598.642,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJCO) + radical(CCOJ)"""),
)

species(
    label = 'CCCC[CH]O[C](O)CC(14055)',
    structure = SMILES('CCCC[CH]O[C](O)CC'),
    E0 = (-183.023,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.211,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.45685,0.125147,-0.00012143,6.66313e-08,-1.53853e-11,-21820.3,42.8603], Tmin=(100,'K'), Tmax=(1023.02,'K')), NASAPolynomial(coeffs=[14.9999,0.0608011,-2.7084e-05,5.14986e-09,-3.60844e-13,-25187.4,-36.9088], Tmin=(1023.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-183.023,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(594.485,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_P) + radical(CCsJOCs)"""),
)

species(
    label = 'CCCCCO[C]([O])CC(14056)',
    structure = SMILES('CCCCCO[C]([O])CC'),
    E0 = (-137.773,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,180,180,180,180,1600,1675.09,2843.34,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154069,'amu*angstrom^2'), symmetry=1, barrier=(3.54236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154069,'amu*angstrom^2'), symmetry=1, barrier=(3.54236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154069,'amu*angstrom^2'), symmetry=1, barrier=(3.54236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154069,'amu*angstrom^2'), symmetry=1, barrier=(3.54236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154069,'amu*angstrom^2'), symmetry=1, barrier=(3.54236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154069,'amu*angstrom^2'), symmetry=1, barrier=(3.54236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154069,'amu*angstrom^2'), symmetry=1, barrier=(3.54236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154069,'amu*angstrom^2'), symmetry=1, barrier=(3.54236,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.211,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.802819,0.11524,-0.000112957,7.47055e-08,-2.28379e-11,-16406,41.9447], Tmin=(100,'K'), Tmax=(757.499,'K')), NASAPolynomial(coeffs=[6.55479,0.0763896,-3.60298e-05,7.00633e-09,-4.96108e-13,-17520.7,8.49161], Tmin=(757.499,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-137.773,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(598.642,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Cs_P)"""),
)

species(
    label = 'CC[CH]CCOC([O])CC(14057)',
    structure = SMILES('CC[CH]CCOC([O])CC'),
    E0 = (-148.561,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,3025,407.5,1350,352.5,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,180,180,180,1600,1620.64,2926.67,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15292,'amu*angstrom^2'), symmetry=1, barrier=(3.51592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15292,'amu*angstrom^2'), symmetry=1, barrier=(3.51592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15292,'amu*angstrom^2'), symmetry=1, barrier=(3.51592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15292,'amu*angstrom^2'), symmetry=1, barrier=(3.51592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15292,'amu*angstrom^2'), symmetry=1, barrier=(3.51592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15292,'amu*angstrom^2'), symmetry=1, barrier=(3.51592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15292,'amu*angstrom^2'), symmetry=1, barrier=(3.51592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15292,'amu*angstrom^2'), symmetry=1, barrier=(3.51592,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.211,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.931381,0.119591,-0.000145384,1.27818e-07,-4.86045e-11,-17701,44.2481], Tmin=(100,'K'), Tmax=(771.959,'K')), NASAPolynomial(coeffs=[3.23458,0.0812796,-3.84421e-05,7.39766e-09,-5.16962e-13,-17845.8,28.4557], Tmin=(771.959,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-148.561,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(598.642,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJCC)"""),
)

species(
    label = 'C[CH]C(O)O[CH]CCCC(14058)',
    structure = SMILES('C[CH]C(O)O[CH]CCCC'),
    E0 = (-188.367,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3615,1277.5,1000,1380,1390,370,380,2900,435,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
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
    molecularWeight = (144.211,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.42062,0.120781,-0.000108509,5.35747e-08,-1.10417e-11,-22461.5,44.5273], Tmin=(100,'K'), Tmax=(1140.9,'K')), NASAPolynomial(coeffs=[16.7364,0.0571228,-2.48148e-05,4.66933e-09,-3.25349e-13,-26604.5,-45.4632], Tmin=(1140.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-188.367,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(594.485,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCsJOCs) + radical(CCJCO)"""),
)

species(
    label = 'C[CH]CCCOC([O])CC(14059)',
    structure = SMILES('C[CH]CCCOC([O])CC'),
    E0 = (-148.574,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,3025,407.5,1350,352.5,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,180,180,180.648,1600,1606.66,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152447,'amu*angstrom^2'), symmetry=1, barrier=(3.50505,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152447,'amu*angstrom^2'), symmetry=1, barrier=(3.50505,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152447,'amu*angstrom^2'), symmetry=1, barrier=(3.50505,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152447,'amu*angstrom^2'), symmetry=1, barrier=(3.50505,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152447,'amu*angstrom^2'), symmetry=1, barrier=(3.50505,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152447,'amu*angstrom^2'), symmetry=1, barrier=(3.50505,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152447,'amu*angstrom^2'), symmetry=1, barrier=(3.50505,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152447,'amu*angstrom^2'), symmetry=1, barrier=(3.50505,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.211,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.939285,0.118731,-0.000138667,1.17192e-07,-4.36917e-11,-17701.1,44.4509], Tmin=(100,'K'), Tmax=(757.199,'K')), NASAPolynomial(coeffs=[4.22823,0.0795076,-3.73405e-05,7.17987e-09,-5.02267e-13,-18141.8,23.2153], Tmin=(757.199,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-148.574,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(598.642,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(CCOJ)"""),
)

species(
    label = '[CH2]CC(O)O[CH]CCCC(14060)',
    structure = SMILES('[CH2]CC(O)O[CH]CCCC'),
    E0 = (-183.023,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,3000,3100,440,815,1455,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
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
    molecularWeight = (144.211,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.45685,0.125147,-0.00012143,6.66313e-08,-1.53853e-11,-21820.3,43.959], Tmin=(100,'K'), Tmax=(1023.02,'K')), NASAPolynomial(coeffs=[14.9999,0.0608011,-2.7084e-05,5.14986e-09,-3.60844e-13,-25187.4,-35.8102], Tmin=(1023.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-183.023,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(594.485,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(CCsJOCs)"""),
)

species(
    label = 'C[CH]C([O])OCCCCC(14061)',
    structure = SMILES('C[CH]C([O])OCCCCC'),
    E0 = (-143.118,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,3025,407.5,1350,352.5,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,180,180,180,1600,1677.84,2833.07,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154008,'amu*angstrom^2'), symmetry=1, barrier=(3.54094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154008,'amu*angstrom^2'), symmetry=1, barrier=(3.54094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154008,'amu*angstrom^2'), symmetry=1, barrier=(3.54094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154008,'amu*angstrom^2'), symmetry=1, barrier=(3.54094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154008,'amu*angstrom^2'), symmetry=1, barrier=(3.54094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154008,'amu*angstrom^2'), symmetry=1, barrier=(3.54094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154008,'amu*angstrom^2'), symmetry=1, barrier=(3.54094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154008,'amu*angstrom^2'), symmetry=1, barrier=(3.54094,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.211,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.404449,0.106517,-8.44031e-05,4.06819e-08,-9.14445e-12,-17062.6,42.3194], Tmin=(100,'K'), Tmax=(980.594,'K')), NASAPolynomial(coeffs=[7.61462,0.073806,-3.43663e-05,6.66414e-09,-4.71777e-13,-18635.3,3.78911], Tmin=(980.594,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-143.118,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(598.642,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CCJCO)"""),
)

species(
    label = 'CCC[CH][CH]OC(O)CC(14062)',
    structure = SMILES('CCC[CH][CH]OC(O)CC'),
    E0 = (-188.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.211,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.42062,0.120781,-0.000108509,5.35747e-08,-1.10417e-11,-22461.5,44.5273], Tmin=(100,'K'), Tmax=(1140.9,'K')), NASAPolynomial(coeffs=[16.7364,0.0571228,-2.48148e-05,4.66933e-09,-3.25349e-13,-26604.5,-45.4632], Tmin=(1140.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-188.367,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(594.485,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJCO) + radical(CCsJOCs)"""),
)

species(
    label = '[CH2]CCCCOC([O])CC(14063)',
    structure = SMILES('[CH2]CCCCOC([O])CC'),
    E0 = (-137.773,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,180,1599.12,1600,2918.37,3200],'cm^-1')),
        HinderedRotor(inertia=(0.151905,'amu*angstrom^2'), symmetry=1, barrier=(3.49259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151905,'amu*angstrom^2'), symmetry=1, barrier=(3.49259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151905,'amu*angstrom^2'), symmetry=1, barrier=(3.49259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151905,'amu*angstrom^2'), symmetry=1, barrier=(3.49259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151905,'amu*angstrom^2'), symmetry=1, barrier=(3.49259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151905,'amu*angstrom^2'), symmetry=1, barrier=(3.49259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151905,'amu*angstrom^2'), symmetry=1, barrier=(3.49259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151905,'amu*angstrom^2'), symmetry=1, barrier=(3.49259,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.211,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.802819,0.11524,-0.000112957,7.47055e-08,-2.28379e-11,-16406,43.0433], Tmin=(100,'K'), Tmax=(757.499,'K')), NASAPolynomial(coeffs=[6.55479,0.0763896,-3.60298e-05,7.00633e-09,-4.96108e-13,-17520.7,9.59023], Tmin=(757.499,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-137.773,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(598.642,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(CCOJ)"""),
)

species(
    label = '[CH2]CC([O])OCCCCC(14064)',
    structure = SMILES('[CH2]CC([O])OCCCCC'),
    E0 = (-137.773,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,180,1599.12,1600,2918.37,3200],'cm^-1')),
        HinderedRotor(inertia=(0.151905,'amu*angstrom^2'), symmetry=1, barrier=(3.49259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151905,'amu*angstrom^2'), symmetry=1, barrier=(3.49259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151905,'amu*angstrom^2'), symmetry=1, barrier=(3.49259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151905,'amu*angstrom^2'), symmetry=1, barrier=(3.49259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151905,'amu*angstrom^2'), symmetry=1, barrier=(3.49259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151905,'amu*angstrom^2'), symmetry=1, barrier=(3.49259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151905,'amu*angstrom^2'), symmetry=1, barrier=(3.49259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151905,'amu*angstrom^2'), symmetry=1, barrier=(3.49259,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.211,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.802819,0.11524,-0.000112957,7.47055e-08,-2.28379e-11,-16406,43.0433], Tmin=(100,'K'), Tmax=(757.499,'K')), NASAPolynomial(coeffs=[6.55479,0.0763896,-3.60298e-05,7.00633e-09,-4.96108e-13,-17520.7,9.59023], Tmin=(757.499,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-137.773,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(598.642,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(CCOJ)"""),
)

species(
    label = 'CC[CH]C[CH]OC(O)CC(14065)',
    structure = SMILES('CC[CH]C[CH]OC(O)CC'),
    E0 = (-193.811,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.211,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.34134,0.125948,-0.000138071,9.40245e-08,-2.75997e-11,-23125.3,44.3303], Tmin=(100,'K'), Tmax=(810.376,'K')), NASAPolynomial(coeffs=[10.382,0.0680799,-3.09536e-05,5.89973e-09,-4.12354e-13,-25025.3,-9.76281], Tmin=(810.376,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-193.811,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(594.485,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(CCsJOCs)"""),
)

species(
    label = 'C[CH]CC[CH]OC(O)CC(14066)',
    structure = SMILES('C[CH]CC[CH]OC(O)CC'),
    E0 = (-193.823,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.211,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.33032,0.124904,-0.000130884,8.29857e-08,-2.25644e-11,-23126.2,44.4627], Tmin=(100,'K'), Tmax=(872.132,'K')), NASAPolynomial(coeffs=[11.4955,0.0660788,-2.97085e-05,5.64601e-09,-3.94545e-13,-25363.4,-15.6595], Tmin=(872.132,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-193.823,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(594.485,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCsJOCs) + radical(RCCJC)"""),
)

species(
    label = '[CH2]CCC[CH]OC(O)CC(14067)',
    structure = SMILES('[CH2]CCC[CH]OC(O)CC'),
    E0 = (-183.023,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.211,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.45685,0.125147,-0.00012143,6.66313e-08,-1.53853e-11,-21820.3,43.959], Tmin=(100,'K'), Tmax=(1023.02,'K')), NASAPolynomial(coeffs=[14.9999,0.0608011,-2.7084e-05,5.14986e-09,-3.60844e-13,-25187.4,-35.8102], Tmin=(1023.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-183.023,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(594.485,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(CCsJOCs)"""),
)

species(
    label = 'CCCC[CH][O](1728)',
    structure = SMILES('CCCC[CH][O]'),
    E0 = (78.0912,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,284.58,286.601,291.211,2078.87],'cm^-1')),
        HinderedRotor(inertia=(0.00200967,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10398,'amu*angstrom^2'), symmetry=1, barrier=(6.1119,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103299,'amu*angstrom^2'), symmetry=1, barrier=(6.13902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0993848,'amu*angstrom^2'), symmetry=1, barrier=(6.07184,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.1323,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20776,0.0668884,-7.57546e-05,6.29928e-08,-2.32085e-11,9487.51,24.9863], Tmin=(100,'K'), Tmax=(767.626,'K')), NASAPolynomial(coeffs=[4.00991,0.045565,-2.09523e-05,3.99079e-09,-2.77682e-13,9255.35,13.4985], Tmin=(767.626,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(78.0912,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = 'CCCCCOC(=O)CC(14068)',
    structure = SMILES('CCCCCOC(=O)CC'),
    E0 = (-560.245,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.211,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.148337,0.0986254,-6.17082e-05,1.98561e-08,-2.76929e-12,-67240.2,36.764], Tmin=(100,'K'), Tmax=(1500.59,'K')), NASAPolynomial(coeffs=[11.9714,0.0663189,-2.94144e-05,5.50902e-09,-3.79052e-13,-70877.5,-26.6254], Tmin=(1500.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-560.245,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(598.642,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsOs)"""),
)

species(
    label = 'CCCC=COC(O)CC(14069)',
    structure = SMILES('CCCC=COC(O)CC'),
    E0 = (-491.97,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.211,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.94657,0.11235,-6.4776e-05,-1.21211e-09,9.76748e-12,-58940,42.3946], Tmin=(100,'K'), Tmax=(1013.55,'K')), NASAPolynomial(coeffs=[25.3555,0.0419202,-1.57733e-05,2.88344e-09,-2.03103e-13,-66391.3,-99.1457], Tmin=(1013.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-491.97,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(598.642,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH)"""),
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
    label = 'CCCC[CH]OC(C)[O](11951)',
    structure = SMILES('CCCC[CH]OC(C)[O]'),
    E0 = (-138.784,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3025,407.5,1350,352.5,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    molecularWeight = (130.185,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4082.6,'J/mol'), sigma=(7.30891,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=637.69 K, Pc=23.73 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.476561,0.10657,-0.000107362,6.74296e-08,-1.86955e-11,-16537.7,37.2826], Tmin=(100,'K'), Tmax=(844.446,'K')), NASAPolynomial(coeffs=[8.82972,0.0624894,-2.90637e-05,5.61726e-09,-3.96475e-13,-18109.5,-6.04165], Tmin=(844.446,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-138.784,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CCsJOCs)"""),
)

species(
    label = 'CCC[CH]OC([O])CC(14070)',
    structure = SMILES('CCC[CH]OC([O])CC'),
    E0 = (-138.784,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3025,407.5,1350,352.5,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    molecularWeight = (130.185,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.476561,0.10657,-0.000107362,6.74296e-08,-1.86955e-11,-16537.7,37.2826], Tmin=(100,'K'), Tmax=(844.446,'K')), NASAPolynomial(coeffs=[8.82972,0.0624894,-2.90637e-05,5.61726e-09,-3.96475e-13,-18109.5,-6.04165], Tmin=(844.446,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-138.784,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCsJOCs) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C(CCC)OC([O])CC(14071)',
    structure = SMILES('[CH2]C(CCC)OC([O])CC'),
    E0 = (-143.061,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    molecularWeight = (144.211,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.3167,0.12665,-0.000141075,1.01138e-07,-3.18179e-11,-17023.6,42.5009], Tmin=(100,'K'), Tmax=(753.812,'K')), NASAPolynomial(coeffs=[8.8055,0.0729365,-3.41875e-05,6.60436e-09,-4.64914e-13,-18549.6,-3.47193], Tmin=(753.812,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-143.061,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(598.642,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CJC(C)OC)"""),
)

species(
    label = 'CCCCC1OC(CC)O1(13631)',
    structure = SMILES('CCCCC1OC(CC)O1'),
    E0 = (-420.193,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.211,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.999446,0.0924285,-2.92335e-05,-1.81291e-08,1.04133e-11,-50342.9,37.7069], Tmin=(100,'K'), Tmax=(1132.15,'K')), NASAPolynomial(coeffs=[19.1164,0.0546903,-2.33972e-05,4.43982e-09,-3.12864e-13,-57034.1,-71.272], Tmin=(1132.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-420.193,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(606.956,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane)"""),
)

species(
    label = '[CH]CCCC(141)',
    structure = SMILES('[CH]CCCC'),
    E0 = (280.131,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49126,0.0476292,-1.77617e-05,-4.75447e-09,3.68063e-12,33788.6,21.8622], Tmin=(100,'K'), Tmax=(1127.24,'K')), NASAPolynomial(coeffs=[10.1002,0.0302049,-1.20402e-05,2.19072e-09,-1.50454e-13,31013.8,-24.4008], Tmin=(1127.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(280.131,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = 'CCC([O])[O](3172)',
    structure = SMILES('CCC([O])[O]'),
    E0 = (-29.172,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,180,2001.39,2003.47],'cm^-1')),
        HinderedRotor(inertia=(0.164747,'amu*angstrom^2'), symmetry=1, barrier=(3.78786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166214,'amu*angstrom^2'), symmetry=1, barrier=(3.82159,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0785,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.47585,0.0406137,-5.60044e-05,6.17463e-08,-2.68606e-11,-3460.67,19.7168], Tmin=(100,'K'), Tmax=(802.392,'K')), NASAPolynomial(coeffs=[-0.627181,0.0384738,-1.9086e-05,3.72267e-09,-2.60878e-13,-2395.84,37.5363], Tmin=(802.392,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-29.172,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CCOJ)"""),
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
    label = 'CC[CH]O[CH]CCCC(10789)',
    structure = SMILES('CC[CH]O[CH]CCCC'),
    E0 = (-7.58188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    molecularWeight = (128.212,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.354,0.119374,-0.000113611,6.0383e-08,-1.3334e-11,-720.404,38.0808], Tmin=(100,'K'), Tmax=(1075.49,'K')), NASAPolynomial(coeffs=[16.087,0.0545067,-2.31391e-05,4.30168e-09,-2.97753e-13,-4471.9,-47.3311], Tmin=(1075.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-7.58188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCsJOCs) + radical(CCsJOCs)"""),
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
    E0 = (-162.564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-47.7361,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-108.02,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-81.8706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-62.2748,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-148.373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-34.3337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-4.5274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (27.7456,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-8.60647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-82.8049,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-102.319,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-99.3426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-96.3824,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-111.686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-83.3117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-91.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-82.8586,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-55.6628,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-55.7169,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (211.218,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-99.1636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-137.591,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (281.078,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (281.078,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (15.5662,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-154.279,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (250.959,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (235.423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['CCCC[CH]OC([O])CC(13628)'],
    products = ['C2H5CHO(70)', 'CCCCC=O(1733)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(3)', 'CCCC=COC([O])CC(14052)'],
    products = ['CCCC[CH]OC([O])CC(13628)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.72e+08,'cm^3/(mol*s)'), n=1.477, Ea=(6.73624,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2825 used for Cds-CsH_Cds-OsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-OsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'CCCC[CH]OC(=O)CC(14053)'],
    products = ['CCCC[CH]OC([O])CC(13628)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.83701,'m^3/(mol*s)'), n=1.71338, Ea=(46.5052,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;HJ] + [CO-NdNd_O;YJ] for rate rule [CO-NdNd_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CC[CH][O](563)', 'CCCCC=O(1733)'],
    products = ['CCCC[CH]OC([O])CC(13628)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4e+09,'cm^3/(mol*s)'), n=1.39, Ea=(35.8862,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Od_CO-CsH;YJ] for rate rule [Od_CO-CsH;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['npropyl(83)', 'C=COC([O])CC(11164)'],
    products = ['CCCC[CH]OC([O])CC(13628)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1071,'cm^3/(mol*s)'), n=2.72, Ea=(34.4762,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2780 used for Cds-HH_Cds-OsH;CsJ-CsHH
Exact match found for rate rule [Cds-HH_Cds-OsH;CsJ-CsHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C2H5(29)', 'CCCC[CH]OC=O(11683)'],
    products = ['CCCC[CH]OC([O])CC(13628)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7.94e+10,'cm^3/(mol*s)'), n=0, Ea=(28.0328,'kJ/mol'), T0=(1,'K'), Tmin=(333,'K'), Tmax=(363,'K'), comment="""Estimated using template [CO_O;CsJ-CsHH] for rate rule [CO-NdH_O;CsJ-CsHH]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CCC[CH]COC([O])CC(14054)'],
    products = ['CCCC[CH]OC([O])CC(13628)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(5.4e-20,'s^-1'), n=9.13, Ea=(108.784,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 341 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CCCC[CH]OC([O])CC(13628)'],
    products = ['CCCC[CH]O[C](O)CC(14055)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.43381e+07,'s^-1'), n=1.70481, Ea=(158.036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_NonDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CCCCCO[C]([O])CC(14056)'],
    products = ['CCCC[CH]OC([O])CC(13628)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.53164e+09,'s^-1'), n=1.09, Ea=(165.519,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_O;Y_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CC[CH]CCOC([O])CC(14057)'],
    products = ['CCCC[CH]OC([O])CC(13628)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.73898e+09,'s^-1'), n=0.87, Ea=(139.955,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_1H] for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C[CH]C(O)O[CH]CCCC(14058)'],
    products = ['CCCC[CH]OC([O])CC(13628)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.71,'s^-1'), n=3.021, Ea=(105.562,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 319 used for R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C[CH]CCCOC([O])CC(14059)'],
    products = ['CCCC[CH]OC([O])CC(13628)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(0.00310895,'s^-1'), n=4.065, Ea=(46.2541,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_1H] for rate rule [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]CC(O)O[CH]CCCC(14060)'],
    products = ['CCCC[CH]OC([O])CC(13628)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C[CH]C([O])OCCCCC(14061)'],
    products = ['CCCC[CH]OC([O])CC(13628)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.00228,'s^-1'), n=3.95, Ea=(46.7353,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC] for rate rule [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CCCC[CH]OC([O])CC(13628)'],
    products = ['CCC[CH][CH]OC(O)CC(14062)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.19599e+09,'s^-1'), n=0.63, Ea=(50.8774,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5HJ_3;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]CCCCOC([O])CC(14063)'],
    products = ['CCCC[CH]OC([O])CC(13628)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(30977.5,'s^-1'), n=1.81667, Ea=(54.4617,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_2H;Cs_H_out_1H] for rate rule [R5H_CCC;C_rad_out_2H;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]CC([O])OCCCCC(14064)'],
    products = ['CCCC[CH]OC([O])CC(13628)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(262000,'s^-1'), n=1.62, Ea=(46.4424,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;C_rad_out_2H;Cs_H_out_H/NonDeC] for rate rule [R5H_SSSS;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CCCC[CH]OC([O])CC(13628)'],
    products = ['CC[CH]C[CH]OC(O)CC(14065)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.12e+09,'s^-1'), n=0, Ea=(79.7052,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6Hall;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R6HJ_3;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CCCC[CH]OC([O])CC(13628)'],
    products = ['C[CH]CC[CH]OC(O)CC(14066)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.9e+08,'s^-1'), n=0, Ea=(106.901,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R7Hall;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R7HJ_3;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['CCCC[CH]OC([O])CC(13628)'],
    products = ['[CH2]CCC[CH]OC(O)CC(14067)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.18263e+09,'s^-1'), n=0.595, Ea=(106.847,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;O_rad_out;Cs_H_out_2H] for rate rule [R8Hall;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CC[CH][O](563)', 'CCCC[CH][O](1728)'],
    products = ['CCCC[CH]OC([O])CC(13628)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['CCCC[CH]OC([O])CC(13628)'],
    products = ['CCCCCOC(=O)CC(14068)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CCCC[CH]OC([O])CC(13628)'],
    products = ['CCCC=COC(O)CC(14069)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CH2(S)(23)', 'CCCC[CH]OC(C)[O](11951)'],
    products = ['CCCC[CH]OC([O])CC(13628)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['CH2(S)(23)', 'CCC[CH]OC([O])CC(14070)'],
    products = ['CCCC[CH]OC([O])CC(13628)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(CCC)OC([O])CC(14071)'],
    products = ['CCCC[CH]OC([O])CC(13628)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CCCC[CH]OC([O])CC(13628)'],
    products = ['CCCCC1OC(CC)O1(13631)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]CCCC(141)', 'CCC([O])[O](3172)'],
    products = ['CCCC[CH]OC([O])CC(13628)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2711.41,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction29',
    reactants = ['O(4)', 'CC[CH]O[CH]CCCC(10789)'],
    products = ['CCCC[CH]OC([O])CC(13628)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/CsO;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

network(
    label = '3030',
    isomers = [
        'CCCC[CH]OC([O])CC(13628)',
    ],
    reactants = [
        ('C2H5CHO(70)', 'CCCCC=O(1733)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '3030',
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

