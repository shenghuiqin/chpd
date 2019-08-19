species(
    label = 'CCCC[CH]OC([O])CCC=O(20308)',
    structure = SMILES('CCCC[CH]OC([O])CCC=O'),
    E0 = (-269.143,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,2782.5,750,1395,475,1775,1000,3025,407.5,1350,352.5,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,475.336,683.338,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156878,'amu*angstrom^2'), symmetry=1, barrier=(3.60693,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156878,'amu*angstrom^2'), symmetry=1, barrier=(3.60693,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156878,'amu*angstrom^2'), symmetry=1, barrier=(3.60693,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156878,'amu*angstrom^2'), symmetry=1, barrier=(3.60693,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156878,'amu*angstrom^2'), symmetry=1, barrier=(3.60693,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156878,'amu*angstrom^2'), symmetry=1, barrier=(3.60693,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156878,'amu*angstrom^2'), symmetry=1, barrier=(3.60693,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156878,'amu*angstrom^2'), symmetry=1, barrier=(3.60693,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156878,'amu*angstrom^2'), symmetry=1, barrier=(3.60693,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.13372,0.132593,-9.493e-05,-4.91092e-08,9.34195e-11,-32205.4,46.4069], Tmin=(100,'K'), Tmax=(500.891,'K')), NASAPolynomial(coeffs=[9.49897,0.0807158,-3.84987e-05,7.44618e-09,-5.22582e-13,-33684.9,-1.67507], Tmin=(500.891,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-269.143,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(644.372,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CCsJOCs)"""),
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
    label = 'O=CCCC=O(5767)',
    structure = SMILES('O=CCCC=O'),
    E0 = (-309.903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,180],'cm^-1')),
        HinderedRotor(inertia=(0.18601,'amu*angstrom^2'), symmetry=1, barrier=(4.27673,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.186498,'amu*angstrom^2'), symmetry=1, barrier=(4.28796,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.186203,'amu*angstrom^2'), symmetry=1, barrier=(4.28117,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3653.08,'J/mol'), sigma=(5.8998,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=570.60 K, Pc=40.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98042,0.0494518,-6.06096e-05,5.52155e-08,-2.1639e-11,-37204.9,20.2227], Tmin=(100,'K'), Tmax=(774.657,'K')), NASAPolynomial(coeffs=[2.99129,0.0355848,-1.7014e-05,3.28724e-09,-2.30123e-13,-37102,17.2786], Tmin=(774.657,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-309.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + group(Cds-OdCsH)"""),
)

species(
    label = 'CCCC[CH]OC1CCC([O])O1(20789)',
    structure = SMILES('CCCC[CH]OC1CCC([O])O1'),
    E0 = (-281.006,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.05946,0.121685,-9.0678e-05,3.66617e-08,-6.13282e-12,-33569.6,45.1046], Tmin=(100,'K'), Tmax=(1397.43,'K')), NASAPolynomial(coeffs=[20.3722,0.0574759,-2.17553e-05,3.78053e-09,-2.50319e-13,-39838.8,-70.6211], Tmin=(1397.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-281.006,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(656.843,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(CCsJOCs) + radical(CCOJ)"""),
)

species(
    label = 'CCCCC1OC([O])CCC1[O](20761)',
    structure = SMILES('CCCCC1OC([O])CCC1[O]'),
    E0 = (-268.562,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.62103,0.108166,-4.93034e-05,-5.41078e-09,8.26275e-12,-32084.9,42.4225], Tmin=(100,'K'), Tmax=(1045.33,'K')), NASAPolynomial(coeffs=[18.107,0.0613052,-2.3143e-05,4.10601e-09,-2.7945e-13,-37773.6,-61.1104], Tmin=(1045.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-268.562,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(665.158,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + ring(Oxane) + radical(CC(C)OJ) + radical(CCOJ)"""),
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
    label = 'CCCC=COC([O])CCC=O(20790)',
    structure = SMILES('CCCC=COC([O])CCC=O'),
    E0 = (-372.844,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,2782.5,750,1395,475,1775,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,180,1582.97,1600,2911.13,3200],'cm^-1')),
        HinderedRotor(inertia=(0.151091,'amu*angstrom^2'), symmetry=1, barrier=(3.47388,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151091,'amu*angstrom^2'), symmetry=1, barrier=(3.47388,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151091,'amu*angstrom^2'), symmetry=1, barrier=(3.47388,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151091,'amu*angstrom^2'), symmetry=1, barrier=(3.47388,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151091,'amu*angstrom^2'), symmetry=1, barrier=(3.47388,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151091,'amu*angstrom^2'), symmetry=1, barrier=(3.47388,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151091,'amu*angstrom^2'), symmetry=1, barrier=(3.47388,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151091,'amu*angstrom^2'), symmetry=1, barrier=(3.47388,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (171.214,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.9054,0.129068,-0.000110355,4.99702e-08,-9.35972e-12,-44629.7,47.385], Tmin=(100,'K'), Tmax=(1249.52,'K')), NASAPolynomial(coeffs=[20.0772,0.0586975,-2.58785e-05,4.89889e-09,-3.42068e-13,-50123.3,-63.5655], Tmin=(1249.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-372.844,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(623.585,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(CCOJ)"""),
)

species(
    label = 'CCCC[CH]OC(=O)CCC=O(20791)',
    structure = SMILES('CCCC[CH]OC(=O)CCC=O'),
    E0 = (-472.897,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,2782.5,750,1395,475,1775,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
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
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (171.214,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.78684,0.138002,-0.000156416,1.12142e-07,-3.50255e-11,-56677.8,46.4522], Tmin=(100,'K'), Tmax=(760.746,'K')), NASAPolynomial(coeffs=[9.7764,0.0772019,-3.65306e-05,7.0805e-09,-4.99261e-13,-58437.1,-6.17157], Tmin=(760.746,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-472.897,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + group(Cds-OdCsH) + radical(CCsJOC(O))"""),
)

species(
    label = '[O][CH]CCC=O(5764)',
    structure = SMILES('[O][CH]CCC=O'),
    E0 = (19.0721,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,3025,407.5,1350,352.5,180,180,2092.49],'cm^-1')),
        HinderedRotor(inertia=(0.166661,'amu*angstrom^2'), symmetry=1, barrier=(3.83187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164989,'amu*angstrom^2'), symmetry=1, barrier=(3.79342,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166144,'amu*angstrom^2'), symmetry=1, barrier=(3.81999,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56305,0.0641173,-0.000107371,1.07453e-07,-4.10375e-11,2371.43,23.494], Tmin=(100,'K'), Tmax=(845.992,'K')), NASAPolynomial(coeffs=[1.74033,0.0387928,-1.90537e-05,3.64333e-09,-2.50242e-13,3217.68,27.8472], Tmin=(845.992,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(19.0721,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + radical(CCsJOH) + radical(CCOJ)"""),
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
    label = 'C=COC([O])CCC=O(18636)',
    structure = SMILES('C=COC([O])CCC=O'),
    E0 = (-290.393,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.123618,0.088702,-7.99792e-05,3.74555e-08,-7.13778e-12,-34776.2,34.0275], Tmin=(100,'K'), Tmax=(1244.26,'K')), NASAPolynomial(coeffs=[16.4867,0.0353033,-1.56048e-05,2.96387e-09,-2.07589e-13,-38909.7,-49.7378], Tmin=(1244.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-290.393,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CCOJ)"""),
)

species(
    label = 'CH2CH2CHO(560)',
    structure = SMILES('[CH2]CC=O'),
    E0 = (11.2619,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000],'cm^-1')),
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
    label = 'CCC[CH]COC([O])CCC=O(20792)',
    structure = SMILES('CCC[CH]COC([O])CCC=O'),
    E0 = (-249.697,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,2782.5,750,1395,475,1775,1000,3025,407.5,1350,352.5,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,439.914,724.105,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156375,'amu*angstrom^2'), symmetry=1, barrier=(3.59538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156375,'amu*angstrom^2'), symmetry=1, barrier=(3.59538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156375,'amu*angstrom^2'), symmetry=1, barrier=(3.59538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156375,'amu*angstrom^2'), symmetry=1, barrier=(3.59538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156375,'amu*angstrom^2'), symmetry=1, barrier=(3.59538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156375,'amu*angstrom^2'), symmetry=1, barrier=(3.59538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156375,'amu*angstrom^2'), symmetry=1, barrier=(3.59538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156375,'amu*angstrom^2'), symmetry=1, barrier=(3.59538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156375,'amu*angstrom^2'), symmetry=1, barrier=(3.59538,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.436862,0.116884,-4.66189e-05,-1.2197e-07,1.42279e-10,-29891.7,47.013], Tmin=(100,'K'), Tmax=(469.866,'K')), NASAPolynomial(coeffs=[6.97981,0.0843473,-4.04404e-05,7.87407e-09,-5.56051e-13,-30926.4,13.239], Tmin=(469.866,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-249.697,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(644.372,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCO) + radical(CCOJ)"""),
)

species(
    label = 'CCCC[CH]O[C](O)CCC=O(20793)',
    structure = SMILES('CCCC[CH]O[C](O)CCC=O'),
    E0 = (-289.602,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.45508,0.15245,-0.000181743,1.27938e-07,-3.78114e-11,-34608.1,50.7321], Tmin=(100,'K'), Tmax=(812.813,'K')), NASAPolynomial(coeffs=[13.8563,0.0721863,-3.36331e-05,6.46964e-09,-4.54066e-13,-37259.9,-24.5817], Tmin=(812.813,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-289.602,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(640.214,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCsJOCs) + radical(Cs_P)"""),
)

species(
    label = 'CCCCCO[C]([O])CCC=O(20794)',
    structure = SMILES('CCCCCO[C]([O])CCC=O'),
    E0 = (-244.353,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2759.09,2768.18,2777.27,2786.36,2795.45,2804.55,2813.64,2822.73,2831.82,2840.91,2850,1425,1430,1435,1440,1445,1450,1225,1235,1245,1255,1265,1275,1270,1284,1298,1312,1326,1340,700,720,740,760,780,800,300,320,340,360,380,400,2782.5,750,1395,475,1775,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,180,180,180,464.425,700.603,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156809,'amu*angstrom^2'), symmetry=1, barrier=(3.60534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156809,'amu*angstrom^2'), symmetry=1, barrier=(3.60534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156809,'amu*angstrom^2'), symmetry=1, barrier=(3.60534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156809,'amu*angstrom^2'), symmetry=1, barrier=(3.60534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156809,'amu*angstrom^2'), symmetry=1, barrier=(3.60534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156809,'amu*angstrom^2'), symmetry=1, barrier=(3.60534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156809,'amu*angstrom^2'), symmetry=1, barrier=(3.60534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156809,'amu*angstrom^2'), symmetry=1, barrier=(3.60534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156809,'amu*angstrom^2'), symmetry=1, barrier=(3.60534,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.07673,0.146513,-0.000190785,1.6445e-07,-6.02483e-11,-29182.4,50.7605], Tmin=(100,'K'), Tmax=(763.169,'K')), NASAPolynomial(coeffs=[6.73701,0.0853349,-4.10908e-05,7.96005e-09,-5.58035e-13,-30091.4,13.4801], Tmin=(763.169,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-244.353,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(644.372,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCOJ) + radical(Cs_P)"""),
)

species(
    label = 'CC[CH]CCOC([O])CCC=O(20795)',
    structure = SMILES('CC[CH]CCOC([O])CCC=O'),
    E0 = (-255.141,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,2782.5,750,1395,475,1775,1000,3025,407.5,1350,352.5,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,328.169,871.411,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155094,'amu*angstrom^2'), symmetry=1, barrier=(3.56592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155094,'amu*angstrom^2'), symmetry=1, barrier=(3.56592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155094,'amu*angstrom^2'), symmetry=1, barrier=(3.56592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155094,'amu*angstrom^2'), symmetry=1, barrier=(3.56592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155094,'amu*angstrom^2'), symmetry=1, barrier=(3.56592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155094,'amu*angstrom^2'), symmetry=1, barrier=(3.56592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155094,'amu*angstrom^2'), symmetry=1, barrier=(3.56592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155094,'amu*angstrom^2'), symmetry=1, barrier=(3.56592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155094,'amu*angstrom^2'), symmetry=1, barrier=(3.56592,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.93812,0.147244,-0.000208009,1.93414e-07,-7.32813e-11,-30488.7,52.1336], Tmin=(100,'K'), Tmax=(804.959,'K')), NASAPolynomial(coeffs=[2.96165,0.0910685,-4.40199e-05,8.47885e-09,-5.89806e-13,-30246.4,35.9631], Tmin=(804.959,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-255.141,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(644.372,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJCC) + radical(CCOJ)"""),
)

species(
    label = 'CCCC[CH]OC(O)[CH]CC=O(20796)',
    structure = SMILES('CCCC[CH]OC(O)[CH]CC=O'),
    E0 = (-294.946,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,2782.5,750,1395,475,1775,1000,180,180,180,180,1600,1662.01,2841.74,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153965,'amu*angstrom^2'), symmetry=1, barrier=(3.53997,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153965,'amu*angstrom^2'), symmetry=1, barrier=(3.53997,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153965,'amu*angstrom^2'), symmetry=1, barrier=(3.53997,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153965,'amu*angstrom^2'), symmetry=1, barrier=(3.53997,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153965,'amu*angstrom^2'), symmetry=1, barrier=(3.53997,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153965,'amu*angstrom^2'), symmetry=1, barrier=(3.53997,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153965,'amu*angstrom^2'), symmetry=1, barrier=(3.53997,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153965,'amu*angstrom^2'), symmetry=1, barrier=(3.53997,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153965,'amu*angstrom^2'), symmetry=1, barrier=(3.53997,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153965,'amu*angstrom^2'), symmetry=1, barrier=(3.53997,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.18668,0.145166,-0.000157829,9.95711e-08,-2.6479e-11,-35258.9,51.5793], Tmin=(100,'K'), Tmax=(896.186,'K')), NASAPolynomial(coeffs=[14.6104,0.070194,-3.23425e-05,6.22135e-09,-4.37899e-13,-38269.5,-27.6157], Tmin=(896.186,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-294.946,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(640.214,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCO) + radical(CCsJOCs)"""),
)

species(
    label = 'C[CH]CCCOC([O])CCC=O(20797)',
    structure = SMILES('C[CH]CCCOC([O])CCC=O'),
    E0 = (-255.153,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,2782.5,750,1395,475,1775,1000,3025,407.5,1350,352.5,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,301.641,889.788,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154548,'amu*angstrom^2'), symmetry=1, barrier=(3.55337,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154548,'amu*angstrom^2'), symmetry=1, barrier=(3.55337,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154548,'amu*angstrom^2'), symmetry=1, barrier=(3.55337,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154548,'amu*angstrom^2'), symmetry=1, barrier=(3.55337,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154548,'amu*angstrom^2'), symmetry=1, barrier=(3.55337,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154548,'amu*angstrom^2'), symmetry=1, barrier=(3.55337,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154548,'amu*angstrom^2'), symmetry=1, barrier=(3.55337,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154548,'amu*angstrom^2'), symmetry=1, barrier=(3.55337,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154548,'amu*angstrom^2'), symmetry=1, barrier=(3.55337,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.96719,0.146663,-0.000202409,1.84453e-07,-6.9178e-11,-30487.9,52.4105], Tmin=(100,'K'), Tmax=(799.252,'K')), NASAPolynomial(coeffs=[4.01212,0.0891925,-4.28551e-05,8.24556e-09,-5.73789e-13,-30563.9,30.4076], Tmin=(799.252,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-255.153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(644.372,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJC) + radical(CCOJ)"""),
)

species(
    label = 'CCCC[CH]OC(O)C[CH]C=O(20798)',
    structure = SMILES('CCCC[CH]OC(O)CC=C[O]'),
    E0 = (-351.597,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.98993,0.146676,-0.000139828,6.93488e-08,-1.37945e-11,-42029.8,51.6476], Tmin=(100,'K'), Tmax=(1210.25,'K')), NASAPolynomial(coeffs=[26.0316,0.0507565,-2.09432e-05,3.86096e-09,-2.66661e-13,-49054.4,-93.9022], Tmin=(1210.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-351.597,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(644.372,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCsJOCs) + radical(C=COJ)"""),
)

species(
    label = 'CCCCCOC([O])[CH]CC=O(20799)',
    structure = SMILES('CCCCCOC([O])[CH]CC=O'),
    E0 = (-249.697,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,2782.5,750,1395,475,1775,1000,3025,407.5,1350,352.5,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,439.914,724.105,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156375,'amu*angstrom^2'), symmetry=1, barrier=(3.59538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156375,'amu*angstrom^2'), symmetry=1, barrier=(3.59538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156375,'amu*angstrom^2'), symmetry=1, barrier=(3.59538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156375,'amu*angstrom^2'), symmetry=1, barrier=(3.59538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156375,'amu*angstrom^2'), symmetry=1, barrier=(3.59538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156375,'amu*angstrom^2'), symmetry=1, barrier=(3.59538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156375,'amu*angstrom^2'), symmetry=1, barrier=(3.59538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156375,'amu*angstrom^2'), symmetry=1, barrier=(3.59538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156375,'amu*angstrom^2'), symmetry=1, barrier=(3.59538,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.436862,0.116884,-4.66189e-05,-1.2197e-07,1.42279e-10,-29891.7,47.013], Tmin=(100,'K'), Tmax=(469.866,'K')), NASAPolynomial(coeffs=[6.97981,0.0843473,-4.04404e-05,7.87407e-09,-5.56051e-13,-30926.4,13.239], Tmin=(469.866,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-249.697,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(644.372,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CCJCO)"""),
)

species(
    label = 'CCC[CH][CH]OC(O)CCC=O(20800)',
    structure = SMILES('CCC[CH][CH]OC(O)CCC=O'),
    E0 = (-294.946,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.18668,0.145166,-0.000157829,9.95711e-08,-2.6479e-11,-35258.9,51.5793], Tmin=(100,'K'), Tmax=(896.186,'K')), NASAPolynomial(coeffs=[14.6104,0.070194,-3.23425e-05,6.22135e-09,-4.37899e-13,-38269.5,-27.6157], Tmin=(896.186,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-294.946,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(640.214,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCsJOCs) + radical(CCJCO)"""),
)

species(
    label = '[CH2]CCCCOC([O])CCC=O(20801)',
    structure = SMILES('[CH2]CCCCOC([O])CCC=O'),
    E0 = (-244.353,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2759.09,2768.18,2777.27,2786.36,2795.45,2804.55,2813.64,2822.73,2831.82,2840.91,2850,1425,1430,1435,1440,1445,1450,1225,1235,1245,1255,1265,1275,1270,1284,1298,1312,1326,1340,700,720,740,760,780,800,300,320,340,360,380,400,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,180,180,180,368.694,798.546,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155284,'amu*angstrom^2'), symmetry=1, barrier=(3.57028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155284,'amu*angstrom^2'), symmetry=1, barrier=(3.57028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155284,'amu*angstrom^2'), symmetry=1, barrier=(3.57028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155284,'amu*angstrom^2'), symmetry=1, barrier=(3.57028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155284,'amu*angstrom^2'), symmetry=1, barrier=(3.57028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155284,'amu*angstrom^2'), symmetry=1, barrier=(3.57028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155284,'amu*angstrom^2'), symmetry=1, barrier=(3.57028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155284,'amu*angstrom^2'), symmetry=1, barrier=(3.57028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155284,'amu*angstrom^2'), symmetry=1, barrier=(3.57028,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.07673,0.146513,-0.000190785,1.6445e-07,-6.02483e-11,-29182.4,51.8591], Tmin=(100,'K'), Tmax=(763.169,'K')), NASAPolynomial(coeffs=[6.73701,0.0853349,-4.10908e-05,7.96005e-09,-5.58035e-13,-30091.4,14.5787], Tmin=(763.169,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-244.353,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(644.372,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCOJ) + radical(RCCJ)"""),
)

species(
    label = 'CCCC[CH]OC(O)CC[C]=O(20802)',
    structure = SMILES('CCCC[CH]OC(O)CC[C]=O'),
    E0 = (-334.888,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.61462,0.155309,-0.000188715,1.33253e-07,-3.90789e-11,-40048.4,51.4324], Tmin=(100,'K'), Tmax=(822.567,'K')), NASAPolynomial(coeffs=[15.082,0.0692532,-3.17869e-05,6.06682e-09,-4.23684e-13,-42959.7,-30.4871], Tmin=(822.567,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-334.888,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(640.214,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(CCsJOCs)"""),
)

species(
    label = 'CCCCCOC([O])C[CH]C=O(20803)',
    structure = SMILES('CCCCCOC([O])CC=C[O]'),
    E0 = (-306.348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.76615,0.130276,-0.000109542,4.99163e-08,-9.60428e-12,-36640.6,48.6723], Tmin=(100,'K'), Tmax=(1200.69,'K')), NASAPolynomial(coeffs=[17.3636,0.0665477,-2.99283e-05,5.71275e-09,-4.00637e-13,-41234.5,-47.1167], Tmin=(1200.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-306.348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(648.529,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(C=COJ)"""),
)

species(
    label = 'CCCCCOC([O])CC[C]=O(20804)',
    structure = SMILES('CCCCCOC([O])CC[C]=O'),
    E0 = (-289.639,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.19751,0.148844,-0.000195529,1.66197e-07,-5.96168e-11,-34624.4,51.326], Tmin=(100,'K'), Tmax=(770.074,'K')), NASAPolynomial(coeffs=[7.90273,0.0825132,-3.93129e-05,7.57409e-09,-5.29097e-13,-35768.8,7.90675], Tmin=(770.074,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-289.639,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(644.372,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(CCOJ)"""),
)

species(
    label = 'CC[CH]C[CH]OC(O)CCC=O(20805)',
    structure = SMILES('CC[CH]C[CH]OC(O)CCC=O'),
    E0 = (-300.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.5459,0.156234,-0.000211527,1.76505e-07,-6.1041e-11,-35904.6,52.9073], Tmin=(100,'K'), Tmax=(787.885,'K')), NASAPolynomial(coeffs=[10.392,0.0773508,-3.62169e-05,6.90383e-09,-4.78627e-13,-37533.7,-3.82631], Tmin=(787.885,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-300.39,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(640.214,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJCC) + radical(CCsJOCs)"""),
)

species(
    label = 'C[CH]CC[CH]OC(O)CCC=O(20806)',
    structure = SMILES('C[CH]CC[CH]OC(O)CCC=O'),
    E0 = (-300.402,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.55274,0.155358,-0.000204724,1.657e-07,-5.60053e-11,-35904.8,53.1065], Tmin=(100,'K'), Tmax=(773.373,'K')), NASAPolynomial(coeffs=[11.3952,0.0755614,-3.51047e-05,6.68345e-09,-4.63711e-13,-37833.2,-9.1196], Tmin=(773.373,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-300.402,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(640.214,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCsJOCs) + radical(RCCJC)"""),
)

species(
    label = '[CH2]CCC[CH]OC(O)CCC=O(20807)',
    structure = SMILES('[CH2]CCC[CH]OC(O)CCC=O'),
    E0 = (-289.602,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.45508,0.15245,-0.000181743,1.27938e-07,-3.78114e-11,-34608.1,51.8307], Tmin=(100,'K'), Tmax=(812.813,'K')), NASAPolynomial(coeffs=[13.8563,0.0721863,-3.36331e-05,6.46964e-09,-4.54066e-13,-37259.9,-23.4831], Tmin=(812.813,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-289.602,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(640.214,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJ) + radical(CCsJOCs)"""),
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
    label = 'CCCC[CH]OC1CC[CH]OO1(20808)',
    structure = SMILES('CCCC[CH]OC1CC[CH]OO1'),
    E0 = (-77.9325,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.40341,0.124898,-7.77565e-05,7.85287e-09,7.793e-12,-9128.38,44.3909], Tmin=(100,'K'), Tmax=(975.848,'K')), NASAPolynomial(coeffs=[23.758,0.0519099,-1.82082e-05,3.13554e-09,-2.12023e-13,-15865,-89.5384], Tmin=(975.848,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-77.9325,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(656.843,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(12dioxane) + radical(CCsJOOC) + radical(CCsJOCs)"""),
)

species(
    label = 'CCCCC1O[CH]CCC([O])O1(20734)',
    structure = SMILES('CCCCC1O[CH]CCC([O])O1'),
    E0 = (-275.651,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.01619,0.0863501,8.09738e-05,-1.95142e-07,8.9197e-11,-32894.1,46.9748], Tmin=(100,'K'), Tmax=(929.352,'K')), NASAPolynomial(coeffs=[41.0233,0.019022,-6.7852e-07,-4.28276e-11,-1.19157e-14,-45986.1,-184.909], Tmin=(929.352,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-275.651,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(665.158,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Cycloheptane) + radical(CCsJOCs) + radical(CCOJ)"""),
)

species(
    label = 'CCCCCOC(=O)CCC=O(20809)',
    structure = SMILES('CCCCCOC(=O)CCC=O'),
    E0 = (-666.825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.2915,0.126967,-0.000123182,8.03362e-08,-2.426e-11,-80019.4,45.2092], Tmin=(100,'K'), Tmax=(765.362,'K')), NASAPolynomial(coeffs=[6.84404,0.0844474,-3.98476e-05,7.74629e-09,-5.48412e-13,-81264.7,8.13565], Tmin=(765.362,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-666.825,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(644.372,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + group(Cds-OdCsH)"""),
)

species(
    label = 'CCCC=COC(O)CCC=O(20810)',
    structure = SMILES('CCCC=COC(O)CCC=O'),
    E0 = (-598.549,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.84629,0.138556,-0.000121297,5.46629e-08,-9.87578e-12,-71732.5,49.9064], Tmin=(100,'K'), Tmax=(1325.7,'K')), NASAPolynomial(coeffs=[26.8965,0.0488144,-1.97569e-05,3.60067e-09,-2.46529e-13,-79618.5,-101.971], Tmin=(1325.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-598.549,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(644.372,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH)"""),
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
    label = 'CCC[CH]OC([O])CCC=O(20811)',
    structure = SMILES('CCC[CH]OC([O])CCC=O'),
    E0 = (-245.363,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,2782.5,750,1395,475,1775,1000,3025,407.5,1350,352.5,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,180,1600,1628.56,2891.51,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152788,'amu*angstrom^2'), symmetry=1, barrier=(3.51289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152788,'amu*angstrom^2'), symmetry=1, barrier=(3.51289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152788,'amu*angstrom^2'), symmetry=1, barrier=(3.51289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152788,'amu*angstrom^2'), symmetry=1, barrier=(3.51289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152788,'amu*angstrom^2'), symmetry=1, barrier=(3.51289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152788,'amu*angstrom^2'), symmetry=1, barrier=(3.51289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152788,'amu*angstrom^2'), symmetry=1, barrier=(3.51289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152788,'amu*angstrom^2'), symmetry=1, barrier=(3.51289,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.368227,0.115383,-6.44007e-05,-1.02426e-07,1.36571e-10,-29373,41.4637], Tmin=(100,'K'), Tmax=(474.9,'K')), NASAPolynomial(coeffs=[8.86209,0.0717513,-3.43353e-05,6.62585e-09,-4.63289e-13,-30634.3,-0.243993], Tmin=(474.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-245.363,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCsJOCs) + radical(CCOJ)"""),
)

species(
    label = 'CO(12)',
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(4223.43,'J/mol'), sigma=(7.62069,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=659.69 K, Pc=21.65 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.07934,0.120801,-0.000119126,7.26186e-08,-1.95438e-11,-19376.8,41.6946], Tmin=(100,'K'), Tmax=(867.569,'K')), NASAPolynomial(coeffs=[9.73791,0.0709261,-3.28929e-05,6.35324e-09,-4.48422e-13,-21253.8,-8.95555], Tmin=(867.569,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-162.564,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(598.642,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCsJOCs) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C(CCC)OC([O])CCC=O(20812)',
    structure = SMILES('[CH2]C(CCC)OC([O])CCC=O'),
    E0 = (-249.64,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,2782.5,750,1395,475,1775,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,180,1600,1623.03,2899.18,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153056,'amu*angstrom^2'), symmetry=1, barrier=(3.51905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153056,'amu*angstrom^2'), symmetry=1, barrier=(3.51905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153056,'amu*angstrom^2'), symmetry=1, barrier=(3.51905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153056,'amu*angstrom^2'), symmetry=1, barrier=(3.51905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153056,'amu*angstrom^2'), symmetry=1, barrier=(3.51905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153056,'amu*angstrom^2'), symmetry=1, barrier=(3.51905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153056,'amu*angstrom^2'), symmetry=1, barrier=(3.51905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153056,'amu*angstrom^2'), symmetry=1, barrier=(3.51905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153056,'amu*angstrom^2'), symmetry=1, barrier=(3.51905,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.605291,0.125101,-3.7344e-05,-2.14976e-07,2.46431e-10,-29884,44.6835], Tmin=(100,'K'), Tmax=(447.974,'K')), NASAPolynomial(coeffs=[9.23149,0.0814268,-3.89682e-05,7.48859e-09,-5.20866e-13,-31208.5,0.180336], Tmin=(447.974,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-249.64,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(644.372,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CJC(C)OC)"""),
)

species(
    label = 'CCCCC1OC(CCC=O)O1(20311)',
    structure = SMILES('CCCCC1OC(CCC=O)O1'),
    E0 = (-526.772,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.6601,0.115753,-7.52932e-05,2.35869e-08,-2.97508e-12,-63146,44.3633], Tmin=(100,'K'), Tmax=(1799.6,'K')), NASAPolynomial(coeffs=[26.578,0.0529882,-2.29774e-05,4.20642e-09,-2.82756e-13,-73309.5,-108.461], Tmin=(1799.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-526.772,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(652.686,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Cyclobutane)"""),
)

species(
    label = 'C=COCC([O])O[CH]CCCC(20813)',
    structure = SMILES('C=COCC([O])O[CH]CCCC'),
    E0 = (-227.547,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
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
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.80924,0.148967,-0.000146408,7.61969e-08,-1.61004e-11,-27121.7,49.2013], Tmin=(100,'K'), Tmax=(1133.53,'K')), NASAPolynomial(coeffs=[23.232,0.0570728,-2.48046e-05,4.67784e-09,-3.26857e-13,-33025.4,-79.6964], Tmin=(1133.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-227.547,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(644.372,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CCsJOCs) + radical(CCOJ)"""),
)

species(
    label = 'CCCC[CH]OC([O])CC=CO(20814)',
    structure = SMILES('CCCC[CH]OC([O])CC=CO'),
    E0 = (-267.354,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1380,1390,370,380,2900,435,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,180,180,180,448.803,688.613,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156071,'amu*angstrom^2'), symmetry=1, barrier=(3.58837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156071,'amu*angstrom^2'), symmetry=1, barrier=(3.58837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156071,'amu*angstrom^2'), symmetry=1, barrier=(3.58837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156071,'amu*angstrom^2'), symmetry=1, barrier=(3.58837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156071,'amu*angstrom^2'), symmetry=1, barrier=(3.58837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156071,'amu*angstrom^2'), symmetry=1, barrier=(3.58837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156071,'amu*angstrom^2'), symmetry=1, barrier=(3.58837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156071,'amu*angstrom^2'), symmetry=1, barrier=(3.58837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156071,'amu*angstrom^2'), symmetry=1, barrier=(3.58837,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.84507,0.148569,-0.000147214,7.72479e-08,-1.63929e-11,-31907.1,50.5818], Tmin=(100,'K'), Tmax=(1132.34,'K')), NASAPolynomial(coeffs=[23.6046,0.0551347,-2.34417e-05,4.37666e-09,-3.04229e-13,-37897.1,-80.3098], Tmin=(1132.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-267.354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(644.372,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCsJOCs) + radical(CCOJ)"""),
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
    label = '[O]C([O])CCC=O(18653)',
    structure = SMILES('[O]C([O])CCC=O'),
    E0 = (-135.752,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,180,180,1940,1940.45],'cm^-1')),
        HinderedRotor(inertia=(0.0847157,'amu*angstrom^2'), symmetry=1, barrier=(1.94778,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0845682,'amu*angstrom^2'), symmetry=1, barrier=(1.94439,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0849048,'amu*angstrom^2'), symmetry=1, barrier=(1.95213,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.53632,0.0352151,-1.37843e-05,9.95017e-10,1.44007e-13,-16394.9,16.3348], Tmin=(100,'K'), Tmax=(2711.71,'K')), NASAPolynomial(coeffs=[49.2619,-0.011801,1.73631e-06,-2.42589e-10,2.04225e-14,-47621.7,-256.909], Tmin=(2711.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-135.752,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CCOJ)"""),
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
    label = 'CCCC[CH]O[CH]CCC=O(18726)',
    structure = SMILES('CCCC[CH]O[CH]CCC=O'),
    E0 = (-114.162,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,2782.5,750,1395,475,1775,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,180,180,180,180,1600,1650.41,2856.84,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153521,'amu*angstrom^2'), symmetry=1, barrier=(3.52975,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153521,'amu*angstrom^2'), symmetry=1, barrier=(3.52975,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153521,'amu*angstrom^2'), symmetry=1, barrier=(3.52975,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153521,'amu*angstrom^2'), symmetry=1, barrier=(3.52975,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153521,'amu*angstrom^2'), symmetry=1, barrier=(3.52975,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153521,'amu*angstrom^2'), symmetry=1, barrier=(3.52975,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153521,'amu*angstrom^2'), symmetry=1, barrier=(3.52975,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153521,'amu*angstrom^2'), symmetry=1, barrier=(3.52975,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153521,'amu*angstrom^2'), symmetry=1, barrier=(3.52975,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.2005,0.144816,-0.000167069,1.12268e-07,-3.14746e-11,-13514.6,45.4134], Tmin=(100,'K'), Tmax=(857.28,'K')), NASAPolynomial(coeffs=[14.5837,0.0665007,-3.00378e-05,5.70367e-09,-3.97769e-13,-16392.3,-32.976], Tmin=(857.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-114.162,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCsJOCs) + radical(CCsJOCs)"""),
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
    E0 = (-269.143,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-242.757,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-247.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-154.316,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-214.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-195.925,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-168.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-244.986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-140.913,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-111.107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-78.8341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-115.186,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-189.384,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-208.899,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-166.998,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-202.962,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-218.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-189.891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-236.406,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-242.261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-208.35,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-189.438,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-162.242,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-162.297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (97.1634,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-77.9325,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-184.313,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-205.743,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-244.17,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (174.499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (75.9493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-91.0135,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (-260.859,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (86.253,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (-123.488,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (144.379,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (128.843,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['CCCC[CH]OC([O])CCC=O(20308)'],
    products = ['CCCCC=O(1733)', 'O=CCCC=O(5767)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CCCC[CH]OC([O])CCC=O(20308)'],
    products = ['CCCC[CH]OC1CCC([O])O1(20789)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(132522,'s^-1'), n=1.48406, Ea=(26.3859,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSS;multiplebond_intra;radadd_intra] for rate rule [R6_SSS_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CCCC[CH]OC([O])CCC=O(20308)'],
    products = ['CCCCC1OC([O])CCC1[O](20761)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2305.07,'s^-1'), n=1.56671, Ea=(21.6023,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_SSSS;multiplebond_intra;radadd_intra_csHNd] for rate rule [R7_SSSS_CO;carbonylbond_intra_H;radadd_intra_csHNd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'CCCC=COC([O])CCC=O(20790)'],
    products = ['CCCC[CH]OC([O])CCC=O(20308)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.72e+08,'cm^3/(mol*s)'), n=1.477, Ea=(6.73624,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2825 used for Cds-CsH_Cds-OsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-OsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'CCCC[CH]OC(=O)CCC=O(20791)'],
    products = ['CCCC[CH]OC([O])CCC=O(20308)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.83701,'m^3/(mol*s)'), n=1.71338, Ea=(46.5052,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;HJ] + [CO-NdNd_O;YJ] for rate rule [CO-NdNd_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CCCCC=O(1733)', '[O][CH]CCC=O(5764)'],
    products = ['CCCC[CH]OC([O])CCC=O(20308)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4e+09,'cm^3/(mol*s)'), n=1.39, Ea=(35.8862,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Od_CO-CsH;YJ] for rate rule [Od_CO-CsH;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['npropyl(83)', 'C=COC([O])CCC=O(18636)'],
    products = ['CCCC[CH]OC([O])CCC=O(20308)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1071,'cm^3/(mol*s)'), n=2.72, Ea=(34.4762,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2780 used for Cds-HH_Cds-OsH;CsJ-CsHH
Exact match found for rate rule [Cds-HH_Cds-OsH;CsJ-CsHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH2CH2CHO(560)', 'CCCC[CH]OC=O(11683)'],
    products = ['CCCC[CH]OC([O])CCC=O(20308)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.94e+10,'cm^3/(mol*s)'), n=0, Ea=(28.0328,'kJ/mol'), T0=(1,'K'), Tmin=(333,'K'), Tmax=(363,'K'), comment="""Estimated using template [CO_O;CsJ-CsHH] for rate rule [CO-NdH_O;CsJ-CsHH]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CCC[CH]COC([O])CCC=O(20792)'],
    products = ['CCCC[CH]OC([O])CCC=O(20308)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.4e-20,'s^-1'), n=9.13, Ea=(108.784,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 341 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CCCC[CH]OC([O])CCC=O(20308)'],
    products = ['CCCC[CH]O[C](O)CCC=O(20793)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.43381e+07,'s^-1'), n=1.70481, Ea=(158.036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_NonDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CCCCCO[C]([O])CCC=O(20794)'],
    products = ['CCCC[CH]OC([O])CCC=O(20308)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.53164e+09,'s^-1'), n=1.09, Ea=(165.519,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_O;Y_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CC[CH]CCOC([O])CCC=O(20795)'],
    products = ['CCCC[CH]OC([O])CCC=O(20308)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.73898e+09,'s^-1'), n=0.87, Ea=(139.955,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_1H] for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CCCC[CH]OC(O)[CH]CC=O(20796)'],
    products = ['CCCC[CH]OC([O])CCC=O(20308)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(5.71,'s^-1'), n=3.021, Ea=(105.562,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 319 used for R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C[CH]CCCOC([O])CCC=O(20797)'],
    products = ['CCCC[CH]OC([O])CCC=O(20308)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.00310895,'s^-1'), n=4.065, Ea=(46.2541,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_1H] for rate rule [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CCCC[CH]OC([O])CCC=O(20308)'],
    products = ['CCCC[CH]OC(O)C[CH]C=O(20798)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(363473,'s^-1'), n=1.92229, Ea=(102.145,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R4H_SSS;O_rad_out;Cs_H_out_H/CO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CCCCCOC([O])[CH]CC=O(20799)'],
    products = ['CCCC[CH]OC([O])CCC=O(20308)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.00228,'s^-1'), n=3.95, Ea=(46.7353,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC] for rate rule [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CCCC[CH]OC([O])CCC=O(20308)'],
    products = ['CCC[CH][CH]OC(O)CCC=O(20800)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.19599e+09,'s^-1'), n=0.63, Ea=(50.8774,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5HJ_3;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]CCCCOC([O])CCC=O(20801)'],
    products = ['CCCC[CH]OC([O])CCC=O(20308)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(30977.5,'s^-1'), n=1.81667, Ea=(54.4617,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_2H;Cs_H_out_1H] for rate rule [R5H_CCC;C_rad_out_2H;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CCCC[CH]OC([O])CCC=O(20308)'],
    products = ['CCCC[CH]OC(O)CC[C]=O(20802)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(214655,'s^-1'), n=1.70206, Ea=(32.737,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;O_rad_out;XH_out] for rate rule [R5H_CCC_O;O_rad_out;CO_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['CCCC[CH]OC([O])CCC=O(20308)'],
    products = ['CCCCCOC([O])C[CH]C=O(20803)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.0756983,'s^-1'), n=3.26, Ea=(26.8822,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;C_rad_out_H/NonDeC;Cs_H_out_1H] for rate rule [R5H_SSSS_OCC;C_rad_out_H/NonDeC;Cs_H_out_H/CO]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CCCC[CH]OC([O])CCC=O(20308)'],
    products = ['CCCCCOC([O])CC[C]=O(20804)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(46.1,'s^-1'), n=3.21, Ea=(60.7935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_SSSSS;C_rad_out_H/NonDeC;XH_out] for rate rule [R6H_SSSSS_OO;C_rad_out_H/NonDeC;CO_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CCCC[CH]OC([O])CCC=O(20308)'],
    products = ['CC[CH]C[CH]OC(O)CCC=O(20805)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.12e+09,'s^-1'), n=0, Ea=(79.7052,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6Hall;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R6HJ_3;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CCCC[CH]OC([O])CCC=O(20308)'],
    products = ['C[CH]CC[CH]OC(O)CCC=O(20806)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.9e+08,'s^-1'), n=0, Ea=(106.901,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R7Hall;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R7HJ_3;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CCCC[CH]OC([O])CCC=O(20308)'],
    products = ['[CH2]CCC[CH]OC(O)CCC=O(20807)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.18263e+09,'s^-1'), n=0.595, Ea=(106.847,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;O_rad_out;Cs_H_out_2H] for rate rule [R8Hall;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CCCC[CH][O](1728)', '[O][CH]CCC=O(5764)'],
    products = ['CCCC[CH]OC([O])CCC=O(20308)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['CCCC[CH]OC([O])CCC=O(20308)'],
    products = ['CCCC[CH]OC1CC[CH]OO1(20808)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(463580,'s^-1'), n=1.14062, Ea=(191.211,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSS;multiplebond_intra;radadd_intra] for rate rule [R6_SSS_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic
Ea raised from 187.7 to 191.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction27',
    reactants = ['CCCC[CH]OC([O])CCC=O(20308)'],
    products = ['CCCCC1O[CH]CCC([O])O1(20734)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.07149e+11,'s^-1'), n=0.2025, Ea=(84.8306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;multiplebond_intra;radadd_intra_csHCs] for rate rule [R7_linear;carbonyl_intra_H;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CCCC[CH]OC([O])CCC=O(20308)'],
    products = ['CCCCCOC(=O)CCC=O(20809)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['CCCC[CH]OC([O])CCC=O(20308)'],
    products = ['CCCC=COC(O)CCC=O(20810)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['CH2(S)(23)', 'CCC[CH]OC([O])CCC=O(20811)'],
    products = ['CCCC[CH]OC([O])CCC=O(20308)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction31',
    reactants = ['CO(12)', 'CCCC[CH]OC([O])CC(13628)'],
    products = ['CCCC[CH]OC([O])CCC=O(20308)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(274200,'cm^3/(mol*s)'), n=2.53, Ea=(357.732,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 4 used for CO;C_pri/NonDeC
Exact match found for rate rule [CO;C_pri/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C(CCC)OC([O])CCC=O(20812)'],
    products = ['CCCC[CH]OC([O])CCC=O(20308)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction33',
    reactants = ['CCCC[CH]OC([O])CCC=O(20308)'],
    products = ['CCCCC1OC(CCC=O)O1(20311)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=COCC([O])O[CH]CCCC(20813)'],
    products = ['CCCC[CH]OC([O])CCC=O(20308)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 7 used for R_ROR;R1_doublebond_CH2;R2_doublebond_H;R_O_C
Exact match found for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction35',
    reactants = ['CCCC[CH]OC([O])CC=CO(20814)'],
    products = ['CCCC[CH]OC([O])CCC=O(20308)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_H;R_O_H]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]CCCC(141)', '[O]C([O])CCC=O(18653)'],
    products = ['CCCC[CH]OC([O])CCC=O(20308)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2711.41,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction37',
    reactants = ['O(4)', 'CCCC[CH]O[CH]CCC=O(18726)'],
    products = ['CCCC[CH]OC([O])CCC=O(20308)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/CsO;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

network(
    label = '3637',
    isomers = [
        'CCCC[CH]OC([O])CCC=O(20308)',
    ],
    reactants = [
        ('CCCCC=O(1733)', 'O=CCCC=O(5767)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '3637',
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

