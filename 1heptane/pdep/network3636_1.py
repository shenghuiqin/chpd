species(
    label = 'CCCCC([O])O[CH]CCC=O(20307)',
    structure = SMILES('CCCCC([O])O[CH]CCC=O'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.13372,0.132593,-9.493e-05,-4.91092e-08,9.34195e-11,-32205.4,46.4069], Tmin=(100,'K'), Tmax=(500.891,'K')), NASAPolynomial(coeffs=[9.49897,0.0807158,-3.84987e-05,7.44618e-09,-5.22582e-13,-33684.9,-1.67507], Tmin=(500.891,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-269.143,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(644.372,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCsJOCs) + radical(CCOJ)"""),
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
    label = 'CCCCC([O])OC1CCC1[O](20733)',
    structure = SMILES('CCCCC([O])OC1CCC1[O]'),
    E0 = (-163.376,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.73137,0.117627,-7.86846e-05,2.59297e-08,-3.47334e-12,-19437.1,47.5115], Tmin=(100,'K'), Tmax=(1696.51,'K')), NASAPolynomial(coeffs=[24.6023,0.0555376,-2.37866e-05,4.35658e-09,-2.94281e-13,-28372.1,-93.4522], Tmin=(1696.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-163.376,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(656.843,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CC(C)OJ) + radical(CCOJ)"""),
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
    label = 'CCCCC([O])OC=CCC=O(20735)',
    structure = SMILES('CCCCC([O])OC=CCC=O'),
    E0 = (-366.662,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,2782.5,750,1395,475,1775,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.37394,0.125618,-9.48095e-05,3.47979e-08,-5.09584e-12,-43858.1,48.6523], Tmin=(100,'K'), Tmax=(1603.12,'K')), NASAPolynomial(coeffs=[29.5649,0.0459265,-2.02441e-05,3.78942e-09,-2.60203e-13,-54098.4,-120.508], Tmin=(1603.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-366.662,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(623.585,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(CCOJ)"""),
)

species(
    label = 'CCCCC(=O)O[CH]CCC=O(20736)',
    structure = SMILES('CCCCC(=O)O[CH]CCC=O'),
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
    label = 'CH2CHO(40)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66874,0.0096233,1.60617e-05,-2.87682e-08,1.2503e-11,219.438,12.5694], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.91637,0.0088465,-3.14955e-06,5.05413e-10,-3.01305e-14,-1047.8,-6.1065], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(1.22925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CH2CHO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C=COC([O])CCCC(11991)',
    structure = SMILES('C=COC([O])CCCC'),
    E0 = (-231.374,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (129.177,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.3056,0.10093,-8.00819e-05,3.21756e-08,-5.1527e-12,-27623.7,38.5053], Tmin=(100,'K'), Tmax=(1492.64,'K')), NASAPolynomial(coeffs=[23.4921,0.034477,-1.33006e-05,2.3485e-09,-1.56992e-13,-35026.4,-91.0619], Tmin=(1492.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-231.374,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CCOJ)"""),
)

species(
    label = 'butyl_1(82)',
    structure = SMILES('[CH2]CCC'),
    E0 = (63.0573,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.0977402,'amu*angstrom^2'), symmetry=1, barrier=(2.24724,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0976865,'amu*angstrom^2'), symmetry=1, barrier=(2.246,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0977534,'amu*angstrom^2'), symmetry=1, barrier=(2.24754,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.1143,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.176,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.25388,0.0316763,2.89994e-06,-1.98049e-08,8.20503e-12,7652.64,17.2725], Tmin=(100,'K'), Tmax=(1050.57,'K')), NASAPolynomial(coeffs=[7.59591,0.0260842,-1.01719e-05,1.85189e-09,-1.28169e-13,5716.37,-12.6366], Tmin=(1050.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(63.0573,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), label="""butyl_1""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'O=CCC[CH]OC=O(18676)',
    structure = SMILES('O=CCC[CH]OC=O'),
    E0 = (-343.299,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,3025,407.5,1350,352.5,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.984844,0.071077,-6.77199e-05,3.62702e-08,-8.30577e-12,-41184.8,28.4337], Tmin=(100,'K'), Tmax=(1018.39,'K')), NASAPolynomial(coeffs=[9.65664,0.0370164,-1.75519e-05,3.42896e-09,-2.43762e-13,-42951,-13.5607], Tmin=(1018.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-343.299,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + group(Cds-OdOsH) + radical(CCsJOC(O)H)"""),
)

species(
    label = 'CCCCC([O])OC[CH]CC=O(20737)',
    structure = SMILES('CCCCC([O])OC[CH]CC=O'),
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
    label = 'CCCC[C](O)O[CH]CCC=O(20738)',
    structure = SMILES('CCCC[C](O)O[CH]CCC=O'),
    E0 = (-289.602,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.45508,0.15245,-0.000181743,1.27938e-07,-3.78114e-11,-34608.1,50.7321], Tmin=(100,'K'), Tmax=(812.813,'K')), NASAPolynomial(coeffs=[13.8563,0.0721863,-3.36331e-05,6.46964e-09,-4.54066e-13,-37259.9,-24.5817], Tmin=(812.813,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-289.602,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(640.214,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Cs_P) + radical(CCsJOCs)"""),
)

species(
    label = 'CCCC[C]([O])OCCCC=O(20739)',
    structure = SMILES('CCCC[C]([O])OCCCC=O'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.07672,0.146513,-0.000190784,1.64449e-07,-6.02479e-11,-29182.4,50.7605], Tmin=(100,'K'), Tmax=(763.171,'K')), NASAPolynomial(coeffs=[6.737,0.0853349,-4.10908e-05,7.96005e-09,-5.58035e-13,-30091.4,13.4802], Tmin=(763.171,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-244.353,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(644.372,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Cs_P) + radical(CCOJ)"""),
)

species(
    label = 'CCCCC([O])OCC[CH]C=O(20740)',
    structure = SMILES('CCCCC([O])OCCC=C[O]'),
    E0 = (-306.348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.76615,0.130276,-0.000109542,4.99163e-08,-9.60428e-12,-36640.6,48.6723], Tmin=(100,'K'), Tmax=(1200.69,'K')), NASAPolynomial(coeffs=[17.3636,0.0665477,-2.99283e-05,5.71275e-09,-4.00637e-13,-41234.5,-47.1167], Tmin=(1200.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-306.348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(648.529,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CCOJ)"""),
)

species(
    label = 'CCC[CH]C(O)O[CH]CCC=O(20741)',
    structure = SMILES('CCC[CH]C(O)O[CH]CCC=O'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.18668,0.145166,-0.000157829,9.95711e-08,-2.6479e-11,-35258.9,51.5793], Tmin=(100,'K'), Tmax=(896.186,'K')), NASAPolynomial(coeffs=[14.6104,0.070194,-3.23425e-05,6.22135e-09,-4.37899e-13,-38269.5,-27.6157], Tmin=(896.186,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-294.946,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(640.214,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCsJOCs) + radical(CCJCO)"""),
)

species(
    label = 'CC[CH]CC(O)O[CH]CCC=O(20742)',
    structure = SMILES('CC[CH]CC(O)O[CH]CCC=O'),
    E0 = (-300.39,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,2782.5,750,1395,475,1775,1000,180,180,180,180,1600,1647.63,2873.32,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15387,'amu*angstrom^2'), symmetry=1, barrier=(3.53777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15387,'amu*angstrom^2'), symmetry=1, barrier=(3.53777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15387,'amu*angstrom^2'), symmetry=1, barrier=(3.53777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15387,'amu*angstrom^2'), symmetry=1, barrier=(3.53777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15387,'amu*angstrom^2'), symmetry=1, barrier=(3.53777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15387,'amu*angstrom^2'), symmetry=1, barrier=(3.53777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15387,'amu*angstrom^2'), symmetry=1, barrier=(3.53777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15387,'amu*angstrom^2'), symmetry=1, barrier=(3.53777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15387,'amu*angstrom^2'), symmetry=1, barrier=(3.53777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15387,'amu*angstrom^2'), symmetry=1, barrier=(3.53777,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.5459,0.156234,-0.000211527,1.76505e-07,-6.1041e-11,-35904.6,52.9073], Tmin=(100,'K'), Tmax=(787.885,'K')), NASAPolynomial(coeffs=[10.392,0.0773508,-3.62169e-05,6.90383e-09,-4.78627e-13,-37533.7,-3.82631], Tmin=(787.885,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-300.39,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(640.214,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJCC) + radical(CCsJOCs)"""),
)

species(
    label = 'CCCCC([O])OCCC[C]=O(20743)',
    structure = SMILES('CCCCC([O])OCCC[C]=O'),
    E0 = (-289.639,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.19751,0.148844,-0.000195529,1.66197e-07,-5.96168e-11,-34624.4,51.326], Tmin=(100,'K'), Tmax=(770.074,'K')), NASAPolynomial(coeffs=[7.90273,0.0825132,-3.93129e-05,7.57409e-09,-5.29097e-13,-35768.8,7.90675], Tmin=(770.074,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-289.639,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(644.372,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CCCJ=O)"""),
)

species(
    label = 'CCC[CH]C([O])OCCCC=O(20744)',
    structure = SMILES('CCC[CH]C([O])OCCCC=O'),
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
    label = 'CCCCC(O)O[CH][CH]CC=O(20745)',
    structure = SMILES('CCCCC(O)O[CH][CH]CC=O'),
    E0 = (-294.946,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.18668,0.145166,-0.000157829,9.95711e-08,-2.6479e-11,-35258.9,51.5793], Tmin=(100,'K'), Tmax=(896.186,'K')), NASAPolynomial(coeffs=[14.6104,0.070194,-3.23425e-05,6.22135e-09,-4.37899e-13,-38269.5,-27.6157], Tmin=(896.186,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-294.946,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(640.214,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCsJOCs) + radical(CCJCO)"""),
)

species(
    label = 'C[CH]CCC(O)O[CH]CCC=O(20746)',
    structure = SMILES('C[CH]CCC(O)O[CH]CCC=O'),
    E0 = (-300.402,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.55274,0.155358,-0.000204724,1.657e-07,-5.60053e-11,-35904.8,53.1065], Tmin=(100,'K'), Tmax=(773.373,'K')), NASAPolynomial(coeffs=[11.3952,0.0755614,-3.51047e-05,6.68345e-09,-4.63711e-13,-37833.2,-9.1196], Tmin=(773.373,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-300.402,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(640.214,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCsJOCs) + radical(RCCJC)"""),
)

species(
    label = 'CC[CH]CC([O])OCCCC=O(20747)',
    structure = SMILES('CC[CH]CC([O])OCCCC=O'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.93812,0.147244,-0.000208009,1.93414e-07,-7.32813e-11,-30488.7,52.1336], Tmin=(100,'K'), Tmax=(804.959,'K')), NASAPolynomial(coeffs=[2.96165,0.0910685,-4.40199e-05,8.47885e-09,-5.89806e-13,-30246.4,35.9631], Tmin=(804.959,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-255.141,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(644.372,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCOJ) + radical(RCCJCC)"""),
)

species(
    label = '[CH2]CCCC(O)O[CH]CCC=O(20748)',
    structure = SMILES('[CH2]CCCC(O)O[CH]CCC=O'),
    E0 = (-289.602,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.45508,0.15245,-0.000181743,1.27938e-07,-3.78114e-11,-34608.1,51.8307], Tmin=(100,'K'), Tmax=(812.813,'K')), NASAPolynomial(coeffs=[13.8563,0.0721863,-3.36331e-05,6.46964e-09,-4.54066e-13,-37259.9,-23.4831], Tmin=(812.813,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-289.602,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(640.214,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJ) + radical(CCsJOCs)"""),
)

species(
    label = 'C[CH]CCC([O])OCCCC=O(20749)',
    structure = SMILES('C[CH]CCC([O])OCCCC=O'),
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
    label = 'CCCCC(O)O[CH]C[CH]C=O(20750)',
    structure = SMILES('CCCCC(O)O[CH]CC=C[O]'),
    E0 = (-351.597,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.98993,0.146676,-0.000139828,6.93488e-08,-1.37945e-11,-42029.8,51.6476], Tmin=(100,'K'), Tmax=(1210.25,'K')), NASAPolynomial(coeffs=[26.0316,0.0507565,-2.09432e-05,3.86096e-09,-2.66661e-13,-49054.4,-93.9022], Tmin=(1210.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-351.597,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(644.372,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CCsJOCs)"""),
)

species(
    label = '[CH2]CCCC([O])OCCCC=O(20751)',
    structure = SMILES('[CH2]CCCC([O])OCCCC=O'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.07672,0.146513,-0.000190784,1.64449e-07,-6.02479e-11,-29182.4,51.8591], Tmin=(100,'K'), Tmax=(763.171,'K')), NASAPolynomial(coeffs=[6.737,0.0853349,-4.10908e-05,7.96005e-09,-5.58035e-13,-30091.4,14.5788], Tmin=(763.171,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-244.353,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(644.372,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJ) + radical(CCOJ)"""),
)

species(
    label = 'CCCCC(O)O[CH]CC[C]=O(20752)',
    structure = SMILES('CCCCC(O)O[CH]CC[C]=O'),
    E0 = (-334.888,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.61462,0.155309,-0.000188715,1.33253e-07,-3.90789e-11,-40048.4,51.4324], Tmin=(100,'K'), Tmax=(822.567,'K')), NASAPolynomial(coeffs=[15.082,0.0692532,-3.17869e-05,6.06682e-09,-4.23684e-13,-42959.7,-30.4871], Tmin=(822.567,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-334.888,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(640.214,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(CCsJOCs)"""),
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
    label = 'CCCCC([O])OC1CC[CH]O1(20753)',
    structure = SMILES('CCCCC([O])OC1CC[CH]O1'),
    E0 = (-281.006,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.05946,0.121685,-9.0678e-05,3.66617e-08,-6.13282e-12,-33569.6,45.7978], Tmin=(100,'K'), Tmax=(1397.43,'K')), NASAPolynomial(coeffs=[20.3722,0.0574759,-2.17553e-05,3.78053e-09,-2.50319e-13,-39838.8,-69.928], Tmin=(1397.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-281.006,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(656.843,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(CCsJOCs) + radical(CCOJ)"""),
)

species(
    label = 'CCCCC1O[CH]CC[CH]OO1(20526)',
    structure = SMILES('CCCCC1O[CH]CC[CH]OO1'),
    E0 = (-48.3707,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.98299,0.113333,-5.75515e-05,3.83782e-09,3.61603e-12,-5587.29,43.0725], Tmin=(100,'K'), Tmax=(1222.59,'K')), NASAPolynomial(coeffs=[22.6848,0.0600259,-2.57641e-05,4.83494e-09,-3.36141e-13,-13666.7,-89.2672], Tmin=(1222.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-48.3707,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(665.158,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Cyclooctane) + radical(CCsJOCs) + radical(CCsJOOC)"""),
)

species(
    label = 'CCCCC(=O)OCCCC=O(20754)',
    structure = SMILES('CCCCC(=O)OCCCC=O'),
    E0 = (-666.825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.2915,0.126967,-0.000123182,8.03362e-08,-2.426e-11,-80019.4,45.2092], Tmin=(100,'K'), Tmax=(765.362,'K')), NASAPolynomial(coeffs=[6.84404,0.0844474,-3.98476e-05,7.74629e-09,-5.48412e-13,-81264.7,8.13565], Tmin=(765.362,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-666.825,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(644.372,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + group(Cds-OdCsH)"""),
)

species(
    label = 'CCCCC(O)OC=CCC=O(20755)',
    structure = SMILES('CCCCC(O)OC=CCC=O'),
    E0 = (-592.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.85134,0.130001,-8.96861e-05,2.12967e-08,1.02378e-12,-70981.2,49.4884], Tmin=(100,'K'), Tmax=(1136.76,'K')), NASAPolynomial(coeffs=[29.2929,0.047074,-2.00848e-05,3.83361e-09,-2.72024e-13,-80239.3,-118.286], Tmin=(1136.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-592.367,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(644.372,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH)"""),
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
    label = 'CCCC([O])O[CH]CCC=O(20756)',
    structure = SMILES('CCCC([O])O[CH]CCC=O'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.368227,0.115383,-6.44007e-05,-1.02426e-07,1.36571e-10,-29373,41.4637], Tmin=(100,'K'), Tmax=(474.9,'K')), NASAPolynomial(coeffs=[8.86209,0.0717513,-3.43353e-05,6.62585e-09,-4.63289e-13,-30634.3,-0.243993], Tmin=(474.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-245.363,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CCsJOCs)"""),
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
    label = 'CC[CH]OC([O])CCCC(13627)',
    structure = SMILES('CC[CH]OC([O])CCCC'),
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
    label = '[CH2]C(CC=O)OC([O])CCCC(20757)',
    structure = SMILES('[CH2]C(CC=O)OC([O])CCCC'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.605291,0.125101,-3.7344e-05,-2.14976e-07,2.46431e-10,-29884,44.6835], Tmin=(100,'K'), Tmax=(447.974,'K')), NASAPolynomial(coeffs=[9.23149,0.0814268,-3.89682e-05,7.48859e-09,-5.20866e-13,-31208.5,0.180336], Tmin=(447.974,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-249.64,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(644.372,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJC(C)OC) + radical(CCOJ)"""),
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
    label = 'C=COC[CH]OC([O])CCCC(20758)',
    structure = SMILES('C=COC[CH]OC([O])CCCC'),
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
    label = 'CCCCC([O])O[CH]CC=CO(20759)',
    structure = SMILES('CCCCC([O])O[CH]CC=CO'),
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
    label = 'CCCCC([O])[O](4945)',
    structure = SMILES('CCCCC([O])[O]'),
    E0 = (-76.7325,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.132,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.91049,0.0514419,-2.35394e-05,4.03792e-09,-2.10455e-13,-9217.86,22.5095], Tmin=(100,'K'), Tmax=(2688.64,'K')), NASAPolynomial(coeffs=[41.7048,0.00470864,-3.59396e-06,6.11596e-10,-3.31305e-14,-34048.1,-210.401], Tmin=(2688.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-76.7325,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[CH]CCC=O(18484)',
    structure = SMILES('[CH]CCC=O'),
    E0 = (221.112,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,338.309,338.317,1520,1520.03],'cm^-1')),
        HinderedRotor(inertia=(0.00147232,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.09193,'amu*angstrom^2'), symmetry=1, barrier=(7.47087,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0919401,'amu*angstrom^2'), symmetry=1, barrier=(7.47082,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.39867,0.0383489,-2.66828e-05,1.05264e-08,-1.89984e-12,26648.5,18.3896], Tmin=(100,'K'), Tmax=(1182.31,'K')), NASAPolynomial(coeffs=[5.79639,0.0268536,-1.20987e-05,2.30283e-09,-1.60963e-13,25845.1,1.42861], Tmin=(1182.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(221.112,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJ2_triplet)"""),
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
    E0 = (-163.376,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-173.263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-148.133,'kJ/mol'),
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
    E0 = (-194.001,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-252.209,'kJ/mol'),
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
    E0 = (-141.531,'kJ/mol'),
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
    E0 = (-216.292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-172.493,'kJ/mol'),
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
    E0 = (-243.412,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-230.414,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-175.631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-212.476,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-185.734,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-189.124,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-192.755,'kJ/mol'),
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
    E0 = (-210.986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-48.3707,'kJ/mol'),
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
    reactants = ['CCCCC([O])O[CH]CCC=O(20307)'],
    products = ['CCCCC=O(1733)', 'O=CCCC=O(5767)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CCCCC([O])O[CH]CCC=O(20307)'],
    products = ['CCCCC([O])OC1CCC1[O](20733)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.13771e+06,'s^-1'), n=1.58803, Ea=(105.768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_csHNd] for rate rule [R5_SS_CO;carbonylbond_intra_H;radadd_intra_csHNd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 101.2 to 105.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['CCCCC([O])O[CH]CCC=O(20307)'],
    products = ['CCCCC1O[CH]CCC([O])O1(20734)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.3251e+11,'s^-1'), n=0.157, Ea=(95.8805,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7plus;multiplebond_intra;radadd_intra_O] + [R8;multiplebond_intra;radadd_intra] for rate rule [R8;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'CCCCC([O])OC=CCC=O(20735)'],
    products = ['CCCCC([O])O[CH]CCC=O(20307)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.72e+08,'cm^3/(mol*s)'), n=1.477, Ea=(6.73624,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2825 used for Cds-CsH_Cds-OsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-OsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'CCCCC(=O)O[CH]CCC=O(20736)'],
    products = ['CCCCC([O])O[CH]CCC=O(20307)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.83701,'m^3/(mol*s)'), n=1.71338, Ea=(46.5052,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;HJ] + [CO-NdNd_O;YJ] for rate rule [CO-NdNd_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CCCC[CH][O](1728)', 'O=CCCC=O(5767)'],
    products = ['CCCCC([O])O[CH]CCC=O(20307)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(8e+09,'cm^3/(mol*s)'), n=1.39, Ea=(35.8862,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Od_CO-CsH;YJ] for rate rule [Od_CO-CsH;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2CHO(40)', 'C=COC([O])CCCC(11991)'],
    products = ['CCCCC([O])O[CH]CCC=O(20307)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.0114756,'m^3/(mol*s)'), n=2.44484, Ea=(36.1431,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-HH_Cds;CsJ-OneDeHH] for rate rule [Cds-HH_Cds-OsH;CsJ-COHH]
Euclidian distance = 1.41421356237
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['butyl_1(82)', 'O=CCC[CH]OC=O(18676)'],
    products = ['CCCCC([O])O[CH]CCC=O(20307)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.94e+10,'cm^3/(mol*s)'), n=0, Ea=(28.0328,'kJ/mol'), T0=(1,'K'), Tmin=(333,'K'), Tmax=(363,'K'), comment="""Estimated using template [CO_O;CsJ-CsHH] for rate rule [CO-NdH_O;CsJ-CsHH]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CCCCC([O])OC[CH]CC=O(20737)'],
    products = ['CCCCC([O])O[CH]CCC=O(20307)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.4e-20,'s^-1'), n=9.13, Ea=(108.784,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 341 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CCCCC([O])O[CH]CCC=O(20307)'],
    products = ['CCCC[C](O)O[CH]CCC=O(20738)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.43381e+07,'s^-1'), n=1.70481, Ea=(158.036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_NonDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CCCC[C]([O])OCCCC=O(20739)'],
    products = ['CCCCC([O])O[CH]CCC=O(20307)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.53164e+09,'s^-1'), n=1.09, Ea=(165.519,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_O;Y_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CCCCC([O])O[CH]CCC=O(20307)'],
    products = ['CCCCC([O])OCC[CH]C=O(20740)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(6.82e+09,'s^-1'), n=0.73, Ea=(127.612,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_1H;Cs_H_out_H/OneDe] for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeO;Cs_H_out_H/CO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CCC[CH]C(O)O[CH]CCC=O(20741)'],
    products = ['CCCCC([O])O[CH]CCC=O(20307)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(5.71,'s^-1'), n=3.021, Ea=(105.562,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 319 used for R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['CC[CH]CC(O)O[CH]CCC=O(20742)'],
    products = ['CCCCC([O])O[CH]CCC=O(20307)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2960,'s^-1'), n=2.11, Ea=(84.0984,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 323 used for R4H_SSS;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CCCCC([O])O[CH]CCC=O(20307)'],
    products = ['CCCCC([O])OCCC[C]=O(20743)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(27900,'s^-1'), n=1.97, Ea=(96.6504,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/NonDeO;XH_out] for rate rule [R4H_SSS;C_rad_out_H/NonDeO;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CCC[CH]C([O])OCCCC=O(20744)'],
    products = ['CCCCC([O])O[CH]CCC=O(20307)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.00228,'s^-1'), n=3.95, Ea=(46.7353,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC] for rate rule [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CCCCC([O])O[CH]CCC=O(20307)'],
    products = ['CCCCC(O)O[CH][CH]CC=O(20745)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.19599e+09,'s^-1'), n=0.63, Ea=(50.8774,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5HJ_3;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CCCCC([O])O[CH]CCC=O(20307)'],
    products = ['C[CH]CCC(O)O[CH]CCC=O(20746)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(8e+10,'s^-1'), n=0, Ea=(25.7316,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(1000,'K'), comment="""From training reaction 307 used for R5H_CCC;O_rad_out;Cs_H_out_H/NonDeC
Exact match found for rate rule [R5H_CCC;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CC[CH]CC([O])OCCCC=O(20747)'],
    products = ['CCCCC([O])O[CH]CCC=O(20307)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.0564,'s^-1'), n=3.28, Ea=(24.7274,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5H_SSSS;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC] for rate rule [R5H_SSSS;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['CCCCC([O])O[CH]CCC=O(20307)'],
    products = ['[CH2]CCCC(O)O[CH]CCC=O(20748)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.68e+09,'s^-1'), n=0, Ea=(93.5124,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 140 used for R6H_SSSSS;O_rad_out;Cs_H_out_2H
Exact match found for rate rule [R6H_SSSSS;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C[CH]CCC([O])OCCCC=O(20749)'],
    products = ['CCCCC([O])O[CH]CCC=O(20307)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(25800,'s^-1'), n=1.67, Ea=(42.6768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_SSSSS;C_rad_out_single;Cs_H_out_H/NonDeC] for rate rule [R6H_SSSSS;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CCCCC([O])O[CH]CCC=O(20307)'],
    products = ['CCCCC(O)O[CH]C[CH]C=O(20750)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(49528.1,'s^-1'), n=1.95205, Ea=(83.4098,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R6HJ_3;O_rad_out;Cs_H_out_H/CO]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]CCCC([O])OCCCC=O(20751)'],
    products = ['CCCCC([O])O[CH]CCC=O(20307)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1062,'s^-1'), n=1.81, Ea=(55.2288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7H;C_rad_out_2H;Cs_H_out_H/NonDeC] for rate rule [R7H;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CCCCC([O])O[CH]CCC=O(20307)'],
    products = ['CCCCC(O)O[CH]CC[C]=O(20752)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(970995,'s^-1'), n=0.905106, Ea=(76.3888,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7Hall;O_rad_out;XH_out] for rate rule [R7HJ_3;O_rad_out;CO_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CCCC[CH][O](1728)', '[O][CH]CCC=O(5764)'],
    products = ['CCCCC([O])O[CH]CCC=O(20307)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['CCCCC([O])O[CH]CCC=O(20307)'],
    products = ['CCCCC([O])OC1CC[CH]O1(20753)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.64e+06,'s^-1'), n=1.15, Ea=(58.1576,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_csHNd] for rate rule [R5_SS_CO;carbonyl_intra_H;radadd_intra_csHO]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CCCCC([O])O[CH]CCC=O(20307)'],
    products = ['CCCCC1O[CH]CC[CH]OO1(20526)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(5.649e+12,'s^-1'), n=0.287, Ea=(220.773,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;carbonyl_intra_H;radadd_intra] for rate rule [R8_linear;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic
Ea raised from 215.4 to 220.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction28',
    reactants = ['CCCCC([O])O[CH]CCC=O(20307)'],
    products = ['CCCCC(=O)OCCCC=O(20754)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['CCCCC([O])O[CH]CCC=O(20307)'],
    products = ['CCCCC(O)OC=CCC=O(20755)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['CH2(S)(23)', 'CCCC([O])O[CH]CCC=O(20756)'],
    products = ['CCCCC([O])O[CH]CCC=O(20307)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction31',
    reactants = ['CO(12)', 'CC[CH]OC([O])CCCC(13627)'],
    products = ['CCCCC([O])O[CH]CCC=O(20307)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(274200,'cm^3/(mol*s)'), n=2.53, Ea=(357.732,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 4 used for CO;C_pri/NonDeC
Exact match found for rate rule [CO;C_pri/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C(CC=O)OC([O])CCCC(20757)'],
    products = ['CCCCC([O])O[CH]CCC=O(20307)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction33',
    reactants = ['CCCCC([O])O[CH]CCC=O(20307)'],
    products = ['CCCCC1OC(CCC=O)O1(20311)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=COC[CH]OC([O])CCCC(20758)'],
    products = ['CCCCC([O])O[CH]CCC=O(20307)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 7 used for R_ROR;R1_doublebond_CH2;R2_doublebond_H;R_O_C
Exact match found for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction35',
    reactants = ['CCCCC([O])O[CH]CC=CO(20759)'],
    products = ['CCCCC([O])O[CH]CCC=O(20307)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_H;R_O_H]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction36',
    reactants = ['CCCCC([O])[O](4945)', '[CH]CCC=O(18484)'],
    products = ['CCCCC([O])O[CH]CCC=O(20307)'],
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
    products = ['CCCCC([O])O[CH]CCC=O(20307)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/CsO;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

network(
    label = '3636',
    isomers = [
        'CCCCC([O])O[CH]CCC=O(20307)',
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
    label = '3636',
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

