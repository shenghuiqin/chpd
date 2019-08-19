species(
    label = 'O=[C]OO[CH]CCC=O(21319)',
    structure = SMILES('O=[C]OO[CH]CCC=O'),
    E0 = (-174.978,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,1855,455,950,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,200,800],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.376243,0.0820076,-7.81698e-05,3.8857e-08,-7.95388e-12,-20916.5,31.8304], Tmin=(100,'K'), Tmax=(1152.53,'K')), NASAPolynomial(coeffs=[14.0013,0.0347201,-1.66259e-05,3.25778e-09,-2.31915e-13,-24057.2,-35.8371], Tmin=(1152.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-174.978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + group(Cds-OdOsH) + radical((O)CJOC) + radical(CCsJOOC)"""),
)

species(
    label = 'CO2(13)',
    structure = SMILES('O=C=O'),
    E0 = (-403.131,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([459.166,1086.67,1086.68,2300.05],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0095,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1622.99,'J/mol'), sigma=(3.941,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.2779,0.00275783,7.12787e-06,-1.07855e-08,4.14228e-12,-48475.6,5.97856], Tmin=(100,'K'), Tmax=(988.185,'K')), NASAPolynomial(coeffs=[4.55071,0.00290728,-1.14643e-06,2.25798e-10,-1.69526e-14,-48986,-1.45662], Tmin=(988.185,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-403.131,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), label="""CO2""", comment="""Thermo library: BurkeH2O2"""),
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
    label = '[O]C1CCC1OO[C]=O(22130)',
    structure = SMILES('[O]C1CCC1OO[C]=O'),
    E0 = (-75.1952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.39685,0.0566461,1.82633e-05,-7.29796e-08,3.36118e-11,-8893.87,30.1192], Tmin=(100,'K'), Tmax=(975.774,'K')), NASAPolynomial(coeffs=[23.0278,0.0183964,-6.75025e-06,1.37227e-09,-1.08563e-13,-15906,-91.8069], Tmin=(975.774,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-75.1952,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-O2s(Cds-O2d)) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-OdOsH) + ring(Cyclobutane) + radical((O)CJOC) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]C1CC[CH]OOC1=O(22131)',
    structure = SMILES('[O]C1CC[CH]OOC1=O'),
    E0 = (-168.337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15391,0.000482196,0.000242873,-3.61543e-07,1.52124e-10,-20085.6,31.6857], Tmin=(100,'K'), Tmax=(909.992,'K')), NASAPolynomial(coeffs=[44.3764,-0.0275894,2.22431e-05,-4.37319e-09,2.80736e-13,-34656.1,-209.597], Tmin=(909.992,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-168.337,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + ring(Cycloheptane) + radical(C=OCOJ) + radical(CCsJOOC)"""),
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
    label = 'O=[C]OOC=CCC=O(22132)',
    structure = SMILES('O=[C]OOC=CCC=O'),
    E0 = (-223.43,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2782.5,750,1395,475,1775,1000,1855,455,950,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,1437.5,1250,1305,750,350,180],'cm^-1')),
        HinderedRotor(inertia=(1.4176,'amu*angstrom^2'), symmetry=1, barrier=(32.5934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41779,'amu*angstrom^2'), symmetry=1, barrier=(32.5978,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41698,'amu*angstrom^2'), symmetry=1, barrier=(32.5791,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41756,'amu*angstrom^2'), symmetry=1, barrier=(32.5925,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41786,'amu*angstrom^2'), symmetry=1, barrier=(32.5995,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (129.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.328857,0.0632627,-1.81674e-05,-2.73027e-08,1.52795e-11,-26725,32.1706], Tmin=(100,'K'), Tmax=(1040.61,'K')), NASAPolynomial(coeffs=[21.6299,0.0192406,-9.28068e-06,1.96395e-09,-1.50544e-13,-33208,-81.2913], Tmin=(1040.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-223.43,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-O2d)) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-OdOsH) + radical((O)CJOC)"""),
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
    label = 'C=COO[C]=O(745)',
    structure = SMILES('C=COO[C]=O'),
    E0 = (-88.1416,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,1855,455,950,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(1.70914,'amu*angstrom^2'), symmetry=1, barrier=(39.2966,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.71763,'amu*angstrom^2'), symmetry=1, barrier=(39.4916,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.7168,'amu*angstrom^2'), symmetry=1, barrier=(39.4727,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.0541,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75162,0.034355,1.16482e-05,-5.01685e-08,2.42041e-11,-10506.3,20.7528], Tmin=(100,'K'), Tmax=(949.643,'K')), NASAPolynomial(coeffs=[17.7534,0.00454878,-6.54467e-07,1.56005e-10,-1.87072e-14,-15240.7,-64.5451], Tmin=(949.643,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-88.1416,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-O2d)) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-OdOsH) + radical((O)CJOC)"""),
)

species(
    label = 'O=[C]OOC[CH]CC=O(22133)',
    structure = SMILES('O=[C]OOC[CH]CC=O'),
    E0 = (-161.004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.394631,0.0821605,-7.97804e-05,4.05426e-08,-8.49047e-12,-19236.9,33.7758], Tmin=(100,'K'), Tmax=(1128,'K')), NASAPolynomial(coeffs=[13.7593,0.0347674,-1.67569e-05,3.2943e-09,-2.34991e-13,-22251.9,-32.3108], Tmin=(1128,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-161.004,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + group(Cds-OdOsH) + radical(CCJCOOH) + radical((O)CJOC)"""),
)

species(
    label = 'O=[C]OOCC[CH]C=O(22134)',
    structure = SMILES('[O]C=CCCOO[C]=O'),
    E0 = (-218.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.278503,0.0766052,-3.67734e-05,-1.99216e-08,1.61218e-11,-26069.7,33.5905], Tmin=(100,'K'), Tmax=(976.145,'K')), NASAPolynomial(coeffs=[24.3102,0.0161105,-5.68518e-06,1.10199e-09,-8.45892e-14,-32788.4,-94.2676], Tmin=(976.145,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-218.167,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdOsH) + radical(C=COJ) + radical((O)CJOC)"""),
)

species(
    label = 'O=[C]CCCOO[C]=O(22135)',
    structure = SMILES('O=[C]CCCOO[C]=O'),
    E0 = (-201.458,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.147287,0.0848953,-8.53542e-05,4.4202e-08,-9.24862e-12,-24091.2,33.1745], Tmin=(100,'K'), Tmax=(1144.73,'K')), NASAPolynomial(coeffs=[15.6951,0.0305671,-1.41653e-05,2.74319e-09,-1.94357e-13,-27650.8,-43.9362], Tmin=(1144.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-201.458,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + group(Cds-OdOsH) + radical(CCCJ=O) + radical((O)CJOC)"""),
)

species(
    label = 'O=CC[CH][CH]OOC=O(22136)',
    structure = SMILES('O=CC[CH][CH]OOC=O'),
    E0 = (-171.012,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.659557,0.0757736,-6.44793e-05,2.79421e-08,-5.023e-12,-20449.9,34.4843], Tmin=(100,'K'), Tmax=(1284.15,'K')), NASAPolynomial(coeffs=[13.666,0.0352594,-1.71547e-05,3.37337e-09,-2.39875e-13,-23790.4,-31.5171], Tmin=(1284.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-171.012,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + group(Cds-OdOsH) + radical(CCsJOOC) + radical(CCJCOOH)"""),
)

species(
    label = 'O=C[CH]C[CH]OOC=O(22137)',
    structure = SMILES('[O]C=CC[CH]OOC=O'),
    E0 = (-228.175,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0468913,0.0705013,-2.20864e-05,-3.20228e-08,1.94219e-11,-27280.9,34.4276], Tmin=(100,'K'), Tmax=(986.515,'K')), NASAPolynomial(coeffs=[23.3753,0.0180222,-6.89817e-06,1.37315e-09,-1.05383e-13,-33969.8,-88.7322], Tmin=(986.515,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-228.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdOsH) + radical(C=COJ) + radical(CCsJOOC)"""),
)

species(
    label = 'O=[C]CC[CH]OOC=O(22138)',
    structure = SMILES('O=[C]CC[CH]OOC=O'),
    E0 = (-211.466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.384112,0.0787797,-7.07696e-05,3.22836e-08,-5.99493e-12,-25302.9,33.9884], Tmin=(100,'K'), Tmax=(1269.84,'K')), NASAPolynomial(coeffs=[15.4446,0.0313389,-1.47299e-05,2.86263e-09,-2.02652e-13,-29127.8,-42.2675], Tmin=(1269.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-211.466,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + group(Cds-OdOsH) + radical(CCCJ=O) + radical(CCsJOOC)"""),
)

species(
    label = '[O][C]=O(669)',
    structure = SMILES('[O][C]=O'),
    E0 = (31.9507,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.90478,-0.000175995,8.15126e-06,-1.13656e-08,4.4768e-12,3848.25,8.04855], Tmin=(100,'K'), Tmax=(975.388,'K')), NASAPolynomial(coeffs=[5.59398,-0.00122084,7.11747e-07,-9.7712e-11,3.97995e-15,3238.91,-1.49318], Tmin=(975.388,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.9507,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(OJC=O) + radical((O)CJOH)"""),
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
    label = 'O=[C]OOC1CC[CH]O1(22139)',
    structure = SMILES('O=[C]OOC1CC[CH]O1'),
    E0 = (-192.825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.450278,0.0559988,2.39669e-05,-8.7405e-08,4.27801e-11,-23043.2,27.0469], Tmin=(100,'K'), Tmax=(914.76,'K')), NASAPolynomial(coeffs=[23.4145,0.0134348,-1.10208e-06,9.82733e-13,-2.17526e-15,-29665,-94.9263], Tmin=(914.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-192.825,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cds-OdOsH) + ring(Tetrahydrofuran) + radical((O)CJOC) + radical(CCsJOCs)"""),
)

species(
    label = 'O=C1O[CH]CC[CH]OO1(22140)',
    structure = SMILES('O=C1O[CH]CC[CH]OO1'),
    E0 = (-208.622,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03561,0.0348038,8.50092e-05,-1.42324e-07,5.78451e-11,-24957.3,25.9714], Tmin=(100,'K'), Tmax=(968.549,'K')), NASAPolynomial(coeffs=[23.5094,0.0178536,-6.23111e-06,1.34912e-09,-1.13781e-13,-32869,-100.103], Tmin=(968.549,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-208.622,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-OdOsOs) + ring(Cyclooctane) + radical(CCsJOC(O)) + radical(CCsJOOC)"""),
)

species(
    label = 'O=CCC=COOC=O(22141)',
    structure = SMILES('O=CCC=COOC=O'),
    E0 = (-419.878,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.546829,0.0553943,1.03159e-05,-5.60849e-08,2.48534e-11,-50357.4,32.761], Tmin=(100,'K'), Tmax=(1029.97,'K')), NASAPolynomial(coeffs=[21.2785,0.0223571,-1.07124e-05,2.27933e-09,-1.7576e-13,-57146.3,-80.0944], Tmin=(1029.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-419.878,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-O2d)) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-OdOsH)"""),
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
    label = 'CC[CH]OO[C]=O(10703)',
    structure = SMILES('CC[CH]OO[C]=O'),
    E0 = (-68.3987,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,1855,455,950,2750,2850,1437.5,1250,1305,750,350,315.234],'cm^-1')),
        HinderedRotor(inertia=(0.559736,'amu*angstrom^2'), symmetry=1, barrier=(39.471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.559737,'amu*angstrom^2'), symmetry=1, barrier=(39.471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.559734,'amu*angstrom^2'), symmetry=1, barrier=(39.471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.559734,'amu*angstrom^2'), symmetry=1, barrier=(39.471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.559736,'amu*angstrom^2'), symmetry=1, barrier=(39.471,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3579.64,'J/mol'), sigma=(6.06667,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=559.13 K, Pc=36.38 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.726087,0.0621671,-4.32071e-05,9.6887e-09,9.17052e-13,-8100.26,26.2973], Tmin=(100,'K'), Tmax=(1105.13,'K')), NASAPolynomial(coeffs=[16.0809,0.0218918,-9.30952e-06,1.76847e-09,-1.25383e-13,-12428.5,-53.5432], Tmin=(1105.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-68.3987,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdOsH) + radical(CCsJOOC) + radical((O)CJOC)"""),
)

species(
    label = '[CH2]C(CC=O)OO[C]=O(22142)',
    structure = SMILES('[CH2]C(CC=O)OO[C]=O'),
    E0 = (-160.099,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,350,500,795,815,2782.5,750,1395,475,1775,1000,200],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0439085,0.0848007,-8.06793e-05,3.82887e-08,-7.24415e-12,-19106.5,34.6899], Tmin=(100,'K'), Tmax=(1268.75,'K')), NASAPolynomial(coeffs=[18.293,0.0269903,-1.23326e-05,2.37603e-09,-1.67821e-13,-23759.6,-58.1398], Tmin=(1268.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-160.099,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Cds-OdOsH) + radical(CJCOOH) + radical((O)CJOC)"""),
)

species(
    label = 'O=CCCC1OOC1=O(22143)',
    structure = SMILES('O=CCCC1OOC1=O'),
    E0 = (-443.471,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1401,0.0378608,6.12254e-05,-1.09455e-07,4.43734e-11,-53211.8,30.9925], Tmin=(100,'K'), Tmax=(983.353,'K')), NASAPolynomial(coeffs=[20.2078,0.0220199,-8.75965e-06,1.81965e-09,-1.43425e-13,-59946,-75.8513], Tmin=(983.353,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-443.471,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsOs) + group(Cds-OdCsH) + ring(Cyclobutane)"""),
)

species(
    label = '[O]O[CH]CCC=O(16180)',
    structure = SMILES('[O]O[CH]CCC=O'),
    E0 = (25.1613,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,492.5,1135,1000,3025,407.5,1350,352.5,180,2535.83],'cm^-1')),
        HinderedRotor(inertia=(0.292549,'amu*angstrom^2'), symmetry=1, barrier=(6.72627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.291727,'amu*angstrom^2'), symmetry=1, barrier=(6.70738,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.292008,'amu*angstrom^2'), symmetry=1, barrier=(6.71385,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.57902,'amu*angstrom^2'), symmetry=1, barrier=(59.2968,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02339,0.0746521,-0.000118135,1.098e-07,-3.98934e-11,3124.52,26.6918], Tmin=(100,'K'), Tmax=(842.043,'K')), NASAPolynomial(coeffs=[4.70518,0.0378582,-1.82034e-05,3.45591e-09,-2.36799e-13,3188.85,13.6261], Tmin=(842.043,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(25.1613,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + radical(CCsJOOH) + radical(ROOJ)"""),
)

species(
    label = 'C=COC[CH]OO[C]=O(22144)',
    structure = SMILES('C=COC[CH]OO[C]=O'),
    E0 = (-133.382,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1855,455,950,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,500,795,815,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.725878,0.0867784,-5.69826e-05,-5.26687e-09,1.25018e-11,-15856.5,32.8277], Tmin=(100,'K'), Tmax=(966.461,'K')), NASAPolynomial(coeffs=[26.6957,0.0129998,-4.11273e-06,7.81785e-10,-6.13145e-14,-23011.6,-108.126], Tmin=(966.461,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-133.382,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-OdOsH) + radical(CCsJOOC) + radical((O)CJOC)"""),
)

species(
    label = 'O=[C]OO[CH]CC=CO(22145)',
    structure = SMILES('O=[C]OO[CH]CC=CO'),
    E0 = (-173.189,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,500,795,815,200,800],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.660423,0.0851711,-5.34874e-05,-9.90724e-09,1.46928e-11,-20646.3,33.8458], Tmin=(100,'K'), Tmax=(953.377,'K')), NASAPolynomial(coeffs=[26.8101,0.0114991,-3.00134e-06,5.39838e-10,-4.35881e-14,-27774.1,-107.283], Tmin=(953.377,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-173.189,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdOsH) + radical((O)CJOC) + radical(CCsJOOC)"""),
)

species(
    label = 'CO3t2(123)',
    structure = SMILES('[O]O[C]=O'),
    E0 = (106.943,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,262.32,262.473],'cm^-1')),
        HinderedRotor(inertia=(0.00244921,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.0089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.24172,0.0164884,-1.87036e-05,1.04487e-08,-2.29623e-12,12889.7,11.1068], Tmin=(100,'K'), Tmax=(1108.89,'K')), NASAPolynomial(coeffs=[6.67163,0.00411578,-1.96698e-06,3.86453e-10,-2.76577e-14,12129,-5.79498], Tmin=(1108.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(106.943,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(78.9875,'J/(mol*K)'), label="""CO3t2""", comment="""Thermo library: primaryThermoLibrary"""),
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
    label = '[C]=O(361)',
    structure = SMILES('[C]=O'),
    E0 = (439.086,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3054.48],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0101,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.08916,0.00200416,-1.61661e-05,2.55058e-08,-1.16424e-11,52802.7,4.52505], Tmin=(100,'K'), Tmax=(856.11,'K')), NASAPolynomial(coeffs=[0.961625,0.00569045,-3.48044e-06,7.19202e-10,-5.08041e-14,53738.7,21.4663], Tmin=(856.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(439.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CdCdJ2_triplet)"""),
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
    E0 = (-174.978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-75.1952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-105.106,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-4.90133,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-50.7693,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-49.4583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-47.3663,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-78.3279,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-131.728,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-91.5685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-95.7427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (51.0228,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-116.821,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-124.093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-150.005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (170.114,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-1.47263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-166.694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-81.5056,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (180.418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-29.3226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (328.054,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (464.247,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=[C]OO[CH]CCC=O(21319)'],
    products = ['CO2(13)', 'O=CCCC=O(5767)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O=[C]OO[CH]CCC=O(21319)'],
    products = ['[O]C1CCC1OO[C]=O(22130)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.13771e+06,'s^-1'), n=1.58803, Ea=(99.7832,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_csHNd] for rate rule [R5_SS_CO;carbonylbond_intra_H;radadd_intra_csHNd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 96.1 to 99.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['O=[C]OO[CH]CCC=O(21319)'],
    products = ['[O]C1CC[CH]OOC1=O(22131)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.19e+11,'s^-1'), n=0.08, Ea=(69.8728,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R8;multiplebond_intra;radadd_intra] for rate rule [R8;carbonylbond_intra_H;radadd_intra_CO]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'O=[C]OOC=CCC=O(22132)'],
    products = ['O=[C]OO[CH]CCC=O(21319)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.72e+08,'cm^3/(mol*s)'), n=1.477, Ea=(6.73624,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2825 used for Cds-CsH_Cds-OsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-OsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2CHO(40)', 'C=COO[C]=O(745)'],
    products = ['O=[C]OO[CH]CCC=O(21319)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.0114756,'m^3/(mol*s)'), n=2.44484, Ea=(36.1431,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-HH_Cds;CsJ-OneDeHH] for rate rule [Cds-HH_Cds-OsH;CsJ-COHH]
Euclidian distance = 1.41421356237
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O=[C]OO[CH]CCC=O(21319)'],
    products = ['O=[C]OOC[CH]CC=O(22133)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(9.4e-19,'s^-1'), n=8.84, Ea=(125.52,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 342 used for R2H_S;C_rad_out_H/NonDeO;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeO;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O=[C]OO[CH]CCC=O(21319)'],
    products = ['O=[C]OOCC[CH]C=O(22134)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6.82e+09,'s^-1'), n=0.73, Ea=(127.612,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_1H;Cs_H_out_H/OneDe] for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeO;Cs_H_out_H/CO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O=[C]OO[CH]CCC=O(21319)'],
    products = ['O=[C]CCCOO[C]=O(22135)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(27900,'s^-1'), n=1.97, Ea=(96.6504,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/NonDeO;XH_out] for rate rule [R4H_SSS;C_rad_out_H/NonDeO;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['O=[C]OO[CH]CCC=O(21319)'],
    products = ['O=CC[CH][CH]OOC=O(22136)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(379583,'s^-1'), n=1.54051, Ea=(43.2505,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5HJ_3;CO_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['O=[C]OO[CH]CCC=O(21319)'],
    products = ['O=C[CH]C[CH]OOC=O(22137)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(49528.1,'s^-1'), n=1.95205, Ea=(83.4098,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R6HJ_3;CO_rad_out;Cs_H_out_H/CO]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['O=[C]OO[CH]CCC=O(21319)'],
    products = ['O=[C]CC[CH]OOC=O(22138)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(364667,'s^-1'), n=1.22214, Ea=(79.2357,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7Hall;Y_rad_out;XH_out] for rate rule [R7HJ_3;CO_rad_out;CO_H_out]
Euclidian distance = 1.73205080757
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O][C]=O(669)', '[O][CH]CCC=O(5764)'],
    products = ['O=[C]OO[CH]CCC=O(21319)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.49215e+07,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['O=[C]OO[CH]CCC=O(21319)'],
    products = ['O=[C]OOC1CC[CH]O1(22139)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.64e+06,'s^-1'), n=1.15, Ea=(58.1576,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_csHNd] for rate rule [R5_SS_CO;carbonyl_intra_H;radadd_intra_csHO]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction14',
    reactants = ['O=[C]OO[CH]CCC=O(21319)'],
    products = ['O=C1O[CH]CC[CH]OO1(22140)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(5.649e+12,'s^-1'), n=0.287, Ea=(50.8858,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;carbonyl_intra_H;radadd_intra_CO] for rate rule [R8_linear;carbonyl_intra_H;radadd_intra_CO]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O=[C]OO[CH]CCC=O(21319)'],
    products = ['O=CCC=COOC=O(22141)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CO(12)', 'CC[CH]OO[C]=O(10703)'],
    products = ['O=[C]OO[CH]CCC=O(21319)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(274200,'cm^3/(mol*s)'), n=2.53, Ea=(357.732,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 4 used for CO;C_pri/NonDeC
Exact match found for rate rule [CO;C_pri/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(CC=O)OO[C]=O(22142)'],
    products = ['O=[C]OO[CH]CCC=O(21319)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction18',
    reactants = ['O=[C]OO[CH]CCC=O(21319)'],
    products = ['O=CCCC1OOC1=O(22143)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CO(12)', '[O]O[CH]CCC=O(16180)'],
    products = ['O=[C]OO[CH]CCC=O(21319)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.41e+07,'cm^3/(mol*s)'), n=0, Ea=(12.552,'kJ/mol'), T0=(1,'K'), Tmin=(250,'K'), Tmax=(2500,'K'), comment="""From training reaction 9 used for COm;O_rad/NonDe
Exact match found for rate rule [COm;O_rad/NonDe]
Euclidian distance = 0
family: R_Addition_COm"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=COC[CH]OO[C]=O(22144)'],
    products = ['O=[C]OO[CH]CCC=O(21319)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 7 used for R_ROR;R1_doublebond_CH2;R2_doublebond_H;R_O_C
Exact match found for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction21',
    reactants = ['O=[C]OO[CH]CC=CO(22145)'],
    products = ['O=[C]OO[CH]CCC=O(21319)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_H;R_O_H]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CO3t2(123)', '[CH]CCC=O(18484)'],
    products = ['O=[C]OO[CH]CCC=O(21319)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['[C]=O(361)', '[O]O[CH]CCC=O(16180)'],
    products = ['O=[C]OO[CH]CCC=O(21319)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

network(
    label = '3726',
    isomers = [
        'O=[C]OO[CH]CCC=O(21319)',
    ],
    reactants = [
        ('CO2(13)', 'O=CCCC=O(5767)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '3726',
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

