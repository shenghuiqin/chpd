species(
    label = 'CCC(=O)[C](O)CO[CH]CCC=O(18458)',
    structure = SMILES('CCC([O])=C(O)CO[CH]CCC=O'),
    E0 = (-463.701,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,2782.5,750,1395,475,1775,1000,180,180,180,180,1600,1754.98,2750.54,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156119,'amu*angstrom^2'), symmetry=1, barrier=(3.58948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156119,'amu*angstrom^2'), symmetry=1, barrier=(3.58948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156119,'amu*angstrom^2'), symmetry=1, barrier=(3.58948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156119,'amu*angstrom^2'), symmetry=1, barrier=(3.58948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156119,'amu*angstrom^2'), symmetry=1, barrier=(3.58948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156119,'amu*angstrom^2'), symmetry=1, barrier=(3.58948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156119,'amu*angstrom^2'), symmetry=1, barrier=(3.58948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156119,'amu*angstrom^2'), symmetry=1, barrier=(3.58948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156119,'amu*angstrom^2'), symmetry=1, barrier=(3.58948,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.12452,0.160107,-0.00018716,1.16833e-07,-2.924e-11,-55516.6,51.9766], Tmin=(100,'K'), Tmax=(972.109,'K')), NASAPolynomial(coeffs=[22.7214,0.0537562,-2.30556e-05,4.29052e-09,-2.96751e-13,-60541.6,-71.9837], Tmin=(972.109,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-463.701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(C=C(C)OJ) + radical(CCsJOCs)"""),
)

species(
    label = 'C=C(O)C(=O)CC(4626)',
    structure = SMILES('C=C(O)C(=O)CC'),
    E0 = (-358.468,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,332.629,332.633],'cm^-1')),
        HinderedRotor(inertia=(0.152178,'amu*angstrom^2'), symmetry=1, barrier=(11.9524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152152,'amu*angstrom^2'), symmetry=1, barrier=(11.9478,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152157,'amu*angstrom^2'), symmetry=1, barrier=(11.9534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152209,'amu*angstrom^2'), symmetry=1, barrier=(11.954,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4071.48,'J/mol'), sigma=(6.49965,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=635.96 K, Pc=33.65 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.785843,0.0702816,-6.54766e-05,3.23543e-08,-6.53762e-12,-42997.7,22.9463], Tmin=(100,'K'), Tmax=(1175.98,'K')), NASAPolynomial(coeffs=[12.9689,0.0288414,-1.26178e-05,2.38822e-09,-1.6711e-13,-45863,-37.8049], Tmin=(1175.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-358.468,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsHH)"""),
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
    label = '[CH2]C(O)=C([O])CC(4557)',
    structure = SMILES('[CH2]C(O)=C([O])CC'),
    E0 = (-184.943,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,402.631,402.631],'cm^-1')),
        HinderedRotor(inertia=(0.100035,'amu*angstrom^2'), symmetry=1, barrier=(11.5078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100034,'amu*angstrom^2'), symmetry=1, barrier=(11.5078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100035,'amu*angstrom^2'), symmetry=1, barrier=(11.5078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100035,'amu*angstrom^2'), symmetry=1, barrier=(11.5078,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4303.36,'J/mol'), sigma=(7.05412,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=672.18 K, Pc=27.82 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.516209,0.0797578,-7.96062e-05,3.90546e-08,-7.19184e-12,-22064.1,30.4076], Tmin=(100,'K'), Tmax=(1528.49,'K')), NASAPolynomial(coeffs=[20.9225,0.0119894,-1.65417e-06,6.23996e-11,2.30461e-15,-27255.3,-77.6604], Tmin=(1528.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-184.943,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(O)CJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CCC([O])=C(O)COC1CCC1[O](19042)',
    structure = SMILES('CCC([O])=C(O)COC1CCC1[O]'),
    E0 = (-357.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.8157,0.129945,-8.7459e-05,1.11882e-08,7.4445e-12,-42786.5,50.0702], Tmin=(100,'K'), Tmax=(1008.63,'K')), NASAPolynomial(coeffs=[29.4985,0.041676,-1.54986e-05,2.82679e-09,-1.99499e-13,-51333.8,-116.162], Tmin=(1008.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-357.933,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(631.9,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(Cyclobutane) + radical(CC(C)OJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CC[C]([O])C1(O)COC1CCC=O(19043)',
    structure = SMILES('CC[C]([O])C1(O)COC1CCC=O'),
    E0 = (-290.379,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.79579,0.149089,-0.000164123,1.01127e-07,-2.51675e-11,-34679.5,48.0454], Tmin=(100,'K'), Tmax=(977.461,'K')), NASAPolynomial(coeffs=[19.9569,0.0559802,-2.12386e-05,3.67398e-09,-2.42581e-13,-39127.5,-61.2043], Tmin=(977.461,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-290.379,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(627.743,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Oxetane) + radical(C2CsJOH) + radical(CC(C)OJ)"""),
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
    label = 'CCC([O])=C(O)COC=CCC=O(19044)',
    structure = SMILES('CCC([O])=C(O)COC=CCC=O'),
    E0 = (-561.22,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2782.5,750,1395,475,1775,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    molecularWeight = (185.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.19151,0.134998,-9.41319e-05,8.76084e-09,1.02893e-11,-67219.7,50.238], Tmin=(100,'K'), Tmax=(999.932,'K')), NASAPolynomial(coeffs=[34.672,0.0318522,-1.1885e-05,2.25081e-09,-1.65143e-13,-77207.5,-144.508], Tmin=(999.932,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-561.22,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(598.642,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=C(C)OJ)"""),
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
    label = 'C=COCC(O)=C([O])CC(6482)',
    structure = SMILES('C=COCC(O)=C([O])CC'),
    E0 = (-425.931,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.16,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.72064,0.105433,-6.15504e-05,-1.84455e-08,2.14357e-11,-51002.9,38.6533], Tmin=(100,'K'), Tmax=(928.112,'K')), NASAPolynomial(coeffs=[30.991,0.0168581,-3.09686e-06,4.06711e-10,-3.04388e-14,-59332.1,-128.881], Tmin=(928.112,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-425.931,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
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
    label = 'O=C=C(O)CO[CH]CCC=O(19045)',
    structure = SMILES('O=C=C(O)CO[CH]CCC=O'),
    E0 = (-364.837,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,3025,407.5,1350,352.5,3615,1277.5,1000,350,440,435,1725,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2782.5,750,1395,475,1775,1000,200,800,1000,1200,1400,1600],'cm^-1')),
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
    molecularWeight = (157.144,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.41081,0.147651,-0.000219371,1.68856e-07,-5.09009e-11,-43655.2,40.8852], Tmin=(100,'K'), Tmax=(851.627,'K')), NASAPolynomial(coeffs=[20.0794,0.0371052,-1.60106e-05,2.89001e-09,-1.92486e-13,-47307.7,-62.9592], Tmin=(851.627,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-364.837,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cdd-O2d)OsHH) + group(Cds-(Cdd-O2d)CsOs) + group(Cds-OdCsH) + radical(CCsJOCs)"""),
)

species(
    label = 'CCC([O])=C(O)COC[CH]CC=O(19046)',
    structure = SMILES('CCC([O])=C(O)COC[CH]CC=O'),
    E0 = (-444.255,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,2782.5,750,1395,475,1775,1000,180,180,180,180,1600,1773.61,2731.34,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156566,'amu*angstrom^2'), symmetry=1, barrier=(3.59976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156566,'amu*angstrom^2'), symmetry=1, barrier=(3.59976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156566,'amu*angstrom^2'), symmetry=1, barrier=(3.59976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156566,'amu*angstrom^2'), symmetry=1, barrier=(3.59976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156566,'amu*angstrom^2'), symmetry=1, barrier=(3.59976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156566,'amu*angstrom^2'), symmetry=1, barrier=(3.59976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156566,'amu*angstrom^2'), symmetry=1, barrier=(3.59976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156566,'amu*angstrom^2'), symmetry=1, barrier=(3.59976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156566,'amu*angstrom^2'), symmetry=1, barrier=(3.59976,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.57241,0.147236,-0.000157224,9.09477e-08,-2.13881e-11,-53197,53.044], Tmin=(100,'K'), Tmax=(1024.63,'K')), NASAPolynomial(coeffs=[20.6488,0.056584,-2.45143e-05,4.60075e-09,-3.20232e-13,-57955.6,-59.5496], Tmin=(1024.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-444.255,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(CCJCO) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CCC([O])=C(O)[CH]OCCCC=O(19047)',
    structure = SMILES('CCC([O])=C(O)[CH]OCCCC=O'),
    E0 = (-533.218,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.43675,0.148528,-0.000164647,1.02222e-07,-2.61999e-11,-63905.4,51.0073], Tmin=(100,'K'), Tmax=(937.047,'K')), NASAPolynomial(coeffs=[17.47,0.0635526,-2.86227e-05,5.44783e-09,-3.81327e-13,-67636.2,-43.7371], Tmin=(937.047,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-533.218,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(C=CCJ(O)C) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CCC([O])=C(O)COCC[CH]C=O(19048)',
    structure = SMILES('CCC([O])=C(O)COCCC=C[O]'),
    E0 = (-500.906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.77,0.153099,-0.000153165,7.72576e-08,-1.52034e-11,-59950.2,54.5478], Tmin=(100,'K'), Tmax=(1246.5,'K')), NASAPolynomial(coeffs=[32.5851,0.0364355,-1.27755e-05,2.17262e-09,-1.44163e-13,-69013.5,-128.855], Tmin=(1246.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-500.906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(623.585,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C[CH]C(O)=C(O)CO[CH]CCC=O(19049)',
    structure = SMILES('C[CH]C(O)=C(O)CO[CH]CCC=O'),
    E0 = (-401.604,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2782.5,750,1395,475,1775,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
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
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.5251,0.161714,-0.000181601,1.04724e-07,-2.38077e-11,-48027.4,54.3987], Tmin=(100,'K'), Tmax=(1075.55,'K')), NASAPolynomial(coeffs=[28.0319,0.0443508,-1.79181e-05,3.26514e-09,-2.24337e-13,-54815.5,-100.143], Tmin=(1075.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-401.604,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(615.271,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(CCJCO) + radical(CCsJOCs)"""),
)

species(
    label = 'CCC([O])=C(O)COCCC[C]=O(19050)',
    structure = SMILES('CCC([O])=C(O)COCCC[C]=O'),
    E0 = (-484.197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.85825,0.155569,-0.000181156,1.14718e-07,-2.93539e-11,-57992.4,52.3966], Tmin=(100,'K'), Tmax=(948.33,'K')), NASAPolynomial(coeffs=[20.6285,0.0565043,-2.44656e-05,4.56776e-09,-3.1621e-13,-62447.1,-59.6676], Tmin=(948.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-484.197,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(C=C(C)OJ) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]CC(O)=C(O)CO[CH]CCC=O(19051)',
    structure = SMILES('[CH2]CC(O)=C(O)CO[CH]CCC=O'),
    E0 = (-396.26,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,325,375,415,465,420,450,1700,1750,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,2782.5,750,1395,475,1775,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
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
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.54163,0.165928,-0.000194301,1.17821e-07,-2.82641e-11,-47387.3,53.7537], Tmin=(100,'K'), Tmax=(1018.9,'K')), NASAPolynomial(coeffs=[26.6476,0.0474098,-1.9821e-05,3.65749e-09,-2.52411e-13,-53539.2,-92.4569], Tmin=(1018.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-396.26,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(615.271,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(RCCJ) + radical(CCsJOCs)"""),
)

species(
    label = 'CCC(O)=C([O])CO[CH]CCC=O(19052)',
    structure = SMILES('CCC(O)=C([O])CO[CH]CCC=O'),
    E0 = (-463.701,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,2782.5,750,1395,475,1775,1000,180,180,180,180,1600,1754.98,2750.54,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156119,'amu*angstrom^2'), symmetry=1, barrier=(3.58948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156119,'amu*angstrom^2'), symmetry=1, barrier=(3.58948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156119,'amu*angstrom^2'), symmetry=1, barrier=(3.58948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156119,'amu*angstrom^2'), symmetry=1, barrier=(3.58948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156119,'amu*angstrom^2'), symmetry=1, barrier=(3.58948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156119,'amu*angstrom^2'), symmetry=1, barrier=(3.58948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156119,'amu*angstrom^2'), symmetry=1, barrier=(3.58948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156119,'amu*angstrom^2'), symmetry=1, barrier=(3.58948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156119,'amu*angstrom^2'), symmetry=1, barrier=(3.58948,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.12452,0.160107,-0.00018716,1.16833e-07,-2.924e-11,-55516.6,51.9766], Tmin=(100,'K'), Tmax=(972.109,'K')), NASAPolynomial(coeffs=[22.7214,0.0537562,-2.30556e-05,4.29052e-09,-2.96751e-13,-60541.6,-71.9837], Tmin=(972.109,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-463.701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(CCsJOCs) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CCC(O)=C(O)[CH]O[CH]CCC=O(19053)',
    structure = SMILES('CCC(O)=C(O)[CH]O[CH]CCC=O'),
    E0 = (-490.567,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.28537,0.161823,-0.000185085,1.1114e-07,-2.66379e-11,-58740.4,51.9858], Tmin=(100,'K'), Tmax=(1015.22,'K')), NASAPolynomial(coeffs=[24.702,0.0515473,-2.21451e-05,4.13791e-09,-2.87414e-13,-64422.9,-83.4589], Tmin=(1015.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-490.567,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(615.271,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(CCsJOCs) + radical(C=CCJ(O)C)"""),
)

species(
    label = 'CCC([O])=C([O])COCCCC=O(19054)',
    structure = SMILES('CCC([O])=C([O])COCCCC=O'),
    E0 = (-506.352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.35929,0.147824,-0.00017039,1.12856e-07,-3.09959e-11,-60678.1,51.2951], Tmin=(100,'K'), Tmax=(876.392,'K')), NASAPolynomial(coeffs=[15.6978,0.0654088,-2.93309e-05,5.55288e-09,-3.86735e-13,-63843.1,-33.4375], Tmin=(876.392,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-506.352,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(623.585,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(C=C(C)OJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C[CH]C([O])=C(O)COCCCC=O(19055)',
    structure = SMILES('C[CH]C([O])=C(O)COCCCC=O'),
    E0 = (-444.255,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,2782.5,750,1395,475,1775,1000,180,180,180,180,1600,1773.61,2731.34,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156566,'amu*angstrom^2'), symmetry=1, barrier=(3.59976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156566,'amu*angstrom^2'), symmetry=1, barrier=(3.59976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156566,'amu*angstrom^2'), symmetry=1, barrier=(3.59976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156566,'amu*angstrom^2'), symmetry=1, barrier=(3.59976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156566,'amu*angstrom^2'), symmetry=1, barrier=(3.59976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156566,'amu*angstrom^2'), symmetry=1, barrier=(3.59976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156566,'amu*angstrom^2'), symmetry=1, barrier=(3.59976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156566,'amu*angstrom^2'), symmetry=1, barrier=(3.59976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156566,'amu*angstrom^2'), symmetry=1, barrier=(3.59976,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.57241,0.147236,-0.000157224,9.09477e-08,-2.13881e-11,-53197,53.044], Tmin=(100,'K'), Tmax=(1024.63,'K')), NASAPolynomial(coeffs=[20.6488,0.056584,-2.45143e-05,4.60075e-09,-3.20232e-13,-57955.6,-59.5496], Tmin=(1024.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-444.255,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(CCJCO) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CCC(O)=C(O)CO[CH][CH]CC=O(19056)',
    structure = SMILES('CCC(O)=C(O)CO[CH][CH]CC=O'),
    E0 = (-401.604,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2782.5,750,1395,475,1775,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
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
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.5251,0.161714,-0.000181601,1.04724e-07,-2.38077e-11,-48027.4,54.3987], Tmin=(100,'K'), Tmax=(1075.55,'K')), NASAPolynomial(coeffs=[28.0319,0.0443508,-1.79181e-05,3.26514e-09,-2.24337e-13,-54815.5,-100.143], Tmin=(1075.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-401.604,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(615.271,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(CCJCO) + radical(CCsJOCs)"""),
)

species(
    label = '[CH2]CC([O])=C(O)COCCCC=O(19057)',
    structure = SMILES('[CH2]CC([O])=C(O)COCCCC=O'),
    E0 = (-438.911,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,2782.5,750,1395,475,1775,1000,3615,1277.5,1000,3000,3100,440,815,1455,1000,325,375,415,465,420,450,1700,1750,180,180,180,180,1600,1733.44,2772.88,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155584,'amu*angstrom^2'), symmetry=1, barrier=(3.57719,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155584,'amu*angstrom^2'), symmetry=1, barrier=(3.57719,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155584,'amu*angstrom^2'), symmetry=1, barrier=(3.57719,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155584,'amu*angstrom^2'), symmetry=1, barrier=(3.57719,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155584,'amu*angstrom^2'), symmetry=1, barrier=(3.57719,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155584,'amu*angstrom^2'), symmetry=1, barrier=(3.57719,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155584,'amu*angstrom^2'), symmetry=1, barrier=(3.57719,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155584,'amu*angstrom^2'), symmetry=1, barrier=(3.57719,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155584,'amu*angstrom^2'), symmetry=1, barrier=(3.57719,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.67105,0.152387,-0.000173066,1.07936e-07,-2.74358e-11,-52553.3,52.6954], Tmin=(100,'K'), Tmax=(950.297,'K')), NASAPolynomial(coeffs=[19.4203,0.0593991,-2.62859e-05,4.96376e-09,-3.45981e-13,-56751.9,-52.7562], Tmin=(950.297,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-438.911,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(C=C(C)OJ) + radical(RCCJ)"""),
)

species(
    label = 'CCC(O)=C(O)CO[CH]C[CH]C=O(19058)',
    structure = SMILES('CCC(O)=C(O)CO[CH]CC=C[O]'),
    E0 = (-458.254,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,180,180,180,180,1600,1614.6,2872.29,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152295,'amu*angstrom^2'), symmetry=1, barrier=(3.50156,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152295,'amu*angstrom^2'), symmetry=1, barrier=(3.50156,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152295,'amu*angstrom^2'), symmetry=1, barrier=(3.50156,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152295,'amu*angstrom^2'), symmetry=1, barrier=(3.50156,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152295,'amu*angstrom^2'), symmetry=1, barrier=(3.50156,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152295,'amu*angstrom^2'), symmetry=1, barrier=(3.50156,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152295,'amu*angstrom^2'), symmetry=1, barrier=(3.50156,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152295,'amu*angstrom^2'), symmetry=1, barrier=(3.50156,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152295,'amu*angstrom^2'), symmetry=1, barrier=(3.50156,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.81503,0.168536,-0.000180371,9.41207e-08,-1.87472e-11,-54776.3,56.2434], Tmin=(100,'K'), Tmax=(1308.92,'K')), NASAPolynomial(coeffs=[39.5705,0.0248913,-6.58352e-06,9.33894e-10,-5.64011e-14,-65710,-167.221], Tmin=(1308.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-458.254,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CCsJOCs)"""),
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
    label = 'CCC([O])=[C]O(4619)',
    structure = SMILES('CCC([O])=[C]O'),
    E0 = (-62.324,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,1685,370,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,297.522,297.522],'cm^-1')),
        HinderedRotor(inertia=(0.185706,'amu*angstrom^2'), symmetry=1, barrier=(11.6652,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185707,'amu*angstrom^2'), symmetry=1, barrier=(11.6652,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185706,'amu*angstrom^2'), symmetry=1, barrier=(11.6652,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.423469,0.0614354,-6.11557e-05,2.97291e-08,-5.38234e-12,-7352.28,27.4708], Tmin=(100,'K'), Tmax=(1581.22,'K')), NASAPolynomial(coeffs=[17.1254,0.00842786,-6.66531e-07,-7.63516e-11,1.03047e-14,-11289.4,-56.5069], Tmin=(1581.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-62.324,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=CJO) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]O[CH]CCC=O(18485)',
    structure = SMILES('[CH2]O[CH]CCC=O'),
    E0 = (2.04251,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.787429,0.0749119,-8.90054e-05,6.27045e-08,-1.84232e-11,357.485,28.3533], Tmin=(100,'K'), Tmax=(820.645,'K')), NASAPolynomial(coeffs=[9.05276,0.0346248,-1.53673e-05,2.88312e-09,-1.9927e-13,-999.093,-9.88836], Tmin=(820.645,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.04251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-OsHHH) + group(Cds-OdCsH) + radical(CCsJOCs) + radical(CsJOCC)"""),
)

species(
    label = 'CC[C]1OC1(O)CO[CH]CCC=O(19059)',
    structure = SMILES('CC[C]1OC1(O)CO[CH]CCC=O'),
    E0 = (-322.23,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.77107,0.163514,-0.000176238,9.09466e-08,-1.56376e-11,-38468.4,50.886], Tmin=(100,'K'), Tmax=(889.081,'K')), NASAPolynomial(coeffs=[29.7061,0.0396478,-1.23879e-05,1.92472e-09,-1.20775e-13,-45478.4,-112.632], Tmin=(889.081,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-322.23,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Ethylene_oxide) + radical(CCsJOCs) + radical(C2CsJO)"""),
)

species(
    label = 'CCC([O])=C(O)COC1CC[CH]O1(19060)',
    structure = SMILES('CCC([O])=C(O)COC1CC[CH]O1'),
    E0 = (-475.563,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.07798,0.144557,-0.000134089,6.3491e-08,-1.1498e-11,-56877.9,51.737], Tmin=(100,'K'), Tmax=(1525.19,'K')), NASAPolynomial(coeffs=[32.1013,0.0327631,-7.51294e-06,8.95071e-10,-4.6196e-14,-65947.3,-131.632], Tmin=(1525.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-475.563,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(631.9,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(Tetrahydrofuran) + radical(CCsJOCs) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CCC1([O])[C](O)COC1CCC=O(19061)',
    structure = SMILES('CCC1([O])[C](O)COC1CCC=O'),
    E0 = (-370.286,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.42201,0.137714,-0.000134242,7.34417e-08,-1.64028e-11,-44300.9,44.7537], Tmin=(100,'K'), Tmax=(1078.78,'K')), NASAPolynomial(coeffs=[19.4585,0.0565849,-2.14378e-05,3.73171e-09,-2.4826e-13,-49021.8,-62.4666], Tmin=(1078.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-370.286,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(631.9,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Tetrahydrofuran) + radical(CC(C)2OJ) + radical(C2CsJOH)"""),
)

species(
    label = 'CCC1OO[CH]CC[CH]OCC=1O(19062)',
    structure = SMILES('CCC1OO[CH]CC[CH]OCC=1O'),
    E0 = (-100.242,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.13522,0.162972,-0.000184302,1.11179e-07,-2.71742e-11,-11804.3,34.1146], Tmin=(100,'K'), Tmax=(987.794,'K')), NASAPolynomial(coeffs=[22.3795,0.0596507,-2.74033e-05,5.28646e-09,-3.73503e-13,-16844.9,-88.6653], Tmin=(987.794,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-100.242,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(644.372,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(Cyclodecane) + radical(CCsJOOC) + radical(CCsJOCs)"""),
)

species(
    label = 'CCC(O)=C(O)COC=CCC=O(19063)',
    structure = SMILES('CCC(O)=C(O)COC=CCC=O'),
    E0 = (-699.024,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.53219,0.13868,-8.56414e-05,-9.51919e-09,1.87911e-11,-83777.9,49.9591], Tmin=(100,'K'), Tmax=(980.487,'K')), NASAPolynomial(coeffs=[38.1968,0.0288494,-1.00309e-05,1.90623e-09,-1.4359e-13,-94864.5,-165.344], Tmin=(980.487,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-699.024,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH)"""),
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
    label = 'CC([O])=C(O)CO[CH]CCC=O(19064)',
    structure = SMILES('CC([O])=C(O)CO[CH]CCC=O'),
    E0 = (-441.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2782.5,750,1395,475,1775,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    molecularWeight = (172.178,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.54489,0.149043,-0.000187526,1.26027e-07,-3.38296e-11,-52815.5,46.9842], Tmin=(100,'K'), Tmax=(909.953,'K')), NASAPolynomial(coeffs=[20.4643,0.0478959,-2.07882e-05,3.86553e-09,-2.66242e-13,-57002.8,-61.8499], Tmin=(909.953,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-441.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(CCsJOCs) + radical(C=C(C)OJ)"""),
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
    label = 'CC[CH]OC[C](O)C(=O)CC(11983)',
    structure = SMILES('CC[CH]OCC(O)=C([O])CC'),
    E0 = (-357.122,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    molecularWeight = (158.195,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4866.2,'J/mol'), sigma=(8.19114,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=760.09 K, Pc=20.09 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.56208,0.137964,-0.000145011,7.93795e-08,-1.72027e-11,-42710,45.6662], Tmin=(100,'K'), Tmax=(1125.41,'K')), NASAPolynomial(coeffs=[24.7296,0.0409615,-1.57209e-05,2.79017e-09,-1.88876e-13,-48852.9,-89.2245], Tmin=(1125.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-357.122,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(CCsJOCs) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C(CC=O)OCC(O)=C([O])CC(19065)',
    structure = SMILES('[CH2]C(CC=O)OCC(O)=C([O])CC'),
    E0 = (-444.198,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3100,440,815,1455,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.16199,0.163426,-0.000199389,1.3131e-07,-3.47727e-11,-53171.9,52.0771], Tmin=(100,'K'), Tmax=(919.262,'K')), NASAPolynomial(coeffs=[21.4412,0.0563711,-2.47038e-05,4.62596e-09,-3.20283e-13,-57695.2,-64.548], Tmin=(919.262,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-444.198,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(C=C(C)OJ) + radical(CJC(C)OC)"""),
)

species(
    label = 'CCC1OC(CCC=O)OCC=1O(19066)',
    structure = SMILES('CCC1OC(CCC=O)OCC=1O'),
    E0 = (-731.194,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.43419,0.118141,-5.01709e-05,-3.24276e-08,2.48187e-11,-87689.5,43.6734], Tmin=(100,'K'), Tmax=(956.465,'K')), NASAPolynomial(coeffs=[30.0149,0.0382111,-1.2289e-05,2.13583e-09,-1.51104e-13,-96448,-124.767], Tmin=(956.465,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-731.194,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(631.9,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + ring(24dihydro13dioxin)"""),
)

species(
    label = 'CCC(=[C]CO[CH]CCC=O)OO(19067)',
    structure = SMILES('CCC(=[C]CO[CH]CCC=O)OO'),
    E0 = (-24.6109,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,3615,1310,387.5,850,1000,350,440,435,1725,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,2782.5,750,1395,475,1775,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
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
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.03482,0.16812,-0.000242399,2.0417e-07,-6.99443e-11,-2719.52,54.5987], Tmin=(100,'K'), Tmax=(794.786,'K')), NASAPolynomial(coeffs=[12.9454,0.0729172,-3.48339e-05,6.67192e-09,-4.62751e-13,-4792.97,-15.8902], Tmin=(794.786,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-24.6109,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(615.271,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCsJOCs) + radical(Cds_S)"""),
)

species(
    label = 'CCC([O])C(=O)CO[CH]CCC=O(19068)',
    structure = SMILES('CCC([O])C(=O)CO[CH]CCC=O'),
    E0 = (-320.082,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.34683,0.148362,-0.000165744,1.04445e-07,-2.73681e-11,-38275.9,52.7337], Tmin=(100,'K'), Tmax=(913.951,'K')), NASAPolynomial(coeffs=[16.4555,0.06607,-3.06799e-05,5.92263e-09,-4.17768e-13,-41712.7,-36.2842], Tmin=(913.951,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-320.082,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCsJOCs)"""),
)

species(
    label = 'C=COC[CH]OCC(O)=C([O])CC(19069)',
    structure = SMILES('C=COC[CH]OCC(O)=C([O])CC'),
    E0 = (-422.105,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,180,180,180,180,1017.27,1190.82,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.159521,'amu*angstrom^2'), symmetry=1, barrier=(3.66769,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159521,'amu*angstrom^2'), symmetry=1, barrier=(3.66769,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159521,'amu*angstrom^2'), symmetry=1, barrier=(3.66769,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159521,'amu*angstrom^2'), symmetry=1, barrier=(3.66769,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159521,'amu*angstrom^2'), symmetry=1, barrier=(3.66769,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159521,'amu*angstrom^2'), symmetry=1, barrier=(3.66769,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159521,'amu*angstrom^2'), symmetry=1, barrier=(3.66769,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159521,'amu*angstrom^2'), symmetry=1, barrier=(3.66769,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159521,'amu*angstrom^2'), symmetry=1, barrier=(3.66769,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.74354,0.171075,-0.000187965,1.01364e-07,-2.09496e-11,-50434.5,54.8195], Tmin=(100,'K'), Tmax=(1231.52,'K')), NASAPolynomial(coeffs=[38.3735,0.0270439,-7.67787e-06,1.14015e-09,-7.03533e-14,-60752.1,-160.947], Tmin=(1231.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-422.105,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CCsJOCs) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CCC([O])=C(O)CO[CH]CC=CO(19070)',
    structure = SMILES('CCC([O])=C(O)CO[CH]CC=CO'),
    E0 = (-461.912,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,180,180,180,180,1600,1674.75,2818.87,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153917,'amu*angstrom^2'), symmetry=1, barrier=(3.53886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153917,'amu*angstrom^2'), symmetry=1, barrier=(3.53886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153917,'amu*angstrom^2'), symmetry=1, barrier=(3.53886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153917,'amu*angstrom^2'), symmetry=1, barrier=(3.53886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153917,'amu*angstrom^2'), symmetry=1, barrier=(3.53886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153917,'amu*angstrom^2'), symmetry=1, barrier=(3.53886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153917,'amu*angstrom^2'), symmetry=1, barrier=(3.53886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153917,'amu*angstrom^2'), symmetry=1, barrier=(3.53886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153917,'amu*angstrom^2'), symmetry=1, barrier=(3.53886,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.76856,0.170553,-0.000188364,1.01938e-07,-2.10637e-11,-55220.4,56.1611], Tmin=(100,'K'), Tmax=(1258.1,'K')), NASAPolynomial(coeffs=[38.5437,0.0254302,-6.4944e-06,8.80066e-10,-5.10563e-14,-65531.7,-160.407], Tmin=(1258.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-461.912,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(CCsJOCs)"""),
)

species(
    label = 'CCC([O])=C(O)C[O](4968)',
    structure = SMILES('CCC([O])=C(O)C[O]'),
    E0 = (-271.29,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,325,375,415,465,420,450,1700,1750,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,419.848,419.884,419.898,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0894139,'amu*angstrom^2'), symmetry=1, barrier=(11.1864,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0894125,'amu*angstrom^2'), symmetry=1, barrier=(11.1858,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0894308,'amu*angstrom^2'), symmetry=1, barrier=(11.1862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0894286,'amu*angstrom^2'), symmetry=1, barrier=(11.1864,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.258431,0.0802067,-8.19403e-05,4.47079e-08,-9.79359e-12,-32492.2,31.5468], Tmin=(100,'K'), Tmax=(1105.45,'K')), NASAPolynomial(coeffs=[14.6794,0.0280257,-1.11354e-05,2.00761e-09,-1.36852e-13,-35680.5,-39.4716], Tmin=(1105.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-271.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(CCOJ) + radical(C=C(C)OJ)"""),
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
    label = 'CC[C]=C(O)CO[CH]CCC=O(19071)',
    structure = SMILES('CC[C]=C(O)CO[CH]CCC=O'),
    E0 = (-149.106,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,350,440,435,1725,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,2782.5,750,1395,475,1775,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    molecularWeight = (170.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.43402,0.14879,-0.000177284,1.17611e-07,-3.18785e-11,-17707.9,48.4214], Tmin=(100,'K'), Tmax=(892.78,'K')), NASAPolynomial(coeffs=[17.5934,0.059056,-2.65135e-05,5.022e-09,-3.49878e-13,-21283.8,-45.9272], Tmin=(892.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-149.106,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(594.485,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Cds_S) + radical(CCsJOCs)"""),
)

species(
    label = 'CCC(=O)C1(O)CO[CH]CCC1[O](19029)',
    structure = SMILES('CCC(=O)C1(O)CO[CH]CCC1[O]'),
    E0 = (-371.832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.90001,0.141365,-0.000125746,5.2508e-08,-5.99592e-12,-44463.2,43.1657], Tmin=(100,'K'), Tmax=(938.832,'K')), NASAPolynomial(coeffs=[24.7829,0.0480867,-1.61265e-05,2.65414e-09,-1.73165e-13,-50748.3,-94.4312], Tmin=(938.832,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-371.832,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(640.214,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(oxepane) + radical(CC(C)OJ) + radical(CCsJOCs)"""),
)

species(
    label = 'CCC(=O)C(=O)CO[CH]CCC=O(19072)',
    structure = SMILES('CCC(=O)C(=O)CO[CH]CCC=O'),
    E0 = (-486.431,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,2782.5,750,1395,475,1775,1000,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,365,385,505,600,445,480,1700,1720,180,180,180,271.083,1476.87,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.149523,'amu*angstrom^2'), symmetry=1, barrier=(3.43782,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149523,'amu*angstrom^2'), symmetry=1, barrier=(3.43782,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149523,'amu*angstrom^2'), symmetry=1, barrier=(3.43782,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149523,'amu*angstrom^2'), symmetry=1, barrier=(3.43782,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149523,'amu*angstrom^2'), symmetry=1, barrier=(3.43782,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149523,'amu*angstrom^2'), symmetry=1, barrier=(3.43782,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149523,'amu*angstrom^2'), symmetry=1, barrier=(3.43782,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149523,'amu*angstrom^2'), symmetry=1, barrier=(3.43782,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149523,'amu*angstrom^2'), symmetry=1, barrier=(3.43782,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (185.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.88187,0.164131,-0.000233394,1.94519e-07,-6.65197e-11,-58268.5,49.4586], Tmin=(100,'K'), Tmax=(772.781,'K')), NASAPolynomial(coeffs=[13.1331,0.0714763,-3.46026e-05,6.68234e-09,-4.66398e-13,-60452.3,-21.7911], Tmin=(772.781,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-486.431,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(594.485,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-OdCsH) + radical(CCsJOCs)"""),
)

species(
    label = 'CCC(=O)C(O)=CO[CH]CCC=O(19073)',
    structure = SMILES('CCC(=O)C(O)=CO[CH]CCC=O'),
    E0 = (-522.556,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,350,440,435,1725,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,375,552.5,462.5,1710,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    molecularWeight = (185.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.48251,0.164993,-0.000192345,1.14459e-07,-2.69226e-11,-62579.7,47.595], Tmin=(100,'K'), Tmax=(1037.82,'K')), NASAPolynomial(coeffs=[27.4385,0.0458144,-2.00891e-05,3.80494e-09,-2.66781e-13,-68997.7,-102.729], Tmin=(1037.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-522.556,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(594.485,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(CCsJOC(O))"""),
)

species(
    label = 'CCC(=O)C([O])CO[CH]CCC=O(19074)',
    structure = SMILES('CCC(=O)C([O])CO[CH]CCC=O'),
    E0 = (-336.293,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,1380,1390,370,380,2900,435,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,2782.5,750,1395,475,1775,1000,180,180,180,519.675,633.347,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.159363,'amu*angstrom^2'), symmetry=1, barrier=(3.66408,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159363,'amu*angstrom^2'), symmetry=1, barrier=(3.66408,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159363,'amu*angstrom^2'), symmetry=1, barrier=(3.66408,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159363,'amu*angstrom^2'), symmetry=1, barrier=(3.66408,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159363,'amu*angstrom^2'), symmetry=1, barrier=(3.66408,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159363,'amu*angstrom^2'), symmetry=1, barrier=(3.66408,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159363,'amu*angstrom^2'), symmetry=1, barrier=(3.66408,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159363,'amu*angstrom^2'), symmetry=1, barrier=(3.66408,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159363,'amu*angstrom^2'), symmetry=1, barrier=(3.66408,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.8031,0.161053,-0.000218444,1.75165e-07,-5.81021e-11,-40212.5,54.1719], Tmin=(100,'K'), Tmax=(760.731,'K')), NASAPolynomial(coeffs=[13.5576,0.0711511,-3.35344e-05,6.42247e-09,-4.47038e-13,-42589.5,-19.5478], Tmin=(760.731,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-336.293,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(CCsJOCs) + radical(C=OCOJ)"""),
)

species(
    label = 'CCC(=O)C(O)[CH]O[CH]CCC=O(19075)',
    structure = SMILES('CCC(=O)C(O)[CH]O[CH]CCC=O'),
    E0 = (-399.57,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,375,552.5,462.5,1710,180,180,180,180,1600,1646.51,2859.1,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153656,'amu*angstrom^2'), symmetry=1, barrier=(3.53285,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153656,'amu*angstrom^2'), symmetry=1, barrier=(3.53285,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153656,'amu*angstrom^2'), symmetry=1, barrier=(3.53285,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153656,'amu*angstrom^2'), symmetry=1, barrier=(3.53285,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153656,'amu*angstrom^2'), symmetry=1, barrier=(3.53285,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153656,'amu*angstrom^2'), symmetry=1, barrier=(3.53285,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153656,'amu*angstrom^2'), symmetry=1, barrier=(3.53285,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153656,'amu*angstrom^2'), symmetry=1, barrier=(3.53285,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153656,'amu*angstrom^2'), symmetry=1, barrier=(3.53285,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153656,'amu*angstrom^2'), symmetry=1, barrier=(3.53285,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.21122,0.169213,-0.000226471,1.69229e-07,-5.14038e-11,-47807.2,54.3142], Tmin=(100,'K'), Tmax=(801.418,'K')), NASAPolynomial(coeffs=[17.8376,0.0641653,-2.98732e-05,5.70281e-09,-3.96937e-13,-51181.3,-42.5764], Tmin=(801.418,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-399.57,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(615.271,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(CCsJOCs) + radical(CCsJOCs)"""),
)

species(
    label = 'C[CH]C(=O)C(O)CO[CH]CCC=O(19076)',
    structure = SMILES('CC=C([O])C(O)CO[CH]CCC=O'),
    E0 = (-428.641,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2782.5,750,1395,475,1775,1000,180,180,180,180,1600,1763.18,2746.55,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156389,'amu*angstrom^2'), symmetry=1, barrier=(3.5957,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156389,'amu*angstrom^2'), symmetry=1, barrier=(3.5957,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156389,'amu*angstrom^2'), symmetry=1, barrier=(3.5957,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156389,'amu*angstrom^2'), symmetry=1, barrier=(3.5957,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156389,'amu*angstrom^2'), symmetry=1, barrier=(3.5957,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156389,'amu*angstrom^2'), symmetry=1, barrier=(3.5957,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156389,'amu*angstrom^2'), symmetry=1, barrier=(3.5957,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156389,'amu*angstrom^2'), symmetry=1, barrier=(3.5957,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156389,'amu*angstrom^2'), symmetry=1, barrier=(3.5957,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.87401,0.159555,-0.000196005,1.33526e-07,-3.70565e-11,-51313.4,53.2452], Tmin=(100,'K'), Tmax=(873.719,'K')), NASAPolynomial(coeffs=[18.5031,0.0616859,-2.798e-05,5.31696e-09,-3.70748e-13,-55048.9,-47.0007], Tmin=(873.719,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-428.641,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCsJOCs) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]CC(=O)C(O)CO[CH]CCC=O(19077)',
    structure = SMILES('[CH2]CC(=O)C(O)CO[CH]CCC=O'),
    E0 = (-368.437,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,3615,1277.5,1000,375,552.5,462.5,1710,1380,1390,370,380,2900,435,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,2782.5,750,1395,475,1775,1000,180,180,180,180,1600,1605.19,2903.23,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152753,'amu*angstrom^2'), symmetry=1, barrier=(3.51209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152753,'amu*angstrom^2'), symmetry=1, barrier=(3.51209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152753,'amu*angstrom^2'), symmetry=1, barrier=(3.51209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152753,'amu*angstrom^2'), symmetry=1, barrier=(3.51209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152753,'amu*angstrom^2'), symmetry=1, barrier=(3.51209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152753,'amu*angstrom^2'), symmetry=1, barrier=(3.51209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152753,'amu*angstrom^2'), symmetry=1, barrier=(3.51209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152753,'amu*angstrom^2'), symmetry=1, barrier=(3.51209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152753,'amu*angstrom^2'), symmetry=1, barrier=(3.51209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152753,'amu*angstrom^2'), symmetry=1, barrier=(3.51209,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.75929,0.162948,-0.000206361,1.37429e-07,-3.21038e-11,-44082.5,53.9291], Tmin=(100,'K'), Tmax=(639.977,'K')), NASAPolynomial(coeffs=[15.7006,0.0681128,-3.22332e-05,6.1971e-09,-4.33041e-13,-46866,-30.1768], Tmin=(639.977,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-368.437,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(615.271,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(CJCC=O) + radical(CCsJOCs)"""),
)

species(
    label = 'CCC(=O)C(O)CO[CH][CH]CC=O(19078)',
    structure = SMILES('CCC(=O)C(O)CO[CH][CH]CC=O'),
    E0 = (-380.124,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,375,552.5,462.5,1710,180,180,180,180,1600,1633.35,2871.24,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153322,'amu*angstrom^2'), symmetry=1, barrier=(3.52517,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153322,'amu*angstrom^2'), symmetry=1, barrier=(3.52517,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153322,'amu*angstrom^2'), symmetry=1, barrier=(3.52517,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153322,'amu*angstrom^2'), symmetry=1, barrier=(3.52517,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153322,'amu*angstrom^2'), symmetry=1, barrier=(3.52517,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153322,'amu*angstrom^2'), symmetry=1, barrier=(3.52517,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153322,'amu*angstrom^2'), symmetry=1, barrier=(3.52517,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153322,'amu*angstrom^2'), symmetry=1, barrier=(3.52517,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153322,'amu*angstrom^2'), symmetry=1, barrier=(3.52517,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153322,'amu*angstrom^2'), symmetry=1, barrier=(3.52517,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.58766,0.155356,-0.000192417,1.37066e-07,-4.04622e-11,-45490.5,55.1349], Tmin=(100,'K'), Tmax=(817.706,'K')), NASAPolynomial(coeffs=[15.3086,0.0678112,-3.18218e-05,6.13204e-09,-4.30497e-13,-48417.2,-27.602], Tmin=(817.706,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-380.124,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(615.271,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(CCJCO) + radical(CCsJOCs)"""),
)

species(
    label = 'CCC(=O)C(O)CO[CH]C[CH]C=O(19079)',
    structure = SMILES('CCC(=O)C(O)CO[CH]CC=C[O]'),
    E0 = (-436.774,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,375,552.5,462.5,1710,180,180,180,435.275,699.411,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155787,'amu*angstrom^2'), symmetry=1, barrier=(3.58185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155787,'amu*angstrom^2'), symmetry=1, barrier=(3.58185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155787,'amu*angstrom^2'), symmetry=1, barrier=(3.58185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155787,'amu*angstrom^2'), symmetry=1, barrier=(3.58185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155787,'amu*angstrom^2'), symmetry=1, barrier=(3.58185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155787,'amu*angstrom^2'), symmetry=1, barrier=(3.58185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155787,'amu*angstrom^2'), symmetry=1, barrier=(3.58185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155787,'amu*angstrom^2'), symmetry=1, barrier=(3.58185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155787,'amu*angstrom^2'), symmetry=1, barrier=(3.58185,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.02944,0.152514,-0.000158822,8.61157e-08,-1.86963e-11,-52276.8,53.9134], Tmin=(100,'K'), Tmax=(1113.87,'K')), NASAPolynomial(coeffs=[25.3101,0.0507425,-2.17696e-05,4.08653e-09,-2.85165e-13,-58590,-85.8645], Tmin=(1113.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-436.774,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CCsJOCs)"""),
)

species(
    label = 'CCC(=O)C(O)CO[CH]CC[C]=O(19080)',
    structure = SMILES('CCC(=O)C(O)CO[CH]CC[C]=O'),
    E0 = (-420.065,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,375,552.5,462.5,1710,180,180,180,180,1600,1717.35,2795.21,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155348,'amu*angstrom^2'), symmetry=1, barrier=(3.57175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155348,'amu*angstrom^2'), symmetry=1, barrier=(3.57175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155348,'amu*angstrom^2'), symmetry=1, barrier=(3.57175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155348,'amu*angstrom^2'), symmetry=1, barrier=(3.57175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155348,'amu*angstrom^2'), symmetry=1, barrier=(3.57175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155348,'amu*angstrom^2'), symmetry=1, barrier=(3.57175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155348,'amu*angstrom^2'), symmetry=1, barrier=(3.57175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155348,'amu*angstrom^2'), symmetry=1, barrier=(3.57175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155348,'amu*angstrom^2'), symmetry=1, barrier=(3.57175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155348,'amu*angstrom^2'), symmetry=1, barrier=(3.57175,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.70612,0.161113,-0.00020349,1.35901e-07,-3.22835e-11,-50293.2,53.9187], Tmin=(100,'K'), Tmax=(647.391,'K')), NASAPolynomial(coeffs=[15.7501,0.0669406,-3.13152e-05,5.99061e-09,-4.17467e-13,-53099.1,-30.3116], Tmin=(647.391,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-420.065,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(615.271,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(CCsJOCs) + radical(CCCJ=O)"""),
)

species(
    label = 'CC[C]1OC(CCC=O)OC[C]1O(19081)',
    structure = SMILES('CC[C]1OC(CCC=O)OC[C]1O'),
    E0 = (-413.895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.96536,0.134554,-0.000130494,6.65186e-08,-1.02542e-11,-49567.9,43.7578], Tmin=(100,'K'), Tmax=(753.144,'K')), NASAPolynomial(coeffs=[15.1183,0.062234,-2.31303e-05,3.94299e-09,-2.57809e-13,-52663.4,-37.2845], Tmin=(753.144,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-413.895,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(631.9,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(1,3-Dioxane) + radical(C2CsJOH) + radical(C2CsJOCs)"""),
)

species(
    label = 'CCC(=O)C1(O)CO[CH]CC[CH]O1(19082)',
    structure = SMILES('CCC(=O)C1(O)CO[CH]CC[CH]O1'),
    E0 = (-403.181,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.99264,0.130304,-7.49623e-05,-8.97672e-09,1.6096e-11,-48218.9,44.5079], Tmin=(100,'K'), Tmax=(987.213,'K')), NASAPolynomial(coeffs=[32.3022,0.0387106,-1.39145e-05,2.55309e-09,-1.83533e-13,-57693,-138.004], Tmin=(987.213,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-403.181,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(640.214,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclooctane) + radical(CCsJOCs) + radical(CCsJOCs)"""),
)

species(
    label = 'CCC(=O)C(O)=COCCCC=O(19083)',
    structure = SMILES('CCC(=O)C(O)=COCCCC=O'),
    E0 = (-716.484,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.08722,0.155304,-0.000164578,9.07755e-08,-2.00646e-11,-85917.3,46.6997], Tmin=(100,'K'), Tmax=(1093.98,'K')), NASAPolynomial(coeffs=[25.2204,0.0518006,-2.26602e-05,4.29099e-09,-3.00796e-13,-92110.8,-92.4107], Tmin=(1093.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-716.484,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsOsH) + group(Cds-OdCsH)"""),
)

species(
    label = 'CCC(=O)C(O)COC=CCC=O(19084)',
    structure = SMILES('CCC(=O)C(O)COC=CCC=O'),
    E0 = (-677.544,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.15953,0.138868,-0.000118499,4.95029e-08,-8.15648e-12,-81216.6,52.7258], Tmin=(100,'K'), Tmax=(1458.13,'K')), NASAPolynomial(coeffs=[32.6976,0.040503,-1.7309e-05,3.23802e-09,-2.2422e-13,-91673.4,-133.788], Tmin=(1458.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-677.544,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH)"""),
)

species(
    label = 'CCC(=O)C(=O)COCCCC=O(19085)',
    structure = SMILES('CCC(=O)C(=O)COCCCC=O'),
    E0 = (-666.887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.13433,0.1342,-9.43958e-05,-6.97916e-08,1.17661e-10,-80044.7,44.0528], Tmin=(100,'K'), Tmax=(488.025,'K')), NASAPolynomial(coeffs=[9.9916,0.0795437,-3.86996e-05,7.52574e-09,-5.28804e-13,-81565.7,-6.09928], Tmin=(488.025,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-666.887,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-OdCsH)"""),
)

species(
    label = 'CCC(O)([C]=O)CO[CH]CCC=O(19086)',
    structure = SMILES('CCC(O)([C]=O)CO[CH]CCC=O'),
    E0 = (-401.094,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,1855,455,950,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,2782.5,750,1395,475,1775,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
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
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.65355,0.180683,-0.000265229,2.16186e-07,-7.05916e-11,-47976.7,53.069], Tmin=(100,'K'), Tmax=(812.737,'K')), NASAPolynomial(coeffs=[17.4555,0.0650783,-3.02468e-05,5.70201e-09,-3.91157e-13,-51021,-42.0129], Tmin=(812.737,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-401.094,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(615.271,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CC(C)(O)CJ=O) + radical(CCsJOCs)"""),
)

species(
    label = 'CCC(=O)C1(O)COC1CCC=O(18462)',
    structure = SMILES('CCC(=O)C1(O)COC1CCC=O'),
    E0 = (-629.979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.50849,0.143842,-0.000155394,9.57147e-08,-2.40355e-11,-75535.2,45.5637], Tmin=(100,'K'), Tmax=(965.076,'K')), NASAPolynomial(coeffs=[18.1152,0.058363,-2.2537e-05,3.93962e-09,-2.61772e-13,-79515.9,-53.2011], Tmin=(965.076,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-629.979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(627.743,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + ring(Oxetane)"""),
)

species(
    label = 'C=C(OC)[C](O)CO[CH]CCC=O(19087)',
    structure = SMILES('[CH2]C(OC)=C(O)CO[CH]CCC=O'),
    E0 = (-401.26,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2782.5,750,1395,475,1775,1000,180,180,180,180,1600,1608,2885.57,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152621,'amu*angstrom^2'), symmetry=1, barrier=(3.50905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152621,'amu*angstrom^2'), symmetry=1, barrier=(3.50905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152621,'amu*angstrom^2'), symmetry=1, barrier=(3.50905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152621,'amu*angstrom^2'), symmetry=1, barrier=(3.50905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152621,'amu*angstrom^2'), symmetry=1, barrier=(3.50905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152621,'amu*angstrom^2'), symmetry=1, barrier=(3.50905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152621,'amu*angstrom^2'), symmetry=1, barrier=(3.50905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152621,'amu*angstrom^2'), symmetry=1, barrier=(3.50905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152621,'amu*angstrom^2'), symmetry=1, barrier=(3.50905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152621,'amu*angstrom^2'), symmetry=1, barrier=(3.50905,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (186.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.49577,0.165736,-0.000194731,1.19321e-07,-2.89973e-11,-47991.1,51.2156], Tmin=(100,'K'), Tmax=(1005.13,'K')), NASAPolynomial(coeffs=[25.8319,0.0490215,-2.05486e-05,3.78913e-09,-2.61138e-13,-53886.6,-90.423], Tmin=(1005.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-401.26,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(615.271,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(CCsJOCs) + radical(C=C(O)CJ)"""),
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
    E0 = (-463.701,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-357.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-290.379,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-342.691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-449.439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-388.559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-202.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-335.471,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-336.089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-336.089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-244.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-367.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-313.207,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-349.98,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-349.98,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-436.819,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-401.578,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-295.289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-383.682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-308.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-165.871,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-59.185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-232.485,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-405.544,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-370.286,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-99.1001,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-438.728,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-21.182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-118.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-285.571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-456.17,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (82.4996,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (-320.082,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (-108.305,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (-318.046,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (-50.1786,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (93.8991,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (-371.832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (-248.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (-303.995,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (-182.384,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (-260.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (-260.082,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (-258.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (-272.516,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (-375.981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (-338.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (-413.895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (-403.181,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (-400.301,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (-455.333,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (-450.626,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (-243.776,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (-455.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (-87.4598,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    products = ['C=C(O)C(=O)CC(4626)', 'O=CCCC=O(5767)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    products = ['CCC([O])=C(O)COC1CCC1[O](19042)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.13771e+06,'s^-1'), n=1.58803, Ea=(105.768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_csHNd] for rate rule [R5_SS_CO;carbonylbond_intra_H;radadd_intra_csHNd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 101.2 to 105.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    products = ['CC[C]([O])C1(O)COC1CCC=O(19043)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.73e+06,'s^-1'), n=1.31, Ea=(173.322,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_csHNd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 171.4 to 173.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'CCC([O])=C(O)COC=CCC=O(19044)'],
    products = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.72e+08,'cm^3/(mol*s)'), n=1.477, Ea=(6.73624,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2825 used for Cds-CsH_Cds-OsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-OsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C(O)=C([O])CC(4557)', 'O=CCCC=O(5767)'],
    products = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.000397355,'m^3/(mol*s)'), n=2.77646, Ea=(45.4073,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_R;CsJ-CdHH] for rate rule [Od_CO-CsH;CsJ-CdHH]
Euclidian distance = 4.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2CHO(40)', 'C=COCC(O)=C([O])CC(6482)'],
    products = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.0114756,'m^3/(mol*s)'), n=2.44484, Ea=(36.1431,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-HH_Cds;CsJ-OneDeHH] for rate rule [Cds-HH_Cds-OsH;CsJ-COHH]
Euclidian distance = 1.41421356237
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C2H5(29)', 'O=C=C(O)CO[CH]CCC=O(19045)'],
    products = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(8.04,'m^3/(mol*s)'), n=1.68, Ea=(54.1828,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cdd_Od;CsJ] for rate rule [Ck_O;CsJ-CsHH]
Euclidian distance = 1.41421356237
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CCC([O])=C(O)COC[CH]CC=O(19046)'],
    products = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(5.4e-20,'s^-1'), n=9.13, Ea=(108.784,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 341 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    products = ['CCC([O])=C(O)[CH]OCCCC=O(19047)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(6.82e+09,'s^-1'), n=0.73, Ea=(127.612,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_H/NonDeC;Cs_H_out_H/Cd] for rate rule [R3H_SS_O;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    products = ['CCC([O])=C(O)COCC[CH]C=O(19048)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6.82e+09,'s^-1'), n=0.73, Ea=(127.612,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_1H;Cs_H_out_H/OneDe] for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeO;Cs_H_out_H/CO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C[CH]C(O)=C(O)CO[CH]CCC=O(19049)'],
    products = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(187863,'s^-1'), n=2.03383, Ea=(157.381,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;C_rad_out_H/NonDeC;O_H_out] + [R3H_SS_2Cd;C_rad_out_H/NonDeC;XH_out] for rate rule [R3H_SS_2Cd;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    products = ['CCC([O])=C(O)COCCC[C]=O(19050)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(27900,'s^-1'), n=1.97, Ea=(96.6504,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/NonDeO;XH_out] for rate rule [R4H_SSS;C_rad_out_H/NonDeO;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]CC(O)=C(O)CO[CH]CCC=O(19051)'],
    products = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.00963743,'s^-1'), n=3.795, Ea=(83.0524,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;C_rad_out_2H;O_H_out] + [R4H_SS(Cd)S;C_rad_out_2H;XH_out] for rate rule [R4H_SS(Cd)S;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['CCC(O)=C([O])CO[CH]CCC=O(19052)'],
    products = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.11e+06,'s^-1','*|/',3), n=1.78, Ea=(113.721,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4H_SDS;O_rad_out;XH_out] for rate rule [R4H_SDS;O_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    products = ['CCC(O)=C(O)[CH]O[CH]CCC=O(19053)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.22e+06,'s^-1','*|/',3), n=1.78, Ea=(113.721,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4H_SDS;O_rad_out;Cs_H_out_1H] for rate rule [R4H_SDS;O_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    products = ['CCC([O])=C([O])COCCCC=O(19054)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.0378492,'s^-1'), n=3.26, Ea=(26.8822,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;C_rad_out_H/NonDeC;XH_out] for rate rule [R5H_SSSS;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C[CH]C([O])=C(O)COCCCC=O(19055)'],
    products = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(25800,'s^-1'), n=1.67, Ea=(42.6768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H;C_rad_out_single;Cs_H_out_H/NonDeC] for rate rule [R6H_SMSSR;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CCC(O)=C(O)CO[CH][CH]CC=O(19056)'],
    products = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(423.689,'s^-1'), n=2.58367, Ea=(106.315,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;C_rad_out_H/NonDeC;O_H_out] for rate rule [R7HJ_1;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]CC([O])=C(O)COCCCC=O(19057)'],
    products = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1062,'s^-1'), n=1.81, Ea=(55.2288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7H;C_rad_out_2H;Cs_H_out_H/NonDeC] for rate rule [R7H;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['CCC(O)=C(O)CO[CH]C[CH]C=O(19058)'],
    products = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(88626.7,'s^-1'), n=2.08905, Ea=(149.401,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;C_rad_out_1H;O_H_out] + [RnH;C_rad_out_H/OneDe;XH_out] for rate rule [R8Hall;C_rad_out_H/OneDe;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(O)=C([O])CC(4557)', '[O][CH]CCC=O(5764)'],
    products = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.56662e+07,'m^3/(mol*s)'), n=0.074875, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [C_rad/H2/Cd;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['CCC([O])=[C]O(4619)', '[CH2]O[CH]CCC=O(18485)'],
    products = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(6.41031e+09,'m^3/(mol*s)'), n=-0.9855, Ea=(1.09644,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_rad;C_pri_rad] for rate rule [Cd_rad;C_rad/H2/O]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    products = ['CC[C]1OC1(O)CO[CH]CCC=O(19059)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd_NdNd;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    products = ['CCC([O])=C(O)COC1CC[CH]O1(19060)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.64e+06,'s^-1'), n=1.15, Ea=(58.1576,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_csHNd] for rate rule [R5_SS_CO;carbonyl_intra_H;radadd_intra_csHO]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    products = ['CCC1([O])[C](O)COC1CCC=O(19061)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.64e+06,'s^-1'), n=1.15, Ea=(93.415,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra;radadd_intra_csHCs] for rate rule [R5_SS_D;doublebond_intra_secNd;radadd_intra_csHCs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 89.9 to 93.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction26',
    reactants = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    products = ['CCC1OO[CH]CC[CH]OCC=1O(19062)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(5.649e+12,'s^-1'), n=0.287, Ea=(364.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;carbonyl_intra_H;radadd_intra] for rate rule [R10_linear;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    products = ['CCC(O)=C(O)COC=CCC=O(19063)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CH2(S)(23)', 'CC([O])=C(O)CO[CH]CCC=O(19064)'],
    products = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.87e+13,'cm^3/(mol*s)','*|/',0.25), n=-0.146, Ea=(0.0118826,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 6 used for carbene;C_pri/Cd
Exact match found for rate rule [carbene;C_pri/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene"""),
)

reaction(
    label = 'reaction29',
    reactants = ['CO(12)', 'CC[CH]OC[C](O)C(=O)CC(11983)'],
    products = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(274200,'cm^3/(mol*s)'), n=2.53, Ea=(357.732,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 4 used for CO;C_pri/NonDeC
Exact match found for rate rule [CO;C_pri/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(CC=O)OCC(O)=C([O])CC(19065)'],
    products = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction31',
    reactants = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    products = ['CCC1OC(CCC=O)OCC=1O(19066)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;C_rad_out_1H;Ypri_rad_out] for rate rule [R6_SSSDS;C_rad_out_H/NonDeC;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['CCC(=[C]CO[CH]CCC=O)OO(19067)'],
    products = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(3.01978e+11,'s^-1'), n=0, Ea=(107.111,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2OOH;Y_rad_out] for rate rule [R2OOH_D;Cd_rad_out_ND]
Euclidian distance = 2.2360679775
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    products = ['CCC([O])C(=O)CO[CH]CCC=O(19068)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1713.16,'s^-1'), n=2.9, Ea=(143.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_Cs;R_O_H]
Euclidian distance = 0
family: ketoenol
Ea raised from 142.3 to 143.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=COC[CH]OCC(O)=C([O])CC(19069)'],
    products = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 7 used for R_ROR;R1_doublebond_CH2;R2_doublebond_H;R_O_C
Exact match found for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction35',
    reactants = ['CCC([O])=C(O)CO[CH]CC=CO(19070)'],
    products = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_H;R_O_H]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction36',
    reactants = ['CCC([O])=C(O)C[O](4968)', '[CH]CCC=O(18484)'],
    products = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction37',
    reactants = ['O(4)', 'CC[C]=C(O)CO[CH]CCC=O(19071)'],
    products = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction38',
    reactants = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    products = ['CCC(=O)C1(O)CO[CH]CCC1[O](19029)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.19e+11,'s^-1'), n=0.08, Ea=(91.869,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R8;multiplebond_intra;radadd_intra_cs] for rate rule [R8;carbonylbond_intra_H;radadd_intra_csNdDe]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 89.1 to 91.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction39',
    reactants = ['H(3)', 'CCC(=O)C(=O)CO[CH]CCC=O(19072)'],
    products = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(10.8137,'m^3/(mol*s)'), n=2.04, Ea=(26.5684,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Od_CO-DeNd;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction40',
    reactants = ['H(3)', 'CCC(=O)C(O)=CO[CH]CCC=O(19073)'],
    products = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(2.182e+10,'cm^3/(mol*s)'), n=0.859, Ea=(6.76971,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [Cds-OsH_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction41',
    reactants = ['CCC(=O)C([O])CO[CH]CCC=O(19074)'],
    products = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.70223e+09,'s^-1'), n=1.15155, Ea=(153.908,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_OneDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_CO]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['CCC(=O)C(O)[CH]O[CH]CCC=O(19075)'],
    products = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(8.1111e+08,'s^-1'), n=1.20078, Ea=(139.26,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_1H;Cs_H_out_OneDe] for rate rule [R2H_S;C_rad_out_H/NonDeO;Cs_H_out_CO]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['C[CH]C(=O)C(O)CO[CH]CCC=O(19076)'],
    products = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(7.66994e+09,'s^-1'), n=0.768667, Ea=(168.559,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_H/NonDeC;Cs_H_out_NonDe] for rate rule [R3H_SS;C_rad_out_H/NonDeC;Cs_H_out_NDMustO]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2]CC(=O)C(O)CO[CH]CCC=O(19077)'],
    products = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.86e+10,'s^-1'), n=0.58, Ea=(109.579,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_2H;Cs_H_out_NonDe] for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_NDMustO]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['CCC(=O)C(O)CO[CH][CH]CC=O(19078)'],
    products = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(0.106062,'s^-1'), n=4.02875, Ea=(107.607,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;C_rad_out_H/NonDeC;Cs_H_out_OneDe] for rate rule [R5HJ_1;C_rad_out_H/NonDeC;Cs_H_out_CO]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['CCC(=O)C(O)CO[CH]C[CH]C=O(19079)'],
    products = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(46.1,'s^-1'), n=3.21, Ea=(60.7935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;C_rad_out_1H;Cs_H_out_OneDe] for rate rule [R6HJ_2;C_rad_out_H/OneDe;Cs_H_out_CO]
Euclidian distance = 1.73205080757
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['CCC(=O)C(O)CO[CH]CC[C]=O(19080)'],
    products = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(20912.2,'s^-1'), n=1.956, Ea=(81.9966,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Y_rad_out;Cs_H_out_OneDe] + [R7Hall;Y_rad_out;Cs_H_out_noH] for rate rule [R7HJ_3;CO_rad_out;Cs_H_out_CO]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    products = ['CC[C]1OC(CCC=O)OC[C]1O(19081)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(487000,'s^-1'), n=1.17, Ea=(49.8062,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_csHCs] for rate rule [R6_linear;carbonyl_intra_Nd;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 46.3 to 49.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction49',
    reactants = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    products = ['CCC(=O)C1(O)CO[CH]CC[CH]O1(19082)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(5.649e+12,'s^-1'), n=0.287, Ea=(60.5199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;carbonyl_intra_H;radadd_intra] for rate rule [R8_linear;carbonyl_intra_H;radadd_intra_csNdDe]
Euclidian distance = 2.82842712475
family: Intra_R_Add_Endocyclic
Ea raised from 56.3 to 60.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction50',
    reactants = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    products = ['CCC(=O)C(O)=COCCCC=O(19083)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 1 used for R3radExo;Y_rad_NDe;XH_Rrad_NDe
Exact match found for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction51',
    reactants = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    products = ['CCC(=O)C(O)COC=CCC=O(19084)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(6.42e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction52',
    reactants = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    products = ['CCC(=O)C(=O)COCCCC=O(19085)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(13.075,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5;Y_rad;XH_Rrad_De] + [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad_De]
Euclidian distance = 1.41421356237
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction53',
    reactants = ['CCC(O)([C]=O)CO[CH]CCC=O(19086)'],
    products = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction54',
    reactants = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    products = ['CCC(=O)C1(O)COC1CCC=O(18462)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_OneDe/O]
Euclidian distance = 3.60555127546
family: Birad_recombination"""),
)

reaction(
    label = 'reaction55',
    reactants = ['C=C(OC)[C](O)CO[CH]CCC=O(19087)'],
    products = ['CCC(=O)[C](O)CO[CH]CCC=O(18458)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

network(
    label = '3511',
    isomers = [
        'CCC(=O)[C](O)CO[CH]CCC=O(18458)',
    ],
    reactants = [
        ('C=C(O)C(=O)CC(4626)', 'O=CCCC=O(5767)'),
        ('[CH2]C(O)=C([O])CC(4557)', 'O=CCCC=O(5767)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '3511',
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

