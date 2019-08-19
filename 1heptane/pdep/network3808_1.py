species(
    label = '[O]C(CCC=O)O[CH]CCC=O(21393)',
    structure = SMILES('[O]C(CCC=O)O[CH]CCC=O'),
    E0 = (-328.163,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,180,180,180,185.102,1592.23,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.151948,'amu*angstrom^2'), symmetry=1, barrier=(3.49359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151948,'amu*angstrom^2'), symmetry=1, barrier=(3.49359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151948,'amu*angstrom^2'), symmetry=1, barrier=(3.49359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151948,'amu*angstrom^2'), symmetry=1, barrier=(3.49359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151948,'amu*angstrom^2'), symmetry=1, barrier=(3.49359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151948,'amu*angstrom^2'), symmetry=1, barrier=(3.49359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151948,'amu*angstrom^2'), symmetry=1, barrier=(3.49359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151948,'amu*angstrom^2'), symmetry=1, barrier=(3.49359,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.178,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.10328,0.15053,-0.00023277,2.12126e-07,-7.69714e-11,-39264.8,49.3964], Tmin=(100,'K'), Tmax=(815.66,'K')), NASAPolynomial(coeffs=[7.66271,0.0731043,-3.60703e-05,6.96523e-09,-4.8353e-13,-39875.5,10.2933], Tmin=(815.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-328.163,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CCsJOCs)"""),
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
    label = '[O]C(CCC=O)OC1CCC1[O](21404)',
    structure = SMILES('[O]C(CCC=O)OC1CCC1[O]'),
    E0 = (-222.395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.178,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.465253,0.104502,-7.59287e-05,2.84482e-08,-4.54748e-12,-26593.3,42.7263], Tmin=(100,'K'), Tmax=(1379.57,'K')), NASAPolynomial(coeffs=[14.615,0.0607773,-2.83862e-05,5.47347e-09,-3.84047e-13,-30754.2,-34.8797], Tmin=(1379.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-222.395,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(561.227,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + ring(Cyclobutane) + radical(CCOJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]C1CCC(O[CH]CCC=O)O1(21405)',
    structure = SMILES('[O]C1CCC(O[CH]CCC=O)O1'),
    E0 = (-340.025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.178,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.04573,0.111155,-9.53622e-05,4.68138e-08,-9.73549e-12,-40714,41.2519], Tmin=(100,'K'), Tmax=(1125.82,'K')), NASAPolynomial(coeffs=[14.1232,0.057261,-2.35557e-05,4.29303e-09,-2.93366e-13,-44129.5,-33.7274], Tmin=(1125.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-340.025,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(561.227,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(Tetrahydrofuran) + radical(CCOJ) + radical(CCsJOCs)"""),
)

species(
    label = '[O]C1CCC([O])C(CCC=O)O1(21406)',
    structure = SMILES('[O]C1CCC([O])C(CCC=O)O1'),
    E0 = (-327.581,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.178,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.11013,0.103426,-7.33416e-05,2.81944e-08,-4.5267e-12,-39207.8,40.3798], Tmin=(100,'K'), Tmax=(1436.69,'K')), NASAPolynomial(coeffs=[17.1822,0.0524971,-2.01693e-05,3.52116e-09,-2.3333e-13,-44463.9,-54.4985], Tmin=(1436.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-327.581,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(569.541,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + ring(Oxane) + radical(CCOJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]C1CC[CH]OC(CCC=O)O1(21407)',
    structure = SMILES('[O]C1CC[CH]OC(CCC=O)O1'),
    E0 = (-334.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.178,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.55029,0.0822138,5.43596e-05,-1.57333e-07,7.41788e-11,-40014.9,45.0903], Tmin=(100,'K'), Tmax=(929.009,'K')), NASAPolynomial(coeffs=[38.2649,0.0131121,7.09484e-07,-2.67332e-10,5.17687e-15,-51828.4,-167.83], Tmin=(929.009,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-334.67,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(569.541,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(Cycloheptane) + radical(CCOJ) + radical(CCsJOCs)"""),
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
    label = '[O]C(CCC=O)OC=CCC=O(21408)',
    structure = SMILES('[O]C(CCC=O)OC=CCC=O'),
    E0 = (-425.681,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    molecularWeight = (171.171,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.08989,0.112362,-9.18441e-05,3.72352e-08,-6.17278e-12,-51015.4,43.7955], Tmin=(100,'K'), Tmax=(1393.97,'K')), NASAPolynomial(coeffs=[20.9696,0.0490619,-2.37288e-05,4.65891e-09,-3.30396e-13,-57165.5,-69.9559], Tmin=(1393.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-425.681,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CCOJ)"""),
)

species(
    label = 'O=CCC[CH]OC(=O)CCC=O(21409)',
    structure = SMILES('O=CCC[CH]OC(=O)CCC=O'),
    E0 = (-531.917,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,3025,407.5,1350,352.5,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
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
    molecularWeight = (171.171,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.71171,0.138976,-0.000203457,1.8052e-07,-6.51246e-11,-63782,45.938], Tmin=(100,'K'), Tmax=(796.675,'K')), NASAPolynomial(coeffs=[8.09438,0.0693518,-3.39762e-05,6.57205e-09,-4.58081e-13,-64697.4,4.91927], Tmin=(796.675,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-531.917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CCsJOC(O))"""),
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
    label = '[O]C(CCC=O)OC[CH]CC=O(21410)',
    structure = SMILES('[O]C(CCC=O)OC[CH]CC=O'),
    E0 = (-308.716,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,180,180,180,292.482,1491.34,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.149707,'amu*angstrom^2'), symmetry=1, barrier=(3.44205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149707,'amu*angstrom^2'), symmetry=1, barrier=(3.44205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149707,'amu*angstrom^2'), symmetry=1, barrier=(3.44205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149707,'amu*angstrom^2'), symmetry=1, barrier=(3.44205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149707,'amu*angstrom^2'), symmetry=1, barrier=(3.44205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149707,'amu*angstrom^2'), symmetry=1, barrier=(3.44205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149707,'amu*angstrom^2'), symmetry=1, barrier=(3.44205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149707,'amu*angstrom^2'), symmetry=1, barrier=(3.44205,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.178,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.58766,0.138091,-0.000204443,1.88717e-07,-7.0481e-11,-36943.6,50.5954], Tmin=(100,'K'), Tmax=(800.013,'K')), NASAPolynomial(coeffs=[5.29428,0.076459,-3.78434e-05,7.35166e-09,-5.13455e-13,-37173.5,24.3742], Tmin=(800.013,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-308.716,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CCJCO)"""),
)

species(
    label = 'O=CCC[CH]O[C](O)CCC=O(21411)',
    structure = SMILES('O=CCC[CH]O[C](O)CCC=O'),
    E0 = (-348.621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.178,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.38744,0.153476,-0.000228774,1.96014e-07,-6.76454e-11,-41711.9,50.247], Tmin=(100,'K'), Tmax=(803.917,'K')), NASAPolynomial(coeffs=[12.1013,0.0644779,-3.11685e-05,5.98383e-09,-4.14856e-13,-43495.1,-13.0923], Tmin=(803.917,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-348.621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(544.598,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CCsJOCs) + radical(Cs_P)"""),
)

species(
    label = '[O][C](CCC=O)OCCCC=O(21412)',
    structure = SMILES('[O][C](CCC=O)OCCCC=O'),
    E0 = (-303.372,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,360,370,350,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,180,180,180,264.539,1524.9,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.150602,'amu*angstrom^2'), symmetry=1, barrier=(3.46263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150602,'amu*angstrom^2'), symmetry=1, barrier=(3.46263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150602,'amu*angstrom^2'), symmetry=1, barrier=(3.46263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150602,'amu*angstrom^2'), symmetry=1, barrier=(3.46263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150602,'amu*angstrom^2'), symmetry=1, barrier=(3.46263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150602,'amu*angstrom^2'), symmetry=1, barrier=(3.46263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150602,'amu*angstrom^2'), symmetry=1, barrier=(3.46263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150602,'amu*angstrom^2'), symmetry=1, barrier=(3.46263,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.178,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.75977,0.144224,-0.000224205,2.11342e-07,-7.91028e-11,-36296.9,49.4036], Tmin=(100,'K'), Tmax=(814.269,'K')), NASAPolynomial(coeffs=[4.6248,0.0782795,-3.90223e-05,7.57125e-09,-5.2709e-13,-36190.2,26.9532], Tmin=(814.269,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-303.372,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(Cs_P) + radical(CCOJ)"""),
)

species(
    label = '[O]C(CCC=O)OCC[CH]C=O(21413)',
    structure = SMILES('[O]C=CCCOC([O])CCC=O'),
    E0 = (-365.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.178,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.27998,0.125491,-0.000132222,8.07904e-08,-2.10772e-11,-43761.3,46.745], Tmin=(100,'K'), Tmax=(905.545,'K')), NASAPolynomial(coeffs=[12.6047,0.0641569,-3.0622e-05,5.98939e-09,-4.25719e-13,-46275.9,-18.8627], Tmin=(905.545,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-365.367,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CCOJ)"""),
)

species(
    label = 'O=CC[CH]C(O)O[CH]CCC=O(21414)',
    structure = SMILES('O=CC[CH]C(O)O[CH]CCC=O'),
    E0 = (-353.966,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,3000,3050,390,425,1340,1360,335,370,3615,1277.5,1000,1380,1390,370,380,2900,435,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    molecularWeight = (172.178,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.17759,0.146845,-0.000207001,1.70344e-07,-5.75052e-11,-42360.2,51.3067], Tmin=(100,'K'), Tmax=(772.127,'K')), NASAPolynomial(coeffs=[12.6845,0.0628149,-3.0085e-05,5.7876e-09,-4.03211e-13,-44445.5,-15.1923], Tmin=(772.127,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-353.966,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(544.598,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CCJCO) + radical(CCsJOCs)"""),
)

species(
    label = 'O=C[CH]CC(O)O[CH]CCC=O(21415)',
    structure = SMILES('[O]C=CCC(O)O[CH]CCC=O'),
    E0 = (-410.616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.178,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.19024,0.138472,-0.000151739,8.76633e-08,-2.03989e-11,-49164.6,48.5761], Tmin=(100,'K'), Tmax=(1038.71,'K')), NASAPolynomial(coeffs=[21.2206,0.0483195,-2.15519e-05,4.10714e-09,-2.88553e-13,-54028.1,-65.2571], Tmin=(1038.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-410.616,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(CCsJOCs) + radical(C=COJ)"""),
)

species(
    label = '[O]C(CCC=O)OCCC[C]=O(21416)',
    structure = SMILES('[O]C(CCC=O)OCCC[C]=O'),
    E0 = (-348.658,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.178,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.8752,0.146491,-0.000228736,2.12867e-07,-7.84271e-11,-41739,49.95], Tmin=(100,'K'), Tmax=(820.721,'K')), NASAPolynomial(coeffs=[5.75398,0.0755238,-3.72841e-05,7.19496e-09,-4.98974e-13,-41853.5,21.5829], Tmin=(820.721,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-348.658,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(CCOJ)"""),
)

species(
    label = '[O]C([CH]CC=O)OCCCC=O(21417)',
    structure = SMILES('[O]C([CH]CC=O)OCCCC=O'),
    E0 = (-308.716,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,180,180,180,292.482,1491.34,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.149707,'amu*angstrom^2'), symmetry=1, barrier=(3.44205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149707,'amu*angstrom^2'), symmetry=1, barrier=(3.44205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149707,'amu*angstrom^2'), symmetry=1, barrier=(3.44205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149707,'amu*angstrom^2'), symmetry=1, barrier=(3.44205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149707,'amu*angstrom^2'), symmetry=1, barrier=(3.44205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149707,'amu*angstrom^2'), symmetry=1, barrier=(3.44205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149707,'amu*angstrom^2'), symmetry=1, barrier=(3.44205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149707,'amu*angstrom^2'), symmetry=1, barrier=(3.44205,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.178,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.58766,0.138091,-0.000204443,1.88717e-07,-7.0481e-11,-36943.6,50.5954], Tmin=(100,'K'), Tmax=(800.013,'K')), NASAPolynomial(coeffs=[5.29428,0.076459,-3.78434e-05,7.35166e-09,-5.13455e-13,-37173.5,24.3742], Tmin=(800.013,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-308.716,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CCJCO)"""),
)

species(
    label = 'O=CC[CH][CH]OC(O)CCC=O(21418)',
    structure = SMILES('O=CC[CH][CH]OC(O)CCC=O'),
    E0 = (-353.966,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.178,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.17759,0.146845,-0.000207001,1.70344e-07,-5.75052e-11,-42360.2,51.3067], Tmin=(100,'K'), Tmax=(772.127,'K')), NASAPolynomial(coeffs=[12.6845,0.0628149,-3.0085e-05,5.7876e-09,-4.03211e-13,-44445.5,-15.1923], Tmin=(772.127,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-353.966,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(544.598,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CCJCO) + radical(CCsJOCs)"""),
)

species(
    label = 'O=[C]CCC(O)O[CH]CCC=O(21419)',
    structure = SMILES('O=[C]CCC(O)O[CH]CCC=O'),
    E0 = (-393.907,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.178,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.51026,0.155839,-0.000233686,1.98103e-07,-6.72442e-11,-47153.8,50.8193], Tmin=(100,'K'), Tmax=(814.156,'K')), NASAPolynomial(coeffs=[13.2492,0.0616884,-2.94099e-05,5.60257e-09,-3.86317e-13,-49165.6,-18.5666], Tmin=(814.156,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-393.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(544.598,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CCsJOCs) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C(C[CH]C=O)OCCCC=O(21420)',
    structure = SMILES('[O]C=CCC([O])OCCCC=O'),
    E0 = (-365.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.178,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.27998,0.125491,-0.000132222,8.07904e-08,-2.10772e-11,-43761.3,46.745], Tmin=(100,'K'), Tmax=(905.545,'K')), NASAPolynomial(coeffs=[12.6047,0.0641569,-3.0622e-05,5.98939e-09,-4.25719e-13,-46275.9,-18.8627], Tmin=(905.545,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-365.367,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CCOJ)"""),
)

species(
    label = '[O]C(CC[C]=O)OCCCC=O(21421)',
    structure = SMILES('[O]C(CC[C]=O)OCCCC=O'),
    E0 = (-348.658,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.178,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.8752,0.146491,-0.000228736,2.12867e-07,-7.84271e-11,-41739,49.95], Tmin=(100,'K'), Tmax=(820.721,'K')), NASAPolynomial(coeffs=[5.75398,0.0755238,-3.72841e-05,7.19496e-09,-4.98974e-13,-41853.5,21.5829], Tmin=(820.721,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-348.658,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CCCJ=O)"""),
)

species(
    label = 'O=C[CH]C[CH]OC(O)CCC=O(21422)',
    structure = SMILES('[O]C=CC[CH]OC(O)CCC=O'),
    E0 = (-410.616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.178,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.19024,0.138472,-0.000151739,8.76633e-08,-2.03989e-11,-49164.6,48.5761], Tmin=(100,'K'), Tmax=(1038.71,'K')), NASAPolynomial(coeffs=[21.2206,0.0483195,-2.15519e-05,4.10714e-09,-2.88553e-13,-54028.1,-65.2571], Tmin=(1038.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-410.616,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(CCsJOCs) + radical(C=COJ)"""),
)

species(
    label = 'O=[C]CC[CH]OC(O)CCC=O(21423)',
    structure = SMILES('O=[C]CC[CH]OC(O)CCC=O'),
    E0 = (-393.907,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.178,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.51026,0.155839,-0.000233686,1.98103e-07,-6.72442e-11,-47153.8,50.8193], Tmin=(100,'K'), Tmax=(814.156,'K')), NASAPolynomial(coeffs=[13.2492,0.0616884,-2.94099e-05,5.60257e-09,-3.86317e-13,-49165.6,-18.5666], Tmin=(814.156,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-393.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(544.598,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CCsJOCs) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C(CCC=O)OC1CC[CH]O1(21424)',
    structure = SMILES('[O]C(CCC=O)OC1CC[CH]O1'),
    E0 = (-340.025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.178,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.04573,0.111155,-9.53622e-05,4.68138e-08,-9.73549e-12,-40714,41.945], Tmin=(100,'K'), Tmax=(1125.82,'K')), NASAPolynomial(coeffs=[14.1232,0.057261,-2.35557e-05,4.29303e-09,-2.93366e-13,-44129.5,-33.0343], Tmin=(1125.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-340.025,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(561.227,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(Tetrahydrofuran) + radical(CCOJ) + radical(CCsJOCs)"""),
)

species(
    label = 'O=CCC[CH]OC1CC[CH]OO1(21425)',
    structure = SMILES('O=CCC[CH]OC1CC[CH]OO1'),
    E0 = (-136.952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.178,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.12173,0.122882,-0.000111479,5.43089e-08,-1.06176e-11,-16241.2,43.1698], Tmin=(100,'K'), Tmax=(1236.94,'K')), NASAPolynomial(coeffs=[22.3335,0.0437999,-1.55788e-05,2.62245e-09,-1.7128e-13,-22291.2,-80.0131], Tmin=(1236.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-136.952,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(561.227,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(12dioxane) + radical(CCsJOCs) + radical(CCsJOOC)"""),
)

species(
    label = 'O=CCCC1O[CH]CC[CH]OO1(21426)',
    structure = SMILES('O=CCCC1O[CH]CC[CH]OO1'),
    E0 = (-107.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.178,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.15006,0.104887,-6.90147e-05,2.1697e-08,-2.73672e-12,-12724.4,39.8675], Tmin=(100,'K'), Tmax=(1805.73,'K')), NASAPolynomial(coeffs=[25.2579,0.0463903,-2.0423e-05,3.75752e-09,-2.53089e-13,-22261.8,-103.142], Tmin=(1805.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-107.39,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(569.541,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(Cyclooctane) + radical(CCsJOCs) + radical(CCsJOOC)"""),
)

species(
    label = 'O=CCCCOC(=O)CCC=O(21427)',
    structure = SMILES('O=CCCCOC(=O)CCC=O'),
    E0 = (-725.844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.178,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.30016,0.12907,-0.000174932,1.56149e-07,-5.82581e-11,-87120.1,44.9872], Tmin=(100,'K'), Tmax=(779.177,'K')), NASAPolynomial(coeffs=[5.29744,0.0763472,-3.71404e-05,7.20022e-09,-5.04013e-13,-87575.9,18.4763], Tmin=(779.177,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-725.844,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + group(Cds-OdCsH) + group(Cds-OdCsH)"""),
)

species(
    label = 'O=CCC=COC(O)CCC=O(21428)',
    structure = SMILES('O=CCC=COC(O)CCC=O'),
    E0 = (-651.386,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.178,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.10968,0.122613,-0.00010481,4.38535e-08,-7.28602e-12,-78114.2,46.613], Tmin=(100,'K'), Tmax=(1435.34,'K')), NASAPolynomial(coeffs=[27.8091,0.0392355,-1.7676e-05,3.38296e-09,-2.37071e-13,-86702.9,-108.54], Tmin=(1435.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-651.386,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-OdCsH)"""),
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
    label = 'CC[CH]OC([O])CCC=O(20314)',
    structure = SMILES('CC[CH]OC([O])CCC=O'),
    E0 = (-221.583,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2782.5,750,1395,475,1775,1000,3025,407.5,1350,352.5,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    molecularWeight = (144.168,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4428.2,'J/mol'), sigma=(7.53396,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=691.67 K, Pc=23.5 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.393274,0.0980395,-3.14934e-05,-1.65986e-07,1.9234e-10,-26540.1,36.543], Tmin=(100,'K'), Tmax=(448.276,'K')), NASAPolynomial(coeffs=[8.26909,0.0626947,-3.01108e-05,5.78963e-09,-4.02579e-13,-27597.2,0.951149], Tmin=(448.276,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-221.583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CCsJOCs)"""),
)

species(
    label = 'CCC([O])O[CH]CCC=O(20313)',
    structure = SMILES('CCC([O])O[CH]CCC=O'),
    E0 = (-221.583,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2782.5,750,1395,475,1775,1000,3025,407.5,1350,352.5,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    molecularWeight = (144.168,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4428.2,'J/mol'), sigma=(7.53396,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=691.67 K, Pc=23.5 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.393274,0.0980395,-3.14934e-05,-1.65986e-07,1.9234e-10,-26540.1,36.543], Tmin=(100,'K'), Tmax=(448.276,'K')), NASAPolynomial(coeffs=[8.26909,0.0626947,-3.01108e-05,5.78963e-09,-4.02579e-13,-27597.2,0.951149], Tmin=(448.276,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-221.583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCsJOCs) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C(CC=O)OC([O])CCC=O(21429)',
    structure = SMILES('[CH2]C(CC=O)OC([O])CCC=O'),
    E0 = (-308.659,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3100,440,815,1455,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.18509,0.154451,-0.000247446,2.30215e-07,-8.42075e-11,-36918.3,49.6506], Tmin=(100,'K'), Tmax=(825.083,'K')), NASAPolynomial(coeffs=[6.68818,0.0751703,-3.73893e-05,7.22066e-09,-5.00283e-13,-37148.2,16.0282], Tmin=(825.083,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-308.659,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CJC(C)OC)"""),
)

species(
    label = 'O=CCCC1OC(CCC=O)O1(21396)',
    structure = SMILES('O=CCCC1OC(CCC=O)O1'),
    E0 = (-585.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.178,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.292577,0.10169,-7.02127e-05,2.40728e-08,-3.47733e-12,-70307.5,39.1936], Tmin=(100,'K'), Tmax=(1505.02,'K')), NASAPolynomial(coeffs=[15.6038,0.0594411,-2.81044e-05,5.42044e-09,-3.78968e-13,-75092.4,-43.9958], Tmin=(1505.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-585.792,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(557.07,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + ring(Cyclobutane)"""),
)

species(
    label = 'C=COC[CH]OC([O])CCC=O(21430)',
    structure = SMILES('C=COC[CH]OC([O])CCC=O'),
    E0 = (-286.566,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2782.5,750,1395,475,1775,1000,180,180,180,180,883.682,1327.67,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157178,'amu*angstrom^2'), symmetry=1, barrier=(3.61382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157178,'amu*angstrom^2'), symmetry=1, barrier=(3.61382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157178,'amu*angstrom^2'), symmetry=1, barrier=(3.61382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157178,'amu*angstrom^2'), symmetry=1, barrier=(3.61382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157178,'amu*angstrom^2'), symmetry=1, barrier=(3.61382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157178,'amu*angstrom^2'), symmetry=1, barrier=(3.61382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157178,'amu*angstrom^2'), symmetry=1, barrier=(3.61382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157178,'amu*angstrom^2'), symmetry=1, barrier=(3.61382,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.178,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.20689,0.14299,-0.000165629,1.03334e-07,-2.62025e-11,-34247.8,46.8441], Tmin=(100,'K'), Tmax=(952.479,'K')), NASAPolynomial(coeffs=[19.0788,0.0535974,-2.48482e-05,4.79649e-09,-3.38546e-13,-38302.6,-54.8104], Tmin=(952.479,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-286.566,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCsJOCs) + radical(CCOJ)"""),
)

species(
    label = '[O]C(CCC=O)O[CH]CC=CO(21431)',
    structure = SMILES('[O]C(CCC=O)O[CH]CC=CO'),
    E0 = (-326.374,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2782.5,750,1395,475,1775,1000,180,180,180,180,1570.68,1600,2925.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.150767,'amu*angstrom^2'), symmetry=1, barrier=(3.46642,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150767,'amu*angstrom^2'), symmetry=1, barrier=(3.46642,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150767,'amu*angstrom^2'), symmetry=1, barrier=(3.46642,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150767,'amu*angstrom^2'), symmetry=1, barrier=(3.46642,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150767,'amu*angstrom^2'), symmetry=1, barrier=(3.46642,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150767,'amu*angstrom^2'), symmetry=1, barrier=(3.46642,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150767,'amu*angstrom^2'), symmetry=1, barrier=(3.46642,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150767,'amu*angstrom^2'), symmetry=1, barrier=(3.46642,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.178,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.22336,0.142388,-0.000165824,1.03711e-07,-2.62457e-11,-39034.2,48.1534], Tmin=(100,'K'), Tmax=(956.679,'K')), NASAPolynomial(coeffs=[19.4964,0.0515729,-2.34312e-05,4.48173e-09,-3.14742e-13,-43189.9,-55.6698], Tmin=(956.679,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-326.374,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CCsJOCs)"""),
)

species(
    label = 'C=COCC([O])O[CH]CCC=O(21432)',
    structure = SMILES('C=COCC([O])O[CH]CCC=O'),
    E0 = (-286.566,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2782.5,750,1395,475,1775,1000,180,180,180,180,883.682,1327.67,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157178,'amu*angstrom^2'), symmetry=1, barrier=(3.61382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157178,'amu*angstrom^2'), symmetry=1, barrier=(3.61382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157178,'amu*angstrom^2'), symmetry=1, barrier=(3.61382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157178,'amu*angstrom^2'), symmetry=1, barrier=(3.61382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157178,'amu*angstrom^2'), symmetry=1, barrier=(3.61382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157178,'amu*angstrom^2'), symmetry=1, barrier=(3.61382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157178,'amu*angstrom^2'), symmetry=1, barrier=(3.61382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157178,'amu*angstrom^2'), symmetry=1, barrier=(3.61382,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.178,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.20689,0.14299,-0.000165629,1.03334e-07,-2.62025e-11,-34247.8,46.8441], Tmin=(100,'K'), Tmax=(952.479,'K')), NASAPolynomial(coeffs=[19.0788,0.0535974,-2.48482e-05,4.79649e-09,-3.38546e-13,-38302.6,-54.8104], Tmin=(952.479,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-286.566,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CCsJOCs) + radical(CCOJ)"""),
)

species(
    label = '[O]C(CC=CO)O[CH]CCC=O(21433)',
    structure = SMILES('[O]C(CC=CO)O[CH]CCC=O'),
    E0 = (-326.374,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2782.5,750,1395,475,1775,1000,180,180,180,180,1570.68,1600,2925.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.150767,'amu*angstrom^2'), symmetry=1, barrier=(3.46642,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150767,'amu*angstrom^2'), symmetry=1, barrier=(3.46642,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150767,'amu*angstrom^2'), symmetry=1, barrier=(3.46642,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150767,'amu*angstrom^2'), symmetry=1, barrier=(3.46642,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150767,'amu*angstrom^2'), symmetry=1, barrier=(3.46642,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150767,'amu*angstrom^2'), symmetry=1, barrier=(3.46642,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150767,'amu*angstrom^2'), symmetry=1, barrier=(3.46642,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150767,'amu*angstrom^2'), symmetry=1, barrier=(3.46642,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.178,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.22336,0.142388,-0.000165824,1.03711e-07,-2.62459e-11,-39034.2,48.1534], Tmin=(100,'K'), Tmax=(956.666,'K')), NASAPolynomial(coeffs=[19.4964,0.0515729,-2.34312e-05,4.48174e-09,-3.14743e-13,-43189.9,-55.6696], Tmin=(956.666,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-326.374,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(CCsJOCs)"""),
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
    label = 'O=CCC[CH]O[CH]CCC=O(21434)',
    structure = SMILES('O=CCC[CH]O[CH]CCC=O'),
    E0 = (-173.181,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,3000,3050,390,425,1340,1360,335,370,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    molecularWeight = (156.179,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.10613,0.145442,-0.000212266,1.77264e-07,-5.96401e-11,-20619.6,44.1443], Tmin=(100,'K'), Tmax=(811.365,'K')), NASAPolynomial(coeffs=[12.6999,0.0590344,-2.77229e-05,5.25507e-09,-3.61763e-13,-22580.7,-21.4702], Tmin=(811.365,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-173.181,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CCsJOCs) + radical(CCsJOCs)"""),
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
    E0 = (-328.163,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-222.395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-301.777,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-306.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-232.282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-207.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-273.619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-254.944,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-253.02,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-304.005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-199.932,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-170.126,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-137.853,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-200.551,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-248.404,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-226.017,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-231.512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-261.981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-277.285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-295.426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-301.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-267.369,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-244.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-251.774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (38.1442,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-270.005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-136.952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-107.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-264.762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-303.189,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (16.9302,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (16.9302,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (-150.033,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (-319.878,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (27.2339,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (-182.507,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (27.2339,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (-182.507,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (85.3601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (69.8241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    products = ['O=CCCC=O(5767)', 'O=CCCC=O(5767)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    products = ['[O]C(CCC=O)OC1CCC1[O](21404)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.13771e+06,'s^-1'), n=1.58803, Ea=(105.768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_csHNd] for rate rule [R5_SS_CO;carbonylbond_intra_H;radadd_intra_csHNd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 101.2 to 105.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    products = ['[O]C1CCC(O[CH]CCC=O)O1(21405)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(132522,'s^-1'), n=1.48406, Ea=(26.3859,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSS;multiplebond_intra;radadd_intra] for rate rule [R6_SSS_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    products = ['[O]C1CCC([O])C(CCC=O)O1(21406)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2305.07,'s^-1'), n=1.56671, Ea=(21.6023,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_SSSS;multiplebond_intra;radadd_intra_csHNd] for rate rule [R7_SSSS_CO;carbonylbond_intra_H;radadd_intra_csHNd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    products = ['[O]C1CC[CH]OC(CCC=O)O1(21407)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.3251e+11,'s^-1'), n=0.157, Ea=(95.8805,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7plus;multiplebond_intra;radadd_intra_O] + [R8;multiplebond_intra;radadd_intra] for rate rule [R8;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[O]C(CCC=O)OC=CCC=O(21408)'],
    products = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.72e+08,'cm^3/(mol*s)'), n=1.477, Ea=(6.73624,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2825 used for Cds-CsH_Cds-OsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-OsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', 'O=CCC[CH]OC(=O)CCC=O(21409)'],
    products = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.83701,'m^3/(mol*s)'), n=1.71338, Ea=(46.5052,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;HJ] + [CO-NdNd_O;YJ] for rate rule [CO-NdNd_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O=CCCC=O(5767)', '[O][CH]CCC=O(5764)'],
    products = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(8e+09,'cm^3/(mol*s)'), n=1.39, Ea=(35.8862,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Od_CO-CsH;YJ] for rate rule [Od_CO-CsH;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CH2CHO(40)', 'C=COC([O])CCC=O(18636)'],
    products = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.0114756,'m^3/(mol*s)'), n=2.44484, Ea=(36.1431,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-HH_Cds;CsJ-OneDeHH] for rate rule [Cds-HH_Cds-OsH;CsJ-COHH]
Euclidian distance = 1.41421356237
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CH2CH2CHO(560)', 'O=CCC[CH]OC=O(18676)'],
    products = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.94e+10,'cm^3/(mol*s)'), n=0, Ea=(28.0328,'kJ/mol'), T0=(1,'K'), Tmin=(333,'K'), Tmax=(363,'K'), comment="""Estimated using template [CO_O;CsJ-CsHH] for rate rule [CO-NdH_O;CsJ-CsHH]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C(CCC=O)OC[CH]CC=O(21410)'],
    products = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.4e-20,'s^-1'), n=9.13, Ea=(108.784,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 341 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    products = ['O=CCC[CH]O[C](O)CCC=O(21411)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.43381e+07,'s^-1'), n=1.70481, Ea=(158.036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_NonDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O][C](CCC=O)OCCCC=O(21412)'],
    products = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.53164e+09,'s^-1'), n=1.09, Ea=(165.519,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_O;Y_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    products = ['[O]C(CCC=O)OCC[CH]C=O(21413)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(6.82e+09,'s^-1'), n=0.73, Ea=(127.612,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_1H;Cs_H_out_H/OneDe] for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeO;Cs_H_out_H/CO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O=CC[CH]C(O)O[CH]CCC=O(21414)'],
    products = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(5.71,'s^-1'), n=3.021, Ea=(105.562,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 319 used for R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    products = ['O=C[CH]CC(O)O[CH]CCC=O(21415)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(363473,'s^-1'), n=1.92229, Ea=(102.145,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R4H_SSS;O_rad_out;Cs_H_out_H/CO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    products = ['[O]C(CCC=O)OCCC[C]=O(21416)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(27900,'s^-1'), n=1.97, Ea=(96.6504,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/NonDeO;XH_out] for rate rule [R4H_SSS;C_rad_out_H/NonDeO;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]C([CH]CC=O)OCCCC=O(21417)'],
    products = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.00228,'s^-1'), n=3.95, Ea=(46.7353,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC] for rate rule [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    products = ['O=CC[CH][CH]OC(O)CCC=O(21418)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.19599e+09,'s^-1'), n=0.63, Ea=(50.8774,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5HJ_3;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    products = ['O=[C]CCC(O)O[CH]CCC=O(21419)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(214655,'s^-1'), n=1.70206, Ea=(32.737,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;O_rad_out;XH_out] for rate rule [R5H_CCC_O;O_rad_out;CO_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    products = ['[O]C(C[CH]C=O)OCCCC=O(21420)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.0756983,'s^-1'), n=3.26, Ea=(26.8822,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;C_rad_out_H/NonDeC;Cs_H_out_1H] for rate rule [R5H_SSSS_OCC;C_rad_out_H/NonDeC;Cs_H_out_H/CO]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    products = ['[O]C(CC[C]=O)OCCCC=O(21421)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(46.1,'s^-1'), n=3.21, Ea=(60.7935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_SSSSS;C_rad_out_H/NonDeC;XH_out] for rate rule [R6H_SSSSS_OO;C_rad_out_H/NonDeC;CO_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    products = ['O=C[CH]C[CH]OC(O)CCC=O(21422)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(49528.1,'s^-1'), n=1.95205, Ea=(83.4098,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R6HJ_3;O_rad_out;Cs_H_out_H/CO]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    products = ['O=[C]CC[CH]OC(O)CCC=O(21423)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(970995,'s^-1'), n=0.905106, Ea=(76.3888,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7Hall;O_rad_out;XH_out] for rate rule [R7HJ_3;O_rad_out;CO_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O][CH]CCC=O(5764)', '[O][CH]CCC=O(5764)'],
    products = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    products = ['[O]C(CCC=O)OC1CC[CH]O1(21424)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.64e+06,'s^-1'), n=1.15, Ea=(58.1576,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_csHNd] for rate rule [R5_SS_CO;carbonyl_intra_H;radadd_intra_csHO]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    products = ['O=CCC[CH]OC1CC[CH]OO1(21425)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(463580,'s^-1'), n=1.14062, Ea=(191.211,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSS;multiplebond_intra;radadd_intra] for rate rule [R6_SSS_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic
Ea raised from 187.7 to 191.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    products = ['O=CCCC1O[CH]CC[CH]OO1(21426)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(5.649e+12,'s^-1'), n=0.287, Ea=(220.773,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;carbonyl_intra_H;radadd_intra] for rate rule [R8_linear;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic
Ea raised from 215.4 to 220.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction29',
    reactants = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    products = ['O=CCCCOC(=O)CCC=O(21427)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    products = ['O=CCC=COC(O)CCC=O(21428)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['CO(12)', 'CC[CH]OC([O])CCC=O(20314)'],
    products = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(274200,'cm^3/(mol*s)'), n=2.53, Ea=(357.732,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 4 used for CO;C_pri/NonDeC
Exact match found for rate rule [CO;C_pri/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction32',
    reactants = ['CO(12)', 'CCC([O])O[CH]CCC=O(20313)'],
    products = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(274200,'cm^3/(mol*s)'), n=2.53, Ea=(357.732,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 4 used for CO;C_pri/NonDeC
Exact match found for rate rule [CO;C_pri/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(CC=O)OC([O])CCC=O(21429)'],
    products = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    products = ['O=CCCC1OC(CCC=O)O1(21396)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=COC[CH]OC([O])CCC=O(21430)'],
    products = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 7 used for R_ROR;R1_doublebond_CH2;R2_doublebond_H;R_O_C
Exact match found for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[O]C(CCC=O)O[CH]CC=CO(21431)'],
    products = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_H;R_O_H]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=COCC([O])O[CH]CCC=O(21432)'],
    products = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 7 used for R_ROR;R1_doublebond_CH2;R2_doublebond_H;R_O_C
Exact match found for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[O]C(CC=CO)O[CH]CCC=O(21433)'],
    products = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_H;R_O_H]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[O]C([O])CCC=O(18653)', '[CH]CCC=O(18484)'],
    products = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(2711.41,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction40',
    reactants = ['O(4)', 'O=CCC[CH]O[CH]CCC=O(21434)'],
    products = ['[O]C(CCC=O)O[CH]CCC=O(21393)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(4171.1,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/CsO;O_birad]
Euclidian distance = 4.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

network(
    label = '3808',
    isomers = [
        '[O]C(CCC=O)O[CH]CCC=O(21393)',
    ],
    reactants = [
        ('O=CCCC=O(5767)', 'O=CCCC=O(5767)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '3808',
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

