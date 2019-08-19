species(
    label = '[CH]=CC([O])(CCCC)CCOO(30981)',
    structure = SMILES('[CH]=CC([O])(CCCC)CCOO'),
    E0 = (54.6067,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.09488,0.140261,-0.000135121,7.31214e-08,-1.66767e-11,6782,49.83], Tmin=(100,'K'), Tmax=(1033.05,'K')), NASAPolynomial(coeffs=[16.4234,0.0685578,-3.10076e-05,5.9336e-09,-4.17176e-13,2955.93,-40.112], Tmin=(1033.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(54.6067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(644.372,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(Cds_P)"""),
)

species(
    label = 'C2H2(32)',
    structure = SMILES('C#C'),
    E0 = (217.784,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,584.389,584.389,2772.01],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.80868,0.0233616,-3.55172e-05,2.80153e-08,-8.50075e-12,26429,13.9397], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.65878,0.00488397,-1.60829e-06,2.46975e-10,-1.38606e-14,25759.4,-3.99838], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(217.784,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(87.302,'J/(mol*K)'), label="""C2H2""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'CCCCC(=O)CCOO(16766)',
    structure = SMILES('CCCCC(=O)CCOO'),
    E0 = (-419.494,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,375,552.5,462.5,1710,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.184,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4468.5,'J/mol'), sigma=(7.51027,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=697.97 K, Pc=23.94 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.021909,0.106478,-5.11605e-05,-1.11163e-07,1.3815e-10,-50327.2,35.0186], Tmin=(100,'K'), Tmax=(471.151,'K')), NASAPolynomial(coeffs=[8.18504,0.0687158,-3.25389e-05,6.25288e-09,-4.36324e-13,-51454.8,-2.15787], Tmin=(471.151,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-419.494,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs)"""),
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
    label = 'C#CC([O])(CCCC)CCOO(31208)',
    structure = SMILES('C#CC([O])(CCCC)CCOO'),
    E0 = (-26.7261,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,750,770,3400,2100,2175,525,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.51878,0.145693,-0.000152125,8.75136e-08,-2.05962e-11,-2981.52,47.4049], Tmin=(100,'K'), Tmax=(1021.56,'K')), NASAPolynomial(coeffs=[19.5621,0.0592348,-2.51768e-05,4.6692e-09,-3.22561e-13,-7493,-59.5941], Tmin=(1021.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-26.7261,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)2OJ)"""),
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
    label = '[CH]=CC(=O)CCOO(31209)',
    structure = SMILES('[CH]=CC(=O)CCOO'),
    E0 = (6.86311,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,3615,1310,387.5,850,1000,375,552.5,462.5,1710,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.669852,'amu*angstrom^2'), symmetry=1, barrier=(15.4012,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182396,'amu*angstrom^2'), symmetry=1, barrier=(4.19364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182946,'amu*angstrom^2'), symmetry=1, barrier=(4.2063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.181274,'amu*angstrom^2'), symmetry=1, barrier=(4.16784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.87648,'amu*angstrom^2'), symmetry=1, barrier=(43.144,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.280826,0.0920268,-0.000144296,1.33403e-07,-4.9286e-11,949.445,29.9539], Tmin=(100,'K'), Tmax=(800.263,'K')), NASAPolynomial(coeffs=[6.0839,0.0451751,-2.30282e-05,4.51361e-09,-3.16247e-13,592.089,6.8208], Tmin=(800.263,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.86311,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'CH2CH2OOH(45)',
    structure = SMILES('[CH2]COO'),
    E0 = (32.3724,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.00650046,'amu*angstrom^2'), symmetry=1, barrier=(11.0723,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2783,'amu*angstrom^2'), symmetry=1, barrier=(29.3906,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.22246,'amu*angstrom^2'), symmetry=1, barrier=(74.0907,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (61.0599,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3284.22,'J/mol'), sigma=(4.037,'angstroms'), dipoleMoment=(1.3,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[5.07456,0.013239,2.53759e-05,-4.3809e-08,1.89062e-11,3710.63,2.7229], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[9.4138,0.0133542,-4.68068e-06,7.43231e-10,-4.39824e-14,1984.6,-22.4517], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(32.3724,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), label="""CH2CH2OOH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[CH]=CC(=O)CCCC(17476)',
    structure = SMILES('[CH]=CC(=O)CCCC'),
    E0 = (47.6251,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.162,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25226,0.0714735,-2.6823e-05,-6.17328e-08,6.9559e-11,5815.79,27.67], Tmin=(100,'K'), Tmax=(478.18,'K')), NASAPolynomial(coeffs=[4.8921,0.0562846,-2.70419e-05,5.30423e-09,-3.77501e-13,5293.24,10.9712], Tmin=(478.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(47.6251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C2H2(T)(175)',
    structure = SMILES('[CH]=[CH]'),
    E0 = (590.758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([597.295,1198.22,1198.22,1198.22,2685.93,3758.47],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88057,-0.000843296,2.04089e-05,-2.51937e-08,9.47543e-12,71059.3,4.72119], Tmin=(100,'K'), Tmax=(935.375,'K')), NASAPolynomial(coeffs=[4.92145,0.00358266,-9.24397e-07,1.57263e-10,-1.19532e-14,70476.2,-2.30676], Tmin=(935.375,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(590.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""C2H2(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CCCC[C]([O])CCOO(16758)',
    structure = SMILES('CCCC[C]([O])CCOO'),
    E0 = (-69.4491,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,360,370,350,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1000,1200,1400,1600],'cm^-1')),
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
    molecularWeight = (146.184,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.38813,0.12719,-0.000146358,1.00264e-07,-2.90472e-11,-8166.49,40.7595], Tmin=(100,'K'), Tmax=(826.048,'K')), NASAPolynomial(coeffs=[11.9547,0.0625799,-2.90356e-05,5.57995e-09,-3.91711e-13,-10370.9,-21.0624], Tmin=(826.048,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-69.4491,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(C2CsJOH)"""),
)

species(
    label = 'C=[C]C([O])(CCCC)CCOO(29324)',
    structure = SMILES('C=[C]C([O])(CCCC)CCOO'),
    E0 = (45.3523,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,2950,3100,1380,975,1025,1650,1685,370,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.0336,0.140973,-0.000141456,8.20703e-08,-2.02623e-11,5664.95,49.6893], Tmin=(100,'K'), Tmax=(955.71,'K')), NASAPolynomial(coeffs=[14.499,0.0717773,-3.28504e-05,6.3104e-09,-4.44301e-13,2504.92,-29.3216], Tmin=(955.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(45.3523,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(644.372,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH]=CC(O)([CH]CCC)CCOO(31210)',
    structure = SMILES('[CH]=CC(O)([CH]CCC)CCOO'),
    E0 = (25.4061,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3615,1310,387.5,850,1000,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,180,180,180,180,850.298,1360.58,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156387,'amu*angstrom^2'), symmetry=1, barrier=(3.59564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156387,'amu*angstrom^2'), symmetry=1, barrier=(3.59564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156387,'amu*angstrom^2'), symmetry=1, barrier=(3.59564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156387,'amu*angstrom^2'), symmetry=1, barrier=(3.59564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156387,'amu*angstrom^2'), symmetry=1, barrier=(3.59564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156387,'amu*angstrom^2'), symmetry=1, barrier=(3.59564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156387,'amu*angstrom^2'), symmetry=1, barrier=(3.59564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156387,'amu*angstrom^2'), symmetry=1, barrier=(3.59564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156387,'amu*angstrom^2'), symmetry=1, barrier=(3.59564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156387,'amu*angstrom^2'), symmetry=1, barrier=(3.59564,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.57084,0.144915,-0.000142917,7.61679e-08,-1.6589e-11,3292.01,52.8088], Tmin=(100,'K'), Tmax=(1098.64,'K')), NASAPolynomial(coeffs=[21.1221,0.0586532,-2.51424e-05,4.70187e-09,-3.26811e-13,-1914.04,-63.7249], Tmin=(1098.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(25.4061,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(640.214,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CCJCO)"""),
)

species(
    label = '[CH]=CC(O)([CH]COO)CCCC(31211)',
    structure = SMILES('[CH]=CC(O)([CH]COO)CCCC'),
    E0 = (25.919,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3615,1310,387.5,850,1000,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,180,180,180,180,841.206,1373.46,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156296,'amu*angstrom^2'), symmetry=1, barrier=(3.59356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156296,'amu*angstrom^2'), symmetry=1, barrier=(3.59356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156296,'amu*angstrom^2'), symmetry=1, barrier=(3.59356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156296,'amu*angstrom^2'), symmetry=1, barrier=(3.59356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156296,'amu*angstrom^2'), symmetry=1, barrier=(3.59356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156296,'amu*angstrom^2'), symmetry=1, barrier=(3.59356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156296,'amu*angstrom^2'), symmetry=1, barrier=(3.59356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156296,'amu*angstrom^2'), symmetry=1, barrier=(3.59356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156296,'amu*angstrom^2'), symmetry=1, barrier=(3.59356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156296,'amu*angstrom^2'), symmetry=1, barrier=(3.59356,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.59622,0.150292,-0.000160239,9.46509e-08,-2.30122e-11,3350.55,52.7205], Tmin=(100,'K'), Tmax=(986.583,'K')), NASAPolynomial(coeffs=[18.8543,0.0633223,-2.80096e-05,5.29802e-09,-3.69984e-13,-881.942,-50.4758], Tmin=(986.583,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(25.919,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(640.214,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCJCOOH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C(O)(CCCC)CCOO(31212)',
    structure = SMILES('[CH]=[C]C(O)(CCCC)CCOO'),
    E0 = (63.3458,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1685,370,3615,1310,387.5,850,1000,2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,180,180,180,180,880.508,1337.5,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15702,'amu*angstrom^2'), symmetry=1, barrier=(3.61021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15702,'amu*angstrom^2'), symmetry=1, barrier=(3.61021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15702,'amu*angstrom^2'), symmetry=1, barrier=(3.61021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15702,'amu*angstrom^2'), symmetry=1, barrier=(3.61021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15702,'amu*angstrom^2'), symmetry=1, barrier=(3.61021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15702,'amu*angstrom^2'), symmetry=1, barrier=(3.61021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15702,'amu*angstrom^2'), symmetry=1, barrier=(3.61021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15702,'amu*angstrom^2'), symmetry=1, barrier=(3.61021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15702,'amu*angstrom^2'), symmetry=1, barrier=(3.61021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15702,'amu*angstrom^2'), symmetry=1, barrier=(3.61021,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.65932,0.152879,-0.000168986,1.04376e-07,-2.65428e-11,7853.14,51.1188], Tmin=(100,'K'), Tmax=(945.536,'K')), NASAPolynomial(coeffs=[18.2541,0.064407,-2.86323e-05,5.41753e-09,-3.77994e-13,3898.27,-48.605], Tmin=(945.536,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(63.3458,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(640.214,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = 'C=CC([O])([CH]CCC)CCOO(31213)',
    structure = SMILES('C=CC([O])([CH]CCC)CCOO'),
    E0 = (7.41258,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.94507,0.133116,-0.00011618,5.53118e-08,-1.10451e-11,1103.51,51.3708], Tmin=(100,'K'), Tmax=(1168.93,'K')), NASAPolynomial(coeffs=[17.9653,0.0649843,-2.87514e-05,5.4491e-09,-3.80922e-13,-3551.26,-47.7932], Tmin=(1168.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(7.41258,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(644.372,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCJCO) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=CC([O])([CH]COO)CCCC(31214)',
    structure = SMILES('C=CC([O])([CH]COO)CCCC'),
    E0 = (7.92548,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.95923,0.138287,-0.000132508,7.22502e-08,-1.67477e-11,1161.81,51.2482], Tmin=(100,'K'), Tmax=(1012.75,'K')), NASAPolynomial(coeffs=[15.203,0.0705026,-3.21122e-05,6.16257e-09,-4.33872e-13,-2314.41,-31.7671], Tmin=(1012.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(7.92548,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(644.372,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCJCOOH)"""),
)

species(
    label = '[CH]=CC(O)(C[CH]CC)CCOO(31215)',
    structure = SMILES('[CH]=CC(O)(C[CH]CC)CCOO'),
    E0 = (19.9623,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3615,1310,387.5,850,1000,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,180,180,180,180,801.676,1420.93,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155741,'amu*angstrom^2'), symmetry=1, barrier=(3.58078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155741,'amu*angstrom^2'), symmetry=1, barrier=(3.58078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155741,'amu*angstrom^2'), symmetry=1, barrier=(3.58078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155741,'amu*angstrom^2'), symmetry=1, barrier=(3.58078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155741,'amu*angstrom^2'), symmetry=1, barrier=(3.58078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155741,'amu*angstrom^2'), symmetry=1, barrier=(3.58078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155741,'amu*angstrom^2'), symmetry=1, barrier=(3.58078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155741,'amu*angstrom^2'), symmetry=1, barrier=(3.58078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155741,'amu*angstrom^2'), symmetry=1, barrier=(3.58078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155741,'amu*angstrom^2'), symmetry=1, barrier=(3.58078,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.39004,0.148896,-0.000168327,1.1109e-07,-3.06581e-11,2623.7,52.2464], Tmin=(100,'K'), Tmax=(869.382,'K')), NASAPolynomial(coeffs=[14.8939,0.0693782,-3.11392e-05,5.89723e-09,-4.1081e-13,-381.77,-28.7206], Tmin=(869.382,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(19.9623,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(640.214,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJCC) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC(O)(C[CH]OO)CCCC(31216)',
    structure = SMILES('[CH]=CC(O)(C[CH]OO)CCCC'),
    E0 = (14.0865,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.62395,0.15107,-0.000161349,9.56483e-08,-2.33549e-11,1928.29,50.6649], Tmin=(100,'K'), Tmax=(982.285,'K')), NASAPolynomial(coeffs=[18.7699,0.0639494,-2.83092e-05,5.35412e-09,-3.73819e-13,-2274.6,-52.1651], Tmin=(982.285,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(14.0865,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(640.214,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CCsJOOH)"""),
)

species(
    label = 'C=CC([O])(C[CH]CC)CCOO(31217)',
    structure = SMILES('C=CC([O])(C[CH]CC)CCOO'),
    E0 = (1.96883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.80764,0.13747,-0.000142345,9.06505e-08,-2.51364e-11,437.424,50.9743], Tmin=(100,'K'), Tmax=(849.998,'K')), NASAPolynomial(coeffs=[11.0648,0.0768953,-3.54514e-05,6.81405e-09,-4.79213e-13,-1750.94,-9.03595], Tmin=(849.998,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1.96883,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(644.372,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJCC) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=CC([O])(C[CH]OO)CCCC(31218)',
    structure = SMILES('C=CC([O])(C[CH]OO)CCCC'),
    E0 = (-3.90698,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.9877,0.139069,-0.00013362,7.32337e-08,-1.70789e-11,-260.415,49.1955], Tmin=(100,'K'), Tmax=(1006.73,'K')), NASAPolynomial(coeffs=[15.1058,0.0711527,-3.24257e-05,6.22204e-09,-4.37994e-13,-3702.12,-33.3852], Tmin=(1006.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-3.90698,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(644.372,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCsJOOH) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH]=CC(O)(CC[CH]C)CCOO(31219)',
    structure = SMILES('[CH]=CC(O)(CC[CH]C)CCOO'),
    E0 = (19.9503,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.42242,0.148406,-0.000163288,1.03175e-07,-2.71272e-11,2624.6,52.5317], Tmin=(100,'K'), Tmax=(912.193,'K')), NASAPolynomial(coeffs=[16.0883,0.0672348,-2.98103e-05,5.6234e-09,-3.91314e-13,-752.443,-35.0706], Tmin=(912.193,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(19.9503,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(640.214,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(Cds_P)"""),
)

species(
    label = 'C=CC([O])(CC[CH]C)CCOO(31220)',
    structure = SMILES('C=CC([O])(CC[CH]C)CCOO'),
    E0 = (1.95681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.81136,0.136652,-0.0001362,8.13344e-08,-2.10081e-11,437.08,51.1563], Tmin=(100,'K'), Tmax=(910.031,'K')), NASAPolynomial(coeffs=[12.2835,0.0747001,-3.40878e-05,6.53113e-09,-4.58907e-13,-2128.33,-15.5147], Tmin=(910.031,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1.95681,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(644.372,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(RCCJC)"""),
)

species(
    label = '[CH]=CC(O)(CCC[CH2])CCOO(31221)',
    structure = SMILES('[CH]=CC(O)(CCC[CH2])CCOO'),
    E0 = (30.7504,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.60549,0.149297,-0.00015603,8.96053e-08,-2.11293e-11,3933.01,52.232], Tmin=(100,'K'), Tmax=(1016.7,'K')), NASAPolynomial(coeffs=[19.5881,0.061981,-2.72067e-05,5.1335e-09,-3.58216e-13,-579.834,-55.207], Tmin=(1016.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(30.7504,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(640.214,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC(O)(CCCC)CCO[O](31222)',
    structure = SMILES('[CH]=CC(O)(CCCC)CCO[O]'),
    E0 = (-22.4913,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.42219,0.145485,-0.0001521,8.90891e-08,-2.15352e-11,-2477.31,50.5465], Tmin=(100,'K'), Tmax=(991.515,'K')), NASAPolynomial(coeffs=[18.0512,0.0628906,-2.71486e-05,5.07527e-09,-3.52072e-13,-6537.26,-48.0515], Tmin=(991.515,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-22.4913,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(644.372,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]CCCC([O])(C=C)CCOO(31223)',
    structure = SMILES('[CH2]CCCC([O])(C=C)CCOO'),
    E0 = (12.7569,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.96608,0.137288,-0.000128384,6.7422e-08,-1.49909e-11,1744.1,50.7492], Tmin=(100,'K'), Tmax=(1054.32,'K')), NASAPolynomial(coeffs=[16.0384,0.0689811,-3.12021e-05,5.97213e-09,-4.19919e-13,-2052.41,-37.0645], Tmin=(1054.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(12.7569,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(644.372,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(RCCJ)"""),
)

species(
    label = 'C=CC([O])(CCCC)CCO[O](31224)',
    structure = SMILES('C=CC([O])(CCCC)CCO[O]'),
    E0 = (-40.4848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.78526,0.133484,-0.000124397,6.67401e-08,-1.52976e-11,-4666.06,49.0742], Tmin=(100,'K'), Tmax=(1021.45,'K')), NASAPolynomial(coeffs=[14.4144,0.0700452,-3.1236e-05,5.93613e-09,-4.15651e-13,-7975.45,-29.4234], Tmin=(1021.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-40.4848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(648.529,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(ROOJ)"""),
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
    label = '[CH]=CC([O])(CCC)CCOO(14889)',
    structure = SMILES('[CH]=CC([O])(CCC)CCOO'),
    E0 = (78.3869,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,180,791.148,1423.38,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155302,'amu*angstrom^2'), symmetry=1, barrier=(3.57069,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155302,'amu*angstrom^2'), symmetry=1, barrier=(3.57069,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155302,'amu*angstrom^2'), symmetry=1, barrier=(3.57069,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155302,'amu*angstrom^2'), symmetry=1, barrier=(3.57069,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155302,'amu*angstrom^2'), symmetry=1, barrier=(3.57069,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155302,'amu*angstrom^2'), symmetry=1, barrier=(3.57069,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155302,'amu*angstrom^2'), symmetry=1, barrier=(3.57069,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155302,'amu*angstrom^2'), symmetry=1, barrier=(3.57069,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.47728,0.125834,-0.000122573,6.67812e-08,-1.52779e-11,9620.57,45.3664], Tmin=(100,'K'), Tmax=(1032.25,'K')), NASAPolynomial(coeffs=[15.4511,0.0602339,-2.72451e-05,5.21369e-09,-3.6658e-13,6125.78,-36.8403], Tmin=(1032.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(78.3869,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(Cds_P)"""),
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
    label = '[CH]=CC1(CCCC)CCOO1(31225)',
    structure = SMILES('[CH]=CC1(CCCC)CCOO1'),
    E0 = (69.7573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (155.214,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.44405,0.0998356,-2.27776e-05,-5.06726e-08,3.04042e-11,8604.24,40.9906], Tmin=(100,'K'), Tmax=(923.714,'K')), NASAPolynomial(coeffs=[23.7655,0.0420276,-1.23038e-05,1.95942e-09,-1.31074e-13,1755.91,-90.4903], Tmin=(923.714,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(69.7573,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(611.114,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(12dioxolane) + radical(Cds_P)"""),
)

species(
    label = 'CCCCC1([O])C=COCC1(31226)',
    structure = SMILES('CCCCC1([O])C=COCC1'),
    E0 = (-218.251,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (155.214,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.0781,0.0874253,9.71087e-06,-7.90936e-08,3.81298e-11,-26044.8,37.4677], Tmin=(100,'K'), Tmax=(955.655,'K')), NASAPolynomial(coeffs=[23.5504,0.0435452,-1.43439e-05,2.51437e-09,-1.77763e-13,-33455.6,-94.3786], Tmin=(955.655,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-218.251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(615.271,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(3,4-Dihydro-2H-pyran) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'CCCCC1(C=CO1)CCOO(31200)',
    structure = SMILES('CCCCC1(C=CO1)CCOO'),
    E0 = (-269.761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (172.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.92022,0.127687,-6.3114e-05,-2.17705e-08,2.07322e-11,-32173.6,45.5857], Tmin=(100,'K'), Tmax=(979.453,'K')), NASAPolynomial(coeffs=[31.9298,0.0405995,-1.43365e-05,2.60792e-09,-1.87023e-13,-41649.9,-135.347], Tmin=(979.453,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-269.761,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(648.529,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene)"""),
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
    label = '[CH]=C[C](CCCC)CCOO(31227)',
    structure = SMILES('[CH]C=C(CCCC)CCOO'),
    E0 = (120.554,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2761.11,2772.22,2783.33,2794.44,2805.56,2816.67,2827.78,2838.89,2850,1425,1431.25,1437.5,1443.75,1450,1225,1237.5,1250,1262.5,1275,1270,1287.5,1305,1322.5,1340,700,725,750,775,800,300,325,350,375,400,3010,987.5,1337.5,450,1655,350,440,435,1725,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
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
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.57183,0.126658,-0.000100989,4.46669e-08,-8.51644e-12,14696.2,44.6186], Tmin=(100,'K'), Tmax=(1193.84,'K')), NASAPolynomial(coeffs=[14.8816,0.0715307,-3.17245e-05,5.98832e-09,-4.16862e-13,10767.6,-37.6748], Tmin=(1193.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(120.554,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(AllylJ2_triplet)"""),
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
    E0 = (54.6067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (205.614,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (97.9532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (108.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (219.422,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (175.806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (160.043,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (130.968,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (131.481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (187.115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (98.9152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (98.9152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (104.061,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (171.759,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (98.012,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (110.456,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (80.3383,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (116.104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (148.119,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (124.149,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (186.884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (206.127,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (521.309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (498.249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (109.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (97.5554,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (62.891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (363.558,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=CC([O])(CCCC)CCOO(30981)'],
    products = ['C2H2(32)', 'CCCCC(=O)CCOO(16766)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(3)', 'C#CC([O])(CCCC)CCOO(31208)'],
    products = ['[CH]=CC([O])(CCCC)CCOO(30981)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.276e+10,'cm^3/(mol*s)'), n=1.103, Ea=(20.5476,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 139 used for Ct-Cs_Ct-H;HJ
Exact match found for rate rule [Ct-Cs_Ct-H;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['butyl_1(82)', '[CH]=CC(=O)CCOO(31209)'],
    products = ['[CH]=CC([O])(CCCC)CCOO(30981)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(7.94e+10,'cm^3/(mol*s)'), n=0, Ea=(28.0328,'kJ/mol'), T0=(1,'K'), Tmin=(333,'K'), Tmax=(363,'K'), comment="""Estimated using template [CO_O;CsJ-CsHH] for rate rule [CO-CdCs_O;CsJ-CsHH]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2CH2OOH(45)', '[CH]=CC(=O)CCCC(17476)'],
    products = ['[CH]=CC([O])(CCCC)CCOO(30981)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7.94e+10,'cm^3/(mol*s)'), n=0, Ea=(28.0328,'kJ/mol'), T0=(1,'K'), Tmin=(333,'K'), Tmax=(363,'K'), comment="""Estimated using template [CO_O;CsJ-CsHH] for rate rule [CO-CdCs_O;CsJ-CsHH]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C2H2(T)(175)', 'CCCCC(=O)CCOO(16766)'],
    products = ['[CH]=CC([O])(CCCC)CCOO(30981)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(6.32e+10,'cm^3/(mol*s)'), n=0, Ea=(48.1578,'kJ/mol'), T0=(1,'K'), Tmin=(413,'K'), Tmax=(563,'K'), comment="""Estimated using an average for rate rule [CO-CsCs_O;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C2H2(32)', 'CCCC[C]([O])CCOO(16758)'],
    products = ['[CH]=CC([O])(CCCC)CCOO(30981)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(46.4627,'m^3/(mol*s)'), n=1.51997, Ea=(27.4714,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-H_Ct-H;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=CC([O])(CCCC)CCOO(30981)'],
    products = ['C=[C]C([O])(CCCC)CCOO(29324)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=CC(O)([CH]CCC)CCOO(31210)'],
    products = ['[CH]=CC([O])(CCCC)CCOO(30981)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(5.71,'s^-1'), n=3.021, Ea=(105.562,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 319 used for R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=CC(O)([CH]COO)CCCC(31211)'],
    products = ['[CH]=CC([O])(CCCC)CCOO(30981)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.71,'s^-1'), n=3.021, Ea=(105.562,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 319 used for R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=[C]C(O)(CCCC)CCOO(31212)'],
    products = ['[CH]=CC([O])(CCCC)CCOO(30981)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(117344,'s^-1'), n=2.01217, Ea=(123.77,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS_Cs;Y_rad_out;O_H_out] + [R3H_SS_Cs;Cd_rad_out;XH_out] for rate rule [R3H_SS_Cs;Cd_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=CC([O])(CCCC)CCOO(30981)'],
    products = ['C=CC([O])([CH]CCC)CCOO(31213)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=CC([O])(CCCC)CCOO(30981)'],
    products = ['C=CC([O])([CH]COO)CCCC(31214)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=CC(O)(C[CH]CC)CCOO(31215)'],
    products = ['[CH]=CC([O])(CCCC)CCOO(30981)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2960,'s^-1'), n=2.11, Ea=(84.0984,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 323 used for R4H_SSS;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=CC([O])(CCCC)CCOO(30981)'],
    products = ['[CH]=CC(O)(C[CH]OO)CCCC(31216)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.94e+12,'s^-1'), n=-0.24, Ea=(117.152,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;O_rad_out;Cs_H_out_OOH/H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=CC([O])(CCCC)CCOO(30981)'],
    products = ['C=CC([O])(C[CH]CC)CCOO(31217)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(366176,'s^-1'), n=1.54456, Ea=(43.4053,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 3.74165738677
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=CC([O])(CCCC)CCOO(30981)'],
    products = ['C=CC([O])(C[CH]OO)CCCC(31218)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(96695.3,'s^-1'), n=1.95578, Ea=(55.8497,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5H;Cd_rad_out_singleH;Cs_H_out] + [R5H_RSSR;Y_rad_out;Cs_H_out_OOH/H] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cs_H_out_OOH/H]
Euclidian distance = 3.60555127546
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=CC([O])(CCCC)CCOO(30981)'],
    products = ['[CH]=CC(O)(CC[CH]C)CCOO(31219)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(8e+10,'s^-1'), n=0, Ea=(25.7316,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(1000,'K'), comment="""From training reaction 307 used for R5H_CCC;O_rad_out;Cs_H_out_H/NonDeC
Exact match found for rate rule [R5H_CCC;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=CC([O])(CCCC)CCOO(30981)'],
    products = ['C=CC([O])(CC[CH]C)CCOO(31220)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(8.63981e+06,'s^-1'), n=0.849259, Ea=(61.4971,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_RSSSR;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R6H_DSSSS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 3.60555127546
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=CC([O])(CCCC)CCOO(30981)'],
    products = ['[CH]=CC(O)(CCC[CH2])CCOO(31221)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.68e+09,'s^-1'), n=0, Ea=(93.5124,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 140 used for R6H_SSSSS;O_rad_out;Cs_H_out_2H
Exact match found for rate rule [R6H_SSSSS;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=CC([O])(CCCC)CCOO(30981)'],
    products = ['[CH]=CC(O)(CCCC)CCO[O](31222)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(123903,'s^-1'), n=1.46258, Ea=(69.5427,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_SSSSS;O_rad_out;XH_out] for rate rule [R6H_SSSSS;O_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=CC([O])(CCCC)CCOO(30981)'],
    products = ['[CH2]CCCC([O])(C=C)CCOO(31223)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(22.7193,'s^-1'), n=3.21897, Ea=(132.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R7H;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=CC([O])(CCCC)CCOO(30981)'],
    products = ['C=CC([O])(CCCC)CCO[O](31224)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R7H;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C2H2(T)(175)', 'CCCC[C]([O])CCOO(16758)'],
    products = ['[CH]=CC([O])(CCCC)CCOO(30981)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.49215e+07,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['CH2(S)(23)', '[CH]=CC([O])(CCC)CCOO(14889)'],
    products = ['[CH]=CC([O])(CCCC)CCOO(30981)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=CC([O])(CCCC)CCOO(30981)'],
    products = ['OH(5)', '[CH]=CC1(CCCC)CCOO1(31225)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.63066e+10,'s^-1'), n=0, Ea=(54.8104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4OO_SSS;Y_rad_intra;OOH]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=CC([O])(CCCC)CCOO(30981)'],
    products = ['OH(5)', 'CCCCC1([O])C=COCC1(31226)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(3.85474e+08,'s^-1'), n=0.385, Ea=(42.9488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5OO;Y_rad_intra;OOH] for rate rule [R5OO_DSSS;Cd_pri_rad_in;OOH]
Euclidian distance = 2.2360679775
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]=CC([O])(CCCC)CCOO(30981)'],
    products = ['CCCCC1(C=CO1)CCOO(31200)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;O_rad;CdsinglepriH_rad_out]
Euclidian distance = 2.44948974278
family: Birad_recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['O(4)', '[CH]=C[C](CCCC)CCOO(31227)'],
    products = ['[CH]=CC([O])(CCCC)CCOO(30981)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/Cs2;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

network(
    label = '4810',
    isomers = [
        '[CH]=CC([O])(CCCC)CCOO(30981)',
    ],
    reactants = [
        ('C2H2(32)', 'CCCCC(=O)CCOO(16766)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4810',
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

