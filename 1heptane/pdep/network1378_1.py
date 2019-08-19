species(
    label = '[CH2]C(CC)C([O])(CC)C(=C)O(6341)',
    structure = SMILES('[CH2]C(CC)C([O])(CC)C(=C)O'),
    E0 = (-115.99,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,180,180,180,180,906.878,1304.46,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157678,'amu*angstrom^2'), symmetry=1, barrier=(3.62533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157678,'amu*angstrom^2'), symmetry=1, barrier=(3.62533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157678,'amu*angstrom^2'), symmetry=1, barrier=(3.62533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157678,'amu*angstrom^2'), symmetry=1, barrier=(3.62533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157678,'amu*angstrom^2'), symmetry=1, barrier=(3.62533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157678,'amu*angstrom^2'), symmetry=1, barrier=(3.62533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157678,'amu*angstrom^2'), symmetry=1, barrier=(3.62533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157678,'amu*angstrom^2'), symmetry=1, barrier=(3.62533,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.78578,0.132175,-0.000117069,5.43523e-08,-1.00362e-11,-13691.8,49.9885], Tmin=(100,'K'), Tmax=(1312.63,'K')), NASAPolynomial(coeffs=[26.5201,0.0428707,-1.5017e-05,2.52182e-09,-1.64734e-13,-21385.4,-99.3679], Tmin=(1312.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-115.99,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(623.585,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(Isobutyl)"""),
)

species(
    label = 'butene1(127)',
    structure = SMILES('C=CCC'),
    E0 = (-16.4325,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,252.555],'cm^-1')),
        HinderedRotor(inertia=(0.178654,'amu*angstrom^2'), symmetry=1, barrier=(7.72883,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185589,'amu*angstrom^2'), symmetry=1, barrier=(7.72103,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.176,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.58773,0.0232778,1.93412e-05,-3.55496e-08,1.36906e-11,-1918.73,14.5751], Tmin=(100,'K'), Tmax=(1007.28,'K')), NASAPolynomial(coeffs=[7.20517,0.0236362,-9.0315e-06,1.65393e-09,-1.16019e-13,-3797.34,-12.4426], Tmin=(1007.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-16.4325,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), label="""butene1""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2]C(CC)C1(CC)OC1([CH2])O(7138)',
    structure = SMILES('[CH2]C(CC)C1(CC)OC1([CH2])O'),
    E0 = (-73.1352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.70737,0.152982,-0.000151089,7.50923e-08,-1.39941e-11,-8449.57,49.813], Tmin=(100,'K'), Tmax=(1535.61,'K')), NASAPolynomial(coeffs=[33.4516,0.0276945,-3.41581e-06,2.55279e-12,1.81829e-14,-17116.4,-140.71], Tmin=(1535.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-73.1352,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(623.585,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(Isobutyl) + radical(CJC(O)2C)"""),
)

species(
    label = '[CH2]C1(O)CC(CC)C1([O])CC(7139)',
    structure = SMILES('[CH2]C1(O)CC(CC)C1([O])CC'),
    E0 = (-54.0457,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.48196,0.125965,-9.37403e-05,2.90062e-08,-1.33114e-12,-6252.68,42.148], Tmin=(100,'K'), Tmax=(1080.09,'K')), NASAPolynomial(coeffs=[25.4387,0.0470211,-1.80695e-05,3.26378e-09,-2.25105e-13,-13710.6,-101.308], Tmin=(1080.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-54.0457,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(631.9,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CJC(C)2O) + radical(CC(C)2OJ)"""),
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
    label = 'C=C(O)C([O])(CC)C(=C)CC(7140)',
    structure = SMILES('C=C(O)C([O])(CC)C(=C)CC'),
    E0 = (-211.289,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,325,375,415,465,420,450,1700,1750,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,180,180,180,964.259,1246.6,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.159311,'amu*angstrom^2'), symmetry=1, barrier=(3.66286,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159311,'amu*angstrom^2'), symmetry=1, barrier=(3.66286,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159311,'amu*angstrom^2'), symmetry=1, barrier=(3.66286,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159311,'amu*angstrom^2'), symmetry=1, barrier=(3.66286,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159311,'amu*angstrom^2'), symmetry=1, barrier=(3.66286,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159311,'amu*angstrom^2'), symmetry=1, barrier=(3.66286,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159311,'amu*angstrom^2'), symmetry=1, barrier=(3.66286,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (155.214,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.43044,0.123092,-9.97008e-05,4.08406e-08,-6.6595e-12,-25165.5,46.1174], Tmin=(100,'K'), Tmax=(1468.56,'K')), NASAPolynomial(coeffs=[27.656,0.0411437,-1.59987e-05,2.84339e-09,-1.91063e-13,-34002.2,-110.594], Tmin=(1468.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-211.289,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(602.799,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH2][CH]CC(130)',
    structure = SMILES('[CH2][CH]CC'),
    E0 = (255.669,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,2031.24],'cm^-1')),
        HinderedRotor(inertia=(0.244974,'amu*angstrom^2'), symmetry=1, barrier=(5.63244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00192352,'amu*angstrom^2'), symmetry=1, barrier=(5.63177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244928,'amu*angstrom^2'), symmetry=1, barrier=(5.63137,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.98997,0.0287412,-9.51469e-06,4.19232e-10,1.90526e-13,30780.1,16.8971], Tmin=(100,'K'), Tmax=(2154.56,'K')), NASAPolynomial(coeffs=[12.4231,0.0182241,-7.06316e-06,1.16769e-09,-7.11818e-14,25091.4,-39.6212], Tmin=(2154.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(255.669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJC) + radical(RCCJ)"""),
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
    label = '[CH2]C(CC)C(=O)C(=C)O(7141)',
    structure = SMILES('[CH2]C(CC)C(=O)C(=C)O'),
    E0 = (-201.278,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.161,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.722636,0.106969,-0.000115436,6.73645e-08,-1.60208e-11,-24040.6,33.0737], Tmin=(100,'K'), Tmax=(1011.67,'K')), NASAPolynomial(coeffs=[15.8,0.0416414,-1.85758e-05,3.53682e-09,-2.48151e-13,-27383.7,-46.8306], Tmin=(1011.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-201.278,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsHH) + radical(CJC(C)C=O)"""),
)

species(
    label = 'CH2COH(99)',
    structure = SMILES('C=[C]O'),
    E0 = (103.269,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3615,1277.5,1000,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.989114,'amu*angstrom^2'), symmetry=1, barrier=(22.7417,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.1624,0.0134245,5.56346e-06,-1.95511e-08,9.36369e-12,12455.2,10.1544], Tmin=(100,'K'), Tmax=(925.618,'K')), NASAPolynomial(coeffs=[8.19875,0.00453462,-8.93448e-07,1.26083e-10,-9.46513e-15,10971.3,-16.733], Tmin=(925.618,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(103.269,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CH2COH""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C(CC)C(=O)CC(7142)',
    structure = SMILES('[CH2]C(CC)C(=O)CC'),
    E0 = (-126.421,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (113.178,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.413552,0.089641,-6.89146e-05,2.23407e-10,2.80351e-11,-15086.2,30.0094], Tmin=(100,'K'), Tmax=(544.798,'K')), NASAPolynomial(coeffs=[7.14274,0.0571757,-2.61725e-05,4.99974e-09,-3.49752e-13,-16070.8,-0.675508], Tmin=(544.798,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-126.421,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJC(C)C=O)"""),
)

species(
    label = 'C=CC([O])(CC)C(=C)O(6516)',
    structure = SMILES('C=CC([O])(CC)C(=C)O'),
    E0 = (-149.589,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,422.355,712.271,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155253,'amu*angstrom^2'), symmetry=1, barrier=(3.56957,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155253,'amu*angstrom^2'), symmetry=1, barrier=(3.56957,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155253,'amu*angstrom^2'), symmetry=1, barrier=(3.56957,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155253,'amu*angstrom^2'), symmetry=1, barrier=(3.56957,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155253,'amu*angstrom^2'), symmetry=1, barrier=(3.56957,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.161,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.623671,0.0878525,-5.65832e-05,8.57472e-09,3.28341e-12,-17813.1,35.9415], Tmin=(100,'K'), Tmax=(1069.66,'K')), NASAPolynomial(coeffs=[20.4152,0.0330793,-1.32915e-05,2.48319e-09,-1.75318e-13,-23681.3,-73.3671], Tmin=(1069.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-149.589,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C(O)C([O])(CC)[C](C)CC(7143)',
    structure = SMILES('C=C(O)C([O])(CC)[C](C)CC'),
    E0 = (-168.499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.23395,0.120436,-8.18703e-05,1.82684e-08,2.20968e-12,-20027,46.415], Tmin=(100,'K'), Tmax=(1054.84,'K')), NASAPolynomial(coeffs=[24.3973,0.0475564,-1.82039e-05,3.29177e-09,-2.27738e-13,-27209,-90.8994], Tmin=(1054.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-168.499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(623.585,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCJ(C)CO) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH2][C](CC)C(O)(CC)C(=C)O(7144)',
    structure = SMILES('[CH2][C](CC)C(O)(CC)C(=C)O'),
    E0 = (-192.519,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.48735,0.127184,-9.36285e-05,2.07246e-08,4.54132e-12,-22907.2,47.8887], Tmin=(100,'K'), Tmax=(966.37,'K')), NASAPolynomial(coeffs=[26.3392,0.0419416,-1.4209e-05,2.42528e-09,-1.6448e-13,-30069.8,-98.4292], Tmin=(966.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-192.519,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Isobutyl) + radical(CCJ(C)CO)"""),
)

species(
    label = '[CH2]C(CC)C(O)([CH]C)C(=C)O(7145)',
    structure = SMILES('[CH2]C(CC)C(O)([CH]C)C(=C)O'),
    E0 = (-145.191,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,412.775,716.937,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155306,'amu*angstrom^2'), symmetry=1, barrier=(3.57078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155306,'amu*angstrom^2'), symmetry=1, barrier=(3.57078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155306,'amu*angstrom^2'), symmetry=1, barrier=(3.57078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155306,'amu*angstrom^2'), symmetry=1, barrier=(3.57078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155306,'amu*angstrom^2'), symmetry=1, barrier=(3.57078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155306,'amu*angstrom^2'), symmetry=1, barrier=(3.57078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155306,'amu*angstrom^2'), symmetry=1, barrier=(3.57078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155306,'amu*angstrom^2'), symmetry=1, barrier=(3.57078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155306,'amu*angstrom^2'), symmetry=1, barrier=(3.57078,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.55152,0.128461,-9.57674e-05,2.01667e-08,5.63054e-12,-17212.4,50.4207], Tmin=(100,'K'), Tmax=(952.449,'K')), NASAPolynomial(coeffs=[27.2692,0.0395228,-1.28693e-05,2.15754e-09,-1.45662e-13,-24539.5,-100.638], Tmin=(952.449,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-145.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCJCO) + radical(Isobutyl)"""),
)

species(
    label = 'C=C(O)C([O])(CC)C(C)[CH]C(7146)',
    structure = SMILES('C=C(O)C([O])(CC)C(C)[CH]C'),
    E0 = (-126.531,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.34966,0.126144,-0.000104865,4.53401e-08,-7.88951e-12,-14978.4,49.1989], Tmin=(100,'K'), Tmax=(1371.92,'K')), NASAPolynomial(coeffs=[24.8232,0.0469177,-1.82416e-05,3.24626e-09,-2.18849e-13,-22434.1,-90.4867], Tmin=(1371.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-126.531,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(623.585,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cs_S) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C(O)C([O])([CH]C)C(C)CC(7147)',
    structure = SMILES('C=C(O)C([O])([CH]C)C(C)CC'),
    E0 = (-121.171,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.29146,0.121666,-8.40146e-05,1.79883e-08,3.04901e-12,-14332.4,48.9211], Tmin=(100,'K'), Tmax=(1031.22,'K')), NASAPolynomial(coeffs=[25.1394,0.045441,-1.70326e-05,3.06268e-09,-2.12057e-13,-21594.4,-92.04], Tmin=(1031.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-121.171,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(623.585,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCJCO) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH2]C([CH]C)C(O)(CC)C(=C)O(7148)',
    structure = SMILES('[CH2]C([CH]C)C(O)(CC)C(=C)O'),
    E0 = (-150.551,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.87685,0.136173,-0.000128361,6.34268e-08,-1.24039e-11,-17846.8,51.6516], Tmin=(100,'K'), Tmax=(1246.12,'K')), NASAPolynomial(coeffs=[26.6579,0.0413656,-1.42369e-05,2.36987e-09,-1.54316e-13,-25207.4,-97.3345], Tmin=(1246.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-150.551,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]CC(O)(C(=C)O)C([CH2])CC(7149)',
    structure = SMILES('[CH2]CC(O)(C(=C)O)C([CH2])CC'),
    E0 = (-139.847,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,180,351.067,782.299,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154392,'amu*angstrom^2'), symmetry=1, barrier=(3.54978,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154392,'amu*angstrom^2'), symmetry=1, barrier=(3.54978,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154392,'amu*angstrom^2'), symmetry=1, barrier=(3.54978,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154392,'amu*angstrom^2'), symmetry=1, barrier=(3.54978,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154392,'amu*angstrom^2'), symmetry=1, barrier=(3.54978,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154392,'amu*angstrom^2'), symmetry=1, barrier=(3.54978,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154392,'amu*angstrom^2'), symmetry=1, barrier=(3.54978,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154392,'amu*angstrom^2'), symmetry=1, barrier=(3.54978,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154392,'amu*angstrom^2'), symmetry=1, barrier=(3.54978,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.22518,0.140403,-0.000135317,6.76615e-08,-1.32638e-11,-16543.8,52.1332], Tmin=(100,'K'), Tmax=(1267.87,'K')), NASAPolynomial(coeffs=[28.9524,0.0374724,-1.18691e-05,1.87151e-09,-1.17927e-13,-24589.7,-110.294], Tmin=(1267.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-139.847,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]CC(C)C([O])(CC)C(=C)O(7150)',
    structure = SMILES('[CH2]CC(C)C([O])(CC)C(=C)O'),
    E0 = (-115.826,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,180,180,180,180,973.971,1237.53,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.159114,'amu*angstrom^2'), symmetry=1, barrier=(3.65834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159114,'amu*angstrom^2'), symmetry=1, barrier=(3.65834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159114,'amu*angstrom^2'), symmetry=1, barrier=(3.65834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159114,'amu*angstrom^2'), symmetry=1, barrier=(3.65834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159114,'amu*angstrom^2'), symmetry=1, barrier=(3.65834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159114,'amu*angstrom^2'), symmetry=1, barrier=(3.65834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159114,'amu*angstrom^2'), symmetry=1, barrier=(3.65834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159114,'amu*angstrom^2'), symmetry=1, barrier=(3.65834,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.70924,0.130478,-0.000112078,4.98031e-08,-8.81552e-12,-13674.9,49.7232], Tmin=(100,'K'), Tmax=(1362.26,'K')), NASAPolynomial(coeffs=[27.0151,0.0431992,-1.59747e-05,2.77178e-09,-1.84445e-13,-21773.4,-102.869], Tmin=(1362.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-115.826,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(623.585,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C(CC)C(O)(CC)C(=C)[O](7151)',
    structure = SMILES('[CH2]C(CC)C(O)(CC)C(=C)[O]'),
    E0 = (-207.288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.66196,0.132979,-0.00012309,6.07238e-08,-1.19414e-11,-24679.7,49.8235], Tmin=(100,'K'), Tmax=(1235.01,'K')), NASAPolynomial(coeffs=[24.7323,0.0442539,-1.53274e-05,2.55263e-09,-1.65955e-13,-31446.2,-88.1201], Tmin=(1235.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-207.288,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(623.585,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=C(O)C(O)(CC)C([CH2])CC(7152)',
    structure = SMILES('[CH]=C(O)C(O)(CC)C([CH2])CC'),
    E0 = (-97.9968,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,180,180,180,180,1581.53,1600,2911.07,3200],'cm^-1')),
        HinderedRotor(inertia=(0.151578,'amu*angstrom^2'), symmetry=1, barrier=(3.48507,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151578,'amu*angstrom^2'), symmetry=1, barrier=(3.48507,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151578,'amu*angstrom^2'), symmetry=1, barrier=(3.48507,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151578,'amu*angstrom^2'), symmetry=1, barrier=(3.48507,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151578,'amu*angstrom^2'), symmetry=1, barrier=(3.48507,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151578,'amu*angstrom^2'), symmetry=1, barrier=(3.48507,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151578,'amu*angstrom^2'), symmetry=1, barrier=(3.48507,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151578,'amu*angstrom^2'), symmetry=1, barrier=(3.48507,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151578,'amu*angstrom^2'), symmetry=1, barrier=(3.48507,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.30227,0.142819,-0.000140343,7.14342e-08,-1.42364e-11,-11508.3,51.0248], Tmin=(100,'K'), Tmax=(1243.07,'K')), NASAPolynomial(coeffs=[29.2095,0.0372343,-1.1768e-05,1.85274e-09,-1.16679e-13,-19516.5,-112.599], Tmin=(1243.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-97.9968,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH2]CC([O])(C(=C)O)C(C)CC(7153)',
    structure = SMILES('[CH2]CC([O])(C(=C)O)C(C)CC'),
    E0 = (-115.826,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,180,180,180,180,973.971,1237.53,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.159114,'amu*angstrom^2'), symmetry=1, barrier=(3.65834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159114,'amu*angstrom^2'), symmetry=1, barrier=(3.65834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159114,'amu*angstrom^2'), symmetry=1, barrier=(3.65834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159114,'amu*angstrom^2'), symmetry=1, barrier=(3.65834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159114,'amu*angstrom^2'), symmetry=1, barrier=(3.65834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159114,'amu*angstrom^2'), symmetry=1, barrier=(3.65834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159114,'amu*angstrom^2'), symmetry=1, barrier=(3.65834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159114,'amu*angstrom^2'), symmetry=1, barrier=(3.65834,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.70924,0.130478,-0.000112078,4.98031e-08,-8.81552e-12,-13674.9,49.7232], Tmin=(100,'K'), Tmax=(1362.26,'K')), NASAPolynomial(coeffs=[27.0151,0.0431992,-1.59747e-05,2.77178e-09,-1.84445e-13,-21773.4,-102.869], Tmin=(1362.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-115.826,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(623.585,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(RCCJ)"""),
)

species(
    label = 'C=C([O])C([O])(CC)C(C)CC(7154)',
    structure = SMILES('C=C([O])C([O])(CC)C(C)CC'),
    E0 = (-183.268,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.12599,0.122863,-9.93456e-05,4.23849e-08,-7.34357e-12,-21811.8,47.3382], Tmin=(100,'K'), Tmax=(1369.06,'K')), NASAPolynomial(coeffs=[22.8614,0.0498564,-1.93568e-05,3.43409e-09,-2.30864e-13,-28653.6,-81.0608], Tmin=(1369.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-183.268,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(627.743,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH]=C(O)C([O])(CC)C(C)CC(7155)',
    structure = SMILES('[CH]=C(O)C([O])(CC)C(C)CC'),
    E0 = (-73.9765,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,3615,1277.5,1000,1380,1390,370,380,2900,435,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,473.566,665.326,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156575,'amu*angstrom^2'), symmetry=1, barrier=(3.59996,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156575,'amu*angstrom^2'), symmetry=1, barrier=(3.59996,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156575,'amu*angstrom^2'), symmetry=1, barrier=(3.59996,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156575,'amu*angstrom^2'), symmetry=1, barrier=(3.59996,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156575,'amu*angstrom^2'), symmetry=1, barrier=(3.59996,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156575,'amu*angstrom^2'), symmetry=1, barrier=(3.59996,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156575,'amu*angstrom^2'), symmetry=1, barrier=(3.59996,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156575,'amu*angstrom^2'), symmetry=1, barrier=(3.59996,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.77032,0.132725,-0.000116597,5.30292e-08,-9.59855e-12,-8640.04,48.5559], Tmin=(100,'K'), Tmax=(1334.38,'K')), NASAPolynomial(coeffs=[27.0717,0.0432696,-1.60389e-05,2.78995e-09,-1.86132e-13,-16604.2,-104.023], Tmin=(1334.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-73.9765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(623.585,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]CC([CH2])C(O)(CC)C(=C)O(7156)',
    structure = SMILES('[CH2]CC([CH2])C(O)(CC)C(=C)O'),
    E0 = (-139.847,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.22518,0.140403,-0.000135317,6.76615e-08,-1.32638e-11,-16543.8,52.1332], Tmin=(100,'K'), Tmax=(1267.87,'K')), NASAPolynomial(coeffs=[28.9524,0.0374724,-1.18691e-05,1.87151e-09,-1.17927e-13,-24589.7,-110.294], Tmin=(1267.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-139.847,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(CC)C1(CC)OC[C]1O(7157)',
    structure = SMILES('[CH2]C(CC)C1(CC)OC[C]1O'),
    E0 = (-66.0821,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.97902,0.1235,-9.53018e-05,2.80352e-08,2.35654e-12,-7724.97,43.1115], Tmin=(100,'K'), Tmax=(882.561,'K')), NASAPolynomial(coeffs=[20.1252,0.0500951,-1.60531e-05,2.55007e-09,-1.62387e-13,-12669.5,-66.6751], Tmin=(882.561,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-66.0821,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(627.743,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(Isobutyl) + radical(C2CsJOH)"""),
)

species(
    label = 'CCC1CC[C](O)C1([O])CC(7158)',
    structure = SMILES('CCC1CC[C](O)C1([O])CC'),
    E0 = (-153.217,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.68072,0.107667,-5.31452e-05,-5.63225e-09,9.32645e-12,-18208.3,41.9741], Tmin=(100,'K'), Tmax=(1040.45,'K')), NASAPolynomial(coeffs=[21.6036,0.0513145,-1.97151e-05,3.58305e-09,-2.48934e-13,-24848.6,-79.9091], Tmin=(1040.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-153.217,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(636.057,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(CC(C)2OJ) + radical(C2CsJOH)"""),
)

species(
    label = 'C=C(O)C(O)(CC)C(=C)CC(7159)',
    structure = SMILES('C=C(O)C(O)(CC)C(=C)CC'),
    E0 = (-440.392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.52682,0.126397,-9.42351e-05,2.66801e-08,2.02107e-13,-52717.2,45.2255], Tmin=(100,'K'), Tmax=(1054.03,'K')), NASAPolynomial(coeffs=[26.6864,0.0437339,-1.67292e-05,3.04126e-09,-2.11619e-13,-60442,-104.679], Tmin=(1054.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-440.392,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(623.585,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
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
    label = '[CH2]C(CC)C(C)([O])C(=C)O(7160)',
    structure = SMILES('[CH2]C(CC)C(C)([O])C(=C)O'),
    E0 = (-92.2101,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,337.07,796.543,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153643,'amu*angstrom^2'), symmetry=1, barrier=(3.53256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153643,'amu*angstrom^2'), symmetry=1, barrier=(3.53256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153643,'amu*angstrom^2'), symmetry=1, barrier=(3.53256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153643,'amu*angstrom^2'), symmetry=1, barrier=(3.53256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153643,'amu*angstrom^2'), symmetry=1, barrier=(3.53256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153643,'amu*angstrom^2'), symmetry=1, barrier=(3.53256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153643,'amu*angstrom^2'), symmetry=1, barrier=(3.53256,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.69604,0.112314,-8.62121e-05,2.53836e-08,4.95704e-13,-10873.8,43.8231], Tmin=(100,'K'), Tmax=(993.19,'K')), NASAPolynomial(coeffs=[22.7469,0.0391271,-1.38224e-05,2.396e-09,-1.62646e-13,-16974.7,-80.2039], Tmin=(993.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-92.2101,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C)C([O])(CC)C(=C)O(6363)',
    structure = SMILES('[CH2]C(C)C([O])(CC)C(=C)O'),
    E0 = (-95.5573,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,337.07,796.543,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153643,'amu*angstrom^2'), symmetry=1, barrier=(3.53256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153643,'amu*angstrom^2'), symmetry=1, barrier=(3.53256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153643,'amu*angstrom^2'), symmetry=1, barrier=(3.53256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153643,'amu*angstrom^2'), symmetry=1, barrier=(3.53256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153643,'amu*angstrom^2'), symmetry=1, barrier=(3.53256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153643,'amu*angstrom^2'), symmetry=1, barrier=(3.53256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153643,'amu*angstrom^2'), symmetry=1, barrier=(3.53256,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4683.75,'J/mol'), sigma=(8.00266,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=731.59 K, Pc=20.74 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.69604,0.112314,-8.62121e-05,2.53836e-08,4.95704e-13,-11276.4,43.8231], Tmin=(100,'K'), Tmax=(993.19,'K')), NASAPolynomial(coeffs=[22.7469,0.0391271,-1.38224e-05,2.396e-09,-1.62646e-13,-17377.3,-80.2039], Tmin=(993.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-95.5573,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Isobutyl) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C(O)C([O])(CC)C[CH]CC(6337)',
    structure = SMILES('C=C(O)C([O])(CC)C[CH]CC'),
    E0 = (-124.267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.7804,0.124386,-0.000107309,5.05006e-08,-9.83902e-12,-14735.6,47.4351], Tmin=(100,'K'), Tmax=(1210.19,'K')), NASAPolynomial(coeffs=[18.6897,0.0567269,-2.34469e-05,4.30284e-09,-2.95513e-13,-19690.1,-55.2267], Tmin=(1210.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-124.267,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(623.585,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(RCCJCC)"""),
)

species(
    label = 'C=C(O)C([O])([CH]CCC)CC(7161)',
    structure = SMILES('C=C(O)C([O])([CH]CCC)CC'),
    E0 = (-118.824,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.70943,0.128694,-0.000108546,4.71945e-08,-8.16164e-12,-14033.8,50.7187], Tmin=(100,'K'), Tmax=(1394.56,'K')), NASAPolynomial(coeffs=[27.4995,0.042047,-1.53478e-05,2.64162e-09,-1.74784e-13,-22459.5,-105.069], Tmin=(1394.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-118.824,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(623.585,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCJCO) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C(O)C1(CC)OCC1CC(6348)',
    structure = SMILES('C=C(O)C1(CC)OCC1CC'),
    E0 = (-372.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.72953,0.102233,-1.13536e-05,-7.43389e-08,4.21204e-11,-44613.5,39.9127], Tmin=(100,'K'), Tmax=(907.173,'K')), NASAPolynomial(coeffs=[27.6378,0.0366171,-8.47379e-06,1.15939e-09,-7.46601e-14,-52570.1,-113.393], Tmin=(907.173,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-372.838,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(631.9,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Oxetane)"""),
)

species(
    label = '[CH2]C(CC)C([C]=C)(CC)OO(7162)',
    structure = SMILES('[CH2]C(CC)C([C]=C)(CC)OO'),
    E0 = (178.812,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.11992,0.132459,-0.000123293,6.29574e-08,-1.32056e-11,21728.3,47.5581], Tmin=(100,'K'), Tmax=(1136.91,'K')), NASAPolynomial(coeffs=[19.4606,0.056532,-2.31181e-05,4.21629e-09,-2.88756e-13,16821.3,-59.3246], Tmin=(1136.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(178.812,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(619.428,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(CC)C([O])(CC)C(C)=O(7163)',
    structure = SMILES('[CH2]C(CC)C([O])(CC)C(C)=O'),
    E0 = (-112.386,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (156.222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.60112,0.128096,-0.000127933,7.58354e-08,-1.90433e-11,-13319.4,46.0433], Tmin=(100,'K'), Tmax=(946.032,'K')), NASAPolynomial(coeffs=[13.3702,0.0647937,-2.75629e-05,5.10442e-09,-3.51748e-13,-16152.1,-25.354], Tmin=(946.032,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-112.386,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(623.585,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(CC(C)(C=O)OJ) + radical(Isobutyl)"""),
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
    label = '[CH2]C(CC)[C](CC)C(=C)O(7164)',
    structure = SMILES('[CH2]C(O)=C(CC)C([CH2])CC'),
    E0 = (-16.1559,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1600],'cm^-1')),
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
    molecularWeight = (140.223,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.97536,0.117705,-8.52051e-05,2.0027e-08,3.22362e-12,-1715.74,43.7712], Tmin=(100,'K'), Tmax=(973.699,'K')), NASAPolynomial(coeffs=[23.3307,0.0428697,-1.47846e-05,2.52967e-09,-1.70729e-13,-8024.4,-84.7308], Tmin=(973.699,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-16.1559,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(598.642,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsCsOs) + radical(C=C(O)CJ) + radical(Isobutyl)"""),
)

species(
    label = 'CH2(19)',
    structure = SMILES('[CH2]'),
    E0 = (381.563,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1032.72,2936.3,3459],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.8328,0.000224446,4.68033e-06,-6.04743e-09,2.59009e-12,45920.8,1.40666], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.16229,0.00281798,-7.56235e-07,5.05446e-11,5.65236e-15,46099.1,4.77656], Tmin=(1000,'K'), Tmax=(3000,'K'))], Tmin=(200,'K'), Tmax=(3000,'K'), E0=(381.563,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C=C(O)C([O])([CH]CC)CC(6800)',
    structure = SMILES('C=C(O)C([O])([CH]CC)CC'),
    E0 = (-95.0432,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,180,1064.56,1147.69,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.60905,0.108824,-7.81789e-05,1.95543e-08,1.52878e-12,-11216.5,44.5071], Tmin=(100,'K'), Tmax=(1033.85,'K')), NASAPolynomial(coeffs=[22.9346,0.0394958,-1.47805e-05,2.65393e-09,-1.83543e-13,-17661.2,-81.3439], Tmin=(1033.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-95.0432,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCJCO) + radical(C=CC(C)2OJ)"""),
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
    E0 = (-115.99,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-71.6544,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-53.2303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (0.502731,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-70.5822,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-65.3711,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (25.006,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-115.99,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-14.3514,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-13.3108,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-40.7097,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-39.6289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (30.8677,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-29.3815,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-62.1004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-56.1666,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-32.9831,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-14.1827,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-53.6883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-63.1079,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-73.3135,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-40.9363,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-84.1501,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (70.7256,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (8.90207,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-0.930328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-52.5902,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (327.652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (324.305,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (43.9448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (43.9448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-107.706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (285.923,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (88.1887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (226.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (286.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    products = ['butene1(127)', 'C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    products = ['[CH2]C(CC)C1(CC)OC1([CH2])O(7138)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.52e+08,'s^-1'), n=0.89, Ea=(44.3359,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H_secNd;radadd_intra] for rate rule [R4_S_D;doublebond_intra_2H_secNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    products = ['[CH2]C1(O)CC(CC)C1([O])CC(7139)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.48e+06,'s^-1'), n=1.25, Ea=(62.76,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 349 used for R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C=C(O)C([O])(CC)C(=C)CC(7140)'],
    products = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.0051739,'m^3/(mol*s)'), n=2.82163, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 102 used for Cds-CsCs_Cds-HH;HJ
Exact match found for rate rule [Cds-CsCs_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -4.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2][CH]CC(130)', 'C=C(O)C(=O)CC(4626)'],
    products = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.97e+07,'cm^3/(mol*s)'), n=1.88, Ea=(32.2168,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO-CdCs_O;YJ] for rate rule [CO-CdCs_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C2H5(29)', '[CH2]C(CC)C(=O)C(=C)O(7141)'],
    products = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7.94e+10,'cm^3/(mol*s)'), n=0, Ea=(28.0328,'kJ/mol'), T0=(1,'K'), Tmin=(333,'K'), Tmax=(363,'K'), comment="""Estimated using template [CO_O;CsJ-CsHH] for rate rule [CO-CdCs_O;CsJ-CsHH]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2COH(99)', '[CH2]C(CC)C(=O)CC(7142)'],
    products = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.16e+10,'cm^3/(mol*s)'), n=0, Ea=(48.1578,'kJ/mol'), T0=(1,'K'), Tmin=(413,'K'), Tmax=(563,'K'), comment="""Estimated using template [CO-CsCs_O;CJ] for rate rule [CO-CsCs_O;CdsJ-O2s]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['butene1(127)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.00168615,'m^3/(mol*s)'), n=2.52599, Ea=(85.3853,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 80.9 to 85.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['C2H5(29)', 'C=CC([O])(CC)C(=C)O(6516)'],
    products = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1020,'cm^3/(mol*s)'), n=2.41, Ea=(27.3634,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 418 used for Cds-CsH_Cds-HH;CsJ-CsHH
Exact match found for rate rule [Cds-CsH_Cds-HH;CsJ-CsHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    products = ['C=C(O)C([O])(CC)[C](C)CC(7143)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.265e-07,'s^-1'), n=5.639, Ea=(102.68,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 38 used for R2H_S;C_rad_out_2H;Cs_H_out_Cs2
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    products = ['[CH2][C](CC)C(O)(CC)C(=C)O(7144)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(CC)C(O)([CH]C)C(=C)O(7145)'],
    products = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(5.71,'s^-1'), n=3.021, Ea=(105.562,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 319 used for R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    products = ['C=C(O)C([O])(CC)C(C)[CH]C(7146)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.18e+10,'s^-1'), n=0.82, Ea=(146.858,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 186 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    products = ['C=C(O)C([O])([CH]C)C(C)CC(7147)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(6.44e+09,'s^-1'), n=0.13, Ea=(86.6088,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 131 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    products = ['[CH2]C([CH]C)C(O)(CC)C(=C)O(7148)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(210000,'s^-1'), n=1.76, Ea=(53.8899,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 326 used for R4H_SSS;O_rad_out;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]CC(O)(C(=C)O)C([CH2])CC(7149)'],
    products = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]CC(C)C([O])(CC)C(=C)O(7150)'],
    products = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(114000,'s^-1'), n=1.74, Ea=(82.8432,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 109 used for R4H_SSS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    products = ['[CH2]C(CC)C(O)(CC)C(=C)[O](7151)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1070.11,'s^-1'), n=2.50856, Ea=(101.808,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;O_H_out] + [R4H_SSS;O_rad_out;XH_out] + [R4H_SS(Cd)S;Y_rad_out;XH_out] for rate rule [R4H_SS(Cd)S;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C(O)C(O)(CC)C([CH2])CC(7152)'],
    products = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]CC([O])(C(=C)O)C(C)CC(7153)'],
    products = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(68850,'s^-1'), n=1.68, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 111 used for R5H_CCC;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R5H_CCC;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    products = ['C=C([O])C([O])(CC)C(C)CC(7154)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2300,'s^-1'), n=1.98, Ea=(42.6768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC(Cd);C_rad_out_2H;XH_out] for rate rule [R5H_CCC(Cd);C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C(O)C([O])(CC)C(C)CC(7155)'],
    products = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5H;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    products = ['[CH2]CC([CH2])C(O)(CC)C(=C)O(7156)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.2e+11,'s^-1'), n=0, Ea=(31.8402,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(1000,'K'), comment="""From training reaction 306 used for R5H_CCC;O_rad_out;Cs_H_out_2H
Exact match found for rate rule [R5H_CCC;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2][CH]CC(130)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    products = ['[CH2]C(CC)C1(CC)OC[C]1O(7157)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.21748e+08,'s^-1'), n=0.95, Ea=(124.892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R4_S_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    products = ['CCC1CC[C](O)C1([O])CC(7158)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.71e+11,'s^-1'), n=0.2, Ea=(115.06,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 4 used for R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    products = ['C=C(O)C(O)(CC)C(=C)CC(7159)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CH2(S)(23)', '[CH2]C(CC)C(C)([O])C(=C)O(7160)'],
    products = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction29',
    reactants = ['CH2(S)(23)', '[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    products = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    products = ['C=C(O)C([O])(CC)C[CH]CC(6337)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    products = ['C=C(O)C([O])([CH]CCC)CC(7161)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    products = ['C=C(O)C1(CC)OCC1CC(6348)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R4_SSS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(CC)C([C]=C)(CC)OO(7162)'],
    products = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(3.01978e+11,'s^-1'), n=0, Ea=(107.111,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2OOH_S;Y_rad_out] for rate rule [R2OOH_S;Cd_rad_out_double]
Euclidian distance = 2.0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    products = ['[CH2]C(CC)C([O])(CC)C(C)=O(7163)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(204.179,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction35',
    reactants = ['O(4)', '[CH2]C(CC)[C](CC)C(=C)O(7164)'],
    products = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/Cs2;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction36',
    reactants = ['CH2(19)', 'C=C(O)C([O])([CH]CC)CC(6800)'],
    products = ['[CH2]C(CC)C([O])(CC)C(=C)O(6341)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '1378',
    isomers = [
        '[CH2]C(CC)C([O])(CC)C(=C)O(6341)',
    ],
    reactants = [
        ('butene1(127)', 'C=C(O)C(=O)CC(4626)'),
        ('butene1(127)', '[CH2]C(O)=C([O])CC(4557)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '1378',
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

