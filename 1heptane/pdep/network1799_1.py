species(
    label = 'CCC(=O)[C](O)CC[C]=O(8853)',
    structure = SMILES('CCC([O])=C(O)CC[C]=O'),
    E0 = (-313.124,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,325,375,415,465,420,450,1700,1750,1855,455,950,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.31561,0.111232,-0.000120018,6.68388e-08,-1.46428e-11,-37463.8,40.5867], Tmin=(100,'K'), Tmax=(1117.59,'K')), NASAPolynomial(coeffs=[21.3204,0.0302145,-1.12782e-05,1.97291e-09,-1.32557e-13,-42523.4,-71.1356], Tmin=(1117.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-313.124,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(C=C(C)OJ) + radical(CCCJ=O)"""),
)

species(
    label = 'CH2CO(28)',
    structure = SMILES('C=C=O'),
    E0 = (-60.8183,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2120,512.5,787.5],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13241,0.0181319,-1.74093e-05,9.35336e-09,-2.01725e-12,-7148.09,13.3808], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.75871,0.00635124,-2.25955e-06,3.62322e-10,-2.15856e-14,-8085.33,-4.9649], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-60.8183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CH2CO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = 'CC[C]([O])C1(O)CCC1=O(8973)',
    structure = SMILES('CC[C]([O])C1(O)CCC1=O'),
    E0 = (-132.761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.909053,0.0989668,-8.81924e-05,4.07324e-08,-7.52729e-12,-15783.3,36.5104], Tmin=(100,'K'), Tmax=(1301.39,'K')), NASAPolynomial(coeffs=[20.2128,0.0340443,-1.33602e-05,2.39698e-09,-1.62809e-13,-21280.7,-70.9541], Tmin=(1301.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-132.761,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(CC(C)OJ) + radical(C2CsJOH)"""),
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
    label = 'CCC([O])=C(O)CC=C=O(8974)',
    structure = SMILES('CCC([O])=C(O)CC=C=O'),
    E0 = (-347.717,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.28302,0.110338,-0.000122334,6.91611e-08,-1.52921e-11,-41625.5,37.3818], Tmin=(100,'K'), Tmax=(1110.7,'K')), NASAPolynomial(coeffs=[21.8123,0.0271637,-1.00062e-05,1.73949e-09,-1.16586e-13,-46755.8,-76.4646], Tmin=(1110.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-347.717,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-Cd(CCO)HH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-(Cdd-O2d)CsH) + radical(C=C(C)OJ)"""),
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
    label = 'O=[C]CCC(O)=C=O(8975)',
    structure = SMILES('O=[C]CCC(O)=C=O'),
    E0 = (-208.907,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2120,512.5,787.5,1855,455,950,3615,1277.5,1000,350,440,435,1725,182.456,182.457],'cm^-1')),
        HinderedRotor(inertia=(0.549328,'amu*angstrom^2'), symmetry=1, barrier=(12.9771,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.549328,'amu*angstrom^2'), symmetry=1, barrier=(12.9771,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.549327,'amu*angstrom^2'), symmetry=1, barrier=(12.9771,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.549328,'amu*angstrom^2'), symmetry=1, barrier=(12.9771,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.173464,0.0840569,-0.000110011,7.26104e-08,-1.86181e-11,-24987.8,28.436], Tmin=(100,'K'), Tmax=(963.096,'K')), NASAPolynomial(coeffs=[16.1263,0.0178027,-6.82501e-06,1.18648e-09,-7.8565e-14,-28060.8,-47.9279], Tmin=(963.096,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-208.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-(Cdd-O2d)CsOs) + group(Cds-OdCsH) + radical(CCCJ=O)"""),
)

species(
    label = 'CCC([O])=C(O)C[CH]C=O(8976)',
    structure = SMILES('CCC([O])=C(O)CC=C[O]'),
    E0 = (-328.564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.34947,0.0981424,-5.42027e-05,-2.02801e-08,2.08194e-11,-39306.5,39.6278], Tmin=(100,'K'), Tmax=(932.438,'K')), NASAPolynomial(coeffs=[29.0888,0.0171212,-3.5811e-06,5.21621e-10,-3.89008e-14,-47137.1,-116.641], Tmin=(932.438,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-328.564,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CCC([O])=C(O)[CH]CC=O(8977)',
    structure = SMILES('CCC([O])=C(O)[CH]CC=O'),
    E0 = (-273.183,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.23731,0.105253,-0.000103868,5.2536e-08,-1.04577e-11,-32659.3,41.9847], Tmin=(100,'K'), Tmax=(1226.29,'K')), NASAPolynomial(coeffs=[22.1305,0.0290292,-1.063e-05,1.84696e-09,-1.23753e-13,-38390.4,-75.518], Tmin=(1226.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-273.183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(CCJCO) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C[CH]C(O)=C(O)CC[C]=O(8978)',
    structure = SMILES('C[CH]C(O)=C(O)CC[C]=O'),
    E0 = (-251.027,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,1855,455,950,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,200,800,1600],'cm^-1')),
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
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.95967,0.115524,-0.000123027,6.47989e-08,-1.3108e-11,-29963.7,43.8952], Tmin=(100,'K'), Tmax=(1265.63,'K')), NASAPolynomial(coeffs=[27.1927,0.019955,-5.69079e-06,8.48741e-10,-5.2422e-14,-37067.9,-102.529], Tmin=(1265.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-251.027,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(CCJCO) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]CC(O)=C(O)CC[C]=O(8979)',
    structure = SMILES('[CH2]CC(O)=C(O)CC[C]=O'),
    E0 = (-245.683,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,1855,455,950,3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,325,375,415,465,420,450,1700,1750,200,800,1600],'cm^-1')),
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
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.84045,0.118256,-0.000131071,7.25178e-08,-1.55218e-11,-29329.8,42.7547], Tmin=(100,'K'), Tmax=(1153.41,'K')), NASAPolynomial(coeffs=[25.4797,0.0235101,-7.85341e-06,1.29783e-09,-8.49133e-14,-35632,-92.9482], Tmin=(1153.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-245.683,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(RCCJ)"""),
)

species(
    label = 'CCC(O)=C([O])CC[C]=O(8980)',
    structure = SMILES('CCC(O)=C([O])CC[C]=O'),
    E0 = (-313.124,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,325,375,415,465,420,450,1700,1750,1855,455,950,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.31561,0.111232,-0.000120018,6.68388e-08,-1.46428e-11,-37463.8,40.5867], Tmin=(100,'K'), Tmax=(1117.59,'K')), NASAPolynomial(coeffs=[21.3204,0.0302145,-1.12782e-05,1.97291e-09,-1.32557e-13,-42523.4,-71.1356], Tmin=(1117.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-313.124,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CCC(O)=C(O)[CH]C[C]=O(8981)',
    structure = SMILES('CCC(O)=C(O)[CH]C[C]=O'),
    E0 = (-251.027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.95967,0.115524,-0.000123027,6.47989e-08,-1.3108e-11,-29963.7,43.8952], Tmin=(100,'K'), Tmax=(1265.63,'K')), NASAPolynomial(coeffs=[27.1927,0.019955,-5.69079e-06,8.48741e-10,-5.2422e-14,-37067.9,-102.529], Tmin=(1265.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-251.027,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(CCJCO) + radical(CCCJ=O)"""),
)

species(
    label = 'CCC([O])=C([O])CCC=O(8982)',
    structure = SMILES('CCC([O])=C([O])CCC=O'),
    E0 = (-335.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.629497,0.101324,-0.000101875,5.56142e-08,-1.23452e-11,-40157.7,38.811], Tmin=(100,'K'), Tmax=(1083.43,'K')), NASAPolynomial(coeffs=[16.0934,0.0395826,-1.6393e-05,3.0137e-09,-2.07492e-13,-43781.2,-43.207], Tmin=(1083.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-335.28,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(C=C(C)OJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CCC(O)=C(O)C[CH][C]=O(8983)',
    structure = SMILES('CCC(O)=C(O)C[CH][C]=O'),
    E0 = (-283.4,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,1855,455,950,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,200,800,1600],'cm^-1')),
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
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.21061,0.116972,-0.000123308,6.37983e-08,-1.25518e-11,-33844.7,42.0611], Tmin=(100,'K'), Tmax=(1361.69,'K')), NASAPolynomial(coeffs=[28.4613,0.0177701,-4.00221e-06,4.78804e-10,-2.53295e-14,-41354,-112.284], Tmin=(1361.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-283.4,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = 'C[CH]C([O])=C(O)CCC=O(8984)',
    structure = SMILES('C[CH]C([O])=C(O)CCC=O'),
    E0 = (-273.183,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.23731,0.105253,-0.000103868,5.2536e-08,-1.04577e-11,-32659.3,41.9847], Tmin=(100,'K'), Tmax=(1226.29,'K')), NASAPolynomial(coeffs=[22.1305,0.0290292,-1.063e-05,1.84696e-09,-1.23753e-13,-38390.4,-75.518], Tmin=(1226.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-273.183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(C=C(C)OJ) + radical(CCJCO)"""),
)

species(
    label = '[CH2]CC([O])=C(O)CCC=O(8985)',
    structure = SMILES('[CH2]CC([O])=C(O)CCC=O'),
    E0 = (-267.838,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2782.5,750,1395,475,1775,1000,3615,1277.5,1000,3000,3100,440,815,1455,1000,325,375,415,465,420,450,1700,1750,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.13133,0.108109,-0.000112225,6.05358e-08,-1.2955e-11,-32024.7,40.8942], Tmin=(100,'K'), Tmax=(1137.17,'K')), NASAPolynomial(coeffs=[20.2673,0.0328403,-1.29409e-05,2.33122e-09,-1.59172e-13,-36891.5,-65.0925], Tmin=(1137.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-267.838,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(C=C(C)OJ) + radical(RCCJ)"""),
)

species(
    label = 'C=[C][O](173)',
    structure = SMILES('[CH2][C]=O'),
    E0 = (160.185,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,539.612,539.669],'cm^-1')),
        HinderedRotor(inertia=(0.000578908,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39563,0.0101365,2.30741e-06,-8.97566e-09,3.68242e-12,19290.3,10.0703], Tmin=(100,'K'), Tmax=(1068.9,'K')), NASAPolynomial(coeffs=[6.35055,0.00638951,-2.69368e-06,5.4221e-10,-4.02476e-14,18240.9,-6.33602], Tmin=(1068.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CsCJ=O) + radical(CJC=O)"""),
)

species(
    label = '[CH2]C[C]=O(2129)',
    structure = SMILES('[CH2]C[C]=O'),
    E0 = (167.212,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(0.182443,'amu*angstrom^2'), symmetry=1, barrier=(4.19471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180479,'amu*angstrom^2'), symmetry=1, barrier=(4.14956,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.63265,0.0331119,-4.63428e-05,4.02599e-08,-1.40723e-11,20157.3,15.987], Tmin=(100,'K'), Tmax=(833.304,'K')), NASAPolynomial(coeffs=[4.91246,0.0168159,-7.37393e-06,1.37543e-09,-9.40246e-14,19963.2,6.51903], Tmin=(833.304,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(167.212,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCCJ=O) + radical(CJCC=O)"""),
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
    label = 'CC[C]1OC1(O)CC[C]=O(8986)',
    structure = SMILES('CC[C]1OC1(O)CC[C]=O'),
    E0 = (-166.096,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.33713,0.118917,-0.000129071,7.01265e-08,-1.43113e-11,-19730.9,41.2252], Tmin=(100,'K'), Tmax=(1379.54,'K')), NASAPolynomial(coeffs=[26.2398,0.0191467,-2.20062e-06,-6.96824e-11,1.99794e-14,-26006.3,-100.004], Tmin=(1379.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-166.096,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Ethylene_oxide) + radical(CCCJ=O) + radical(C2CsJO)"""),
)

species(
    label = 'CCC1([O])[C](O)CCC1=O(8987)',
    structure = SMILES('CCC1([O])[C](O)CCC1=O'),
    E0 = (-196.149,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.417197,0.0926637,-7.88445e-05,3.60863e-08,-6.77457e-12,-23428.6,35.3701], Tmin=(100,'K'), Tmax=(1260.94,'K')), NASAPolynomial(coeffs=[16.2064,0.0399305,-1.61146e-05,2.92115e-09,-1.99179e-13,-27621,-48.6837], Tmin=(1260.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-196.149,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(469.768,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclopentanone) + radical(CC(C)(C=O)OJ) + radical(C2CsJOH)"""),
)

species(
    label = 'CCC(O)=C(O)CC=C=O(8988)',
    structure = SMILES('CCC(O)=C(O)CC=C=O'),
    E0 = (-485.522,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.70916,0.115035,-0.000117434,5.5611e-08,-8.84651e-12,-58179.9,37.4091], Tmin=(100,'K'), Tmax=(986.423,'K')), NASAPolynomial(coeffs=[25.498,0.0238874,-7.99469e-06,1.35781e-09,-9.19618e-14,-64480.6,-98.2075], Tmin=(986.423,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-485.522,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-Cd(CCO)HH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-(Cdd-O2d)CsH)"""),
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
    label = 'CC([O])=C(O)CC[C]=O(8989)',
    structure = SMILES('CC([O])=C(O)CC[C]=O'),
    E0 = (-290.479,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,325,375,415,465,420,450,1700,1750,1855,455,950,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.601199,0.0985743,-0.000114811,6.88269e-08,-1.6172e-11,-34768.5,35.1114], Tmin=(100,'K'), Tmax=(1045.62,'K')), NASAPolynomial(coeffs=[18.5969,0.0251313,-9.45186e-06,1.65086e-09,-1.10507e-13,-38783.2,-58.3645], Tmin=(1045.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-290.479,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(C=C(C)OJ) + radical(CCCJ=O)"""),
)

species(
    label = 'C=C([O])C[C](O)C(=O)CC(8845)',
    structure = SMILES('C=C([O])CC(O)=C([O])CC'),
    E0 = (-337.987,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.16894,0.113535,-0.000115587,5.78171e-08,-1.09746e-11,-40409.5,42.8655], Tmin=(100,'K'), Tmax=(1430.46,'K')), NASAPolynomial(coeffs=[28.2173,0.0177505,-3.80482e-06,4.35562e-10,-2.23884e-14,-47996.2,-110.741], Tmin=(1430.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-337.987,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CCC1OC(=O)CCC=1O(8990)',
    structure = SMILES('CCC1OC(=O)CCC=1O'),
    E0 = (-615.595,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.161637,0.0700749,7.58316e-06,-5.80365e-08,2.57351e-11,-73870.5,28.4439], Tmin=(100,'K'), Tmax=(1031.17,'K')), NASAPolynomial(coeffs=[20.7982,0.0375699,-1.61212e-05,3.18362e-09,-2.34215e-13,-80787.6,-85.899], Tmin=(1031.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-615.595,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(469.768,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsOs) + ring(Cyclohexane)"""),
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
    label = '[CH2]CC(O)=C([O])CC(6008)',
    structure = SMILES('[CH2]CC(O)=C([O])CC'),
    E0 = (-161.259,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,325,375,415,465,420,450,1700,1750,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.288716,0.0824972,-5.7298e-05,5.85097e-09,6.62095e-12,-19229.8,33.5912], Tmin=(100,'K'), Tmax=(955.389,'K')), NASAPolynomial(coeffs=[20.4852,0.0231146,-7.38732e-06,1.25368e-09,-8.65032e-14,-24458.6,-72.2735], Tmin=(955.389,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-161.259,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(RCCJ)"""),
)

species(
    label = 'CCC(=[C]CC[C]=O)OO(8991)',
    structure = SMILES('CCC(=[C]CC[C]=O)OO'),
    E0 = (125.966,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,3615,1310,387.5,850,1000,350,440,435,1725,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,200,800,1600],'cm^-1')),
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
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.803605,0.113997,-0.000155609,1.26732e-07,-4.26973e-11,15315.3,41.7116], Tmin=(100,'K'), Tmax=(769.849,'K')), NASAPolynomial(coeffs=[10.2395,0.0516303,-2.43717e-05,4.6676e-09,-3.24707e-13,13762.8,-7.71626], Tmin=(769.849,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(125.966,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(Cds_S)"""),
)

species(
    label = 'CCC([O])C(=O)CC[C]=O(8992)',
    structure = SMILES('CCC([O])C(=O)CC[C]=O'),
    E0 = (-180.159,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.452062,0.105316,-0.000130566,9.8517e-08,-3.13545e-11,-21514.6,41.2408], Tmin=(100,'K'), Tmax=(756.765,'K')), NASAPolynomial(coeffs=[9.72076,0.0515388,-2.39575e-05,4.58814e-09,-3.20502e-13,-23054.1,-5.0006], Tmin=(756.765,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-180.159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(C=OCOJ)"""),
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
    label = 'CC[C]=C(O)CC[C]=O(8993)',
    structure = SMILES('CC[C]=C(O)CC[C]=O'),
    E0 = (1.47133,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1855,455,950,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.438362,0.0977483,-0.000102734,5.82141e-08,-1.33379e-11,336.756,36.3593], Tmin=(100,'K'), Tmax=(1054.15,'K')), NASAPolynomial(coeffs=[15.7935,0.0361552,-1.50891e-05,2.7848e-09,-1.92162e-13,-3085.36,-42.8057], Tmin=(1054.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1.47133,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = 'CCC(=O)C(=O)CC[C]=O(8994)',
    structure = SMILES('CCC(=O)C(=O)CC[C]=O'),
    E0 = (-346.507,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,1855,455,950,2750,2800,2850,1350,1500,750,1050,1375,1000,365,385,505,600,445,480,1700,1720,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.750656,0.118125,-0.000186879,1.71977e-07,-6.22806e-11,-41517.3,37.1277], Tmin=(100,'K'), Tmax=(830.077,'K')), NASAPolynomial(coeffs=[6.37381,0.0569663,-2.78826e-05,5.34655e-09,-3.689e-13,-41775.9,9.65028], Tmin=(830.077,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-346.507,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-OdCsH) + radical(CCCJ=O)"""),
)

species(
    label = 'CCC(=O)C(O)=CC[C]=O(8995)',
    structure = SMILES('CCC(=O)C(O)=CC[C]=O'),
    E0 = (-333.795,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,3615,1277.5,1000,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.200191,0.0961926,-8.96967e-05,4.29888e-08,-8.48729e-12,-39998.4,33.613], Tmin=(100,'K'), Tmax=(1187.74,'K')), NASAPolynomial(coeffs=[16.1118,0.0412578,-2.0319e-05,4.04756e-09,-2.90722e-13,-43873.2,-47.8891], Tmin=(1187.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-333.795,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCCJ=O)"""),
)

species(
    label = 'CCC(=O)C([O])CC[C]=O(8996)',
    structure = SMILES('CCC(=O)C([O])CC[C]=O'),
    E0 = (-180.159,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,375,552.5,462.5,1710,1855,455,950,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.452062,0.105316,-0.000130566,9.8517e-08,-3.13545e-11,-21514.6,41.2408], Tmin=(100,'K'), Tmax=(756.765,'K')), NASAPolynomial(coeffs=[9.72076,0.0515388,-2.39575e-05,4.58814e-09,-3.20502e-13,-23054.1,-5.0006], Tmin=(756.765,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-180.159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCCJ=O)"""),
)

species(
    label = 'CCC(=O)C(O)[CH]C[C]=O(8997)',
    structure = SMILES('CCC(=O)C(O)[CH]C[C]=O'),
    E0 = (-223.99,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,200,800,1200,1600],'cm^-1')),
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
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.298231,0.100551,-0.000108851,6.77033e-08,-1.76802e-11,-26790.1,42.4123], Tmin=(100,'K'), Tmax=(913.132,'K')), NASAPolynomial(coeffs=[11.7885,0.0476053,-2.18785e-05,4.20679e-09,-2.96141e-13,-28997.5,-14.8009], Tmin=(913.132,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-223.99,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(CCJCO)"""),
)

species(
    label = 'CCC(=O)C(O)C[CH][C]=O(8998)',
    structure = SMILES('CCC(=O)C(O)C[CH][C]=O'),
    E0 = (-256.362,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,200,800,1200,1600],'cm^-1')),
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
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.383825,0.100169,-0.000103236,5.96411e-08,-1.43117e-11,-30678.5,39.9764], Tmin=(100,'K'), Tmax=(994.284,'K')), NASAPolynomial(coeffs=[13.3448,0.0449377,-1.99109e-05,3.77047e-09,-2.63496e-13,-33408.5,-26.1773], Tmin=(994.284,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-256.362,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(CCJCHO)"""),
)

species(
    label = 'C[CH]C(=O)C(O)CC[C]=O(8999)',
    structure = SMILES('CC=C([O])C(O)CC[C]=O'),
    E0 = (-272.507,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.686369,0.106031,-0.000117314,7.11008e-08,-1.75369e-11,-32608.7,40.8823], Tmin=(100,'K'), Tmax=(978.974,'K')), NASAPolynomial(coeffs=[15.2074,0.0410895,-1.78079e-05,3.33707e-09,-2.31824e-13,-35720.6,-35.4579], Tmin=(978.974,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-272.507,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]CC(=O)C(O)CC[C]=O(9000)',
    structure = SMILES('[CH2]CC(=O)C(O)CC[C]=O'),
    E0 = (-212.303,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1855,455,950,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,375,552.5,462.5,1710,200,800,1200,1600],'cm^-1')),
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
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.788215,0.112704,-0.000143676,1.05307e-07,-3.18304e-11,-25368.4,42.3039], Tmin=(100,'K'), Tmax=(801.146,'K')), NASAPolynomial(coeffs=[12.2006,0.0478512,-2.22482e-05,4.26024e-09,-2.97608e-13,-27449.6,-17.4797], Tmin=(801.146,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-212.303,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(CJCC=O)"""),
)

species(
    label = 'CC[C]1OC(=O)CC[C]1O(9001)',
    structure = SMILES('CC[C]1OC(=O)CC[C]1O'),
    E0 = (-289.181,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.153746,0.0821025,-5.02633e-05,1.34325e-08,-1.38086e-12,-34640.5,31.6404], Tmin=(100,'K'), Tmax=(2248.5,'K')), NASAPolynomial(coeffs=[31.2252,0.0268281,-1.33894e-05,2.4997e-09,-1.65302e-13,-48613.4,-143.438], Tmin=(2248.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-289.181,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(469.768,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + ring(Cyclohexanone) + radical(C2CsJOC(O)C) + radical(C2CsJOH)"""),
)

species(
    label = 'CCC(=O)C(O)=CCC=O(9002)',
    structure = SMILES('CCC(=O)C(O)=CCC=O'),
    E0 = (-493.756,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.153925,0.0936557,-7.64698e-05,3.11818e-08,-5.23783e-12,-59238.1,32.6176], Tmin=(100,'K'), Tmax=(1367.55,'K')), NASAPolynomial(coeffs=[17.1362,0.0430832,-2.09993e-05,4.14053e-09,-2.9445e-13,-63967.1,-56.2095], Tmin=(1367.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-493.756,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-OdCsH)"""),
)

species(
    label = 'CCC(=O)C(O)CC=C=O(9003)',
    structure = SMILES('CCC(=O)C(O)CC=C=O'),
    E0 = (-458.432,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.133277,0.0952994,-8.97404e-05,4.69251e-08,-1.03483e-11,-54991.5,37.6744], Tmin=(100,'K'), Tmax=(1063.72,'K')), NASAPolynomial(coeffs=[12.7416,0.0468848,-2.14686e-05,4.1369e-09,-2.92032e-13,-57730.5,-25.2347], Tmin=(1063.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-458.432,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-(Cdd-O2d)CsH)"""),
)

species(
    label = 'CCC(=O)C(=O)CCC=O(9004)',
    structure = SMILES('CCC(=O)C(=O)CCC=O'),
    E0 = (-506.468,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.648821,0.11484,-0.000170896,1.56957e-07,-5.80497e-11,-60758.8,35.9429], Tmin=(100,'K'), Tmax=(808.076,'K')), NASAPolynomial(coeffs=[5.32966,0.0622302,-3.05158e-05,5.89595e-09,-4.10192e-13,-60973.6,13.0238], Tmin=(808.076,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-506.468,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-OdCsH)"""),
)

species(
    label = '[CH2]C(O)(C[C]=O)C(=O)CC(8854)',
    structure = SMILES('[CH2]C(O)(C[C]=O)C(=O)CC'),
    E0 = (-222.408,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1855,455,950,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,180,180,180,180,1600,1777.09,2739.24,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157068,'amu*angstrom^2'), symmetry=1, barrier=(3.6113,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157068,'amu*angstrom^2'), symmetry=1, barrier=(3.6113,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157068,'amu*angstrom^2'), symmetry=1, barrier=(3.6113,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157068,'amu*angstrom^2'), symmetry=1, barrier=(3.6113,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157068,'amu*angstrom^2'), symmetry=1, barrier=(3.6113,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157068,'amu*angstrom^2'), symmetry=1, barrier=(3.6113,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157068,'amu*angstrom^2'), symmetry=1, barrier=(3.6113,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.88299,0.13861,-0.000212102,1.72889e-07,-5.52672e-11,-26546.4,40.186], Tmin=(100,'K'), Tmax=(855.013,'K')), NASAPolynomial(coeffs=[15.7712,0.0422671,-1.89583e-05,3.48228e-09,-2.33891e-13,-29062.7,-39.2807], Tmin=(855.013,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-222.408,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(CJC(C)(C=O)O) + radical(CCCJ=O)"""),
)

species(
    label = 'CCC(O)([C]=O)CC[C]=O(9005)',
    structure = SMILES('CCC(O)([C]=O)CC[C]=O'),
    E0 = (-244.96,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,1850,1860,440,470,900,1000,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,180,1600,1766.53,2749.5,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15672,'amu*angstrom^2'), symmetry=1, barrier=(3.6033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15672,'amu*angstrom^2'), symmetry=1, barrier=(3.6033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15672,'amu*angstrom^2'), symmetry=1, barrier=(3.6033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15672,'amu*angstrom^2'), symmetry=1, barrier=(3.6033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15672,'amu*angstrom^2'), symmetry=1, barrier=(3.6033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15672,'amu*angstrom^2'), symmetry=1, barrier=(3.6033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15672,'amu*angstrom^2'), symmetry=1, barrier=(3.6033,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.35307,0.125633,-0.000180238,1.4413e-07,-4.62674e-11,-29276.7,40.3138], Tmin=(100,'K'), Tmax=(814.554,'K')), NASAPolynomial(coeffs=[13.7076,0.0453009,-2.05687e-05,3.84271e-09,-2.62482e-13,-31518.8,-27.9582], Tmin=(814.554,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-244.96,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(CC(C)(O)CJ=O)"""),
)

species(
    label = 'CCC(=O)C1(O)CCC1=O(8855)',
    structure = SMILES('CCC(=O)C1(O)CCC1=O'),
    E0 = (-489.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.497262,0.0859157,-6.20011e-05,2.25653e-08,-3.30352e-12,-58721,36.0152], Tmin=(100,'K'), Tmax=(1605.09,'K')), NASAPolynomial(coeffs=[20.3845,0.0338765,-1.33687e-05,2.36597e-09,-1.5736e-13,-65424.4,-74.608], Tmin=(1605.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-489.67,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-OdCsCs) + ring(Cyclobutanone)"""),
)

species(
    label = 'C=C(OC)[C](O)CC[C]=O(9006)',
    structure = SMILES('[CH2]C(OC)=C(O)CC[C]=O'),
    E0 = (-250.683,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
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
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.76278,0.11771,-0.000130348,7.26408e-08,-1.57143e-11,-29934.9,40.1013], Tmin=(100,'K'), Tmax=(1139.52,'K')), NASAPolynomial(coeffs=[24.5731,0.0252644,-8.65814e-06,1.44682e-09,-9.50219e-14,-35936.9,-90.3935], Tmin=(1139.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-250.683,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(C=C(O)CJ) + radical(CCCJ=O)"""),
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
    E0 = (-313.124,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-132.761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-124.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-196.635,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-46.8503,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-156.251,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-80.1611,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-93.646,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-162.63,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-199.403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-199.403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-270.447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-188.787,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-212.389,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-182.903,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-24.758,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (104.888,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-81.908,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-196.149,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-288.151,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (129.395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-67.5234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-305.593,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-254.859,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (233.077,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-170.011,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (277.827,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (244.476,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-108.147,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-121.502,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-196.103,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-26.2505,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (-102.235,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (-191.385,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (-103.948,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (-102.724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (-254.496,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (-249.724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (-249.724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (-295.342,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (-65.0894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (-87.6419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (-304.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (63.1172,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    products = ['CH2CO(28)', 'C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    products = ['CC[C]([O])C1(O)CCC1=O(8973)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.48517e+07,'s^-1'), n=1.03851, Ea=(180.363,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra;radadd_intra] for rate rule [R5_SS_D;doublebond_intra;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 178.7 to 180.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'CCC([O])=C(O)CC=C=O(8974)'],
    products = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.82e-16,'cm^3/(molecule*s)'), n=1.61, Ea=(10.992,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Ck;HJ] for rate rule [Cds-CsH_Ck;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2CO(28)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.000154228,'m^3/(mol*s)'), n=2.70143, Ea=(49.1261,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cdd;CsJ-CdHH] for rate rule [Cds-HH_Ck;CsJ-CdHH]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C2H5(29)', 'O=[C]CCC(O)=C=O(8975)'],
    products = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(8.04,'m^3/(mol*s)'), n=1.68, Ea=(54.1828,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cdd_Od;CsJ] for rate rule [Ck_O;CsJ-CsHH]
Euclidian distance = 1.41421356237
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    products = ['CCC([O])=C(O)C[CH]C=O(8976)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(791180,'s^-1'), n=2.19286, Ea=(156.873,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CCC([O])=C(O)[CH]CC=O(8977)'],
    products = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.32587e+10,'s^-1'), n=0.723333, Ea=(193.022,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_H/Cd;XH_out] for rate rule [R3H_SS_Cs;C_rad_out_H/Cd;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C[CH]C(O)=C(O)CC[C]=O(8978)'],
    products = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(187863,'s^-1'), n=2.03383, Ea=(157.381,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;C_rad_out_H/NonDeC;O_H_out] + [R3H_SS_2Cd;C_rad_out_H/NonDeC;XH_out] for rate rule [R3H_SS_2Cd;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]CC(O)=C(O)CC[C]=O(8979)'],
    products = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.00963743,'s^-1'), n=3.795, Ea=(83.0524,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;C_rad_out_2H;O_H_out] + [R4H_SS(Cd)S;C_rad_out_2H;XH_out] for rate rule [R4H_SS(Cd)S;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CCC(O)=C([O])CC[C]=O(8980)'],
    products = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.11e+06,'s^-1','*|/',3), n=1.78, Ea=(113.721,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4H_SDS;O_rad_out;XH_out] for rate rule [R4H_SDS;O_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    products = ['CCC(O)=C(O)[CH]C[C]=O(8981)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.22e+06,'s^-1','*|/',3), n=1.78, Ea=(113.721,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 289 used for R4H_SDS;O_rad_out;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SDS;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    products = ['CCC([O])=C([O])CCC=O(8982)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2300,'s^-1'), n=1.98, Ea=(42.6768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC(Cd);Y_rad_out;XH_out] for rate rule [R5H_CCC(Cd);CO_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CCC(O)=C(O)C[CH][C]=O(8983)'],
    products = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(722272,'s^-1'), n=1.6737, Ea=(94.6126,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSMS;Y_rad_out;XH_out] for rate rule [R5H_SSMS;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C[CH]C([O])=C(O)CCC=O(8984)'],
    products = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(46.1,'s^-1'), n=3.21, Ea=(60.7935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H;C_rad_out_H/NonDeC;XH_out] for rate rule [R6H_SMSSR;C_rad_out_H/NonDeC;CO_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]CC([O])=C(O)CCC=O(8985)'],
    products = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(342721,'s^-1'), n=1.49389, Ea=(84.9352,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7H;C_rad_out_2H;XH_out] for rate rule [R7H;C_rad_out_2H;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=[C][O](173)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.56662e+07,'m^3/(mol*s)'), n=0.074875, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [C_rad/H2/Cd;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C[C]=O(2129)', 'CCC([O])=[C]O(4619)'],
    products = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.83556e+08,'m^3/(mol*s)'), n=-0.424942, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;C_rad/H2/Cs] + [Cd_rad;C_pri_rad] for rate rule [Cd_rad;C_rad/H2/Cs]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    products = ['CC[C]1OC1(O)CC[C]=O(8986)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd_NdNd;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    products = ['CCC1([O])[C](O)CCC1=O(8987)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.31572e+09,'s^-1'), n=0.656667, Ea=(116.975,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_secNd;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_secNd;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 114.6 to 117.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction20',
    reactants = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    products = ['CCC(O)=C(O)CC=C=O(8988)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CH2(S)(23)', 'CC([O])=C(O)CC[C]=O(8989)'],
    products = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.87e+13,'cm^3/(mol*s)','*|/',0.25), n=-0.146, Ea=(0.0118826,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 6 used for carbene;C_pri/Cd
Exact match found for rate rule [carbene;C_pri/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    products = ['C=C([O])C[C](O)C(=O)CC(8845)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    products = ['CCC1OC(=O)CCC=1O(8990)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Ypri_rad_out] for rate rule [R6_SSSDS;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CO(12)', '[CH2]CC(O)=C([O])CC(6008)'],
    products = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1230.59,'m^3/(mol*s)'), n=1.0523, Ea=(25.6182,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [COm;C_rad/H2/Cs]
Euclidian distance = 0
family: R_Addition_COm"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CCC(=[C]CC[C]=O)OO(8991)'],
    products = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.01978e+11,'s^-1'), n=0, Ea=(107.111,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2OOH;Y_rad_out] for rate rule [R2OOH_D;Cd_rad_out_ND]
Euclidian distance = 2.2360679775
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    products = ['CCC([O])C(=O)CC[C]=O(8992)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[C]=O(361)', '[CH2]CC(O)=C([O])CC(6008)'],
    products = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['O(4)', 'CC[C]=C(O)CC[C]=O(8993)'],
    products = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction29',
    reactants = ['H(3)', 'CCC(=O)C(=O)CC[C]=O(8994)'],
    products = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(10.8137,'m^3/(mol*s)'), n=2.04, Ea=(26.5684,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Od_CO-DeNd;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction30',
    reactants = ['H(3)', 'CCC(=O)C(O)=CC[C]=O(8995)'],
    products = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(25.3411,'m^3/(mol*s)'), n=1.79231, Ea=(0.501101,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsH_Cds-OneDe;HJ] for rate rule [Cds-CsH_Cds-COOs;HJ]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=[C][O](173)', 'C=C(O)C(=O)CC(4626)'],
    products = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(0.0105038,'m^3/(mol*s)'), n=2.40545, Ea=(2.17933,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-HH_Cds-OneDe;CJ] for rate rule [Cds-HH_Cds-COOs;CJ]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction32',
    reactants = ['CCC(=O)C([O])CC[C]=O(8996)'],
    products = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.70223e+09,'s^-1'), n=1.15155, Ea=(153.908,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_OneDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_CO]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['CCC(=O)C(O)[CH]C[C]=O(8997)'],
    products = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(8.83e+10,'s^-1'), n=0.3, Ea=(121.754,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_OneDe] for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_CO]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['CCC(=O)C(O)C[CH][C]=O(8998)'],
    products = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(2.3012e-19,'s^-1'), n=9.03667, Ea=(64.9775,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Y_rad_out;Cs_H_out_OneDe] for rate rule [R3H_SS_Cs;Y_rad_out;Cs_H_out_CO]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C[CH]C(=O)C(O)CC[C]=O(8999)'],
    products = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(7.66994e+09,'s^-1'), n=0.768667, Ea=(168.559,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_H/NonDeC;Cs_H_out_NonDe] for rate rule [R3H_SS;C_rad_out_H/NonDeC;Cs_H_out_NDMustO]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]CC(=O)C(O)CC[C]=O(9000)'],
    products = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.86e+10,'s^-1'), n=0.58, Ea=(109.579,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_2H;Cs_H_out_NonDe] for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_NDMustO]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    products = ['CC[C]1OC(=O)CC[C]1O(9001)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(4.39512e+11,'s^-1'), n=0.277081, Ea=(58.6284,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6plus;carbonyl_intra;radadd_intra_CO] + [R6_linear;multiplebond_intra;radadd_intra] for rate rule [R6_linear;carbonyl_intra_Nd;radadd_intra_CO]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction38',
    reactants = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    products = ['CCC(=O)C(O)=CCC=O(9002)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction39',
    reactants = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    products = ['CCC(=O)C(O)CC=C=O(9003)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction40',
    reactants = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    products = ['CCC(=O)C(=O)CCC=O(9004)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.07e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_De] for rate rule [R5radEndo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]C(O)(C[C]=O)C(=O)CC(8854)'],
    products = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction42',
    reactants = ['CCC(O)([C]=O)CC[C]=O(9005)'],
    products = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction43',
    reactants = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    products = ['CCC(=O)C1(O)CCC1=O(8855)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_OneDe/O;Ypri_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction44',
    reactants = ['C=C(OC)[C](O)CC[C]=O(9006)'],
    products = ['CCC(=O)[C](O)CC[C]=O(8853)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

network(
    label = '1799',
    isomers = [
        'CCC(=O)[C](O)CC[C]=O(8853)',
    ],
    reactants = [
        ('CH2CO(28)', 'C=C(O)C(=O)CC(4626)'),
        ('CH2CO(28)', '[CH2]C(O)=C([O])CC(4557)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '1799',
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

