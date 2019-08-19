species(
    label = 'C=C(O)C([O])(CC)C([O])O(13385)',
    structure = SMILES('C=C(O)C([O])(CC)C([O])O'),
    E0 = (-381.738,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,180,1075.76,1144,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.163243,'amu*angstrom^2'), symmetry=1, barrier=(3.75328,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163243,'amu*angstrom^2'), symmetry=1, barrier=(3.75328,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163243,'amu*angstrom^2'), symmetry=1, barrier=(3.75328,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163243,'amu*angstrom^2'), symmetry=1, barrier=(3.75328,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163243,'amu*angstrom^2'), symmetry=1, barrier=(3.75328,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163243,'amu*angstrom^2'), symmetry=1, barrier=(3.75328,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (146.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.873821,0.108296,-0.000114859,6.47338e-08,-1.47436e-11,-45737.8,41.6147], Tmin=(100,'K'), Tmax=(1058.86,'K')), NASAPolynomial(coeffs=[17.374,0.0393624,-1.72067e-05,3.25102e-09,-2.27374e-13,-49602.2,-47.4639], Tmin=(1058.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-381.738,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCOJ) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'HOCHO(59)',
    structure = SMILES('O=CO'),
    E0 = (-389.211,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1051.41,1051.94,1052.55,2787.25,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00861396,'amu*angstrom^2'), symmetry=1, barrier=(6.76824,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4005.91,'J/mol'), sigma=(3.626,'angstroms'), dipoleMoment=(1.7,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.89836,-0.00355878,3.55205e-05,-4.385e-08,1.71078e-11,-46778.6,7.34954], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.61383,0.00644964,-2.29083e-06,3.6716e-10,-2.18737e-14,-45330.3,0.847884], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-389.211,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""HOCHO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = '[CH2]C1(O)OC1(CC)C([O])O(14965)',
    structure = SMILES('[CH2]C1(O)OC1(CC)C([O])O'),
    E0 = (-335.536,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.6455,0.127577,-0.000144454,8.06713e-08,-1.69481e-11,-40100.4,40.8824], Tmin=(100,'K'), Tmax=(1318.6,'K')), NASAPolynomial(coeffs=[28.5234,0.0177705,-2.18643e-06,-3.10471e-11,1.61775e-14,-46994,-113.081], Tmin=(1318.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-335.536,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(CJC(O)2C) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C1(O)OC(O)C1([O])CC(14717)',
    structure = SMILES('[CH2]C1(O)OC(O)C1([O])CC'),
    E0 = (-343.729,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.54341,0.127922,-0.000136728,6.97847e-08,-1.30347e-11,-41036.8,42.4824], Tmin=(100,'K'), Tmax=(1569.79,'K')), NASAPolynomial(coeffs=[31.2606,0.0106832,2.58257e-06,-9.65719e-10,7.81624e-14,-48445.6,-129.916], Tmin=(1569.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-343.729,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CJC(C)OC) + radical(CC(C)2OJ)"""),
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
    label = 'C=C(O)C([O])(CC)C(=O)O(14966)',
    structure = SMILES('C=C(O)C([O])(CC)C(=O)O'),
    E0 = (-563.984,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,350,440,435,1725,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (145.133,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.00041,0.113675,-0.000137331,8.86714e-08,-2.29734e-11,-67654.7,36.0476], Tmin=(100,'K'), Tmax=(939.493,'K')), NASAPolynomial(coeffs=[16.7264,0.0382015,-1.68295e-05,3.16348e-09,-2.1972e-13,-70985.5,-48.3679], Tmin=(939.493,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-563.984,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(C=CC(C)(C=O)OJ)"""),
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
    label = 'C=C(O)C(=O)C([O])O(14967)',
    structure = SMILES('C=C(O)C(=O)C([O])O'),
    E0 = (-449.222,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,375,552.5,462.5,1710,3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,219.362,219.436,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00349823,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.490796,'amu*angstrom^2'), symmetry=1, barrier=(16.7653,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.490613,'amu*angstrom^2'), symmetry=1, barrier=(16.7651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.490721,'amu*angstrom^2'), symmetry=1, barrier=(16.7655,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.510219,0.0786043,-9.68102e-05,6.05617e-08,-1.4967e-11,-53904.8,29.2188], Tmin=(100,'K'), Tmax=(988.812,'K')), NASAPolynomial(coeffs=[14.6897,0.0212442,-9.79566e-06,1.89511e-09,-1.34213e-13,-56708.9,-39.0292], Tmin=(988.812,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-449.222,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsHH) + radical(C=OCOJ)"""),
)

species(
    label = '[O][CH]O(171)',
    structure = SMILES('[O][CH]O'),
    E0 = (5.37818,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,1989.68],'cm^-1')),
        HinderedRotor(inertia=(0.0937094,'amu*angstrom^2'), symmetry=1, barrier=(2.15456,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.40204,0.0209571,-4.96202e-05,5.95355e-08,-2.44634e-11,660.882,10.3084], Tmin=(100,'K'), Tmax=(868.739,'K')), NASAPolynomial(coeffs=[-0.343229,0.0179325,-9.40003e-06,1.81361e-09,-1.23835e-13,2076.48,32.2523], Tmin=(868.739,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.37818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsHH) + radical(OCJO) + radical(OCOJ)"""),
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
    label = 'CCC(=O)C([O])O(11574)',
    structure = SMILES('CCC(=O)C([O])O'),
    E0 = (-374.366,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,3615,1277.5,1000,1380,1390,370,380,2900,435,214.497,1592.09,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00366396,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.310634,'amu*angstrom^2'), symmetry=1, barrier=(10.1426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.310642,'amu*angstrom^2'), symmetry=1, barrier=(10.1426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.310674,'amu*angstrom^2'), symmetry=1, barrier=(10.1427,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06367,0.0702375,-9.47739e-05,7.96619e-08,-2.77969e-11,-44925.5,28.1307], Tmin=(100,'K'), Tmax=(787.846,'K')), NASAPolynomial(coeffs=[6.58624,0.0357236,-1.67338e-05,3.19334e-09,-2.21563e-13,-45594.7,4.07975], Tmin=(787.846,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-374.366,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(C=OCOJ)"""),
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
    label = 'C=C(O)C([O])(C=O)CC(6168)',
    structure = SMILES('C=C(O)C([O])(C=O)CC'),
    E0 = (-298.298,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,180,180,180,180,1600,1818.45,2694.67,3200],'cm^-1')),
        HinderedRotor(inertia=(0.159255,'amu*angstrom^2'), symmetry=1, barrier=(3.66158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159255,'amu*angstrom^2'), symmetry=1, barrier=(3.66158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159255,'amu*angstrom^2'), symmetry=1, barrier=(3.66158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159255,'amu*angstrom^2'), symmetry=1, barrier=(3.66158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159255,'amu*angstrom^2'), symmetry=1, barrier=(3.66158,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.608054,0.10466,-0.000128562,8.50037e-08,-2.25179e-11,-35713.8,33.3289], Tmin=(100,'K'), Tmax=(920.683,'K')), NASAPolynomial(coeffs=[15.433,0.0349675,-1.50146e-05,2.78299e-09,-1.91503e-13,-38667.5,-42.7339], Tmin=(920.683,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-298.298,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=CC(C)(C=O)OJ)"""),
)

species(
    label = 'C=C(O)C([O])(CC)[C](O)O(14968)',
    structure = SMILES('C=C(O)C([O])(CC)[C](O)O'),
    E0 = (-402.197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.50923,0.11555,-0.000126742,7.0472e-08,-1.53551e-11,-48169.9,43.0214], Tmin=(100,'K'), Tmax=(1124.45,'K')), NASAPolynomial(coeffs=[22.8156,0.0290186,-1.13097e-05,2.03372e-09,-1.39084e-13,-53640.3,-77.1851], Tmin=(1124.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-402.197,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cs_P) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C(O)C(O)([CH]C)C([O])O(14969)',
    structure = SMILES('C=C(O)C(O)([CH]C)C([O])O'),
    E0 = (-410.939,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3615,3650,1210,1277.5,1345,900,1000,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,180,180,180,539.477,601.155,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15934,'amu*angstrom^2'), symmetry=1, barrier=(3.66353,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15934,'amu*angstrom^2'), symmetry=1, barrier=(3.66353,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15934,'amu*angstrom^2'), symmetry=1, barrier=(3.66353,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15934,'amu*angstrom^2'), symmetry=1, barrier=(3.66353,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15934,'amu*angstrom^2'), symmetry=1, barrier=(3.66353,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15934,'amu*angstrom^2'), symmetry=1, barrier=(3.66353,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15934,'amu*angstrom^2'), symmetry=1, barrier=(3.66353,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (146.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.39385,0.113411,-0.000124017,6.92577e-08,-1.51874e-11,-49225.8,44.7558], Tmin=(100,'K'), Tmax=(1116.65,'K')), NASAPolynomial(coeffs=[22.0027,0.0296007,-1.14341e-05,2.04302e-09,-1.39095e-13,-54450.9,-70.7005], Tmin=(1116.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-410.939,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCOJ) + radical(CCJCO)"""),
)

species(
    label = 'C=C(O)C(O)(CC)[C]([O])O(14970)',
    structure = SMILES('C=C(O)C(O)(CC)[C]([O])O'),
    E0 = (-405.594,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.37925,0.117273,-0.000135564,8.09619e-08,-1.90915e-11,-48587,42.8994], Tmin=(100,'K'), Tmax=(1037.2,'K')), NASAPolynomial(coeffs=[20.5208,0.032815,-1.34222e-05,2.45468e-09,-1.6872e-13,-53130,-63.5555], Tmin=(1037.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-405.594,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cs_P)"""),
)

species(
    label = 'C=C(O)C([O])([CH]C)C(O)O(14971)',
    structure = SMILES('C=C(O)C([O])([CH]C)C(O)O'),
    E0 = (-407.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.60432,0.112572,-0.000117995,6.2018e-08,-1.26868e-11,-48805,45.1709], Tmin=(100,'K'), Tmax=(1201.17,'K')), NASAPolynomial(coeffs=[24.5853,0.025357,-9.08127e-06,1.5683e-09,-1.05189e-13,-55096.6,-85.9791], Tmin=(1201.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-407.541,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCJCO) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH2]CC(O)(C(=C)O)C([O])O(14972)',
    structure = SMILES('[CH2]CC(O)(C(=C)O)C([O])O'),
    E0 = (-405.594,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3615,3650,1210,1277.5,1345,900,1000,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,556.083,586.23,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.158443,'amu*angstrom^2'), symmetry=1, barrier=(3.64292,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158443,'amu*angstrom^2'), symmetry=1, barrier=(3.64292,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158443,'amu*angstrom^2'), symmetry=1, barrier=(3.64292,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158443,'amu*angstrom^2'), symmetry=1, barrier=(3.64292,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158443,'amu*angstrom^2'), symmetry=1, barrier=(3.64292,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158443,'amu*angstrom^2'), symmetry=1, barrier=(3.64292,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158443,'amu*angstrom^2'), symmetry=1, barrier=(3.64292,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (146.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.37925,0.117273,-0.000135564,8.09619e-08,-1.90915e-11,-48587,43.9981], Tmin=(100,'K'), Tmax=(1037.2,'K')), NASAPolynomial(coeffs=[20.5208,0.032815,-1.34222e-05,2.45468e-09,-1.6872e-13,-53130,-62.4569], Tmin=(1037.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-405.594,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCOJ) + radical(RCCJ)"""),
)

species(
    label = 'C=C(O)C(O)(CC)C([O])[O](14973)',
    structure = SMILES('C=C(O)C(O)(CC)C([O])[O]'),
    E0 = (-385.136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.859904,0.111332,-0.000128021,8.05231e-08,-2.06144e-11,-46149.9,40.5265], Tmin=(100,'K'), Tmax=(944.201,'K')), NASAPolynomial(coeffs=[15.3201,0.0427846,-1.91187e-05,3.62737e-09,-2.53485e-13,-49205.2,-36.603], Tmin=(944.201,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-385.136,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = 'C=C([O])C(O)(CC)C([O])O(14974)',
    structure = SMILES('C=C([O])C(O)(CC)C([O])O'),
    E0 = (-473.036,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.953091,0.111347,-0.000128065,7.95229e-08,-1.98799e-11,-56716.7,42.1884], Tmin=(100,'K'), Tmax=(971.594,'K')), NASAPolynomial(coeffs=[16.5713,0.0391992,-1.66778e-05,3.09255e-09,-2.13452e-13,-60122,-41.8515], Tmin=(971.594,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-473.036,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(CCOJ)"""),
)

species(
    label = '[CH]=C(O)C(O)(CC)C([O])O(14975)',
    structure = SMILES('[CH]=C(O)C(O)(CC)C([O])O'),
    E0 = (-363.745,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3615,3650,1210,1277.5,1345,900,1000,1100,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,180,1600,1721.37,2783.41,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155049,'amu*angstrom^2'), symmetry=1, barrier=(3.56487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155049,'amu*angstrom^2'), symmetry=1, barrier=(3.56487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155049,'amu*angstrom^2'), symmetry=1, barrier=(3.56487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155049,'amu*angstrom^2'), symmetry=1, barrier=(3.56487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155049,'amu*angstrom^2'), symmetry=1, barrier=(3.56487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155049,'amu*angstrom^2'), symmetry=1, barrier=(3.56487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155049,'amu*angstrom^2'), symmetry=1, barrier=(3.56487,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (146.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.49706,0.120137,-0.000142016,8.63939e-08,-2.06963e-11,-43549.7,43.0379], Tmin=(100,'K'), Tmax=(1022.72,'K')), NASAPolynomial(coeffs=[20.9669,0.0322796,-1.31595e-05,2.39941e-09,-1.64549e-13,-48144.6,-65.8427], Tmin=(1022.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-363.745,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]CC([O])(C(=C)O)C(O)O(14976)',
    structure = SMILES('[CH2]CC([O])(C(=C)O)C(O)O'),
    E0 = (-402.197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.50923,0.11555,-0.000126742,7.0472e-08,-1.53551e-11,-48169.9,44.12], Tmin=(100,'K'), Tmax=(1124.45,'K')), NASAPolynomial(coeffs=[22.8156,0.0290186,-1.13097e-05,2.03372e-09,-1.39084e-13,-53640.3,-76.0865], Tmin=(1124.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-402.197,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(RCCJ)"""),
)

species(
    label = 'C=C([O])C([O])(CC)C(O)O(14977)',
    structure = SMILES('C=C([O])C([O])(CC)C(O)O'),
    E0 = (-469.639,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.0095,0.108792,-0.000116492,6.56821e-08,-1.48005e-11,-56302.8,42.0442], Tmin=(100,'K'), Tmax=(1077.31,'K')), NASAPolynomial(coeffs=[18.6741,0.0357065,-1.47307e-05,2.70891e-09,-1.86801e-13,-60543.8,-54.3836], Tmin=(1077.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-469.639,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C(O)C([O])(CC)C(O)O(14978)',
    structure = SMILES('[CH]=C(O)C([O])(CC)C(O)O'),
    E0 = (-360.347,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.6121,0.118242,-0.000132613,7.51864e-08,-1.66708e-11,-43133.2,43.106], Tmin=(100,'K'), Tmax=(1106.66,'K')), NASAPolynomial(coeffs=[23.187,0.0286057,-1.11158e-05,1.99438e-09,-1.36214e-13,-48621.9,-79.0486], Tmin=(1106.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-360.347,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'CCC1(OC[C]1O)C([O])O(14979)',
    structure = SMILES('CCC1(OC[C]1O)C([O])O'),
    E0 = (-328.483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.705351,0.106752,-0.000115801,6.38313e-08,-1.09389e-11,-39340.6,37.0504], Tmin=(100,'K'), Tmax=(741.302,'K')), NASAPolynomial(coeffs=[14.5061,0.0408021,-1.49941e-05,2.52688e-09,-1.63554e-13,-42039.1,-34.7721], Tmin=(741.302,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-328.483,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CCOJ) + radical(C2CsJOH)"""),
)

species(
    label = 'CCC1([O])[C](O)COC1O(14696)',
    structure = SMILES('CCC1([O])[C](O)COC1O'),
    E0 = (-403.734,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.20742,0.10257,-9.97973e-05,5.22272e-08,-1.0764e-11,-48360.3,35.1279], Tmin=(100,'K'), Tmax=(1248.83,'K')), NASAPolynomial(coeffs=[19.9485,0.0309378,-9.11017e-06,1.33407e-09,-7.9158e-14,-53342.6,-70.4298], Tmin=(1248.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-403.734,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(CC(C)2OJ) + radical(C2CsJOH)"""),
)

species(
    label = 'C=C(O)C(O)(CC)C(=O)O(14980)',
    structure = SMILES('C=C(O)C(O)(CC)C(=O)O'),
    E0 = (-808.27,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.62292,0.115848,-0.000124019,6.7147e-08,-1.4222e-11,-97003.3,37.2772], Tmin=(100,'K'), Tmax=(1157.44,'K')), NASAPolynomial(coeffs=[23.6714,0.0284346,-1.07347e-05,1.89761e-09,-1.28609e-13,-102859,-88.4516], Tmin=(1157.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-808.27,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-OdCsOs) + group(Cds-CdsHH)"""),
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
    label = 'C=C(O)C(C)([O])C([O])O(14981)',
    structure = SMILES('C=C(O)C(C)([O])C([O])O'),
    E0 = (-357.958,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,180,180,180,535.122,609.206,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.158859,'amu*angstrom^2'), symmetry=1, barrier=(3.65247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158859,'amu*angstrom^2'), symmetry=1, barrier=(3.65247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158859,'amu*angstrom^2'), symmetry=1, barrier=(3.65247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158859,'amu*angstrom^2'), symmetry=1, barrier=(3.65247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158859,'amu*angstrom^2'), symmetry=1, barrier=(3.65247,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.25776,0.0938847,-0.000102356,5.84392e-08,-1.33604e-11,-42899.2,37.1569], Tmin=(100,'K'), Tmax=(1060.17,'K')), NASAPolynomial(coeffs=[16.3963,0.0310489,-1.34506e-05,2.53268e-09,-1.76913e-13,-46430.4,-44.1622], Tmin=(1060.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-357.958,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCOJ)"""),
)

species(
    label = 'C=C(O)C1(CC)OOC1O(13395)',
    structure = SMILES('C=C(O)C1(CC)OOC1O'),
    E0 = (-435.859,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.24889,0.102746,-8.68279e-05,3.00177e-08,-1.62252e-12,-52221.6,34.7418], Tmin=(100,'K'), Tmax=(1014.91,'K')), NASAPolynomial(coeffs=[23.4558,0.0281577,-1.02561e-05,1.83432e-09,-1.27586e-13,-58409.3,-90.5893], Tmin=(1014.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-435.859,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(12dioxetane)"""),
)

species(
    label = 'C=C(O)C([O])([CH]OO)CC(14982)',
    structure = SMILES('C=C(O)C([O])([CH]OO)CC'),
    E0 = (-147.125,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
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
    molecularWeight = (146.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.37785,0.117517,-0.000129668,7.42502e-08,-1.69269e-11,-17500.6,40.3128], Tmin=(100,'K'), Tmax=(1066.51,'K')), NASAPolynomial(coeffs=[20.262,0.0363564,-1.55209e-05,2.89867e-09,-2.01602e-13,-22116.5,-65.4808], Tmin=(1066.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-147.125,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCsJOOH) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C(O)C([CH][O])(CC)OO(14983)',
    structure = SMILES('C=C(O)C([CH][O])(CC)OO'),
    E0 = (-158.807,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,227.442,910.116,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.151459,'amu*angstrom^2'), symmetry=1, barrier=(3.48235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151459,'amu*angstrom^2'), symmetry=1, barrier=(3.48235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151459,'amu*angstrom^2'), symmetry=1, barrier=(3.48235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151459,'amu*angstrom^2'), symmetry=1, barrier=(3.48235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151459,'amu*angstrom^2'), symmetry=1, barrier=(3.48235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151459,'amu*angstrom^2'), symmetry=1, barrier=(3.48235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151459,'amu*angstrom^2'), symmetry=1, barrier=(3.48235,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (146.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.6006,0.126148,-0.000155888,1.00778e-07,-2.58715e-11,-18900.9,40.1655], Tmin=(100,'K'), Tmax=(952.215,'K')), NASAPolynomial(coeffs=[19.4967,0.0375238,-1.62815e-05,3.03711e-09,-2.10074e-13,-22918.7,-60.5839], Tmin=(952.215,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-158.807,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = 'C=[C]C(CC)(OO)C([O])O(14984)',
    structure = SMILES('C=[C]C(CC)(OO)C([O])O'),
    E0 = (-86.9355,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3615,1310,387.5,850,1000,180,180,180,226.977,920.171,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.151679,'amu*angstrom^2'), symmetry=1, barrier=(3.4874,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151679,'amu*angstrom^2'), symmetry=1, barrier=(3.4874,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151679,'amu*angstrom^2'), symmetry=1, barrier=(3.4874,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151679,'amu*angstrom^2'), symmetry=1, barrier=(3.4874,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151679,'amu*angstrom^2'), symmetry=1, barrier=(3.4874,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151679,'amu*angstrom^2'), symmetry=1, barrier=(3.4874,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151679,'amu*angstrom^2'), symmetry=1, barrier=(3.4874,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (146.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.974847,0.117384,-0.000150838,1.10732e-07,-3.35e-11,-10284,41.9519], Tmin=(100,'K'), Tmax=(800.405,'K')), NASAPolynomial(coeffs=[12.645,0.0493185,-2.32786e-05,4.48491e-09,-3.14369e-13,-12464.3,-20.7237], Tmin=(800.405,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-86.9355,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CCOJ)"""),
)

species(
    label = 'C=[C]C([O])(CC)C(O)OO(14985)',
    structure = SMILES('C=[C]C([O])(CC)C(O)OO'),
    E0 = (-83.538,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3615,1310,387.5,850,1000,180,180,180,203.16,936.833,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.151009,'amu*angstrom^2'), symmetry=1, barrier=(3.47198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151009,'amu*angstrom^2'), symmetry=1, barrier=(3.47198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151009,'amu*angstrom^2'), symmetry=1, barrier=(3.47198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151009,'amu*angstrom^2'), symmetry=1, barrier=(3.47198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151009,'amu*angstrom^2'), symmetry=1, barrier=(3.47198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151009,'amu*angstrom^2'), symmetry=1, barrier=(3.47198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151009,'amu*angstrom^2'), symmetry=1, barrier=(3.47198,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (146.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.791464,0.111856,-0.000128195,8.15255e-08,-2.13837e-11,-9880.26,41.6505], Tmin=(100,'K'), Tmax=(916.185,'K')), NASAPolynomial(coeffs=[14.0632,0.0470021,-2.20145e-05,4.26336e-09,-3.01207e-13,-12602.2,-28.7143], Tmin=(916.185,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-83.538,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'CCC([O])(C(C)=O)C([O])O(14986)',
    structure = SMILES('CCC([O])(C(C)=O)C([O])O'),
    E0 = (-378.134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.884537,0.118579,-0.000177295,1.55317e-07,-5.44134e-11,-45313.9,41.9405], Tmin=(100,'K'), Tmax=(823.506,'K')), NASAPolynomial(coeffs=[8.5179,0.0541261,-2.56815e-05,4.88176e-09,-3.36044e-13,-46225.6,2.27191], Tmin=(823.506,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-378.134,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(CCOJ) + radical(CC(C)(C=O)OJ)"""),
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
    label = 'C=C(O)[C](CC)C([O])O(14987)',
    structure = SMILES('[CH2]C(O)=C(CC)C([O])O'),
    E0 = (-281.598,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.567599,0.10281,-0.000115286,7.00903e-08,-1.72386e-11,-33705.8,36.5165], Tmin=(100,'K'), Tmax=(984.465,'K')), NASAPolynomial(coeffs=[15.4261,0.0378245,-1.62696e-05,3.03671e-09,-2.10503e-13,-36854.8,-40.3933], Tmin=(984.465,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-281.598,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(C=C(O)CJ) + radical(CCOJ)"""),
)

species(
    label = 'C=C(O)C([O])([CH]O)CC(6171)',
    structure = SMILES('C=C(O)C([O])([CH]O)CC'),
    E0 = (-226.914,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,504.12,637.597,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157562,'amu*angstrom^2'), symmetry=1, barrier=(3.62267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157562,'amu*angstrom^2'), symmetry=1, barrier=(3.62267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157562,'amu*angstrom^2'), symmetry=1, barrier=(3.62267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157562,'amu*angstrom^2'), symmetry=1, barrier=(3.62267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157562,'amu*angstrom^2'), symmetry=1, barrier=(3.62267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157562,'amu*angstrom^2'), symmetry=1, barrier=(3.62267,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.14686,0.108696,-0.000120647,6.88117e-08,-1.5414e-11,-27102.2,38.1158], Tmin=(100,'K'), Tmax=(1094.21,'K')), NASAPolynomial(coeffs=[20.7662,0.0285919,-1.08375e-05,1.90964e-09,-1.28799e-13,-31897.8,-69.5751], Tmin=(1094.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-226.914,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCsJOH) + radical(C=CC(C)2OJ)"""),
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
    E0 = (-381.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-334.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-289.337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-305.687,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-313.316,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-320.873,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-222.938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-381.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-175.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-223.702,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-305.377,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-306.458,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-327.848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-321.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-276.097,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-279.931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-319.436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-349.898,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-344.031,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-270.293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-179.565,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-256.846,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-319.954,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-318.338,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (61.904,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-373.454,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-54.3797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-51.6964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (20.175,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (29.0114,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-177.559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-38.5929,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (16.0904,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C(O)C([O])(CC)C([O])O(13385)'],
    products = ['HOCHO(59)', 'C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C(O)C([O])(CC)C([O])O(13385)'],
    products = ['[CH2]C1(O)OC1(CC)C([O])O(14965)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.52e+08,'s^-1'), n=0.89, Ea=(47.6831,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H_secNd;radadd_intra] for rate rule [R4_S_D;doublebond_intra_2H_secNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=C(O)C([O])(CC)C([O])O(13385)'],
    products = ['[CH2]C1(O)OC(O)C1([O])CC(14717)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.20137e+08,'s^-1'), n=0.864, Ea=(92.4015,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5_SS_D;doublebond_intra_2H;radadd_intra_O] + [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C=C(O)C([O])(CC)C(=O)O(14966)'],
    products = ['C=C(O)C([O])(CC)C([O])O(13385)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.83701,'m^3/(mol*s)'), n=1.71338, Ea=(46.5052,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;HJ] + [CO-NdNd_O;YJ] for rate rule [CO-NdNd_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C2H5(29)', 'C=C(O)C(=O)C([O])O(14967)'],
    products = ['C=C(O)C([O])(CC)C([O])O(13385)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7.94e+10,'cm^3/(mol*s)'), n=0, Ea=(28.0328,'kJ/mol'), T0=(1,'K'), Tmin=(333,'K'), Tmax=(363,'K'), comment="""Estimated using template [CO_O;CsJ-CsHH] for rate rule [CO-CdCs_O;CsJ-CsHH]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O][CH]O(171)', 'C=C(O)C(=O)CC(4626)'],
    products = ['C=C(O)C([O])(CC)C([O])O(13385)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.97e+07,'cm^3/(mol*s)'), n=1.88, Ea=(32.2168,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO-CdCs_O;YJ] for rate rule [CO-CdCs_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2COH(99)', 'CCC(=O)C([O])O(11574)'],
    products = ['C=C(O)C([O])(CC)C([O])O(13385)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.16e+10,'cm^3/(mol*s)'), n=0, Ea=(48.1578,'kJ/mol'), T0=(1,'K'), Tmin=(413,'K'), Tmax=(563,'K'), comment="""Estimated using template [CO-CsCs_O;CJ] for rate rule [CO-CsCs_O;CdsJ-O2s]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['HOCHO(59)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['C=C(O)C([O])(CC)C([O])O(13385)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6.87291,'m^3/(mol*s)'), n=1.39198, Ea=(192.416,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;CJ] + [CO-NdH_O;YJ] for rate rule [CO-NdH_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 189.3 to 192.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['OH(5)', 'C=C(O)C([O])(C=O)CC(6168)'],
    products = ['C=C(O)C([O])(CC)C([O])O(13385)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.04245,'cm^3/(mol*s)'), n=3.486, Ea=(94.7258,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [CO-CsH_O;OJ] for rate rule [CO-CsH_O;OJ_pri]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=C(O)C([O])(CC)C([O])O(13385)'],
    products = ['C=C(O)C([O])(CC)[C](O)O(14968)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.43381e+07,'s^-1'), n=1.70481, Ea=(158.036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_NonDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=C(O)C(O)([CH]C)C([O])O(14969)'],
    products = ['C=C(O)C([O])(CC)C([O])O(13385)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.71,'s^-1'), n=3.021, Ea=(105.562,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 319 used for R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=C(O)C([O])(CC)C([O])O(13385)'],
    products = ['C=C(O)C(O)(CC)[C]([O])O(14970)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C(O)C([O])(CC)C([O])O(13385)'],
    products = ['C=C(O)C([O])([CH]C)C(O)O(14971)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(210000,'s^-1'), n=1.76, Ea=(53.8899,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 326 used for R4H_SSS;O_rad_out;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]CC(O)(C(=C)O)C([O])O(14972)'],
    products = ['C=C(O)C([O])(CC)C([O])O(13385)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=C(O)C([O])(CC)C([O])O(13385)'],
    products = ['C=C(O)C(O)(CC)C([O])[O](14973)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(94.0113,'s^-1'), n=2.81534, Ea=(105.641,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;O_H_out] + [R4H_SSS;O_rad_out;XH_out] for rate rule [R4H_SSS;O_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=C(O)C([O])(CC)C([O])O(13385)'],
    products = ['C=C([O])C(O)(CC)C([O])O(14974)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1070.11,'s^-1'), n=2.50856, Ea=(101.808,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;O_H_out] + [R4H_SSS;O_rad_out;XH_out] + [R4H_SS(Cd)S;Y_rad_out;XH_out] for rate rule [R4H_SS(Cd)S;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C(O)C(O)(CC)C([O])O(14975)'],
    products = ['C=C(O)C([O])(CC)C([O])O(13385)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=C(O)C([O])(CC)C([O])O(13385)'],
    products = ['[CH2]CC([O])(C(=C)O)C(O)O(14976)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.2e+11,'s^-1'), n=0, Ea=(31.8402,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(1000,'K'), comment="""From training reaction 306 used for R5H_CCC;O_rad_out;Cs_H_out_2H
Exact match found for rate rule [R5H_CCC;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=C(O)C([O])(CC)C([O])O(13385)'],
    products = ['C=C([O])C([O])(CC)C(O)O(14977)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(22219.5,'s^-1'), n=1.84103, Ea=(37.7069,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5H_CCC;O_rad_out;XH_out] + [R5H_CCC(Cd);Y_rad_out;XH_out] for rate rule [R5H_CCC(Cd);O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C(O)C([O])(CC)C([O])O(13385)'],
    products = ['[CH]=C(O)C([O])(CC)C(O)O(14978)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.468e+06,'s^-1','*|/',3), n=1.554, Ea=(111.445,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 290 used for R5H_SSSD;O_rad_out;Cd_H_out_singleH
Exact match found for rate rule [R5H_SSSD;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O][CH]O(171)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['C=C(O)C([O])(CC)C([O])O(13385)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=C(O)C([O])(CC)C([O])O(13385)'],
    products = ['CCC1(OC[C]1O)C([O])O(14979)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.21748e+08,'s^-1'), n=0.95, Ea=(124.892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R4_S_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C(O)C([O])(CC)C([O])O(13385)'],
    products = ['CCC1([O])[C](O)COC1O(14696)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.90996e+09,'s^-1'), n=0.623333, Ea=(61.7837,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=C(O)C([O])(CC)C([O])O(13385)'],
    products = ['C=C(O)C(O)(CC)C(=O)O(14980)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CH2(S)(23)', 'C=C(O)C(C)([O])C([O])O(14981)'],
    products = ['C=C(O)C([O])(CC)C([O])O(13385)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=C(O)C([O])(CC)C([O])O(13385)'],
    products = ['C=C(O)C1(CC)OOC1O(13395)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=C(O)C([O])([CH]OO)CC(14982)'],
    products = ['C=C(O)C([O])(CC)C([O])O(13385)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.263e+10,'s^-1'), n=0, Ea=(92.7451,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnOOH;C_rad_out_H/NonDeC] for rate rule [ROOH;C_rad_out_H/NonDeC]
Euclidian distance = 1.0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=C(O)C([CH][O])(CC)OO(14983)'],
    products = ['C=C(O)C([O])(CC)C([O])O(13385)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(3.01978e+11,'s^-1'), n=0, Ea=(107.111,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OOH_S;Y_rad_out]
Euclidian distance = 0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=[C]C(CC)(OO)C([O])O(14984)'],
    products = ['C=C(O)C([O])(CC)C([O])O(13385)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(3.01978e+11,'s^-1'), n=0, Ea=(107.111,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2OOH_S;Y_rad_out] for rate rule [R2OOH_S;Cd_rad_out_double]
Euclidian distance = 2.0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=[C]C([O])(CC)C(O)OO(14985)'],
    products = ['C=C(O)C([O])(CC)C([O])O(13385)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(3.95074e+10,'s^-1'), n=0, Ea=(112.549,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OOH_SS;Y_rad_out] for rate rule [R3OOH_SS;Cd_rad_out_double]
Euclidian distance = 2.0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=C(O)C([O])(CC)C([O])O(13385)'],
    products = ['CCC([O])(C(C)=O)C([O])O(14986)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(204.179,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction32',
    reactants = ['O(4)', 'C=C(O)[C](CC)C([O])O(14987)'],
    products = ['C=C(O)C([O])(CC)C([O])O(13385)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/Cs2;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction33',
    reactants = ['O(4)', 'C=C(O)C([O])([CH]O)CC(6171)'],
    products = ['C=C(O)C([O])(CC)C([O])O(13385)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/CsO;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

network(
    label = '2751',
    isomers = [
        'C=C(O)C([O])(CC)C([O])O(13385)',
    ],
    reactants = [
        ('HOCHO(59)', 'C=C(O)C(=O)CC(4626)'),
        ('HOCHO(59)', '[CH2]C(O)=C([O])CC(4557)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '2751',
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

