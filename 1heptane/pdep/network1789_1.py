species(
    label = 'C=C([O])C([O])(CC)C(=C)O(8843)',
    structure = SMILES('C=C([O])C([O])(CC)C(=C)O'),
    E0 = (-226.343,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,325,375,415,465,420,450,1700,1750,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,180,1067.68,1145.06,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.164189,'amu*angstrom^2'), symmetry=1, barrier=(3.77502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164189,'amu*angstrom^2'), symmetry=1, barrier=(3.77502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164189,'amu*angstrom^2'), symmetry=1, barrier=(3.77502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164189,'amu*angstrom^2'), symmetry=1, barrier=(3.77502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164189,'amu*angstrom^2'), symmetry=1, barrier=(3.77502,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.27039,0.1028,-8.76054e-05,3.21232e-08,-2.89632e-12,-27021.6,39.8182], Tmin=(100,'K'), Tmax=(1050.58,'K')), NASAPolynomial(coeffs=[23.6019,0.0285773,-1.08689e-05,1.98135e-09,-1.38613e-13,-33377.7,-86.7815], Tmin=(1050.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-226.343,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CC(C)2OJ)"""),
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
    label = '[CH2]C1(O)OC1(CC)C(=C)[O](9616)',
    structure = SMILES('[CH2]C1(O)OC1(CC)C(=C)[O]'),
    E0 = (-181.327,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.98361,0.1303,-0.000139679,7.10851e-08,-1.31029e-11,-21482.1,44.6874], Tmin=(100,'K'), Tmax=(1621.06,'K')), NASAPolynomial(coeffs=[31.0325,0.00922409,4.43806e-06,-1.37774e-09,1.07198e-13,-28279,-127.109], Tmin=(1621.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-181.327,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CJC(O)2C) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C1([O])OC1(CC)C(=C)O(9597)',
    structure = SMILES('[CH2]C1([O])OC1(CC)C(=C)O'),
    E0 = (-87.1414,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.93112,0.122311,-0.00012784,6.56092e-08,-1.25052e-11,-10205.3,40.6655], Tmin=(100,'K'), Tmax=(1509.63,'K')), NASAPolynomial(coeffs=[28.6685,0.0163524,-4.68808e-07,-3.93871e-10,4.05628e-14,-17212.9,-116.409], Tmin=(1509.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-87.1414,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CC(C)(O)OJ) + radical(CJC(O)2C)"""),
)

species(
    label = '[CH2]C1(O)OC(=C)C1([O])CC(9297)',
    structure = SMILES('[CH2]C1(O)OC(=C)C1([O])CC'),
    E0 = (-109.622,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.68042,0.0957367,-1.62082e-05,-8.44731e-08,5.10013e-11,-12952.4,33.3488], Tmin=(100,'K'), Tmax=(897.673,'K')), NASAPolynomial(coeffs=[37.5035,0.00159783,6.64338e-06,-1.59097e-09,1.09864e-13,-23229.2,-169.518], Tmin=(897.673,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-109.622,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(CJC(C)OC) + radical(C=CC(C)2OJ)"""),
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
    label = 'C=C([O])C(=O)C(=C)O(9617)',
    structure = SMILES('[CH2]C(=O)C(=O)C(=C)O'),
    E0 = (-240.684,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,3615,1277.5,1000,3000,3100,440,815,1455,1000,365,385,505,600,445,480,1700,1720,263.446,263.811],'cm^-1')),
        HinderedRotor(inertia=(0.00240092,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.408454,'amu*angstrom^2'), symmetry=1, barrier=(20.2151,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.409306,'amu*angstrom^2'), symmetry=1, barrier=(20.2075,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.411838,'amu*angstrom^2'), symmetry=1, barrier=(20.2117,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.225253,0.0801935,-9.26942e-05,5.31482e-08,-1.18754e-11,-28809.2,27.7474], Tmin=(100,'K'), Tmax=(1099.04,'K')), NASAPolynomial(coeffs=[17.414,0.0176348,-7.31286e-06,1.35711e-09,-9.45168e-14,-32587.5,-56.8019], Tmin=(1099.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-240.684,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)HHH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsHH) + radical(CJCC=O)"""),
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
    label = '[CH2]C(=O)C(=O)CC(4620)',
    structure = SMILES('[CH2]C(=O)C(=O)CC'),
    E0 = (-162.882,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,365,385,505,600,445,480,1700,1720,352.088,2238.1],'cm^-1')),
        HinderedRotor(inertia=(0.0859868,'amu*angstrom^2'), symmetry=1, barrier=(7.56458,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0859642,'amu*angstrom^2'), symmetry=1, barrier=(7.56459,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0859615,'amu*angstrom^2'), symmetry=1, barrier=(7.56392,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.37635,'amu*angstrom^2'), symmetry=1, barrier=(33.1076,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15233,0.0680038,-7.99912e-05,5.95539e-08,-1.92307e-11,-19492.6,25.214], Tmin=(100,'K'), Tmax=(738.439,'K')), NASAPolynomial(coeffs=[6.798,0.0374228,-1.78728e-05,3.47419e-09,-2.45165e-13,-20326.4,-0.311287], Tmin=(738.439,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-162.882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)Cs) + radical(CJCC=O)"""),
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
    label = 'C=C([O])C(O)([CH]C)C(=C)O(9618)',
    structure = SMILES('C=C([O])C(O)([CH]C)C(=C)O'),
    E0 = (-255.543,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,325,375,415,465,420,450,1700,1750,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,545.63,587.683,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.159338,'amu*angstrom^2'), symmetry=1, barrier=(3.6635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159338,'amu*angstrom^2'), symmetry=1, barrier=(3.6635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159338,'amu*angstrom^2'), symmetry=1, barrier=(3.6635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159338,'amu*angstrom^2'), symmetry=1, barrier=(3.6635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159338,'amu*angstrom^2'), symmetry=1, barrier=(3.6635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159338,'amu*angstrom^2'), symmetry=1, barrier=(3.6635,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.36897,0.102886,-7.89782e-05,1.34318e-08,6.58124e-12,-30527.7,41.4522], Tmin=(100,'K'), Tmax=(960.245,'K')), NASAPolynomial(coeffs=[26.4848,0.0217662,-6.79037e-06,1.1722e-09,-8.33311e-14,-37486.4,-100.176], Tmin=(960.245,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-255.543,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(CCJCO)"""),
)

species(
    label = '[CH]=C(O)C([O])(CC)C(=C)O(9619)',
    structure = SMILES('[CH]=C(O)C([O])(CC)C(=C)O'),
    E0 = (-117.051,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,180,180,180,180,1600,1602.56,2883.22,3200],'cm^-1')),
        HinderedRotor(inertia=(0.149965,'amu*angstrom^2'), symmetry=1, barrier=(3.44799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149965,'amu*angstrom^2'), symmetry=1, barrier=(3.44799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149965,'amu*angstrom^2'), symmetry=1, barrier=(3.44799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149965,'amu*angstrom^2'), symmetry=1, barrier=(3.44799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149965,'amu*angstrom^2'), symmetry=1, barrier=(3.44799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149965,'amu*angstrom^2'), symmetry=1, barrier=(3.44799,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.65503,0.109672,-9.47166e-05,3.00374e-08,1.05864e-13,-13861.4,40.0992], Tmin=(100,'K'), Tmax=(995.43,'K')), NASAPolynomial(coeffs=[27.1287,0.023134,-8.20176e-06,1.48931e-09,-1.06392e-13,-21034.8,-105.882], Tmin=(995.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-117.051,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C(O)C([O])([CH]C)C(=C)O(9620)',
    structure = SMILES('C=C(O)C([O])([CH]C)C(=C)O'),
    E0 = (-164.245,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,325,375,415,465,420,450,1700,1750,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,432.059,692.006,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155236,'amu*angstrom^2'), symmetry=1, barrier=(3.56917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155236,'amu*angstrom^2'), symmetry=1, barrier=(3.56917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155236,'amu*angstrom^2'), symmetry=1, barrier=(3.56917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155236,'amu*angstrom^2'), symmetry=1, barrier=(3.56917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155236,'amu*angstrom^2'), symmetry=1, barrier=(3.56917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155236,'amu*angstrom^2'), symmetry=1, barrier=(3.56917,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.36876,0.10073,-6.8752e-05,2.35936e-09,1.01682e-11,-19545.3,40.4717], Tmin=(100,'K'), Tmax=(973.903,'K')), NASAPolynomial(coeffs=[27.2817,0.0219584,-7.34451e-06,1.33818e-09,-9.79758e-14,-26970.7,-106.464], Tmin=(973.903,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-164.245,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCJCO) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH2]CC(O)(C(=C)[O])C(=C)O(9621)',
    structure = SMILES('[CH2]CC(O)(C(=C)[O])C(=C)O'),
    E0 = (-250.199,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,325,375,415,465,420,450,1700,1750,180,180,180,524.075,611.672,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157874,'amu*angstrom^2'), symmetry=1, barrier=(3.62982,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157874,'amu*angstrom^2'), symmetry=1, barrier=(3.62982,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157874,'amu*angstrom^2'), symmetry=1, barrier=(3.62982,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157874,'amu*angstrom^2'), symmetry=1, barrier=(3.62982,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157874,'amu*angstrom^2'), symmetry=1, barrier=(3.62982,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157874,'amu*angstrom^2'), symmetry=1, barrier=(3.62982,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.02343,0.114555,-0.000117366,5.91569e-08,-1.14536e-11,-29859.9,43.0993], Tmin=(100,'K'), Tmax=(1327.37,'K')), NASAPolynomial(coeffs=[28.0619,0.0199075,-5.90539e-06,9.14239e-10,-5.79777e-14,-37495.7,-109.243], Tmin=(1327.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-250.199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(RCCJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])C(O)(CC)C(=C)[O](9622)',
    structure = SMILES('C=C([O])C(O)(CC)C(=C)[O]'),
    E0 = (-317.641,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.46082,0.107143,-0.000105193,5.22925e-08,-1.01604e-11,-37995.8,40.0982], Tmin=(100,'K'), Tmax=(1260.91,'K')), NASAPolynomial(coeffs=[23.9604,0.0265,-9.25963e-06,1.57157e-09,-1.0408e-13,-44406.6,-88.4382], Tmin=(1260.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-317.641,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C(O)C(O)(CC)C(=C)[O](9623)',
    structure = SMILES('[CH]=C(O)C(O)(CC)C(=C)[O]'),
    E0 = (-208.349,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,180,180,180,180,1600,1694.47,2803.25,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153709,'amu*angstrom^2'), symmetry=1, barrier=(3.53407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153709,'amu*angstrom^2'), symmetry=1, barrier=(3.53407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153709,'amu*angstrom^2'), symmetry=1, barrier=(3.53407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153709,'amu*angstrom^2'), symmetry=1, barrier=(3.53407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153709,'amu*angstrom^2'), symmetry=1, barrier=(3.53407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153709,'amu*angstrom^2'), symmetry=1, barrier=(3.53407,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.09715,0.116936,-0.000122291,6.28222e-08,-1.23884e-11,-24824.5,41.9784], Tmin=(100,'K'), Tmax=(1296.96,'K')), NASAPolynomial(coeffs=[28.3301,0.0196515,-5.79435e-06,8.93154e-10,-5.65394e-14,-32427.5,-111.611], Tmin=(1296.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-208.349,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])C(O)(CC)C(=C)O(9624)',
    structure = SMILES('[CH]=C([O])C(O)(CC)C(=C)O'),
    E0 = (-208.349,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,180,180,180,180,1600,1694.47,2803.25,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153709,'amu*angstrom^2'), symmetry=1, barrier=(3.53407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153709,'amu*angstrom^2'), symmetry=1, barrier=(3.53407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153709,'amu*angstrom^2'), symmetry=1, barrier=(3.53407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153709,'amu*angstrom^2'), symmetry=1, barrier=(3.53407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153709,'amu*angstrom^2'), symmetry=1, barrier=(3.53407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153709,'amu*angstrom^2'), symmetry=1, barrier=(3.53407,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.09715,0.116936,-0.000122291,6.28222e-08,-1.23884e-11,-24824.5,41.9784], Tmin=(100,'K'), Tmax=(1296.96,'K')), NASAPolynomial(coeffs=[28.3301,0.0196515,-5.79435e-06,8.93154e-10,-5.65394e-14,-32427.5,-111.611], Tmin=(1296.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-208.349,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]CC([O])(C(=C)O)C(=C)O(9625)',
    structure = SMILES('[CH2]CC([O])(C(=C)O)C(=C)O'),
    E0 = (-158.901,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,325,375,415,465,420,450,1700,1750,180,180,180,395.279,730.412,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154415,'amu*angstrom^2'), symmetry=1, barrier=(3.55031,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154415,'amu*angstrom^2'), symmetry=1, barrier=(3.55031,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154415,'amu*angstrom^2'), symmetry=1, barrier=(3.55031,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154415,'amu*angstrom^2'), symmetry=1, barrier=(3.55031,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154415,'amu*angstrom^2'), symmetry=1, barrier=(3.55031,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154415,'amu*angstrom^2'), symmetry=1, barrier=(3.55031,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.54057,0.10684,-8.83518e-05,2.46877e-08,1.68512e-12,-18898.6,40.3788], Tmin=(100,'K'), Tmax=(996.419,'K')), NASAPolynomial(coeffs=[26.6608,0.0237089,-8.48837e-06,1.55043e-09,-1.11061e-13,-26011.9,-103.068], Tmin=(996.419,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-158.901,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(RCCJ) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C(O)C([O])(CC)[C]1CO1(9626)',
    structure = SMILES('C=C(O)C([O])(CC)[C]1CO1'),
    E0 = (-73.9671,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.67758,0.111298,-0.00011727,6.47195e-08,-1.37529e-11,-8680.11,40.5155], Tmin=(100,'K'), Tmax=(1275.1,'K')), NASAPolynomial(coeffs=[22.1807,0.026727,-6.33914e-06,7.38119e-10,-3.54512e-14,-13973.7,-77.2845], Tmin=(1275.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-73.9671,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(C2CsJO) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C([O])C1(CC)OC[C]1O(9627)',
    structure = SMILES('C=C([O])C1(CC)OC[C]1O'),
    E0 = (-174.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.77027,0.107148,-0.000107194,5.58831e-08,-1.10774e-11,-20735.6,39.8153], Tmin=(100,'K'), Tmax=(1418.57,'K')), NASAPolynomial(coeffs=[22.3119,0.0246162,-4.4591e-06,3.34006e-10,-6.83395e-15,-26096.4,-79.6003], Tmin=(1418.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-174.274,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Oxetane) + radical(C=C(C)OJ) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2][C](O)C1(CC)OCC1=O(9007)',
    structure = SMILES('[CH2][C](O)C1(CC)OCC1=O'),
    E0 = (-90.2878,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.18656,0.10783,-0.000106116,5.35459e-08,-1.07367e-11,-10667.1,36.0584], Tmin=(100,'K'), Tmax=(1207.83,'K')), NASAPolynomial(coeffs=[21.2922,0.0333852,-1.36609e-05,2.51369e-09,-1.73684e-13,-16097.1,-76.6325], Tmin=(1207.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-90.2878,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutane) + radical(CJCO) + radical(C2CsJOH)"""),
)

species(
    label = 'C=C1OC[C](O)C1([O])CC(9354)',
    structure = SMILES('C=C1OC[C](O)C1([O])CC'),
    E0 = (-161.9,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.537705,0.0799695,-1.73053e-05,-4.44115e-08,2.54545e-11,-19290.8,35.7453], Tmin=(100,'K'), Tmax=(958.284,'K')), NASAPolynomial(coeffs=[23.7289,0.0264972,-8.45681e-06,1.50605e-09,-1.09637e-13,-26137.3,-91.7489], Tmin=(958.284,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-161.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(469.768,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(C=CC(C)2OJ) + radical(C2CsJOH)"""),
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
    label = 'C=C([O])C(C)([O])C(=C)O(9628)',
    structure = SMILES('C=C([O])C(C)([O])C(=C)O'),
    E0 = (-202.562,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,529.646,606.383,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.159015,'amu*angstrom^2'), symmetry=1, barrier=(3.65606,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159015,'amu*angstrom^2'), symmetry=1, barrier=(3.65606,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159015,'amu*angstrom^2'), symmetry=1, barrier=(3.65606,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159015,'amu*angstrom^2'), symmetry=1, barrier=(3.65606,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.494516,0.0865026,-6.8532e-05,1.73993e-08,2.02287e-12,-24189.9,34.7875], Tmin=(100,'K'), Tmax=(1000.91,'K')), NASAPolynomial(coeffs=[21.9231,0.0214408,-7.7852e-06,1.42075e-09,-1.01166e-13,-29906.1,-79.5224], Tmin=(1000.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-202.562,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C(O)C1(CC)OOC1=C(9485)',
    structure = SMILES('C=C(O)C1(CC)OOC1=C'),
    E0 = (-140.814,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.570867,0.0770441,-2.53517e-06,-5.96323e-08,2.97477e-11,-16750.5,34.4034], Tmin=(100,'K'), Tmax=(987.2,'K')), NASAPolynomial(coeffs=[25.3135,0.0271508,-1.02746e-05,2.01615e-09,-1.5258e-13,-24540.5,-103.711], Tmin=(987.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-140.814,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutane)"""),
)

species(
    label = 'C=[C]C(CC)(OO)C(=C)[O](9629)',
    structure = SMILES('C=[C]C(CC)(OO)C(=C)[O]'),
    E0 = (68.46,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3615,1310,387.5,850,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,275.536,858.892,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15179,'amu*angstrom^2'), symmetry=1, barrier=(3.48996,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15179,'amu*angstrom^2'), symmetry=1, barrier=(3.48996,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15179,'amu*angstrom^2'), symmetry=1, barrier=(3.48996,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15179,'amu*angstrom^2'), symmetry=1, barrier=(3.48996,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15179,'amu*angstrom^2'), symmetry=1, barrier=(3.48996,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15179,'amu*angstrom^2'), symmetry=1, barrier=(3.48996,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.89301,0.106354,-0.000104606,5.36812e-08,-1.11285e-11,8411.11,38.431], Tmin=(100,'K'), Tmax=(1156.7,'K')), NASAPolynomial(coeffs=[18.6213,0.0388716,-1.70955e-05,3.24445e-09,-2.27553e-13,3896.65,-58.5551], Tmin=(1156.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(68.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C([O])(CC)C(=C)OO(9630)',
    structure = SMILES('C=[C]C([O])(CC)C(=C)OO'),
    E0 = (212.748,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3615,1310,387.5,850,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,289.596,851.282,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152317,'amu*angstrom^2'), symmetry=1, barrier=(3.50206,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152317,'amu*angstrom^2'), symmetry=1, barrier=(3.50206,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152317,'amu*angstrom^2'), symmetry=1, barrier=(3.50206,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152317,'amu*angstrom^2'), symmetry=1, barrier=(3.50206,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152317,'amu*angstrom^2'), symmetry=1, barrier=(3.50206,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152317,'amu*angstrom^2'), symmetry=1, barrier=(3.50206,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0714265,0.0971948,-9.24923e-05,4.96655e-08,-1.14909e-11,25727.9,38.4954], Tmin=(100,'K'), Tmax=(1003.33,'K')), NASAPolynomial(coeffs=[11.2143,0.052203,-2.52303e-05,4.97421e-09,-3.55401e-13,23463.2,-15.9894], Tmin=(1003.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(212.748,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C([O])C([O])(CC)C(C)=O(9631)',
    structure = SMILES('C=C([O])C([O])(CC)C(C)=O'),
    E0 = (-215.364,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.67477,0.112483,-0.000135215,8.05256e-08,-1.32612e-11,-25742.8,36.7104], Tmin=(100,'K'), Tmax=(633.693,'K')), NASAPolynomial(coeffs=[12.5379,0.0466262,-2.08543e-05,3.90288e-09,-2.6839e-13,-27769.6,-23.785], Tmin=(633.693,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-215.364,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)(C=O)OJ) + radical(C=C(C)OJ)"""),
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
    label = 'C=C([O])[C](CC)C(=C)O(9632)',
    structure = SMILES('[CH2]C(O)=C(CC)C(=C)[O]'),
    E0 = (-134.07,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,325,350,375,415,440,465,420,435,450,1700,1725,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.90116,0.111026,-0.000117768,6.16102e-08,-1.22243e-11,-15896.4,35.9856], Tmin=(100,'K'), Tmax=(1367.16,'K')), NASAPolynomial(coeffs=[26.7136,0.0170337,-3.37366e-06,3.32883e-10,-1.4186e-14,-22760.7,-107.502], Tmin=(1367.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-134.07,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(C=C(O)CJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=[C]C([O])(CC)C(=C)O(9633)',
    structure = SMILES('C=[C]C([O])(CC)C(=C)O'),
    E0 = (88.2529,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,420.019,717.559,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155273,'amu*angstrom^2'), symmetry=1, barrier=(3.57004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155273,'amu*angstrom^2'), symmetry=1, barrier=(3.57004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155273,'amu*angstrom^2'), symmetry=1, barrier=(3.57004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155273,'amu*angstrom^2'), symmetry=1, barrier=(3.57004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155273,'amu*angstrom^2'), symmetry=1, barrier=(3.57004,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.628126,0.092044,-7.95451e-05,3.48078e-08,-6.07079e-12,10789,36.4352], Tmin=(100,'K'), Tmax=(1374.89,'K')), NASAPolynomial(coeffs=[20.5059,0.0305585,-1.24648e-05,2.28149e-09,-1.56468e-13,4977.61,-72.2528], Tmin=(1374.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(88.2529,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH2]C1(O)CC(=O)C1([O])CC(9037)',
    structure = SMILES('[CH2]C1(O)CC(=O)C1([O])CC'),
    E0 = (-103.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.848824,0.103598,-0.000102798,5.42637e-08,-1.15305e-11,-12299.2,36.3886], Tmin=(100,'K'), Tmax=(1135.75,'K')), NASAPolynomial(coeffs=[18.0241,0.0371284,-1.50091e-05,2.73212e-09,-1.87284e-13,-16586.1,-57.0645], Tmin=(1135.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-103.735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(CJC(C)2O) + radical(CC(C)(C=O)OJ)"""),
)

species(
    label = 'C=C(O)C([O])([CH]C)C(C)=O(9634)',
    structure = SMILES('C=C(O)C([O])([CH]C)C(C)=O'),
    E0 = (-153.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.0931,0.115012,-0.000137396,8.80849e-08,-2.2636e-11,-18252.8,39.1586], Tmin=(100,'K'), Tmax=(947.831,'K')), NASAPolynomial(coeffs=[17.0065,0.0386293,-1.6518e-05,3.06466e-09,-2.11362e-13,-21683.9,-47.192], Tmin=(947.831,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-153.266,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-CdsCsOs) + group(Cds-OdCsCs) + group(Cds-CdsHH) + radical(CCJCO) + radical(C=CC(C)(C=O)OJ)"""),
)

species(
    label = '[CH2]CC([O])(C(=C)O)C(C)=O(9635)',
    structure = SMILES('[CH2]CC([O])(C(=C)O)C(C)=O'),
    E0 = (-147.922,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,378.685,769.526,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.164237,'amu*angstrom^2'), symmetry=1, barrier=(3.77613,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164237,'amu*angstrom^2'), symmetry=1, barrier=(3.77613,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164237,'amu*angstrom^2'), symmetry=1, barrier=(3.77613,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164237,'amu*angstrom^2'), symmetry=1, barrier=(3.77613,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164237,'amu*angstrom^2'), symmetry=1, barrier=(3.77613,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164237,'amu*angstrom^2'), symmetry=1, barrier=(3.77613,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.26226,0.121073,-0.000156768,1.10107e-07,-3.10175e-11,-17606.2,39.0575], Tmin=(100,'K'), Tmax=(866.971,'K')), NASAPolynomial(coeffs=[16.1305,0.0408257,-1.79245e-05,3.33999e-09,-2.29745e-13,-20621.9,-42.3696], Tmin=(866.971,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-147.922,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-CdsCsOs) + group(Cds-OdCsCs) + group(Cds-CdsHH) + radical(C=CC(C)(C=O)OJ) + radical(RCCJ)"""),
)

species(
    label = '[CH]=C(O)C([O])(CC)C(C)=O(9636)',
    structure = SMILES('[CH]=C(O)C([O])(CC)C(C)=O'),
    E0 = (-106.072,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1863.6,2650.18,3200],'cm^-1')),
        HinderedRotor(inertia=(0.160436,'amu*angstrom^2'), symmetry=1, barrier=(3.68874,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160436,'amu*angstrom^2'), symmetry=1, barrier=(3.68874,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160436,'amu*angstrom^2'), symmetry=1, barrier=(3.68874,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160436,'amu*angstrom^2'), symmetry=1, barrier=(3.68874,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160436,'amu*angstrom^2'), symmetry=1, barrier=(3.68874,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160436,'amu*angstrom^2'), symmetry=1, barrier=(3.68874,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.40026,0.124196,-0.000164211,1.16943e-07,-3.32706e-11,-12568,38.1684], Tmin=(100,'K'), Tmax=(859.976,'K')), NASAPolynomial(coeffs=[16.6777,0.0401135,-1.75578e-05,3.25976e-09,-2.2348e-13,-15677.4,-46.321], Tmin=(859.976,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-106.072,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CC(C)(C=O)OJ)"""),
)

species(
    label = '[CH2][C](O)C1(CC)OOC1=C(9458)',
    structure = SMILES('[CH2][C](O)C1(CC)OOC1=C'),
    E0 = (167.915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.28222,0.102658,-9.13752e-05,4.09378e-08,-7.24113e-12,20397,39.6178], Tmin=(100,'K'), Tmax=(1368.23,'K')), NASAPolynomial(coeffs=[23.601,0.0299125,-1.16243e-05,2.07958e-09,-1.4107e-13,13587.7,-88.231], Tmin=(1368.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(167.915,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(CJCO) + radical(C2CsJOH)"""),
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
    label = 'C=C(O)C([O])(CC)C[C]=O(8851)',
    structure = SMILES('C=C(O)C([O])(CC)C[C]=O'),
    E0 = (-194.004,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,559.641,584.928,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.159624,'amu*angstrom^2'), symmetry=1, barrier=(3.67006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159624,'amu*angstrom^2'), symmetry=1, barrier=(3.67006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159624,'amu*angstrom^2'), symmetry=1, barrier=(3.67006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159624,'amu*angstrom^2'), symmetry=1, barrier=(3.67006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159624,'amu*angstrom^2'), symmetry=1, barrier=(3.67006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159624,'amu*angstrom^2'), symmetry=1, barrier=(3.67006,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4868.81,'J/mol'), sigma=(7.87367,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=760.50 K, Pc=22.63 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.09161,0.109079,-0.000115091,6.35948e-08,-1.39994e-11,-23147.4,40.7623], Tmin=(100,'K'), Tmax=(1103.94,'K')), NASAPolynomial(coeffs=[19.4402,0.0346838,-1.40053e-05,2.54948e-09,-1.74915e-13,-27680.6,-60.3221], Tmin=(1103.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-194.004,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCCJ=O)"""),
)

species(
    label = 'C=C(O)C1(CC)OCC1=O(8862)',
    structure = SMILES('C=C(O)C1(CC)OCC1=O'),
    E0 = (-393.945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.04615,0.0926444,-4.81736e-05,-1.54823e-08,1.55951e-11,-47182.6,32.9995], Tmin=(100,'K'), Tmax=(974.152,'K')), NASAPolynomial(coeffs=[25.6358,0.0249969,-8.54674e-06,1.56485e-09,-1.14232e-13,-54369.7,-105.234], Tmin=(974.152,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-393.945,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclobutane)"""),
)

species(
    label = '[CH2]C(OCC)=C([O])C(=C)O(9637)',
    structure = SMILES('[CH2]C(OCC)=C([O])C(=C)O'),
    E0 = (-306.591,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,325,350,375,415,440,465,420,435,450,1700,1725,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.80466,0.130923,-0.000149723,8.25595e-08,-1.71653e-11,-36613.3,38.419], Tmin=(100,'K'), Tmax=(1297.71,'K')), NASAPolynomial(coeffs=[31.3458,0.0144576,-2.15496e-06,9.82895e-11,2.07924e-15,-44533.5,-131.603], Tmin=(1297.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-306.591,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]C(OC(=C)O)=C([O])CC(9638)',
    structure = SMILES('[CH2]C(OC(=C)O)=C([O])CC'),
    E0 = (-205.623,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,325,350,375,415,440,465,420,435,450,1700,1725,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.3844,0.115851,-0.000131711,7.78879e-08,-1.81667e-11,-24534.5,39.1258], Tmin=(100,'K'), Tmax=(1049.9,'K')), NASAPolynomial(coeffs=[20.5061,0.0324507,-1.25552e-05,2.22633e-09,-1.50292e-13,-29131.1,-67.5491], Tmin=(1049.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-205.623,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(C=C(O)CJ) + radical(C=C(C)OJ)"""),
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
    label = 'C=C(O)C([O])([C]=O)CC(9639)',
    structure = SMILES('C=C(O)C([O])([C]=O)CC'),
    E0 = (-141.286,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,180,1600,1777.87,2737.22,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157467,'amu*angstrom^2'), symmetry=1, barrier=(3.62047,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157467,'amu*angstrom^2'), symmetry=1, barrier=(3.62047,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157467,'amu*angstrom^2'), symmetry=1, barrier=(3.62047,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157467,'amu*angstrom^2'), symmetry=1, barrier=(3.62047,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157467,'amu*angstrom^2'), symmetry=1, barrier=(3.62047,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.530013,0.107824,-0.000145397,9.90387e-08,-2.37928e-11,-16837.2,33.614], Tmin=(100,'K'), Tmax=(668.051,'K')), NASAPolynomial(coeffs=[14.0045,0.0355226,-1.61173e-05,3.02008e-09,-2.07201e-13,-19107.7,-33.1025], Tmin=(668.051,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-141.286,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=CC(C)(C=O)OJ) + radical(C=CC(C)(O)CJ=O)"""),
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
    E0 = (-226.343,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-180.075,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-85.8745,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-109.452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-104.777,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-1.38918,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-166.066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-202.423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-149.981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (75.0775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-80.147,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-166.519,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-124.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-164.041,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-164.041,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-116.224,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-24.758,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (4.87354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-101.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-89.3355,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-161.9,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (217.3,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-140.814,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (175.571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (325.297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-22.1637,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (108.934,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (331.258,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-103.379,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-139.734,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-95.2036,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-73.032,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (167.915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (-111.283,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (51.597,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (-218.058,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (7.20906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (108.177,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (240.277,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    products = ['CH2CO(28)', 'C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    products = ['[CH2]C1(O)OC1(CC)C(=C)[O](9616)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.52e+08,'s^-1'), n=0.89, Ea=(46.2682,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H_secNd;radadd_intra] for rate rule [R4_S_D;doublebond_intra_2H_secNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    products = ['[CH2]C1([O])OC1(CC)C(=C)O(9597)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.62453e+10,'s^-1'), n=0.5217, Ea=(140.468,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H;radadd_intra] for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    products = ['[CH2]C1(O)OC(=C)C1([O])CC(9297)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.20137e+08,'s^-1'), n=0.864, Ea=(116.891,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5_SS_D;doublebond_intra_2H;radadd_intra_O] + [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C2H5(29)', 'C=C([O])C(=O)C(=C)O(9617)'],
    products = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7.94e+10,'cm^3/(mol*s)'), n=0, Ea=(28.0328,'kJ/mol'), T0=(1,'K'), Tmin=(333,'K'), Tmax=(363,'K'), comment="""Estimated using template [CO_O;CsJ-CsHH] for rate rule [CO-DeDe_O;CsJ-CsHH]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2COH(99)', '[CH2]C(=O)C(=O)CC(4620)'],
    products = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2339.95,'m^3/(mol*s)'), n=0.573452, Ea=(58.2237,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;CJ] for rate rule [CO-CdCs_O;CdsJ-O2s]
Euclidian distance = 2.82842712475
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=[C][O](173)', 'C=C(O)C(=O)CC(4626)'],
    products = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.97e+07,'cm^3/(mol*s)'), n=1.88, Ea=(32.2168,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO-CdCs_O;YJ] for rate rule [CO-CdCs_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH2CO(28)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3525.68,'m^3/(mol*s)'), n=1.07323, Ea=(43.3386,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cdd_Od;CJ] + [Ck_O;YJ] for rate rule [Ck_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=C([O])C(O)([CH]C)C(=C)O(9618)'],
    products = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.71,'s^-1'), n=3.021, Ea=(105.562,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 319 used for R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=C(O)C([O])(CC)C(=C)O(9619)'],
    products = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=C(O)C([O])([CH]C)C(=C)O(9620)'],
    products = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5920,'s^-1'), n=2.11, Ea=(84.0984,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/NonDeC;O_H_out] for rate rule [R4H_SS(Cd)S;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]CC(O)(C(=C)[O])C(=C)O(9621)'],
    products = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    products = ['C=C([O])C(O)(CC)C(=C)[O](9622)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1070.11,'s^-1'), n=2.50856, Ea=(101.808,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;O_H_out] + [R4H_SSS;O_rad_out;XH_out] + [R4H_SS(Cd)S;Y_rad_out;XH_out] for rate rule [R4H_SS(Cd)S;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=C(O)C(O)(CC)C(=C)[O](9623)'],
    products = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=C([O])C(O)(CC)C(=C)O(9624)'],
    products = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]CC([O])(C(=C)O)C(=C)O(9625)'],
    products = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4600,'s^-1'), n=1.98, Ea=(42.6768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC(Cd);C_rad_out_2H;XH_out] for rate rule [R5H_CCC(Cd);C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=[C][O](173)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    products = ['C=C(O)C([O])(CC)[C]1CO1(9626)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    products = ['C=C([O])C1(CC)OC[C]1O(9627)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.21748e+08,'s^-1'), n=0.95, Ea=(124.892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R4_S_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    products = ['[CH2][C](O)C1(CC)OCC1=O(9007)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(137.007,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    products = ['C=C1OC[C](O)C1([O])CC(9354)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.90996e+09,'s^-1'), n=0.623333, Ea=(64.4423,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 61.8 to 64.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction22',
    reactants = ['CH2(S)(23)', 'C=C([O])C(C)([O])C(=C)O(9628)'],
    products = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    products = ['C=C(O)C1(CC)OOC1=C(9485)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(85.5288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination
Ea raised from 82.5 to 85.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=[C]C(CC)(OO)C(=C)[O](9629)'],
    products = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.01978e+11,'s^-1'), n=0, Ea=(107.111,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2OOH_S;Y_rad_out] for rate rule [R2OOH_S;Cd_rad_out_double]
Euclidian distance = 2.0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=[C]C([O])(CC)C(=C)OO(9630)'],
    products = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.95074e+10,'s^-1'), n=0, Ea=(112.549,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OOH_SS;Y_rad_out] for rate rule [R3OOH_SS;Cd_rad_out_double]
Euclidian distance = 2.0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    products = ['C=C([O])C([O])(CC)C(C)=O(9631)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(204.179,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction27',
    reactants = ['O(4)', 'C=C([O])[C](CC)C(=C)O(9632)'],
    products = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/Cs;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['O(4)', 'C=[C]C([O])(CC)C(=C)O(9633)'],
    products = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    products = ['[CH2]C1(O)CC(=O)C1([O])CC(9037)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(6.48e+06,'s^-1'), n=1.25, Ea=(122.964,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 349 used for R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    products = ['C=C(O)C([O])([CH]C)C(C)=O(9634)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(6.44e+09,'s^-1'), n=0.13, Ea=(86.6088,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 131 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]CC([O])(C(=C)O)C(C)=O(9635)'],
    products = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(68850,'s^-1'), n=1.68, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_2H;Cs_H_out_2H] for rate rule [R5H_CC(O2d)CC;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]=C(O)C([O])(CC)C(C)=O(9636)'],
    products = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5H;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    products = ['[CH2][C](O)C1(CC)OOC1=C(9458)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(394.258,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 394.0 to 394.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    products = ['CCC1([O])[C](O)CCC1=O(8987)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.71e+11,'s^-1'), n=0.2, Ea=(115.06,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 4 used for R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=C(O)C([O])(CC)C[C]=O(8851)'],
    products = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    products = ['C=C(O)C1(CC)OCC1=O(8862)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R4_SSS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C(OCC)=C([O])C(=C)O(9637)'],
    products = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C(OC(=C)O)=C([O])CC(9638)'],
    products = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction39',
    reactants = ['CH2(19)', 'C=C(O)C([O])([C]=O)CC(9639)'],
    products = ['C=C([O])C([O])(CC)C(=C)O(8843)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '1789',
    isomers = [
        'C=C([O])C([O])(CC)C(=C)O(8843)',
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
    label = '1789',
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

