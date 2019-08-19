species(
    label = '[CH2]C(O)C([CH2])(O)C(=O)CC(6869)',
    structure = SMILES('[CH2]C(O)C([CH2])(O)C(=O)CC'),
    E0 = (-260.45,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1705.44,2804.31,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154818,'amu*angstrom^2'), symmetry=1, barrier=(3.55956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154818,'amu*angstrom^2'), symmetry=1, barrier=(3.55956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154818,'amu*angstrom^2'), symmetry=1, barrier=(3.55956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154818,'amu*angstrom^2'), symmetry=1, barrier=(3.55956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154818,'amu*angstrom^2'), symmetry=1, barrier=(3.55956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154818,'amu*angstrom^2'), symmetry=1, barrier=(3.55956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154818,'amu*angstrom^2'), symmetry=1, barrier=(3.55956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154818,'amu*angstrom^2'), symmetry=1, barrier=(3.55956,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.35196,0.14577,-0.00020325,1.49875e-07,-4.37062e-11,-31101.8,42.0175], Tmin=(100,'K'), Tmax=(842.814,'K')), NASAPolynomial(coeffs=[19.5929,0.0416236,-1.79014e-05,3.26972e-09,-2.21019e-13,-34801,-60.1018], Tmin=(842.814,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-260.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJCO) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = 'CH2CHOH(42)',
    structure = SMILES('C=CO'),
    E0 = (-138.725,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3615,1277.5,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(1.72808,'amu*angstrom^2'), symmetry=1, barrier=(39.7321,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.11,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28758,0.0197013,1.96383e-06,-1.9439e-08,1.02617e-11,-16537.3,14.1333], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[7.49818,0.0103957,-3.66891e-06,5.85206e-10,-3.47374e-14,-18164.3,-13.8388], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-138.725,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CH2CHOH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = '[CH2]C(O)C1(O)CC1([O])CC(12800)',
    structure = SMILES('[CH2]C(O)C1(O)CC1([O])CC'),
    E0 = (-165.922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.07192,0.12027,-0.000120423,6.11667e-08,-1.21419e-11,-19725.9,40.2982], Tmin=(100,'K'), Tmax=(1235.08,'K')), NASAPolynomial(coeffs=[25.9097,0.0296463,-1.03599e-05,1.75711e-09,-1.16342e-13,-26637.8,-100.605], Tmin=(1235.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-165.922,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CJCO) + radical(CC(C)2OJ)"""),
)

species(
    label = '[CH2]C1(O)C(O)CC1([O])CC(12801)',
    structure = SMILES('[CH2]C1(O)C(O)CC1([O])CC'),
    E0 = (-169.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.93739,0.116512,-0.000110333,5.31472e-08,-1.00627e-11,-20210.3,38.3451], Tmin=(100,'K'), Tmax=(1287.94,'K')), NASAPolynomial(coeffs=[25.3857,0.0316549,-1.15049e-05,1.9919e-09,-1.33127e-13,-27248.5,-100.387], Tmin=(1287.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-169.914,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CC(C)2OJ) + radical(CJC(C)2O)"""),
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
    label = '[CH2]C(O)(C(=C)O)C(=O)CC(12802)',
    structure = SMILES('[CH2]C(O)(C(=C)O)C(=O)CC'),
    E0 = (-385.868,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1713.58,2794.83,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15488,'amu*angstrom^2'), symmetry=1, barrier=(3.56099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15488,'amu*angstrom^2'), symmetry=1, barrier=(3.56099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15488,'amu*angstrom^2'), symmetry=1, barrier=(3.56099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15488,'amu*angstrom^2'), symmetry=1, barrier=(3.56099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15488,'amu*angstrom^2'), symmetry=1, barrier=(3.56099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15488,'amu*angstrom^2'), symmetry=1, barrier=(3.56099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15488,'amu*angstrom^2'), symmetry=1, barrier=(3.56099,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.16,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.02186,0.136388,-0.000180575,1.24801e-07,-3.40231e-11,-46195.8,38.1431], Tmin=(100,'K'), Tmax=(901.029,'K')), NASAPolynomial(coeffs=[20.2733,0.0374116,-1.58019e-05,2.88614e-09,-1.96226e-13,-50213.5,-67.0942], Tmin=(901.029,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-385.868,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(O)(C=O)CJ)"""),
)

species(
    label = '[CH2][CH]O(284)',
    structure = SMILES('[CH2][CH]O'),
    E0 = (143.484,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.217851,'amu*angstrom^2'), symmetry=1, barrier=(5.00882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0197382,'amu*angstrom^2'), symmetry=1, barrier=(14.867,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.84769,0.0229973,-1.85068e-05,7.77211e-09,-1.30725e-12,17300.5,13.0245], Tmin=(100,'K'), Tmax=(1418.34,'K')), NASAPolynomial(coeffs=[7.97636,0.00853337,-3.21015e-06,5.82141e-10,-3.99268e-14,15845.7,-13.5107], Tmin=(1418.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(143.484,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CJCO) + radical(CCsJOH)"""),
)

species(
    label = 'C2H5CO(71)',
    structure = SMILES('CC[C]=O'),
    E0 = (-44.1874,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(0.192008,'amu*angstrom^2'), symmetry=1, barrier=(4.41465,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191563,'amu*angstrom^2'), symmetry=1, barrier=(4.40441,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0712,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3133.67,'J/mol'), sigma=(5.35118,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=489.47 K, Pc=46.4 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.91313,0.0225012,-8.59328e-06,4.80319e-10,2.48743e-13,-5274.24,14.655], Tmin=(100,'K'), Tmax=(1673.18,'K')), NASAPolynomial(coeffs=[7.46306,0.0158512,-6.42141e-06,1.12499e-09,-7.32055e-14,-7388.54,-11.406], Tmin=(1673.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-44.1874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""CH3CH2CO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C(O)C(=C)O(5812)',
    structure = SMILES('[CH2]C(O)C(=C)O'),
    E0 = (-192.298,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,380.141],'cm^-1')),
        HinderedRotor(inertia=(0.137777,'amu*angstrom^2'), symmetry=1, barrier=(14.0659,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139215,'amu*angstrom^2'), symmetry=1, barrier=(14.0622,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137815,'amu*angstrom^2'), symmetry=1, barrier=(14.0638,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136883,'amu*angstrom^2'), symmetry=1, barrier=(14.0691,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.0972,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.809015,0.0586654,-3.40006e-05,-1.04469e-08,1.19044e-11,-23002.6,25.4896], Tmin=(100,'K'), Tmax=(931.884,'K')), NASAPolynomial(coeffs=[19.1312,0.00933009,-1.76755e-06,2.45497e-10,-1.87827e-14,-27690.1,-68.4401], Tmin=(931.884,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-192.298,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CJCO)"""),
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
    label = '[CH2]C(O)C(=C)C(=O)CC(12803)',
    structure = SMILES('[CH2]C(O)C(=C)C(=O)CC'),
    E0 = (-179.117,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.10158,0.0932452,-9.20851e-05,5.39415e-08,-1.37874e-11,-21408.8,34.927], Tmin=(100,'K'), Tmax=(912.645,'K')), NASAPolynomial(coeffs=[9.42597,0.0523776,-2.4916e-05,4.87584e-09,-3.46871e-13,-23110.7,-9.20543], Tmin=(912.645,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-179.117,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsHH) + radical(CJCO)"""),
)

species(
    label = '[CH2]C(O)(C=C)C(=O)CC(6377)',
    structure = SMILES('[CH2]C(O)(C=C)C(=O)CC'),
    E0 = (-171.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1706.04,2810.13,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154561,'amu*angstrom^2'), symmetry=1, barrier=(3.55366,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154561,'amu*angstrom^2'), symmetry=1, barrier=(3.55366,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154561,'amu*angstrom^2'), symmetry=1, barrier=(3.55366,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154561,'amu*angstrom^2'), symmetry=1, barrier=(3.55366,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154561,'amu*angstrom^2'), symmetry=1, barrier=(3.55366,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154561,'amu*angstrom^2'), symmetry=1, barrier=(3.55366,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.161,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.749509,0.113763,-0.000140295,8.89869e-08,-1.85042e-11,-20441.3,33.559], Tmin=(100,'K'), Tmax=(647.718,'K')), NASAPolynomial(coeffs=[12.7354,0.0463508,-2.09198e-05,3.93392e-09,-2.71362e-13,-22521,-28.2104], Tmin=(647.718,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-171.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(O)(C=O)CJ)"""),
)

species(
    label = '[CH2]C(O)([C](C)O)C(=O)CC(6866)',
    structure = SMILES('[CH2]C(O)([C](C)O)C(=O)CC'),
    E0 = (-295.412,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.43779,0.150406,-0.000224877,1.79374e-07,-5.63098e-11,-35306.4,41.1097], Tmin=(100,'K'), Tmax=(851.937,'K')), NASAPolynomial(coeffs=[17.6282,0.0450777,-1.98575e-05,3.62676e-09,-2.43221e-13,-38322,-50.1146], Tmin=(851.937,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-295.412,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJC(C)(C=O)O) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C(O)C(C)([O])C(=O)CC(12804)',
    structure = SMILES('[CH2]C(O)C(C)([O])C(=O)CC'),
    E0 = (-230.079,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,315.595,834.782,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.164482,'amu*angstrom^2'), symmetry=1, barrier=(3.78177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164482,'amu*angstrom^2'), symmetry=1, barrier=(3.78177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164482,'amu*angstrom^2'), symmetry=1, barrier=(3.78177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164482,'amu*angstrom^2'), symmetry=1, barrier=(3.78177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164482,'amu*angstrom^2'), symmetry=1, barrier=(3.78177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164482,'amu*angstrom^2'), symmetry=1, barrier=(3.78177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164482,'amu*angstrom^2'), symmetry=1, barrier=(3.78177,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.62487,0.132968,-0.000190542,1.55678e-07,-5.10546e-11,-27478.4,42.1994], Tmin=(100,'K'), Tmax=(829.88,'K')), NASAPolynomial(coeffs=[12.6191,0.0524671,-2.36289e-05,4.39259e-09,-2.9886e-13,-29434.7,-21.4058], Tmin=(829.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-230.079,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CC(C)(C=O)OJ) + radical(CJCO)"""),
)

species(
    label = '[CH2][C](O)C(C)(O)C(=O)CC(12805)',
    structure = SMILES('[CH2][C](O)C(C)(O)C(=O)CC'),
    E0 = (-295.434,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.14975,0.144607,-0.000214977,1.73983e-07,-5.56776e-11,-35320,41.0923], Tmin=(100,'K'), Tmax=(846.493,'K')), NASAPolynomial(coeffs=[15.6418,0.047942,-2.1371e-05,3.93217e-09,-2.65095e-13,-37880.9,-39.1115], Tmin=(846.493,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-295.434,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJCO) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C(O)(C(=O)CC)C(C)[O](6328)',
    structure = SMILES('[CH2]C(O)(C(=O)CC)C(C)[O]'),
    E0 = (-241.679,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,566.073,582.003,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15914,'amu*angstrom^2'), symmetry=1, barrier=(3.65895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15914,'amu*angstrom^2'), symmetry=1, barrier=(3.65895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15914,'amu*angstrom^2'), symmetry=1, barrier=(3.65895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15914,'amu*angstrom^2'), symmetry=1, barrier=(3.65895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15914,'amu*angstrom^2'), symmetry=1, barrier=(3.65895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15914,'amu*angstrom^2'), symmetry=1, barrier=(3.65895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15914,'amu*angstrom^2'), symmetry=1, barrier=(3.65895,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4929.37,'J/mol'), sigma=(8.10058,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=769.96 K, Pc=21.04 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.90639,0.134962,-0.000171496,1.17148e-07,-3.20109e-11,-28859,39.5985], Tmin=(100,'K'), Tmax=(893.955,'K')), NASAPolynomial(coeffs=[18.3767,0.0442012,-1.91969e-05,3.56532e-09,-2.45143e-13,-32485.2,-55.9808], Tmin=(893.955,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-241.679,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CC(C)OJ) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = '[CH2]C([O])(C(=O)CC)C(C)O(6870)',
    structure = SMILES('[CH2]C([O])(C(=O)CC)C(C)O'),
    E0 = (-230.057,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,379.326,771.708,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.163177,'amu*angstrom^2'), symmetry=1, barrier=(3.75177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163177,'amu*angstrom^2'), symmetry=1, barrier=(3.75177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163177,'amu*angstrom^2'), symmetry=1, barrier=(3.75177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163177,'amu*angstrom^2'), symmetry=1, barrier=(3.75177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163177,'amu*angstrom^2'), symmetry=1, barrier=(3.75177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163177,'amu*angstrom^2'), symmetry=1, barrier=(3.75177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163177,'amu*angstrom^2'), symmetry=1, barrier=(3.75177,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.90963,0.138724,-0.000200266,1.60796e-07,-5.15458e-11,-27465,42.2054], Tmin=(100,'K'), Tmax=(833.394,'K')), NASAPolynomial(coeffs=[14.6024,0.0496085,-2.21188e-05,4.08804e-09,-2.77059e-13,-29874.7,-32.3916], Tmin=(833.394,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-230.057,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJC(C)(C=O)O) + radical(CC(C)(C=O)OJ)"""),
)

species(
    label = '[CH2]C([O])C(C)(O)C(=O)CC(6871)',
    structure = SMILES('[CH2]C([O])C(C)(O)C(=O)CC'),
    E0 = (-241.701,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.59665,0.128874,-0.000160422,1.09978e-07,-3.04938e-11,-28873.5,39.5053], Tmin=(100,'K'), Tmax=(877.496,'K')), NASAPolynomial(coeffs=[16.3233,0.047187,-2.07838e-05,3.88866e-09,-2.68543e-13,-32018.4,-44.6063], Tmin=(877.496,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-241.701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJCO) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C(O)C(C)(O)C(=O)[CH]C(12806)',
    structure = SMILES('[CH2]C(O)C(C)(O)C([O])=CC'),
    E0 = (-335.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.77419,0.121837,-0.00013311,7.5929e-08,-1.70762e-11,-40149.7,43.1987], Tmin=(100,'K'), Tmax=(1088.82,'K')), NASAPolynomial(coeffs=[22.0208,0.0344218,-1.26836e-05,2.19384e-09,-1.4613e-13,-45331.4,-73.6231], Tmin=(1088.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-335.59,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C(O)(C(=O)[CH]C)C(C)O(6873)',
    structure = SMILES('[CH2]C(O)(C([O])=CC)C(C)O'),
    E0 = (-333.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.63024,0.121226,-0.000133891,7.82001e-08,-1.81314e-11,-39934,42.9264], Tmin=(100,'K'), Tmax=(1053.16,'K')), NASAPolynomial(coeffs=[20.4633,0.0373109,-1.43708e-05,2.54099e-09,-1.71135e-13,-44587.6,-64.8063], Tmin=(1053.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-333.735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]CC(=O)C(C)(O)C([CH2])O(12807)',
    structure = SMILES('[CH2]CC(=O)C(C)(O)C([CH2])O'),
    E0 = (-260.473,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.05982,0.139915,-0.000193112,1.44113e-07,-4.28863e-11,-31115.6,41.9859], Tmin=(100,'K'), Tmax=(823.839,'K')), NASAPolynomial(coeffs=[17.5889,0.0445202,-1.94346e-05,3.57995e-09,-2.43306e-13,-34353.3,-49.0016], Tmin=(823.839,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-260.473,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJCO) + radical(CJCC=O)"""),
)

species(
    label = '[CH2]CC(=O)C([CH2])(O)C(C)O(6875)',
    structure = SMILES('[CH2]CC(=O)C([CH2])(O)C(C)O'),
    E0 = (-260.45,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1705.44,2804.31,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154818,'amu*angstrom^2'), symmetry=1, barrier=(3.55956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154818,'amu*angstrom^2'), symmetry=1, barrier=(3.55956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154818,'amu*angstrom^2'), symmetry=1, barrier=(3.55956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154818,'amu*angstrom^2'), symmetry=1, barrier=(3.55956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154818,'amu*angstrom^2'), symmetry=1, barrier=(3.55956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154818,'amu*angstrom^2'), symmetry=1, barrier=(3.55956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154818,'amu*angstrom^2'), symmetry=1, barrier=(3.55956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154818,'amu*angstrom^2'), symmetry=1, barrier=(3.55956,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.35196,0.14577,-0.00020325,1.49875e-07,-4.37062e-11,-31101.8,42.0175], Tmin=(100,'K'), Tmax=(842.814,'K')), NASAPolynomial(coeffs=[19.5929,0.0416236,-1.79014e-05,3.26972e-09,-2.21019e-13,-34801,-60.1018], Tmin=(842.814,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-260.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJC(C)(C=O)O) + radical(CJCC=O)"""),
)

species(
    label = '[CH2]C(O)C1(O)CO[C]1CC(9063)',
    structure = SMILES('[CH2]C(O)C1(O)CO[C]1CC'),
    E0 = (-174.715,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.66966,0.121156,-0.000139545,8.84856e-08,-2.2146e-11,-20806.2,39.5129], Tmin=(100,'K'), Tmax=(1032.75,'K')), NASAPolynomial(coeffs=[18.5728,0.0380817,-1.20983e-05,1.83512e-09,-1.09922e-13,-24738.1,-57.5915], Tmin=(1032.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-174.715,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CJCO) + radical(C2CsJOCs)"""),
)

species(
    label = '[CH2]C1(O)[C](CC)OCC1O(12808)',
    structure = SMILES('[CH2]C1(O)[C](CC)OCC1O'),
    E0 = (-252.75,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.16631,0.110821,-0.000106363,4.95494e-08,-6.33738e-12,-30209.9,34.0337], Tmin=(100,'K'), Tmax=(839.714,'K')), NASAPolynomial(coeffs=[17.6917,0.0397727,-1.29983e-05,2.06228e-09,-1.2986e-13,-34039.2,-57.594], Tmin=(839.714,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-252.75,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(C2CsJOCs) + radical(CJC(C)2O)"""),
)

species(
    label = 'CH3CH2OH(54)',
    structure = SMILES('CCO'),
    E0 = (-248.759,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
        HinderedRotor(inertia=(0.00248522,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168417,'amu*angstrom^2'), symmetry=1, barrier=(8.09857,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (46.0684,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3014.83,'J/mol'), sigma=(4.53,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.8587,-0.00374017,6.95554e-05,-8.86548e-08,3.51688e-11,-29996.1,4.80185], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[6.56244,0.0152042,-5.38968e-06,8.6225e-10,-5.12898e-14,-31525.6,-9.47302], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-248.759,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""CH3CH2OH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[CH2]C(=C=O)C([CH2])O(12809)',
    structure = SMILES('[CH2]C(O)C(=C)[C]=O'),
    E0 = (55.0586,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,287.659],'cm^-1')),
        HinderedRotor(inertia=(0.00203721,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.24041,'amu*angstrom^2'), symmetry=1, barrier=(14.1172,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.240419,'amu*angstrom^2'), symmetry=1, barrier=(14.1172,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.240426,'amu*angstrom^2'), symmetry=1, barrier=(14.1172,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18709,0.0654261,-6.67631e-05,3.55692e-08,-7.7954e-12,6720.37,27.2714], Tmin=(100,'K'), Tmax=(1081.31,'K')), NASAPolynomial(coeffs=[11.579,0.0269846,-1.34371e-05,2.6919e-09,-1.94169e-13,4473,-23.6757], Tmin=(1081.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(55.0586,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cd-CdCs(CO)) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)CJ=O) + radical(CJCO)"""),
)

species(
    label = 'C=C(O)C(C)(O)C(=O)CC(12810)',
    structure = SMILES('C=C(O)C(C)(O)C(=O)CC'),
    E0 = (-599.091,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.73839,0.126097,-0.000144626,8.73263e-08,-2.09577e-11,-71847.2,37.0682], Tmin=(100,'K'), Tmax=(1016.9,'K')), NASAPolynomial(coeffs=[20.4615,0.0387733,-1.58171e-05,2.88102e-09,-1.97247e-13,-76362.2,-70.4057], Tmin=(1016.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-599.091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH)"""),
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
    label = '[CH2]C(O)C([CH2])(O)C(C)=O(12811)',
    structure = SMILES('[CH2]C(O)C([CH2])(O)C(C)=O'),
    E0 = (-235.033,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,3580,3650,1210,1345,900,1100,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1390,370,380,2900,435,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.49513,0.122679,-0.000160898,1.08354e-07,-2.85728e-11,-28071.7,37.6098], Tmin=(100,'K'), Tmax=(933.729,'K')), NASAPolynomial(coeffs=[20.1126,0.0301147,-1.21985e-05,2.18562e-09,-1.47188e-13,-32106.9,-65.1536], Tmin=(933.729,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-235.033,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(CJC(C)(C=O)O) + radical(CJCO)"""),
)

species(
    label = '[CH2]C(O)CC(O)=C([O])CC(6910)',
    structure = SMILES('[CH2]C(O)CC(O)=C([O])CC'),
    E0 = (-351.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.22026,0.123634,-0.000130007,6.9053e-08,-1.42411e-11,-42000.4,43.9743], Tmin=(100,'K'), Tmax=(1194.95,'K')), NASAPolynomial(coeffs=[26.6275,0.0270677,-8.78671e-06,1.42304e-09,-9.18035e-14,-48894.6,-100.337], Tmin=(1194.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-351.167,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C(O)[C](O)CC(=O)CC(12812)',
    structure = SMILES('[CH2]C(O)[C](O)CC(=O)CC'),
    E0 = (-287.895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.84659,0.142063,-0.000222257,1.93863e-07,-6.63381e-11,-34428.2,42.5615], Tmin=(100,'K'), Tmax=(847.8,'K')), NASAPolynomial(coeffs=[10.9752,0.056473,-2.64231e-05,4.95386e-09,-3.36807e-13,-35700.4,-11.8603], Tmin=(847.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-287.895,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(C2CsJOH) + radical(CJCO)"""),
)

species(
    label = '[CH2]C(O)(C[CH]O)C(=O)CC(12316)',
    structure = SMILES('[CH2]C(O)(C[CH]O)C(=O)CC'),
    E0 = (-279.099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.51804,0.1522,-0.000228567,1.82211e-07,-5.70618e-11,-33341.5,41.8678], Tmin=(100,'K'), Tmax=(854.885,'K')), NASAPolynomial(coeffs=[18.0441,0.0445413,-1.95769e-05,3.56792e-09,-2.38827e-13,-36438.8,-51.6619], Tmin=(854.885,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-279.099,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CCsJOH) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = 'CCC(=O)C1(O)CCC1O(12318)',
    structure = SMILES('CCC(=O)C1(O)CCC1O'),
    E0 = (-525.474,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.22013,0.103776,-8.82909e-05,3.89614e-08,-6.89576e-12,-63002.8,36.5771], Tmin=(100,'K'), Tmax=(1353.56,'K')), NASAPolynomial(coeffs=[21.3468,0.0370878,-1.43883e-05,2.5626e-09,-1.73041e-13,-69112,-79.1272], Tmin=(1353.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-525.474,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutane)"""),
)

species(
    label = '[CH2]C(=C(CC)OO)C([CH2])O(12813)',
    structure = SMILES('[CH2]C(=C(CC)OO)C([CH2])O'),
    E0 = (-14.0301,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3615,1310,387.5,850,1000,200,800],'cm^-1')),
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
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.28365,0.114889,-0.00011357,5.89241e-08,-1.2355e-11,-1496.02,42.804], Tmin=(100,'K'), Tmax=(1144.72,'K')), NASAPolynomial(coeffs=[19.5156,0.0422097,-1.83343e-05,3.46017e-09,-2.41975e-13,-6257.88,-60.3514], Tmin=(1144.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-14.0301,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(Allyl_P) + radical(CJCO)"""),
)

species(
    label = '[CH2]C(O)=C(CC)OC([CH2])O(7292)',
    structure = SMILES('[CH2]C(O)=C(CC)OC([CH2])O'),
    E0 = (-330.253,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1200,1600],'cm^-1')),
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
    molecularWeight = (144.168,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4739.66,'J/mol'), sigma=(7.90457,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=740.32 K, Pc=21.77 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.88799,0.130798,-0.000138172,7.12636e-08,-1.40256e-11,-39454.7,44.6742], Tmin=(100,'K'), Tmax=(1332.1,'K')), NASAPolynomial(coeffs=[31.7765,0.0196796,-5.13306e-06,7.21273e-10,-4.33192e-14,-48066.3,-130.161], Tmin=(1332.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-330.253,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(CJCO) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]C(O)C([CH2])(O)C(=C)OC(12814)',
    structure = SMILES('[CH2]C(O)C([CH2])(O)C(=C)OC'),
    E0 = (-205.241,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,180,478.847,658.617,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156632,'amu*angstrom^2'), symmetry=1, barrier=(3.60128,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156632,'amu*angstrom^2'), symmetry=1, barrier=(3.60128,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156632,'amu*angstrom^2'), symmetry=1, barrier=(3.60128,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156632,'amu*angstrom^2'), symmetry=1, barrier=(3.60128,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156632,'amu*angstrom^2'), symmetry=1, barrier=(3.60128,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156632,'amu*angstrom^2'), symmetry=1, barrier=(3.60128,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156632,'amu*angstrom^2'), symmetry=1, barrier=(3.60128,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156632,'amu*angstrom^2'), symmetry=1, barrier=(3.60128,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.2059,0.128785,-0.000143188,8.07523e-08,-1.77606e-11,-24454.7,43.9508], Tmin=(100,'K'), Tmax=(1118.18,'K')), NASAPolynomial(coeffs=[25.3641,0.0301598,-1.08861e-05,1.87265e-09,-1.24865e-13,-30620.3,-92.1382], Tmin=(1118.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-205.241,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)(O)CJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C(O)C([CH2])(O)C(O)=CC(12815)',
    structure = SMILES('[CH2]C(O)C([CH2])(O)C(O)=CC'),
    E0 = (-259.951,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3615,3650,1210,1277.5,1345,900,1000,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.32239,0.132483,-0.000153646,9.02674e-08,-2.06172e-11,-31031.5,44.5396], Tmin=(100,'K'), Tmax=(1079.55,'K')), NASAPolynomial(coeffs=[25.525,0.0293006,-1.0276e-05,1.72981e-09,-1.13633e-13,-37044,-91.9392], Tmin=(1079.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-259.951,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CJCO) + radical(C=CC(C)(O)CJ)"""),
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
    label = '[CH2]C(O)[C](O)C(=O)CC(12816)',
    structure = SMILES('[CH2]C(O)C(O)=C([O])CC'),
    E0 = (-327.723,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.77594,0.108516,-0.000113488,5.81074e-08,-1.13218e-11,-39192.2,41.0486], Tmin=(100,'K'), Tmax=(1368.2,'K')), NASAPolynomial(coeffs=[27.0573,0.0163672,-3.85349e-06,4.91862e-10,-2.76311e-14,-46347.1,-104.408], Tmin=(1368.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-327.723,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(CJCO) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C(O)([CH]O)C(=O)CC(5972)',
    structure = SMILES('[CH2]C(O)([CH]O)C(=O)CC'),
    E0 = (-255.318,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1720.81,2793.18,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15521,'amu*angstrom^2'), symmetry=1, barrier=(3.56859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15521,'amu*angstrom^2'), symmetry=1, barrier=(3.56859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15521,'amu*angstrom^2'), symmetry=1, barrier=(3.56859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15521,'amu*angstrom^2'), symmetry=1, barrier=(3.56859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15521,'amu*angstrom^2'), symmetry=1, barrier=(3.56859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15521,'amu*angstrom^2'), symmetry=1, barrier=(3.56859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15521,'amu*angstrom^2'), symmetry=1, barrier=(3.56859,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.84691,0.137121,-0.000213608,1.72496e-07,-5.408e-11,-30505.3,37.2135], Tmin=(100,'K'), Tmax=(874.621,'K')), NASAPolynomial(coeffs=[17.0471,0.036257,-1.5836e-05,2.85283e-09,-1.88614e-13,-33257.4,-48.2475], Tmin=(874.621,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-255.318,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJC(C)(C=O)O) + radical(CCsJOH)"""),
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
    E0 = (-260.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-165.922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-169.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-162.786,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-199.377,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-142.555,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-260.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-136.34,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-58.7345,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-157.771,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-91.401,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-142.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-103.001,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-176.77,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-176.77,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-173.842,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-214.008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-207.732,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-210.661,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-41.4592,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-74.0252,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-222.582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (135.58,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-197.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (184.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-103.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-103.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-101.824,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-252.166,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (129.084,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-187.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (53.7483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (-116.837,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (53.8402,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (126.245,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    products = ['CH2CHOH(42)', 'C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    products = ['[CH2]C(O)C1(O)CC1([O])CC(12800)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.34238e+09,'s^-1'), n=0.889391, Ea=(94.5285,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonylbond_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 90.6 to 94.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    products = ['[CH2]C1(O)C(O)CC1([O])CC(12801)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.46159e+06,'s^-1'), n=1.55572, Ea=(90.5362,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5_SS_CO;carbonylbond_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 86.1 to 90.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[CH2]C(O)(C(=C)O)C(=O)CC(12802)'],
    products = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(170.641,'m^3/(mol*s)'), n=1.56204, Ea=(11.2897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;HJ] for rate rule [Cds-OsCs_Cds;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2][CH]O(284)', 'C=C(O)C(=O)CC(4626)'],
    products = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00392164,'m^3/(mol*s)'), n=2.41519, Ea=(15.6067,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;CJ] for rate rule [Cds-COOs_Cds;CJ]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C2H5CO(71)', '[CH2]C(O)C(=C)O(5812)'],
    products = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=0, Ea=(93.9308,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_R;CO_rad] for rate rule [Cds-OsCs_Cds;CO_rad/NonDe]
Euclidian distance = 2.82842712475
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2CHOH(42)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.30847,'m^3/(mol*s)'), n=2.06429, Ea=(63.218,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cds_Cds;CJ] + [Cds-OsH_Cds;YJ] for rate rule [Cds-OsH_Cds;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 62.0 to 63.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['OH(5)', '[CH2]C(O)C(=C)C(=O)CC(12803)'],
    products = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.52832,'m^3/(mol*s)'), n=2.02802, Ea=(14.4047,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDeCs_Cds;YJ] for rate rule [Cds-COCs_Cds;OJ_pri]
Euclidian distance = 2.2360679775
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['OH(5)', '[CH2]C(O)(C=C)C(=O)CC(6377)'],
    products = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.00674586,'m^3/(mol*s)'), n=2.3625, Ea=(84.2031,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;OJ] for rate rule [Cds-Cs\O2s/H_Cds-HH;OJ_pri]
Euclidian distance = 1.41421356237
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    products = ['[CH2]C(O)([C](C)O)C(=O)CC(6866)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.265e-07,'s^-1'), n=5.639, Ea=(102.68,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_2H;Cs_H_out_NonDe] for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_NDMustO]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(O)C(C)([O])C(=O)CC(12804)'],
    products = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    products = ['[CH2][C](O)C(C)(O)C(=O)CC(12805)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(O)(C(=O)CC)C(C)[O](6328)'],
    products = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    products = ['[CH2]C([O])(C(=O)CC)C(C)O(6870)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    products = ['[CH2]C([O])C(C)(O)C(=O)CC(6871)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    products = ['[CH2]C(O)C(C)(O)C(=O)[CH]C(12806)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(6.44e+09,'s^-1'), n=0.13, Ea=(86.6088,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 131 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    products = ['[CH2]C(O)(C(=O)[CH]C)C(C)O(6873)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(262000,'s^-1'), n=1.62, Ea=(46.4424,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_2H;Cs_H_out_H/NonDeC] for rate rule [R5H_CC(O2d)CC;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    products = ['[CH2]CC(=O)C(C)(O)C([CH2])O(12807)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(68850,'s^-1'), n=1.68, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_2H;Cs_H_out_2H] for rate rule [R5H_CCC(O2d)C;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]CC(=O)C([CH2])(O)C(C)O(6875)'],
    products = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3690,'s^-1'), n=1.79, Ea=(49.7896,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 113 used for R6H_SSSSS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R6H_SSSSS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][CH]O(284)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    products = ['[CH2]C(O)C1(O)CO[C]1CC(9063)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(7.01137e+09,'s^-1'), n=0.572544, Ea=(186.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4_S;multiplebond_intra;radadd_intra_cs2H] + [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    products = ['[CH2]C1(O)[C](CC)OCC1O(12808)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.64784e+07,'s^-1'), n=0.990488, Ea=(37.8683,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5_SS_CO;carbonyl_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CH3CH2OH(54)', '[CH2]C(=C=O)C([CH2])O(12809)'],
    products = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.79e-11,'m^3/(mol*s)'), n=3.97, Ea=(329.281,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [doublebond;Cs_OH] for rate rule [Cdd_Cd;C_pri_OH]
Euclidian distance = 1.41421356237
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    products = ['C=C(O)C(C)(O)C(=O)CC(12810)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CH2(S)(23)', '[CH2]C(O)C([CH2])(O)C(C)=O(12811)'],
    products = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;C_pri/De]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    products = ['[CH2]C(O)CC(O)=C([O])CC(6910)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    products = ['[CH2]C(O)[C](O)CC(=O)CC(12812)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ-HH;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;CO]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    products = ['[CH2]C(O)(C[CH]O)C(=O)CC(12316)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    products = ['CCC(=O)C1(O)CCC1O(12318)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(=C(CC)OO)C([CH2])O(12813)'],
    products = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O] for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O_R]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(O)=C(CC)OC([CH2])O(7292)'],
    products = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O] for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O_C]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C(O)C([CH2])(O)C(=C)OC(12814)'],
    products = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(37989.5,'s^-1'), n=2.515, Ea=(258.99,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C] + [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O] for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_C]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(O)C([CH2])(O)C(O)=CC(12815)'],
    products = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O_H] for rate rule [R_ROR;R1_doublebond_CHCH3;R2_doublebond_CsC;R_O_H]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction34',
    reactants = ['CH2(19)', '[CH2]C(O)[C](O)C(=O)CC(12816)'],
    products = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/ODMustO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction35',
    reactants = ['CH2(19)', '[CH2]C(O)([CH]O)C(=O)CC(5972)'],
    products = ['[CH2]C(O)C([CH2])(O)C(=O)CC(6869)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/CsO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '2248',
    isomers = [
        '[CH2]C(O)C([CH2])(O)C(=O)CC(6869)',
    ],
    reactants = [
        ('CH2CHOH(42)', 'C=C(O)C(=O)CC(4626)'),
        ('CH2CHOH(42)', '[CH2]C(O)=C([O])CC(4557)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '2248',
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

