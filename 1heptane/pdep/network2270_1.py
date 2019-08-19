species(
    label = '[CH2]C(O)(C(=O)CC)C([O])C=O(12333)',
    structure = SMILES('[CH2]C(O)(C(=O)CC)C([O])C=O'),
    E0 = (-308.517,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,461.37,685.843,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.161293,'amu*angstrom^2'), symmetry=1, barrier=(3.70845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161293,'amu*angstrom^2'), symmetry=1, barrier=(3.70845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161293,'amu*angstrom^2'), symmetry=1, barrier=(3.70845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161293,'amu*angstrom^2'), symmetry=1, barrier=(3.70845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161293,'amu*angstrom^2'), symmetry=1, barrier=(3.70845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161293,'amu*angstrom^2'), symmetry=1, barrier=(3.70845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161293,'amu*angstrom^2'), symmetry=1, barrier=(3.70845,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.82822,0.133213,-0.000171644,1.16877e-07,-3.16847e-11,-36900.5,44.7783], Tmin=(100,'K'), Tmax=(902.123,'K')), NASAPolynomial(coeffs=[19.0005,0.0408576,-1.80782e-05,3.39084e-09,-2.34602e-13,-40658.5,-53.5624], Tmin=(902.123,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-308.517,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = 'OCHCHO(48)',
    structure = SMILES('O=CC=O'),
    E0 = (-225.3,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,180,857.339,857.946,861.419,1717.87],'cm^-1')),
        HinderedRotor(inertia=(0.00964493,'amu*angstrom^2'), symmetry=1, barrier=(5.0669,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3660.03,'J/mol'), sigma=(4.01,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.68412,0.000478013,4.26391e-05,-5.79018e-08,2.31669e-11,-27198.5,4.51187], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[8.72507,0.00633097,-2.35575e-06,3.89783e-10,-2.37487e-14,-29102.4,-20.3904], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-225.3,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""OCHCHO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = 'CCC1([O])CC1(O)C([O])C=O(12614)',
    structure = SMILES('CCC1([O])CC1(O)C([O])C=O'),
    E0 = (-213.988,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.52095,0.107533,-8.88972e-05,2.93785e-08,-1.15833e-12,-25526,42.9523], Tmin=(100,'K'), Tmax=(1025.26,'K')), NASAPolynomial(coeffs=[24.5921,0.0300281,-1.11639e-05,2.02053e-09,-1.41364e-13,-32161.5,-89.9267], Tmin=(1025.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-213.988,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(486.397,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(C=OCOJ) + radical(CC(C)2OJ)"""),
)

species(
    label = '[CH2]C(O)(C(=O)CC)C1OC1[O](12615)',
    structure = SMILES('[CH2]C(O)(C(=O)CC)C1OC1[O]'),
    E0 = (-257.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.17487,0.139637,-0.000192604,1.41917e-07,-4.08728e-11,-30691.5,41.1073], Tmin=(100,'K'), Tmax=(928.341,'K')), NASAPolynomial(coeffs=[18.6933,0.0396212,-1.46808e-05,2.42597e-09,-1.52405e-13,-34130.8,-55.6739], Tmin=(928.341,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-257.002,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Ethylene_oxide) + radical(CCOJ) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = '[CH2]C1(O)C(C=O)OC1([O])CC(12568)',
    structure = SMILES('[CH2]C1(O)C(C=O)OC1([O])CC'),
    E0 = (-255.552,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.79488,0.118395,-0.000112881,5.37501e-08,-9.65683e-12,-30464.5,44.3302], Tmin=(100,'K'), Tmax=(1570.93,'K')), NASAPolynomial(coeffs=[28.1003,0.0219971,-3.90599e-06,3.19337e-10,-1.04219e-14,-37983.6,-111.712], Tmin=(1570.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-255.552,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(486.397,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsOs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Oxetane) + radical(CJC(C)2O) + radical(CC(C)(O)OJ)"""),
)

species(
    label = 'CCC(=O)C1(O)CC([O])C1[O](12600)',
    structure = SMILES('CCC(=O)C1(O)CC([O])C1[O]'),
    E0 = (-237.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.46894,0.11162,-0.000105236,5.07823e-08,-9.74288e-12,-28327.2,38.7471], Tmin=(100,'K'), Tmax=(1259.7,'K')), NASAPolynomial(coeffs=[22.6167,0.0351397,-1.41655e-05,2.58531e-09,-1.77672e-13,-34395.3,-83.0129], Tmin=(1259.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-237.223,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(490.554,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutane) + radical(CC(C)OJ) + radical(CC(C)OJ)"""),
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
    label = '[CH2]C(O)(C(=O)C=O)C(=O)CC(12616)',
    structure = SMILES('[CH2]C(O)(C(=O)C=O)C(=O)CC'),
    E0 = (-473.444,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3615,1277.5,1000,365,385,505,600,445,480,1700,1720,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,180,180,180,245.278,1501.73,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.148306,'amu*angstrom^2'), symmetry=1, barrier=(3.40985,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148306,'amu*angstrom^2'), symmetry=1, barrier=(3.40985,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148306,'amu*angstrom^2'), symmetry=1, barrier=(3.40985,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148306,'amu*angstrom^2'), symmetry=1, barrier=(3.40985,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148306,'amu*angstrom^2'), symmetry=1, barrier=(3.40985,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148306,'amu*angstrom^2'), symmetry=1, barrier=(3.40985,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148306,'amu*angstrom^2'), symmetry=1, barrier=(3.40985,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (157.144,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.58779,0.131267,-0.000179887,1.34712e-07,-4.07378e-11,-56748.6,40.596], Tmin=(100,'K'), Tmax=(806.665,'K')), NASAPolynomial(coeffs=[15.527,0.0463882,-2.20337e-05,4.238e-09,-2.961e-13,-59509.4,-38.2938], Tmin=(806.665,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-473.444,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = '[O]C=C[O](2537)',
    structure = SMILES('[O]C=C[O]'),
    E0 = (-18.8461,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,714.947,715.113],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.7018,-0.0027881,5.18222e-05,-5.96735e-08,2.03775e-11,-2247.83,13.3237], Tmin=(100,'K'), Tmax=(1035.3,'K')), NASAPolynomial(coeffs=[6.65087,0.0113883,-5.76543e-06,1.26602e-09,-9.87778e-14,-4228.83,-7.62441], Tmin=(1035.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-18.8461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=COJ) + radical(C=COJ)"""),
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
    label = 'C=C(O)C([O])C=O(6279)',
    structure = SMILES('C=C(O)C([O])C=O'),
    E0 = (-241.163,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,3615,1277.5,1000,1380,1390,370,380,2900,435,246.261,4000],'cm^-1')),
        HinderedRotor(inertia=(0.435262,'amu*angstrom^2'), symmetry=1, barrier=(18.3212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.419167,'amu*angstrom^2'), symmetry=1, barrier=(18.3263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.416358,'amu*angstrom^2'), symmetry=1, barrier=(18.3383,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06631,0.0549406,-3.62242e-05,-9.0911e-10,6.6365e-12,-28890.9,27.6228], Tmin=(100,'K'), Tmax=(961.285,'K')), NASAPolynomial(coeffs=[16.8081,0.0117502,-3.64622e-06,6.43589e-10,-4.69304e-14,-32948.3,-53.0627], Tmin=(961.285,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-241.163,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ)"""),
)

species(
    label = 'HCO(16)',
    structure = SMILES('[CH]=O'),
    E0 = (32.4782,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1131.19,1955.83,1955.83],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.23755,-0.00332075,1.4003e-05,-1.3424e-08,4.37416e-12,3872.41,3.30835], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.92002,0.00252279,-6.71004e-07,1.05616e-10,-7.43798e-15,3653.43,3.58077], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(32.4782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HCO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[CH2]C(O)(C=O)C(=O)CC(5970)',
    structure = SMILES('[CH2]C(O)(C=O)C(=O)CC'),
    E0 = (-365.452,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3615,1277.5,1000,375,552.5,462.5,1710,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.500114,0.10405,-0.000128411,8.70522e-08,-2.39327e-11,-43796.2,34.2575], Tmin=(100,'K'), Tmax=(883.269,'K')), NASAPolynomial(coeffs=[13.9292,0.0387082,-1.74495e-05,3.30532e-09,-2.30069e-13,-46345.3,-33.5653], Tmin=(883.269,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-365.452,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(CJC(C)(C=O)O)"""),
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
    label = 'C=C(C(=O)CC)C([O])C=O(12617)',
    structure = SMILES('C=C(C(=O)CC)C([O])C=O'),
    E0 = (-227.982,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.232014,0.0910325,-9.97511e-05,7.10352e-08,-2.26574e-11,-27291.6,37.5145], Tmin=(100,'K'), Tmax=(737.124,'K')), NASAPolynomial(coeffs=[6.8079,0.055348,-2.71341e-05,5.35817e-09,-3.82271e-13,-28261.1,7.79542], Tmin=(737.124,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-227.982,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ)"""),
)

species(
    label = '[CH2]C(O)([C](O)C=O)C(=O)CC(12618)',
    structure = SMILES('[CH2]C(O)(C(=O)CC)C(O)=C[O]'),
    E0 = (-453.198,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.81299,0.147247,-0.000189438,1.20626e-07,-2.96685e-11,-54259.5,42.2439], Tmin=(100,'K'), Tmax=(1005.61,'K')), NASAPolynomial(coeffs=[27.4061,0.0270441,-1.01369e-05,1.758e-09,-1.17013e-13,-60337.1,-103.715], Tmin=(1005.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-453.198,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=CC(O)(C=O)CJ) + radical(C=COJ)"""),
)

species(
    label = 'CCC(=O)C(C)([O])C([O])C=O(12619)',
    structure = SMILES('CCC(=O)C(C)([O])C([O])C=O'),
    E0 = (-278.146,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,375,552.5,462.5,1710,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,180,180,180,823.023,1399.79,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.172544,'amu*angstrom^2'), symmetry=1, barrier=(3.96712,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172544,'amu*angstrom^2'), symmetry=1, barrier=(3.96712,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172544,'amu*angstrom^2'), symmetry=1, barrier=(3.96712,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172544,'amu*angstrom^2'), symmetry=1, barrier=(3.96712,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172544,'amu*angstrom^2'), symmetry=1, barrier=(3.96712,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172544,'amu*angstrom^2'), symmetry=1, barrier=(3.96712,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.78382,0.115719,-0.000136883,8.27866e-08,-1.479e-11,-33290.7,43.8754], Tmin=(100,'K'), Tmax=(625.027,'K')), NASAPolynomial(coeffs=[11.6782,0.0523711,-2.42262e-05,4.61927e-09,-3.216e-13,-35168.9,-12.9536], Tmin=(625.027,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-278.146,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CC(C)(C=O)OJ)"""),
)

species(
    label = 'CCC(=O)C(C)(O)[C]([O])C=O(12620)',
    structure = SMILES('CCC(=O)C(C)(O)C([O])=C[O]'),
    E0 = (-528.617,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.12827,0.132509,-0.000159094,9.7481e-08,-2.33653e-11,-63355.2,41.2343], Tmin=(100,'K'), Tmax=(1025.89,'K')), NASAPolynomial(coeffs=[23.5845,0.0322509,-1.24997e-05,2.21545e-09,-1.49465e-13,-68630.7,-83.4714], Tmin=(1025.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-528.617,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C(O)(C(=O)CC)C(O)[C]=O(12621)',
    structure = SMILES('[CH2]C(O)(C(=O)CC)C(O)[C]=O'),
    E0 = (-392.289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.18367,0.139643,-0.000184989,1.25994e-07,-3.37221e-11,-46962,46.0888], Tmin=(100,'K'), Tmax=(918.39,'K')), NASAPolynomial(coeffs=[21.6705,0.03575,-1.53044e-05,2.82129e-09,-1.93182e-13,-51343.6,-66.9633], Tmin=(918.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-392.289,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(CJC(C)(C=O)O) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C([O])(C(=O)CC)C(O)C=O(12622)',
    structure = SMILES('[CH2]C([O])(C(=O)CC)C(O)C=O'),
    E0 = (-310.267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.74402,0.13379,-0.000182504,1.35474e-07,-4.03881e-11,-37116.5,46.1096], Tmin=(100,'K'), Tmax=(819.776,'K')), NASAPolynomial(coeffs=[16.3804,0.0453619,-2.07148e-05,3.91333e-09,-2.70518e-13,-40088.4,-37.73], Tmin=(819.776,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-310.267,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(CJC(C)(C=O)O) + radical(CC(C)(C=O)OJ)"""),
)

species(
    label = 'CCC(=O)C(C)(O)C([O])[C]=O(12623)',
    structure = SMILES('CCC(=O)C(C)(O)C([O])[C]=O'),
    E0 = (-360.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.46503,0.125232,-0.000157121,1.06557e-07,-2.9018e-11,-43125.8,44.6763], Tmin=(100,'K'), Tmax=(895.348,'K')), NASAPolynomial(coeffs=[17.0259,0.0426207,-1.87177e-05,3.50075e-09,-2.41851e-13,-46436.9,-42.4874], Tmin=(895.348,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-360.168,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCCJ=O)"""),
)

species(
    label = 'C[CH]C(=O)C(C)(O)C([O])C=O(12624)',
    structure = SMILES('CC=C([O])C(C)(O)C([O])C=O'),
    E0 = (-383.656,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.52718,0.11262,-0.000113513,5.89231e-08,-1.20536e-11,-45936.6,46.9466], Tmin=(100,'K'), Tmax=(1193.37,'K')), NASAPolynomial(coeffs=[22.6419,0.0316094,-1.1687e-05,2.03892e-09,-1.36888e-13,-51705.1,-73.9279], Tmin=(1193.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-383.656,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C(O)(C(=O)[CH]C)C(O)C=O(12625)',
    structure = SMILES('[CH2]C(O)(C([O])=CC)C(O)C=O'),
    E0 = (-413.946,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.74411,0.119782,-0.000129245,7.12068e-08,-1.53975e-11,-49573.7,47.8203], Tmin=(100,'K'), Tmax=(1132.53,'K')), NASAPolynomial(coeffs=[23.3526,0.0311422,-1.18442e-05,2.09857e-09,-1.42227e-13,-55258.3,-76.3802], Tmin=(1132.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-413.946,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(C=CC(C)(O)CJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]CC(=O)C(C)(O)C([O])C=O(12626)',
    structure = SMILES('[CH2]CC(=O)C(C)(O)C([O])C=O'),
    E0 = (-308.539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.51577,0.127091,-0.00016044,1.09522e-07,-3.00812e-11,-36915.1,44.6755], Tmin=(100,'K'), Tmax=(886.971,'K')), NASAPolynomial(coeffs=[16.9411,0.0438539,-1.96713e-05,3.71566e-09,-2.58126e-13,-40189.2,-42.1543], Tmin=(886.971,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-308.539,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(CJCC=O) + radical(C=OCOJ)"""),
)

species(
    label = '[CH2]CC(=O)C([CH2])(O)C(O)C=O(12627)',
    structure = SMILES('[CH2]CC(=O)C([CH2])(O)C(O)C=O'),
    E0 = (-340.661,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.22888,0.141437,-0.000188075,1.2865e-07,-3.46478e-11,-40751.6,46.0681], Tmin=(100,'K'), Tmax=(911.873,'K')), NASAPolynomial(coeffs=[21.5819,0.0369888,-1.62607e-05,3.03679e-09,-2.09501e-13,-45094.1,-66.6085], Tmin=(911.873,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-340.661,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(CJCC=O) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = 'CC[C]1OCC1(O)C([O])C=O(12628)',
    structure = SMILES('CC[C]1OCC1(O)C([O])C=O'),
    E0 = (-222.781,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.06113,0.107533,-0.000103873,4.96841e-08,-7.38776e-12,-26608.5,41.9738], Tmin=(100,'K'), Tmax=(869.051,'K')), NASAPolynomial(coeffs=[17.8125,0.0376147,-1.2453e-05,1.99924e-09,-1.27132e-13,-30529.1,-50.1144], Tmin=(869.051,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-222.781,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(486.397,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Oxetane) + radical(C2CsJOCs) + radical(C=OCOJ)"""),
)

species(
    label = '[CH2]C(O)(C(=O)CC)C1[CH]OO1(12629)',
    structure = SMILES('[CH2]C(O)(C(=O)CC)C1[CH]OO1'),
    E0 = (-43.7323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.08479,0.142473,-0.000205783,1.62989e-07,-5.17303e-11,-5048.85,37.9058], Tmin=(100,'K'), Tmax=(802.194,'K')), NASAPolynomial(coeffs=[15.9899,0.0483217,-2.22047e-05,4.1698e-09,-2.85761e-13,-7819.22,-44.5033], Tmin=(802.194,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-43.7323,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(12dioxetane) + radical(CJC(C)(C=O)O) + radical(CCsJOO)"""),
)

species(
    label = '[CH2]C1(O)[C](CC)OOC1C=O(12516)',
    structure = SMILES('[CH2]C1(O)[C](CC)OOC1C=O'),
    E0 = (-117.456,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.04014,0.113527,-0.000108114,5.32895e-08,-1.01847e-11,-13892.7,41.9167], Tmin=(100,'K'), Tmax=(1390.9,'K')), NASAPolynomial(coeffs=[24.7998,0.0283978,-7.74373e-06,1.07614e-09,-6.20518e-14,-20590.9,-93.6646], Tmin=(1390.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-117.456,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(486.397,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(12dioxolane) + radical(C2CsJOOC) + radical(CJC(C)2O)"""),
)

species(
    label = 'CCC(=O)C1(O)CO[CH]C1[O](12486)',
    structure = SMILES('CCC(=O)C1(O)CO[CH]C1[O]'),
    E0 = (-327.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.44815,0.113338,-0.000105866,4.17471e-08,-1.94513e-12,-39137.6,34.994], Tmin=(100,'K'), Tmax=(866.275,'K')), NASAPolynomial(coeffs=[21.6507,0.0321342,-9.33344e-06,1.37812e-09,-8.39908e-14,-44094.7,-78.6415], Tmin=(866.275,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-327.092,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(490.554,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Tetrahydrofuran) + radical(CC(C)OJ) + radical(CCsJOCs)"""),
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
    label = '[CH2]C(=C=O)C([O])C=O(12630)',
    structure = SMILES('C=C([C]=O)C([O])C=O'),
    E0 = (6.19354,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,1855,455,950,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,369.036,370.777],'cm^-1')),
        HinderedRotor(inertia=(0.00124316,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0998564,'amu*angstrom^2'), symmetry=1, barrier=(9.68043,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100022,'amu*angstrom^2'), symmetry=1, barrier=(9.67691,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5422,0.0603682,-6.35453e-05,3.70842e-08,-9.27655e-12,828.135,29.0664], Tmin=(100,'K'), Tmax=(934.737,'K')), NASAPolynomial(coeffs=[8.44027,0.0308487,-1.61733e-05,3.297e-09,-2.39769e-13,-461.409,-3.74707], Tmin=(934.737,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.19354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cd-CdCs(CO)) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)CJ=O) + radical(C=OCOJ)"""),
)

species(
    label = 'CCC(=O)C(C)(O)C(=O)C=O(12631)',
    structure = SMILES('CCC(=O)C(C)(O)C(=O)C=O'),
    E0 = (-685.056,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.00823,0.118506,-0.000143401,1.00316e-07,-2.92773e-11,-82220.3,38.9071], Tmin=(100,'K'), Tmax=(824.169,'K')), NASAPolynomial(coeffs=[12.3424,0.0537146,-2.54865e-05,4.9411e-09,-3.48355e-13,-84421,-22.9212], Tmin=(824.169,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-685.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H)"""),
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
    label = '[CH2]C(O)(C(C)=O)C([O])C=O(12632)',
    structure = SMILES('[CH2]C(O)(C(C)=O)C([O])C=O'),
    E0 = (-283.1,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,375,552.5,462.5,1710,180,180,180,180,1600,1802.98,2704.99,3200],'cm^-1')),
        HinderedRotor(inertia=(0.158024,'amu*angstrom^2'), symmetry=1, barrier=(3.63328,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158024,'amu*angstrom^2'), symmetry=1, barrier=(3.63328,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158024,'amu*angstrom^2'), symmetry=1, barrier=(3.63328,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158024,'amu*angstrom^2'), symmetry=1, barrier=(3.63328,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158024,'amu*angstrom^2'), symmetry=1, barrier=(3.63328,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158024,'amu*angstrom^2'), symmetry=1, barrier=(3.63328,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.125,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.07514,0.11145,-0.000134423,8.27005e-08,-1.99924e-11,-33866,40.7358], Tmin=(100,'K'), Tmax=(1014.3,'K')), NASAPolynomial(coeffs=[19.9243,0.0286363,-1.19536e-05,2.20523e-09,-1.52231e-13,-38126,-60.8723], Tmin=(1014.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-283.1,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(CJC(C)(C=O)O) + radical(C=OCOJ)"""),
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
    label = '[CH2]C(O)(C[O])C(=O)CC(5530)',
    structure = SMILES('[CH2]C(O)(C[O])C(=O)CC'),
    E0 = (-209.911,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,544.926,612.52,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.160633,'amu*angstrom^2'), symmetry=1, barrier=(3.69326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160633,'amu*angstrom^2'), symmetry=1, barrier=(3.69326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160633,'amu*angstrom^2'), symmetry=1, barrier=(3.69326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160633,'amu*angstrom^2'), symmetry=1, barrier=(3.69326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160633,'amu*angstrom^2'), symmetry=1, barrier=(3.69326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160633,'amu*angstrom^2'), symmetry=1, barrier=(3.69326,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4774.5,'J/mol'), sigma=(7.75786,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=745.77 K, Pc=23.2 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.23669,0.1259,-0.000195729,1.6687e-07,-5.58239e-11,-25068.2,36.334], Tmin=(100,'K'), Tmax=(850.788,'K')), NASAPolynomial(coeffs=[11.7735,0.0458102,-2.11646e-05,3.94237e-09,-2.66959e-13,-26597.1,-20.3056], Tmin=(850.788,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-209.911,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJC(C)(C=O)O) + radical(CCOJ)"""),
)

species(
    label = 'CCC(=O)[C](O)CC([O])C=O(12332)',
    structure = SMILES('CCC([O])=C(O)CC([O])C=O'),
    E0 = (-399.233,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(5510.93,'J/mol'), sigma=(8.61691,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=860.79 K, Pc=19.54 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.60869,0.110162,-9.58124e-05,3.36685e-08,-1.66773e-12,-47803.1,46.4125], Tmin=(100,'K'), Tmax=(994.333,'K')), NASAPolynomial(coeffs=[25.2744,0.0275297,-9.6449e-06,1.70054e-09,-1.18072e-13,-54410.4,-89.4722], Tmin=(994.333,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-399.233,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(C=C(C)OJ) + radical(C=OCOJ)"""),
)

species(
    label = 'CCC(=O)C[C](O)C([O])C=O(12633)',
    structure = SMILES('CCC(=O)C[C](O)C([O])C=O'),
    E0 = (-335.962,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.43414,0.130919,-0.000196137,1.68934e-07,-5.82818e-11,-40222.2,45.7153], Tmin=(100,'K'), Tmax=(817.151,'K')), NASAPolynomial(coeffs=[10.4962,0.0555084,-2.64831e-05,5.04705e-09,-3.48053e-13,-41604,-5.95735], Tmin=(817.151,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-335.962,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(C2CsJOH)"""),
)

species(
    label = 'CCC(=O)C1(O)COC1C=O(12339)',
    structure = SMILES('CCC(=O)C1(O)COC1C=O'),
    E0 = (-579.83,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.59758,0.107471,-0.000101005,5.0413e-08,-9.87102e-12,-69522.6,40.2543], Tmin=(100,'K'), Tmax=(1330.8,'K')), NASAPolynomial(coeffs=[21.9772,0.0312443,-9.03711e-06,1.31098e-09,-7.75247e-14,-75321.9,-78.432], Tmin=(1330.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-579.83,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(486.397,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + ring(Oxetane)"""),
)

species(
    label = '[CH2][C](C(=O)CC)C(C=O)OO(12634)',
    structure = SMILES('[CH2]C(=C([O])CC)C(C=O)OO'),
    E0 = (-221.813,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,3615,1310,387.5,850,1000,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,200,800,1600],'cm^-1')),
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
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.48337,0.114055,-0.000109611,5.36823e-08,-1.04753e-11,-26475.1,44.0158], Tmin=(100,'K'), Tmax=(1236.16,'K')), NASAPolynomial(coeffs=[22.5099,0.0364164,-1.54013e-05,2.87458e-09,-1.9995e-13,-32407,-76.8249], Tmin=(1236.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-221.813,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(Allyl_P) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C(=C(CC)OO)C([O])C=O(12635)',
    structure = SMILES('[CH2]C(=C(CC)OO)C([O])C=O'),
    E0 = (-62.8951,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,3615,1310,387.5,850,1000,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,200,800,1600],'cm^-1')),
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
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.811081,0.108551,-0.00010634,5.57966e-08,-1.20557e-11,-7393.53,44.1704], Tmin=(100,'K'), Tmax=(1099.15,'K')), NASAPolynomial(coeffs=[16.3678,0.0460343,-2.10237e-05,4.04993e-09,-2.86029e-13,-11170,-40.332], Tmin=(1099.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-62.8951,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(Allyl_P)"""),
)

species(
    label = 'C=C(O)[C](CC)OC([O])C=O(12331)',
    structure = SMILES('[CH2]C(O)=C(CC)OC([O])C=O'),
    E0 = (-379.118,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.3716,0.124026,-0.000129744,6.69238e-08,-1.33164e-11,-45354.3,45.8771], Tmin=(100,'K'), Tmax=(1240.93,'K')), NASAPolynomial(coeffs=[29.0136,0.0228579,-7.4518e-06,1.22364e-09,-8.01289e-14,-53143.5,-112.312], Tmin=(1240.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-379.118,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(C=C(O)CJ) + radical(C=OCOJ)"""),
)

species(
    label = '[CH2]C(O)(C(=C)OC)C([O])C=O(12636)',
    structure = SMILES('[CH2]C(O)(C(=C)OC)C([O])C=O'),
    E0 = (-253.308,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2782.5,750,1395,475,1775,1000,180,180,180,180,1099.73,1114.6,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.162677,'amu*angstrom^2'), symmetry=1, barrier=(3.74026,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162677,'amu*angstrom^2'), symmetry=1, barrier=(3.74026,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162677,'amu*angstrom^2'), symmetry=1, barrier=(3.74026,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162677,'amu*angstrom^2'), symmetry=1, barrier=(3.74026,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162677,'amu*angstrom^2'), symmetry=1, barrier=(3.74026,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162677,'amu*angstrom^2'), symmetry=1, barrier=(3.74026,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162677,'amu*angstrom^2'), symmetry=1, barrier=(3.74026,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.98736,0.119876,-0.000124545,6.48267e-08,-1.31384e-11,-30240.3,47.8028], Tmin=(100,'K'), Tmax=(1211.75,'K')), NASAPolynomial(coeffs=[26.0925,0.027183,-9.80203e-06,1.69837e-09,-1.14096e-13,-37045.4,-93.0593], Tmin=(1211.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-253.308,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]C(O)(C(O)=CC)C([O])C=O(12637)',
    structure = SMILES('[CH2]C(O)(C(O)=CC)C([O])C=O'),
    E0 = (-308.017,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,180,180,180,180,1600,1744.65,2757.55,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155746,'amu*angstrom^2'), symmetry=1, barrier=(3.5809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155746,'amu*angstrom^2'), symmetry=1, barrier=(3.5809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155746,'amu*angstrom^2'), symmetry=1, barrier=(3.5809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155746,'amu*angstrom^2'), symmetry=1, barrier=(3.5809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155746,'amu*angstrom^2'), symmetry=1, barrier=(3.5809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155746,'amu*angstrom^2'), symmetry=1, barrier=(3.5809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155746,'amu*angstrom^2'), symmetry=1, barrier=(3.5809,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.06472,0.123138,-0.000133601,7.26995e-08,-1.53687e-11,-36818.8,48.2496], Tmin=(100,'K'), Tmax=(1164.58,'K')), NASAPolynomial(coeffs=[26.008,0.0267165,-9.40851e-06,1.60497e-09,-1.0686e-13,-43357.4,-91.462], Tmin=(1164.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-308.017,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]C(O)(O[CH]C=O)C(=O)CC(12337)',
    structure = SMILES('[CH2]C(O)(OC=C[O])C(=O)CC'),
    E0 = (-452.214,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
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
    molecularWeight = (158.152,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(5120.42,'J/mol'), sigma=(8.11082,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=799.80 K, Pc=21.77 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.88148,0.151253,-0.000167649,8.46852e-08,-1.56275e-11,-54031.1,49.1652], Tmin=(100,'K'), Tmax=(1561.77,'K')), NASAPolynomial(coeffs=[41.7805,0.0017608,4.72585e-06,-1.18748e-09,8.61319e-14,-64949.7,-185.049], Tmin=(1561.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-452.214,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(CJC(C)OC) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C(O)(C(=O)CC)C([O])=CO(12638)',
    structure = SMILES('[CH2]C(O)(C(=O)CC)C([O])=CO'),
    E0 = (-456.856,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,542.196,602.52,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.159493,'amu*angstrom^2'), symmetry=1, barrier=(3.66705,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159493,'amu*angstrom^2'), symmetry=1, barrier=(3.66705,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159493,'amu*angstrom^2'), symmetry=1, barrier=(3.66705,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159493,'amu*angstrom^2'), symmetry=1, barrier=(3.66705,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159493,'amu*angstrom^2'), symmetry=1, barrier=(3.66705,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159493,'amu*angstrom^2'), symmetry=1, barrier=(3.66705,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159493,'amu*angstrom^2'), symmetry=1, barrier=(3.66705,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.93096,0.151129,-0.000203597,1.35945e-07,-3.49739e-11,-54696.3,42.7566], Tmin=(100,'K'), Tmax=(963.367,'K')), NASAPolynomial(coeffs=[27.1133,0.0263876,-9.38025e-06,1.55026e-09,-9.91381e-14,-60485.3,-101.07], Tmin=(963.367,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-456.856,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=CC(O)(C=O)CJ)"""),
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
    label = 'CCC(=O)[C](O)C([O])C=O(12639)',
    structure = SMILES('CCC([O])=C(O)C([O])C=O'),
    E0 = (-376.588,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.125,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.26422,0.101795,-0.000105214,5.39321e-08,-1.06673e-11,-45091.6,42.2684], Tmin=(100,'K'), Tmax=(1267.19,'K')), NASAPolynomial(coeffs=[24.4867,0.0192503,-6.0141e-06,9.5887e-10,-6.16271e-14,-51516.8,-87.6637], Tmin=(1267.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-376.588,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(C=C(C)OJ)"""),
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
    label = '[CH2]C(O)([CH]C=O)C(=O)CC(9040)',
    structure = SMILES('[CH2]C(O)(C=C[O])C(=O)CC'),
    E0 = (-238.64,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,506.646,639.981,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157779,'amu*angstrom^2'), symmetry=1, barrier=(3.62765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157779,'amu*angstrom^2'), symmetry=1, barrier=(3.62765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157779,'amu*angstrom^2'), symmetry=1, barrier=(3.62765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157779,'amu*angstrom^2'), symmetry=1, barrier=(3.62765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157779,'amu*angstrom^2'), symmetry=1, barrier=(3.62765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157779,'amu*angstrom^2'), symmetry=1, barrier=(3.62765,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.60403,0.125948,-0.000157352,1.02784e-07,-2.65938e-11,-28502.2,37.8558], Tmin=(100,'K'), Tmax=(946.303,'K')), NASAPolynomial(coeffs=[19.5637,0.0364722,-1.55234e-05,2.86671e-09,-1.97038e-13,-32508.4,-63.0979], Tmin=(946.303,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-238.64,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CC(O)(C=O)CJ) + radical(C=COJ)"""),
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
    E0 = (-308.517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-213.988,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-206.226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-186.474,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-237.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-229.226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-308.517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-191.42,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-308.517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-239.043,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-185.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-154.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-139.468,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-190.578,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-233.236,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-202.875,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-232.492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-221.908,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-282.785,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-255.799,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-215.005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-203.789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-122.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-42.7305,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-117.456,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-270.649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (86.7153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-245.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (136.762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-23.8237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-151.199,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-151.199,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (-300.233,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (-114.703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (80.2185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (-236.005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (5.68176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (-164.904,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (-138.414,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (-308.517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (4.97515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (4.36523,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    products = ['OCHCHO(48)', 'C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    products = ['CCC1([O])CC1(O)C([O])C=O(12614)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.34238e+09,'s^-1'), n=0.889391, Ea=(94.5285,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonylbond_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 90.6 to 94.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    products = ['[CH2]C(O)(C(=O)CC)C1OC1[O](12615)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(7.785e+11,'s^-1'), n=0.342, Ea=(102.29,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    products = ['[CH2]C1(O)C(C=O)OC1([O])CC(12568)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_O] for rate rule [R5_SS_CO;carbonylbond_intra_Nd;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    products = ['CCC(=O)C1(O)CC([O])C1[O](12600)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.46159e+06,'s^-1'), n=1.55572, Ea=(71.2943,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5_SS_CO;carbonylbond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 68.1 to 71.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[CH2]C(O)(C(=O)C=O)C(=O)CC(12616)'],
    products = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(25.1243,'m^3/(mol*s)'), n=1.86, Ea=(32.426,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO-DeNd_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C=C[O](2537)', 'C=C(O)C(=O)CC(4626)'],
    products = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00784329,'m^3/(mol*s)'), n=2.41519, Ea=(68.7968,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;CJ] for rate rule [Cds-COOs_Cds;CJ]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 68.6 to 68.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['C2H5CO(71)', 'C=C(O)C([O])C=O(6279)'],
    products = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=0, Ea=(93.9308,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_R;CO_rad] for rate rule [Cds-OsCs_Cds;CO_rad/NonDe]
Euclidian distance = 2.82842712475
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['OCHCHO(48)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(210.531,'m^3/(mol*s)'), n=1.35673, Ea=(101.726,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;CJ] + [CO-DeH_O;YJ] for rate rule [CO-DeH_O;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 99.6 to 101.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['HCO(16)', '[CH2]C(O)(C=O)C(=O)CC(5970)'],
    products = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.2e+11,'cm^3/(mol*s)'), n=0, Ea=(93.9308,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO_O;CO_pri_rad] for rate rule [CO-CsH_O;CO_pri_rad]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['OH(5)', 'C=C(C(=O)CC)C([O])C=O(12617)'],
    products = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(0.52832,'m^3/(mol*s)'), n=2.02802, Ea=(14.4047,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDeCs_Cds;YJ] for rate rule [Cds-COCs_Cds;OJ_pri]
Euclidian distance = 2.2360679775
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    products = ['[CH2]C(O)([C](O)C=O)C(=O)CC(12618)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.70223e+09,'s^-1'), n=1.15155, Ea=(153.908,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_OneDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_CO]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CCC(=O)C(C)([O])C([O])C=O(12619)'],
    products = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    products = ['CCC(=O)C(C)(O)[C]([O])C=O(12620)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    products = ['[CH2]C(O)(C(=O)CC)C(O)[C]=O(12621)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;O_rad_out;XH_out] for rate rule [R3H_SS_Cs;O_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    products = ['[CH2]C([O])(C(=O)CC)C(O)C=O(12622)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(94.0113,'s^-1'), n=2.81534, Ea=(105.641,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;O_H_out] + [R4H_SSS;O_rad_out;XH_out] for rate rule [R4H_SSS;O_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    products = ['CCC(=O)C(C)(O)C([O])[C]=O(12623)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(463.959,'s^-1'), n=2.50105, Ea=(76.0245,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_2H;XH_out] for rate rule [R4H_SSS;C_rad_out_2H;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    products = ['C[CH]C(=O)C(C)(O)C([O])C=O(12624)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(6.44e+09,'s^-1'), n=0.13, Ea=(86.6088,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 131 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    products = ['[CH2]C(O)(C(=O)[CH]C)C(O)C=O(12625)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(8e+10,'s^-1'), n=0, Ea=(25.7316,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(1000,'K'), comment="""Estimated using template [R5H_CCC;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5H_CC(O2d)CC;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    products = ['[CH2]CC(=O)C(C)(O)C([O])C=O(12626)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(68850,'s^-1'), n=1.68, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_2H;Cs_H_out_2H] for rate rule [R5H_CCC(O2d)C;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    products = ['[CH2]CC(=O)C([CH2])(O)C(O)C=O(12627)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.68e+09,'s^-1'), n=0, Ea=(93.5124,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 140 used for R6H_SSSSS;O_rad_out;Cs_H_out_2H
Exact match found for rate rule [R6H_SSSSS;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O]C=C[O](2537)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.49215e+07,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    products = ['CC[C]1OCC1(O)C([O])C=O(12628)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(7.01137e+09,'s^-1'), n=0.572544, Ea=(186.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4_S;multiplebond_intra;radadd_intra_cs2H] + [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    products = ['[CH2]C(O)(C(=O)CC)C1[CH]OO1(12629)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(265.786,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    products = ['[CH2]C1(O)[C](CC)OOC1C=O(12516)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.81184e+09,'s^-1'), n=0.551229, Ea=(191.061,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra] for rate rule [R5_SS_CO;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic
Ea raised from 187.2 to 191.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    products = ['CCC(=O)C1(O)CO[CH]C1[O](12486)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.64784e+07,'s^-1'), n=0.990488, Ea=(37.8683,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5_SS_CO;carbonyl_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CH3CH2OH(54)', '[CH2]C(=C=O)C([O])C=O(12630)'],
    products = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.79e-11,'m^3/(mol*s)'), n=3.97, Ea=(329.281,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [doublebond;Cs_OH] for rate rule [Cdd_Cd;C_pri_OH]
Euclidian distance = 1.41421356237
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    products = ['CCC(=O)C(C)(O)C(=O)C=O(12631)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['CH2(S)(23)', '[CH2]C(O)(C(C)=O)C([O])C=O(12632)'],
    products = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;C_pri/De]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction30',
    reactants = ['CO(12)', '[CH2]C(O)(C[O])C(=O)CC(5530)'],
    products = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(0.118397,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    products = ['CCC(=O)[C](O)CC([O])C=O(12332)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    products = ['CCC(=O)C[C](O)C([O])C=O(12633)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ-HH;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;CO]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    products = ['CCC(=O)C1(O)COC1C=O(12339)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2][C](C(=O)CC)C(C=O)OO(12634)'],
    products = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(3.01978e+11,'s^-1'), n=0, Ea=(107.111,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OOH_S;Y_rad_out]
Euclidian distance = 0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C(=C(CC)OO)C([O])C=O(12635)'],
    products = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O] for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O_R]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C=C(O)[C](CC)OC([O])C=O(12331)'],
    products = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O] for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O_C]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C(O)(C(=C)OC)C([O])C=O(12636)'],
    products = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(37989.5,'s^-1'), n=2.515, Ea=(258.99,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C] + [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O] for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_C]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C(O)(C(O)=CC)C([O])C=O(12637)'],
    products = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O_H] for rate rule [R_ROR;R1_doublebond_CHCH3;R2_doublebond_CsC;R_O_H]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(O)(O[CH]C=O)C(=O)CC(12337)'],
    products = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]C(O)(C(=O)CC)C([O])=CO(12638)'],
    products = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(148.339,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol
Ea raised from 146.1 to 148.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction41',
    reactants = ['CH2(19)', 'CCC(=O)[C](O)C([O])C=O(12639)'],
    products = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/ODMustO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction42',
    reactants = ['O(4)', '[CH2]C(O)([CH]C=O)C(=O)CC(9040)'],
    products = ['[CH2]C(O)(C(=O)CC)C([O])C=O(12333)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/OneDeC;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

network(
    label = '2270',
    isomers = [
        '[CH2]C(O)(C(=O)CC)C([O])C=O(12333)',
    ],
    reactants = [
        ('OCHCHO(48)', 'C=C(O)C(=O)CC(4626)'),
        ('OCHCHO(48)', '[CH2]C(O)=C([O])CC(4557)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '2270',
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

