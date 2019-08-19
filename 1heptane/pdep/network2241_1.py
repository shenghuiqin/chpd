species(
    label = 'C=C(O)C([O])(CC)C[CH]O(12313)',
    structure = SMILES('C=C(O)C([O])(CC)C[CH]O'),
    E0 = (-250.695,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,180,1050.67,1166.56,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.16155,'amu*angstrom^2'), symmetry=1, barrier=(3.71436,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16155,'amu*angstrom^2'), symmetry=1, barrier=(3.71436,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16155,'amu*angstrom^2'), symmetry=1, barrier=(3.71436,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16155,'amu*angstrom^2'), symmetry=1, barrier=(3.71436,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16155,'amu*angstrom^2'), symmetry=1, barrier=(3.71436,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16155,'amu*angstrom^2'), symmetry=1, barrier=(3.71436,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16155,'amu*angstrom^2'), symmetry=1, barrier=(3.71436,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.75922,0.123069,-0.000133041,7.49908e-08,-1.6757e-11,-29941.1,42.5599], Tmin=(100,'K'), Tmax=(1091.3,'K')), NASAPolynomial(coeffs=[21.7533,0.0368872,-1.45822e-05,2.62513e-09,-1.79012e-13,-35072.9,-72.9286], Tmin=(1091.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-250.695,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCsJOH)"""),
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
    label = '[CH2]C1(O)OC1(CC)C[CH]O(12934)',
    structure = SMILES('[CH2]C1(O)OC1(CC)C[CH]O'),
    E0 = (-204.492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.526,0.142276,-0.00016234,9.05754e-08,-1.88462e-11,-24303.7,41.8117], Tmin=(100,'K'), Tmax=(1355.77,'K')), NASAPolynomial(coeffs=[31.7807,0.0170684,-5.3317e-07,-4.35841e-10,4.66968e-14,-31943.5,-132.137], Tmin=(1355.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-204.492,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(CCsJOH) + radical(CJC(O)2C)"""),
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
    label = 'C=C(O)C([O])(CC)CC=O(12935)',
    structure = SMILES('C=C(O)C([O])(CC)CC=O'),
    E0 = (-353.964,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,180,180,180,493.254,649.952,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157352,'amu*angstrom^2'), symmetry=1, barrier=(3.61783,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157352,'amu*angstrom^2'), symmetry=1, barrier=(3.61783,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157352,'amu*angstrom^2'), symmetry=1, barrier=(3.61783,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157352,'amu*angstrom^2'), symmetry=1, barrier=(3.61783,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157352,'amu*angstrom^2'), symmetry=1, barrier=(3.61783,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157352,'amu*angstrom^2'), symmetry=1, barrier=(3.61783,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.16,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.03991,0.106371,-0.000100934,5.03546e-08,-1.01188e-11,-42386.9,39.7567], Tmin=(100,'K'), Tmax=(1196.02,'K')), NASAPolynomial(coeffs=[19.2573,0.0384894,-1.57994e-05,2.90091e-09,-1.99798e-13,-47242.1,-61.7991], Tmin=(1196.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-353.964,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C(O)C([O])(C=CO)CC(12936)',
    structure = SMILES('C=C(O)C([O])(C=CO)CC'),
    E0 = (-358.382,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,407.894,719.951,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154809,'amu*angstrom^2'), symmetry=1, barrier=(3.55936,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154809,'amu*angstrom^2'), symmetry=1, barrier=(3.55936,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154809,'amu*angstrom^2'), symmetry=1, barrier=(3.55936,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154809,'amu*angstrom^2'), symmetry=1, barrier=(3.55936,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154809,'amu*angstrom^2'), symmetry=1, barrier=(3.55936,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154809,'amu*angstrom^2'), symmetry=1, barrier=(3.55936,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.16,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.5931,0.102908,-5.91317e-05,-1.42453e-08,1.7527e-11,-42883.7,39.2752], Tmin=(100,'K'), Tmax=(957.484,'K')), NASAPolynomial(coeffs=[29.3064,0.0215534,-6.45663e-06,1.1429e-09,-8.48825e-14,-50988.8,-119.88], Tmin=(957.484,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-358.382,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ)"""),
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
    label = 'C=C(O)C(=O)C[CH]O(12937)',
    structure = SMILES('C=C(O)C(=O)C[CH]O'),
    E0 = (-337.997,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.254221,0.0992394,-0.000143581,1.10109e-07,-3.35567e-11,-40503.8,28.9371], Tmin=(100,'K'), Tmax=(805.31,'K')), NASAPolynomial(coeffs=[13.7803,0.0295251,-1.37207e-05,2.59852e-09,-1.79278e-13,-42764,-35.7317], Tmin=(805.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-337.997,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsHH) + radical(CCsJOH)"""),
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
    label = 'CCC(=O)C[CH]O(11226)',
    structure = SMILES('CCC(=O)C[CH]O'),
    E0 = (-263.14,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,3025,407.5,1350,352.5,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.124,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.352567,0.0902447,-0.00013928,1.25963e-07,-4.47479e-11,-31526.8,27.657], Tmin=(100,'K'), Tmax=(846.618,'K')), NASAPolynomial(coeffs=[5.89743,0.043601,-2.04134e-05,3.83647e-09,-2.6148e-13,-31732.9,6.15701], Tmin=(846.618,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-263.14,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CCsJOH)"""),
)

species(
    label = 'C=C(O)C([O])(CC)CC[O](12938)',
    structure = SMILES('C=C(O)C([O])(CC)CC[O]'),
    E0 = (-205.287,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2950,3100,1380,975,1025,1650,350,440,435,1725,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.913064,0.108983,-0.000104679,5.4926e-08,-1.18756e-11,-24514.2,40.8391], Tmin=(100,'K'), Tmax=(1100.9,'K')), NASAPolynomial(coeffs=[16.1442,0.047008,-2.02364e-05,3.7912e-09,-2.63679e-13,-28269.9,-43.0924], Tmin=(1100.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-205.287,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCOJ)"""),
)

species(
    label = 'C=C(O)C([O])([CH]CO)CC(12939)',
    structure = SMILES('C=C(O)C([O])([CH]CO)CC'),
    E0 = (-231.09,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,180,1065.37,1149.63,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.161848,'amu*angstrom^2'), symmetry=1, barrier=(3.72121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161848,'amu*angstrom^2'), symmetry=1, barrier=(3.72121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161848,'amu*angstrom^2'), symmetry=1, barrier=(3.72121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161848,'amu*angstrom^2'), symmetry=1, barrier=(3.72121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161848,'amu*angstrom^2'), symmetry=1, barrier=(3.72121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161848,'amu*angstrom^2'), symmetry=1, barrier=(3.72121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161848,'amu*angstrom^2'), symmetry=1, barrier=(3.72121,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.70038,0.113879,-0.000109754,5.44261e-08,-1.06459e-11,-27578.8,45.2957], Tmin=(100,'K'), Tmax=(1245.83,'K')), NASAPolynomial(coeffs=[23.6136,0.0326029,-1.18964e-05,2.06054e-09,-1.37689e-13,-33886.2,-82.3937], Tmin=(1245.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-231.09,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCJCO) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C(O)C(O)([CH]C)C[CH]O(12940)',
    structure = SMILES('C=C(O)C(O)([CH]C)C[CH]O'),
    E0 = (-279.895,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,3580,3615,3650,1210,1277.5,1345,900,1000,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,560.039,578.863,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15814,'amu*angstrom^2'), symmetry=1, barrier=(3.63595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15814,'amu*angstrom^2'), symmetry=1, barrier=(3.63595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15814,'amu*angstrom^2'), symmetry=1, barrier=(3.63595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15814,'amu*angstrom^2'), symmetry=1, barrier=(3.63595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15814,'amu*angstrom^2'), symmetry=1, barrier=(3.63595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15814,'amu*angstrom^2'), symmetry=1, barrier=(3.63595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15814,'amu*angstrom^2'), symmetry=1, barrier=(3.63595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15814,'amu*angstrom^2'), symmetry=1, barrier=(3.63595,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.30747,0.12848,-0.000143081,8.04802e-08,-1.75499e-11,-33427.7,45.8049], Tmin=(100,'K'), Tmax=(1132.03,'K')), NASAPolynomial(coeffs=[26.364,0.0271712,-8.84238e-06,1.42601e-09,-9.15414e-14,-39919.2,-96.0742], Tmin=(1132.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-279.895,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCsJOH) + radical(CCJCO)"""),
)

species(
    label = 'C=C(O)C(O)([CH][CH]O)CC(12941)',
    structure = SMILES('C=C(O)C(O)([CH][CH]O)CC'),
    E0 = (-279.895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.30747,0.12848,-0.000143081,8.04802e-08,-1.75499e-11,-33427.7,45.8049], Tmin=(100,'K'), Tmax=(1132.03,'K')), NASAPolynomial(coeffs=[26.364,0.0271712,-8.84238e-06,1.42601e-09,-9.15414e-14,-39919.2,-96.0742], Tmin=(1132.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-279.895,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCsJOH) + radical(CCJCO)"""),
)

species(
    label = 'C=C(O)C([O])([CH]C)CCO(12942)',
    structure = SMILES('C=C(O)C([O])([CH]C)CCO'),
    E0 = (-231.09,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,180,1065.37,1149.63,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.161848,'amu*angstrom^2'), symmetry=1, barrier=(3.72121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161848,'amu*angstrom^2'), symmetry=1, barrier=(3.72121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161848,'amu*angstrom^2'), symmetry=1, barrier=(3.72121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161848,'amu*angstrom^2'), symmetry=1, barrier=(3.72121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161848,'amu*angstrom^2'), symmetry=1, barrier=(3.72121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161848,'amu*angstrom^2'), symmetry=1, barrier=(3.72121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161848,'amu*angstrom^2'), symmetry=1, barrier=(3.72121,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.70038,0.113879,-0.000109754,5.44261e-08,-1.06459e-11,-27578.8,45.2957], Tmin=(100,'K'), Tmax=(1245.83,'K')), NASAPolynomial(coeffs=[23.6136,0.0326029,-1.18964e-05,2.06054e-09,-1.37689e-13,-33886.2,-82.3937], Tmin=(1245.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-231.09,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCJCO)"""),
)

species(
    label = '[CH2]CC(O)(C[CH]O)C(=C)O(12943)',
    structure = SMILES('[CH2]CC(O)(C[CH]O)C(=C)O'),
    E0 = (-274.551,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3580,3615,3650,1210,1277.5,1345,900,1000,1100,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,180,180,180,522.962,617.176,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157493,'amu*angstrom^2'), symmetry=1, barrier=(3.62108,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157493,'amu*angstrom^2'), symmetry=1, barrier=(3.62108,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157493,'amu*angstrom^2'), symmetry=1, barrier=(3.62108,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157493,'amu*angstrom^2'), symmetry=1, barrier=(3.62108,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157493,'amu*angstrom^2'), symmetry=1, barrier=(3.62108,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157493,'amu*angstrom^2'), symmetry=1, barrier=(3.62108,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157493,'amu*angstrom^2'), symmetry=1, barrier=(3.62108,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157493,'amu*angstrom^2'), symmetry=1, barrier=(3.62108,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.25834,0.131971,-0.000153482,9.08817e-08,-2.09652e-11,-32790.5,44.9207], Tmin=(100,'K'), Tmax=(1067.97,'K')), NASAPolynomial(coeffs=[24.8571,0.0304118,-1.08386e-05,1.83839e-09,-1.2115e-13,-38582.2,-87.6787], Tmin=(1067.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-274.551,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJ) + radical(CCsJOH)"""),
)

species(
    label = 'C=C([O])C(O)(CC)C[CH]O(12944)',
    structure = SMILES('C=C([O])C(O)(CC)C[CH]O'),
    E0 = (-341.993,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.78879,0.125558,-0.000144391,8.75157e-08,-2.09819e-11,-40922.2,42.9538], Tmin=(100,'K'), Tmax=(1021.26,'K')), NASAPolynomial(coeffs=[20.8533,0.0368748,-1.41338e-05,2.48454e-09,-1.66501e-13,-45546.8,-66.7574], Tmin=(1021.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-341.993,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCsJOH) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C(O)C(O)(CC)C[CH]O(12945)',
    structure = SMILES('[CH]=C(O)C(O)(CC)C[CH]O'),
    E0 = (-232.701,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,3580,3615,3650,1210,1277.5,1345,900,1000,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,180,1600,1693.6,2808.17,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154342,'amu*angstrom^2'), symmetry=1, barrier=(3.54862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154342,'amu*angstrom^2'), symmetry=1, barrier=(3.54862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154342,'amu*angstrom^2'), symmetry=1, barrier=(3.54862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154342,'amu*angstrom^2'), symmetry=1, barrier=(3.54862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154342,'amu*angstrom^2'), symmetry=1, barrier=(3.54862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154342,'amu*angstrom^2'), symmetry=1, barrier=(3.54862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154342,'amu*angstrom^2'), symmetry=1, barrier=(3.54862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154342,'amu*angstrom^2'), symmetry=1, barrier=(3.54862,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.36844,0.134749,-0.000159651,9.59734e-08,-2.24349e-11,-27753.5,43.9326], Tmin=(100,'K'), Tmax=(1055.08,'K')), NASAPolynomial(coeffs=[25.2871,0.0299012,-1.05892e-05,1.78608e-09,-1.17211e-13,-33589.3,-90.9723], Tmin=(1055.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-232.701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCsJOH) + radical(Cds_P)"""),
)

species(
    label = '[CH2]CC([O])(CCO)C(=C)O(12946)',
    structure = SMILES('[CH2]CC([O])(CCO)C(=C)O'),
    E0 = (-225.746,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2950,3100,1380,975,1025,1650,3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,350,440,435,1725,180,180,180,180,1034.73,1181.3,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.161134,'amu*angstrom^2'), symmetry=1, barrier=(3.70479,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161134,'amu*angstrom^2'), symmetry=1, barrier=(3.70479,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161134,'amu*angstrom^2'), symmetry=1, barrier=(3.70479,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161134,'amu*angstrom^2'), symmetry=1, barrier=(3.70479,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161134,'amu*angstrom^2'), symmetry=1, barrier=(3.70479,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161134,'amu*angstrom^2'), symmetry=1, barrier=(3.70479,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161134,'amu*angstrom^2'), symmetry=1, barrier=(3.70479,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.57723,0.116552,-0.000117554,6.181e-08,-1.2921e-11,-26945,44.1423], Tmin=(100,'K'), Tmax=(1162.1,'K')), NASAPolynomial(coeffs=[21.6676,0.0365401,-1.42741e-05,2.55957e-09,-1.74269e-13,-32347.4,-71.4922], Tmin=(1162.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-225.746,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJ) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C([O])C([O])(CC)CCO(12947)',
    structure = SMILES('C=C([O])C([O])(CC)CCO'),
    E0 = (-293.188,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.05914,0.109591,-0.00010666,5.62723e-08,-1.20821e-11,-35078.7,41.9996], Tmin=(100,'K'), Tmax=(1117.8,'K')), NASAPolynomial(coeffs=[17.4551,0.043339,-1.7755e-05,3.2482e-09,-2.23057e-13,-39217.7,-49.3825], Tmin=(1117.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-293.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C(O)C([O])(CC)CCO(12948)',
    structure = SMILES('[CH]=C(O)C([O])(CC)CCO'),
    E0 = (-183.896,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,350,440,435,1725,3120,650,792.5,1650,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,540.894,601.583,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.158144,'amu*angstrom^2'), symmetry=1, barrier=(3.63604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158144,'amu*angstrom^2'), symmetry=1, barrier=(3.63604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158144,'amu*angstrom^2'), symmetry=1, barrier=(3.63604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158144,'amu*angstrom^2'), symmetry=1, barrier=(3.63604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158144,'amu*angstrom^2'), symmetry=1, barrier=(3.63604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158144,'amu*angstrom^2'), symmetry=1, barrier=(3.63604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158144,'amu*angstrom^2'), symmetry=1, barrier=(3.63604,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.6751,0.119187,-0.00012323,6.62868e-08,-1.41433e-11,-21908.4,43.1103], Tmin=(100,'K'), Tmax=(1140.86,'K')), NASAPolynomial(coeffs=[21.9939,0.0362011,-1.41216e-05,2.52981e-09,-1.72181e-13,-27309.1,-74.1981], Tmin=(1140.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-183.896,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C(O)C(O)(CC)C[CH][O](12949)',
    structure = SMILES('C=C(O)C(O)(CC)C[CH][O]'),
    E0 = (-254.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.6684,0.125236,-0.000143326,8.72578e-08,-2.12014e-11,-30356.5,41.8866], Tmin=(100,'K'), Tmax=(1003.05,'K')), NASAPolynomial(coeffs=[19.5854,0.0404807,-1.6583e-05,3.02068e-09,-2.06601e-13,-34620.3,-60.7156], Tmin=(1003.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-254.092,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = 'CCC1(C[CH]O)OC[C]1O(12950)',
    structure = SMILES('CCC1(C[CH]O)OC[C]1O'),
    E0 = (-197.439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.78443,0.124248,-0.000145915,9.38647e-08,-2.37979e-11,-23535.7,38.6643], Tmin=(100,'K'), Tmax=(1019.72,'K')), NASAPolynomial(coeffs=[18.9288,0.0381954,-1.2269e-05,1.87283e-09,-1.12578e-13,-27510.3,-60.4453], Tmin=(1019.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-197.439,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CCsJOH) + radical(C2CsJOH)"""),
)

species(
    label = 'CCC1([O])CC(O)C[C]1O(12826)',
    structure = SMILES('CCC1([O])CC(O)C[C]1O'),
    E0 = (-269.085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.00306,0.0967158,-6.4858e-05,1.27322e-08,2.81094e-12,-32171.8,37.6894], Tmin=(100,'K'), Tmax=(1033.08,'K')), NASAPolynomial(coeffs=[20.44,0.037739,-1.41441e-05,2.53926e-09,-1.75471e-13,-37885.6,-72.6702], Tmin=(1033.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-269.085,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(515.497,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(CC(C)2OJ) + radical(C2CsJOH)"""),
)

species(
    label = 'C=C(O)C(O)(C=CO)CC(12951)',
    structure = SMILES('C=C(O)C(O)(C=CO)CC'),
    E0 = (-587.484,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.94807,0.108886,-6.12097e-05,-2.1268e-08,2.24989e-11,-70423.7,39.3357], Tmin=(100,'K'), Tmax=(937.144,'K')), NASAPolynomial(coeffs=[32.2267,0.0180401,-3.86843e-06,5.90454e-10,-4.52493e-14,-79245.2,-136.21], Tmin=(937.144,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-587.484,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C(O)C(O)(CC)CC=O(12952)',
    structure = SMILES('C=C(O)C(O)(CC)CC=O'),
    E0 = (-583.067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.68376,0.115773,-0.000115087,5.91923e-08,-1.20272e-11,-69914.3,40.8526], Tmin=(100,'K'), Tmax=(1199.65,'K')), NASAPolynomial(coeffs=[22.9233,0.0337256,-1.24978e-05,2.1813e-09,-1.46388e-13,-75818.2,-82.3412], Tmin=(1199.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-583.067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH)"""),
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
    label = 'C=C(O)C(C)([O])C[CH]O(12953)',
    structure = SMILES('C=C(O)C(C)([O])C[CH]O'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.14686,0.108696,-0.000120647,6.88117e-08,-1.5414e-11,-27102.2,38.1158], Tmin=(100,'K'), Tmax=(1094.21,'K')), NASAPolynomial(coeffs=[20.7662,0.0285919,-1.08375e-05,1.90964e-09,-1.28799e-13,-31897.8,-69.5751], Tmin=(1094.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-226.914,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C(O)C([O])(CC)C(=C)O(7272)',
    structure = SMILES('[CH2]C(O)C([O])(CC)C(=C)O'),
    E0 = (-232.046,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,476.214,661.706,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15667,'amu*angstrom^2'), symmetry=1, barrier=(3.60216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15667,'amu*angstrom^2'), symmetry=1, barrier=(3.60216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15667,'amu*angstrom^2'), symmetry=1, barrier=(3.60216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15667,'amu*angstrom^2'), symmetry=1, barrier=(3.60216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15667,'amu*angstrom^2'), symmetry=1, barrier=(3.60216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15667,'amu*angstrom^2'), symmetry=1, barrier=(3.60216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15667,'amu*angstrom^2'), symmetry=1, barrier=(3.60216,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(5105.78,'J/mol'), sigma=(8.38548,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=797.51 K, Pc=19.65 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.97881,0.121303,-0.00012457,6.52765e-08,-1.3414e-11,-27684.8,44.0852], Tmin=(100,'K'), Tmax=(1191.33,'K')), NASAPolynomial(coeffs=[24.7714,0.0314856,-1.14799e-05,1.99084e-09,-1.33386e-13,-34058.4,-89.6517], Tmin=(1191.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-232.046,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CJCO)"""),
)

species(
    label = 'C=C(O)C1(CC)CC(O)O1(12319)',
    structure = SMILES('C=C(O)C1(CC)CC(O)O1'),
    E0 = (-516.467,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.983496,0.0884102,-1.45898e-05,-6.4649e-08,3.8616e-11,-61917.2,36.1867], Tmin=(100,'K'), Tmax=(895.094,'K')), NASAPolynomial(coeffs=[26.9196,0.0225866,-2.9369e-06,1.49883e-10,-4.79607e-15,-69270.7,-108.511], Tmin=(895.094,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-516.467,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Oxetane)"""),
)

species(
    label = 'C=[C]C(CC)(C[CH]O)OO(12954)',
    structure = SMILES('C=[C]C(CC)(C[CH]O)OO'),
    E0 = (44.108,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,250.137,892.606,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152557,'amu*angstrom^2'), symmetry=1, barrier=(3.5076,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152557,'amu*angstrom^2'), symmetry=1, barrier=(3.5076,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152557,'amu*angstrom^2'), symmetry=1, barrier=(3.5076,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152557,'amu*angstrom^2'), symmetry=1, barrier=(3.5076,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152557,'amu*angstrom^2'), symmetry=1, barrier=(3.5076,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152557,'amu*angstrom^2'), symmetry=1, barrier=(3.5076,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152557,'amu*angstrom^2'), symmetry=1, barrier=(3.5076,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152557,'amu*angstrom^2'), symmetry=1, barrier=(3.5076,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.61877,0.129243,-0.00015849,1.06733e-07,-2.91002e-11,5502.36,42.0353], Tmin=(100,'K'), Tmax=(891.393,'K')), NASAPolynomial(coeffs=[16.577,0.047588,-2.10761e-05,3.95728e-09,-2.74064e-13,2258.6,-43.6561], Tmin=(891.393,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(44.108,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CCsJOH)"""),
)

species(
    label = 'CCC([O])(C[CH]O)C(C)=O(12955)',
    structure = SMILES('CCC([O])(C[CH]O)C(C)=O'),
    E0 = (-247.091,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.51782,0.130238,-0.000183869,1.49361e-07,-4.89175e-11,-29527.9,41.9903], Tmin=(100,'K'), Tmax=(825.575,'K')), NASAPolynomial(coeffs=[12.2956,0.0526855,-2.36583e-05,4.39869e-09,-2.99576e-13,-31446.7,-19.811], Tmin=(825.575,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-247.091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(CC(C)(C=O)OJ) + radical(CCsJOH)"""),
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
    label = 'C=C(O)[C](CC)C[CH]O(12956)',
    structure = SMILES('[CH2]C(O)=C(CC)C[CH]O'),
    E0 = (-149.419,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
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
    molecularWeight = (128.169,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.50292,0.115284,-0.000124383,6.99019e-08,-1.5478e-11,-17767.9,38.3037], Tmin=(100,'K'), Tmax=(1105.86,'K')), NASAPolynomial(coeffs=[21.4425,0.0322881,-1.18059e-05,2.0347e-09,-1.35326e-13,-22842.7,-74.7032], Tmin=(1105.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-149.419,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(CCsJOH) + radical(C=C(O)CJ)"""),
)

species(
    label = 'HCOH(T)(285)',
    structure = SMILES('[CH]O'),
    E0 = (205.906,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,403.876,3308.82],'cm^-1')),
        HinderedRotor(inertia=(0.0103144,'amu*angstrom^2'), symmetry=1, barrier=(22.7121,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.75938,0.0029613,8.90411e-06,-1.35016e-08,5.39816e-12,24775.6,6.76286], Tmin=(100,'K'), Tmax=(940.429,'K')), NASAPolynomial(coeffs=[5.09112,0.00321239,-9.31686e-07,1.59615e-10,-1.15729e-14,24263.5,-0.971], Tmin=(940.429,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(78.9875,'J/(mol*K)'), label="""HCOH(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C([O])(CC)C(=C)O(6187)',
    structure = SMILES('[CH2]C([O])(CC)C(=C)O'),
    E0 = (-33.9413,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,180,1600,1602.72,2893.45,3200],'cm^-1')),
        HinderedRotor(inertia=(0.149077,'amu*angstrom^2'), symmetry=1, barrier=(3.42758,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149077,'amu*angstrom^2'), symmetry=1, barrier=(3.42758,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149077,'amu*angstrom^2'), symmetry=1, barrier=(3.42758,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149077,'amu*angstrom^2'), symmetry=1, barrier=(3.42758,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149077,'amu*angstrom^2'), symmetry=1, barrier=(3.42758,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.676536,0.0908249,-8.54045e-05,4.10682e-08,-7.75587e-12,-3903.68,35.0552], Tmin=(100,'K'), Tmax=(1293.32,'K')), NASAPolynomial(coeffs=[20.7584,0.0245293,-8.51315e-06,1.43237e-09,-9.40947e-14,-9448.05,-73.8692], Tmin=(1293.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-33.9413,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)(O)CJ) + radical(C=CC(C)2OJ)"""),
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
    E0 = (-250.695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-203.012,'kJ/mol'),
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
    E0 = (-106.286,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-139.853,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-202.09,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-182.767,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-111.713,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-250.695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-116.294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-122.306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-174.333,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-175.414,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-184.836,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-190.871,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-148.887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-188.393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-171.284,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-221.658,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-135.025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-157.118,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-41.4592,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-125.802,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-164.923,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-187.295,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-225.721,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (192.948,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-73.4196,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-242.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (151.219,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-46.5157,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (93.5854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (171.965,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C(O)C([O])(CC)C[CH]O(12313)'],
    products = ['CH2CHOH(42)', 'C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C(O)C([O])(CC)C[CH]O(12313)'],
    products = ['[CH2]C1(O)OC1(CC)C[CH]O(12934)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.52e+08,'s^-1'), n=0.89, Ea=(47.6831,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H_secNd;radadd_intra] for rate rule [R4_S_D;doublebond_intra_2H_secNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=C(O)C([O])(CC)C[CH]O(12313)'],
    products = ['[CH2]C1(O)C(O)CC1([O])CC(12801)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.53628e+06,'s^-1'), n=1.28, Ea=(80.7805,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5_SS_D;doublebond_intra_2H;radadd_intra_csHNd] + [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_cs] for rate rule [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_csHNd]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 79.4 to 80.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C=C(O)C([O])(CC)CC=O(12935)'],
    products = ['C=C(O)C([O])(CC)C[CH]O(12313)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4e+09,'cm^3/(mol*s)'), n=1.39, Ea=(35.8862,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2818 used for Od_CO-CsH;HJ
Exact match found for rate rule [Od_CO-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C=C(O)C([O])(C=CO)CC(12936)'],
    products = ['C=C(O)C([O])(CC)C[CH]O(12313)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.72e+08,'cm^3/(mol*s)'), n=1.477, Ea=(6.73624,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2825 used for Cds-CsH_Cds-OsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-OsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C2H5(29)', 'C=C(O)C(=O)C[CH]O(12937)'],
    products = ['C=C(O)C([O])(CC)C[CH]O(12313)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7.94e+10,'cm^3/(mol*s)'), n=0, Ea=(28.0328,'kJ/mol'), T0=(1,'K'), Tmin=(333,'K'), Tmax=(363,'K'), comment="""Estimated using template [CO_O;CsJ-CsHH] for rate rule [CO-CdCs_O;CsJ-CsHH]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2][CH]O(284)', 'C=C(O)C(=O)CC(4626)'],
    products = ['C=C(O)C([O])(CC)C[CH]O(12313)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.97e+07,'cm^3/(mol*s)'), n=1.88, Ea=(32.2168,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO-CdCs_O;YJ] for rate rule [CO-CdCs_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH2COH(99)', 'CCC(=O)C[CH]O(11226)'],
    products = ['C=C(O)C([O])(CC)C[CH]O(12313)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.16e+10,'cm^3/(mol*s)'), n=0, Ea=(48.1578,'kJ/mol'), T0=(1,'K'), Tmin=(413,'K'), Tmax=(563,'K'), comment="""Estimated using template [CO-CsCs_O;CJ] for rate rule [CO-CsCs_O;CdsJ-O2s]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CH2CHOH(42)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['C=C(O)C([O])(CC)C[CH]O(12313)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.00668046,'m^3/(mol*s)'), n=2.5095, Ea=(72.9737,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-OsH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 68.7 to 73.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=C(O)C([O])(CC)CC[O](12938)'],
    products = ['C=C(O)C([O])(CC)C[CH]O(12313)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(153000,'s^-1'), n=2.26, Ea=(88.9937,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 325 used for R2H_S;O_rad_out;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R2H_S;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=C(O)C([O])([CH]CO)CC(12939)'],
    products = ['C=C(O)C([O])(CC)C[CH]O(12313)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.4e-20,'s^-1'), n=9.13, Ea=(108.784,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 341 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=C(O)C(O)([CH]C)C[CH]O(12940)'],
    products = ['C=C(O)C([O])(CC)C[CH]O(12313)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(5.71,'s^-1'), n=3.021, Ea=(105.562,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 319 used for R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C(O)C([O])(CC)C[CH]O(12313)'],
    products = ['C=C(O)C(O)([CH][CH]O)CC(12941)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(223829,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=C(O)C([O])([CH]C)CCO(12942)'],
    products = ['C=C(O)C([O])(CC)C[CH]O(12313)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.00310895,'s^-1'), n=4.065, Ea=(46.2541,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_1H] for rate rule [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]CC(O)(C[CH]O)C(=C)O(12943)'],
    products = ['C=C(O)C([O])(CC)C[CH]O(12313)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=C(O)C([O])(CC)C[CH]O(12313)'],
    products = ['C=C([O])C(O)(CC)C[CH]O(12944)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1070.11,'s^-1'), n=2.50856, Ea=(101.808,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;O_H_out] + [R4H_SSS;O_rad_out;XH_out] + [R4H_SS(Cd)S;Y_rad_out;XH_out] for rate rule [R4H_SS(Cd)S;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C(O)C(O)(CC)C[CH]O(12945)'],
    products = ['C=C(O)C([O])(CC)C[CH]O(12313)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]CC([O])(CCO)C(=C)O(12946)'],
    products = ['C=C(O)C([O])(CC)C[CH]O(12313)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(30977.5,'s^-1'), n=1.81667, Ea=(54.4617,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_2H;Cs_H_out_1H] for rate rule [R5H_CCC;C_rad_out_2H;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=C(O)C([O])(CC)C[CH]O(12313)'],
    products = ['C=C([O])C([O])(CC)CCO(12947)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.0508,'s^-1'), n=3.24, Ea=(29.037,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_1H;XH_out] for rate rule [R5H_CCC(Cd);C_rad_out_H/NonDeO;O_H_out]
Euclidian distance = 1.73205080757
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C(O)C([O])(CC)CCO(12948)'],
    products = ['C=C(O)C([O])(CC)C[CH]O(12313)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(107770,'s^-1'), n=1.98401, Ea=(48.8716,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5H;Cd_rad_out_singleH;Cs_H_out] + [R5H_RSSR;Y_rad_out;Cs_H_out_H/NonDeO] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cs_H_out_H/NonDeO]
Euclidian distance = 3.60555127546
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=C(O)C([O])(CC)C[CH]O(12313)'],
    products = ['C=C(O)C(O)(CC)C[CH][O](12949)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(36450.5,'s^-1'), n=1.97932, Ea=(93.5766,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;XH_out] for rate rule [R5HJ_3;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][CH]O(284)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['C=C(O)C([O])(CC)C[CH]O(12313)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C(O)C([O])(CC)C[CH]O(12313)'],
    products = ['CCC1(C[CH]O)OC[C]1O(12950)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.21748e+08,'s^-1'), n=0.95, Ea=(124.892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R4_S_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=C(O)C([O])(CC)C[CH]O(12313)'],
    products = ['CCC1([O])CC(O)C[C]1O(12826)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.95277e+09,'s^-1'), n=0.6, Ea=(85.772,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_cs] for rate rule [R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_csHO]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=C(O)C([O])(CC)C[CH]O(12313)'],
    products = ['C=C(O)C(O)(C=CO)CC(12951)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=C(O)C([O])(CC)C[CH]O(12313)'],
    products = ['C=C(O)C(O)(CC)CC=O(12952)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CH2(S)(23)', 'C=C(O)C(C)([O])C[CH]O(12953)'],
    products = ['C=C(O)C([O])(CC)C[CH]O(12313)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(O)C([O])(CC)C(=C)O(7272)'],
    products = ['C=C(O)C([O])(CC)C[CH]O(12313)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=C(O)C([O])(CC)C[CH]O(12313)'],
    products = ['C=C(O)C1(CC)CC(O)O1(12319)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;O_rad;Cpri_rad_out_H/NonDeO]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=[C]C(CC)(C[CH]O)OO(12954)'],
    products = ['C=C(O)C([O])(CC)C[CH]O(12313)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(3.01978e+11,'s^-1'), n=0, Ea=(107.111,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2OOH_S;Y_rad_out] for rate rule [R2OOH_S;Cd_rad_out_double]
Euclidian distance = 2.0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=C(O)C([O])(CC)C[CH]O(12313)'],
    products = ['CCC([O])(C[CH]O)C(C)=O(12955)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(204.179,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction32',
    reactants = ['O(4)', 'C=C(O)[C](CC)C[CH]O(12956)'],
    products = ['C=C(O)C([O])(CC)C[CH]O(12313)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/Cs2;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction33',
    reactants = ['HCOH(T)(285)', '[CH2]C([O])(CC)C(=C)O(6187)'],
    products = ['C=C(O)C([O])(CC)C[CH]O(12313)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '2241',
    isomers = [
        'C=C(O)C([O])(CC)C[CH]O(12313)',
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
    label = '2241',
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

