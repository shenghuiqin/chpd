species(
    label = '[CH2]OO[C](C)OO(5260)',
    structure = SMILES('[CH2]OO[C](C)OO'),
    E0 = (84.9753,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,360,370,350,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
        HinderedRotor(inertia=(0.765514,'amu*angstrom^2'), symmetry=1, barrier=(38.803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.765514,'amu*angstrom^2'), symmetry=1, barrier=(38.803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0323312,'amu*angstrom^2'), symmetry=1, barrier=(38.803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167222,'amu*angstrom^2'), symmetry=1, barrier=(8.47629,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.765512,'amu*angstrom^2'), symmetry=1, barrier=(38.803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167221,'amu*angstrom^2'), symmetry=1, barrier=(8.47629,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.944732,0.0735817,-8.49389e-05,5.46166e-08,-1.47529e-11,10324.7,26.1307], Tmin=(100,'K'), Tmax=(881.621,'K')), NASAPolynomial(coeffs=[9.70466,0.0338355,-1.73114e-05,3.47582e-09,-2.50353e-13,8780.21,-15.0269], Tmin=(881.621,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.9753,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-OsHHH) + radical(Cs_P) + radical(CsJOOC)"""),
)

species(
    label = 'CH3C(O)OOH(67)',
    structure = SMILES('CC(=O)OO'),
    E0 = (-367.861,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,367.717,367.721,367.767,1420.36],'cm^-1')),
        HinderedRotor(inertia=(0.00493127,'amu*angstrom^2'), symmetry=1, barrier=(7.05973,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00124622,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.712235,'amu*angstrom^2'), symmetry=1, barrier=(68.3461,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (76.0514,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3635.24,'J/mol'), sigma=(5.76225,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=567.82 K, Pc=43.11 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43748,0.0339384,-2.25467e-05,5.77205e-09,1.28452e-13,-44070.1,19.8194], Tmin=(298,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.64132,0.0217231,-1.11238e-05,2.75682e-09,-2.67754e-13,-45140.9,-1.61172], Tmin=(1000,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-367.861,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), label="""CH3C(O)OOH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'CH2O(15)',
    structure = SMILES('C=O'),
    E0 = (-119.055,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.79372,-0.00990833,3.7322e-05,-3.79285e-08,1.31773e-11,-14379.2,0.602798], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.16953,0.00619321,-2.25056e-06,3.65976e-10,-2.20149e-14,-14548.7,6.04208], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-119.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH2O""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = '[CH2]OOC(=C)OO(5657)',
    structure = SMILES('[CH2]OOC(=C)OO'),
    E0 = (118.731,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,350,440,435,1725],'cm^-1')),
        HinderedRotor(inertia=(2.09533,'amu*angstrom^2'), symmetry=1, barrier=(48.1758,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.192849,'amu*angstrom^2'), symmetry=1, barrier=(4.43397,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.192152,'amu*angstrom^2'), symmetry=1, barrier=(4.41796,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.192398,'amu*angstrom^2'), symmetry=1, barrier=(4.42361,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.09589,'amu*angstrom^2'), symmetry=1, barrier=(48.1886,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.069,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.866263,0.0782395,-0.000129425,1.20918e-07,-4.43688e-11,14383.8,26.2879], Tmin=(100,'K'), Tmax=(814.539,'K')), NASAPolynomial(coeffs=[6.00611,0.0354671,-1.83725e-05,3.60189e-09,-2.51572e-13,14128.1,6.11555], Tmin=(814.539,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(118.731,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-OsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CsJOOC)"""),
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
    label = '[CH2]OOC(C)=O(5658)',
    structure = SMILES('[CH2]OOC(C)=O'),
    E0 = (-252.233,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,539.971,542.137,543.072,543.417],'cm^-1')),
        HinderedRotor(inertia=(0.000569639,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.231558,'amu*angstrom^2'), symmetry=1, barrier=(48.2345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.232615,'amu*angstrom^2'), symmetry=1, barrier=(48.287,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.232013,'amu*angstrom^2'), symmetry=1, barrier=(48.2104,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (89.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90384,0.0326432,1.8558e-05,-4.8508e-08,2.07201e-11,-30249.2,21.0761], Tmin=(100,'K'), Tmax=(992.12,'K')), NASAPolynomial(coeffs=[13.7244,0.0172339,-6.90137e-06,1.36257e-09,-1.02352e-13,-34181.8,-43.8564], Tmin=(992.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-252.233,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(Cs-(Cds-O2d)HHH) + group(Cs-OsHHH) + group(Cds-OdCsOs) + radical(CsJOOC)"""),
)

species(
    label = '[CH2]OOC([CH2])OO(5659)',
    structure = SMILES('[CH2]OOC([CH2])OO'),
    E0 = (93.6915,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1310,387.5,850,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1390,370,380,2900,435],'cm^-1')),
        HinderedRotor(inertia=(0.675229,'amu*angstrom^2'), symmetry=1, barrier=(33.222,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0388919,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.675234,'amu*angstrom^2'), symmetry=1, barrier=(33.222,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172151,'amu*angstrom^2'), symmetry=1, barrier=(8.47006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.675237,'amu*angstrom^2'), symmetry=1, barrier=(33.222,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.675233,'amu*angstrom^2'), symmetry=1, barrier=(33.222,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.823363,0.0758167,-8.86377e-05,5.63692e-08,-1.48538e-11,11377.8,27.8717], Tmin=(100,'K'), Tmax=(908.099,'K')), NASAPolynomial(coeffs=[10.7651,0.0320259,-1.6305e-05,3.26799e-09,-2.35207e-13,9572.2,-19.1331], Tmin=(908.099,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(93.6915,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-OsHHH) + radical(CsJOOC) + radical(CJCOOH)"""),
)

species(
    label = '[CH2]OOC(C)O[O](5660)',
    structure = SMILES('[CH2]OOC(C)O[O]'),
    E0 = (31.7336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05078,0.0707429,-8.47061e-05,5.93289e-08,-1.75959e-11,3917.65,25.8167], Tmin=(100,'K'), Tmax=(806.116,'K')), NASAPolynomial(coeffs=[8.37889,0.0343794,-1.704e-05,3.36676e-09,-2.39968e-13,2736.22,-7.95765], Tmin=(806.116,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.7336,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-OsHHH) + radical(CsJOOC) + radical(ROOJ)"""),
)

species(
    label = '[CH2][C](OO)OOC(5661)',
    structure = SMILES('[CH2][C](OO)OOC'),
    E0 = (105.341,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,360,370,350,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
        HinderedRotor(inertia=(0.220809,'amu*angstrom^2'), symmetry=1, barrier=(7.44,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0073913,'amu*angstrom^2'), symmetry=1, barrier=(34.5917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220853,'amu*angstrom^2'), symmetry=1, barrier=(7.43993,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02672,'amu*angstrom^2'), symmetry=1, barrier=(34.5917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02678,'amu*angstrom^2'), symmetry=1, barrier=(34.5916,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02689,'amu*angstrom^2'), symmetry=1, barrier=(34.5916,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.755579,0.0781494,-0.000104237,7.95033e-08,-2.52559e-11,12780.2,29.1018], Tmin=(100,'K'), Tmax=(759.589,'K')), NASAPolynomial(coeffs=[9.08049,0.0343121,-1.76723e-05,3.53155e-09,-2.52608e-13,11515.5,-8.77225], Tmin=(759.589,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(105.341,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-OsHHH) + radical(CJCOOH) + radical(Cs_P)"""),
)

species(
    label = 'COO[C](C)O[O](5662)',
    structure = SMILES('COO[C](C)O[O]'),
    E0 = (43.3832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60368,0.0630656,-4.65138e-05,-3.38566e-08,5.8986e-11,5293.48,24.9598], Tmin=(100,'K'), Tmax=(481.828,'K')), NASAPolynomial(coeffs=[6.99781,0.0361064,-1.80662e-05,3.5464e-09,-2.50196e-13,4566.8,0.728113], Tmin=(481.828,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(43.3832,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-OsHHH) + radical(Cs_P) + radical(ROOJ)"""),
)

species(
    label = '[CH2][O](167)',
    structure = SMILES('[CH2][O]'),
    E0 = (192.903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88409,-0.00363885,3.28543e-05,-4.13611e-08,1.59631e-11,23210.8,7.47983], Tmin=(100,'K'), Tmax=(933.06,'K')), NASAPolynomial(coeffs=[6.69335,0.000289989,8.61416e-07,-1.56351e-10,7.33778e-15,21991.3,-9.6043], Tmin=(933.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(H3COJ) + radical(CsJOH)"""),
)

species(
    label = 'C[C]([O])OO(5646)',
    structure = SMILES('C[C]([O])OO'),
    E0 = (45.6545,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.204777,'amu*angstrom^2'), symmetry=1, barrier=(4.70822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203884,'amu*angstrom^2'), symmetry=1, barrier=(4.68769,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0414863,'amu*angstrom^2'), symmetry=1, barrier=(47.3212,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (76.0514,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.16949,0.0456319,-6.91049e-05,6.57388e-08,-2.5341e-11,5551.65,19.8558], Tmin=(100,'K'), Tmax=(783.975,'K')), NASAPolynomial(coeffs=[4.05265,0.0263061,-1.35357e-05,2.67394e-09,-1.88616e-13,5555,13.1335], Tmin=(783.975,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(45.6545,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Cs_P)"""),
)

species(
    label = 'C=C(OO)OOC(5663)',
    structure = SMILES('C=C(OO)OOC'),
    E0 = (-74.8662,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.917167,0.0779416,-0.00012762,1.22174e-07,-4.58683e-11,-8903.16,26.4081], Tmin=(100,'K'), Tmax=(816.702,'K')), NASAPolynomial(coeffs=[4.43348,0.0401917,-2.05836e-05,4.02439e-09,-2.8084e-13,-8792.91,14.3472], Tmin=(816.702,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-74.8662,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-OsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH)"""),
)

species(
    label = 'CC1(COO1)OO(5261)',
    structure = SMILES('CC1(COO1)OO'),
    E0 = (-192.796,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.532005,0.0607655,-2.24784e-05,-2.74915e-08,1.83322e-11,-23049.1,20.0167], Tmin=(100,'K'), Tmax=(950.946,'K')), NASAPolynomial(coeffs=[21.561,0.0102213,-2.551e-06,4.61375e-10,-3.79309e-14,-28762.7,-89.3913], Tmin=(950.946,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-192.796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(12dioxetane)"""),
)

species(
    label = '[CH2]OOC(C)([O])O(5664)',
    structure = SMILES('[CH2]OOC(C)([O])O'),
    E0 = (-146.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.077,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.962369,0.0755981,-0.000117413,1.08585e-07,-3.9759e-11,-17512.9,27.3101], Tmin=(100,'K'), Tmax=(824.727,'K')), NASAPolynomial(coeffs=[4.99219,0.038778,-1.90252e-05,3.65509e-09,-2.52729e-13,-17590.1,12.2069], Tmin=(824.727,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-146.449,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsOs) + group(Cs-CsHHH) + group(Cs-OsHHH) + radical(CCOJ) + radical(CsJOOC)"""),
)

species(
    label = 'C[C]OO(178)',
    structure = SMILES('C[C]OO'),
    E0 = (283.255,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,427.45],'cm^-1')),
        HinderedRotor(inertia=(0.110494,'amu*angstrom^2'), symmetry=1, barrier=(14.3257,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110496,'amu*angstrom^2'), symmetry=1, barrier=(14.3257,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1105,'amu*angstrom^2'), symmetry=1, barrier=(14.3257,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.4285,0.0278106,-1.11983e-05,-7.95465e-09,5.41923e-12,34130.2,15.5328], Tmin=(100,'K'), Tmax=(1026.15,'K')), NASAPolynomial(coeffs=[11.3253,0.00771013,-3.12825e-06,6.48476e-10,-5.00284e-14,31536.7,-31.359], Tmin=(1026.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(283.255,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CH2_triplet)"""),
)

species(
    label = '[CH2]O[O](92)',
    structure = SMILES('[CH2]O[O]'),
    E0 = (206.805,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2219.2,2221.31],'cm^-1')),
        HinderedRotor(inertia=(0.156598,'amu*angstrom^2'), symmetry=1, barrier=(14.2064,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.45042,0.0115293,-8.64543e-06,3.98163e-09,-7.73155e-13,24893.2,10.7919], Tmin=(100,'K'), Tmax=(1212.82,'K')), NASAPolynomial(coeffs=[5.07483,0.00617178,-2.01911e-06,3.39169e-10,-2.23136e-14,24499.1,2.64172], Tmin=(1212.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(206.805,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CsJOOH) + radical(ROOJ)"""),
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
    label = 'C[C](O[O])OO(5665)',
    structure = SMILES('C[C](O[O])OO'),
    E0 = (43.459,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3615,1310,387.5,850,1000,492.5,1135,1000],'cm^-1')),
        HinderedRotor(inertia=(0.283797,'amu*angstrom^2'), symmetry=1, barrier=(6.52506,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.720352,'amu*angstrom^2'), symmetry=1, barrier=(16.5623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284519,'amu*angstrom^2'), symmetry=1, barrier=(6.54165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.19905,'amu*angstrom^2'), symmetry=1, barrier=(50.5605,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0508,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66765,0.057692,-7.36115e-05,4.31886e-08,-4.98778e-12,5304.81,22.363], Tmin=(100,'K'), Tmax=(585.451,'K')), NASAPolynomial(coeffs=[8.17275,0.0233898,-1.1712e-05,2.2943e-09,-1.6156e-13,4369.3,-7.02217], Tmin=(585.451,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(43.459,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(ROOJ) + radical(Cs_P)"""),
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
    E0 = (84.9753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (333.507,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (84.9753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (225.976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (258.681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (159.524,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (158.624,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (238.558,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (119.493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (93.2596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (171.166,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (490.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (425.022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]OO[C](C)OO(5260)'],
    products = ['CH3C(O)OOH(67)', 'CH2O(15)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(3)', '[CH2]OOC(=C)OO(5657)'],
    products = ['[CH2]OO[C](C)OO(5260)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1999.59,'m^3/(mol*s)'), n=1.33685, Ea=(2.9843,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-HH_Cds;HJ] for rate rule [Cds-HH_Cds-OsOs;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['OH(5)', '[CH2]OOC(C)=O(5658)'],
    products = ['[CH2]OO[C](C)OO(5260)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(0.0323,'m^3/(mol*s)'), n=2.98, Ea=(308.836,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Od_CO-NdNd;YJ] for rate rule [Od_CO-NdNd;OJ_pri]
Euclidian distance = 2.0
family: R_Addition_MultipleBond
Ea raised from 306.8 to 308.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]OOC([CH2])OO(5659)'],
    products = ['[CH2]OO[C](C)OO(5260)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(309968,'s^-1'), n=2.08546, Ea=(132.285,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_2H;Cs_H_out] for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_OOH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]OO[C](C)OO(5260)'],
    products = ['[CH2]OOC(C)O[O](5660)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(9.31743e+07,'s^-1'), n=1.36595, Ea=(173.706,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_NonDe;XH_out] for rate rule [R3H_SS_O;C_rad_out_NDMustO;O_H_out]
Euclidian distance = 1.73205080757
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2][C](OO)OOC(5661)'],
    products = ['[CH2]OO[C](C)OO(5260)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(91273.5,'s^-1'), n=1.79, Ea=(54.1828,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;Cs_H_out_2H] for rate rule [R5HJ_1;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]OO[C](C)OO(5260)'],
    products = ['COO[C](C)O[O](5662)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1389.4,'s^-1'), n=2.18252, Ea=(73.6483,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;C_rad_out_2H;XH_out] for rate rule [R6HJ_3;C_rad_out_2H;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2][O](167)', 'C[C]([O])OO(5646)'],
    products = ['[CH2]OO[C](C)OO(5260)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]OO[C](C)OO(5260)'],
    products = ['C=C(OO)OOC(5663)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(9.63e+09,'s^-1'), n=0.137, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_NDe] for rate rule [R5radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]OO[C](C)OO(5260)'],
    products = ['CC1(COO1)OO(5261)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_NDMustO;Cpri_rad_out_2H]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]OO[C](C)OO(5260)'],
    products = ['[CH2]OOC(C)([O])O(5664)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.04257e+10,'s^-1'), n=0, Ea=(86.1903,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnOOH;C_rad_out_NonDe] for rate rule [ROOH;C_rad_out_NDMustO]
Euclidian distance = 1.41421356237
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C[C]OO(178)', '[CH2]O[O](92)'],
    products = ['[CH2]OO[C](C)OO(5260)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['CH2(19)', 'C[C](O[O])OO(5665)'],
    products = ['[CH2]OO[C](C)OO(5260)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

network(
    label = '1246',
    isomers = [
        '[CH2]OO[C](C)OO(5260)',
    ],
    reactants = [
        ('CH3C(O)OOH(67)', 'CH2O(15)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '1246',
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

