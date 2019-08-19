species(
    label = '[CH2]C(C)C([CH2])(O)C(=O)CC(6366)',
    structure = SMILES('[CH2]C(C)C([CH2])(O)C(=O)CC'),
    E0 = (-123.961,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1615.53,2892.47,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152251,'amu*angstrom^2'), symmetry=1, barrier=(3.50056,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152251,'amu*angstrom^2'), symmetry=1, barrier=(3.50056,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152251,'amu*angstrom^2'), symmetry=1, barrier=(3.50056,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152251,'amu*angstrom^2'), symmetry=1, barrier=(3.50056,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152251,'amu*angstrom^2'), symmetry=1, barrier=(3.50056,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152251,'amu*angstrom^2'), symmetry=1, barrier=(3.50056,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152251,'amu*angstrom^2'), symmetry=1, barrier=(3.50056,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152251,'amu*angstrom^2'), symmetry=1, barrier=(3.50056,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.07377,0.136772,-0.000164485,1.08758e-07,-2.88992e-11,-14693.2,41.7753], Tmin=(100,'K'), Tmax=(917.919,'K')), NASAPolynomial(coeffs=[18.2742,0.0481018,-1.95875e-05,3.5222e-09,-2.3776e-13,-18428.7,-54.6492], Tmin=(917.919,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-123.961,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(Isobutyl) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = 'C3H6(72)',
    structure = SMILES('C=CC'),
    E0 = (5.9763,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.497558,'amu*angstrom^2'), symmetry=1, barrier=(11.4398,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2218.31,'J/mol'), sigma=(4.982,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31912,0.00817959,3.34736e-05,-4.36194e-08,1.58213e-11,749.325,9.54025], Tmin=(100,'K'), Tmax=(983.754,'K')), NASAPolynomial(coeffs=[5.36755,0.0170743,-6.35108e-06,1.1662e-09,-8.2762e-14,-487.138,-4.54468], Tmin=(983.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.9763,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""C3H6""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2]C(C)C1(O)CC1([O])CC(6724)',
    structure = SMILES('[CH2]C(C)C1(O)CC1([O])CC'),
    E0 = (-32.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.61862,0.109306,-7.5336e-05,1.28159e-08,5.25687e-12,-3727.52,39.422], Tmin=(100,'K'), Tmax=(978.759,'K')), NASAPolynomial(coeffs=[23.1286,0.0385309,-1.34014e-05,2.32435e-09,-1.58856e-13,-10026.1,-86.8669], Tmin=(978.759,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-32.78,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(557.07,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(Isobutyl) + radical(CC(C)2OJ)"""),
)

species(
    label = '[CH2]C1(O)C(C)CC1([O])CC(6725)',
    structure = SMILES('[CH2]C1(O)C(C)CC1([O])CC'),
    E0 = (-30.2655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.68847,0.109483,-7.41333e-05,1.37695e-08,3.72372e-12,-3421.71,37.0526], Tmin=(100,'K'), Tmax=(1029.75,'K')), NASAPolynomial(coeffs=[23.5548,0.0402161,-1.517e-05,2.74552e-09,-1.91097e-13,-10146.9,-92.8828], Tmin=(1029.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-30.2655,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(561.227,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CJC(C)2O) + radical(CC(C)2OJ)"""),
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
    label = '[CH2]C(O)(C(=C)C)C(=O)CC(6726)',
    structure = SMILES('[CH2]C(O)(C(=C)C)C(=O)CC'),
    E0 = (-210.365,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1721.64,2794.66,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155285,'amu*angstrom^2'), symmetry=1, barrier=(3.57031,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155285,'amu*angstrom^2'), symmetry=1, barrier=(3.57031,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155285,'amu*angstrom^2'), symmetry=1, barrier=(3.57031,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155285,'amu*angstrom^2'), symmetry=1, barrier=(3.57031,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155285,'amu*angstrom^2'), symmetry=1, barrier=(3.57031,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155285,'amu*angstrom^2'), symmetry=1, barrier=(3.57031,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155285,'amu*angstrom^2'), symmetry=1, barrier=(3.57031,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.188,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.88352,0.137987,-0.000190306,1.50529e-07,-4.83743e-11,-25097.3,38.3402], Tmin=(100,'K'), Tmax=(801.384,'K')), NASAPolynomial(coeffs=[13.7476,0.0547376,-2.46972e-05,4.61906e-09,-3.16441e-13,-27434.7,-32.5625], Tmin=(801.384,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-210.365,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsHH) + radical(C=CC(O)(C=O)CJ)"""),
)

species(
    label = 'C3H6(T)(143)',
    structure = SMILES('[CH2][CH]C'),
    E0 = (284.865,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.238389,'amu*angstrom^2'), symmetry=1, barrier=(5.48103,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00909639,'amu*angstrom^2'), symmetry=1, barrier=(22.1005,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.93778,0.0190991,4.26842e-06,-1.44873e-08,5.74941e-12,34303.2,12.9695], Tmin=(100,'K'), Tmax=(1046.81,'K')), NASAPolynomial(coeffs=[5.93909,0.0171892,-6.69152e-06,1.21546e-09,-8.39795e-14,33151.2,-4.14888], Tmin=(1046.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(284.865,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""C3H6(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.91313,0.0225012,-8.59328e-06,4.80319e-10,2.48743e-13,-5274.24,14.655], Tmin=(100,'K'), Tmax=(1673.18,'K')), NASAPolynomial(coeffs=[7.46306,0.0158512,-6.42141e-06,1.12499e-09,-7.32055e-14,-7388.54,-11.406], Tmin=(1673.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-44.1874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""CH3CH2CO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C(C)C(=C)O(2940)',
    structure = SMILES('[CH2]C(C)C(=C)O'),
    E0 = (-57.7502,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,487.358],'cm^-1')),
        HinderedRotor(inertia=(0.0819617,'amu*angstrom^2'), symmetry=1, barrier=(13.8096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0818807,'amu*angstrom^2'), symmetry=1, barrier=(13.808,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.081887,'amu*angstrom^2'), symmetry=1, barrier=(13.8087,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.08198,'amu*angstrom^2'), symmetry=1, barrier=(13.8114,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00656,0.0529817,-7.49544e-06,-3.6877e-08,2.10514e-11,-6826.03,24.6616], Tmin=(100,'K'), Tmax=(917.958,'K')), NASAPolynomial(coeffs=[16.9638,0.0171451,-3.99897e-06,5.73081e-10,-3.87537e-14,-11175.4,-58.6898], Tmin=(917.958,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-57.7502,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Isobutyl)"""),
)

species(
    label = 'CH3(17)',
    structure = SMILES('[CH3]'),
    E0 = (136.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([604.263,1333.71,1492.19,2836.77,2836.77,3806.92],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (15.0345,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.65718,0.0021266,5.45839e-06,-6.6181e-09,2.46571e-12,16422.7,1.67354], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.97812,0.00579785,-1.97558e-06,3.07298e-10,-1.79174e-14,16509.5,4.72248], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(136.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH3""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = '[CH2]C(C)C(=C)C(=O)CC(6727)',
    structure = SMILES('[CH2]C(C)C(=C)C(=O)CC'),
    E0 = (-44.5692,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.188,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.289552,0.0878621,-6.74557e-05,3.09751e-08,-6.50866e-12,-5232.11,34.1205], Tmin=(100,'K'), Tmax=(1053.29,'K')), NASAPolynomial(coeffs=[7.86556,0.0590909,-2.64817e-05,5.0407e-09,-3.52996e-13,-6828.03,-2.82256], Tmin=(1053.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-44.5692,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsHH) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(O)([C](C)C)C(=O)CC(6728)',
    structure = SMILES('[CH2]C(O)([C](C)C)C(=O)CC'),
    E0 = (-176.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.655,0.126489,-0.000133821,7.76892e-08,-1.83973e-11,-21022.5,37.9932], Tmin=(100,'K'), Tmax=(1016.96,'K')), NASAPolynomial(coeffs=[17.7288,0.0502464,-2.13647e-05,3.96852e-09,-2.74489e-13,-24965,-55.8486], Tmin=(1016.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-176.47,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CCJ(C)CO) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = '[CH2]C(C)C(C)([O])C(=O)CC(6729)',
    structure = SMILES('[CH2]C(C)C(C)([O])C(=O)CC'),
    E0 = (-93.5899,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,180,180,180,412.961,743.296,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.16268,'amu*angstrom^2'), symmetry=1, barrier=(3.74033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16268,'amu*angstrom^2'), symmetry=1, barrier=(3.74033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16268,'amu*angstrom^2'), symmetry=1, barrier=(3.74033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16268,'amu*angstrom^2'), symmetry=1, barrier=(3.74033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16268,'amu*angstrom^2'), symmetry=1, barrier=(3.74033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16268,'amu*angstrom^2'), symmetry=1, barrier=(3.74033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16268,'amu*angstrom^2'), symmetry=1, barrier=(3.74033,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.803147,0.116097,-0.000115515,4.99153e-08,2.7056e-12,-11093,40.0897], Tmin=(100,'K'), Tmax=(611.632,'K')), NASAPolynomial(coeffs=[10.7948,0.0599046,-2.59119e-05,4.79399e-09,-3.28462e-13,-12879.4,-13.1674], Tmin=(611.632,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-93.5899,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CC(C)(C=O)OJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C](C)C(C)(O)C(=O)CC(6730)',
    structure = SMILES('[CH2][C](C)C(C)(O)C(=O)CC'),
    E0 = (-182.999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.26006,0.117927,-0.000119206,6.87668e-08,-1.64499e-11,-21822.1,39.4331], Tmin=(100,'K'), Tmax=(1000.52,'K')), NASAPolynomial(coeffs=[14.947,0.0531329,-2.20661e-05,4.04054e-09,-2.7686e-13,-25065.2,-38.7649], Tmin=(1000.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-182.999,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(Isobutyl) + radical(CCJ(C)CO)"""),
)

species(
    label = '[CH2]C([O])(C(=O)CC)C(C)C(6731)',
    structure = SMILES('[CH2]C([O])(C(=O)CC)C(C)C'),
    E0 = (-87.0609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.63466,0.131083,-0.00016042,1.14368e-07,-3.36604e-11,-10274.6,40.1456], Tmin=(100,'K'), Tmax=(822.505,'K')), NASAPolynomial(coeffs=[13.6192,0.0568968,-2.51187e-05,4.69632e-09,-3.23706e-13,-12783.7,-30.464], Tmin=(822.505,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-87.0609,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CC(C)(C=O)OJ) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = '[CH2]C([CH2])C(C)(O)C(=O)CC(6732)',
    structure = SMILES('[CH2]C([CH2])C(C)(O)C(=O)CC'),
    E0 = (-130.49,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.69967,0.128483,-0.000150948,1.01389e-07,-2.76768e-11,-15492,41.9017], Tmin=(100,'K'), Tmax=(891.288,'K')), NASAPolynomial(coeffs=[15.6172,0.0507682,-2.01586e-05,3.56284e-09,-2.37491e-13,-18578.9,-39.6493], Tmin=(891.288,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-130.49,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C)C(C)(O)C(=O)[CH]C(6733)',
    structure = SMILES('[CH2]C(C)C(C)(O)C([O])=CC'),
    E0 = (-199.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.87209,0.117289,-0.000109974,5.5182e-08,-1.10236e-11,-23724.8,44.3043], Tmin=(100,'K'), Tmax=(1218.14,'K')), NASAPolynomial(coeffs=[22.0967,0.0385828,-1.30558e-05,2.13989e-09,-1.37701e-13,-29564.2,-76.0607], Tmin=(1218.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-199.1,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(O)(C(=O)[CH]C)C(C)C(6734)',
    structure = SMILES('[CH2]C(O)(C([O])=CC)C(C)C'),
    E0 = (-190.739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.74163,0.118371,-0.000111836,5.63168e-08,-1.13736e-11,-22727.2,42.2369], Tmin=(100,'K'), Tmax=(1197.05,'K')), NASAPolynomial(coeffs=[21.235,0.0415937,-1.56275e-05,2.73562e-09,-1.83351e-13,-28228.1,-72.7447], Tmin=(1197.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-190.739,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=CC(C)(O)CJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]CC(=O)C(C)(O)C([CH2])C(6735)',
    structure = SMILES('[CH2]CC(=O)C(C)(O)C([CH2])C'),
    E0 = (-123.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.75452,0.130569,-0.000152994,1.01021e-07,-2.71258e-11,-14708.1,41.6481], Tmin=(100,'K'), Tmax=(904.297,'K')), NASAPolynomial(coeffs=[16.21,0.0511053,-2.11842e-05,3.84779e-09,-2.61343e-13,-17957.1,-43.2132], Tmin=(904.297,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-123.984,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJCC=O) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]CC(=O)C([CH2])(O)C(C)C(6736)',
    structure = SMILES('[CH2]CC(=O)C([CH2])(O)C(C)C'),
    E0 = (-117.455,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1689.58,2817.88,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154335,'amu*angstrom^2'), symmetry=1, barrier=(3.54847,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154335,'amu*angstrom^2'), symmetry=1, barrier=(3.54847,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154335,'amu*angstrom^2'), symmetry=1, barrier=(3.54847,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154335,'amu*angstrom^2'), symmetry=1, barrier=(3.54847,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154335,'amu*angstrom^2'), symmetry=1, barrier=(3.54847,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154335,'amu*angstrom^2'), symmetry=1, barrier=(3.54847,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154335,'amu*angstrom^2'), symmetry=1, barrier=(3.54847,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154335,'amu*angstrom^2'), symmetry=1, barrier=(3.54847,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.12945,0.138875,-0.000166634,1.08574e-07,-2.84473e-11,-13909.3,40.1379], Tmin=(100,'K'), Tmax=(929.275,'K')), NASAPolynomial(coeffs=[18.8915,0.0483941,-2.05859e-05,3.80045e-09,-2.6104e-13,-17816.2,-59.7345], Tmin=(929.275,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-117.455,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJCC=O) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = '[CH2]C(C)C1(O)CO[C]1CC(6737)',
    structure = SMILES('[CH2]C(C)C1(O)CO[C]1CC'),
    E0 = (-41.573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.11872,0.108664,-8.7149e-05,2.73978e-08,2.34302e-12,-4811.74,38.3096], Tmin=(100,'K'), Tmax=(828.504,'K')), NASAPolynomial(coeffs=[16.702,0.0455473,-1.43742e-05,2.23052e-09,-1.38744e-13,-8551.31,-49.0597], Tmin=(828.504,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-41.573,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(557.07,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(Isobutyl) + radical(C2CsJOCs)"""),
)

species(
    label = '[CH2]C1(O)[C](CC)OCC1C(6738)',
    structure = SMILES('[CH2]C1(O)[C](CC)OCC1C'),
    E0 = (-113.102,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.20557,0.107225,-8.23648e-05,2.63738e-08,2.9167e-13,-13408.8,33.7729], Tmin=(100,'K'), Tmax=(906.031,'K')), NASAPolynomial(coeffs=[17.4854,0.0456722,-1.51685e-05,2.46947e-09,-1.59507e-13,-17656.2,-59.3048], Tmin=(906.031,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-113.102,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(561.227,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(C2CsJOCs) + radical(CJC(C)2O)"""),
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
    label = '[CH2]C(=C=O)C([CH2])C(6739)',
    structure = SMILES('[CH2]C(C)C(=C)[C]=O'),
    E0 = (189.606,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.0339455,'amu*angstrom^2'), symmetry=1, barrier=(18.4667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179988,'amu*angstrom^2'), symmetry=1, barrier=(4.13829,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180374,'amu*angstrom^2'), symmetry=1, barrier=(4.14715,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.803183,'amu*angstrom^2'), symmetry=1, barrier=(18.4667,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22132,0.0618244,-4.8147e-05,2.00338e-08,-3.51982e-12,22903.6,27.0174], Tmin=(100,'K'), Tmax=(1299.81,'K')), NASAPolynomial(coeffs=[10.7294,0.0325648,-1.43813e-05,2.7156e-09,-1.88941e-13,20431.8,-21.3467], Tmin=(1299.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(189.606,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cd-CdCs(CO)) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)CJ=O) + radical(Isobutyl)"""),
)

species(
    label = 'C=C(C)C(C)(O)C(=O)CC(6740)',
    structure = SMILES('C=C(C)C(C)(O)C(=O)CC'),
    E0 = (-423.588,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.35053,0.12445,-0.000141515,9.41101e-08,-2.61199e-11,-50759.1,36.3897], Tmin=(100,'K'), Tmax=(865.691,'K')), NASAPolynomial(coeffs=[13.1324,0.0575317,-2.55668e-05,4.82081e-09,-3.34939e-13,-53266.8,-31.3936], Tmin=(865.691,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-423.588,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsHH)"""),
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
    label = '[CH2]CC([CH2])(O)C(=O)CC(6318)',
    structure = SMILES('[CH2]CC([CH2])(O)C(=O)CC'),
    E0 = (-94.3227,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,3615,1277.5,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,180,1600,1714.53,2796.77,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154963,'amu*angstrom^2'), symmetry=1, barrier=(3.5629,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154963,'amu*angstrom^2'), symmetry=1, barrier=(3.5629,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154963,'amu*angstrom^2'), symmetry=1, barrier=(3.5629,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154963,'amu*angstrom^2'), symmetry=1, barrier=(3.5629,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154963,'amu*angstrom^2'), symmetry=1, barrier=(3.5629,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154963,'amu*angstrom^2'), symmetry=1, barrier=(3.5629,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154963,'amu*angstrom^2'), symmetry=1, barrier=(3.5629,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (128.169,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4369.31,'J/mol'), sigma=(7.38606,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=682.48 K, Pc=24.6 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.33096,0.122451,-0.000154591,1.07306e-07,-3.00235e-11,-11157.1,37.1524], Tmin=(100,'K'), Tmax=(871.218,'K')), NASAPolynomial(coeffs=[15.83,0.0436639,-1.89463e-05,3.51334e-09,-2.41133e-13,-14147.4,-43.2742], Tmin=(871.218,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-94.3227,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(RCCJ) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = '[CH2]C(C)C([CH2])(O)C(C)=O(6741)',
    structure = SMILES('[CH2]C(C)C([CH2])(O)C(C)=O'),
    E0 = (-98.5444,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,375,552.5,462.5,1710,3615,1277.5,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1390,370,380,2900,435,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.36321,0.115504,-0.000128974,7.67776e-08,-1.8142e-11,-11656.9,37.8857], Tmin=(100,'K'), Tmax=(1034.81,'K')), NASAPolynomial(coeffs=[19.261,0.0357807,-1.34086e-05,2.32438e-09,-1.54416e-13,-15925.2,-62.3193], Tmin=(1034.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-98.5444,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(Isobutyl) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = '[CH2]C(C)C[C](O)C(=O)CC(6365)',
    structure = SMILES('[CH2]C(C)CC(O)=C([O])CC'),
    E0 = (-218.025,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
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
    molecularWeight = (142.196,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4737.91,'J/mol'), sigma=(8.02363,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=740.05 K, Pc=20.81 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.70705,0.11192,-8.20779e-05,1.6683e-08,5.02528e-12,-26004.6,42.8861], Tmin=(100,'K'), Tmax=(954.452,'K')), NASAPolynomial(coeffs=[23.9193,0.0358577,-1.17854e-05,1.98211e-09,-1.3376e-13,-32323.7,-87.0284], Tmin=(954.452,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-218.025,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(Isobutyl) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C(C)[C](O)CC(=O)CC(6742)',
    structure = SMILES('[CH2]C(C)[C](O)CC(=O)CC'),
    E0 = (-151.406,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.58463,0.133249,-0.000184121,1.53578e-07,-5.19165e-11,-18018.9,42.378], Tmin=(100,'K'), Tmax=(834.73,'K')), NASAPolynomial(coeffs=[9.60319,0.0630511,-2.81709e-05,5.22169e-09,-3.54871e-13,-19308.8,-6.11448], Tmin=(834.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-151.406,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(Isobutyl) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C(O)(C[CH]C)C(=O)CC(6362)',
    structure = SMILES('[CH2]C(O)(C[CH]C)C(=O)CC'),
    E0 = (-128.903,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.4624,0.131445,-0.000152452,8.59456e-08,-1.13099e-11,-15317.2,40.8806], Tmin=(100,'K'), Tmax=(630.021,'K')), NASAPolynomial(coeffs=[13.6745,0.056618,-2.49574e-05,4.6402e-09,-3.18008e-13,-17646.8,-28.5042], Tmin=(630.021,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-128.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(RCCJC) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = '[CH2]C(O)([CH]CC)C(=O)CC(6743)',
    structure = SMILES('[CH2]C(O)([CH]CC)C(=O)CC'),
    E0 = (-123.447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.78751,0.130909,-0.000148081,9.20269e-08,-2.32161e-11,-14641.9,41.7464], Tmin=(100,'K'), Tmax=(959.638,'K')), NASAPolynomial(coeffs=[17.6694,0.0498083,-2.13129e-05,3.96063e-09,-2.7358e-13,-18376.2,-51.3204], Tmin=(959.638,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-123.447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJC(C)(C=O)O) + radical(CCJCO)"""),
)

species(
    label = 'CCC(=O)C1(O)CCC1C(6368)',
    structure = SMILES('CCC(=O)C1(O)CCC1C'),
    E0 = (-385.825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.15869,0.0989075,-5.94141e-05,8.74462e-09,3.13339e-12,-46206,35.9606], Tmin=(100,'K'), Tmax=(1086.57,'K')), NASAPolynomial(coeffs=[19.9739,0.0448652,-1.76006e-05,3.20923e-09,-2.22144e-13,-52200.6,-74.1987], Tmin=(1086.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-385.825,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(561.227,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutane)"""),
)

species(
    label = '[CH2]C(=C(CC)OO)C([CH2])C(6744)',
    structure = SMILES('[CH2]C(=C(CC)OO)C([CH2])C'),
    E0 = (124.702,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1310,387.5,850,1000,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1390,370,380,2900,435,200,800],'cm^-1')),
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
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.30219,0.11181,-9.63914e-05,4.48276e-08,-8.55674e-12,15193,42.7468], Tmin=(100,'K'), Tmax=(1241.26,'K')), NASAPolynomial(coeffs=[18.4347,0.0482061,-1.95285e-05,3.54471e-09,-2.41907e-13,10293.3,-56.7375], Tmin=(1241.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(124.702,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsCsOs) + radical(Allyl_P) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C)O[C](CC)C(=C)O(6364)',
    structure = SMILES('[CH2]C(O)=C(CC)OC([CH2])C'),
    E0 = (-167.525,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1200,1600],'cm^-1')),
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
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.12854,0.136861,-0.000140679,7.21528e-08,-1.4257e-11,-19874.9,43.3956], Tmin=(100,'K'), Tmax=(1301.77,'K')), NASAPolynomial(coeffs=[31.3026,0.026218,-7.60391e-06,1.14284e-09,-7.06509e-14,-28428.7,-130.218], Tmin=(1301.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-167.525,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(CJC(C)OC) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]C(C)C([CH2])(O)C(=C)OC(6745)',
    structure = SMILES('[CH2]C(C)C([CH2])(O)C(=C)OC'),
    E0 = (-68.7521,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,180,342.353,791.203,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154039,'amu*angstrom^2'), symmetry=1, barrier=(3.54167,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154039,'amu*angstrom^2'), symmetry=1, barrier=(3.54167,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154039,'amu*angstrom^2'), symmetry=1, barrier=(3.54167,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154039,'amu*angstrom^2'), symmetry=1, barrier=(3.54167,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154039,'amu*angstrom^2'), symmetry=1, barrier=(3.54167,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154039,'amu*angstrom^2'), symmetry=1, barrier=(3.54167,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154039,'amu*angstrom^2'), symmetry=1, barrier=(3.54167,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154039,'amu*angstrom^2'), symmetry=1, barrier=(3.54167,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.0039,0.120729,-0.000107984,4.47454e-08,-5.40532e-12,-8042.78,43.9793], Tmin=(100,'K'), Tmax=(985.713,'K')), NASAPolynomial(coeffs=[23.8818,0.0369083,-1.27253e-05,2.16091e-09,-1.44552e-13,-14177,-85.7615], Tmin=(985.713,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-68.7521,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Isobutyl) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]C(C)C([CH2])(O)C(O)=CC(6746)',
    structure = SMILES('[CH2]C(C)C([CH2])(O)C(O)=CC'),
    E0 = (-123.462,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.10142,0.124193,-0.000117561,5.30125e-08,-7.67636e-12,-14620.4,44.5008], Tmin=(100,'K'), Tmax=(964.677,'K')), NASAPolynomial(coeffs=[24.0863,0.0359813,-1.20785e-05,2.00981e-09,-1.3266e-13,-20621,-85.8115], Tmin=(964.677,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-123.462,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=CC(C)(O)CJ) + radical(Isobutyl)"""),
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
    label = '[CH2]C(C)[C](O)C(=O)CC(6747)',
    structure = SMILES('[CH2]C(C)C(O)=C([O])CC'),
    E0 = (-193.175,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.8946,0.106443,-9.91443e-05,4.69297e-08,-8.51577e-12,-23001.6,41.3639], Tmin=(100,'K'), Tmax=(1504.4,'K')), NASAPolynomial(coeffs=[25.1414,0.0236738,-5.76464e-06,7.39832e-10,-4.07833e-14,-29904.5,-96.0169], Tmin=(1504.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-193.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(O)([CH]C)C(=O)CC(6379)',
    structure = SMILES('[CH2]C(O)([CH]C)C(=O)CC'),
    E0 = (-99.667,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1754.16,2755.36,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156202,'amu*angstrom^2'), symmetry=1, barrier=(3.5914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156202,'amu*angstrom^2'), symmetry=1, barrier=(3.5914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156202,'amu*angstrom^2'), symmetry=1, barrier=(3.5914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156202,'amu*angstrom^2'), symmetry=1, barrier=(3.5914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156202,'amu*angstrom^2'), symmetry=1, barrier=(3.5914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156202,'amu*angstrom^2'), symmetry=1, barrier=(3.5914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156202,'amu*angstrom^2'), symmetry=1, barrier=(3.5914,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (128.169,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.16691,0.116456,-0.000135482,8.5667e-08,-2.18244e-11,-11803.5,37.2714], Tmin=(100,'K'), Tmax=(953.997,'K')), NASAPolynomial(coeffs=[16.727,0.0414301,-1.75175e-05,3.23263e-09,-2.22295e-13,-15217.7,-48.2139], Tmin=(953.997,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-99.667,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CCJCO) + radical(CJC(C)(C=O)O)"""),
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
    E0 = (-123.961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-32.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-30.2655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (1.4276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-57.9958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-8.00671,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-123.961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-5.37377,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-1.7926,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-21.2817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (45.0881,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-6.0221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-40.2813,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-41.1181,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-37.3525,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-77.5189,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-71.2429,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-67.6649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (99.9221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (62.4639,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-86.093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (270.128,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-60.5611,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (325.539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (321.318,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (33.3571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (33.3571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (35.9739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (70.8039,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-115.677,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (267.815,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-24.411,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (190.237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (19.6518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (188.388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (281.896,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    products = ['C3H6(72)', 'C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    products = ['[CH2]C(C)C1(O)CC1([O])CC(6724)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.34238e+09,'s^-1'), n=0.889391, Ea=(91.1813,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonylbond_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 87.3 to 91.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    products = ['[CH2]C1(O)C(C)CC1([O])CC(6725)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.46159e+06,'s^-1'), n=1.55572, Ea=(93.6958,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5_SS_CO;carbonylbond_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 89.8 to 93.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[CH2]C(O)(C(=C)C)C(=O)CC(6726)'],
    products = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.0051739,'m^3/(mol*s)'), n=2.82163, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 102 used for Cds-CsCs_Cds-HH;HJ
Exact match found for rate rule [Cds-CsCs_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -4.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['C3H6(T)(143)', 'C=C(O)C(=O)CC(4626)'],
    products = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00392164,'m^3/(mol*s)'), n=2.41519, Ea=(15.6067,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;CJ] for rate rule [Cds-COOs_Cds;CJ]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C2H5CO(71)', '[CH2]C(C)C(=C)O(2940)'],
    products = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=0, Ea=(93.9308,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_R;CO_rad] for rate rule [Cds-OsCs_Cds;CO_rad/NonDe]
Euclidian distance = 2.82842712475
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C3H6(72)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00168615,'m^3/(mol*s)'), n=2.52599, Ea=(55.0055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 54.0 to 55.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH3(17)', '[CH2]C(O)(C=C)C(=O)CC(6377)'],
    products = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(10000,'cm^3/(mol*s)'), n=2.41, Ea=(29.7482,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;CsJ-HHH] for rate rule [Cds-Cs\O2s/H_Cds-HH;CsJ-HHH]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['OH(5)', '[CH2]C(C)C(=C)C(=O)CC(6727)'],
    products = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.52832,'m^3/(mol*s)'), n=2.02802, Ea=(14.4047,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDeCs_Cds;YJ] for rate rule [Cds-COCs_Cds;OJ_pri]
Euclidian distance = 2.2360679775
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    products = ['[CH2]C(O)([C](C)C)C(=O)CC(6728)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.265e-07,'s^-1'), n=5.639, Ea=(102.68,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 38 used for R2H_S;C_rad_out_2H;Cs_H_out_Cs2
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(C)C(C)([O])C(=O)CC(6729)'],
    products = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    products = ['[CH2][C](C)C(C)(O)C(=O)CC(6730)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    products = ['[CH2]C([O])(C(=O)CC)C(C)C(6731)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    products = ['[CH2]C([CH2])C(C)(O)C(=O)CC(6732)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(114000,'s^-1'), n=1.74, Ea=(82.8432,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 109 used for R4H_SSS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    products = ['[CH2]C(C)C(C)(O)C(=O)[CH]C(6733)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(6.44e+09,'s^-1'), n=0.13, Ea=(86.6088,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 131 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    products = ['[CH2]C(O)(C(=O)[CH]C)C(C)C(6734)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(262000,'s^-1'), n=1.62, Ea=(46.4424,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_2H;Cs_H_out_H/NonDeC] for rate rule [R5H_CC(O2d)CC;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    products = ['[CH2]CC(=O)C(C)(O)C([CH2])C(6735)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(68850,'s^-1'), n=1.68, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_2H;Cs_H_out_2H] for rate rule [R5H_CCC(O2d)C;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]CC(=O)C([CH2])(O)C(C)C(6736)'],
    products = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7380,'s^-1'), n=1.79, Ea=(49.7896,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 113 used for R6H_SSSSS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R6H_SSSSS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C3H6(T)(143)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    products = ['[CH2]C(C)C1(O)CO[C]1CC(6737)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.01137e+09,'s^-1'), n=0.572544, Ea=(186.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4_S;multiplebond_intra;radadd_intra_cs2H] + [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    products = ['[CH2]C1(O)[C](CC)OCC1C(6738)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.64784e+07,'s^-1'), n=0.990488, Ea=(37.8683,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5_SS_CO;carbonyl_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CH3CH2OH(54)', '[CH2]C(=C=O)C([CH2])C(6739)'],
    products = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.79e-11,'m^3/(mol*s)'), n=3.97, Ea=(329.281,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [doublebond;Cs_OH] for rate rule [Cdd_Cd;C_pri_OH]
Euclidian distance = 1.41421356237
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    products = ['C=C(C)C(C)(O)C(=O)CC(6740)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CH2(S)(23)', '[CH2]CC([CH2])(O)C(=O)CC(6318)'],
    products = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['CH2(S)(23)', '[CH2]C(C)C([CH2])(O)C(C)=O(6741)'],
    products = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;C_pri/De]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    products = ['[CH2]C(C)C[C](O)C(=O)CC(6365)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    products = ['[CH2]C(C)[C](O)CC(=O)CC(6742)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ-HH;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;CO]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    products = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    products = ['[CH2]C(O)([CH]CC)C(=O)CC(6743)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(5.59192e+09,'s^-1'), n=1.025, Ea=(194.765,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CH3] for rate rule [cCs(-HC)CJ;CsJ-HH;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    products = ['CCC(=O)C1(O)CCC1C(6368)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(=C(CC)OO)C([CH2])C(6744)'],
    products = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O] for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O_R]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C(C)O[C](CC)C(=C)O(6364)'],
    products = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O] for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O_C]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(C)C([CH2])(O)C(=C)OC(6745)'],
    products = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(37989.5,'s^-1'), n=2.515, Ea=(258.99,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C] + [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O] for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_C]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C(C)C([CH2])(O)C(O)=CC(6746)'],
    products = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O_H] for rate rule [R_ROR;R1_doublebond_CHCH3;R2_doublebond_CsC;R_O_H]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction35',
    reactants = ['CH2(19)', '[CH2]C(C)[C](O)C(=O)CC(6747)'],
    products = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/ODMustO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction36',
    reactants = ['CH2(19)', '[CH2]C(O)([CH]C)C(=O)CC(6379)'],
    products = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '1403',
    isomers = [
        '[CH2]C(C)C([CH2])(O)C(=O)CC(6366)',
    ],
    reactants = [
        ('C3H6(72)', 'C=C(O)C(=O)CC(4626)'),
        ('C3H6(72)', '[CH2]C(O)=C([O])CC(4557)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '1403',
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

