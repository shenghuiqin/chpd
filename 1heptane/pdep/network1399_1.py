species(
    label = '[CH2]C(O)(C[CH]C)C(=O)CC(6362)',
    structure = SMILES('[CH2]C(O)(C[CH]C)C(=O)CC'),
    E0 = (-128.903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,180,180,180,560.729,594.512,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.159113,'amu*angstrom^2'), symmetry=1, barrier=(3.65831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159113,'amu*angstrom^2'), symmetry=1, barrier=(3.65831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159113,'amu*angstrom^2'), symmetry=1, barrier=(3.65831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159113,'amu*angstrom^2'), symmetry=1, barrier=(3.65831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159113,'amu*angstrom^2'), symmetry=1, barrier=(3.65831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159113,'amu*angstrom^2'), symmetry=1, barrier=(3.65831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159113,'amu*angstrom^2'), symmetry=1, barrier=(3.65831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159113,'amu*angstrom^2'), symmetry=1, barrier=(3.65831,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.4624,0.131445,-0.000152452,8.59456e-08,-1.13099e-11,-15317.2,40.8806], Tmin=(100,'K'), Tmax=(630.021,'K')), NASAPolynomial(coeffs=[13.6745,0.056618,-2.49574e-05,4.6402e-09,-3.18008e-13,-17646.8,-28.5042], Tmin=(630.021,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-128.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(RCCJC) + radical(CJC(C)(C=O)O)"""),
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
    label = 'C[CH]CC1(O)CC1([O])CC(6840)',
    structure = SMILES('C[CH]CC1(O)CC1([O])CC'),
    E0 = (-34.3745,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.49153,0.110748,-9.3279e-05,4.1654e-08,-7.53047e-12,-3928.31,40.2045], Tmin=(100,'K'), Tmax=(1319.9,'K')), NASAPolynomial(coeffs=[20.7639,0.0433016,-1.66295e-05,2.93897e-09,-1.97487e-13,-9803.26,-73.342], Tmin=(1319.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-34.3745,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(557.07,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(RCCJC) + radical(CC(C)2OJ)"""),
)

species(
    label = '[CH2]C1(O)CC(C)C1([O])CC(6781)',
    structure = SMILES('[CH2]C1(O)CC(C)C1([O])CC'),
    E0 = (-30.2655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.68847,0.109483,-7.41333e-05,1.37695e-08,3.72372e-12,-3421.71,37.0526], Tmin=(100,'K'), Tmax=(1029.75,'K')), NASAPolynomial(coeffs=[23.5548,0.0402161,-1.517e-05,2.74552e-09,-1.91097e-13,-10146.9,-92.8828], Tmin=(1029.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-30.2655,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(561.227,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CC(C)2OJ) + radical(CJC(C)2O)"""),
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
    label = '[CH2]C(O)(C=CC)C(=O)CC(6841)',
    structure = SMILES('[CH2]C(O)(C=CC)C(=O)CC'),
    E0 = (-207.335,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1745.98,2771.42,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15609,'amu*angstrom^2'), symmetry=1, barrier=(3.58882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15609,'amu*angstrom^2'), symmetry=1, barrier=(3.58882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15609,'amu*angstrom^2'), symmetry=1, barrier=(3.58882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15609,'amu*angstrom^2'), symmetry=1, barrier=(3.58882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15609,'amu*angstrom^2'), symmetry=1, barrier=(3.58882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15609,'amu*angstrom^2'), symmetry=1, barrier=(3.58882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15609,'amu*angstrom^2'), symmetry=1, barrier=(3.58882,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.188,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.69108,0.133604,-0.000180436,1.42028e-07,-4.57921e-11,-24739.7,38.6998], Tmin=(100,'K'), Tmax=(788.098,'K')), NASAPolynomial(coeffs=[12.8864,0.0558587,-2.53097e-05,4.754e-09,-3.26933e-13,-26920.7,-27.4162], Tmin=(788.098,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-207.335,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C=CC(O)(C=O)CJ)"""),
)

species(
    label = '[CH2]C(O)(CC=C)C(=O)CC(6842)',
    structure = SMILES('[CH2]C(O)(CC=C)C(=O)CC'),
    E0 = (-195.567,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,180,180,180,560.97,586.724,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.159235,'amu*angstrom^2'), symmetry=1, barrier=(3.66114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159235,'amu*angstrom^2'), symmetry=1, barrier=(3.66114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159235,'amu*angstrom^2'), symmetry=1, barrier=(3.66114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159235,'amu*angstrom^2'), symmetry=1, barrier=(3.66114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159235,'amu*angstrom^2'), symmetry=1, barrier=(3.66114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159235,'amu*angstrom^2'), symmetry=1, barrier=(3.66114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159235,'amu*angstrom^2'), symmetry=1, barrier=(3.66114,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.188,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.67157,0.128304,-0.000147003,9.1779e-08,-2.31834e-11,-23320,39.3079], Tmin=(100,'K'), Tmax=(959.639,'K')), NASAPolynomial(coeffs=[17.7703,0.0472694,-2.03438e-05,3.79167e-09,-2.62367e-13,-27051.6,-53.6875], Tmin=(959.639,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-195.567,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJC(C)(C=O)O)"""),
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
    label = 'C=C(O)C[CH]C(2383)',
    structure = SMILES('C=C(O)C[CH]C'),
    E0 = (-60.4143,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,433.312,433.344],'cm^-1')),
        HinderedRotor(inertia=(0.000897782,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0761217,'amu*angstrom^2'), symmetry=1, barrier=(10.1393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0760956,'amu*angstrom^2'), symmetry=1, barrier=(10.1388,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0761163,'amu*angstrom^2'), symmetry=1, barrier=(10.1386,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13709,0.0547856,-2.79212e-05,-3.31604e-09,5.58123e-12,-7155.87,25.3906], Tmin=(100,'K'), Tmax=(991.19,'K')), NASAPolynomial(coeffs=[12.8569,0.024853,-8.89992e-06,1.56377e-09,-1.07184e-13,-10332.1,-35.3496], Tmin=(991.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-60.4143,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJC)"""),
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
    label = 'C=C(C[CH]C)C(=O)CC(6843)',
    structure = SMILES('C=C(C[CH]C)C(=O)CC'),
    E0 = (-47.2334,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.117597,0.0934345,-0.000102401,8.63132e-08,-3.31485e-11,-5548.73,35.9244], Tmin=(100,'K'), Tmax=(740.57,'K')), NASAPolynomial(coeffs=[3.13984,0.0679046,-3.20447e-05,6.19227e-09,-4.35064e-13,-5743.92,23.9559], Tmin=(740.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-47.2334,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsHH) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C(O)([CH]CC)C(=O)CC(6743)',
    structure = SMILES('[CH2]C(O)([CH]CC)C(=O)CC'),
    E0 = (-123.447,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,180,180,180,508.593,638.262,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15976,'amu*angstrom^2'), symmetry=1, barrier=(3.6732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15976,'amu*angstrom^2'), symmetry=1, barrier=(3.6732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15976,'amu*angstrom^2'), symmetry=1, barrier=(3.6732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15976,'amu*angstrom^2'), symmetry=1, barrier=(3.6732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15976,'amu*angstrom^2'), symmetry=1, barrier=(3.6732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15976,'amu*angstrom^2'), symmetry=1, barrier=(3.6732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15976,'amu*angstrom^2'), symmetry=1, barrier=(3.6732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15976,'amu*angstrom^2'), symmetry=1, barrier=(3.6732,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.78751,0.130909,-0.000148081,9.20269e-08,-2.32161e-11,-14641.9,41.7464], Tmin=(100,'K'), Tmax=(959.638,'K')), NASAPolynomial(coeffs=[17.6694,0.0498083,-2.13129e-05,3.96063e-09,-2.7358e-13,-18376.2,-51.3204], Tmin=(959.638,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-123.447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJC(C)(C=O)O) + radical(CCJCO)"""),
)

species(
    label = '[CH2]CCC([CH2])(O)C(=O)CC(6844)',
    structure = SMILES('[CH2]CCC([CH2])(O)C(=O)CC'),
    E0 = (-118.103,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,375,552.5,462.5,1710,3615,1277.5,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,560.879,588.379,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.158897,'amu*angstrom^2'), symmetry=1, barrier=(3.65335,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158897,'amu*angstrom^2'), symmetry=1, barrier=(3.65335,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158897,'amu*angstrom^2'), symmetry=1, barrier=(3.65335,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158897,'amu*angstrom^2'), symmetry=1, barrier=(3.65335,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158897,'amu*angstrom^2'), symmetry=1, barrier=(3.65335,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158897,'amu*angstrom^2'), symmetry=1, barrier=(3.65335,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158897,'amu*angstrom^2'), symmetry=1, barrier=(3.65335,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158897,'amu*angstrom^2'), symmetry=1, barrier=(3.65335,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.94869,0.136855,-0.000166946,1.13256e-07,-3.12016e-11,-13995.6,41.6181], Tmin=(100,'K'), Tmax=(881.671,'K')), NASAPolynomial(coeffs=[16.7271,0.0521248,-2.27918e-05,4.25361e-09,-2.93464e-13,-17288.7,-46.1295], Tmin=(881.671,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-118.103,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(RCCJ) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = 'C[CH]CC(C)([O])C(=O)CC(6845)',
    structure = SMILES('C[CH]CC(C)([O])C(=O)CC'),
    E0 = (-98.5317,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,375,552.5,462.5,1710,180,180,180,180,806.536,1417.48,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.170644,'amu*angstrom^2'), symmetry=1, barrier=(3.92345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170644,'amu*angstrom^2'), symmetry=1, barrier=(3.92345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170644,'amu*angstrom^2'), symmetry=1, barrier=(3.92345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170644,'amu*angstrom^2'), symmetry=1, barrier=(3.92345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170644,'amu*angstrom^2'), symmetry=1, barrier=(3.92345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170644,'amu*angstrom^2'), symmetry=1, barrier=(3.92345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170644,'amu*angstrom^2'), symmetry=1, barrier=(3.92345,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.199,0.125326,-0.000170372,1.46026e-07,-5.1124e-11,-11674,42.6575], Tmin=(100,'K'), Tmax=(826.462,'K')), NASAPolynomial(coeffs=[7.02297,0.0668401,-3.02942e-05,5.66497e-09,-3.87338e-13,-12394.6,8.42052], Tmin=(826.462,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-98.5317,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(RCCJC) + radical(CC(C)(C=O)OJ)"""),
)

species(
    label = 'C[CH][CH]C(C)(O)C(=O)CC(6846)',
    structure = SMILES('C[CH][CH]C(C)(O)C(=O)CC'),
    E0 = (-140.612,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.15384,0.120143,-0.000137738,9.57158e-08,-2.8033e-11,-16732.3,42.5935], Tmin=(100,'K'), Tmax=(820.069,'K')), NASAPolynomial(coeffs=[11.3734,0.0590398,-2.59744e-05,4.8598e-09,-3.35577e-13,-18787,-15.3585], Tmin=(820.069,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-140.612,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CCJCO) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C([O])(CCC)C(=O)CC(6847)',
    structure = SMILES('[CH2]C([O])(CCC)C(=O)CC'),
    E0 = (-81.3666,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3100,440,815,1455,1000,375,552.5,462.5,1710,180,180,180,180,950.951,1276.47,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.167278,'amu*angstrom^2'), symmetry=1, barrier=(3.84604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167278,'amu*angstrom^2'), symmetry=1, barrier=(3.84604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167278,'amu*angstrom^2'), symmetry=1, barrier=(3.84604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167278,'amu*angstrom^2'), symmetry=1, barrier=(3.84604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167278,'amu*angstrom^2'), symmetry=1, barrier=(3.84604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167278,'amu*angstrom^2'), symmetry=1, barrier=(3.84604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167278,'amu*angstrom^2'), symmetry=1, barrier=(3.84604,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.29619,0.128279,-0.000144632,7.81997e-08,-7.91218e-12,-9606.39,39.9705], Tmin=(100,'K'), Tmax=(616.332,'K')), NASAPolynomial(coeffs=[12.5696,0.0590156,-2.65023e-05,4.98175e-09,-3.4393e-13,-11709.2,-23.4065], Tmin=(616.332,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-81.3666,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJC(C)(C=O)O) + radical(CC(C)(C=O)OJ)"""),
)

species(
    label = 'C[CH]CC(C)(O)C(=O)[CH]C(6848)',
    structure = SMILES('C[CH]CC(C)(O)C([O])=CC'),
    E0 = (-204.042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.16292,0.111939,-0.000104644,5.477e-08,-1.18171e-11,-24353.3,42.9958], Tmin=(100,'K'), Tmax=(1106.68,'K')), NASAPolynomial(coeffs=[16.2034,0.0491689,-1.95633e-05,3.51657e-09,-2.38764e-13,-28197,-42.5468], Tmin=(1106.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-204.042,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(RCCJC) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C(O)(CCC)C(=O)[CH]C(6849)',
    structure = SMILES('[CH2]C(O)(CCC)C([O])=CC'),
    E0 = (-185.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.60088,0.119082,-0.000115604,6.05966e-08,-1.28377e-11,-22050.6,42.7061], Tmin=(100,'K'), Tmax=(1137.75,'K')), NASAPolynomial(coeffs=[19.5539,0.044708,-1.75511e-05,3.14248e-09,-2.13253e-13,-26864.4,-62.0837], Tmin=(1137.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-185.045,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2][CH]CC(C)(O)C(=O)CC(6850)',
    structure = SMILES('[CH2][CH]CC(C)(O)C(=O)CC'),
    E0 = (-135.268,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.44759,0.127883,-0.000164039,1.284e-07,-4.17973e-11,-16080.5,42.9262], Tmin=(100,'K'), Tmax=(793.059,'K')), NASAPolynomial(coeffs=[10.8975,0.0605077,-2.694e-05,5.02722e-09,-3.44772e-13,-17877.9,-12.7562], Tmin=(793.059,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-135.268,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(RCCJ) + radical(RCCJC)"""),
)

species(
    label = '[CH2]CC(=O)C(C)(O)C[CH]C(6851)',
    structure = SMILES('[CH2]CC(=O)C(C)(O)C[CH]C'),
    E0 = (-128.925,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.978902,0.122509,-0.000126091,4.66244e-08,1.30121e-11,-15339.2,40.2061], Tmin=(100,'K'), Tmax=(580.861,'K')), NASAPolynomial(coeffs=[11.7518,0.0593748,-2.64096e-05,4.93129e-09,-3.38707e-13,-17232,-17.8588], Tmin=(580.861,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-128.925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(RCCJC) + radical(CJCC=O)"""),
)

species(
    label = '[CH2]CC(=O)C([CH2])(O)CCC(6852)',
    structure = SMILES('[CH2]CC(=O)C([CH2])(O)CCC'),
    E0 = (-111.76,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,375,552.5,462.5,1710,3615,1277.5,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,554.277,594.502,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.159001,'amu*angstrom^2'), symmetry=1, barrier=(3.65575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159001,'amu*angstrom^2'), symmetry=1, barrier=(3.65575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159001,'amu*angstrom^2'), symmetry=1, barrier=(3.65575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159001,'amu*angstrom^2'), symmetry=1, barrier=(3.65575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159001,'amu*angstrom^2'), symmetry=1, barrier=(3.65575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159001,'amu*angstrom^2'), symmetry=1, barrier=(3.65575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159001,'amu*angstrom^2'), symmetry=1, barrier=(3.65575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159001,'amu*angstrom^2'), symmetry=1, barrier=(3.65575,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.17854,0.141809,-0.000178098,1.22721e-07,-3.40797e-11,-13224.4,41.2887], Tmin=(100,'K'), Tmax=(877.612,'K')), NASAPolynomial(coeffs=[17.8726,0.0504202,-2.1898e-05,4.06569e-09,-2.79374e-13,-16743.9,-52.829], Tmin=(877.612,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-111.76,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJC(C)(C=O)O) + radical(CJCC=O)"""),
)

species(
    label = 'C[CH]CC1(O)CO[C]1CC(6853)',
    structure = SMILES('C[CH]CC1(O)CO[C]1CC'),
    E0 = (-43.1676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.772498,0.107365,-9.46893e-05,4.13082e-08,-3.20052e-12,-5021.83,38.3163], Tmin=(100,'K'), Tmax=(746.582,'K')), NASAPolynomial(coeffs=[12.2625,0.0537431,-1.95379e-05,3.2959e-09,-2.14373e-13,-7420.13,-23.7876], Tmin=(746.582,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-43.1676,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(557.07,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(RCCJC) + radical(C2CsJOCs)"""),
)

species(
    label = '[CH2]C1(O)CC(C)O[C]1CC(6761)',
    structure = SMILES('[CH2]C1(O)CC(C)O[C]1CC'),
    E0 = (-116.703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.26672,0.108127,-8.35458e-05,2.6869e-08,1.89202e-13,-13839.4,33.9266], Tmin=(100,'K'), Tmax=(912.263,'K')), NASAPolynomial(coeffs=[18.0023,0.0450236,-1.49512e-05,2.43828e-09,-1.57841e-13,-18245,-62.1433], Tmin=(912.263,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-116.703,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(561.227,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(C2CsJOCs) + radical(CJC(C)2O)"""),
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
    label = '[CH2]C(=C=O)C[CH]C(6854)',
    structure = SMILES('C=C([C]=O)C[CH]C'),
    E0 = (186.942,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1855,455,950,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,232.701,232.702],'cm^-1')),
        HinderedRotor(inertia=(0.00311324,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00873855,'amu*angstrom^2'), symmetry=1, barrier=(15.4641,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162232,'amu*angstrom^2'), symmetry=1, barrier=(6.23377,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162231,'amu*angstrom^2'), symmetry=1, barrier=(6.23376,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.12554,0.0506886,3.05529e-06,-1.05541e-07,1.04525e-10,22541.6,25.1744], Tmin=(100,'K'), Tmax=(432.849,'K')), NASAPolynomial(coeffs=[4.25805,0.0444102,-2.1722e-05,4.29255e-09,-3.06648e-13,22231.2,15.2188], Tmin=(432.849,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(186.942,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cd-CdCs(CO)) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)CJ=O) + radical(RCCJC)"""),
)

species(
    label = 'CC=CC(C)(O)C(=O)CC(6855)',
    structure = SMILES('CC=CC(C)(O)C(=O)CC'),
    E0 = (-420.558,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.14284,0.119882,-0.000130963,8.4641e-08,-2.30686e-11,-50402.2,36.6949], Tmin=(100,'K'), Tmax=(876.751,'K')), NASAPolynomial(coeffs=[12.2827,0.0586293,-2.61641e-05,4.95183e-09,-3.45085e-13,-52756.3,-26.3094], Tmin=(876.751,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-420.558,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH)"""),
)

species(
    label = 'C=CCC(C)(O)C(=O)CC(6856)',
    structure = SMILES('C=CCC(C)(O)C(=O)CC'),
    E0 = (-407.178,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.14395,0.116306,-0.00011387,6.26633e-08,-1.43711e-11,-48789.6,37.7952], Tmin=(100,'K'), Tmax=(1035.95,'K')), NASAPolynomial(coeffs=[15.0601,0.0537416,-2.32834e-05,4.36984e-09,-3.04024e-13,-52147.1,-40.953], Tmin=(1035.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-407.178,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH)"""),
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
    label = '[CH2]C(O)(C[CH]C)C(C)=O(6857)',
    structure = SMILES('[CH2]C(O)(C[CH]C)C(C)=O'),
    E0 = (-103.486,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1744.31,2769.16,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155952,'amu*angstrom^2'), symmetry=1, barrier=(3.58564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155952,'amu*angstrom^2'), symmetry=1, barrier=(3.58564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155952,'amu*angstrom^2'), symmetry=1, barrier=(3.58564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155952,'amu*angstrom^2'), symmetry=1, barrier=(3.58564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155952,'amu*angstrom^2'), symmetry=1, barrier=(3.58564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155952,'amu*angstrom^2'), symmetry=1, barrier=(3.58564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155952,'amu*angstrom^2'), symmetry=1, barrier=(3.58564,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (128.169,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.971287,0.113768,-0.000135702,9.12551e-08,-2.50007e-11,-12271.4,37.7237], Tmin=(100,'K'), Tmax=(885.876,'K')), NASAPolynomial(coeffs=[14.2764,0.0449181,-1.91215e-05,3.52055e-09,-2.4096e-13,-14972.9,-33.9898], Tmin=(885.876,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-103.486,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(RCCJC) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = 'C[CH]CC[C](O)C(=O)CC(6361)',
    structure = SMILES('C[CH]CCC(O)=C([O])CC'),
    E0 = (-219.619,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(4722.71,'J/mol'), sigma=(7.99553,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=737.68 K, Pc=20.96 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.62738,0.113965,-0.000102364,4.89442e-08,-9.40486e-12,-26203.3,43.8362], Tmin=(100,'K'), Tmax=(1253.9,'K')), NASAPolynomial(coeffs=[21.1812,0.0412039,-1.53219e-05,2.66568e-09,-1.77871e-13,-31923.2,-71.3624], Tmin=(1253.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-219.619,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(RCCJC) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C[CH]C[C](O)CC(=O)CC(6858)',
    structure = SMILES('C[CH]C[C](O)CC(=O)CC'),
    E0 = (-156.348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.41858,0.134392,-0.000201959,1.83999e-07,-6.62907e-11,-18623.9,43.0121], Tmin=(100,'K'), Tmax=(842.04,'K')), NASAPolynomial(coeffs=[5.3811,0.0708426,-3.30865e-05,6.22581e-09,-4.25249e-13,-18661.2,17.9548], Tmin=(842.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-156.348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(C2CsJOH) + radical(RCCJC)"""),
)

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
    collisionModel = TransportData(shapeIndex=2, epsilon=(4524.8,'J/mol'), sigma=(7.72749,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=706.76 K, Pc=22.25 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.07377,0.136772,-0.000164485,1.08758e-07,-2.88992e-11,-14693.2,41.7753], Tmin=(100,'K'), Tmax=(917.919,'K')), NASAPolynomial(coeffs=[18.2742,0.0481018,-1.95875e-05,3.5222e-09,-2.3776e-13,-18428.7,-54.6492], Tmin=(917.919,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-123.961,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(Isobutyl) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = 'CCC(=O)C1(O)CC(C)C1(6367)',
    structure = SMILES('CCC(=O)C1(O)CC(C)C1'),
    E0 = (-385.825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.15869,0.0989075,-5.94141e-05,8.74462e-09,3.13339e-12,-46206,35.9606], Tmin=(100,'K'), Tmax=(1086.57,'K')), NASAPolynomial(coeffs=[19.9739,0.0448652,-1.76006e-05,3.20923e-09,-2.22144e-13,-52200.6,-74.1987], Tmin=(1086.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-385.825,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(561.227,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutane)"""),
)

species(
    label = '[CH2]C(C[CH]C)=C(CC)OO(6859)',
    structure = SMILES('[CH2]C(C[CH]C)=C(CC)OO'),
    E0 = (117.854,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1310,387.5,850,1000,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.624817,0.10698,-9.28611e-05,4.64272e-08,-1.00561e-11,14336.6,41.528], Tmin=(100,'K'), Tmax=(1065.24,'K')), NASAPolynomial(coeffs=[11.8631,0.0600866,-2.68285e-05,5.10091e-09,-3.5717e-13,11676.2,-19.508], Tmin=(1065.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.854,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(Allyl_P) + radical(RCCJC)"""),
)

species(
    label = 'C=C(O)[C](CC)OC[CH]C(6360)',
    structure = SMILES('[CH2]C(O)=C(CC)OC[CH]C'),
    E0 = (-165.49,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.00063,0.114336,-7.64967e-05,3.60551e-09,1.09278e-11,-19671.8,42.4312], Tmin=(100,'K'), Tmax=(959.873,'K')), NASAPolynomial(coeffs=[27.8003,0.0309665,-1.00012e-05,1.72384e-09,-1.20585e-13,-27273.2,-109.916], Tmin=(959.873,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-165.49,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(O)CJ) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C(O)(C[CH]C)C(=C)OC(6860)',
    structure = SMILES('[CH2]C(O)(C[CH]C)C(=C)OC'),
    E0 = (-73.6938,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,180,1031.07,1186.85,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.160458,'amu*angstrom^2'), symmetry=1, barrier=(3.68925,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160458,'amu*angstrom^2'), symmetry=1, barrier=(3.68925,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160458,'amu*angstrom^2'), symmetry=1, barrier=(3.68925,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160458,'amu*angstrom^2'), symmetry=1, barrier=(3.68925,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160458,'amu*angstrom^2'), symmetry=1, barrier=(3.68925,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160458,'amu*angstrom^2'), symmetry=1, barrier=(3.68925,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160458,'amu*angstrom^2'), symmetry=1, barrier=(3.68925,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160458,'amu*angstrom^2'), symmetry=1, barrier=(3.68925,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.5946,0.118896,-0.000114791,5.97165e-08,-1.25616e-11,-8658.34,43.747], Tmin=(100,'K'), Tmax=(1144.89,'K')), NASAPolynomial(coeffs=[19.618,0.044785,-1.76951e-05,3.17862e-09,-2.16102e-13,-13515.6,-61.4622], Tmin=(1144.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-73.6938,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)(O)CJ) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C(O)(C[CH]C)C(O)=CC(6861)',
    structure = SMILES('[CH2]C(O)(C[CH]C)C(O)=CC'),
    E0 = (-128.404,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,180,1600,1691.85,2813.61,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154354,'amu*angstrom^2'), symmetry=1, barrier=(3.54889,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154354,'amu*angstrom^2'), symmetry=1, barrier=(3.54889,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154354,'amu*angstrom^2'), symmetry=1, barrier=(3.54889,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154354,'amu*angstrom^2'), symmetry=1, barrier=(3.54889,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154354,'amu*angstrom^2'), symmetry=1, barrier=(3.54889,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154354,'amu*angstrom^2'), symmetry=1, barrier=(3.54889,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154354,'amu*angstrom^2'), symmetry=1, barrier=(3.54889,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154354,'amu*angstrom^2'), symmetry=1, barrier=(3.54889,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.714,0.122611,-0.000125239,6.91531e-08,-1.53672e-11,-15235,44.3476], Tmin=(100,'K'), Tmax=(1089.91,'K')), NASAPolynomial(coeffs=[19.6784,0.0440995,-1.71866e-05,3.06001e-09,-2.06896e-13,-19898.1,-60.7], Tmin=(1089.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-128.404,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=CC(C)(O)CJ) + radical(RCCJC)"""),
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
    label = 'C[CH]CC(O)=C([O])CC(6536)',
    structure = SMILES('C[CH]CC(O)=C([O])CC'),
    E0 = (-195.839,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.03946,0.0998398,-9.06771e-05,4.34978e-08,-8.3118e-12,-23363.3,39.4826], Tmin=(100,'K'), Tmax=(1269.03,'K')), NASAPolynomial(coeffs=[20.1845,0.0329417,-1.16034e-05,1.95759e-09,-1.28347e-13,-28750.1,-67.9673], Tmin=(1269.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-195.839,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(RCCJC) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CHCH3(T)(95)',
    structure = SMILES('[CH]C'),
    E0 = (343.893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,592.414,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00438699,'amu*angstrom^2'), symmetry=1, barrier=(26.7685,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.82363,-0.000909515,3.2138e-05,-3.7348e-08,1.3309e-11,41371.4,7.10948], Tmin=(100,'K'), Tmax=(960.812,'K')), NASAPolynomial(coeffs=[4.30487,0.00943069,-3.27559e-06,5.95121e-10,-4.27307e-14,40709.1,1.84202], Tmin=(960.812,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CHCH3(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C([CH2])(O)C(=O)CC(5988)',
    structure = SMILES('[CH2]C([CH2])(O)C(=O)CC'),
    E0 = (-64.1774,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,3615,1277.5,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.24882,0.119041,-0.00016428,1.17739e-07,-3.31244e-11,-7533.04,31.7431], Tmin=(100,'K'), Tmax=(875.483,'K')), NASAPolynomial(coeffs=[18.1268,0.0305211,-1.26247e-05,2.26293e-09,-1.51404e-13,-10925.9,-59.1577], Tmin=(875.483,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-64.1774,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJC(C)(C=O)O) + radical(CJC(C)(C=O)O)"""),
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
    E0 = (-128.903,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-34.3745,'kJ/mol'),
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
    E0 = (10.1891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (18.7355,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-57.9958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-10.6709,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-128.903,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-4.45674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (35.5448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (29.5923,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (40.1463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-10.9638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-27.4767,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-42.2942,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-99.866,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-74.7202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-76.1846,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-69.0834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (99.9221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (57.5222,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-70.7454,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (267.464,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-65.5029,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-103.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (316.376,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (28.4154,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (28.4154,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (35.9739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-120.619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (260.967,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-22.3761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (185.296,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (14.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (185.724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (279.715,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    products = ['C3H6(72)', 'C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    products = ['C[CH]CC1(O)CC1([O])CC(6840)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.34238e+09,'s^-1'), n=0.889391, Ea=(94.5285,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonylbond_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 90.6 to 94.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    products = ['[CH2]C1(O)CC(C)C1([O])CC(6781)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.13771e+06,'s^-1'), n=1.58803, Ea=(98.6375,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_csHNd] for rate rule [R5_SS_CO;carbonylbond_intra_Nd;radadd_intra_csHNd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 94.9 to 98.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[CH2]C(O)(C=CC)C(=O)CC(6841)'],
    products = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', '[CH2]C(O)(CC=C)C(=O)CC(6842)'],
    products = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.36e+08,'cm^3/(mol*s)'), n=1.56, Ea=(2.5104,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 9 used for Cds-HH_Cds-CsH;HJ
Exact match found for rate rule [Cds-HH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C3H6(T)(143)', 'C=C(O)C(=O)CC(4626)'],
    products = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.00392164,'m^3/(mol*s)'), n=2.41519, Ea=(15.6067,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;CJ] for rate rule [Cds-COOs_Cds;CJ]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C2H5CO(71)', 'C=C(O)C[CH]C(2383)'],
    products = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=0, Ea=(93.9308,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_R;CO_rad] for rate rule [Cds-OsCs_Cds;CO_rad/NonDe]
Euclidian distance = 2.82842712475
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C3H6(72)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.00620445,'m^3/(mol*s)'), n=2.46568, Ea=(50.0638,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-Cs\H3/H;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 49.0 to 50.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['OH(5)', 'C=C(C[CH]C)C(=O)CC(6843)'],
    products = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.52832,'m^3/(mol*s)'), n=2.02802, Ea=(14.4047,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDeCs_Cds;YJ] for rate rule [Cds-COCs_Cds;OJ_pri]
Euclidian distance = 2.2360679775
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(O)([CH]CC)C(=O)CC(6743)'],
    products = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6.76e+09,'s^-1'), n=0.88, Ea=(158.992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 357 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]CCC([CH2])(O)C(=O)CC(6844)'],
    products = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6.48e+07,'s^-1'), n=1.57, Ea=(147.695,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 106 used for R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C[CH]CC(C)([O])C(=O)CC(6845)'],
    products = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    products = ['C[CH][CH]C(C)(O)C(=O)CC(6846)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C([O])(CCC)C(=O)CC(6847)'],
    products = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(210000,'s^-1'), n=1.76, Ea=(53.8899,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 326 used for R4H_SSS;O_rad_out;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    products = ['C[CH]CC(C)(O)C(=O)[CH]C(6848)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(6.44e+09,'s^-1'), n=0.13, Ea=(86.6088,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 131 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    products = ['[CH2]C(O)(CCC)C(=O)[CH]C(6849)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.1016,'s^-1'), n=3.24, Ea=(29.037,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_H/NonDeC;Cs_H_out_1H] for rate rule [R5H_CC(O2d)CC;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    products = ['[CH2][CH]CC(C)(O)C(=O)CC(6850)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(91273.5,'s^-1'), n=1.79, Ea=(54.1828,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;Cs_H_out_2H] for rate rule [R5HJ_3;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    products = ['[CH2]CC(=O)C(C)(O)C[CH]C(6851)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(68850,'s^-1'), n=1.68, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_2H;Cs_H_out_2H] for rate rule [R5H_CCC(O2d)C;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]CC(=O)C([CH2])(O)CCC(6852)'],
    products = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(25800,'s^-1'), n=1.67, Ea=(42.6768,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 114 used for R6H_SSSSS;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R6H_SSSSS;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C3H6(T)(143)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    products = ['C[CH]CC1(O)CO[C]1CC(6853)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(7.01137e+09,'s^-1'), n=0.572544, Ea=(186.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4_S;multiplebond_intra;radadd_intra_cs2H] + [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    products = ['[CH2]C1(O)CC(C)O[C]1CC(6761)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.64e+06,'s^-1'), n=1.15, Ea=(58.1576,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_csHCs] for rate rule [R5_SS_CO;carbonyl_intra_Nd;radadd_intra_csHCs]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CH3CH2OH(54)', '[CH2]C(=C=O)C[CH]C(6854)'],
    products = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.79e-11,'m^3/(mol*s)'), n=3.97, Ea=(329.281,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [doublebond;Cs_OH] for rate rule [Cdd_Cd;C_pri_OH]
Euclidian distance = 1.41421356237
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    products = ['CC=CC(C)(O)C(=O)CC(6855)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    products = ['C=CCC(C)(O)C(=O)CC(6856)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CH2(S)(23)', '[CH2]C(O)(C[CH]C)C(C)=O(6857)'],
    products = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;C_pri/De]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    products = ['C[CH]CC[C](O)C(=O)CC(6361)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    products = ['C[CH]C[C](O)CC(=O)CC(6858)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ-HH;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;CO]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C(C)C([CH2])(O)C(=O)CC(6366)'],
    products = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    products = ['CCC(=O)C1(O)CC(C)C1(6367)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(C[CH]C)=C(CC)OO(6859)'],
    products = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O] for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O_R]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=C(O)[C](CC)OC[CH]C(6360)'],
    products = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O] for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O_C]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(O)(C[CH]C)C(=C)OC(6860)'],
    products = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(37989.5,'s^-1'), n=2.515, Ea=(258.99,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C] + [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O] for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_C]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C(O)(C[CH]C)C(O)=CC(6861)'],
    products = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O_H] for rate rule [R_ROR;R1_doublebond_CHCH3;R2_doublebond_CsC;R_O_H]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction35',
    reactants = ['CH2(19)', 'C[CH]CC(O)=C([O])CC(6536)'],
    products = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/ODMustO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction36',
    reactants = ['CHCH3(T)(95)', '[CH2]C([CH2])(O)C(=O)CC(5988)'],
    products = ['[CH2]C(O)(C[CH]C)C(=O)CC(6362)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2.13464e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '1399',
    isomers = [
        '[CH2]C(O)(C[CH]C)C(=O)CC(6362)',
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
    label = '1399',
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

