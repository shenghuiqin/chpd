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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.69604,0.112314,-8.62121e-05,2.53836e-08,4.95704e-13,-11276.4,43.8231], Tmin=(100,'K'), Tmax=(993.19,'K')), NASAPolynomial(coeffs=[22.7469,0.0391271,-1.38224e-05,2.396e-09,-1.62646e-13,-17377.3,-80.2039], Tmin=(993.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-95.5573,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Isobutyl) + radical(C=CC(C)2OJ)"""),
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
    label = '[CH2]C(C)C1(CC)OC1([CH2])O(6780)',
    structure = SMILES('[CH2]C(C)C1(CC)OC1([CH2])O'),
    E0 = (-52.7022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.11452,0.138799,-0.000139222,6.94714e-08,-1.28542e-11,-6012.33,45.4421], Tmin=(100,'K'), Tmax=(1575.24,'K')), NASAPolynomial(coeffs=[30.8422,0.021798,-9.2137e-07,-4.38826e-10,4.69072e-14,-13522.2,-127.969], Tmin=(1575.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-52.7022,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(CJC(O)2C) + radical(Isobutyl)"""),
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
    label = 'C=C(C)C([O])(CC)C(=C)O(6782)',
    structure = SMILES('C=C(C)C([O])(CC)C(=C)O'),
    E0 = (-188.644,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,325,375,415,465,420,450,1700,1750,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,180,180,421.612,715.777,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155443,'amu*angstrom^2'), symmetry=1, barrier=(3.57395,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155443,'amu*angstrom^2'), symmetry=1, barrier=(3.57395,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155443,'amu*angstrom^2'), symmetry=1, barrier=(3.57395,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155443,'amu*angstrom^2'), symmetry=1, barrier=(3.57395,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155443,'amu*angstrom^2'), symmetry=1, barrier=(3.57395,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155443,'amu*angstrom^2'), symmetry=1, barrier=(3.57395,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.188,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.56788,0.108855,-8.96929e-05,3.75133e-08,-6.26654e-12,-22476.9,40.0985], Tmin=(100,'K'), Tmax=(1430.57,'K')), NASAPolynomial(coeffs=[23.9089,0.0376187,-1.4999e-05,2.70457e-09,-1.83477e-13,-29766.1,-91.9345], Tmin=(1430.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-188.644,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ)"""),
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
    label = '[CH2]C(C)C(=O)C(=C)O(6783)',
    structure = SMILES('[CH2]C(C)C(=O)C(=C)O'),
    E0 = (-177.498,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,375,552.5,462.5,1710,348.913,353.5],'cm^-1')),
        HinderedRotor(inertia=(0.141425,'amu*angstrom^2'), symmetry=1, barrier=(12.1817,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141268,'amu*angstrom^2'), symmetry=1, barrier=(12.1886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142762,'amu*angstrom^2'), symmetry=1, barrier=(12.1706,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141086,'amu*angstrom^2'), symmetry=1, barrier=(12.171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142088,'amu*angstrom^2'), symmetry=1, barrier=(12.1827,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.103763,0.0925295,-0.000102855,6.09942e-08,-1.46132e-11,-21202.1,28.6053], Tmin=(100,'K'), Tmax=(1008.95,'K')), NASAPolynomial(coeffs=[14.835,0.0333044,-1.48052e-05,2.81492e-09,-1.97385e-13,-24216.6,-43.5988], Tmin=(1008.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-177.498,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsHH) + radical(CJC(C)C=O)"""),
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
    label = '[CH2]C(C)C(=O)CC(6784)',
    structure = SMILES('[CH2]C(C)C(=O)CC'),
    E0 = (-102.641,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,292.285,2269.5],'cm^-1')),
        HinderedRotor(inertia=(0.149715,'amu*angstrom^2'), symmetry=1, barrier=(9.08411,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149809,'amu*angstrom^2'), symmetry=1, barrier=(9.08268,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149801,'amu*angstrom^2'), symmetry=1, barrier=(9.08044,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149863,'amu*angstrom^2'), symmetry=1, barrier=(9.08348,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149789,'amu*angstrom^2'), symmetry=1, barrier=(9.08409,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1509,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20579,0.0721848,-3.85334e-05,-4.85428e-08,6.40128e-11,-12255.1,24.9679], Tmin=(100,'K'), Tmax=(496.027,'K')), NASAPolynomial(coeffs=[6.44739,0.0483312,-2.20879e-05,4.19976e-09,-2.9226e-13,-13001.6,1.0716], Tmin=(496.027,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-102.641,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJC(C)C=O)"""),
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
    label = 'C=C(O)C([O])(CC)[C](C)C(6785)',
    structure = SMILES('C=C(O)C([O])(CC)[C](C)C'),
    E0 = (-148.066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.45794,0.104144,-6.28528e-05,3.67344e-09,7.04785e-12,-17597.8,40.6902], Tmin=(100,'K'), Tmax=(1016.18,'K')), NASAPolynomial(coeffs=[22.6787,0.040483,-1.51548e-05,2.73905e-09,-1.90927e-13,-24121.8,-84.1066], Tmin=(1016.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-148.066,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCJ(C)CO)"""),
)

species(
    label = '[CH2][C](C)C(O)(CC)C(=C)O(6786)',
    structure = SMILES('[CH2][C](C)C(O)(CC)C(=C)O'),
    E0 = (-172.086,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.73754,0.111154,-7.52775e-05,6.57773e-09,9.39472e-12,-20476.9,42.9538], Tmin=(100,'K'), Tmax=(944.185,'K')), NASAPolynomial(coeffs=[24.9751,0.0342894,-1.08359e-05,1.79767e-09,-1.2156e-13,-27139.4,-92.9538], Tmin=(944.185,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-172.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Isobutyl) + radical(CCJ(C)CO)"""),
)

species(
    label = '[CH2]C(C)C(O)([CH]C)C(=C)O(6787)',
    structure = SMILES('[CH2]C(C)C(O)([CH]C)C(=C)O'),
    E0 = (-124.758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.80199,0.124005,-0.000116902,5.59015e-08,-1.02595e-11,-14738.2,49.0896], Tmin=(100,'K'), Tmax=(1476.78,'K')), NASAPolynomial(coeffs=[28.6857,0.0271586,-6.79237e-06,8.94335e-10,-5.0253e-14,-22777.8,-110.829], Tmin=(1476.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-124.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCJCO) + radical(Isobutyl)"""),
)

species(
    label = 'C=C(O)C([O])([CH]C)C(C)C(6788)',
    structure = SMILES('C=C(O)C([O])([CH]C)C(C)C'),
    E0 = (-100.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.52394,0.105461,-6.52292e-05,3.5826e-09,7.85926e-12,-11902.9,43.2276], Tmin=(100,'K'), Tmax=(998.632,'K')), NASAPolynomial(coeffs=[23.5338,0.0381847,-1.38818e-05,2.48657e-09,-1.73344e-13,-18557.7,-85.8888], Tmin=(998.632,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-100.738,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C([CH2])C(O)(CC)C(=C)O(6789)',
    structure = SMILES('[CH2]C([CH2])C(O)(CC)C(=C)O'),
    E0 = (-119.578,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.69514,0.127776,-0.000128026,6.61265e-08,-1.31714e-11,-14124.2,47.2834], Tmin=(100,'K'), Tmax=(1340.98,'K')), NASAPolynomial(coeffs=[27.1686,0.0293327,-7.435e-06,9.67713e-10,-5.30147e-14,-21291.6,-102.415], Tmin=(1340.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-119.578,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]CC(O)(C(=C)O)C([CH2])C(6790)',
    structure = SMILES('[CH2]CC(O)(C(=C)O)C([CH2])C'),
    E0 = (-119.414,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.62234,0.126114,-0.000123118,6.16467e-08,-1.19698e-11,-14107.1,47.7255], Tmin=(100,'K'), Tmax=(1343.44,'K')), NASAPolynomial(coeffs=[27.622,0.029731,-8.43277e-06,1.2271e-09,-7.35046e-14,-21662,-104.988], Tmin=(1343.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-119.414,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C(C)C(O)(CC)C(=C)[O](6791)',
    structure = SMILES('[CH2]C(C)C(O)(CC)C(=C)[O]'),
    E0 = (-186.855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.6813,0.114295,-9.58338e-05,3.56964e-08,-2.78241e-12,-22259.4,44.0568], Tmin=(100,'K'), Tmax=(970.256,'K')), NASAPolynomial(coeffs=[21.8701,0.0390695,-1.33453e-05,2.24813e-09,-1.49498e-13,-27859,-74.1586], Tmin=(970.256,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-186.855,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=C(O)C(O)(CC)C([CH2])C(6792)',
    structure = SMILES('[CH]=C(O)C(O)(CC)C([CH2])C'),
    E0 = (-77.5638,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.69669,0.128502,-0.000128063,6.53317e-08,-1.29108e-11,-9071.68,46.6069], Tmin=(100,'K'), Tmax=(1316.98,'K')), NASAPolynomial(coeffs=[27.9228,0.0294253,-8.2952e-06,1.20011e-09,-7.15985e-14,-16609.7,-107.544], Tmin=(1316.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-77.5638,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH2]CC([O])(C(=C)O)C(C)C(6793)',
    structure = SMILES('[CH2]CC([O])(C(=C)O)C(C)C'),
    E0 = (-95.3932,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,429.526,704.561,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155597,'amu*angstrom^2'), symmetry=1, barrier=(3.57748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155597,'amu*angstrom^2'), symmetry=1, barrier=(3.57748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155597,'amu*angstrom^2'), symmetry=1, barrier=(3.57748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155597,'amu*angstrom^2'), symmetry=1, barrier=(3.57748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155597,'amu*angstrom^2'), symmetry=1, barrier=(3.57748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155597,'amu*angstrom^2'), symmetry=1, barrier=(3.57748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155597,'amu*angstrom^2'), symmetry=1, barrier=(3.57748,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.72969,0.111952,-8.60501e-05,2.73055e-08,-1.12341e-12,-11254.7,43.2575], Tmin=(100,'K'), Tmax=(1046.92,'K')), NASAPolynomial(coeffs=[23.1775,0.0394985,-1.47792e-05,2.64153e-09,-1.81735e-13,-17714.5,-83.9914], Tmin=(1046.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-95.3932,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJ) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C([O])C([O])(CC)C(C)C(6794)',
    structure = SMILES('C=C([O])C([O])(CC)C(C)C'),
    E0 = (-162.835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.55703,0.108925,-8.81685e-05,3.7439e-08,-6.41051e-12,-19373.4,42.3623], Tmin=(100,'K'), Tmax=(1393.29,'K')), NASAPolynomial(coeffs=[21.9192,0.041527,-1.56095e-05,2.72088e-09,-1.81017e-13,-25915.3,-78.6831], Tmin=(1393.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-162.835,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(557.07,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C(O)C([O])(CC)C(C)C(6795)',
    structure = SMILES('[CH]=C(O)C([O])(CC)C(C)C'),
    E0 = (-53.5435,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,3615,1277.5,1000,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,180,1600,1633.01,2861.9,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152016,'amu*angstrom^2'), symmetry=1, barrier=(3.49515,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152016,'amu*angstrom^2'), symmetry=1, barrier=(3.49515,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152016,'amu*angstrom^2'), symmetry=1, barrier=(3.49515,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152016,'amu*angstrom^2'), symmetry=1, barrier=(3.49515,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152016,'amu*angstrom^2'), symmetry=1, barrier=(3.49515,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152016,'amu*angstrom^2'), symmetry=1, barrier=(3.49515,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152016,'amu*angstrom^2'), symmetry=1, barrier=(3.49515,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.85119,0.114867,-9.26999e-05,3.30164e-08,-2.85201e-12,-6217.21,42.3101], Tmin=(100,'K'), Tmax=(1050.73,'K')), NASAPolynomial(coeffs=[23.6812,0.0388637,-1.44586e-05,2.57244e-09,-1.7641e-13,-12752.8,-87.7004], Tmin=(1050.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-53.5435,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C(C)C1(CC)OC[C]1O(6796)',
    structure = SMILES('[CH2]C(C)C1(CC)OC[C]1O'),
    E0 = (-45.649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.95346,0.116196,-0.000108372,5.60439e-08,-1.14541e-11,-5263.4,40.7627], Tmin=(100,'K'), Tmax=(1298.07,'K')), NASAPolynomial(coeffs=[20.4885,0.0394633,-1.09462e-05,1.51044e-09,-8.51504e-14,-10451.2,-70.9022], Tmin=(1298.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-45.649,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(557.07,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(C2CsJOH) + radical(Isobutyl)"""),
)

species(
    label = 'CCC1([O])[C](O)CCC1C(6797)',
    structure = SMILES('CCC1([O])[C](O)CCC1C'),
    E0 = (-129.437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.919339,0.0915417,-3.46749e-05,-1.95757e-08,1.39152e-11,-15376,36.9953], Tmin=(100,'K'), Tmax=(1010.74,'K')), NASAPolynomial(coeffs=[19.9882,0.0440714,-1.65705e-05,3.00813e-09,-2.10306e-13,-21404,-73.0075], Tmin=(1010.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-129.437,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(565.384,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(CC(C)2OJ) + radical(C2CsJOH)"""),
)

species(
    label = 'C=C(C)C(O)(CC)C(=C)O(6798)',
    structure = SMILES('C=C(C)C(O)(CC)C(=C)O'),
    E0 = (-417.747,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.87259,0.114456,-9.15444e-05,3.18863e-08,-2.66008e-12,-50019.3,39.9653], Tmin=(100,'K'), Tmax=(1072.58,'K')), NASAPolynomial(coeffs=[24.397,0.0379216,-1.44859e-05,2.62134e-09,-1.81487e-13,-56887.4,-94.3584], Tmin=(1072.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-417.747,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
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
    label = '[CH2]CC([O])(CC)C(=C)O(6315)',
    structure = SMILES('[CH2]CC([O])(CC)C(=C)O'),
    E0 = (-65.9187,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,180,180,180,472.217,665.457,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156676,'amu*angstrom^2'), symmetry=1, barrier=(3.6023,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156676,'amu*angstrom^2'), symmetry=1, barrier=(3.6023,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156676,'amu*angstrom^2'), symmetry=1, barrier=(3.6023,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156676,'amu*angstrom^2'), symmetry=1, barrier=(3.6023,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156676,'amu*angstrom^2'), symmetry=1, barrier=(3.6023,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156676,'amu*angstrom^2'), symmetry=1, barrier=(3.6023,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (128.169,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4527.98,'J/mol'), sigma=(7.66109,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=707.26 K, Pc=22.85 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.0521,0.0991652,-8.03971e-05,2.90918e-08,-2.72955e-12,-7736.02,39.5541], Tmin=(100,'K'), Tmax=(1056.62,'K')), NASAPolynomial(coeffs=[21.0046,0.0335009,-1.24976e-05,2.22591e-09,-1.52659e-13,-13392.7,-72.7829], Tmin=(1056.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-65.9187,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C(C)C(C)([O])C(=C)O(6799)',
    structure = SMILES('[CH2]C(C)C(C)([O])C(=C)O'),
    E0 = (-71.777,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.933302,0.0961403,-6.74066e-05,1.07386e-08,5.51236e-12,-8444.03,38.8414], Tmin=(100,'K'), Tmax=(959.329,'K')), NASAPolynomial(coeffs=[21.282,0.031643,-1.05449e-05,1.79074e-09,-1.21566e-13,-14000.8,-74.1588], Tmin=(959.329,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-71.777,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Isobutyl) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C(O)C([O])(CC)C[CH]C(6359)',
    structure = SMILES('C=C(O)C([O])(CC)C[CH]C'),
    E0 = (-100.499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.37847,0.111562,-9.67436e-05,4.49979e-08,-8.52632e-12,-11888.1,43.919], Tmin=(100,'K'), Tmax=(1257.13,'K')), NASAPolynomial(coeffs=[19.371,0.0455405,-1.79667e-05,3.2218e-09,-2.18469e-13,-17105,-60.9332], Tmin=(1257.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-100.499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJC) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C(O)C([O])([CH]CC)CC(6800)',
    structure = SMILES('C=C(O)C([O])([CH]CC)CC'),
    E0 = (-95.0432,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.60905,0.108824,-7.81789e-05,1.95543e-08,1.52878e-12,-11216.5,44.5071], Tmin=(100,'K'), Tmax=(1033.85,'K')), NASAPolynomial(coeffs=[22.9346,0.0394958,-1.47805e-05,2.65393e-09,-1.83543e-13,-17661.2,-81.3439], Tmin=(1033.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-95.0432,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCJCO) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C(O)C1(CC)OCC1C(6370)',
    structure = SMILES('C=C(O)C1(CC)OCC1C'),
    E0 = (-349.058,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.00796,0.0865248,5.95177e-06,-8.72994e-08,4.65699e-11,-41779.4,35.0797], Tmin=(100,'K'), Tmax=(898.721,'K')), NASAPolynomial(coeffs=[26.4665,0.0286383,-4.91288e-06,4.87515e-10,-2.80744e-14,-49318.4,-109.003], Tmin=(898.721,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-349.058,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(561.227,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Oxetane)"""),
)

species(
    label = '[CH2]C(C)C([C]=C)(CC)OO(6801)',
    structure = SMILES('[CH2]C(C)C([C]=C)(CC)OO'),
    E0 = (199.245,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.51195,0.118132,-0.000111036,5.69274e-08,-1.19162e-11,24164.8,43.1301], Tmin=(100,'K'), Tmax=(1143.76,'K')), NASAPolynomial(coeffs=[18.4729,0.0482401,-1.93764e-05,3.50171e-09,-2.38629e-13,19593.2,-55.9697], Tmin=(1143.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(199.245,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(C)C([O])(CC)C(C)=O(6802)',
    structure = SMILES('[CH2]C(C)C([O])(CC)C(C)=O'),
    E0 = (-91.9533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.982168,0.113662,-0.000115404,6.95678e-08,-1.76919e-11,-10883.5,41.5742], Tmin=(100,'K'), Tmax=(937.124,'K')), NASAPolynomial(coeffs=[12.4309,0.0564104,-2.37646e-05,4.37576e-09,-3.0041e-13,-13397.5,-22.2651], Tmin=(937.124,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-91.9533,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(Isobutyl) + radical(CC(C)(C=O)OJ)"""),
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
    label = '[CH2]C(C)[C](CC)C(=C)O(6803)',
    structure = SMILES('[CH2]C(O)=C(CC)C([CH2])C'),
    E0 = (7.62438,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800],'cm^-1')),
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
    molecularWeight = (126.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.21982,0.101609,-6.6625e-05,5.59552e-09,8.1907e-12,1116.92,38.8158], Tmin=(100,'K'), Tmax=(946.492,'K')), NASAPolynomial(coeffs=[21.9351,0.0352707,-1.1442e-05,1.90923e-09,-1.28402e-13,-4678.04,-79.078], Tmin=(946.492,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(7.62438,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsCsOs) + radical(Isobutyl) + radical(C=C(O)CJ)"""),
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
    label = 'C=C(O)C([O])([CH]C)CC(6518)',
    structure = SMILES('C=C(O)C([O])([CH]C)CC'),
    E0 = (-71.263,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,524.731,610.588,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157874,'amu*angstrom^2'), symmetry=1, barrier=(3.62983,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157874,'amu*angstrom^2'), symmetry=1, barrier=(3.62983,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157874,'amu*angstrom^2'), symmetry=1, barrier=(3.62983,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157874,'amu*angstrom^2'), symmetry=1, barrier=(3.62983,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157874,'amu*angstrom^2'), symmetry=1, barrier=(3.62983,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157874,'amu*angstrom^2'), symmetry=1, barrier=(3.62983,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (128.169,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.840076,0.0926017,-5.93304e-05,5.06435e-09,6.37595e-12,-8384.46,39.5015], Tmin=(100,'K'), Tmax=(995.569,'K')), NASAPolynomial(coeffs=[21.324,0.0322481,-1.16347e-05,2.07897e-09,-1.44926e-13,-14219.8,-74.4714], Tmin=(995.569,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-71.263,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCJCO) + radical(C=CC(C)2OJ)"""),
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
    E0 = (-95.5573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-51.2214,'kJ/mol'),
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
    E0 = (23.1482,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-41.3857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-41.5909,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (48.7863,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-95.5573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (16.3468,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (7.12227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-20.2766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-19.1959,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-8.94847,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-35.8976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-35.7336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (6.25033,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-33.2552,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-42.6748,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-52.8805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-20.5033,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (99.9221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (29.3351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (19.5027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-32.1571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (353.943,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (348.085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (64.3779,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (99.2079,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-87.273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (306.356,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (108.622,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (250.629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (310.3,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    products = ['C3H6(72)', 'C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    products = ['[CH2]C(C)C1(CC)OC1([CH2])O(6780)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.52e+08,'s^-1'), n=0.89, Ea=(44.3359,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H_secNd;radadd_intra] for rate rule [R4_S_D;doublebond_intra_2H_secNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    products = ['[CH2]C1(O)CC(C)C1([O])CC(6781)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.48e+06,'s^-1'), n=1.25, Ea=(65.2918,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 349 used for R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 65.0 to 65.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C=C(C)C([O])(CC)C(=C)O(6782)'],
    products = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
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
    products = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.97e+07,'cm^3/(mol*s)'), n=1.88, Ea=(32.2168,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO-CdCs_O;YJ] for rate rule [CO-CdCs_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C2H5(29)', '[CH2]C(C)C(=O)C(=C)O(6783)'],
    products = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7.94e+10,'cm^3/(mol*s)'), n=0, Ea=(28.0328,'kJ/mol'), T0=(1,'K'), Tmin=(333,'K'), Tmax=(363,'K'), comment="""Estimated using template [CO_O;CsJ-CsHH] for rate rule [CO-CdCs_O;CsJ-CsHH]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2COH(99)', '[CH2]C(C)C(=O)CC(6784)'],
    products = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.16e+10,'cm^3/(mol*s)'), n=0, Ea=(48.1578,'kJ/mol'), T0=(1,'K'), Tmin=(413,'K'), Tmax=(563,'K'), comment="""Estimated using template [CO-CsCs_O;CJ] for rate rule [CO-CsCs_O;CdsJ-O2s]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C3H6(72)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.00168615,'m^3/(mol*s)'), n=2.52599, Ea=(83.4095,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 78.8 to 83.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['CH3(17)', 'C=CC([O])(CC)C(=C)O(6516)'],
    products = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(10000,'cm^3/(mol*s)'), n=2.41, Ea=(29.7482,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 417 used for Cds-CsH_Cds-HH;CsJ-HHH
Exact match found for rate rule [Cds-CsH_Cds-HH;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    products = ['C=C(O)C([O])(CC)[C](C)C(6785)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.265e-07,'s^-1'), n=5.639, Ea=(102.68,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 38 used for R2H_S;C_rad_out_2H;Cs_H_out_Cs2
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    products = ['[CH2][C](C)C(O)(CC)C(=C)O(6786)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(C)C(O)([CH]C)C(=C)O(6787)'],
    products = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(5.71,'s^-1'), n=3.021, Ea=(105.562,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 319 used for R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    products = ['C=C(O)C([O])([CH]C)C(C)C(6788)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(6.44e+09,'s^-1'), n=0.13, Ea=(86.6088,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 131 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C([CH2])C(O)(CC)C(=C)O(6789)'],
    products = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.72e-08,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]CC(O)(C(=C)O)C([CH2])C(6790)'],
    products = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    products = ['[CH2]C(C)C(O)(CC)C(=C)[O](6791)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1070.11,'s^-1'), n=2.50856, Ea=(101.808,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;O_H_out] + [R4H_SSS;O_rad_out;XH_out] + [R4H_SS(Cd)S;Y_rad_out;XH_out] for rate rule [R4H_SS(Cd)S;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C(O)C(O)(CC)C([CH2])C(6792)'],
    products = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]CC([O])(C(=C)O)C(C)C(6793)'],
    products = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(137700,'s^-1'), n=1.68, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 111 used for R5H_CCC;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R5H_CCC;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    products = ['C=C([O])C([O])(CC)C(C)C(6794)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2300,'s^-1'), n=1.98, Ea=(42.6768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC(Cd);C_rad_out_2H;XH_out] for rate rule [R5H_CCC(Cd);C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C(O)C([O])(CC)C(C)C(6795)'],
    products = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(816000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5H;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C3H6(T)(143)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    products = ['[CH2]C(C)C1(CC)OC[C]1O(6796)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.21748e+08,'s^-1'), n=0.95, Ea=(124.892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R4_S_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    products = ['CCC1([O])[C](O)CCC1C(6797)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.71e+11,'s^-1'), n=0.2, Ea=(115.06,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 4 used for R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    products = ['C=C(C)C(O)(CC)C(=C)O(6798)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CH2(S)(23)', '[CH2]CC([O])(CC)C(=C)O(6315)'],
    products = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['CH2(S)(23)', '[CH2]C(C)C(C)([O])C(=C)O(6799)'],
    products = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    products = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    products = ['C=C(O)C([O])([CH]CC)CC(6800)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(5.59192e+09,'s^-1'), n=1.025, Ea=(194.765,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CH3] for rate rule [cCs(-HC)CJ;CsJ-HH;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    products = ['C=C(O)C1(CC)OCC1C(6370)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R4_SSS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(C)C([C]=C)(CC)OO(6801)'],
    products = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(3.01978e+11,'s^-1'), n=0, Ea=(107.111,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2OOH_S;Y_rad_out] for rate rule [R2OOH_S;Cd_rad_out_double]
Euclidian distance = 2.0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    products = ['[CH2]C(C)C([O])(CC)C(C)=O(6802)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(204.179,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction32',
    reactants = ['O(4)', '[CH2]C(C)[C](CC)C(=C)O(6803)'],
    products = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/Cs2;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction33',
    reactants = ['CH2(19)', 'C=C(O)C([O])([CH]C)CC(6518)'],
    products = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '1400',
    isomers = [
        '[CH2]C(C)C([O])(CC)C(=C)O(6363)',
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
    label = '1400',
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

