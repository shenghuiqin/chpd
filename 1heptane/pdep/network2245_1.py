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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.97881,0.121303,-0.00012457,6.52765e-08,-1.3414e-11,-27684.8,44.0852], Tmin=(100,'K'), Tmax=(1191.33,'K')), NASAPolynomial(coeffs=[24.7714,0.0314856,-1.14799e-05,1.99084e-09,-1.33386e-13,-34058.4,-89.6517], Tmin=(1191.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-232.046,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CJCO)"""),
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
    label = '[CH2]C(O)C1(CC)OC1([CH2])O(12907)',
    structure = SMILES('[CH2]C(O)C1(CC)OC1([CH2])O'),
    E0 = (-185.844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.8231,0.141322,-0.000156321,8.36223e-08,-1.65409e-11,-22043.8,43.6229], Tmin=(100,'K'), Tmax=(1451.62,'K')), NASAPolynomial(coeffs=[33.0292,0.014361,1.13265e-06,-7.49571e-10,6.684e-14,-30065.3,-138.68], Tmin=(1451.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-185.844,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(CJCO) + radical(CJC(O)2C)"""),
)

species(
    label = '[CH2]C1(O)CC(O)C1([O])CC(12846)',
    structure = SMILES('[CH2]C1(O)CC(O)C1([O])CC'),
    E0 = (-169.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.93739,0.116512,-0.000110333,5.31472e-08,-1.00627e-11,-20210.3,38.3451], Tmin=(100,'K'), Tmax=(1287.94,'K')), NASAPolynomial(coeffs=[25.3857,0.0316549,-1.15049e-05,1.9919e-09,-1.33127e-13,-27248.5,-100.387], Tmin=(1287.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-169.914,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CJC(C)2O) + radical(CC(C)2OJ)"""),
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
    label = 'C=C(O)C([O])(CC)C(=C)O(12908)',
    structure = SMILES('C=C(O)C([O])(CC)C(=C)O'),
    E0 = (-364.147,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,325,375,415,465,420,450,1700,1750,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,454.217,676.019,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156007,'amu*angstrom^2'), symmetry=1, barrier=(3.58692,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156007,'amu*angstrom^2'), symmetry=1, barrier=(3.58692,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156007,'amu*angstrom^2'), symmetry=1, barrier=(3.58692,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156007,'amu*angstrom^2'), symmetry=1, barrier=(3.58692,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156007,'amu*angstrom^2'), symmetry=1, barrier=(3.58692,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156007,'amu*angstrom^2'), symmetry=1, barrier=(3.58692,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.16,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.57548,0.10608,-7.78101e-05,1.23252e-08,6.16488e-12,-43581.4,38.7175], Tmin=(100,'K'), Tmax=(990.041,'K')), NASAPolynomial(coeffs=[26.8606,0.0260142,-9.26308e-06,1.69453e-09,-1.21794e-13,-50918.6,-106.805], Tmin=(990.041,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-364.147,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ)"""),
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
    label = '[CH2]C(O)C(=O)C(=C)O(12909)',
    structure = SMILES('[CH2]C(O)C(=O)C(=C)O'),
    E0 = (-316.76,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,375,552.5,462.5,1710,380.512,380.516],'cm^-1')),
        HinderedRotor(inertia=(0.00116428,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146833,'amu*angstrom^2'), symmetry=1, barrier=(15.0864,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146833,'amu*angstrom^2'), symmetry=1, barrier=(15.0865,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146829,'amu*angstrom^2'), symmetry=1, barrier=(15.0864,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146829,'amu*angstrom^2'), symmetry=1, barrier=(15.0864,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.369417,0.0848859,-8.7961e-05,4.42745e-08,-8.58833e-12,-37930.6,32.5428], Tmin=(100,'K'), Tmax=(1270.88,'K')), NASAPolynomial(coeffs=[21.8672,0.0148981,-5.35601e-06,9.42496e-10,-6.43626e-14,-43582.7,-80.0664], Tmin=(1270.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-316.76,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsHH) + radical(CJCO)"""),
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
    label = '[CH2]C(O)C(=O)CC(5567)',
    structure = SMILES('[CH2]C(O)C(=O)CC'),
    E0 = (-241.903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,294.511,294.529],'cm^-1')),
        HinderedRotor(inertia=(0.00194327,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00194344,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00194348,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152949,'amu*angstrom^2'), symmetry=1, barrier=(9.41625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152967,'amu*angstrom^2'), symmetry=1, barrier=(9.41629,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.124,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.836082,0.0687768,-5.84164e-05,2.66337e-08,-5.03381e-12,-28979.8,29.1175], Tmin=(100,'K'), Tmax=(1239.65,'K')), NASAPolynomial(coeffs=[12.3054,0.0317684,-1.36353e-05,2.55091e-09,-1.77024e-13,-31823.3,-28.6792], Tmin=(1239.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-241.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJCO)"""),
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
    label = 'C=C(O)C([O])(CC)[C](C)O(7269)',
    structure = SMILES('C=C(O)C([O])(CC)[C](C)O'),
    E0 = (-267.008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.67948,0.121281,-0.00012937,7.21554e-08,-1.59866e-11,-31905.9,41.8036], Tmin=(100,'K'), Tmax=(1099.12,'K')), NASAPolynomial(coeffs=[21.3859,0.0373407,-1.48146e-05,2.67254e-09,-1.82455e-13,-36976.3,-71.6536], Tmin=(1099.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-267.008,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C(O)C(O)([CH]C)C(=C)O(12910)',
    structure = SMILES('[CH2]C(O)C(O)([CH]C)C(=C)O'),
    E0 = (-261.247,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3615,3650,1210,1277.5,1345,900,1000,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,180,180,180,180,1600,1659.78,2832.8,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153199,'amu*angstrom^2'), symmetry=1, barrier=(3.52234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153199,'amu*angstrom^2'), symmetry=1, barrier=(3.52234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153199,'amu*angstrom^2'), symmetry=1, barrier=(3.52234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153199,'amu*angstrom^2'), symmetry=1, barrier=(3.52234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153199,'amu*angstrom^2'), symmetry=1, barrier=(3.52234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153199,'amu*angstrom^2'), symmetry=1, barrier=(3.52234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153199,'amu*angstrom^2'), symmetry=1, barrier=(3.52234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153199,'amu*angstrom^2'), symmetry=1, barrier=(3.52234,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.56749,0.127126,-0.000135796,7.2019e-08,-1.46473e-11,-31169.5,47.4801], Tmin=(100,'K'), Tmax=(1276.75,'K')), NASAPolynomial(coeffs=[29.2182,0.0220574,-5.91059e-06,8.32846e-10,-4.93825e-14,-38838.9,-111.883], Tmin=(1276.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-261.247,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCJCO) + radical(CJCO)"""),
)

species(
    label = '[CH2][C](O)C(O)(CC)C(=C)O(12911)',
    structure = SMILES('[CH2][C](O)C(O)(CC)C(=C)O'),
    E0 = (-284.521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.38539,0.134848,-0.000159873,9.60084e-08,-2.24e-11,-33985.2,43.7535], Tmin=(100,'K'), Tmax=(1057.56,'K')), NASAPolynomial(coeffs=[25.4907,0.0294102,-1.03227e-05,1.73256e-09,-1.13397e-13,-39881.2,-92.2919], Tmin=(1057.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-284.521,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C2CsJOH) + radical(CJCO)"""),
)

species(
    label = 'C=C(O)C([O])(CC)C(C)[O](6325)',
    structure = SMILES('C=C(O)C([O])(CC)C(C)[O]'),
    E0 = (-213.275,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,180,1017.44,1196.65,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.161442,'amu*angstrom^2'), symmetry=1, barrier=(3.71187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161442,'amu*angstrom^2'), symmetry=1, barrier=(3.71187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161442,'amu*angstrom^2'), symmetry=1, barrier=(3.71187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161442,'amu*angstrom^2'), symmetry=1, barrier=(3.71187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161442,'amu*angstrom^2'), symmetry=1, barrier=(3.71187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161442,'amu*angstrom^2'), symmetry=1, barrier=(3.71187,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(5105.78,'J/mol'), sigma=(8.38548,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=797.51 K, Pc=19.65 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.88961,0.114667,-0.000107262,5.109e-08,-9.56476e-12,-25426.5,42.9466], Tmin=(100,'K'), Tmax=(1302.58,'K')), NASAPolynomial(coeffs=[25.2787,0.0312377,-1.11881e-05,1.91867e-09,-1.27447e-13,-32504.2,-95.3067], Tmin=(1302.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-213.275,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=C(O)C([O])([CH]C)C(C)O(7273)',
    structure = SMILES('C=C(O)C([O])([CH]C)C(C)O'),
    E0 = (-243.733,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.60492,0.110608,-9.53424e-05,3.58415e-08,-3.24401e-12,-29101.5,44.6047], Tmin=(100,'K'), Tmax=(1016.18,'K')), NASAPolynomial(coeffs=[23.947,0.031915,-1.14903e-05,2.02574e-09,-1.39163e-13,-35424.6,-84.6392], Tmin=(1016.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-243.733,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCJCO) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH2]CC(O)(C(=C)O)C([CH2])O(12912)',
    structure = SMILES('[CH2]CC(O)(C(=C)O)C([CH2])O'),
    E0 = (-255.903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3615,3650,1210,1277.5,1345,900,1000,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,180,180,1600,1639.49,2853,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152636,'amu*angstrom^2'), symmetry=1, barrier=(3.5094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152636,'amu*angstrom^2'), symmetry=1, barrier=(3.5094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152636,'amu*angstrom^2'), symmetry=1, barrier=(3.5094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152636,'amu*angstrom^2'), symmetry=1, barrier=(3.5094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152636,'amu*angstrom^2'), symmetry=1, barrier=(3.5094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152636,'amu*angstrom^2'), symmetry=1, barrier=(3.5094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152636,'amu*angstrom^2'), symmetry=1, barrier=(3.5094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152636,'amu*angstrom^2'), symmetry=1, barrier=(3.5094,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.23995,0.1274,-0.000135266,6.87119e-08,-1.24211e-11,-30544.5,45.5928], Tmin=(100,'K'), Tmax=(974.965,'K')), NASAPolynomial(coeffs=[26.6704,0.0270208,-8.88035e-06,1.4717e-09,-9.7555e-14,-37048.4,-97.5939], Tmin=(974.965,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-255.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C([O])C(O)(CC)C(=C)O(7275)',
    structure = SMILES('[CH2]C([O])C(O)(CC)C(=C)O'),
    E0 = (-230.788,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,502.136,635.28,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157202,'amu*angstrom^2'), symmetry=1, barrier=(3.61437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157202,'amu*angstrom^2'), symmetry=1, barrier=(3.61437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157202,'amu*angstrom^2'), symmetry=1, barrier=(3.61437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157202,'amu*angstrom^2'), symmetry=1, barrier=(3.61437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157202,'amu*angstrom^2'), symmetry=1, barrier=(3.61437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157202,'amu*angstrom^2'), symmetry=1, barrier=(3.61437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157202,'amu*angstrom^2'), symmetry=1, barrier=(3.61437,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.50119,0.127179,-0.000134364,7.09601e-08,-1.44666e-11,-27509.9,44.5545], Tmin=(100,'K'), Tmax=(1233.53,'K')), NASAPolynomial(coeffs=[28.5523,0.0246375,-7.43039e-06,1.14654e-09,-7.19708e-14,-35030.7,-111.21], Tmin=(1233.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-230.788,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CJCO) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C(O)C(O)(CC)C(=C)[O](12913)',
    structure = SMILES('[CH2]C(O)C(O)(CC)C(=C)[O]'),
    E0 = (-323.344,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.93054,0.122917,-0.000133063,7.43858e-08,-1.63053e-11,-38669.3,44.1971], Tmin=(100,'K'), Tmax=(1119.71,'K')), NASAPolynomial(coeffs=[23.5001,0.0320695,-1.1361e-05,1.92558e-09,-1.26967e-13,-44364.3,-81.3665], Tmin=(1119.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-323.344,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CJCO) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C(O)C(O)(CC)C([CH2])O(12914)',
    structure = SMILES('[CH]=C(O)C(O)(CC)C([CH2])O'),
    E0 = (-214.053,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3615,3650,1210,1277.5,1345,900,1000,1100,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.34821,0.130156,-0.000141353,7.36849e-08,-1.38324e-11,-25507.6,44.5981], Tmin=(100,'K'), Tmax=(969.962,'K')), NASAPolynomial(coeffs=[27.1189,0.0264794,-8.61332e-06,1.41525e-09,-9.32749e-14,-32063.3,-100.991], Tmin=(969.962,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-214.053,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cds_P) + radical(CJCO)"""),
)

species(
    label = '[CH2]CC([O])(C(=C)O)C(C)O(7278)',
    structure = SMILES('[CH2]CC([O])(C(=C)O)C(C)O'),
    E0 = (-238.389,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.80077,0.116962,-0.00011556,5.85622e-08,-1.16894e-11,-28453.8,44.6002], Tmin=(100,'K'), Tmax=(1222.27,'K')), NASAPolynomial(coeffs=[23.9726,0.0326151,-1.20478e-05,2.10274e-09,-1.41229e-13,-34754.1,-84.9144], Tmin=(1222.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-238.389,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(RCCJ)"""),
)

species(
    label = 'C=C([O])C([O])(CC)C(C)O(7279)',
    structure = SMILES('C=C([O])C([O])(CC)C(C)O'),
    E0 = (-305.831,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.25016,0.109663,-0.000103675,5.19621e-08,-1.04768e-11,-36589,42.3375], Tmin=(100,'K'), Tmax=(1194.97,'K')), NASAPolynomial(coeffs=[19.7249,0.0394515,-1.55411e-05,2.7927e-09,-1.90027e-13,-41601.9,-62.5912], Tmin=(1194.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-305.831,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C(O)C([O])(CC)C(C)O(7280)',
    structure = SMILES('[CH]=C(O)C([O])(CC)C(C)O'),
    E0 = (-196.539,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,180,1600,1650.8,2847.06,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152643,'amu*angstrom^2'), symmetry=1, barrier=(3.50955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152643,'amu*angstrom^2'), symmetry=1, barrier=(3.50955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152643,'amu*angstrom^2'), symmetry=1, barrier=(3.50955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152643,'amu*angstrom^2'), symmetry=1, barrier=(3.50955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152643,'amu*angstrom^2'), symmetry=1, barrier=(3.50955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152643,'amu*angstrom^2'), symmetry=1, barrier=(3.50955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152643,'amu*angstrom^2'), symmetry=1, barrier=(3.50955,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.88608,0.119461,-0.000120824,6.25801e-08,-1.2745e-11,-23417.8,43.5222], Tmin=(100,'K'), Tmax=(1199.85,'K')), NASAPolynomial(coeffs=[24.2384,0.0323687,-1.19446e-05,2.08392e-09,-1.40005e-13,-29686.9,-87.2732], Tmin=(1199.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-196.539,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C(O)C1(CC)OC[C]1O(8926)',
    structure = SMILES('[CH2]C(O)C1(CC)OC[C]1O'),
    E0 = (-178.791,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.4446,0.115627,-0.000112207,4.94363e-08,-4.65776e-12,-21303.4,38.1995], Tmin=(100,'K'), Tmax=(841.89,'K')), NASAPolynomial(coeffs=[20.2061,0.0357575,-1.08775e-05,1.64298e-09,-1.00517e-13,-25763.9,-67.3666], Tmin=(841.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-178.791,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(C2CsJOH) + radical(CJCO)"""),
)

species(
    label = 'CCC1([O])[C](O)CCC1O(12886)',
    structure = SMILES('CCC1([O])[C](O)CCC1O'),
    E0 = (-269.085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.00306,0.0967158,-6.4858e-05,1.27322e-08,2.81094e-12,-32171.8,37.6894], Tmin=(100,'K'), Tmax=(1033.08,'K')), NASAPolynomial(coeffs=[20.44,0.037739,-1.41441e-05,2.53926e-09,-1.75471e-13,-37885.6,-72.6702], Tmin=(1033.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-269.085,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(515.497,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(C2CsJOH) + radical(CC(C)2OJ)"""),
)

species(
    label = 'C=C(O)C(O)(CC)C(=C)O(12915)',
    structure = SMILES('C=C(O)C(O)(CC)C(=C)O'),
    E0 = (-593.25,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.91777,0.111925,-7.95199e-05,4.99915e-09,1.11732e-11,-71121.9,38.7316], Tmin=(100,'K'), Tmax=(956.038,'K')), NASAPolynomial(coeffs=[29.6404,0.0227345,-6.80747e-06,1.17303e-09,-8.47046e-14,-79114.2,-122.34], Tmin=(956.038,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-593.25,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
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
    label = '[CH2]C(O)C(C)([O])C(=C)O(12916)',
    structure = SMILES('[CH2]C(O)C(C)([O])C(=C)O'),
    E0 = (-208.266,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,180,180,180,180,1600,1623.74,2870.13,3200],'cm^-1')),
        HinderedRotor(inertia=(0.150955,'amu*angstrom^2'), symmetry=1, barrier=(3.47076,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150955,'amu*angstrom^2'), symmetry=1, barrier=(3.47076,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150955,'amu*angstrom^2'), symmetry=1, barrier=(3.47076,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150955,'amu*angstrom^2'), symmetry=1, barrier=(3.47076,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150955,'amu*angstrom^2'), symmetry=1, barrier=(3.47076,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150955,'amu*angstrom^2'), symmetry=1, barrier=(3.47076,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.22773,0.105305,-0.000106596,5.20852e-08,-9.20502e-12,-24851.9,39.1432], Tmin=(100,'K'), Tmax=(1019.76,'K')), NASAPolynomial(coeffs=[22.9078,0.0246454,-8.56016e-06,1.46783e-09,-9.89884e-14,-30503,-81.3409], Tmin=(1019.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-208.266,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CJCO)"""),
)

species(
    label = 'C=C(O)C([O])(CC)C[CH]O(12313)',
    structure = SMILES('C=C(O)C([O])(CC)C[CH]O'),
    E0 = (-250.695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.75922,0.123069,-0.000133041,7.49908e-08,-1.6757e-11,-29941.1,42.5599], Tmin=(100,'K'), Tmax=(1091.3,'K')), NASAPolynomial(coeffs=[21.7533,0.0368872,-1.45822e-05,2.62513e-09,-1.79012e-13,-35072.9,-72.9286], Tmin=(1091.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-250.695,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCsJOH)"""),
)

species(
    label = 'C=C(O)C1(CC)OCC1O(12320)',
    structure = SMILES('C=C(O)C1(CC)OCC1O'),
    E0 = (-488.706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.0554,0.0912429,-2.24867e-05,-5.75134e-08,3.66814e-11,-58576.8,35.6453], Tmin=(100,'K'), Tmax=(886.665,'K')), NASAPolynomial(coeffs=[26.8325,0.0224515,-2.57028e-06,3.84193e-11,5.12224e-15,-65763.6,-108.182], Tmin=(886.665,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-488.706,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Oxetane)"""),
)

species(
    label = '[CH2][CH]C(CC)(OO)C(=C)O(12917)',
    structure = SMILES('[CH2][CH]C(CC)(OO)C(=C)O'),
    E0 = (-23.1014,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.02366,0.123695,-0.000128077,6.75872e-08,-1.40173e-11,-2554.08,44.8943], Tmin=(100,'K'), Tmax=(1178.43,'K')), NASAPolynomial(coeffs=[24.7173,0.0329266,-1.25403e-05,2.22494e-09,-1.50936e-13,-8856.53,-88.5056], Tmin=(1178.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-23.1014,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJ) + radical(CCJCOOH)"""),
)

species(
    label = '[CH2]C(O)C([C]=C)(CC)OO(12918)',
    structure = SMILES('[CH2]C(O)C([C]=C)(CC)OO'),
    E0 = (62.7563,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3615,1310,387.5,850,1000,180,180,180,296.152,1433.79,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.147328,'amu*angstrom^2'), symmetry=1, barrier=(3.38736,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147328,'amu*angstrom^2'), symmetry=1, barrier=(3.38736,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147328,'amu*angstrom^2'), symmetry=1, barrier=(3.38736,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147328,'amu*angstrom^2'), symmetry=1, barrier=(3.38736,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147328,'amu*angstrom^2'), symmetry=1, barrier=(3.38736,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147328,'amu*angstrom^2'), symmetry=1, barrier=(3.38736,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147328,'amu*angstrom^2'), symmetry=1, barrier=(3.38736,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147328,'amu*angstrom^2'), symmetry=1, barrier=(3.38736,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.55674,0.124229,-0.00013901,8.32228e-08,-2.00628e-11,7746.35,42.5458], Tmin=(100,'K'), Tmax=(1005.52,'K')), NASAPolynomial(coeffs=[18.7045,0.0436293,-1.87754e-05,3.50696e-09,-2.43356e-13,3671.71,-55.3145], Tmin=(1005.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(62.7563,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CJCO)"""),
)

species(
    label = '[CH2]C(O)C([O])(CC)C(C)=O(12919)',
    structure = SMILES('[CH2]C(O)C([O])(CC)C(C)=O'),
    E0 = (-228.442,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.29178,0.12302,-0.000155376,1.12223e-07,-3.31642e-11,-27290.8,41.93], Tmin=(100,'K'), Tmax=(821.85,'K')), NASAPolynomial(coeffs=[13.7074,0.0500166,-2.21329e-05,4.13715e-09,-2.84885e-13,-29756.1,-27.4892], Tmin=(821.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-228.442,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(CC(C)(C=O)OJ) + radical(CJCO)"""),
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
    label = '[CH2]C(O)[C](CC)C(=C)O(12920)',
    structure = SMILES('[CH2]C(O)=C(CC)C([CH2])O'),
    E0 = (-131.107,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.92039,0.113098,-0.000112801,5.65514e-08,-1.09641e-11,-15540.9,41.4561], Tmin=(100,'K'), Tmax=(1321.31,'K')), NASAPolynomial(coeffs=[26.389,0.0236872,-7.08644e-06,1.08833e-09,-6.80396e-14,-22698.1,-101.782], Tmin=(1321.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-131.107,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(CJCO) + radical(C=C(O)CJ)"""),
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
    E0 = (-232.046,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-184.363,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-169.286,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-141.066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-180.853,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-182.767,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-90.4762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-232.046,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-37.0138,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-129.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-155.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-156.766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-74.5967,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-145.438,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-172.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-125.147,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-130.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-169.744,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-179.328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-189.37,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-163.499,'kJ/mol'),
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
    E0 = (-107.154,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-116.986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-168.646,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (211.596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-73.4196,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-223.762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (84.0091,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (169.867,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-27.8674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (111.897,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (154.649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(O)C([O])(CC)C(=C)O(7272)'],
    products = ['CH2CHOH(42)', 'C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(O)C([O])(CC)C(=C)O(7272)'],
    products = ['[CH2]C(O)C1(CC)OC1([CH2])O(12907)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.52e+08,'s^-1'), n=0.89, Ea=(47.6831,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H_secNd;radadd_intra] for rate rule [R4_S_D;doublebond_intra_2H_secNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(O)C([O])(CC)C(=C)O(7272)'],
    products = ['[CH2]C1(O)CC(O)C1([O])CC(12846)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.48e+06,'s^-1'), n=1.25, Ea=(62.76,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 349 used for R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C=C(O)C([O])(CC)C(=C)O(12908)'],
    products = ['[CH2]C(O)C([O])(CC)C(=C)O(7272)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(341.282,'m^3/(mol*s)'), n=1.56204, Ea=(11.2897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;HJ] for rate rule [Cds-OsCs_Cds;HJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C2H5(29)', '[CH2]C(O)C(=O)C(=C)O(12909)'],
    products = ['[CH2]C(O)C([O])(CC)C(=C)O(7272)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7.94e+10,'cm^3/(mol*s)'), n=0, Ea=(28.0328,'kJ/mol'), T0=(1,'K'), Tmin=(333,'K'), Tmax=(363,'K'), comment="""Estimated using template [CO_O;CsJ-CsHH] for rate rule [CO-CdCs_O;CsJ-CsHH]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2][CH]O(284)', 'C=C(O)C(=O)CC(4626)'],
    products = ['[CH2]C(O)C([O])(CC)C(=C)O(7272)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.97e+07,'cm^3/(mol*s)'), n=1.88, Ea=(32.2168,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO-CdCs_O;YJ] for rate rule [CO-CdCs_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2COH(99)', '[CH2]C(O)C(=O)CC(5567)'],
    products = ['[CH2]C(O)C([O])(CC)C(=C)O(7272)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.16e+10,'cm^3/(mol*s)'), n=0, Ea=(48.1578,'kJ/mol'), T0=(1,'K'), Tmin=(413,'K'), Tmax=(563,'K'), comment="""Estimated using template [CO-CsCs_O;CJ] for rate rule [CO-CsCs_O;CdsJ-O2s]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH2CHOH(42)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['[CH2]C(O)C([O])(CC)C(=C)O(7272)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.30847,'m^3/(mol*s)'), n=2.06429, Ea=(91.622,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cds_Cds;CJ] + [Cds-OsH_Cds;YJ] for rate rule [Cds-OsH_Cds;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 86.8 to 91.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['OH(5)', 'C=CC([O])(CC)C(=C)O(6516)'],
    products = ['[CH2]C(O)C([O])(CC)C(=C)O(7272)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.00674586,'m^3/(mol*s)'), n=2.3625, Ea=(84.2031,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;OJ] for rate rule [Cds-CsH_Cds-HH;OJ_pri]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(O)C([O])(CC)C(=C)O(7272)'],
    products = ['C=C(O)C([O])(CC)[C](C)O(7269)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.265e-07,'s^-1'), n=5.639, Ea=(102.68,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_2H;Cs_H_out_NonDe] for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_NDMustO]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(O)C(O)([CH]C)C(=C)O(12910)'],
    products = ['[CH2]C(O)C([O])(CC)C(=C)O(7272)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.71,'s^-1'), n=3.021, Ea=(105.562,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 319 used for R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(O)C([O])(CC)C(=C)O(7272)'],
    products = ['[CH2][C](O)C(O)(CC)C(=C)O(12911)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C(O)C([O])(CC)C(C)[O](6325)'],
    products = ['[CH2]C(O)C([O])(CC)C(=C)O(7272)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(O)C([O])(CC)C(=C)O(7272)'],
    products = ['C=C(O)C([O])([CH]C)C(C)O(7273)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(6.44e+09,'s^-1'), n=0.13, Ea=(86.6088,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 131 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]CC(O)(C(=C)O)C([CH2])O(12912)'],
    products = ['[CH2]C(O)C([O])(CC)C(=C)O(7272)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C([O])C(O)(CC)C(=C)O(7275)'],
    products = ['[CH2]C(O)C([O])(CC)C(=C)O(7272)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(94.0113,'s^-1'), n=2.81534, Ea=(105.641,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;O_H_out] + [R4H_SSS;O_rad_out;XH_out] for rate rule [R4H_SSS;O_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(O)C([O])(CC)C(=C)O(7272)'],
    products = ['[CH2]C(O)C(O)(CC)C(=C)[O](12913)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1070.11,'s^-1'), n=2.50856, Ea=(101.808,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;O_H_out] + [R4H_SSS;O_rad_out;XH_out] + [R4H_SS(Cd)S;Y_rad_out;XH_out] for rate rule [R4H_SS(Cd)S;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=C(O)C(O)(CC)C([CH2])O(12914)'],
    products = ['[CH2]C(O)C([O])(CC)C(=C)O(7272)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(O)C([O])(CC)C(=C)O(7272)'],
    products = ['[CH2]CC([O])(C(=C)O)C(C)O(7278)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(68850,'s^-1'), n=1.68, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 111 used for R5H_CCC;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R5H_CCC;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(O)C([O])(CC)C(=C)O(7272)'],
    products = ['C=C([O])C([O])(CC)C(C)O(7279)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2300,'s^-1'), n=1.98, Ea=(42.6768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC(Cd);C_rad_out_2H;XH_out] for rate rule [R5H_CCC(Cd);C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=C(O)C([O])(CC)C(C)O(7280)'],
    products = ['[CH2]C(O)C([O])(CC)C(=C)O(7272)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5H;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][CH]O(284)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['[CH2]C(O)C([O])(CC)C(=C)O(7272)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C(O)C([O])(CC)C(=C)O(7272)'],
    products = ['[CH2]C(O)C1(CC)OC[C]1O(8926)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.21748e+08,'s^-1'), n=0.95, Ea=(124.892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R4_S_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(O)C([O])(CC)C(=C)O(7272)'],
    products = ['CCC1([O])[C](O)CCC1O(12886)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.71e+11,'s^-1'), n=0.2, Ea=(115.06,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 4 used for R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(O)C([O])(CC)C(=C)O(7272)'],
    products = ['C=C(O)C(O)(CC)C(=C)O(12915)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CH2(S)(23)', '[CH2]C(O)C(C)([O])C(=C)O(12916)'],
    products = ['[CH2]C(O)C([O])(CC)C(=C)O(7272)'],
    transitionState = 'TS26',
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
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C(O)C([O])(CC)C(=C)O(7272)'],
    products = ['C=C(O)C1(CC)OCC1O(12320)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R4_SSS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][CH]C(CC)(OO)C(=C)O(12917)'],
    products = ['[CH2]C(O)C([O])(CC)C(=C)O(7272)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(3.01978e+11,'s^-1'), n=0, Ea=(107.111,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OOH_S;Y_rad_out]
Euclidian distance = 0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(O)C([C]=C)(CC)OO(12918)'],
    products = ['[CH2]C(O)C([O])(CC)C(=C)O(7272)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(3.01978e+11,'s^-1'), n=0, Ea=(107.111,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2OOH_S;Y_rad_out] for rate rule [R2OOH_S;Cd_rad_out_double]
Euclidian distance = 2.0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(O)C([O])(CC)C(=C)O(7272)'],
    products = ['[CH2]C(O)C([O])(CC)C(C)=O(12919)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(204.179,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction32',
    reactants = ['O(4)', '[CH2]C(O)[C](CC)C(=C)O(12920)'],
    products = ['[CH2]C(O)C([O])(CC)C(=C)O(7272)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/Cs2;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction33',
    reactants = ['CH2(19)', 'C=C(O)C([O])([CH]O)CC(6171)'],
    products = ['[CH2]C(O)C([O])(CC)C(=C)O(7272)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/CsO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '2245',
    isomers = [
        '[CH2]C(O)C([O])(CC)C(=C)O(7272)',
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
    label = '2245',
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

