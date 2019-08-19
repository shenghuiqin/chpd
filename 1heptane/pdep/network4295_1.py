species(
    label = 'C=C(O)C([O])(CC)C[C]=CC(20179)',
    structure = SMILES('C=C(O)C([O])(CC)C[C]=CC'),
    E0 = (34.6531,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,180,1085.24,1133.56,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.162437,'amu*angstrom^2'), symmetry=1, barrier=(3.73475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162437,'amu*angstrom^2'), symmetry=1, barrier=(3.73475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162437,'amu*angstrom^2'), symmetry=1, barrier=(3.73475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162437,'amu*angstrom^2'), symmetry=1, barrier=(3.73475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162437,'amu*angstrom^2'), symmetry=1, barrier=(3.73475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162437,'amu*angstrom^2'), symmetry=1, barrier=(3.73475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162437,'amu*angstrom^2'), symmetry=1, barrier=(3.73475,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.94699,0.123714,-0.000112923,5.42935e-08,-1.0529e-11,4387.81,45.9784], Tmin=(100,'K'), Tmax=(1236.57,'K')), NASAPolynomial(coeffs=[22.1353,0.0458134,-1.8426e-05,3.3475e-09,-2.29054e-13,-1568.05,-75.3184], Tmin=(1236.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(34.6531,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(Cds_S)"""),
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
    label = 'CH3CHCCH2(18175)',
    structure = SMILES('C=C=CC'),
    E0 = (145.615,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.759584,'amu*angstrom^2'), symmetry=1, barrier=(17.4643,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2996.71,'J/mol'), sigma=(5.18551,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=468.08 K, Pc=48.77 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.74635,0.0218189,8.22353e-06,-2.14768e-08,8.55624e-12,17563.6,12.7381], Tmin=(100,'K'), Tmax=(1025.6,'K')), NASAPolynomial(coeffs=[6.82078,0.0192338,-7.45622e-06,1.36536e-09,-9.53195e-14,16028,-10.4333], Tmin=(1025.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.615,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), label="""CH3CHCCH2""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2]C1(O)OC1(CC)C[C]=CC(26125)',
    structure = SMILES('[CH2]C1(O)OC1(CC)C[C]=CC'),
    E0 = (80.8555,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.88024,0.144669,-0.000147438,7.5563e-08,-1.4649e-11,10032.9,45.8432], Tmin=(100,'K'), Tmax=(1437.8,'K')), NASAPolynomial(coeffs=[32.0494,0.0261756,-4.48103e-06,3.1152e-10,-5.4654e-15,1616.9,-133.88], Tmin=(1437.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(80.8555,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Ethylene_oxide) + radical(Cds_S) + radical(CJC(O)2C)"""),
)

species(
    label = '[CH2]C1(O)C(=CC)CC1([O])CC(25640)',
    structure = SMILES('[CH2]C1(O)C(=CC)CC1([O])CC'),
    E0 = (47.0461,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[8.85565,0.0465354,6.31884e-05,-8.80389e-08,2.10357e-11,5348.73,-4.61778], Tmin=(100,'K'), Tmax=(1804.37,'K')), NASAPolynomial(coeffs=[102.358,0.0217628,-6.79402e-05,1.64672e-08,-1.21087e-12,-58103.9,-593.229], Tmin=(1804.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(47.0461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(methylenecyclobutane) + radical(C=CC(C)(O)CJ) + radical(CC(C)2OJ)"""),
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
    label = 'C=C(O)C([O])(C=C=CC)CC(26126)',
    structure = SMILES('C=C(O)C([O])(C=C=CC)CC'),
    E0 = (-45.0358,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,442.227,697.955,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156045,'amu*angstrom^2'), symmetry=1, barrier=(3.58778,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156045,'amu*angstrom^2'), symmetry=1, barrier=(3.58778,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156045,'amu*angstrom^2'), symmetry=1, barrier=(3.58778,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156045,'amu*angstrom^2'), symmetry=1, barrier=(3.58778,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156045,'amu*angstrom^2'), symmetry=1, barrier=(3.58778,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156045,'amu*angstrom^2'), symmetry=1, barrier=(3.58778,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (153.198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.6972,0.116485,-9.88895e-05,4.29918e-08,-7.53207e-12,-5204.33,42.4616], Tmin=(100,'K'), Tmax=(1356.48,'K')), NASAPolynomial(coeffs=[23.0108,0.0436253,-1.83201e-05,3.39408e-09,-2.34121e-13,-11907.4,-84.273], Tmin=(1356.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-45.0358,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(557.07,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C(O)C([O])(CC)CC#CC(26127)',
    structure = SMILES('C=C(O)C([O])(CC)CC#CC'),
    E0 = (-43.2408,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2100,2250,500,550,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,180,724.805,1482.02,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153382,'amu*angstrom^2'), symmetry=1, barrier=(3.52656,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153382,'amu*angstrom^2'), symmetry=1, barrier=(3.52656,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153382,'amu*angstrom^2'), symmetry=1, barrier=(3.52656,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153382,'amu*angstrom^2'), symmetry=1, barrier=(3.52656,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153382,'amu*angstrom^2'), symmetry=1, barrier=(3.52656,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153382,'amu*angstrom^2'), symmetry=1, barrier=(3.52656,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153382,'amu*angstrom^2'), symmetry=1, barrier=(3.52656,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (153.198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.26575,0.121741,-0.000113565,5.59337e-08,-1.08783e-11,-4961.48,46.1076], Tmin=(100,'K'), Tmax=(1276.38,'K')), NASAPolynomial(coeffs=[24.5545,0.0363189,-1.15658e-05,1.81674e-09,-1.13745e-13,-11696.4,-89.3921], Tmin=(1276.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-43.2408,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CtHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(C=CC(C)2OJ)"""),
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
    label = 'C=C(O)C(=O)C[C]=CC(26128)',
    structure = SMILES('C=C(O)C(=O)C[C]=CC'),
    E0 = (-46.4671,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1277.5,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.101428,0.0850763,-7.06743e-05,2.89774e-08,-4.82747e-12,-5448.33,29.952], Tmin=(100,'K'), Tmax=(1396.19,'K')), NASAPolynomial(coeffs=[17.5431,0.0351072,-1.69901e-05,3.34386e-09,-2.37563e-13,-10318.7,-60.015], Tmin=(1396.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-46.4671,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C=[C][CH]C(18176)',
    structure = SMILES('[CH2][C]=CC'),
    E0 = (361.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.352622,'amu*angstrom^2'), symmetry=1, barrier=(8.10748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.828631,'amu*angstrom^2'), symmetry=1, barrier=(19.0519,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.42015,0.030446,-1.69076e-05,4.64684e-09,-5.12013e-13,43485.7,14.8304], Tmin=(100,'K'), Tmax=(2065.83,'K')), NASAPolynomial(coeffs=[10.7464,0.014324,-5.20136e-06,8.69079e-10,-5.48385e-14,40045.6,-31.3799], Tmin=(2065.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Cds_S) + radical(Allyl_P)"""),
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
    label = 'CC=[C]CC(=O)CC(24960)',
    structure = SMILES('CC=[C]CC(=O)CC'),
    E0 = (28.3898,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,375,552.5,462.5,1710,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.162,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47644,0.0676262,-3.85299e-05,9.66829e-09,-9.50797e-13,3491.76,25.8592], Tmin=(100,'K'), Tmax=(2184.48,'K')), NASAPolynomial(coeffs=[18.1228,0.037145,-1.75996e-05,3.28071e-09,-2.19779e-13,-3780.98,-67.4571], Tmin=(2184.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.3898,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S)"""),
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
    label = 'C#CCC([O])(CC)C(=C)O(26129)',
    structure = SMILES('C#CCC([O])(CC)C(=C)O'),
    E0 = (-1.05972,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,750,770,3400,2100,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,217.119,913.389,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.150245,'amu*angstrom^2'), symmetry=1, barrier=(3.45443,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150245,'amu*angstrom^2'), symmetry=1, barrier=(3.45443,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150245,'amu*angstrom^2'), symmetry=1, barrier=(3.45443,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150245,'amu*angstrom^2'), symmetry=1, barrier=(3.45443,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150245,'amu*angstrom^2'), symmetry=1, barrier=(3.45443,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150245,'amu*angstrom^2'), symmetry=1, barrier=(3.45443,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (139.172,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.41169,0.1085,-9.69197e-05,3.86469e-08,-3.90927e-12,76.7642,39.6173], Tmin=(100,'K'), Tmax=(976.107,'K')), NASAPolynomial(coeffs=[22.5256,0.0314078,-1.07239e-05,1.81807e-09,-1.21968e-13,-5596.77,-80.4119], Tmin=(976.107,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1.05972,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C(O)C([O])([CH]C=CC)CC(20184)',
    structure = SMILES('C=C(O)C([O])([CH]C=CC)CC'),
    E0 = (-133.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.84236,0.107144,-4.616e-05,-2.43036e-08,1.88708e-11,-15772.3,45.9259], Tmin=(100,'K'), Tmax=(989.156,'K')), NASAPolynomial(coeffs=[26.8585,0.0385403,-1.4095e-05,2.59014e-09,-1.85441e-13,-23772,-103.962], Tmin=(989.156,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-133.045,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CCJC(O)C=C) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C(O)C([O])(CC)CC=[C]C(20176)',
    structure = SMILES('C=C(O)C([O])(CC)CC=[C]C'),
    E0 = (34.6531,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,180,1085.24,1133.56,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.162437,'amu*angstrom^2'), symmetry=1, barrier=(3.73475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162437,'amu*angstrom^2'), symmetry=1, barrier=(3.73475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162437,'amu*angstrom^2'), symmetry=1, barrier=(3.73475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162437,'amu*angstrom^2'), symmetry=1, barrier=(3.73475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162437,'amu*angstrom^2'), symmetry=1, barrier=(3.73475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162437,'amu*angstrom^2'), symmetry=1, barrier=(3.73475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162437,'amu*angstrom^2'), symmetry=1, barrier=(3.73475,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.94699,0.123714,-0.000112923,5.42936e-08,-1.0529e-11,4387.81,45.9785], Tmin=(100,'K'), Tmax=(1236.56,'K')), NASAPolynomial(coeffs=[22.1353,0.0458134,-1.84261e-05,3.3475e-09,-2.29054e-13,-1568.04,-75.3182], Tmin=(1236.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(34.6531,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C(O)C(O)([CH]C)C[C]=CC(26130)',
    structure = SMILES('C=C(O)C(O)([CH]C)C[C]=CC'),
    E0 = (5.45255,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,548.136,592.708,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.158793,'amu*angstrom^2'), symmetry=1, barrier=(3.65096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158793,'amu*angstrom^2'), symmetry=1, barrier=(3.65096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158793,'amu*angstrom^2'), symmetry=1, barrier=(3.65096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158793,'amu*angstrom^2'), symmetry=1, barrier=(3.65096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158793,'amu*angstrom^2'), symmetry=1, barrier=(3.65096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158793,'amu*angstrom^2'), symmetry=1, barrier=(3.65096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158793,'amu*angstrom^2'), symmetry=1, barrier=(3.65096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158793,'amu*angstrom^2'), symmetry=1, barrier=(3.65096,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.56184,0.129814,-0.000124977,6.19313e-08,-1.20766e-11,904.291,49.4695], Tmin=(100,'K'), Tmax=(1252.22,'K')), NASAPolynomial(coeffs=[26.7831,0.0360772,-1.26925e-05,2.15294e-09,-1.42155e-13,-6445.03,-98.7035], Tmin=(1252.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.45255,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCJCO) + radical(Cds_S)"""),
)

species(
    label = 'C=C(O)C(O)([CH][C]=CC)CC(26131)',
    structure = SMILES('C=C(O)C(O)([CH][C]=CC)CC'),
    E0 = (-124.306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.23224,0.117769,-7.33342e-05,-1.41096e-09,1.24727e-11,-14708.9,46.584], Tmin=(100,'K'), Tmax=(973.301,'K')), NASAPolynomial(coeffs=[28.5061,0.0346618,-1.18595e-05,2.10394e-09,-1.48536e-13,-22739.5,-111.395], Tmin=(973.301,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-124.306,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = 'C=C[CH]CC([O])(CC)C(=C)O(19803)',
    structure = SMILES('[CH2]C=CCC([O])(CC)C(=C)O'),
    E0 = (-51.6895,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,180,180,180,180,902.467,1303.41,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157691,'amu*angstrom^2'), symmetry=1, barrier=(3.62563,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157691,'amu*angstrom^2'), symmetry=1, barrier=(3.62563,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157691,'amu*angstrom^2'), symmetry=1, barrier=(3.62563,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157691,'amu*angstrom^2'), symmetry=1, barrier=(3.62563,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157691,'amu*angstrom^2'), symmetry=1, barrier=(3.62563,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157691,'amu*angstrom^2'), symmetry=1, barrier=(3.62563,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157691,'amu*angstrom^2'), symmetry=1, barrier=(3.62563,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4873.87,'J/mol'), sigma=(8.13902,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=761.29 K, Pc=20.51 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.15688,0.119489,-8.97057e-05,2.5025e-08,6.07937e-13,-5981.34,46.5382], Tmin=(100,'K'), Tmax=(1037.02,'K')), NASAPolynomial(coeffs=[25.4516,0.0405121,-1.52672e-05,2.75569e-09,-1.91319e-13,-13186.9,-94.7938], Tmin=(1037.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-51.6895,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(Allyl_P)"""),
)

species(
    label = 'C=C(O)C([O])([CH]C)CC=CC(20186)',
    structure = SMILES('C=C(O)C([O])([CH]C)CC=CC'),
    E0 = (-3.28658,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,180,1069.08,1145.33,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.161939,'amu*angstrom^2'), symmetry=1, barrier=(3.7233,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161939,'amu*angstrom^2'), symmetry=1, barrier=(3.7233,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161939,'amu*angstrom^2'), symmetry=1, barrier=(3.7233,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161939,'amu*angstrom^2'), symmetry=1, barrier=(3.7233,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161939,'amu*angstrom^2'), symmetry=1, barrier=(3.7233,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161939,'amu*angstrom^2'), symmetry=1, barrier=(3.7233,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161939,'amu*angstrom^2'), symmetry=1, barrier=(3.7233,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.01577,0.11738,-9.16594e-05,3.13671e-08,-2.54215e-12,-165.946,48.2492], Tmin=(100,'K'), Tmax=(1077.28,'K')), NASAPolynomial(coeffs=[24.1627,0.041545,-1.58187e-05,2.84514e-09,-1.95877e-13,-7046.12,-85.7502], Tmin=(1077.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-3.28658,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCJCO) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH2]CC(O)(C[C]=CC)C(=C)O(26132)',
    structure = SMILES('[CH2]CC(O)(C[C]=CC)C(=C)O'),
    E0 = (10.7969,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,180,180,180,565.818,576.026,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.158327,'amu*angstrom^2'), symmetry=1, barrier=(3.64025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158327,'amu*angstrom^2'), symmetry=1, barrier=(3.64025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158327,'amu*angstrom^2'), symmetry=1, barrier=(3.64025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158327,'amu*angstrom^2'), symmetry=1, barrier=(3.64025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158327,'amu*angstrom^2'), symmetry=1, barrier=(3.64025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158327,'amu*angstrom^2'), symmetry=1, barrier=(3.64025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158327,'amu*angstrom^2'), symmetry=1, barrier=(3.64025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158327,'amu*angstrom^2'), symmetry=1, barrier=(3.64025,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.42569,0.132354,-0.000132389,6.89073e-08,-1.42109e-11,1537.54,48.268], Tmin=(100,'K'), Tmax=(1180.41,'K')), NASAPolynomial(coeffs=[24.8157,0.0400407,-1.50811e-05,2.65385e-09,-1.78847e-13,-4893.57,-87.6739], Tmin=(1180.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(10.7969,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(RCCJ)"""),
)

species(
    label = 'C=C([O])C(O)(CC)C[C]=CC(26133)',
    structure = SMILES('C=C([O])C(O)(CC)C[C]=CC'),
    E0 = (-56.6447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.89368,0.125247,-0.000121065,6.29043e-08,-1.32066e-11,-6596.84,46.074], Tmin=(100,'K'), Tmax=(1147.28,'K')), NASAPolynomial(coeffs=[20.5926,0.0468487,-1.85637e-05,3.34222e-09,-2.27573e-13,-11756.5,-65.4988], Tmin=(1147.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-56.6447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C(O)C(O)(CC)C[C]=CC(26134)',
    structure = SMILES('[CH]=C(O)C(O)(CC)C[C]=CC'),
    E0 = (52.6466,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,180,1600,1720.41,2783.48,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155143,'amu*angstrom^2'), symmetry=1, barrier=(3.56704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155143,'amu*angstrom^2'), symmetry=1, barrier=(3.56704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155143,'amu*angstrom^2'), symmetry=1, barrier=(3.56704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155143,'amu*angstrom^2'), symmetry=1, barrier=(3.56704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155143,'amu*angstrom^2'), symmetry=1, barrier=(3.56704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155143,'amu*angstrom^2'), symmetry=1, barrier=(3.56704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155143,'amu*angstrom^2'), symmetry=1, barrier=(3.56704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155143,'amu*angstrom^2'), symmetry=1, barrier=(3.56704,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.51903,0.134939,-0.000137916,7.32182e-08,-1.53726e-11,6573.85,47.2194], Tmin=(100,'K'), Tmax=(1161.38,'K')), NASAPolynomial(coeffs=[25.1255,0.0397265,-1.49417e-05,2.62694e-09,-1.76981e-13,152.724,-90.2851], Tmin=(1161.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(52.6466,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]CC([O])(CC=CC)C(=C)O(20187)',
    structure = SMILES('[CH2]CC([O])(CC=CC)C(=C)O'),
    E0 = (2.05772,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,180,180,180,180,1028.3,1187.76,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.160992,'amu*angstrom^2'), symmetry=1, barrier=(3.70151,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160992,'amu*angstrom^2'), symmetry=1, barrier=(3.70151,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160992,'amu*angstrom^2'), symmetry=1, barrier=(3.70151,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160992,'amu*angstrom^2'), symmetry=1, barrier=(3.70151,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160992,'amu*angstrom^2'), symmetry=1, barrier=(3.70151,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160992,'amu*angstrom^2'), symmetry=1, barrier=(3.70151,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160992,'amu*angstrom^2'), symmetry=1, barrier=(3.70151,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.13795,0.122853,-0.000108728,4.9893e-08,-9.15153e-12,478.564,47.981], Tmin=(100,'K'), Tmax=(1312.08,'K')), NASAPolynomial(coeffs=[24.4056,0.0419319,-1.62173e-05,2.88821e-09,-1.95345e-13,-6486.87,-87.2858], Tmin=(1312.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.05772,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJ) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C([O])C([O])(CC)CC=CC(20188)',
    structure = SMILES('C=C([O])C([O])(CC)CC=CC'),
    E0 = (-65.3838,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.56269,0.11531,-9.61711e-05,4.26251e-08,-7.72124e-12,-7657.9,45.6264], Tmin=(100,'K'), Tmax=(1305.68,'K')), NASAPolynomial(coeffs=[20.1729,0.0487234,-1.96768e-05,3.56884e-09,-2.43285e-13,-13334,-65.0333], Tmin=(1305.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-65.3838,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH]=C(O)C([O])(CC)CC=CC(20189)',
    structure = SMILES('[CH]=C(O)C([O])(CC)CC=CC'),
    E0 = (43.9075,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,528.184,614.111,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157884,'amu*angstrom^2'), symmetry=1, barrier=(3.63007,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157884,'amu*angstrom^2'), symmetry=1, barrier=(3.63007,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157884,'amu*angstrom^2'), symmetry=1, barrier=(3.63007,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157884,'amu*angstrom^2'), symmetry=1, barrier=(3.63007,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157884,'amu*angstrom^2'), symmetry=1, barrier=(3.63007,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157884,'amu*angstrom^2'), symmetry=1, barrier=(3.63007,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157884,'amu*angstrom^2'), symmetry=1, barrier=(3.63007,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.21001,0.12521,-0.000113555,5.34292e-08,-1.00358e-11,5513.91,46.8546], Tmin=(100,'K'), Tmax=(1283.91,'K')), NASAPolynomial(coeffs=[24.5284,0.0419078,-1.62341e-05,2.89636e-09,-1.96273e-13,-1352.12,-88.8253], Tmin=(1283.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(43.9075,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(Cds_P)"""),
)

species(
    label = 'C=C(O)C(O)(CC)C[C]=[C]C(26135)',
    structure = SMILES('C=C(O)C(O)(CC)C[C]=[C]C'),
    E0 = (43.3923,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,3580,3650,1210,1345,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,525.574,618.252,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15933,'amu*angstrom^2'), symmetry=1, barrier=(3.66331,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15933,'amu*angstrom^2'), symmetry=1, barrier=(3.66331,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15933,'amu*angstrom^2'), symmetry=1, barrier=(3.66331,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15933,'amu*angstrom^2'), symmetry=1, barrier=(3.66331,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15933,'amu*angstrom^2'), symmetry=1, barrier=(3.66331,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15933,'amu*angstrom^2'), symmetry=1, barrier=(3.66331,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15933,'amu*angstrom^2'), symmetry=1, barrier=(3.66331,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15933,'amu*angstrom^2'), symmetry=1, barrier=(3.66331,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.31136,0.134034,-0.000139078,7.6063e-08,-1.6579e-11,5450.28,46.5462], Tmin=(100,'K'), Tmax=(1115.15,'K')), NASAPolynomial(coeffs=[23.0092,0.0432112,-1.6911e-05,3.02894e-09,-2.05898e-13,-196.994,-78.3705], Tmin=(1115.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(43.3923,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=[C]CC(O)(CC)C(=C)O(20181)',
    structure = SMILES('[CH2]C=[C]CC(O)(CC)C(=C)O'),
    E0 = (-42.9503,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,180,180,180,399.146,736.022,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155125,'amu*angstrom^2'), symmetry=1, barrier=(3.56663,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155125,'amu*angstrom^2'), symmetry=1, barrier=(3.56663,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155125,'amu*angstrom^2'), symmetry=1, barrier=(3.56663,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155125,'amu*angstrom^2'), symmetry=1, barrier=(3.56663,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155125,'amu*angstrom^2'), symmetry=1, barrier=(3.56663,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155125,'amu*angstrom^2'), symmetry=1, barrier=(3.56663,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155125,'amu*angstrom^2'), symmetry=1, barrier=(3.56663,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155125,'amu*angstrom^2'), symmetry=1, barrier=(3.56663,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.89003,0.134059,-0.000130153,6.43323e-08,-1.24324e-11,-4902.88,48.4341], Tmin=(100,'K'), Tmax=(1268.41,'K')), NASAPolynomial(coeffs=[28.9781,0.0335639,-1.13119e-05,1.87187e-09,-1.21965e-13,-12987.4,-112.889], Tmin=(1268.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-42.9503,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'CC=[C]CC1(CC)OC[C]1O(23147)',
    structure = SMILES('CC=[C]CC1(CC)OC[C]1O'),
    E0 = (87.9086,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3010,987.5,1337.5,450,1655,1685,370,2750,3150,900,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.86548,0.123605,-0.00012124,6.7337e-08,-1.51701e-11,10788.8,41.7031], Tmin=(100,'K'), Tmax=(1075.88,'K')), NASAPolynomial(coeffs=[18.2156,0.0489455,-1.71497e-05,2.83764e-09,-1.82578e-13,6467.79,-56.6453], Tmin=(1075.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(87.9086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Oxetane) + radical(Cds_S) + radical(C2CsJOH)"""),
)

species(
    label = 'CC=C1C[C](O)C([O])(CC)C1(25796)',
    structure = SMILES('CC=C1C[C](O)C([O])(CC)C1'),
    E0 = (-47.8312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.64501,0.108449,-7.36e-05,2.16752e-08,-1.52882e-12,-5536.5,40.2178], Tmin=(100,'K'), Tmax=(1205.98,'K')), NASAPolynomial(coeffs=[22.101,0.0471145,-1.89867e-05,3.46663e-09,-2.38011e-13,-12531.2,-84.0443], Tmin=(1205.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-47.8312,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(590.328,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(methylenecyclopentane) + radical(C2CsJOH) + radical(CC(C)2OJ)"""),
)

species(
    label = 'C=C(O)C(O)(C=C=CC)CC(26136)',
    structure = SMILES('C=C(O)C(O)(C=C=CC)CC'),
    E0 = (-274.138,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.35603,0.126014,-0.000113311,5.20113e-08,-9.47482e-12,-32730.9,43.6152], Tmin=(100,'K'), Tmax=(1326.07,'K')), NASAPolynomial(coeffs=[26.3649,0.0393788,-1.5312e-05,2.74308e-09,-1.86359e-13,-40348.1,-103.052], Tmin=(1326.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-274.138,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
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
    label = 'C=C(O)C(C)([O])C[C]=CC(26137)',
    structure = SMILES('C=C(O)C(C)([O])C[C]=CC'),
    E0 = (58.4334,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,547.816,595.626,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.158691,'amu*angstrom^2'), symmetry=1, barrier=(3.64862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158691,'amu*angstrom^2'), symmetry=1, barrier=(3.64862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158691,'amu*angstrom^2'), symmetry=1, barrier=(3.64862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158691,'amu*angstrom^2'), symmetry=1, barrier=(3.64862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158691,'amu*angstrom^2'), symmetry=1, barrier=(3.64862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158691,'amu*angstrom^2'), symmetry=1, barrier=(3.64862,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.35444,0.109544,-0.000101114,4.87288e-08,-9.39823e-12,7227.56,41.6076], Tmin=(100,'K'), Tmax=(1249.95,'K')), NASAPolynomial(coeffs=[21.1484,0.0375313,-1.46947e-05,2.63615e-09,-1.79245e-13,1602.12,-71.9759], Tmin=(1249.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(58.4334,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=[C]CC([O])(CC)C(=C)O(26119)',
    structure = SMILES('C=[C]CC([O])(CC)C(=C)O'),
    E0 = (70.6787,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,180,1036.49,1179.57,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.162038,'amu*angstrom^2'), symmetry=1, barrier=(3.72556,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162038,'amu*angstrom^2'), symmetry=1, barrier=(3.72556,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162038,'amu*angstrom^2'), symmetry=1, barrier=(3.72556,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162038,'amu*angstrom^2'), symmetry=1, barrier=(3.72556,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162038,'amu*angstrom^2'), symmetry=1, barrier=(3.72556,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162038,'amu*angstrom^2'), symmetry=1, barrier=(3.72556,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.5578,0.111128,-0.000102613,4.89134e-08,-9.25811e-12,8710.09,42.7781], Tmin=(100,'K'), Tmax=(1279.66,'K')), NASAPolynomial(coeffs=[22.8455,0.0348461,-1.31952e-05,2.32871e-09,-1.56992e-13,2464.6,-80.9709], Tmin=(1279.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(70.6787,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(Cds_S)"""),
)

species(
    label = 'C=C(O)C([O])(CC)C(=C)[CH]C(25632)',
    structure = SMILES('[CH2]C(=CC)C([O])(CC)C(=C)O'),
    E0 = (-73.1703,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,262.04,866.396,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15197,'amu*angstrom^2'), symmetry=1, barrier=(3.49409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15197,'amu*angstrom^2'), symmetry=1, barrier=(3.49409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15197,'amu*angstrom^2'), symmetry=1, barrier=(3.49409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15197,'amu*angstrom^2'), symmetry=1, barrier=(3.49409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15197,'amu*angstrom^2'), symmetry=1, barrier=(3.49409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15197,'amu*angstrom^2'), symmetry=1, barrier=(3.49409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15197,'amu*angstrom^2'), symmetry=1, barrier=(3.49409,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4889.28,'J/mol'), sigma=(8.15989,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=763.69 K, Pc=20.42 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.24133,0.122242,-0.00010271,4.37536e-08,-7.43059e-12,-8563.24,44.6026], Tmin=(100,'K'), Tmax=(1410.71,'K')), NASAPolynomial(coeffs=[26.3925,0.0410524,-1.63822e-05,2.95775e-09,-2.00988e-13,-16642.1,-103.392], Tmin=(1410.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-73.1703,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(Allyl_P)"""),
)

species(
    label = 'C=C(O)C1(CC)CC(=CC)O1(26047)',
    structure = SMILES('C=C(O)C1(CC)CC(=CC)O1'),
    E0 = (-311.179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.86538,0.100295,-5.96529e-06,-8.80266e-08,4.88019e-11,-37188.3,36.8255], Tmin=(100,'K'), Tmax=(915.082,'K')), NASAPolynomial(coeffs=[33.1866,0.0227904,-3.03084e-06,2.54342e-10,-1.87432e-14,-46773.5,-146.491], Tmin=(915.082,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-311.179,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(2methyleneoxetane)"""),
)

species(
    label = 'C=[C]C(CC)(C[C]=CC)OO(26138)',
    structure = SMILES('C=[C]C(CC)(C[C]=CC)OO'),
    E0 = (329.456,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,317.953,825.807,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153831,'amu*angstrom^2'), symmetry=1, barrier=(3.53689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153831,'amu*angstrom^2'), symmetry=1, barrier=(3.53689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153831,'amu*angstrom^2'), symmetry=1, barrier=(3.53689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153831,'amu*angstrom^2'), symmetry=1, barrier=(3.53689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153831,'amu*angstrom^2'), symmetry=1, barrier=(3.53689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153831,'amu*angstrom^2'), symmetry=1, barrier=(3.53689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153831,'amu*angstrom^2'), symmetry=1, barrier=(3.53689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153831,'amu*angstrom^2'), symmetry=1, barrier=(3.53689,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.53229,0.126656,-0.000127158,7.17294e-08,-1.69004e-11,39819.5,44.471], Tmin=(100,'K'), Tmax=(1008.19,'K')), NASAPolynomial(coeffs=[15.5681,0.0588102,-2.62153e-05,4.98067e-09,-3.48762e-13,36371.4,-38.1681], Tmin=(1008.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.456,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'CC=[C]CC([O])(CC)C(C)=O(26139)',
    structure = SMILES('CC=[C]CC([O])(CC)C(C)=O'),
    E0 = (38.2571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.389,0.126908,-0.000148802,1.07958e-07,-3.33341e-11,4787.63,44.2901], Tmin=(100,'K'), Tmax=(776.668,'K')), NASAPolynomial(coeffs=[10.5991,0.0651682,-2.95645e-05,5.61045e-09,-3.9036e-13,2925.44,-10.5158], Tmin=(776.668,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(38.2571,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)(C=O)OJ) + radical(Cds_S)"""),
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
    label = 'C=C(O)[C](CC)C[C]=CC(26140)',
    structure = SMILES('[CH2]C(O)=C(CC)C[C]=CC'),
    E0 = (137.198,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1277.5,1000,1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
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
    molecularWeight = (138.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.74205,0.116371,-0.000105478,5.04117e-08,-9.6573e-12,16716.1,41.9372], Tmin=(100,'K'), Tmax=(1258.51,'K')), NASAPolynomial(coeffs=[22.0275,0.0408234,-1.54348e-05,2.71368e-09,-1.82285e-13,10733.2,-78.2027], Tmin=(1258.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(137.198,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(C=C(O)CJ)"""),
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
    label = '[C]=CC(24199)',
    structure = SMILES('[C]=CC'),
    E0 = (564.071,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.40488,'amu*angstrom^2'), symmetry=1, barrier=(9.30899,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28451,0.0144722,-2.9151e-06,-2.04635e-09,7.91551e-13,67868.8,10.0416], Tmin=(100,'K'), Tmax=(1455.86,'K')), NASAPolynomial(coeffs=[5.67821,0.0121697,-4.94658e-06,9.00486e-10,-6.07653e-14,66718.9,-3.96135], Tmin=(1455.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(564.071,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CdCdJ2_triplet)"""),
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
    E0 = (34.6531,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (82.3363,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (97.4131,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (175.25,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (187.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (89.4398,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (34.8055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (179.817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (34.6531,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (165.253,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (234.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (268.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (111.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (109.934,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (196.574,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (78.9617,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (94.4769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (136.461,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (96.9552,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (70.4556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (77.3299,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (76.9477,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (220.634,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (128.166,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (176.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (159.546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (77.7483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (112.9,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (478.295,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (490.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (204.691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (42.9375,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (436.566,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (238.832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (380.203,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (530.13,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    products = ['C=C(O)C(=O)CC(4626)', 'CH3CHCCH2(18175)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    products = ['[CH2]C1(O)OC1(CC)C[C]=CC(26125)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.52e+08,'s^-1'), n=0.89, Ea=(47.6831,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H_secNd;radadd_intra] for rate rule [R4_S_D;doublebond_intra_2H_secNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    products = ['[CH2]C1(O)C(=CC)CC1([O])CC(25640)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.48e+06,'s^-1'), n=1.25, Ea=(62.76,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C=C(O)C([O])(C=C=CC)CC(26126)'],
    products = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.23e+08,'cm^3/(mol*s)'), n=1.64, Ea=(8.49352,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2566 used for Cds-CsH_Ca;HJ
Exact match found for rate rule [Cds-CsH_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C=C(O)C([O])(CC)CC#CC(26127)'],
    products = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.02e+09,'cm^3/(mol*s)'), n=1.64, Ea=(18.4933,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2702 used for Ct-Cs_Ct-Cs;HJ
Exact match found for rate rule [Ct-Cs_Ct-Cs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C2H5(29)', 'C=C(O)C(=O)C[C]=CC(26128)'],
    products = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7.94e+10,'cm^3/(mol*s)'), n=0, Ea=(28.0328,'kJ/mol'), T0=(1,'K'), Tmin=(333,'K'), Tmax=(363,'K'), comment="""Estimated using template [CO_O;CsJ-CsHH] for rate rule [CO-CdCs_O;CsJ-CsHH]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=C(O)C(=O)CC(4626)', 'C=[C][CH]C(18176)'],
    products = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.97e+07,'cm^3/(mol*s)'), n=1.88, Ea=(32.2168,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO-CdCs_O;YJ] for rate rule [CO-CdCs_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH2COH(99)', 'CC=[C]CC(=O)CC(24960)'],
    products = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.16e+10,'cm^3/(mol*s)'), n=0, Ea=(48.1578,'kJ/mol'), T0=(1,'K'), Tmin=(413,'K'), Tmax=(563,'K'), comment="""Estimated using template [CO-CsCs_O;CJ] for rate rule [CO-CsCs_O;CdsJ-O2s]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(O)=C([O])CC(4557)', 'CH3CHCCH2(18175)'],
    products = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.00429749,'m^3/(mol*s)'), n=2.45395, Ea=(73.9811,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Ca;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 69.2 to 74.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['CH3(17)', 'C#CCC([O])(CC)C(=C)O(26129)'],
    products = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(178000,'cm^3/(mol*s)'), n=2.41, Ea=(30.1248,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2271 used for Ct-H_Ct-Cs;CsJ-HHH
Exact match found for rate rule [Ct-H_Ct-Cs;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    products = ['C=C(O)C([O])([CH]C=CC)CC(20184)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.9054e+11,'s^-1'), n=0.853, Ea=(200.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=C(O)C([O])(CC)CC=[C]C(20176)'],
    products = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.231e+11,'s^-1'), n=0.765, Ea=(234.057,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 82 used for R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C(O)C(O)([CH]C)C[C]=CC(26130)'],
    products = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(5.71,'s^-1'), n=3.021, Ea=(105.562,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 319 used for R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    products = ['C=C(O)C(O)([CH][C]=CC)CC(26131)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(223829,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    products = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    products = ['C=C(O)C([O])([CH]C)CC=CC(20186)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out_1H] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]CC(O)(C[C]=CC)C(=C)O(26132)'],
    products = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    products = ['C=C([O])C(O)(CC)C[C]=CC(26133)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1070.11,'s^-1'), n=2.50856, Ea=(101.808,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;O_H_out] + [R4H_SSS;O_rad_out;XH_out] + [R4H_SS(Cd)S;Y_rad_out;XH_out] for rate rule [R4H_SS(Cd)S;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C(O)C(O)(CC)C[C]=CC(26134)'],
    products = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    products = ['[CH2]CC([O])(CC=CC)C(=C)O(20187)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(561575,'s^-1'), n=1.6076, Ea=(35.8025,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;Y_rad_out;Cs_H_out_2H] for rate rule [R5H_CCC;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    products = ['C=C([O])C([O])(CC)CC=CC(20188)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2300,'s^-1'), n=1.98, Ea=(42.6768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC(Cd);Y_rad_out;XH_out] for rate rule [R5H_CCC(Cd);Cd_rad_out_Cd;O_H_out]
Euclidian distance = 3.16227766017
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C(O)C([O])(CC)CC=CC(20189)'],
    products = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_singleH;XH_out] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 3.60555127546
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C(O)C(O)(CC)C[C]=[C]C(26135)'],
    products = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(8.44313e+09,'s^-1'), n=0.985167, Ea=(177.241,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cs;XH_out] for rate rule [R5HJ_1;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    products = ['[CH2]C=[C]CC(O)(CC)C(=C)O(20181)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.68e+09,'s^-1'), n=0, Ea=(93.5124,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6Hall;O_rad_out;Cs_H_out_2H] for rate rule [R6HJ_3;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(O)=C([O])CC(4557)', 'C=[C][CH]C(18176)'],
    products = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    products = ['CC=[C]CC1(CC)OC[C]1O(23147)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(3.21748e+08,'s^-1'), n=0.95, Ea=(124.892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R4_S_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    products = ['CC=C1C[C](O)C([O])(CC)C1(25796)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.6e+11,'s^-1'), n=0.27, Ea=(43.0952,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 8 used for R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_cddouble
Exact match found for rate rule [R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    products = ['C=C(O)C(O)(C=C=CC)CC(26136)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['CH2(S)(23)', 'C=C(O)C(C)([O])C[C]=CC(26137)'],
    products = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction30',
    reactants = ['CH2(S)(23)', 'C=[C]CC([O])(CC)C(=C)O(26119)'],
    products = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    products = ['C=C(O)C([O])(CC)C(=C)[CH]C(25632)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    products = ['C=C(O)C1(CC)CC(=CC)O1(26047)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=[C]C(CC)(C[C]=CC)OO(26138)'],
    products = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(3.01978e+11,'s^-1'), n=0, Ea=(107.111,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2OOH_S;Y_rad_out] for rate rule [R2OOH_S;Cd_rad_out_double]
Euclidian distance = 2.0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    products = ['CC=[C]CC([O])(CC)C(C)=O(26139)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(204.179,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction35',
    reactants = ['O(4)', 'C=C(O)[C](CC)C[C]=CC(26140)'],
    products = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/Cs2;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C([O])(CC)C(=C)O(6187)', '[C]=CC(24199)'],
    products = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4295',
    isomers = [
        'C=C(O)C([O])(CC)C[C]=CC(20179)',
    ],
    reactants = [
        ('C=C(O)C(=O)CC(4626)', 'CH3CHCCH2(18175)'),
        ('[CH2]C(O)=C([O])CC(4557)', 'CH3CHCCH2(18175)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4295',
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

