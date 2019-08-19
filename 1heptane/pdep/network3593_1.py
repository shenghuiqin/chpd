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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.15688,0.119489,-8.97057e-05,2.5025e-08,6.07937e-13,-5981.34,46.5382], Tmin=(100,'K'), Tmax=(1037.02,'K')), NASAPolynomial(coeffs=[25.4516,0.0405121,-1.52672e-05,2.75569e-09,-1.91319e-13,-13186.9,-94.7938], Tmin=(1037.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-51.6895,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(Allyl_P)"""),
)

species(
    label = 'butadiene13(2459)',
    structure = SMILES('C=CC=C'),
    E0 = (96.4553,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(1.30712,'amu*angstrom^2'), symmetry=1, barrier=(30.0532,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.18,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.80599,0.0102584,6.1726e-05,-9.01643e-08,3.59117e-11,11658.5,12.0621], Tmin=(100,'K'), Tmax=(946.047,'K')), NASAPolynomial(coeffs=[12.4694,0.0100554,-2.41207e-06,4.57077e-10,-3.93161e-14,8010.78,-43.6375], Tmin=(946.047,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(96.4553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), label="""butadiene13""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2]C=CCC1(CC)OC1([CH2])O(20172)',
    structure = SMILES('[CH2]C=CCC1(CC)OC1([CH2])O'),
    E0 = (-5.48714,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.64037,0.146544,-0.000143857,6.95146e-08,-1.25006e-11,-311.452,48.4046], Tmin=(100,'K'), Tmax=(1605.37,'K')), NASAPolynomial(coeffs=[34.1367,0.022247,-1.85688e-06,-1.93888e-10,2.73939e-14,-9195.08,-145.919], Tmin=(1605.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-5.48714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Ethylene_oxide) + radical(Allyl_P) + radical(CJC(O)2C)"""),
)

species(
    label = '[CH2][C](O)C1(CC)CC(C=C)O1(19941)',
    structure = SMILES('[CH2][C](O)C1(CC)CC(C=C)O1'),
    E0 = (60.922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.7203,0.127144,-0.000117733,5.7631e-08,-1.10276e-11,7586.75,45.2462], Tmin=(100,'K'), Tmax=(1379.13,'K')), NASAPolynomial(coeffs=[25.9116,0.0361551,-1.01271e-05,1.43685e-09,-8.37833e-14,444.978,-99.3504], Tmin=(1379.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(60.922,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Oxetane) + radical(CJCO) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C1(O)CC=CCC1([O])CC(19823)',
    structure = SMILES('[CH2]C1(O)CC=CCC1([O])CC'),
    E0 = (-32.5857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.91028,0.112845,-7.13855e-05,8.97001e-09,5.3574e-12,-3691.54,37.9232], Tmin=(100,'K'), Tmax=(1040.64,'K')), NASAPolynomial(coeffs=[24.3451,0.0434348,-1.67545e-05,3.06809e-09,-2.1474e-13,-10862.2,-97.9871], Tmin=(1040.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-32.5857,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(590.328,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclohexene) + radical(CC(C)2OJ) + radical(CJC(C)2O)"""),
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
    label = 'C=C=CCC([O])(CC)C(=C)O(20173)',
    structure = SMILES('C=C=CCC([O])(CC)C(=C)O'),
    E0 = (-26.5843,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,180,997.07,1217.64,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.1609,'amu*angstrom^2'), symmetry=1, barrier=(3.6994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1609,'amu*angstrom^2'), symmetry=1, barrier=(3.6994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1609,'amu*angstrom^2'), symmetry=1, barrier=(3.6994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1609,'amu*angstrom^2'), symmetry=1, barrier=(3.6994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1609,'amu*angstrom^2'), symmetry=1, barrier=(3.6994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1609,'amu*angstrom^2'), symmetry=1, barrier=(3.6994,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (153.198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.26623,0.123263,-0.000112409,5.24633e-08,-9.68717e-12,-2959.56,45.7108], Tmin=(100,'K'), Tmax=(1313.41,'K')), NASAPolynomial(coeffs=[26.0746,0.0369501,-1.38337e-05,2.42775e-09,-1.63136e-13,-10404.1,-98.7436], Tmin=(1313.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-26.5843,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(557.07,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CC(C)2OJ)"""),
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
    label = '[CH2]C=CCC(=O)C(=C)O(20174)',
    structure = SMILES('[CH2]C=CCC(=O)C(=C)O'),
    E0 = (-132.81,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.528291,0.0855178,-6.25504e-05,1.76891e-08,-7.22871e-13,-15798.6,32.0382], Tmin=(100,'K'), Tmax=(1196.03,'K')), NASAPolynomial(coeffs=[22.1383,0.0277317,-1.26775e-05,2.48695e-09,-1.78307e-13,-22509.4,-86.7608], Tmin=(1196.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-132.81,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]C=C(2458)',
    structure = SMILES('[CH2]C=C[CH2]'),
    E0 = (274.714,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.210311,'amu*angstrom^2'), symmetry=1, barrier=(25.2351,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.779031,'amu*angstrom^2'), symmetry=1, barrier=(93.4717,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.56318,0.0223429,1.87067e-05,-3.93099e-08,1.63982e-11,33100.5,13.4097], Tmin=(100,'K'), Tmax=(974.264,'K')), NASAPolynomial(coeffs=[9.82995,0.0151966,-5.22272e-06,9.67656e-10,-7.0786e-14,30607.7,-26.985], Tmin=(974.264,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Allyl_P) + radical(Allyl_P)"""),
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
    label = '[CH2]C=CCC(=O)CC(20175)',
    structure = SMILES('[CH2]C=CCC(=O)CC'),
    E0 = (-57.9528,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,375,552.5,462.5,1710,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.268571,0.0736198,-4.51714e-05,1.26555e-08,-1.38487e-12,-6828.48,30.1234], Tmin=(100,'K'), Tmax=(2112.75,'K')), NASAPolynomial(coeffs=[25.0929,0.026621,-1.18035e-05,2.1265e-09,-1.38984e-13,-17318,-108.208], Tmin=(2112.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-57.9528,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P)"""),
)

species(
    label = 'C=C(O)C([O])(CC)CC=[C]C(20176)',
    structure = SMILES('C=C(O)C([O])(CC)CC=[C]C'),
    E0 = (34.6531,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.94699,0.123714,-0.000112923,5.42936e-08,-1.0529e-11,4387.81,45.9785], Tmin=(100,'K'), Tmax=(1236.56,'K')), NASAPolynomial(coeffs=[22.1353,0.0458134,-1.84261e-05,3.3475e-09,-2.29054e-13,-1568.04,-75.3182], Tmin=(1236.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(34.6531,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH2]C=CCC(O)([CH]C)C(=C)O(20177)',
    structure = SMILES('[CH2]C=CCC(O)([CH]C)C(=C)O'),
    E0 = (-80.8901,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,379.414,746.693,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154537,'amu*angstrom^2'), symmetry=1, barrier=(3.55312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154537,'amu*angstrom^2'), symmetry=1, barrier=(3.55312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154537,'amu*angstrom^2'), symmetry=1, barrier=(3.55312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154537,'amu*angstrom^2'), symmetry=1, barrier=(3.55312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154537,'amu*angstrom^2'), symmetry=1, barrier=(3.55312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154537,'amu*angstrom^2'), symmetry=1, barrier=(3.55312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154537,'amu*angstrom^2'), symmetry=1, barrier=(3.55312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154537,'amu*angstrom^2'), symmetry=1, barrier=(3.55312,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.25836,0.119608,-8.11805e-05,6.45143e-09,1.00411e-11,-9487.27,48.1827], Tmin=(100,'K'), Tmax=(965.735,'K')), NASAPolynomial(coeffs=[28.3598,0.0336602,-1.11661e-05,1.94135e-09,-1.35617e-13,-17307,-108.332], Tmin=(965.735,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-80.8901,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCJCO) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=C[CH]C(O)(CC)C(=C)O(20178)',
    structure = SMILES('[CH2]C=C[CH]C(O)(CC)C(=C)O'),
    E0 = (-210.649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.17433,0.110367,-3.88391e-05,-4.56817e-08,3.0227e-11,-25089.8,46.1835], Tmin=(100,'K'), Tmax=(951.665,'K')), NASAPolynomial(coeffs=[31.792,0.0294413,-8.75776e-06,1.52726e-09,-1.1215e-13,-34355,-130.715], Tmin=(951.665,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-210.649,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CCJC(O)C=C) + radical(Allyl_P)"""),
)

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
    label = '[CH2]C=CCC(O)(C[CH2])C(=C)O(20180)',
    structure = SMILES('[CH2]C=CCC(O)(C[CH2])C(=C)O'),
    E0 = (-75.5457,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,180,341.064,786.592,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153887,'amu*angstrom^2'), symmetry=1, barrier=(3.53816,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153887,'amu*angstrom^2'), symmetry=1, barrier=(3.53816,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153887,'amu*angstrom^2'), symmetry=1, barrier=(3.53816,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153887,'amu*angstrom^2'), symmetry=1, barrier=(3.53816,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153887,'amu*angstrom^2'), symmetry=1, barrier=(3.53816,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153887,'amu*angstrom^2'), symmetry=1, barrier=(3.53816,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153887,'amu*angstrom^2'), symmetry=1, barrier=(3.53816,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153887,'amu*angstrom^2'), symmetry=1, barrier=(3.53816,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.43273,0.125754,-0.000100935,2.90284e-08,1.42862e-12,-8840.47,48.0986], Tmin=(100,'K'), Tmax=(985.23,'K')), NASAPolynomial(coeffs=[27.7255,0.0354316,-1.23212e-05,2.15612e-09,-1.48903e-13,-16341.9,-104.86], Tmin=(985.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-75.5457,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJ) + radical(Allyl_P)"""),
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
    label = '[CH2]C=CCC(O)(CC)C(=C)[O](20182)',
    structure = SMILES('[CH2]C=CCC(O)(CC)C(=C)[O]'),
    E0 = (-142.987,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.1714,0.121831,-0.000100661,3.717e-08,-3.49678e-12,-16963.2,46.8758], Tmin=(100,'K'), Tmax=(1033.61,'K')), NASAPolynomial(coeffs=[24.6457,0.0403297,-1.47164e-05,2.59009e-09,-1.76681e-13,-23696.9,-89.1445], Tmin=(1033.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-142.987,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C(O)C(O)(CC)CC=C[CH2](20183)',
    structure = SMILES('[CH]=C(O)C(O)(CC)CC=C[CH2]'),
    E0 = (-33.696,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,180,180,180,180,1566.46,1600,2920.19,3200],'cm^-1')),
        HinderedRotor(inertia=(0.150498,'amu*angstrom^2'), symmetry=1, barrier=(3.46025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150498,'amu*angstrom^2'), symmetry=1, barrier=(3.46025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150498,'amu*angstrom^2'), symmetry=1, barrier=(3.46025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150498,'amu*angstrom^2'), symmetry=1, barrier=(3.46025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150498,'amu*angstrom^2'), symmetry=1, barrier=(3.46025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150498,'amu*angstrom^2'), symmetry=1, barrier=(3.46025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150498,'amu*angstrom^2'), symmetry=1, barrier=(3.46025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150498,'amu*angstrom^2'), symmetry=1, barrier=(3.46025,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.54577,0.128569,-0.000107236,3.42915e-08,-1.11935e-13,-3803.31,47.1209], Tmin=(100,'K'), Tmax=(983.672,'K')), NASAPolynomial(coeffs=[28.1893,0.0348638,-1.20388e-05,2.096e-09,-1.44316e-13,-11363.1,-108.343], Tmin=(983.672,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-33.696,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Allyl_P)"""),
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
    label = '[CH2][C]=CCC(O)(CC)C(=C)O(20185)',
    structure = SMILES('[CH2][C]=CCC(O)(CC)C(=C)O'),
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
    label = '[CH]=C[CH2](2461)',
    structure = SMILES('[CH]C=C'),
    E0 = (376.654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,229.711,230.18,230.787],'cm^-1')),
        HinderedRotor(inertia=(1.33306,'amu*angstrom^2'), symmetry=1, barrier=(50.5153,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31912,0.00817959,3.34736e-05,-4.36194e-08,1.58213e-11,45331.5,10.6389], Tmin=(100,'K'), Tmax=(983.754,'K')), NASAPolynomial(coeffs=[5.36755,0.0170743,-6.35108e-06,1.1662e-09,-8.2762e-14,44095,-3.44606], Tmin=(983.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(376.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(AllylJ2_triplet)"""),
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
    label = 'C=C(O)C([O])(CC)CC1[CH]C1(20190)',
    structure = SMILES('C=C(O)C([O])(CC)CC1[CH]C1'),
    E0 = (60.7422,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.62938,0.105512,-5.39593e-05,-9.46804e-09,1.22352e-11,7524.35,45.9489], Tmin=(100,'K'), Tmax=(1006.66,'K')), NASAPolynomial(coeffs=[24.1355,0.0414125,-1.54848e-05,2.82621e-09,-1.9918e-13,397.564,-88.1562], Tmin=(1006.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(60.7422,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(cyclopropane) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH2]C=CCC1(CC)OC[C]1O(20191)',
    structure = SMILES('[CH2]C=CCC1(CC)OC[C]1O'),
    E0 = (1.566,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.41242,0.123258,-0.000111058,5.40826e-08,-1.04338e-11,434.283,43.4774], Tmin=(100,'K'), Tmax=(1333.26,'K')), NASAPolynomial(coeffs=[23.6251,0.0401983,-1.20497e-05,1.79518e-09,-1.07952e-13,-6069.38,-87.9803], Tmin=(1333.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1.566,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Oxetane) + radical(Allyl_P) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C1[CH]CC(CC)(O1)C(=C)O(19976)',
    structure = SMILES('[CH2]C1[CH]CC(CC)(O1)C(=C)O'),
    E0 = (-44.4895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.59546,0.0988202,-1.18323e-05,-7.6257e-08,4.40877e-11,-5126.6,39.4161], Tmin=(100,'K'), Tmax=(899.626,'K')), NASAPolynomial(coeffs=[29.2549,0.0280601,-4.57848e-06,4.22883e-10,-2.37534e-14,-13364.7,-121.092], Tmin=(899.626,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-44.4895,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Tetrahydrofuran) + radical(CJC(C)OC) + radical(CCJCO)"""),
)

species(
    label = 'CCC1([O])CC=CCC[C]1O(20048)',
    structure = SMILES('CCC1([O])CC=CCC[C]1O'),
    E0 = (-27.9323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.28373,0.0644414,0.000128747,-2.53025e-07,1.13983e-10,-3120.84,40.8789], Tmin=(100,'K'), Tmax=(910.716,'K')), NASAPolynomial(coeffs=[44.2014,0.000599487,1.0005e-05,-2.20755e-09,1.4092e-13,-17042.9,-205.256], Tmin=(910.716,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-27.9323,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(594.485,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(C2CsJOH) + radical(CC(C)2OJ)"""),
)

species(
    label = 'C=C=CCC(O)(CC)C(=C)O(20192)',
    structure = SMILES('C=C=CCC(O)(CC)C(=C)O'),
    E0 = (-255.687,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.45625,0.127461,-0.000109125,3.98861e-08,-3.01247e-12,-30506.8,45.1696], Tmin=(100,'K'), Tmax=(1005.47,'K')), NASAPolynomial(coeffs=[27.0102,0.0366396,-1.30237e-05,2.28376e-09,-1.5667e-13,-37767,-103.786], Tmin=(1005.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-255.687,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
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
    label = '[CH2]C=CCC(C)([O])C(=C)O(20193)',
    structure = SMILES('[CH2]C=CCC(C)([O])C(=C)O'),
    E0 = (-27.9092,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,317.206,811.817,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152646,'amu*angstrom^2'), symmetry=1, barrier=(3.50964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152646,'amu*angstrom^2'), symmetry=1, barrier=(3.50964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152646,'amu*angstrom^2'), symmetry=1, barrier=(3.50964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152646,'amu*angstrom^2'), symmetry=1, barrier=(3.50964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152646,'amu*angstrom^2'), symmetry=1, barrier=(3.50964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152646,'amu*angstrom^2'), symmetry=1, barrier=(3.50964,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.38621,0.103249,-7.08102e-05,1.04966e-08,5.46064e-12,-3149.36,41.5263], Tmin=(100,'K'), Tmax=(1000.57,'K')), NASAPolynomial(coeffs=[23.8191,0.0332999,-1.21412e-05,2.1853e-09,-1.53074e-13,-9735.8,-87.7973], Tmin=(1000.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-27.9092,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(Allyl_P)"""),
)

species(
    label = 'C=C(O)C1(CC)CC=CCO1(19811)',
    structure = SMILES('C=C(O)C1(CC)CC=CCO1'),
    E0 = (-340.162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.25633,0.0898311,6.62812e-06,-8.53891e-08,4.3503e-11,-40699.1,35.9898], Tmin=(100,'K'), Tmax=(931.399,'K')), NASAPolynomial(coeffs=[27.0306,0.0339523,-9.03229e-06,1.44253e-09,-1.01807e-13,-48814,-113.744], Tmin=(931.399,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-340.162,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(590.328,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(36dihydro2hpyran)"""),
)

species(
    label = '[CH2]C=CCC([C]=C)(CC)OO(20194)',
    structure = SMILES('[CH2]C=CCC([C]=C)(CC)OO'),
    E0 = (243.113,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
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
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.8925,0.124392,-0.000111421,5.25674e-08,-1.0077e-11,29456.2,45.5543], Tmin=(100,'K'), Tmax=(1242.25,'K')), NASAPolynomial(coeffs=[21.4812,0.0491297,-2.05432e-05,3.79694e-09,-2.62026e-13,23648.9,-72.2808], Tmin=(1242.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.113,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=CCC([O])(CC)C(C)=O(20195)',
    structure = SMILES('[CH2]C=CCC([O])(CC)C(C)=O'),
    E0 = (-48.0855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.19462,0.117988,-0.000109295,5.72195e-08,-1.26425e-11,-5599.39,43.3929], Tmin=(100,'K'), Tmax=(1064.07,'K')), NASAPolynomial(coeffs=[14.5804,0.0586887,-2.57026e-05,4.84799e-09,-3.38163e-13,-8956.61,-33.6926], Tmin=(1064.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-48.0855,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(CC(C)(C=O)OJ)"""),
)

species(
    label = '[CH2]C(C=C)C[C](O)C(=O)CC(19808)',
    structure = SMILES('[CH2]C(C=C)CC(O)=C([O])CC'),
    E0 = (-112.953,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    molecularWeight = (154.206,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4873.63,'J/mol'), sigma=(8.12413,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=761.25 K, Pc=20.62 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.12093,0.119539,-8.80957e-05,1.73632e-08,5.69966e-12,-13351.1,47.1247], Tmin=(100,'K'), Tmax=(957.103,'K')), NASAPolynomial(coeffs=[26.0227,0.036442,-1.1968e-05,2.02345e-09,-1.37465e-13,-20319.6,-95.6787], Tmin=(957.103,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-112.953,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(C=C(C)OJ)"""),
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
    label = '[CH2]C=CC[C](CC)C(=C)O(20196)',
    structure = SMILES('[CH2]C=CCC(CC)=C([CH2])O'),
    E0 = (50.8553,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.67643,0.108941,-7.12715e-05,7.21923e-09,7.27242e-12,6334.96,41.5057], Tmin=(100,'K'), Tmax=(988.411,'K')), NASAPolynomial(coeffs=[24.3198,0.0372179,-1.3235e-05,2.34511e-09,-1.6286e-13,-439.493,-91.8814], Tmin=(988.411,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(50.8553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(C=C(O)CJ)"""),
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
    label = '[CH]=CCC([O])(CC)C(=C)O(20197)',
    structure = SMILES('[CH]=CCC([O])(CC)C(=C)O'),
    E0 = (79.9331,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,467.89,671.497,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156629,'amu*angstrom^2'), symmetry=1, barrier=(3.60122,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156629,'amu*angstrom^2'), symmetry=1, barrier=(3.60122,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156629,'amu*angstrom^2'), symmetry=1, barrier=(3.60122,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156629,'amu*angstrom^2'), symmetry=1, barrier=(3.60122,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156629,'amu*angstrom^2'), symmetry=1, barrier=(3.60122,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156629,'amu*angstrom^2'), symmetry=1, barrier=(3.60122,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.84314,0.112844,-0.000103852,4.8651e-08,-8.95984e-12,9837.29,43.7376], Tmin=(100,'K'), Tmax=(1323.42,'K')), NASAPolynomial(coeffs=[25.2596,0.0309267,-1.10042e-05,1.87929e-09,-1.24445e-13,2663.63,-94.6124], Tmin=(1323.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(79.9331,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH2]C1(O)C(C=C)CC1([O])CC(19911)',
    structure = SMILES('[CH2]C1(O)C(C=C)CC1([O])CC'),
    E0 = (74.8061,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.09621,0.117031,-7.99055e-05,1.41415e-08,4.52408e-12,9231.45,41.2692], Tmin=(100,'K'), Tmax=(1023.66,'K')), NASAPolynomial(coeffs=[25.6225,0.0408599,-1.53864e-05,2.79477e-09,-1.95453e-13,1872.53,-101.331], Tmin=(1023.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.8061,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(CJC(C)2O) + radical(CC(C)2OJ)"""),
)

species(
    label = 'C=CC=CC([O])(CC)C(=C)O(20198)',
    structure = SMILES('C=CC=CC([O])(CC)C(=C)O'),
    E0 = (-97.817,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,180,980.761,1228.03,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.160275,'amu*angstrom^2'), symmetry=1, barrier=(3.68503,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160275,'amu*angstrom^2'), symmetry=1, barrier=(3.68503,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160275,'amu*angstrom^2'), symmetry=1, barrier=(3.68503,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160275,'amu*angstrom^2'), symmetry=1, barrier=(3.68503,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160275,'amu*angstrom^2'), symmetry=1, barrier=(3.68503,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160275,'amu*angstrom^2'), symmetry=1, barrier=(3.68503,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (153.198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.91272,0.112295,-7.44718e-05,9.79102e-09,5.77459e-12,-11536.3,43.2013], Tmin=(100,'K'), Tmax=(1030.71,'K')), NASAPolynomial(coeffs=[26.0619,0.0374035,-1.44869e-05,2.68897e-09,-1.90784e-13,-19091.7,-101.283], Tmin=(1030.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-97.817,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(557.07,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=CC[CH]C([O])(CC)C(=C)O(20167)',
    structure = SMILES('C=CC[CH]C([O])(CC)C(=C)O'),
    E0 = (8.95876,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
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
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.02309,0.11672,-8.56518e-05,2.23318e-08,1.29885e-12,1308,48.7121], Tmin=(100,'K'), Tmax=(1035.48,'K')), NASAPolynomial(coeffs=[24.8559,0.0405006,-1.52385e-05,2.75023e-09,-1.90965e-13,-5738.86,-89.0487], Tmin=(1035.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(8.95876,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCJCO)"""),
)

species(
    label = 'C=[C]CCC([O])(CC)C(=C)O(20199)',
    structure = SMILES('C=[C]CCC([O])(CC)C(=C)O'),
    E0 = (46.8985,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,350,440,435,1725,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
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
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.14469,0.125242,-0.000114268,5.43259e-08,-1.03395e-11,5870.06,47.1279], Tmin=(100,'K'), Tmax=(1265.37,'K')), NASAPolynomial(coeffs=[23.834,0.0431206,-1.69203e-05,3.03825e-09,-2.0663e-13,-704.543,-84.3192], Tmin=(1265.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(46.8985,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CCCC([O])(CC)C(=C)O(20200)',
    structure = SMILES('[CH]=CCCC([O])(CC)C(=C)O'),
    E0 = (56.1529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.42131,0.126873,-0.000115278,5.38409e-08,-9.97062e-12,6996.83,48.0547], Tmin=(100,'K'), Tmax=(1307.95,'K')), NASAPolynomial(coeffs=[26.25,0.0391901,-1.47196e-05,2.58597e-09,-1.73813e-13,-503.302,-97.9651], Tmin=(1307.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(56.1529,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=CCCC([O])([CH]C)C(=C)O(20201)',
    structure = SMILES('C=CCCC([O])([CH]C)C(=C)O'),
    E0 = (8.95876,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
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
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.02309,0.11672,-8.56518e-05,2.23318e-08,1.29885e-12,1308,48.7121], Tmin=(100,'K'), Tmax=(1035.48,'K')), NASAPolynomial(coeffs=[24.8559,0.0405006,-1.52385e-05,2.75023e-09,-1.90965e-13,-5738.86,-89.0487], Tmin=(1035.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(8.95876,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCJCO)"""),
)

species(
    label = '[CH2]CC([O])(CCC=C)C(=C)O(20202)',
    structure = SMILES('[CH2]CC([O])(CCC=C)C(=C)O'),
    E0 = (14.3031,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,350,440,435,1725,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3100,440,815,1455,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
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
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.35527,0.124576,-0.000110617,5.04703e-08,-9.14012e-12,1961.77,49.2035], Tmin=(100,'K'), Tmax=(1335.25,'K')), NASAPolynomial(coeffs=[26.15,0.0391834,-1.46883e-05,2.57494e-09,-1.72679e-13,-5650.6,-96.5597], Tmin=(1335.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(14.3031,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(RCCJ)"""),
)

species(
    label = 'C=CCCC([O])(CC)C(=C)[O](20203)',
    structure = SMILES('C=CCCC([O])(CC)C(=C)[O]'),
    E0 = (-53.1385,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.77548,0.116991,-9.79559e-05,4.31103e-08,-7.68334e-12,-6174.93,46.8317], Tmin=(100,'K'), Tmax=(1335.24,'K')), NASAPolynomial(coeffs=[21.9503,0.0459175,-1.81141e-05,3.24752e-09,-2.19945e-13,-12511,-74.4917], Tmin=(1335.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-53.1385,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C(O)C([O])(CC)CCC=C(20204)',
    structure = SMILES('[CH]=C(O)C([O])(CC)CCC=C'),
    E0 = (56.1529,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,180,180,180,180,1023.87,1191.08,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.160858,'amu*angstrom^2'), symmetry=1, barrier=(3.69845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160858,'amu*angstrom^2'), symmetry=1, barrier=(3.69845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160858,'amu*angstrom^2'), symmetry=1, barrier=(3.69845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160858,'amu*angstrom^2'), symmetry=1, barrier=(3.69845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160858,'amu*angstrom^2'), symmetry=1, barrier=(3.69845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160858,'amu*angstrom^2'), symmetry=1, barrier=(3.69845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160858,'amu*angstrom^2'), symmetry=1, barrier=(3.69845,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.42131,0.126873,-0.000115278,5.38408e-08,-9.97061e-12,6996.83,48.0547], Tmin=(100,'K'), Tmax=(1307.96,'K')), NASAPolynomial(coeffs=[26.25,0.03919,-1.47195e-05,2.58597e-09,-1.73813e-13,-503.307,-97.9652], Tmin=(1307.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(56.1529,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C[CH]CC(O)(CC)C(=C)O(20205)',
    structure = SMILES('[CH]C=CCC(O)(CC)C(=C)O'),
    E0 = (-61.6068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.50547,0.127852,-0.000100038,3.1248e-08,-5.46084e-13,-7162.04,47.4017], Tmin=(100,'K'), Tmax=(1011.3,'K')), NASAPolynomial(coeffs=[25.8086,0.043491,-1.58924e-05,2.79445e-09,-1.90842e-13,-14301.7,-96.5011], Tmin=(1011.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-61.6068,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=CC1C[C](O)C([O])(CC)C1(19884)',
    structure = SMILES('C=CC1C[C](O)C([O])(CC)C1'),
    E0 = (-24.3652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.32994,0.0991218,-4.05522e-05,-1.90808e-08,1.46697e-11,-2722.67,41.2222], Tmin=(100,'K'), Tmax=(1007.81,'K')), NASAPolynomial(coeffs=[22.0771,0.0446803,-1.67672e-05,3.05281e-09,-2.14288e-13,-9393.83,-81.5754], Tmin=(1007.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-24.3652,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(590.328,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(CC(C)2OJ) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2][C](O)C1(CC)CC=CCO1(20004)',
    structure = SMILES('[CH2][C](O)C1(CC)CC=CCO1'),
    E0 = (-30.2471,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.15198,0.120146,-0.000103588,4.77198e-08,-8.81075e-12,-3403.67,39.6699], Tmin=(100,'K'), Tmax=(1308.11,'K')), NASAPolynomial(coeffs=[23.1825,0.0426764,-1.47541e-05,2.44599e-09,-1.58179e-13,-10031.7,-89.3587], Tmin=(1308.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-30.2471,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(36dihydro2hpyran) + radical(CJCO) + radical(C2CsJOH)"""),
)

species(
    label = 'C=CC=CC(O)(CC)C(=C)O(20206)',
    structure = SMILES('C=CC=CC(O)(CC)C(=C)O'),
    E0 = (-326.92,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.24506,0.118056,-7.60657e-05,2.61012e-09,1.0574e-11,-39077.3,43.1777], Tmin=(100,'K'), Tmax=(986.516,'K')), NASAPolynomial(coeffs=[28.6293,0.0344701,-1.22248e-05,2.21216e-09,-1.57337e-13,-47193.2,-115.614], Tmin=(986.516,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-326.92,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2]C(C=C)C([O])(CC)C(=C)O(19804)',
    structure = SMILES('[CH2]C(C=C)C([O])(CC)C(=C)O'),
    E0 = (6.16714,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,180,891.493,1315.34,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15743,'amu*angstrom^2'), symmetry=1, barrier=(3.61962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15743,'amu*angstrom^2'), symmetry=1, barrier=(3.61962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15743,'amu*angstrom^2'), symmetry=1, barrier=(3.61962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15743,'amu*angstrom^2'), symmetry=1, barrier=(3.61962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15743,'amu*angstrom^2'), symmetry=1, barrier=(3.61962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15743,'amu*angstrom^2'), symmetry=1, barrier=(3.61962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15743,'amu*angstrom^2'), symmetry=1, barrier=(3.61962,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4819.15,'J/mol'), sigma=(8.10355,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=752.74 K, Pc=20.55 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.10588,0.119884,-9.2049e-05,2.58179e-08,1.27997e-12,974.314,48.0474], Tmin=(100,'K'), Tmax=(991.066,'K')), NASAPolynomial(coeffs=[24.8373,0.0397338,-1.4018e-05,2.44045e-09,-1.66611e-13,-5770.46,-88.7809], Tmin=(991.066,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.16714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=CC1CC(CC)(O1)C(=C)O(19815)',
    structure = SMILES('C=CC1CC(CC)(O1)C(=C)O'),
    E0 = (-248.993,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (154.206,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.43733,0.0921676,9.3464e-06,-9.83822e-08,5.16117e-11,-29725.6,40.1818], Tmin=(100,'K'), Tmax=(908.003,'K')), NASAPolynomial(coeffs=[29.893,0.0273119,-4.3771e-06,4.33278e-10,-2.77891e-14,-38431.2,-124.553], Tmin=(908.003,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-248.993,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Oxetane)"""),
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
    E0 = (-51.6895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-4.00634,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (70.3535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (80.1065,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (197.007,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (3.09718,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-51.6895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (93.4743,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (156.042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (24.6719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (90.2874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (196.574,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (8.13425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (58.2272,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (50.1181,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (10.6126,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (366.711,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (1.3475,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (57.5069,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (65.2361,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (41.5172,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (176.184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (89.7708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (342.712,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (174.247,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (73.2029,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (5.66314,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-24.9119,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-26.7162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (391.953,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-44.1583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (350.224,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (152.49,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (95.3599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (293.86,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (461.496,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (74.8061,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (115.565,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (-51.6895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (134.061,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (247.095,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (216.505,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (50.5896,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (44.8463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (-22.6525,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (93.6781,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (77.0313,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (69.4091,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (9.19689,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (26.5576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (166.102,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (-43.4051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    products = ['butadiene13(2459)', 'C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    products = ['[CH2]C=CCC1(CC)OC1([CH2])O(20172)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.52e+08,'s^-1'), n=0.89, Ea=(47.6831,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H_secNd;radadd_intra] for rate rule [R4_S_D;doublebond_intra_2H_secNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    products = ['[CH2][C](O)C1(CC)CC(C=C)O1(19941)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    products = ['[CH2]C1(O)CC=CCC1([O])CC(19823)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.95e+10,'s^-1'), n=0.53, Ea=(131.796,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_SMSR;doublebond_intra_2H;radadd_intra_cs2H] for rate rule [R7_SMSS_D;doublebond_intra_2H_secNd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C=C=CCC([O])(CC)C(=C)O(20173)'],
    products = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.42e+08,'cm^3/(mol*s)'), n=1.64, Ea=(11.7989,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2713 used for Ca_Cds-HH;HJ
Exact match found for rate rule [Ca_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C2H5(29)', '[CH2]C=CCC(=O)C(=C)O(20174)'],
    products = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7.94e+10,'cm^3/(mol*s)'), n=0, Ea=(28.0328,'kJ/mol'), T0=(1,'K'), Tmin=(333,'K'), Tmax=(363,'K'), comment="""Estimated using template [CO_O;CsJ-CsHH] for rate rule [CO-CdCs_O;CsJ-CsHH]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2][CH]C=C(2458)', 'C=C(O)C(=O)CC(4626)'],
    products = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.000108643,'m^3/(mol*s)'), n=3.00879, Ea=(32.0644,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;CsJ-OneDeHH] for rate rule [CO-CdCs_O;CsJ-CdHH]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 27.6 to 32.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH2COH(99)', '[CH2]C=CCC(=O)CC(20175)'],
    products = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.16e+10,'cm^3/(mol*s)'), n=0, Ea=(48.1578,'kJ/mol'), T0=(1,'K'), Tmin=(413,'K'), Tmax=(563,'K'), comment="""Estimated using template [CO-CsCs_O;CJ] for rate rule [CO-CsCs_O;CdsJ-O2s]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    products = ['C=C(O)C([O])(CC)CC=[C]C(20176)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C=CCC(O)([CH]C)C(=C)O(20177)'],
    products = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.71,'s^-1'), n=3.021, Ea=(105.562,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 319 used for R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    products = ['[CH2]C=C[CH]C(O)(CC)C(=C)O(20178)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.07519e+07,'s^-1'), n=1.60667, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_H/Cd] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=C(O)C([O])(CC)C[C]=CC(20179)'],
    products = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C=CCC(O)(C[CH2])C(=C)O(20180)'],
    products = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C=[C]CC(O)(CC)C(=C)O(20181)'],
    products = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    products = ['[CH2]C=CCC(O)(CC)C(=C)[O](20182)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1070.11,'s^-1'), n=2.50856, Ea=(101.808,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;O_H_out] + [R4H_SSS;O_rad_out;XH_out] + [R4H_SS(Cd)S;Y_rad_out;XH_out] for rate rule [R4H_SS(Cd)S;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=C(O)C(O)(CC)CC=C[CH2](20183)'],
    products = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    products = ['C=C(O)C([O])([CH]C=CC)CC(20184)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2e+10,'s^-1'), n=0, Ea=(418.4,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4H_SDS;C_rad_out_single;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R4H_SDS;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][C]=CCC(O)(CC)C(=C)O(20185)'],
    products = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(380071,'s^-1'), n=1.62386, Ea=(44.2978,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;XH_out] for rate rule [R5H_DSSS;Cd_rad_out;O_H_out]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=C(O)C([O])([CH]C)CC=CC(20186)'],
    products = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(138.3,'s^-1'), n=3.21, Ea=(60.7935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H;C_rad_out_H/NonDeC;Cs_H_out] for rate rule [R6H_RSSMS;C_rad_out_H/NonDeC;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]CC([O])(CC=CC)C(=C)O(20187)'],
    products = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(64.2,'s^-1'), n=2.1, Ea=(63.1784,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 115 used for R7H;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R7H;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=C([O])C([O])(CC)CC=CC(20188)'],
    products = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(5.85e+08,'s^-1'), n=0, Ea=(106.901,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 143 used for R7H;O_rad_out;Cs_H_out_2H
Exact match found for rate rule [R7H;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C(O)C([O])(CC)CC=CC(20189)'],
    products = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(22.7193,'s^-1'), n=3.21897, Ea=(132.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R7H;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2][CH]C=C(2458)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.13324e+07,'m^3/(mol*s)'), n=0.074875, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_rad/H2/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -0.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=C[CH2](2461)', '[CH2]C([O])(CC)C(=C)O(6187)'],
    products = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.83556e+08,'m^3/(mol*s)'), n=-0.424942, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [C_pri_rad;Cd_rad] + [C_rad/H2/Cs;Y_rad] for rate rule [C_rad/H2/Cs;Cd_rad]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    products = ['C=C(O)C([O])(CC)CC1[CH]C1(20190)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra_cs2H] for rate rule [R3_D;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    products = ['[CH2]C=CCC1(CC)OC[C]1O(20191)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(3.21748e+08,'s^-1'), n=0.95, Ea=(124.892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R4_S_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    products = ['[CH2]C1[CH]CC(CC)(O1)C(=C)O(19976)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.66591e+07,'s^-1'), n=1.01661, Ea=(57.3526,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    products = ['CCC1([O])CC=CCC[C]1O(20048)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(185000,'s^-1'), n=1.07, Ea=(26.7776,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 50 used for R7_linear;doublebond_intra_secNd_2H;radadd_intra_cs2H
Exact match found for rate rule [R7_linear;doublebond_intra_secNd_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    products = ['C=C=CCC(O)(CC)C(=C)O(20192)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['CH2(S)(23)', '[CH2]C=CCC(C)([O])C(=C)O(20193)'],
    products = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    products = ['C=C(O)C1(CC)CC=CCO1(19811)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), Tmin=(550,'K'), Tmax=(650,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R6_SSSDS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C=CCC([C]=C)(CC)OO(20194)'],
    products = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(3.01978e+11,'s^-1'), n=0, Ea=(107.111,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2OOH_S;Y_rad_out] for rate rule [R2OOH_S;Cd_rad_out_double]
Euclidian distance = 2.0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    products = ['[CH2]C=CCC([O])(CC)C(C)=O(20195)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(204.179,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    products = ['[CH2]C(C=C)C[C](O)C(=O)CC(19808)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(6.21184e+10,'s^-1'), n=0.288169, Ea=(147.049,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_5_unsaturated_hexane] for rate rule [1_5_hexadiene]
Euclidian distance = 1.0
family: 6_membered_central_C-C_shift"""),
)

reaction(
    label = 'reaction35',
    reactants = ['O(4)', '[CH2]C=CC[C](CC)C(=C)O(20196)'],
    products = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/Cs2;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction36',
    reactants = ['CH2(19)', '[CH]=CCC([O])(CC)C(=C)O(20197)'],
    products = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    products = ['[CH2]C1(O)C(C=C)CC1([O])CC(19911)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(6.48e+06,'s^-1'), n=1.25, Ea=(126.496,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_cs] for rate rule [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_csHCd]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 126.3 to 126.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction38',
    reactants = ['H(3)', 'C=CC=CC([O])(CC)C(=C)O(20198)'],
    products = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.35e+08,'cm^3/(mol*s)'), n=1.64, Ea=(1.58992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2557 used for Cds-CsH_Cds-CdH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction39',
    reactants = ['butadiene13(2459)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(0.0534234,'m^3/(mol*s)'), n=2.459, Ea=(36.7983,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CdH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 32.5 to 36.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction40',
    reactants = ['C=CC[CH]C([O])(CC)C(=C)O(20167)'],
    products = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.682e+10,'s^-1'), n=0.35, Ea=(125.102,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 160 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C=[C]CCC([O])(CC)C(=C)O(20199)'],
    products = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.9054e+11,'s^-1'), n=0.853, Ea=(200.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    products = ['[CH]=CCCC([O])(CC)C(=C)O(20200)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(8.32e+10,'s^-1'), n=0.77, Ea=(268.194,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 195 used for R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleH
Exact match found for rate rule [R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['C=CCCC([O])([CH]C)C(=C)O(20201)'],
    products = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(0.502,'s^-1'), n=3.86, Ea=(41.6308,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 332 used for R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/Cd
Exact match found for rate rule [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2]CC([O])(CCC=C)C(=C)O(20202)'],
    products = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(15400,'s^-1'), n=1.87, Ea=(30.5432,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 87 used for R5H_CCC;C_rad_out_2H;Cs_H_out_H/Cd
Exact match found for rate rule [R5H_CCC;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    products = ['C=CCCC([O])(CC)C(=C)[O](20203)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(0.0508,'s^-1'), n=3.24, Ea=(29.037,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_1H;XH_out] for rate rule [R5H_CCC(Cd);C_rad_out_H/Cd;O_H_out]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH]=C(O)C([O])(CC)CCC=C(20204)'],
    products = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(14044.1,'s^-1'), n=1.88125, Ea=(37.5252,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cs_H_out_H/Cd]
Euclidian distance = 3.60555127546
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    products = ['[CH]=C[CH]CC(O)(CC)C(=C)O(20205)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(3.427,'s^-1'), n=3.311, Ea=(128.721,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;O_rad_out;Cd_H_out_singleH] for rate rule [R6HJ_3;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    products = ['C=CC1C[C](O)C([O])(CC)C1(19884)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(3.65565e+12,'s^-1'), n=-0.288384, Ea=(121.099,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5_SS_D;doublebond_intra;radadd_intra_csHCd] + [R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_cs] for rate rule [R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_csHCd]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction49',
    reactants = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    products = ['[CH2][C](O)C1(CC)CC=CCO1(20004)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(9.91671e+09,'s^-1'), n=0.30082, Ea=(60.8864,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra_pri_2H;radadd_intra] for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction50',
    reactants = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    products = ['C=CC=CC(O)(CC)C(=C)O(20206)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(C=C)C([O])(CC)C(=C)O(19804)'],
    products = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction52',
    reactants = ['C=C[CH]CC([O])(CC)C(=C)O(19803)'],
    products = ['C=CC1CC(CC)(O1)C(=C)O(19815)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;O_rad;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

network(
    label = '3593',
    isomers = [
        'C=C[CH]CC([O])(CC)C(=C)O(19803)',
    ],
    reactants = [
        ('butadiene13(2459)', 'C=C(O)C(=O)CC(4626)'),
        ('butadiene13(2459)', '[CH2]C(O)=C([O])CC(4557)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '3593',
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

