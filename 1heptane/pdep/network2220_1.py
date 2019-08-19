species(
    label = 'C=C(O)C([O])(CC)C([O])CC(11978)',
    structure = SMILES('C=C(O)C([O])(CC)C([O])CC'),
    E0 = (-237.055,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
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
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.47139,0.128732,-0.000118784,5.63731e-08,-1.06052e-11,-28266.7,47.2772], Tmin=(100,'K'), Tmax=(1288.49,'K')), NASAPolynomial(coeffs=[26.2726,0.0394988,-1.49038e-05,2.62568e-09,-1.76857e-13,-35674.1,-98.682], Tmin=(1288.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-237.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C2H5CHO(70)',
    structure = SMILES('CCC=O'),
    E0 = (-204.33,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.207559,'amu*angstrom^2'), symmetry=1, barrier=(4.77219,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208362,'amu*angstrom^2'), symmetry=1, barrier=(4.79065,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0791,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3133.67,'J/mol'), sigma=(5.35118,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=489.47 K, Pc=46.4 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.90578,0.0240644,-7.06356e-06,-9.81837e-10,5.55825e-13,-24535.9,13.5806], Tmin=(100,'K'), Tmax=(1712.49,'K')), NASAPolynomial(coeffs=[7.69109,0.0189242,-7.84934e-06,1.38273e-09,-8.99057e-14,-27060.1,-14.6647], Tmin=(1712.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-204.33,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""propanal""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2]C1(O)OC1(CC)C([O])CC(12253)',
    structure = SMILES('[CH2]C1(O)OC1(CC)C([O])CC'),
    E0 = (-190.853,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.3995,0.149613,-0.000153031,7.73429e-08,-1.46313e-11,-22621.7,47.125], Tmin=(100,'K'), Tmax=(1500.14,'K')), NASAPolynomial(coeffs=[34.399,0.0225549,-2.38205e-06,-9.51214e-11,2.18365e-14,-31606.3,-146.937], Tmin=(1500.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-190.853,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(CC(C)OJ) + radical(CJC(O)2C)"""),
)

species(
    label = '[CH2]C1(O)OC(CC)C1([O])CC(12126)',
    structure = SMILES('[CH2]C1(O)OC(CC)C1([O])CC'),
    E0 = (-203.702,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.59877,0.145275,-0.000142525,6.89071e-08,-1.23829e-11,-24152.3,47.2731], Tmin=(100,'K'), Tmax=(1611.25,'K')), NASAPolynomial(coeffs=[33.6611,0.0221112,-1.62847e-06,-2.45362e-10,3.10865e-14,-32823.4,-144.207], Tmin=(1611.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-203.702,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CC(C)2OJ) + radical(CJC(C)OC)"""),
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
    label = 'C=C(O)C([O])(CC)C(=O)CC(12254)',
    structure = SMILES('C=C(O)C([O])(CC)C(=O)CC'),
    E0 = (-378.585,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,180,180,180,180,971.001,1254.74,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.166718,'amu*angstrom^2'), symmetry=1, barrier=(3.83317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166718,'amu*angstrom^2'), symmetry=1, barrier=(3.83317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166718,'amu*angstrom^2'), symmetry=1, barrier=(3.83317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166718,'amu*angstrom^2'), symmetry=1, barrier=(3.83317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166718,'amu*angstrom^2'), symmetry=1, barrier=(3.83317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166718,'amu*angstrom^2'), symmetry=1, barrier=(3.83317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166718,'amu*angstrom^2'), symmetry=1, barrier=(3.83317,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (157.187,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.16852,0.1436,-0.000189391,1.40541e-07,-4.23338e-11,-45318.5,41.8545], Tmin=(100,'K'), Tmax=(808.996,'K')), NASAPolynomial(coeffs=[15.8314,0.0546025,-2.43798e-05,4.56275e-09,-3.13848e-13,-48230.9,-41.1696], Tmin=(808.996,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-378.585,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)(C=O)OJ)"""),
)

species(
    label = 'CC[CH][O](563)',
    structure = SMILES('CC[CH][O]'),
    E0 = (133.127,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,298.357,1774.23],'cm^-1')),
        HinderedRotor(inertia=(0.129074,'amu*angstrom^2'), symmetry=1, barrier=(8.14273,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00364816,'amu*angstrom^2'), symmetry=1, barrier=(8.14268,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0791,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.1585,0.0245341,-8.42945e-06,1.83944e-10,2.32791e-13,16036.2,14.3859], Tmin=(100,'K'), Tmax=(2077.96,'K')), NASAPolynomial(coeffs=[11.8474,0.0146996,-6.30487e-06,1.09829e-09,-6.9226e-14,10937.4,-37.4679], Tmin=(2077.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(133.127,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCOJ) + radical(CCsJOH)"""),
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
    label = 'C=C(O)C(=O)C([O])CC(12255)',
    structure = SMILES('C=C(O)C(=O)C([O])CC'),
    E0 = (-308.396,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.587807,0.0911445,-8.4399e-05,3.93407e-08,-7.23762e-12,-36918.2,35.7185], Tmin=(100,'K'), Tmax=(1317.09,'K')), NASAPolynomial(coeffs=[20.7988,0.0261923,-1.04254e-05,1.89705e-09,-1.30224e-13,-42551.7,-73.3494], Tmin=(1317.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-308.396,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsHH) + radical(C=OCOJ)"""),
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
    label = 'CCC(=O)C([O])CC(11179)',
    structure = SMILES('CCC(=O)C([O])CC'),
    E0 = (-233.54,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.15,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.647945,0.074755,-5.41576e-05,2.1084e-08,-3.50558e-12,-27968.9,32.1785], Tmin=(100,'K'), Tmax=(1349.41,'K')), NASAPolynomial(coeffs=[11.4824,0.0426388,-1.84574e-05,3.44652e-09,-2.37962e-13,-30893,-23.3384], Tmin=(1349.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-233.54,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(C=OCOJ)"""),
)

species(
    label = 'C=C(O)C([O])(C=O)CC(6168)',
    structure = SMILES('C=C(O)C([O])(C=O)CC'),
    E0 = (-298.298,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,180,180,180,180,1600,1818.45,2694.67,3200],'cm^-1')),
        HinderedRotor(inertia=(0.159255,'amu*angstrom^2'), symmetry=1, barrier=(3.66158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159255,'amu*angstrom^2'), symmetry=1, barrier=(3.66158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159255,'amu*angstrom^2'), symmetry=1, barrier=(3.66158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159255,'amu*angstrom^2'), symmetry=1, barrier=(3.66158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159255,'amu*angstrom^2'), symmetry=1, barrier=(3.66158,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.608054,0.10466,-0.000128562,8.50037e-08,-2.25179e-11,-35713.8,33.3289], Tmin=(100,'K'), Tmax=(920.683,'K')), NASAPolynomial(coeffs=[15.433,0.0349675,-1.50146e-05,2.78299e-09,-1.91503e-13,-38667.5,-42.7339], Tmin=(920.683,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-298.298,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=CC(C)(C=O)OJ)"""),
)

species(
    label = 'C=C(O)C([O])(CC)[C](O)CC(12256)',
    structure = SMILES('C=C(O)C([O])(CC)[C](O)CC'),
    E0 = (-290.788,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.29158,0.135652,-0.000141756,7.8325e-08,-1.73259e-11,-34744.7,46.2468], Tmin=(100,'K'), Tmax=(1096.03,'K')), NASAPolynomial(coeffs=[22.3727,0.0456365,-1.85595e-05,3.38806e-09,-2.3267e-13,-40151.2,-75.0052], Tmin=(1096.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-290.788,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(C2CsJOH)"""),
)

species(
    label = 'C=C(O)C(O)(CC)[C]([O])CC(12257)',
    structure = SMILES('C=C(O)C(O)(CC)[C]([O])CC'),
    E0 = (-289.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.76878,0.141059,-0.000150174,8.25304e-08,-1.78558e-11,-34571.9,46.5493], Tmin=(100,'K'), Tmax=(1130.8,'K')), NASAPolynomial(coeffs=[26.1478,0.0387716,-1.44891e-05,2.53687e-09,-1.70561e-13,-41111.7,-96.5108], Tmin=(1130.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-289.53,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C2CsJOH) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=C(O)C(O)([CH]C)C([O])CC(12258)',
    structure = SMILES('C=C(O)C(O)([CH]C)C([O])CC'),
    E0 = (-266.256,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,180,1079.9,1130.97,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.161347,'amu*angstrom^2'), symmetry=1, barrier=(3.70968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161347,'amu*angstrom^2'), symmetry=1, barrier=(3.70968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161347,'amu*angstrom^2'), symmetry=1, barrier=(3.70968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161347,'amu*angstrom^2'), symmetry=1, barrier=(3.70968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161347,'amu*angstrom^2'), symmetry=1, barrier=(3.70968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161347,'amu*angstrom^2'), symmetry=1, barrier=(3.70968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161347,'amu*angstrom^2'), symmetry=1, barrier=(3.70968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161347,'amu*angstrom^2'), symmetry=1, barrier=(3.70968,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.50845,0.128179,-0.000108391,3.61678e-08,-8.52239e-13,-31775.5,48.6854], Tmin=(100,'K'), Tmax=(983.237,'K')), NASAPolynomial(coeffs=[27.9116,0.0346884,-1.19336e-05,2.0707e-09,-1.42191e-13,-39220.4,-104.999], Tmin=(983.237,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-266.256,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CCJCO)"""),
)

species(
    label = 'C=C(O)C([O])(CC)C(O)[CH]C(12259)',
    structure = SMILES('C=C(O)C([O])(CC)C(O)[CH]C'),
    E0 = (-267.514,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,180,1039.39,1172.77,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.160509,'amu*angstrom^2'), symmetry=1, barrier=(3.69042,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160509,'amu*angstrom^2'), symmetry=1, barrier=(3.69042,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160509,'amu*angstrom^2'), symmetry=1, barrier=(3.69042,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160509,'amu*angstrom^2'), symmetry=1, barrier=(3.69042,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160509,'amu*angstrom^2'), symmetry=1, barrier=(3.69042,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160509,'amu*angstrom^2'), symmetry=1, barrier=(3.69042,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160509,'amu*angstrom^2'), symmetry=1, barrier=(3.69042,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160509,'amu*angstrom^2'), symmetry=1, barrier=(3.69042,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.55359,0.128847,-0.000120717,5.79733e-08,-1.09751e-11,-31925.6,50.2612], Tmin=(100,'K'), Tmax=(1286.32,'K')), NASAPolynomial(coeffs=[27.0775,0.036705,-1.32696e-05,2.28646e-09,-1.52242e-13,-39548.7,-100.152], Tmin=(1286.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-267.514,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCJCO) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C(O)C([O])([CH]C)C(O)CC(12260)',
    structure = SMILES('C=C(O)C([O])([CH]C)C(O)CC'),
    E0 = (-267.514,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.55359,0.128847,-0.000120717,5.79733e-08,-1.09751e-11,-31925.6,50.2612], Tmin=(100,'K'), Tmax=(1286.32,'K')), NASAPolynomial(coeffs=[27.0775,0.036705,-1.32696e-05,2.28646e-09,-1.52242e-13,-39548.7,-100.152], Tmin=(1286.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-267.514,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCJCO)"""),
)

species(
    label = 'C=C(O)C(O)(CC)C([O])[CH]C(12261)',
    structure = SMILES('C=C(O)C(O)(CC)C([O])[CH]C'),
    E0 = (-266.256,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.50845,0.128179,-0.000108391,3.61678e-08,-8.52239e-13,-31775.5,48.6854], Tmin=(100,'K'), Tmax=(983.237,'K')), NASAPolynomial(coeffs=[27.9116,0.0346884,-1.19336e-05,2.0707e-09,-1.42191e-13,-39220.4,-104.999], Tmin=(983.237,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-266.256,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CCJCO)"""),
)

species(
    label = '[CH2]CC(O)(C(=C)O)C([O])CC(12262)',
    structure = SMILES('[CH2]CC(O)(C(=C)O)C([O])CC'),
    E0 = (-260.911,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,180,180,180,180,1037.46,1175.28,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.16049,'amu*angstrom^2'), symmetry=1, barrier=(3.68998,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16049,'amu*angstrom^2'), symmetry=1, barrier=(3.68998,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16049,'amu*angstrom^2'), symmetry=1, barrier=(3.68998,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16049,'amu*angstrom^2'), symmetry=1, barrier=(3.68998,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16049,'amu*angstrom^2'), symmetry=1, barrier=(3.68998,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16049,'amu*angstrom^2'), symmetry=1, barrier=(3.68998,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16049,'amu*angstrom^2'), symmetry=1, barrier=(3.68998,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16049,'amu*angstrom^2'), symmetry=1, barrier=(3.68998,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.92486,0.137104,-0.000137451,7.01206e-08,-1.39818e-11,-31118.2,49.474], Tmin=(100,'K'), Tmax=(1230.06,'K')), NASAPolynomial(coeffs=[28.7879,0.0339777,-1.16923e-05,1.96157e-09,-1.2898e-13,-38919.8,-110.088], Tmin=(1230.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-260.911,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]CC(O)C([O])(CC)C(=C)O(12263)',
    structure = SMILES('[CH2]CC(O)C([O])(CC)C(=C)O'),
    E0 = (-262.169,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,180,180,180,180,1015.12,1198.03,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.160021,'amu*angstrom^2'), symmetry=1, barrier=(3.67921,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160021,'amu*angstrom^2'), symmetry=1, barrier=(3.67921,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160021,'amu*angstrom^2'), symmetry=1, barrier=(3.67921,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160021,'amu*angstrom^2'), symmetry=1, barrier=(3.67921,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160021,'amu*angstrom^2'), symmetry=1, barrier=(3.67921,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160021,'amu*angstrom^2'), symmetry=1, barrier=(3.67921,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160021,'amu*angstrom^2'), symmetry=1, barrier=(3.67921,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160021,'amu*angstrom^2'), symmetry=1, barrier=(3.67921,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.39597,0.13116,-0.000127451,6.4214e-08,-1.28504e-11,-31293.4,48.9808], Tmin=(100,'K'), Tmax=(1212.58,'K')), NASAPolynomial(coeffs=[24.9718,0.0408792,-1.57699e-05,2.81208e-09,-1.90883e-13,-37930.4,-88.3277], Tmin=(1212.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-262.169,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(RCCJ)"""),
)

species(
    label = 'C=C([O])C(O)(CC)C([O])CC(12264)',
    structure = SMILES('C=C([O])C(O)(CC)C([O])CC'),
    E0 = (-328.353,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.36826,0.129746,-0.000125402,6.33572e-08,-1.27158e-11,-39253.7,47.189], Tmin=(100,'K'), Tmax=(1210.46,'K')), NASAPolynomial(coeffs=[24.5509,0.0407918,-1.51711e-05,2.64782e-09,-1.77451e-13,-45770.7,-87.8224], Tmin=(1210.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-328.353,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=C(O)C(O)(CC)C([O])CC(12265)',
    structure = SMILES('[CH]=C(O)C(O)(CC)C([O])CC'),
    E0 = (-219.061,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,535.133,604.051,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157707,'amu*angstrom^2'), symmetry=1, barrier=(3.62598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157707,'amu*angstrom^2'), symmetry=1, barrier=(3.62598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157707,'amu*angstrom^2'), symmetry=1, barrier=(3.62598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157707,'amu*angstrom^2'), symmetry=1, barrier=(3.62598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157707,'amu*angstrom^2'), symmetry=1, barrier=(3.62598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157707,'amu*angstrom^2'), symmetry=1, barrier=(3.62598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157707,'amu*angstrom^2'), symmetry=1, barrier=(3.62598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157707,'amu*angstrom^2'), symmetry=1, barrier=(3.62598,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.00667,0.139569,-0.000142615,7.40361e-08,-1.5003e-11,-26082.4,48.383], Tmin=(100,'K'), Tmax=(1211.63,'K')), NASAPolynomial(coeffs=[29.0527,0.0337303,-1.15874e-05,1.94218e-09,-1.27696e-13,-33851.2,-112.439], Tmin=(1211.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-219.061,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cds_P) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]CC([O])(C(=C)O)C(O)CC(12266)',
    structure = SMILES('[CH2]CC([O])(C(=C)O)C(O)CC'),
    E0 = (-262.169,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.39598,0.13116,-0.000127451,6.42141e-08,-1.28504e-11,-31293.4,48.9808], Tmin=(100,'K'), Tmax=(1212.56,'K')), NASAPolynomial(coeffs=[24.9718,0.0408792,-1.577e-05,2.81209e-09,-1.90884e-13,-37930.4,-88.3275], Tmin=(1212.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-262.169,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJ) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C([O])C([O])(CC)C(O)CC(12267)',
    structure = SMILES('C=C([O])C([O])(CC)C(O)CC'),
    E0 = (-329.611,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.85099,0.123917,-0.000115721,5.77711e-08,-1.16898e-11,-39428.4,46.7391], Tmin=(100,'K'), Tmax=(1184.85,'K')), NASAPolynomial(coeffs=[20.7157,0.0477335,-1.92752e-05,3.50514e-09,-2.39955e-13,-44776,-65.96], Tmin=(1184.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-329.611,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C(O)C([O])(CC)C(O)CC(12268)',
    structure = SMILES('[CH]=C(O)C([O])(CC)C(O)CC'),
    E0 = (-220.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.48457,0.133693,-0.000132809,6.83292e-08,-1.39389e-11,-26257.3,47.9151], Tmin=(100,'K'), Tmax=(1191.64,'K')), NASAPolynomial(coeffs=[25.2382,0.0406342,-1.56686e-05,2.79387e-09,-1.89719e-13,-32864.3,-90.6915], Tmin=(1191.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-220.32,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH2]CC([O])C(O)(CC)C(=C)O(12269)',
    structure = SMILES('[CH2]CC([O])C(O)(CC)C(=C)O'),
    E0 = (-260.911,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.92486,0.137104,-0.000137451,7.01205e-08,-1.39818e-11,-31118.2,49.474], Tmin=(100,'K'), Tmax=(1230.07,'K')), NASAPolynomial(coeffs=[28.7879,0.0339777,-1.16923e-05,1.96157e-09,-1.2898e-13,-38919.8,-110.088], Tmin=(1230.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-260.911,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(RCCJ)"""),
)

species(
    label = 'CCC([O])C1(CC)OC[C]1O(12270)',
    structure = SMILES('CCC([O])C1(CC)OC[C]1O'),
    E0 = (-183.799,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.96011,0.123573,-0.000109533,4.6753e-08,-5.61574e-12,-21884.5,41.4588], Tmin=(100,'K'), Tmax=(917.817,'K')), NASAPolynomial(coeffs=[20.6525,0.0454432,-1.52171e-05,2.48581e-09,-1.60728e-13,-26895.4,-70.3799], Tmin=(917.817,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-183.799,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(C2CsJOH) + radical(CC(C)OJ)"""),
)

species(
    label = 'CCC1OC[C](O)C1([O])CC(12076)',
    structure = SMILES('CCC1OC[C](O)C1([O])CC'),
    E0 = (-263.707,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.59136,0.112238,-7.97434e-05,1.91945e-08,3.05605e-12,-31505.5,37.4937], Tmin=(100,'K'), Tmax=(942.081,'K')), NASAPolynomial(coeffs=[19.9076,0.0464744,-1.56656e-05,2.60309e-09,-1.71387e-13,-36688.7,-70.9537], Tmin=(942.081,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-263.707,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(586.17,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(CC(C)2OJ) + radical(C2CsJOH)"""),
)

species(
    label = 'C=C(O)C(O)(CC)C(=O)CC(12271)',
    structure = SMILES('C=C(O)C(O)(CC)C(=O)CC'),
    E0 = (-622.872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.35737,0.140538,-0.000157212,9.37037e-08,-2.23684e-11,-74685.7,41.5369], Tmin=(100,'K'), Tmax=(1018.29,'K')), NASAPolynomial(coeffs=[21.4273,0.047109,-1.95868e-05,3.60274e-09,-2.47998e-13,-79529.7,-73.6418], Tmin=(1018.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-622.872,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH)"""),
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
    label = 'C=C(O)C(C)([O])C([O])CC(12272)',
    structure = SMILES('C=C(O)C(C)([O])C([O])CC'),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.88961,0.114667,-0.000107262,5.109e-08,-9.56476e-12,-25426.5,42.9466], Tmin=(100,'K'), Tmax=(1302.58,'K')), NASAPolynomial(coeffs=[25.2787,0.0312377,-1.11881e-05,1.91867e-09,-1.27447e-13,-32504.2,-95.3067], Tmin=(1302.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-213.275,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C(O)C1(CC)OOC1CC(11988)',
    structure = SMILES('C=C(O)C1(CC)OOC1CC'),
    E0 = (-295.831,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.06893,0.117573,-8.48605e-05,2.04498e-08,2.1996e-12,-35347.9,38.6705], Tmin=(100,'K'), Tmax=(1025.62,'K')), NASAPolynomial(coeffs=[24.9573,0.0410987,-1.53265e-05,2.75463e-09,-1.91025e-13,-42413.2,-99.817], Tmin=(1025.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-295.831,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(12dioxetane)"""),
)

species(
    label = 'C=[C]C(CC)(OO)C([O])CC(12273)',
    structure = SMILES('C=[C]C(CC)(OO)C([O])CC'),
    E0 = (57.7478,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1380,1390,370,380,2900,435,180,180,180,180,836.201,1375.13,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156185,'amu*angstrom^2'), symmetry=1, barrier=(3.591,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156185,'amu*angstrom^2'), symmetry=1, barrier=(3.591,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156185,'amu*angstrom^2'), symmetry=1, barrier=(3.591,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156185,'amu*angstrom^2'), symmetry=1, barrier=(3.591,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156185,'amu*angstrom^2'), symmetry=1, barrier=(3.591,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156185,'amu*angstrom^2'), symmetry=1, barrier=(3.591,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156185,'amu*angstrom^2'), symmetry=1, barrier=(3.591,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156185,'amu*angstrom^2'), symmetry=1, barrier=(3.591,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.86925,0.129672,-0.000126913,6.69895e-08,-1.44713e-11,7156.37,45.0825], Tmin=(100,'K'), Tmax=(1105.75,'K')), NASAPolynomial(coeffs=[19.315,0.0530384,-2.29557e-05,4.31218e-09,-3.00441e-13,2471.5,-59.2488], Tmin=(1105.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(57.7478,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C([O])(CC)C(CC)OO(12274)',
    structure = SMILES('C=[C]C([O])(CC)C(CC)OO'),
    E0 = (56.4896,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1380,1390,370,380,2900,435,180,180,180,180,789.681,1422.7,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155216,'amu*angstrom^2'), symmetry=1, barrier=(3.56871,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155216,'amu*angstrom^2'), symmetry=1, barrier=(3.56871,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155216,'amu*angstrom^2'), symmetry=1, barrier=(3.56871,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155216,'amu*angstrom^2'), symmetry=1, barrier=(3.56871,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155216,'amu*angstrom^2'), symmetry=1, barrier=(3.56871,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155216,'amu*angstrom^2'), symmetry=1, barrier=(3.56871,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155216,'amu*angstrom^2'), symmetry=1, barrier=(3.56871,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155216,'amu*angstrom^2'), symmetry=1, barrier=(3.56871,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.44115,0.124778,-0.000120013,6.44299e-08,-1.45323e-11,6985.85,44.9609], Tmin=(100,'K'), Tmax=(1045.41,'K')), NASAPolynomial(coeffs=[15.4805,0.0600311,-2.71114e-05,5.18556e-09,-3.64522e-13,3447.85,-37.4277], Tmin=(1045.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(56.4896,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(573.699,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'CCC([O])C([O])(CC)C(C)=O(12275)',
    structure = SMILES('CCC([O])C([O])(CC)C(C)=O'),
    E0 = (-233.451,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.43599,0.126294,-0.000134882,8.40345e-08,-2.20308e-11,-27887.7,43.8759], Tmin=(100,'K'), Tmax=(909.918,'K')), NASAPolynomial(coeffs=[13.4181,0.0609939,-2.7234e-05,5.16275e-09,-3.60434e-13,-30590.8,-26.384], Tmin=(909.918,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-233.451,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(CC(C)OJ) + radical(CC(C)(C=O)OJ)"""),
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
    label = 'C=C(O)[C](CC)C([O])CC(12276)',
    structure = SMILES('[CH2]C(O)=C(CC)C([O])CC'),
    E0 = (-136.116,'kJ/mol'),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.8825,0.114467,-8.6794e-05,2.29118e-08,1.72584e-12,-16146.1,42.7328], Tmin=(100,'K'), Tmax=(1006.43,'K')), NASAPolynomial(coeffs=[24.886,0.0365994,-1.32506e-05,2.356e-09,-1.63105e-13,-22978.7,-93.7578], Tmin=(1006.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-136.116,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(C=C(O)CJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=C(O)C([O])([CH]CC)CC(6800)',
    structure = SMILES('C=C(O)C([O])([CH]CC)CC'),
    E0 = (-95.0432,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,180,1064.56,1147.69,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.60905,0.108824,-7.81789e-05,1.95543e-08,1.52878e-12,-11216.5,44.5071], Tmin=(100,'K'), Tmax=(1033.85,'K')), NASAPolynomial(coeffs=[22.9346,0.0394958,-1.47805e-05,2.65393e-09,-1.83543e-13,-17661.2,-81.3439], Tmin=(1033.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-95.0432,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCJCO) + radical(C=CC(C)2OJ)"""),
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
    E0 = (-237.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-189.372,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-144.653,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-95.3513,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-193.124,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-172.49,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-82.1124,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-237.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-162.391,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-71.3846,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-161.774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-160.694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-161.952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-183.165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-183.165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-177.231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-178.489,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-135.247,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-174.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-205.215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-199.348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-125.61,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-205.215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-51.8162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-112.163,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-175.271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-173.655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (206.587,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (206.587,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-228.771,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (164.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (169.039,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (-32.8759,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (106.889,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (147.962,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    products = ['C2H5CHO(70)', 'C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    products = ['[CH2]C1(O)OC1(CC)C([O])CC(12253)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.52e+08,'s^-1'), n=0.89, Ea=(47.6831,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H_secNd;radadd_intra] for rate rule [R4_S_D;doublebond_intra_2H_secNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    products = ['[CH2]C1(O)OC(CC)C1([O])CC(12126)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.20137e+08,'s^-1'), n=0.864, Ea=(92.4015,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5_SS_D;doublebond_intra_2H;radadd_intra_O] + [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C=C(O)C([O])(CC)C(=O)CC(12254)'],
    products = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.0366254,'m^3/(mol*s)'), n=1.743, Ea=(71.4418,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsCs_O;YJ] for rate rule [CO-CsCs_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CC[CH][O](563)', 'C=C(O)C(=O)CC(4626)'],
    products = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.97e+07,'cm^3/(mol*s)'), n=1.88, Ea=(32.2168,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO-CdCs_O;YJ] for rate rule [CO-CdCs_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C2H5(29)', 'C=C(O)C(=O)C([O])CC(12255)'],
    products = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7.94e+10,'cm^3/(mol*s)'), n=0, Ea=(28.0328,'kJ/mol'), T0=(1,'K'), Tmin=(333,'K'), Tmax=(363,'K'), comment="""Estimated using template [CO_O;CsJ-CsHH] for rate rule [CO-CdCs_O;CsJ-CsHH]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2COH(99)', 'CCC(=O)C([O])CC(11179)'],
    products = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.16e+10,'cm^3/(mol*s)'), n=0, Ea=(48.1578,'kJ/mol'), T0=(1,'K'), Tmin=(413,'K'), Tmax=(563,'K'), comment="""Estimated using template [CO-CsCs_O;CJ] for rate rule [CO-CsCs_O;CdsJ-O2s]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C2H5CHO(70)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.0201871,'m^3/(mol*s)'), n=2.2105, Ea=(152.218,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsH_O;YJ] for rate rule [CO-CsH_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 147.5 to 152.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['C2H5(29)', 'C=C(O)C([O])(C=O)CC(6168)'],
    products = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.94e+10,'cm^3/(mol*s)'), n=0, Ea=(28.0328,'kJ/mol'), T0=(1,'K'), Tmin=(333,'K'), Tmax=(363,'K'), comment="""Estimated using template [CO_O;CsJ-CsHH] for rate rule [CO-CsH_O;CsJ-CsHH]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    products = ['C=C(O)C([O])(CC)[C](O)CC(12256)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.56178e+08,'s^-1'), n=1.25272, Ea=(165.67,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_Cs2] for rate rule [R2H_S;O_rad_out;Cs_H_out_Cs2]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    products = ['C=C(O)C(O)(CC)[C]([O])CC(12257)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=C(O)C(O)([CH]C)C([O])CC(12258)'],
    products = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(5.71,'s^-1'), n=3.021, Ea=(105.562,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 319 used for R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C(O)C([O])(CC)C(O)[CH]C(12259)'],
    products = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(5.71,'s^-1'), n=3.021, Ea=(105.562,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 319 used for R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    products = ['C=C(O)C([O])([CH]C)C(O)CC(12260)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(210000,'s^-1'), n=1.76, Ea=(53.8899,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 326 used for R4H_SSS;O_rad_out;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    products = ['C=C(O)C(O)(CC)C([O])[CH]C(12261)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(210000,'s^-1'), n=1.76, Ea=(53.8899,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 326 used for R4H_SSS;O_rad_out;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]CC(O)(C(=C)O)C([O])CC(12262)'],
    products = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]CC(O)C([O])(CC)C(=C)O(12263)'],
    products = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    products = ['C=C([O])C(O)(CC)C([O])CC(12264)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1070.11,'s^-1'), n=2.50856, Ea=(101.808,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;O_H_out] + [R4H_SSS;O_rad_out;XH_out] + [R4H_SS(Cd)S;Y_rad_out;XH_out] for rate rule [R4H_SS(Cd)S;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C(O)C(O)(CC)C([O])CC(12265)'],
    products = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    products = ['[CH2]CC([O])(C(=C)O)C(O)CC(12266)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.2e+11,'s^-1'), n=0, Ea=(31.8402,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(1000,'K'), comment="""From training reaction 306 used for R5H_CCC;O_rad_out;Cs_H_out_2H
Exact match found for rate rule [R5H_CCC;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    products = ['C=C([O])C([O])(CC)C(O)CC(12267)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(22219.5,'s^-1'), n=1.84103, Ea=(37.7069,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5H_CCC;O_rad_out;XH_out] + [R5H_CCC(Cd);Y_rad_out;XH_out] for rate rule [R5H_CCC(Cd);O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    products = ['[CH]=C(O)C([O])(CC)C(O)CC(12268)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.468e+06,'s^-1','*|/',3), n=1.554, Ea=(111.445,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 290 used for R5H_SSSD;O_rad_out;Cd_H_out_singleH
Exact match found for rate rule [R5H_SSSD;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    products = ['[CH2]CC([O])C(O)(CC)C(=C)O(12269)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.2e+11,'s^-1'), n=0, Ea=(31.8402,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(1000,'K'), comment="""From training reaction 306 used for R5H_CCC;O_rad_out;Cs_H_out_2H
Exact match found for rate rule [R5H_CCC;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CC[CH][O](563)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    products = ['CCC([O])C1(CC)OC[C]1O(12270)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.21748e+08,'s^-1'), n=0.95, Ea=(124.892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R4_S_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    products = ['CCC1OC[C](O)C1([O])CC(12076)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.90996e+09,'s^-1'), n=0.623333, Ea=(61.7837,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    products = ['C=C(O)C(O)(CC)C(=O)CC(12271)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CH2(S)(23)', 'C=C(O)C([O])(CC)C(C)[O](6325)'],
    products = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction29',
    reactants = ['CH2(S)(23)', 'C=C(O)C(C)([O])C([O])CC(12272)'],
    products = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    products = ['C=C(O)C1(CC)OOC1CC(11988)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=[C]C(CC)(OO)C([O])CC(12273)'],
    products = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(3.01978e+11,'s^-1'), n=0, Ea=(107.111,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2OOH_S;Y_rad_out] for rate rule [R2OOH_S;Cd_rad_out_double]
Euclidian distance = 2.0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=[C]C([O])(CC)C(CC)OO(12274)'],
    products = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(3.95074e+10,'s^-1'), n=0, Ea=(112.549,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OOH_SS;Y_rad_out] for rate rule [R3OOH_SS;Cd_rad_out_double]
Euclidian distance = 2.0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    products = ['CCC([O])C([O])(CC)C(C)=O(12275)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(204.179,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction34',
    reactants = ['O(4)', 'C=C(O)[C](CC)C([O])CC(12276)'],
    products = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/Cs2;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction35',
    reactants = ['O(4)', 'C=C(O)C([O])([CH]CC)CC(6800)'],
    products = ['C=C(O)C([O])(CC)C([O])CC(11978)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/NonDeC;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

network(
    label = '2220',
    isomers = [
        'C=C(O)C([O])(CC)C([O])CC(11978)',
    ],
    reactants = [
        ('C2H5CHO(70)', 'C=C(O)C(=O)CC(4626)'),
        ('C2H5CHO(70)', '[CH2]C(O)=C([O])CC(4557)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '2220',
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

