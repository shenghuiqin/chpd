species(
    label = '[CH2]C(CC)O[C](C)OO(11852)',
    structure = SMILES('[CH2]C(CC)O[C](C)OO'),
    E0 = (-44.4539,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3615,1310,387.5,850,1000,200,800,1200,1600],'cm^-1')),
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
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.92285,0.116268,-0.000138181,9.40864e-08,-2.66874e-11,-5176.32,37.1866], Tmin=(100,'K'), Tmax=(847.207,'K')), NASAPolynomial(coeffs=[12.6524,0.052175,-2.47041e-05,4.79267e-09,-3.38309e-13,-7476.57,-26.0555], Tmin=(847.207,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-44.4539,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CJC(C)OC) + radical(Cs_P)"""),
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
    label = 'butene1(127)',
    structure = SMILES('C=CCC'),
    E0 = (-16.4325,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,252.555],'cm^-1')),
        HinderedRotor(inertia=(0.178654,'amu*angstrom^2'), symmetry=1, barrier=(7.72883,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185589,'amu*angstrom^2'), symmetry=1, barrier=(7.72103,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.176,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.58773,0.0232778,1.93412e-05,-3.55496e-08,1.36906e-11,-1918.73,14.5751], Tmin=(100,'K'), Tmax=(1007.28,'K')), NASAPolynomial(coeffs=[7.20517,0.0236362,-9.0315e-06,1.65393e-09,-1.16019e-13,-3797.34,-12.4426], Tmin=(1007.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-16.4325,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), label="""butene1""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = 'C=C(CC)O[C](C)OO(16305)',
    structure = SMILES('C=C(CC)O[C](C)OO'),
    E0 = (-173.424,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,360,370,350,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3615,1310,387.5,850,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (131.15,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.19173,0.102453,-9.22804e-05,4.15832e-08,-7.4049e-12,-20661.4,36.8541], Tmin=(100,'K'), Tmax=(1356.41,'K')), NASAPolynomial(coeffs=[23.3516,0.0300754,-1.22413e-05,2.24463e-09,-1.54411e-13,-27319.6,-89.0351], Tmin=(1356.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-173.424,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cs_P)"""),
)

species(
    label = '[CH2]C(CC)OC(=C)OO(16306)',
    structure = SMILES('[CH2]C(CC)OC(=C)OO'),
    E0 = (-67.8422,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (131.15,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.02117,0.115583,-0.000137647,8.8966e-08,-2.33126e-11,-7983.09,34.4506], Tmin=(100,'K'), Tmax=(924.724,'K')), NASAPolynomial(coeffs=[15.8321,0.0426782,-1.93804e-05,3.6984e-09,-2.5908e-13,-11099.8,-45.5368], Tmin=(924.724,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-67.8422,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2][CH]CC(130)',
    structure = SMILES('[CH2][CH]CC'),
    E0 = (255.669,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,2031.24],'cm^-1')),
        HinderedRotor(inertia=(0.244974,'amu*angstrom^2'), symmetry=1, barrier=(5.63244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00192352,'amu*angstrom^2'), symmetry=1, barrier=(5.63177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244928,'amu*angstrom^2'), symmetry=1, barrier=(5.63137,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.98997,0.0287412,-9.51469e-06,4.19232e-10,1.90526e-13,30780.1,16.8971], Tmin=(100,'K'), Tmax=(2154.56,'K')), NASAPolynomial(coeffs=[12.4231,0.0182241,-7.06316e-06,1.16769e-09,-7.11818e-14,25091.4,-39.6212], Tmin=(2154.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(255.669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJC) + radical(RCCJ)"""),
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
    label = 'C=CO[C](C)OO(10134)',
    structure = SMILES('C=CO[C](C)OO'),
    E0 = (-108.987,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,360,370,350,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.194598,0.0707156,-4.75626e-05,-3.78975e-10,8.08278e-12,-12959.4,26.9515], Tmin=(100,'K'), Tmax=(976.074,'K')), NASAPolynomial(coeffs=[21.2204,0.013669,-4.64328e-06,8.70141e-10,-6.52942e-14,-18451,-81.0823], Tmin=(976.074,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-108.987,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cs_P)"""),
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
    label = '[CH2]C(CC)OC(C)=O(16307)',
    structure = SMILES('[CH2]C(CC)OC(C)=O'),
    E0 = (-313.181,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.562159,0.0762144,-5.41187e-05,2.03217e-08,-3.24242e-12,-37544.2,30.0671], Tmin=(100,'K'), Tmax=(1401.42,'K')), NASAPolynomial(coeffs=[12.1406,0.0431669,-1.87467e-05,3.49497e-09,-2.40709e-13,-40789.5,-29.6996], Tmin=(1401.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-313.181,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsOs) + radical(CJC(C)OC)"""),
)

species(
    label = 'CC[C](C)O[C](C)OO(16308)',
    structure = SMILES('CC[C](C)O[C](C)OO'),
    E0 = (-74.2602,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.046992,0.10382,-8.84203e-05,-1.19157e-08,5.48107e-11,-8800.47,34.3978], Tmin=(100,'K'), Tmax=(513.551,'K')), NASAPolynomial(coeffs=[9.12865,0.0576931,-2.77106e-05,5.36132e-09,-3.75779e-13,-10077.1,-7.00826], Tmin=(513.551,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-74.2602,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_P) + radical(C2CsJOCs)"""),
)

species(
    label = '[CH2]C(CC)OC([CH2])OO(16309)',
    structure = SMILES('[CH2]C(CC)OC([CH2])OO'),
    E0 = (-35.7377,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1310,387.5,850,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,200,800,1200,1600],'cm^-1')),
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
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.03897,0.118427,-0.000141545,9.53064e-08,-2.65166e-11,-4123.43,38.9097], Tmin=(100,'K'), Tmax=(865.823,'K')), NASAPolynomial(coeffs=[13.6805,0.0504244,-2.37334e-05,4.59361e-09,-3.2391e-13,-6672.32,-29.9827], Tmin=(865.823,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-35.7377,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CJC(C)OC) + radical(CJCOOH)"""),
)

species(
    label = '[CH2][C](CC)OC(C)OO(16310)',
    structure = SMILES('[CH2][C](CC)OC(C)OO'),
    E0 = (-68.9964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.42197,0.130934,-0.000196908,1.70632e-07,-5.92112e-11,-8114.36,38.2099], Tmin=(100,'K'), Tmax=(816.494,'K')), NASAPolynomial(coeffs=[10.1781,0.0562755,-2.69976e-05,5.15552e-09,-3.55898e-13,-9414.34,-11.7627], Tmin=(816.494,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-68.9964,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CJC(C)OC) + radical(C2CsJOCs)"""),
)

species(
    label = '[CH2]C(CC)OC(C)O[O](16311)',
    structure = SMILES('[CH2]C(CC)OC(C)O[O]'),
    E0 = (-97.6956,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,492.5,1135,1000,3000,3100,440,815,1455,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,200,800,1200,1600],'cm^-1')),
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
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.809383,0.113372,-0.000137893,9.89341e-08,-2.96783e-11,-11583.8,36.8441], Tmin=(100,'K'), Tmax=(802.925,'K')), NASAPolynomial(coeffs=[11.3898,0.052598,-2.43569e-05,4.66465e-09,-3.2628e-13,-13542.8,-19.3323], Tmin=(802.925,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-97.6956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CJC(C)OC)"""),
)

species(
    label = 'C[CH]C(C)O[C](C)OO(16312)',
    structure = SMILES('C[CH]C(C)O[C](C)OO'),
    E0 = (-55.0619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.301945,0.0977021,-8.51796e-05,3.99801e-08,-7.88473e-12,-6470.22,37.6812], Tmin=(100,'K'), Tmax=(1177.7,'K')), NASAPolynomial(coeffs=[14.2898,0.0481417,-2.2056e-05,4.24737e-09,-2.99436e-13,-9907.16,-35.1024], Tmin=(1177.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-55.0619,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_P) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C([CH]C)OC(C)OO(16313)',
    structure = SMILES('[CH2]C([CH]C)OC(C)OO'),
    E0 = (-49.7982,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.703012,0.109642,-0.000116975,6.9832e-08,-1.73901e-11,-5825.17,38.2024], Tmin=(100,'K'), Tmax=(956.837,'K')), NASAPolynomial(coeffs=[13.6055,0.0498243,-2.31982e-05,4.49199e-09,-3.17699e-13,-8563.27,-30.196], Tmin=(956.837,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-49.7982,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJCO) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]CC(C)O[C](C)OO(16314)',
    structure = SMILES('[CH2]CC(C)O[C](C)OO'),
    E0 = (-49.7176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.355701,0.102225,-9.84597e-05,5.33055e-08,-1.22847e-11,-5828.13,37.1795], Tmin=(100,'K'), Tmax=(1015.6,'K')), NASAPolynomial(coeffs=[12.3528,0.0521718,-2.45328e-05,4.77781e-09,-3.39126e-13,-8409.48,-24.3283], Tmin=(1015.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-49.7176,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_P) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C](OO)OC(C)CC(16315)',
    structure = SMILES('[CH2][C](OO)OC(C)CC'),
    E0 = (-41.0014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.504691,0.104783,-0.000103273,5.64939e-08,-1.3001e-11,-4773.84,37.9212], Tmin=(100,'K'), Tmax=(1024.02,'K')), NASAPolynomial(coeffs=[13.4351,0.0503313,-2.35114e-05,4.56698e-09,-3.23765e-13,-7628.77,-29.6613], Tmin=(1024.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-41.0014,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_P) + radical(CJCOOH)"""),
)

species(
    label = '[CH2]CC([CH2])OC(C)OO(16316)',
    structure = SMILES('[CH2]CC([CH2])OC(C)OO'),
    E0 = (-44.4539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.92285,0.116268,-0.000138181,9.40864e-08,-2.66874e-11,-5176.32,38.2852], Tmin=(100,'K'), Tmax=(847.207,'K')), NASAPolynomial(coeffs=[12.6524,0.052175,-2.47041e-05,4.79267e-09,-3.38309e-13,-7476.57,-24.9569], Tmin=(847.207,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-44.4539,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CJC(C)OC) + radical(RCCJ)"""),
)

species(
    label = 'CCC(C)O[C](C)O[O](16317)',
    structure = SMILES('CCC(C)O[C](C)O[O]'),
    E0 = (-102.959,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.19501,0.0986545,-9.52701e-05,5.36263e-08,-1.30061e-11,-12237.4,35.5769], Tmin=(100,'K'), Tmax=(964.584,'K')), NASAPolynomial(coeffs=[10.7756,0.0531614,-2.45261e-05,4.73264e-09,-3.34125e-13,-14353.9,-16.9547], Tmin=(964.584,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-102.959,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_P) + radical(ROOJ)"""),
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
    label = 'C=C(CC)OC(C)OO(16318)',
    structure = SMILES('C=C(CC)OC(C)OO'),
    E0 = (-378.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.44131,0.104079,-8.94446e-05,3.83694e-08,-6.49316e-12,-45334.6,37.0717], Tmin=(100,'K'), Tmax=(1427.06,'K')), NASAPolynomial(coeffs=[24.8863,0.0302816,-1.18733e-05,2.13013e-09,-1.44388e-13,-52848.5,-99.3054], Tmin=(1427.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-378.67,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C(OO)OC(C)CC(16319)',
    structure = SMILES('C=C(OO)OC(C)CC'),
    E0 = (-278.352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.693424,0.10329,-9.65496e-05,4.81519e-08,-9.85394e-12,-33309.1,33.5081], Tmin=(100,'K'), Tmax=(1159.8,'K')), NASAPolynomial(coeffs=[16.6774,0.043379,-1.90631e-05,3.61067e-09,-2.52669e-13,-37338.4,-52.8711], Tmin=(1159.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-278.352,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH)"""),
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
    label = '[CH2]C(C)O[C](C)OO(8913)',
    structure = SMILES('[CH2]C(C)O[C](C)OO'),
    E0 = (-20.6736,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,200,800,1600],'cm^-1')),
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
    molecularWeight = (118.131,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3963.94,'J/mol'), sigma=(6.91174,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=619.16 K, Pc=27.24 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.314811,0.101988,-0.000126306,8.88367e-08,-2.58477e-11,-2337.39,32.755], Tmin=(100,'K'), Tmax=(829.223,'K')), NASAPolynomial(coeffs=[11.7722,0.0436844,-2.0841e-05,4.04819e-09,-2.85624e-13,-4342,-23.2946], Tmin=(829.223,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-20.6736,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CJC(C)OC) + radical(Cs_P)"""),
)

species(
    label = 'CCC[CH]O[C](C)OO(11854)',
    structure = SMILES('CCC[CH]O[C](C)OO'),
    E0 = (-63.9571,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,360,370,350,3615,1310,387.5,850,1000,3025,407.5,1350,352.5,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.800774,0.111935,-0.000122322,7.47414e-08,-1.89835e-11,-7524.72,36.7838], Tmin=(100,'K'), Tmax=(941.057,'K')), NASAPolynomial(coeffs=[13.8375,0.0497124,-2.31392e-05,4.47606e-09,-3.1628e-13,-10279.7,-32.9475], Tmin=(941.057,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-63.9571,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCsJOCs) + radical(Cs_P)"""),
)

species(
    label = 'CCC1CC(C)(OO)O1(13485)',
    structure = SMILES('CCC1CC(C)(OO)O1'),
    E0 = (-333.608,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.731421,0.0766185,2.29627e-05,-1.0728e-07,5.43954e-11,-39927.7,30.3982], Tmin=(100,'K'), Tmax=(904.739,'K')), NASAPolynomial(coeffs=[29.3159,0.0172553,-4.4362e-07,-2.63803e-10,1.88771e-14,-48372.1,-128.175], Tmin=(904.739,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-333.608,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(486.397,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane)"""),
)

species(
    label = '[CH2]C(CC)OC(C)([O])O(16320)',
    structure = SMILES('[CH2]C(CC)OC(C)([O])O'),
    E0 = (-275.878,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.880559,0.118001,-0.000169693,1.46815e-07,-5.1151e-11,-33015,38.277], Tmin=(100,'K'), Tmax=(828.461,'K')), NASAPolynomial(coeffs=[7.97198,0.057053,-2.63761e-05,4.96126e-09,-3.39745e-13,-33857,1.00541], Tmin=(828.461,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-275.878,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsOs) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CJC(C)OC)"""),
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
    label = '[CH2]C([O])CC(1666)',
    structure = SMILES('[CH2]C([O])CC'),
    E0 = (125.894,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,326.431,326.435],'cm^-1')),
        HinderedRotor(inertia=(0.00158206,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102503,'amu*angstrom^2'), symmetry=1, barrier=(7.75102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102506,'amu*angstrom^2'), symmetry=1, barrier=(7.75103,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69126,0.0496706,-3.78506e-05,1.56223e-08,-2.69593e-12,15225.4,19.2246], Tmin=(100,'K'), Tmax=(1336.85,'K')), NASAPolynomial(coeffs=[9.88284,0.0251604,-1.03492e-05,1.90768e-09,-1.31205e-13,13035.2,-22.6731], Tmin=(1336.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(125.894,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo library: CBS_QB3_1dHR + radical(CJCO) + radical(CC(C)OJ)"""),
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
    label = 'CC[CH]O[C](C)OO(11606)',
    structure = SMILES('CC[CH]O[C](C)OO'),
    E0 = (-40.1768,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,360,370,350,3615,1310,387.5,850,1000,3025,407.5,1350,352.5,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,200,800,1200,1600],'cm^-1')),
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
    molecularWeight = (118.131,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.183322,0.0975184,-0.000109849,6.85421e-08,-1.76601e-11,-4686.16,32.3201], Tmin=(100,'K'), Tmax=(931.221,'K')), NASAPolynomial(coeffs=[12.8973,0.041331,-1.93423e-05,3.74777e-09,-2.64976e-13,-7122.34,-29.8542], Tmin=(931.221,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-40.1768,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCsJOCs) + radical(Cs_P)"""),
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
    E0 = (-44.4539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (49.6582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (146.934,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-44.4539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (21.4395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-44.4539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (58.2256,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (96.5469,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (129.252,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (56.4844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (102.404,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (19.2893,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (38.3893,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (9.7289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (5.7541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (29.1944,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (301.323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (18.9463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-9.9359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (399.188,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (114.173,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-36.1696,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (41.7364,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (409.149,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (341.386,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(CC)O[C](C)OO(11852)'],
    products = ['CH3C(O)OOH(67)', 'butene1(127)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(3)', 'C=C(CC)O[C](C)OO(16305)'],
    products = ['[CH2]C(CC)O[C](C)OO(11852)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(170.641,'m^3/(mol*s)'), n=1.56204, Ea=(11.2897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;HJ] for rate rule [Cds-OsCs_Cds;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', '[CH2]C(CC)OC(=C)OO(16306)'],
    products = ['[CH2]C(CC)O[C](C)OO(11852)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1999.59,'m^3/(mol*s)'), n=1.33685, Ea=(2.9843,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-HH_Cds;HJ] for rate rule [Cds-HH_Cds-OsOs;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2][CH]CC(130)', 'CH3C(O)OOH(67)'],
    products = ['[CH2]C(CC)O[C](C)OO(11852)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.0323,'m^3/(mol*s)'), n=2.98, Ea=(67.7388,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Od_CO-NdNd;YJ] for rate rule [Od_CO-NdNd;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 67.5 to 67.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction5',
    reactants = ['C2H5(29)', 'C=CO[C](C)OO(10134)'],
    products = ['[CH2]C(CC)O[C](C)OO(11852)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00139222,'m^3/(mol*s)'), n=2.42243, Ea=(22.5521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;CsJ-CsHH] for rate rule [Cds-OsH_Cds;CsJ-CsHH]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['OH(5)', '[CH2]C(CC)OC(C)=O(16307)'],
    products = ['[CH2]C(CC)O[C](C)OO(11852)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.0323,'m^3/(mol*s)'), n=2.98, Ea=(240.355,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Od_CO-NdNd;YJ] for rate rule [Od_CO-NdNd;OJ_pri]
Euclidian distance = 2.0
family: R_Addition_MultipleBond
Ea raised from 237.9 to 240.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C(CC)O[C](C)OO(11852)'],
    products = ['CC[C](C)O[C](C)OO(16308)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(5.265e-07,'s^-1'), n=5.639, Ea=(102.68,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_2H;Cs_H_out_NonDe] for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_NDMustO]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C(CC)OC([CH2])OO(16309)'],
    products = ['[CH2]C(CC)O[C](C)OO(11852)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(309968,'s^-1'), n=2.08546, Ea=(132.285,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_2H;Cs_H_out] for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_OOH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(CC)O[C](C)OO(11852)'],
    products = ['[CH2][C](CC)OC(C)OO(16310)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(9.31743e+07,'s^-1'), n=1.36595, Ea=(173.706,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_NonDe;XH_out] for rate rule [R3H_SS_O;C_rad_out_NDMustO;XH_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(CC)OC(C)O[O](16311)'],
    products = ['[CH2]C(CC)O[C](C)OO(11852)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3e+08,'s^-1'), n=1.23, Ea=(154.18,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 270 used for R3H_SS_O;O_rad_out;Cs_H_out_NDMustO
Exact match found for rate rule [R3H_SS_O;O_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(CC)O[C](C)OO(11852)'],
    products = ['C[CH]C(C)O[C](C)OO(16312)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.18e+10,'s^-1'), n=0.82, Ea=(146.858,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 186 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(CC)O[C](C)OO(11852)'],
    products = ['[CH2]C([CH]C)OC(C)OO(16313)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(398.595,'s^-1'), n=2.415, Ea=(63.7433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_single;Cs_H_out_H/NonDeC] for rate rule [R4H_SSS_OCs;C_rad_out_NDMustO;Cs_H_out_H/NonDeC]
Euclidian distance = 3.16227766017
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(CC)O[C](C)OO(11852)'],
    products = ['[CH2]CC(C)O[C](C)OO(16314)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(114000,'s^-1'), n=1.74, Ea=(82.8432,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 109 used for R4H_SSS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(CC)O[C](C)OO(11852)'],
    products = ['[CH2][C](OO)OC(C)CC(16315)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(91273.5,'s^-1'), n=1.79, Ea=(54.1828,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;Cs_H_out_2H] for rate rule [R5HJ_3;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(CC)O[C](C)OO(11852)'],
    products = ['[CH2]CC([CH2])OC(C)OO(16316)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(38738.3,'s^-1'), n=1.755, Ea=(50.208,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;C_rad_out_single;Cs_H_out_2H] for rate rule [R5H_SSSS_OCC;C_rad_out_NDMustO;Cs_H_out_2H]
Euclidian distance = 3.16227766017
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C(CC)O[C](C)OO(11852)'],
    products = ['CCC(C)O[C](C)O[O](16317)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1389.4,'s^-1'), n=2.18252, Ea=(73.6483,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;C_rad_out_2H;XH_out] for rate rule [R6HJ_3;C_rad_out_2H;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][CH]CC(130)', 'C[C]([O])OO(5646)'],
    products = ['[CH2]C(CC)O[C](C)OO(11852)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(CC)O[C](C)OO(11852)'],
    products = ['C=C(CC)OC(C)OO(16318)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(CC)O[C](C)OO(11852)'],
    products = ['C=C(OO)OC(C)CC(16319)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(9.63e+09,'s^-1'), n=0.137, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_NDe] for rate rule [R5radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['CH2(S)(23)', '[CH2]C(C)O[C](C)OO(8913)'],
    products = ['[CH2]C(CC)O[C](C)OO(11852)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(CC)O[C](C)OO(11852)'],
    products = ['CCC[CH]O[C](C)OO(11854)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(CC)O[C](C)OO(11852)'],
    products = ['CCC1CC(C)(OO)O1(13485)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_NDMustO;Cpri_rad_out_2H]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C(CC)O[C](C)OO(11852)'],
    products = ['[CH2]C(CC)OC(C)([O])O(16320)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.04257e+10,'s^-1'), n=0, Ea=(86.1903,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnOOH;C_rad_out_NonDe] for rate rule [ROOH;C_rad_out_NDMustO]
Euclidian distance = 1.41421356237
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C[C]OO(178)', '[CH2]C([O])CC(1666)'],
    products = ['[CH2]C(CC)O[C](C)OO(11852)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['CH2(19)', 'CC[CH]O[C](C)OO(11606)'],
    products = ['[CH2]C(CC)O[C](C)OO(11852)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/CsO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '2862',
    isomers = [
        '[CH2]C(CC)O[C](C)OO(11852)',
    ],
    reactants = [
        ('CH3C(O)OOH(67)', 'butene1(127)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '2862',
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

