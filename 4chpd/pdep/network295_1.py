species(
    label = '[CH2]CCCC=CC=O(1870)',
    structure = SMILES('[CH2]CCCC=CC=O'),
    E0 = (17.676,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.162,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3796.55,'J/mol'), sigma=(6.47699,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=593.01 K, Pc=31.7 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.917115,0.0717441,-4.53172e-05,1.37194e-08,-1.70697e-12,2232.04,30.4737], Tmin=(100,'K'), Tmax=(1738.36,'K')), NASAPolynomial(coeffs=[14.2937,0.0409643,-1.87579e-05,3.5339e-09,-2.42159e-13,-2418.65,-41.4573], Tmin=(1738.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(17.676,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(RCCJ)"""),
)

species(
    label = 'H(19)',
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
    label = 'C=CCCC=CC=O(1889)',
    structure = SMILES('C=CCCC=CC=O'),
    E0 = (-59.7881,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,228.993,229.051,229.059],'cm^-1')),
        HinderedRotor(inertia=(0.00321524,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.37984,'amu*angstrom^2'), symmetry=1, barrier=(14.136,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.379796,'amu*angstrom^2'), symmetry=1, barrier=(14.136,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.37974,'amu*angstrom^2'), symmetry=1, barrier=(14.136,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.154,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.64578,0.0690372,-4.33348e-05,1.27352e-08,-1.48766e-12,-7066.73,30.1787], Tmin=(100,'K'), Tmax=(1942.67,'K')), NASAPolynomial(coeffs=[19.5141,0.0301875,-1.33382e-05,2.44137e-09,-1.62985e-13,-14397.8,-73.3802], Tmin=(1942.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-59.7881,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'C2H4(45)',
    structure = SMILES('C=C'),
    E0 = (42.0619,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2334.71,'J/mol'), sigma=(3.971,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.9592,-0.00757051,5.7099e-05,-6.91588e-08,2.69884e-11,5089.78,4.0973], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.99183,0.0104834,-3.71721e-06,5.94628e-10,-3.5363e-14,4268.66,-0.269082], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(42.0619,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""C2H4""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH2]CC=CC=O(1890)',
    structure = SMILES('[CH2]CC=CC=O'),
    E0 = (65.2365,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.566367,'amu*angstrom^2'), symmetry=1, barrier=(13.0219,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0103939,'amu*angstrom^2'), symmetry=1, barrier=(5.37868,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.567537,'amu*angstrom^2'), symmetry=1, barrier=(13.0488,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90082,0.0452877,-2.64603e-05,6.86669e-09,-6.89619e-13,7922.04,22.494], Tmin=(100,'K'), Tmax=(2280.95,'K')), NASAPolynomial(coeffs=[18.0657,0.0169393,-7.81719e-06,1.4176e-09,-9.23615e-14,548.021,-68.8211], Tmin=(2280.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(65.2365,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[CH2](171)',
    structure = SMILES('[CH2]C[CH2]'),
    E0 = (290.499,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.000703322,'amu*angstrom^2'), symmetry=1, barrier=(7.98556,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0167612,'amu*angstrom^2'), symmetry=1, barrier=(29.9599,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.12863,0.0145638,1.31003e-05,-2.16419e-08,7.8594e-12,34974.2,14.113], Tmin=(100,'K'), Tmax=(1043.7,'K')), NASAPolynomial(coeffs=[5.43609,0.0169627,-6.50459e-06,1.20126e-09,-8.4342e-14,33880.3,-0.0512239], Tmin=(1043.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(290.499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=CC=O(1897)',
    structure = SMILES('C=CC=C[O]'),
    E0 = (26.5912,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.32203,'amu*angstrom^2'), symmetry=1, barrier=(30.3961,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05417,0.0257395,3.62618e-05,-7.76905e-08,3.52495e-11,3284.1,16.2493], Tmin=(100,'K'), Tmax=(919.876,'K')), NASAPolynomial(coeffs=[17.5824,0.00251509,1.89717e-06,-4.33504e-10,2.49804e-14,-1446.91,-67.5556], Tmin=(919.876,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(26.5912,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ)"""),
)

species(
    label = 'C2H4(T)(169)',
    structure = SMILES('[CH2][CH2]'),
    E0 = (318.146,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,180,1436.54,1437.15,2688.96,2689.16],'cm^-1')),
        HinderedRotor(inertia=(0.0257549,'amu*angstrom^2'), symmetry=1, barrier=(17.2441,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.40736,0.0100312,6.40927e-06,-1.41291e-08,5.92671e-12,38288.2,6.11703], Tmin=(100,'K'), Tmax=(954.26,'K')), NASAPolynomial(coeffs=[5.52249,0.00856173,-2.90743e-06,5.02353e-10,-3.44572e-14,37547.8,-5.75276], Tmin=(954.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.146,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""C2H4(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C[CH]CC=CC=O(1898)',
    structure = SMILES('[CH2]C[CH]CC=CC=O'),
    E0 = (212.134,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.154,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.873816,0.0754692,-8.05145e-05,6.64511e-08,-2.59671e-11,25620.2,33.2216], Tmin=(100,'K'), Tmax=(689.248,'K')), NASAPolynomial(coeffs=[3.56081,0.0553957,-2.708e-05,5.33754e-09,-3.80094e-13,25356.2,22.0303], Tmin=(689.248,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(212.134,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(RCCJCC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC[CH2](390)',
    structure = SMILES('[CH2]CC[CH2]'),
    E0 = (266.469,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2188.14],'cm^-1')),
        HinderedRotor(inertia=(0.165367,'amu*angstrom^2'), symmetry=1, barrier=(4.74005,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00483329,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.176255,'amu*angstrom^2'), symmetry=1, barrier=(4.7527,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2896.88,'J/mol'), sigma=(5.41451,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=452.49 K, Pc=41.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.35532,0.0344795,-1.7272e-05,4.08056e-09,-3.78294e-13,32109.4,17.5603], Tmin=(100,'K'), Tmax=(2451.48,'K')), NASAPolynomial(coeffs=[13.9856,0.0155026,-5.66049e-06,9.22849e-10,-5.6271e-14,26407.1,-48.9781], Tmin=(2451.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(266.469,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = 'CHCHCHO(1899)',
    structure = SMILES('[CH]=CC=O'),
    E0 = (171.79,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3120,650,792.5,1650,691.226,691.35,691.975,692.297,692.622],'cm^-1')),
        HinderedRotor(inertia=(0.00997499,'amu*angstrom^2'), symmetry=1, barrier=(3.38584,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (55.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.95468,0.0145151,1.89382e-05,-3.59857e-08,1.47962e-11,20706.8,12.171], Tmin=(100,'K'), Tmax=(981.049,'K')), NASAPolynomial(coeffs=[9.63601,0.00810119,-3.09994e-06,6.30252e-10,-4.9152e-14,18393.6,-25.043], Tmin=(981.049,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(171.79,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CHCHCHO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]CC[CH]C=CC=O(1900)',
    structure = SMILES('[CH2]CCC=CC=C[O]'),
    E0 = (149.386,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,329.598,329.651,329.675,329.69],'cm^-1')),
        HinderedRotor(inertia=(0.224928,'amu*angstrom^2'), symmetry=1, barrier=(17.3402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224793,'amu*angstrom^2'), symmetry=1, barrier=(17.3387,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.225011,'amu*angstrom^2'), symmetry=1, barrier=(17.3391,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224852,'amu*angstrom^2'), symmetry=1, barrier=(17.3394,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.154,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0364764,0.070078,-1.60729e-05,-3.75396e-08,2.17413e-11,18125.2,32.2394], Tmin=(100,'K'), Tmax=(958.329,'K')), NASAPolynomial(coeffs=[21.0825,0.0234145,-7.49215e-06,1.33158e-09,-9.66323e-14,12200.3,-78.2661], Tmin=(958.329,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(149.386,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(RCCJ) + radical(C=COJ)"""),
)

species(
    label = '[CH2][CH]CCC=CC=O(1901)',
    structure = SMILES('[CH2][CH]CCC=CC=O'),
    E0 = (212.122,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.154,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01577,0.072524,-6.47125e-05,4.073e-08,-1.268e-11,25613.7,32.9053], Tmin=(100,'K'), Tmax=(722.115,'K')), NASAPolynomial(coeffs=[4.36227,0.0539859,-2.62026e-05,5.17551e-09,-3.70201e-13,25130.4,17.8501], Tmin=(722.115,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(212.122,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CCC[C]=CC=O(1902)',
    structure = SMILES('[CH2]CCC[C]=CC=O'),
    E0 = (255.518,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,1685,370,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.154,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.970414,0.0747362,-6.24412e-05,3.12902e-08,-7.28979e-12,30834,30.8056], Tmin=(100,'K'), Tmax=(950.057,'K')), NASAPolynomial(coeffs=[6.62182,0.0509424,-2.48747e-05,4.92949e-09,-3.53243e-13,29760.2,3.83033], Tmin=(950.057,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(255.518,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_S) + radical(RCCJ)"""),
)

species(
    label = 'HCO(30)',
    structure = SMILES('[CH]=O'),
    E0 = (32.4782,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1131.19,1955.83,1955.83],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.23755,-0.00332075,1.4003e-05,-1.3424e-08,4.37416e-12,3872.41,3.30835], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.92002,0.00252279,-6.71004e-07,1.05616e-10,-7.43798e-15,3653.43,3.58077], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(32.4782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HCO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH]=CCCC[CH2](1903)',
    structure = SMILES('[CH]=CCCC[CH2]'),
    E0 = (388.26,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3120,650,792.5,1650,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3100,440,815,1455,1000,180,803.677],'cm^-1')),
        HinderedRotor(inertia=(0.132333,'amu*angstrom^2'), symmetry=1, barrier=(3.04261,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130856,'amu*angstrom^2'), symmetry=1, barrier=(3.00864,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133479,'amu*angstrom^2'), symmetry=1, barrier=(3.06894,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.676879,'amu*angstrom^2'), symmetry=1, barrier=(15.5628,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.924506,0.0584054,-3.73455e-05,1.20693e-08,-1.58147e-12,46815.4,27.8224], Tmin=(100,'K'), Tmax=(1759.65,'K')), NASAPolynomial(coeffs=[15.0124,0.0263809,-1.00462e-05,1.72652e-09,-1.12027e-13,41857.5,-48.1045], Tmin=(1759.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(388.26,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CCCC=[C]C=O(1904)',
    structure = SMILES('[CH2]CCCC=C=C[O]'),
    E0 = (214.413,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,540,610,2055,3000,3100,440,815,1455,1000,303.301,303.302,303.302,303.303],'cm^-1')),
        HinderedRotor(inertia=(0.274535,'amu*angstrom^2'), symmetry=1, barrier=(17.9205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.274521,'amu*angstrom^2'), symmetry=1, barrier=(17.9205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.274525,'amu*angstrom^2'), symmetry=1, barrier=(17.9205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00183254,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.154,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0520094,0.0771262,-4.69808e-05,3.29624e-09,4.90661e-12,25944.2,33.0271], Tmin=(100,'K'), Tmax=(1026.08,'K')), NASAPolynomial(coeffs=[18.086,0.0295514,-1.12501e-05,2.05313e-09,-1.43819e-13,21004.2,-60.8798], Tmin=(1026.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.413,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CCCC=C[C]=O(1905)',
    structure = SMILES('[CH2]CCCC=C[C]=O'),
    E0 = (178.297,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,1855,455,950,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.154,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.56413,0.082779,-7.91206e-05,4.4157e-08,-1.08022e-11,21561.9,31.3126], Tmin=(100,'K'), Tmax=(945.759,'K')), NASAPolynomial(coeffs=[8.94691,0.0473233,-2.28842e-05,4.51417e-09,-3.22641e-13,19976.4,-8.66158], Tmin=(945.759,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(178.297,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(RCCJ) + radical(C=CCJ=O)"""),
)

species(
    label = 'O=C[CH]C1CCCC1(1888)',
    structure = SMILES('[O]C=CC1CCCC1'),
    E0 = (-86.754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.162,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19083,0.0312223,0.000105547,-1.65125e-07,6.71728e-11,-10305,27.6782], Tmin=(100,'K'), Tmax=(942.328,'K')), NASAPolynomial(coeffs=[21.2456,0.0236543,-5.86818e-06,1.04256e-09,-8.41787e-14,-17528.3,-86.1554], Tmin=(942.328,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-86.754,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclopentane) + radical(C=COJ)"""),
)

species(
    label = '[O]C1C=CCCCC1(160)',
    structure = SMILES('[O]C1C=CCCCC1'),
    E0 = (34.5629,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2790,2830,2870,2910,2950,2990,3030,3070,3110,3150,900,920,940,960,980,1000,1020,1040,1060,1080,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.162,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4335.08,'J/mol'), sigma=(7.34722,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=677.13 K, Pc=24.8 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24659,-0.000730058,0.000255851,-3.7063e-07,1.54011e-10,4313.44,28.414], Tmin=(100,'K'), Tmax=(909.883,'K')), NASAPolynomial(coeffs=[40.9949,-0.015947,1.79534e-05,-3.63655e-09,2.32493e-13,-9523.18,-195.881], Tmin=(909.883,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(34.5629,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptene) + radical(CC(C)OJ)"""),
)

species(
    label = 'C[CH]CCC=CC=O(1891)',
    structure = SMILES('C[CH]CCC=CC=O'),
    E0 = (6.87596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.162,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62685,0.0653899,-3.62736e-05,9.12991e-09,-9.25014e-13,898.057,28.8239], Tmin=(100,'K'), Tmax=(2009.3,'K')), NASAPolynomial(coeffs=[12.6316,0.0434822,-1.99188e-05,3.70353e-09,-2.49854e-13,-3524.28,-31.9464], Tmin=(2009.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.87596,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(RCCJC)"""),
)

species(
    label = 'CC[CH]CC=CC=O(1892)',
    structure = SMILES('CC[CH]CC=CC=O'),
    E0 = (6.88799,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.162,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92208,0.0631567,-3.30215e-05,7.39017e-09,-6.30396e-13,883.842,27.562], Tmin=(100,'K'), Tmax=(2606.67,'K')), NASAPolynomial(coeffs=[23.6116,0.0298738,-1.38692e-05,2.49191e-09,-1.60618e-13,-10423.7,-97.8577], Tmin=(2606.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.88799,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(RCCJCC)"""),
)

species(
    label = 'CCC[CH]C=CC=O(1893)',
    structure = SMILES('CCCC=CC=C[O]'),
    E0 = (-55.8601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.162,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.00419796,0.0693848,-5.75947e-06,-4.96203e-08,2.6109e-11,-6557.34,30.5989], Tmin=(100,'K'), Tmax=(957.959,'K')), NASAPolynomial(coeffs=[21.3162,0.0256629,-8.2344e-06,1.46807e-09,-1.06737e-13,-12720.8,-82.1934], Tmin=(957.959,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-55.8601,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(440.667,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ)"""),
)

species(
    label = 'CCCC[C]=CC=O(1894)',
    structure = SMILES('CCCC[C]=CC=O'),
    E0 = (50.2715,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3010,987.5,1337.5,450,1655,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11006,0.0724492,-4.86246e-05,1.68422e-08,-2.57388e-12,6142.09,28.4774], Tmin=(100,'K'), Tmax=(1362.49,'K')), NASAPolynomial(coeffs=[8.87262,0.0496599,-2.35352e-05,4.56596e-09,-3.21341e-13,4026.81,-11.3736], Tmin=(1362.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(50.2715,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_S)"""),
)

species(
    label = 'CCCCC=[C]C=O(1895)',
    structure = SMILES('CCCCC=C=C[O]'),
    E0 = (9.1664,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.162,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0825568,0.0763186,-3.62957e-05,-9.21512e-09,9.43148e-12,1261.25,31.3499], Tmin=(100,'K'), Tmax=(1011.97,'K')), NASAPolynomial(coeffs=[18.2434,0.0319263,-1.2064e-05,2.2063e-09,-1.55294e-13,-3883.83,-64.3755], Tmin=(1011.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(9.1664,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(440.667,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ)"""),
)

species(
    label = 'CCCCC=C[C]=O(1896)',
    structure = SMILES('CCCCC=C[C]=O'),
    E0 = (-26.9489,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.162,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.597216,0.0814088,-6.72378e-05,3.09421e-08,-6.24606e-12,-3124.24,29.3936], Tmin=(100,'K'), Tmax=(1115.79,'K')), NASAPolynomial(coeffs=[9.7845,0.0484735,-2.29618e-05,4.48804e-09,-3.18897e-13,-5174.46,-15.9364], Tmin=(1115.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-26.9489,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJ=O)"""),
)

species(
    label = 'O=CC1[CH]CCCC1(1906)',
    structure = SMILES('O=CC1[CH]CCCC1'),
    E0 = (-54.3059,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.162,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34751,0.0392988,5.655e-05,-8.80975e-08,3.21092e-11,-6419.63,25.0912], Tmin=(100,'K'), Tmax=(1044.26,'K')), NASAPolynomial(coeffs=[12.2238,0.0449317,-1.94756e-05,3.80812e-09,-2.76395e-13,-11269.8,-40.1983], Tmin=(1044.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-54.3059,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-OdCsH) + ring(Cyclohexane) + radical(CCJCC=O)"""),
)

species(
    label = '[CH]1C=CCCCCO1(1907)',
    structure = SMILES('[CH]1C=COCCCC1'),
    E0 = (-31.0492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.162,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49064,0.0203434,0.000140834,-2.00511e-07,7.88311e-11,-3612.25,20.3424], Tmin=(100,'K'), Tmax=(949.752,'K')), NASAPolynomial(coeffs=[21.2505,0.0255571,-7.06994e-06,1.34875e-09,-1.10424e-13,-11354.2,-94.9664], Tmin=(949.752,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-31.0492,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclooctane) + radical(Allyl_S)"""),
)

species(
    label = 'CO(28)',
    structure = SMILES('[C-]#[O+]'),
    E0 = (-119.219,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2084.51],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0101,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(762.44,'J/mol'), sigma=(3.69,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.5971,-0.00102424,2.83336e-06,-1.75825e-09,3.42587e-13,-14343.2,3.45822], Tmin=(100,'K'), Tmax=(1669.93,'K')), NASAPolynomial(coeffs=[2.92796,0.00181931,-8.35308e-07,1.51269e-10,-9.88872e-15,-14292.7,6.51157], Tmin=(1669.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-119.219,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""CO""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = '[CH2]CCCC=C(1908)',
    structure = SMILES('[CH2]CCCC=C'),
    E0 = (141.164,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3100,440,815,1455,1000,180,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.08215,'amu*angstrom^2'), symmetry=1, barrier=(27.5211,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.521949,'amu*angstrom^2'), symmetry=1, barrier=(12.0006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19687,'amu*angstrom^2'), symmetry=1, barrier=(27.5184,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00242347,'amu*angstrom^2'), symmetry=1, barrier=(27.5163,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1515,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17125,0.05314,-1.58866e-05,-9.87118e-09,5.64078e-12,17087.5,26.5128], Tmin=(100,'K'), Tmax=(1116.79,'K')), NASAPolynomial(coeffs=[10.9831,0.0348882,-1.40594e-05,2.58114e-09,-1.78433e-13,13842.6,-26.6236], Tmin=(1116.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(141.164,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CCCC=C=CO(1909)',
    structure = SMILES('[CH2]CCCC=C=CO'),
    E0 = (72.9501,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,3615,1277.5,1000,540,610,2055,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.434622,0.0837899,-4.95015e-05,-3.73617e-09,9.82974e-12,8945.88,32.9884], Tmin=(100,'K'), Tmax=(969.032,'K')), NASAPolynomial(coeffs=[21.0346,0.0263683,-8.9115e-06,1.56465e-09,-1.09657e-13,3320.15,-77.4711], Tmin=(969.032,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(72.9501,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(RCCJ)"""),
)

species(
    label = 'CH2(T)(33)',
    structure = SMILES('[CH2]'),
    E0 = (381.08,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([971.045,2816.03,3444.23],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.71758,0.00127391,2.17347e-06,-3.48858e-09,1.65209e-12,45872.4,1.75298], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.14632,0.00303671,-9.96474e-07,1.50484e-10,-8.57336e-15,46041.3,4.72342], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(381.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH2]CCC=CC=O(1910)',
    structure = SMILES('[CH2]CCC=CC=O'),
    E0 = (41.4563,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,180,820.333],'cm^-1')),
        HinderedRotor(inertia=(0.174138,'amu*angstrom^2'), symmetry=1, barrier=(4.00377,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.513411,'amu*angstrom^2'), symmetry=1, barrier=(11.8043,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0247682,'amu*angstrom^2'), symmetry=1, barrier=(11.8121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.513542,'amu*angstrom^2'), symmetry=1, barrier=(11.8073,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43773,0.0582771,-3.53817e-05,9.92115e-09,-1.11161e-12,5075.3,26.3718], Tmin=(100,'K'), Tmax=(1963.5,'K')), NASAPolynomial(coeffs=[15.4265,0.029779,-1.36106e-05,2.52913e-09,-1.70418e-13,-418.039,-50.5546], Tmin=(1963.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(41.4563,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(RCCJ)"""),
)

species(
    label = '[CH]CCCC=CC=O(1911)',
    structure = SMILES('[CH]CCCC=CC=O'),
    E0 = (260.645,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2782.5,750,1395,475,1775,1000,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.154,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.85578,0.0736413,-5.27425e-05,1.91519e-08,-2.94465e-12,31457,29.9217], Tmin=(100,'K'), Tmax=(1431.01,'K')), NASAPolynomial(coeffs=[12.0754,0.0422797,-1.98688e-05,3.83695e-09,-2.69096e-13,28245.9,-28.2275], Tmin=(1431.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.645,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(CCJ2_triplet)"""),
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
    E0 = (162.046,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (128.469,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (317.09,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (383.383,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (423.926,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (439.355,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (366.942,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (423.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (467.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (420.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (430.401,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (390.09,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (38.7704,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (87.5488,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (165.371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (164.534,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (74.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (86.0739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (68.1246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (102.611,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (46.5456,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (89.1115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (327.252,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (216.817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (422.536,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (472.437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction22',
    reactants = ['H(19)', 'C=CCCC=CC=O(1889)'],
    products = ['[CH2]CCCC=CC=O(1870)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(3.01e+08,'cm^3/(mol*s)'), n=1.6, Ea=(10.0416,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 10 used for Cds-CsH_Cds-HH;HJ
Exact match found for rate rule [Cds-CsH_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C2H4(45)', '[CH2]CC=CC=O(1890)'],
    products = ['[CH2]CCCC=CC=O(1870)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(4240,'cm^3/(mol*s)'), n=2.41, Ea=(21.171,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 220 used for Cds-HH_Cds-HH;CsJ-CsHH
Exact match found for rate rule [Cds-HH_Cds-HH;CsJ-CsHH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C[CH2](171)', '[CH2]C=CC=O(1897)'],
    products = ['[CH2]CCCC=CC=O(1870)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 126 used for C_rad/H2/Cd;C_rad/H2/Cs
Exact match found for rate rule [C_rad/H2/Cd;C_rad/H2/Cs]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -0.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction31',
    reactants = ['C2H4(T)(169)', '[CH2]CC=CC=O(1890)'],
    products = ['[CH2]CCCC=CC=O(1870)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.94704e+07,'m^3/(mol*s)'), n=0.0484062, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [C_rad/H2/Cs;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -0.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction32',
    reactants = ['H(19)', '[CH2]C[CH]CC=CC=O(1898)'],
    products = ['[CH2]CCCC=CC=O(1870)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]CC[CH2](390)', 'CHCHCHO(1899)'],
    products = ['[CH2]CCCC=CC=O(1870)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.28206e+10,'m^3/(mol*s)'), n=-0.9855, Ea=(1.09644,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_rad;C_pri_rad] for rate rule [Cd_pri_rad;C_rad/H2/Cs]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['H(19)', '[CH2]CC[CH]C=CC=O(1900)'],
    products = ['[CH2]CCCC=CC=O(1870)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.71464e+07,'m^3/(mol*s)'), n=0.107721, Ea=(5.76381,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 36 used for C_rad/H/CdCs;H_rad
Exact match found for rate rule [C_rad/H/CdCs;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['H(19)', '[CH2][CH]CCC=CC=O(1901)'],
    products = ['[CH2]CCCC=CC=O(1870)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.34601e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction36',
    reactants = ['H(19)', '[CH2]CCC[C]=CC=O(1902)'],
    products = ['[CH2]CCCC=CC=O(1870)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['HCO(30)', '[CH]=CCCC[CH2](1903)'],
    products = ['[CH2]CCCC=CC=O(1870)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.81e+13,'cm^3/(mol*s)','*|/',3), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 91 used for Cd_pri_rad;CO_pri_rad
Exact match found for rate rule [CO_pri_rad;Cd_pri_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['H(19)', '[CH2]CCCC=[C]C=O(1904)'],
    products = ['[CH2]CCCC=CC=O(1870)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad/OneDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['H(19)', '[CH2]CCCC=C[C]=O(1905)'],
    products = ['[CH2]CCCC=CC=O(1870)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.53107e+07,'m^3/(mol*s)'), n=-0.133333, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;CO_sec_rad] for rate rule [H_rad;CO_rad/OneDe]
Euclidian distance = 1.41421356237
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]CCCC=CC=O(1870)'],
    products = ['O=C[CH]C1CCCC1(1888)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.3294e+07,'s^-1'), n=0.933333, Ea=(21.0943,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSS_D;doublebond_intra;radadd_intra_cs2H] for rate rule [R6_SSS_D;doublebond_intra_HDe_pri;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]CCCC=CC=O(1870)'],
    products = ['[O]C1C=CCCCC1(160)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.19e+11,'s^-1'), n=0.08, Ea=(69.8728,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R8;multiplebond_intra;radadd_intra_cs2H] for rate rule [R8;carbonylbond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]CCCC=CC=O(1870)'],
    products = ['C[CH]CCC=CC=O(1891)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(6.48e+07,'s^-1'), n=1.57, Ea=(147.695,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 106 used for R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]CCCC=CC=O(1870)'],
    products = ['CC[CH]CC=CC=O(1892)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.18e+10,'s^-1'), n=0.82, Ea=(146.858,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 186 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]CCCC=CC=O(1870)'],
    products = ['CCC[CH]C=CC=O(1893)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(42200,'s^-1'), n=1.93, Ea=(56.484,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 86 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/Cd
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CCCC[C]=CC=O(1894)'],
    products = ['[CH2]CCCC=CC=O(1870)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(561575,'s^-1'), n=1.6076, Ea=(35.8025,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;Y_rad_out;Cs_H_out_2H] for rate rule [R5H_CCC;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]CCCC=CC=O(1870)'],
    products = ['CCCCC=[C]C=O(1895)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(456.249,'s^-1'), n=2.07385, Ea=(50.4486,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_SSSSR;C_rad_out_2H;XH_out] for rate rule [R6H_SSSSD;C_rad_out_2H;Cd_H_out_singleDe]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]CCCC=CC=O(1870)'],
    products = ['CCCCC=C[C]=O(1896)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(342721,'s^-1'), n=1.49389, Ea=(84.9352,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7H;C_rad_out_2H;XH_out] for rate rule [R7H;C_rad_out_2H;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]CCCC=CC=O(1870)'],
    products = ['O=CC1[CH]CCCC1(1906)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(815229,'s^-1'), n=1.06667, Ea=(28.8696,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_CsCsCs_HH_D;doublebond_intra_pri;radadd_intra_cs2H] for rate rule [R6_CsCsCs_HH_D;doublebond_intra_pri_HCO;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]CCCC=CC=O(1870)'],
    products = ['[CH]1C=CCCCCO1(1907)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.18057e+11,'s^-1'), n=0.420859, Ea=(71.4354,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6plus;multiplebond_intra;radadd_intra_cs2H] + [R6plus;carbonyl_intra_H;radadd_intra] for rate rule [R8_linear;carbonyl_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.82842712475
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction42',
    reactants = ['CO(28)', '[CH2]CCCC=C(1908)'],
    products = ['[CH2]CCCC=CC=O(1870)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(0.118397,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO;R_H] for rate rule [CO;Cd_pri]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2]CCCC=C=CO(1909)'],
    products = ['[CH2]CCCC=CC=O(1870)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction44',
    reactants = ['CH2(T)(33)', '[CH2]CCC=CC=O(1910)'],
    products = ['[CH2]CCCC=CC=O(1870)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.23625e+06,'m^3/(mol*s)'), n=0.36814, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction45',
    reactants = ['H(19)', '[CH]CCCC=CC=O(1911)'],
    products = ['[CH2]CCCC=CC=O(1870)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '295',
    isomers = [
        '[CH2]CCCC=CC=O(1870)',
    ],
    reactants = [
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '295',
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

