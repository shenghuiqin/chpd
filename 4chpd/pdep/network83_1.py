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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.35532,0.0344795,-1.7272e-05,4.08056e-09,-3.78294e-13,32109.4,17.5603], Tmin=(100,'K'), Tmax=(2451.48,'K')), NASAPolynomial(coeffs=[13.9856,0.0155026,-5.66049e-06,9.22849e-10,-5.6271e-14,26407.1,-48.9781], Tmin=(2451.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(266.469,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJ) + radical(RCCJ)"""),
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
    label = 'buten3yl1(651)',
    structure = SMILES('[CH2]CC=C'),
    E0 = (191.701,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,193.365],'cm^-1')),
        HinderedRotor(inertia=(0.286961,'amu*angstrom^2'), symmetry=1, barrier=(7.61231,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00450981,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.48601,0.0277403,-7.5101e-07,-1.39971e-08,6.14188e-12,23115.6,15.6915], Tmin=(100,'K'), Tmax=(1051,'K')), NASAPolynomial(coeffs=[7.36564,0.0210619,-8.1931e-06,1.49014e-09,-1.03098e-13,21433,-11.2175], Tmin=(1051,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), label="""buten3yl1""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2]C[CH]C(652)',
    structure = SMILES('[CH2]C[CH]C'),
    E0 = (255.669,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.98997,0.0287412,-9.51469e-06,4.19232e-10,1.90526e-13,30780.1,16.8971], Tmin=(100,'K'), Tmax=(2154.56,'K')), NASAPolynomial(coeffs=[12.4231,0.0182241,-7.06316e-06,1.16769e-09,-7.11818e-14,25091.4,-39.6212], Tmin=(2154.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(255.669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CC(653)',
    structure = SMILES('[CH2][CH]CC'),
    E0 = (255.669,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.98997,0.0287412,-9.51469e-06,4.19232e-10,1.90526e-13,30780.1,16.8971], Tmin=(100,'K'), Tmax=(2154.56,'K')), NASAPolynomial(coeffs=[12.4231,0.0182241,-7.06316e-06,1.16769e-09,-7.11818e-14,25091.4,-39.6212], Tmin=(2154.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(255.669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'butene1(654)',
    structure = SMILES('C=CCC'),
    E0 = (-16.4325,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.58773,0.0232778,1.93412e-05,-3.55496e-08,1.36906e-11,-1918.73,14.5751], Tmin=(100,'K'), Tmax=(1007.28,'K')), NASAPolynomial(coeffs=[7.20517,0.0236362,-9.0315e-06,1.65393e-09,-1.16019e-13,-3797.34,-12.4426], Tmin=(1007.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-16.4325,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), label="""butene1""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C1CCC1(391)',
    structure = SMILES('C1CCC1'),
    E0 = (13.873,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.2702,-0.00230965,0.00010246,-1.25882e-07,4.64171e-11,1711.44,10.4049], Tmin=(100,'K'), Tmax=(958.659,'K')), NASAPolynomial(coeffs=[9.08685,0.0204702,-6.80115e-06,1.29477e-09,-9.89856e-14,-1565.8,-28.6877], Tmin=(958.659,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(13.873,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + ring(Cyclobutane)"""),
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
    E0 = (266.469,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (413.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (381.666,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (418.348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (384.408,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (636.292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (329.869,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (274.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (671.578,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]CC[CH2](390)'],
    products = ['C2H4(45)', 'C2H4(45)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(19)', 'buten3yl1(651)'],
    products = ['[CH2]CC[CH2](390)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.01e+08,'cm^3/(mol*s)'), n=1.6, Ea=(10.0416,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 10 used for Cds-CsH_Cds-HH;HJ
Exact match found for rate rule [Cds-CsH_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C2H4(T)(169)', 'C2H4(45)'],
    products = ['[CH2]CC[CH2](390)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(0.00620033,'m^3/(mol*s)'), n=2.49247, Ea=(21.458,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-HH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]CC[CH2](390)'],
    products = ['[CH2]C[CH]C(652)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.436e+06,'s^-1'), n=2.05, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(500,'K'), Tmax=(2000,'K'), comment="""From training reaction 147 used for R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]CC[CH2](390)'],
    products = ['[CH2][CH]CC(653)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(333380,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C2H4(T)(169)', 'C2H4(T)(169)'],
    products = ['[CH2]CC[CH2](390)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.47003e+07,'m^3/(mol*s)'), n=0.0284742, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]CC[CH2](390)'],
    products = ['butene1(654)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.9748e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]CC[CH2](390)'],
    products = ['C1CCC1(391)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Matched reaction 2 C4H8 <=> C4H8-2 in Birad_recombination/training
This reaction matched rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C[CH2](171)', 'CH2(T)(33)'],
    products = ['[CH2]CC[CH2](390)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.4725e+06,'m^3/(mol*s)'), n=0.36814, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

network(
    label = '83',
    isomers = [
        '[CH2]CC[CH2](390)',
    ],
    reactants = [
        ('C2H4(45)', 'C2H4(45)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '83',
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

