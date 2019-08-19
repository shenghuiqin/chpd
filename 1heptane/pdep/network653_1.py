species(
    label = 'C=[C]OO[C]=C(1076)',
    structure = SMILES('C=[C]OO[C]=C'),
    E0 = (542.634,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,1670,1700,300,440,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180],'cm^-1')),
        HinderedRotor(inertia=(0.150286,'amu*angstrom^2'), symmetry=1, barrier=(3.45536,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149149,'amu*angstrom^2'), symmetry=1, barrier=(3.42922,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.925521,'amu*angstrom^2'), symmetry=1, barrier=(21.2796,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.24651,0.0434753,-4.27674e-05,1.75398e-08,2.32461e-12,65322.3,26.6208], Tmin=(100,'K'), Tmax=(567.426,'K')), NASAPolynomial(coeffs=[5.79377,0.0251761,-1.21227e-05,2.36572e-09,-1.67388e-13,64811.7,10.566], Tmin=(567.426,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(542.634,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = 'CH2CO(28)',
    structure = SMILES('C=C=O'),
    E0 = (-60.8183,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2120,512.5,787.5],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13241,0.0181319,-1.74093e-05,9.35336e-09,-2.01725e-12,-7148.09,13.3808], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.75871,0.00635124,-2.25955e-06,3.62322e-10,-2.15856e-14,-8085.33,-4.9649], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-60.8183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CH2CO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = 'C#COO[C]=C(4570)',
    structure = SMILES('C#COO[C]=C'),
    E0 = (534.406,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2175,525,1685,370,750,770,3400,2100,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(1.34534,'amu*angstrom^2'), symmetry=1, barrier=(30.9319,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34576,'amu*angstrom^2'), symmetry=1, barrier=(30.9416,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34545,'amu*angstrom^2'), symmetry=1, barrier=(30.9345,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88133,0.0503015,-6.29738e-05,4.33222e-08,-1.22468e-11,64347.3,22.3829], Tmin=(100,'K'), Tmax=(853.695,'K')), NASAPolynomial(coeffs=[8.31092,0.0201755,-1.00404e-05,1.98546e-09,-1.41559e-13,63249.5,-7.61905], Tmin=(853.695,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(534.406,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsCt) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(C=CJO)"""),
)

species(
    label = '[CH]=COO[C]=C(4571)',
    structure = SMILES('[CH]=COO[C]=C'),
    E0 = (549.986,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,1685,370,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.918356,'amu*angstrom^2'), symmetry=1, barrier=(21.1148,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.720452,'amu*angstrom^2'), symmetry=1, barrier=(16.5646,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.720819,'amu*angstrom^2'), symmetry=1, barrier=(16.5731,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73695,0.0485798,-4.82298e-05,2.46636e-08,-5.05784e-12,66230.5,25.6851], Tmin=(100,'K'), Tmax=(1173.77,'K')), NASAPolynomial(coeffs=[11.1557,0.0164823,-7.21151e-06,1.36645e-09,-9.58106e-14,64019.4,-21.2641], Tmin=(1173.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(549.986,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]OOC=C(4572)',
    structure = SMILES('[CH]=[C]OOC=C'),
    E0 = (549.986,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,1685,370,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.918356,'amu*angstrom^2'), symmetry=1, barrier=(21.1148,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.720452,'amu*angstrom^2'), symmetry=1, barrier=(16.5646,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.720819,'amu*angstrom^2'), symmetry=1, barrier=(16.5731,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73695,0.0485798,-4.82298e-05,2.46636e-08,-5.05784e-12,66230.5,25.6851], Tmin=(100,'K'), Tmax=(1173.77,'K')), NASAPolynomial(coeffs=[11.1557,0.0164823,-7.21151e-06,1.36645e-09,-9.58106e-14,64019.4,-21.2641], Tmin=(1173.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(549.986,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CJO)"""),
)

species(
    label = 'C=[C][O](173)',
    structure = SMILES('[CH2][C]=O'),
    E0 = (160.185,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,539.612,539.669],'cm^-1')),
        HinderedRotor(inertia=(0.000578908,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39563,0.0101365,2.30741e-06,-8.97566e-09,3.68242e-12,19290.3,10.0703], Tmin=(100,'K'), Tmax=(1068.9,'K')), NASAPolynomial(coeffs=[6.35055,0.00638951,-2.69368e-06,5.4221e-10,-4.02476e-14,18240.9,-6.33602], Tmin=(1068.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CsCJ=O) + radical(CJC=O)"""),
)

species(
    label = 'C=C1OOC1=C(4573)',
    structure = SMILES('C=C1OOC1=C'),
    E0 = (178.533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26274,0.0337843,-1.66016e-05,1.99158e-09,3.87466e-13,21538.7,15.006], Tmin=(100,'K'), Tmax=(1462.36,'K')), NASAPolynomial(coeffs=[10.3315,0.0187829,-8.46519e-06,1.58796e-09,-1.08648e-13,18423,-29.5719], Tmin=(1462.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(178.533,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(12methylenecyclobutane)"""),
)

species(
    label = 'C=[C]O[O](110)',
    structure = SMILES('C=[C]O[O]'),
    E0 = (339.001,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,492.5,1135,1000,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.0942477,'amu*angstrom^2'), symmetry=1, barrier=(2.16694,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.10787,0.026676,-5.52143e-05,6.05409e-08,-2.36332e-11,40797.8,12.0288], Tmin=(100,'K'), Tmax=(872.797,'K')), NASAPolynomial(coeffs=[1.34334,0.0171001,-8.40162e-06,1.59781e-09,-1.08413e-13,41778.6,24.1556], Tmin=(872.797,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(339.001,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(ROOJ) + radical(C=CJO)"""),
)

species(
    label = 'H2CC(41)',
    structure = SMILES('[C]=C'),
    E0 = (401.202,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2480.69,'J/mol'), sigma=(4.48499,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=387.48 K, Pc=62.39 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28155,0.00697643,-2.38528e-06,-1.21078e-09,9.82042e-13,48319.2,5.92036], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.27807,0.00475623,-1.63007e-06,2.54623e-10,-1.4886e-14,48014,0.639979], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(401.202,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""H2CC""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    E0 = (542.634,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (761.308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (655.423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (859.072,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (542.634,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (550.918,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (740.203,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C]OO[C]=C(1076)'],
    products = ['CH2CO(28)', 'CH2CO(28)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(3)', 'C#COO[C]=C(4570)'],
    products = ['C=[C]OO[C]=C(1076)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=COO[C]=C(4571)'],
    products = ['C=[C]OO[C]=C(1076)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=[C]OOC=C(4572)'],
    products = ['C=[C]OO[C]=C(1076)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.24929e+07,'s^-1'), n=1.719, Ea=(309.087,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_doubleC] for rate rule [R5HJ_1;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=[C][O](173)', 'C=[C][O](173)'],
    products = ['C=[C]OO[C]=C(1076)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.73038e+06,'m^3/(mol*s)'), n=0.027223, Ea=(222.264,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 222.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=[C]OO[C]=C(1076)'],
    products = ['C=C1OOC1=C(4573)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=[C]O[O](110)', 'H2CC(41)'],
    products = ['C=[C]OO[C]=C(1076)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

network(
    label = '653',
    isomers = [
        'C=[C]OO[C]=C(1076)',
    ],
    reactants = [
        ('CH2CO(28)', 'CH2CO(28)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '653',
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

