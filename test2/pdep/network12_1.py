species(
    label = 'C1=CC2CCCC12(59)',
    structure = SMILES('C1=CC2CCCC12'),
    E0 = (108.528,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3725.22,'J/mol'), sigma=(6.44341,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=581.87 K, Pc=31.6 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.30701,0.00901903,0.000132148,-1.74471e-07,6.6169e-11,13139.4,18.0522], Tmin=(100,'K'), Tmax=(958.829,'K')), NASAPolynomial(coeffs=[15.1122,0.0274004,-8.93448e-06,1.7229e-09,-1.34172e-13,7383.31,-60.3982], Tmin=(958.829,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(108.528,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
)

species(
    label = '[CH]1C[C]2CCCC12(73)',
    structure = SMILES('[CH]1C[C]2CCCC12'),
    E0 = (338.717,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,303.444,805.837,805.837,805.837,805.837,805.837,805.837,805.837,805.837,805.837,805.837,805.837,805.837,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.58889,0.012735,9.48241e-05,-1.16716e-07,4.09595e-11,40805.1,22.7077], Tmin=(100,'K'), Tmax=(1010.55,'K')), NASAPolynomial(coeffs=[7.78141,0.039824,-1.61028e-05,3.11628e-09,-2.27358e-13,37323,-14.4339], Tmin=(1010.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.717,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(bicyclo[3.2.0]heptane-tertiary) + radical(bicyclo[3.2.0]heptane-C5-6)"""),
)

species(
    label = '[CH]1CC2CCC[C]12(74)',
    structure = SMILES('[CH]1CC2CCC[C]12'),
    E0 = (338.717,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,303.444,805.837,805.837,805.837,805.837,805.837,805.837,805.837,805.837,805.837,805.837,805.837,805.837,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57,1609.57],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.58891,0.0127347,9.48252e-05,-1.16717e-07,4.09601e-11,40805.1,22.7076], Tmin=(100,'K'), Tmax=(1010.54,'K')), NASAPolynomial(coeffs=[7.7813,0.0398242,-1.61029e-05,3.11631e-09,-2.2736e-13,37323,-14.4332], Tmin=(1010.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.717,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(bicyclo[3.2.0]heptane-C5-6) + radical(bicyclo[3.2.0]heptane-tertiary)"""),
)

species(
    label = '[CH]1CCC2[CH]CC12(75)',
    structure = SMILES('[CH]1CCC2[CH]CC12'),
    E0 = (319.794,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66775,0.00592624,0.000123203,-1.51785e-07,5.46311e-11,38531.1,22.9424], Tmin=(100,'K'), Tmax=(986.716,'K')), NASAPolynomial(coeffs=[10.1915,0.0357275,-1.37706e-05,2.6963e-09,-2.01715e-13,34110.9,-28.1299], Tmin=(986.716,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(319.794,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(bicyclo[3.2.0]heptane-C5-6) + radical(bicyclo[3.2.0]heptane-C5-2)"""),
)

species(
    label = '[CH]1CCC2C[CH]C12(76)',
    structure = SMILES('[CH]1CCC2C[CH]C12'),
    E0 = (319.794,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66775,0.00592624,0.000123203,-1.51785e-07,5.46311e-11,38531.1,22.9424], Tmin=(100,'K'), Tmax=(986.716,'K')), NASAPolynomial(coeffs=[10.1915,0.0357275,-1.37706e-05,2.6963e-09,-2.01715e-13,34110.9,-28.1299], Tmin=(986.716,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(319.794,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(bicyclo[3.2.0]heptane-C5-6) + radical(bicyclo[3.2.0]heptane-C5-2)"""),
)

species(
    label = '[CH]1CC2[CH]CC2C1(77)',
    structure = SMILES('[CH]1CC2[CH]CC2C1'),
    E0 = (326.489,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66775,0.00592624,0.000123203,-1.51785e-07,5.46311e-11,39336.3,22.9424], Tmin=(100,'K'), Tmax=(986.716,'K')), NASAPolynomial(coeffs=[10.1915,0.0357275,-1.37706e-05,2.6963e-09,-2.01715e-13,34916,-28.1299], Tmin=(986.716,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(326.489,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(bicyclo[3.2.0]heptane-C5-6) + radical(bicyclo[3.2.0]heptane-C5-3)"""),
)

species(
    label = '[CH]=CC1[CH]CCC1(78)',
    structure = SMILES('[CH]=CC1[CH]CCC1'),
    E0 = (403.14,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,3120,650,792.5,1650,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88819,0.0268691,7.11195e-05,-1.03575e-07,3.90066e-11,48579.9,25.2513], Tmin=(100,'K'), Tmax=(993.488,'K')), NASAPolynomial(coeffs=[12.4137,0.033516,-1.29354e-05,2.4985e-09,-1.84363e-13,44069.1,-37.6356], Tmin=(993.488,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(403.14,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cs_S) + radical(Cds_P)"""),
)

species(
    label = '[CH2]CC1C=CC1[CH2](79)',
    structure = SMILES('[CH2]CC1C=CC1[CH2]'),
    E0 = (468.452,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,1082.03,1082.03,1082.03,1082.03,1082.03,1082.03,1082.03,1082.03,1082.03,1082.03,1082.03,1082.03,1082.03,2231.58],'cm^-1')),
        HinderedRotor(inertia=(0.0951468,'amu*angstrom^2'), symmetry=1, barrier=(2.18761,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0951468,'amu*angstrom^2'), symmetry=1, barrier=(2.18761,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0951468,'amu*angstrom^2'), symmetry=1, barrier=(2.18761,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44248,0.0408177,3.32491e-05,-6.77816e-08,2.77761e-11,56447.6,26.1745], Tmin=(100,'K'), Tmax=(977.723,'K')), NASAPolynomial(coeffs=[12.6761,0.0324508,-1.15862e-05,2.11346e-09,-1.50784e-13,52454.1,-36.9565], Tmin=(977.723,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(468.452,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CCC1[CH]C=C1(80)',
    structure = SMILES('[CH2]CCC1[CH]C=C1'),
    E0 = (434.123,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,363.369,1018.96,1018.96,1018.96,1018.96,1018.96,1018.96,1018.96,1018.96,1018.96,1018.96,1018.96,1018.96,2256.16],'cm^-1')),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00975557,'amu*angstrom^2'), symmetry=1, barrier=(0.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79043,0.0277409,7.43956e-05,-1.1054e-07,4.23935e-11,52311,22.78], Tmin=(100,'K'), Tmax=(981.018,'K')), NASAPolynomial(coeffs=[13.5491,0.0319656,-1.18326e-05,2.26558e-09,-1.67938e-13,47493.6,-46.5185], Tmin=(981.018,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(434.123,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(cyclobutene-allyl) + radical(RCCJ)"""),
)

species(
    label = '[C]1CC2CCCC12(81)',
    structure = SMILES('[C]1CC2CCCC12'),
    E0 = (337.886,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2794.44,2838.89,2883.33,2927.78,2972.22,3016.67,3061.11,3105.56,3150,900,922.222,944.444,966.667,988.889,1011.11,1033.33,1055.56,1077.78,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2288,0.0110772,0.000127348,-1.69745e-07,6.44766e-11,40727.3,19.2909], Tmin=(100,'K'), Tmax=(960.68,'K')), NASAPolynomial(coeffs=[15.2612,0.027668,-9.18817e-06,1.77724e-09,-1.37976e-13,34953.7,-60.0772], Tmin=(960.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(337.886,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""QM MopacMolPM3 calculation attempt 1"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.61263,-0.00100893,2.49899e-06,-1.43376e-09,2.58637e-13,-1051.1,2.6527], Tmin=(100,'K'), Tmax=(1817.03,'K')), NASAPolynomial(coeffs=[2.97589,0.00164143,-7.1973e-07,1.25379e-10,-7.91538e-15,-1025.84,5.53764], Tmin=(1817.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.69489,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""N2""", comment="""Thermo library: BurkeH2O2"""),
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
    E0 = (361.863,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (402.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (383.195,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (328.162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (334.857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (411.048,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (475.564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (441.445,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (429.22,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction39',
    reactants = ['[CH]1C[C]2CCCC12(73)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5.10299e+10,'s^-1'), n=0.2847, Ea=(23.1459,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad_NDe] + [R2radExo;Y_rad;XH_Rrad_NDe] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]1CC2CCC[C]12(74)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH]1CCC2[CH]CC12(75)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 1 used for R3radExo;Y_rad_NDe;XH_Rrad_NDe
Exact match found for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH]1CCC2C[CH]C12(76)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH]1CC2[CH]CC2C1(77)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad_NDe] for rate rule [R4radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH]=CC1[CH]CCC1(78)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/NonDeC;CdsinglepriH_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]CC1C=CC1[CH2](79)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for R5_SSSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R5_SSSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2]CCC1[CH]C=C1(80)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.49159e+11,'s^-1'), n=0.1555, Ea=(7.322,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Cpri_rad_out_2H] + [R5_SSSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R5_SSSS;C_rad_out_H/OneDe;Cpri_rad_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[C]1CC2CCCC12(81)'],
    products = ['C1=CC2CCCC12(59)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.83659e+12,'s^-1'), n=0.345437, Ea=(91.3338,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [CsJ2-C;CsJ2(CsC);CH] for rate rule [CsJ2-C;CsJ2(CsC);CH2(C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

network(
    label = '12',
    isomers = [
        'C1=CC2CCCC12(59)',
    ],
    reactants = [
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '12',
    Tmin = (300,'K'),
    Tmax = (3000,'K'),
    Tcount = 8,
    Tlist = ([302.617,324.619,374.997,470.374,649.057,1000.02,1706.11,2761.25],'K'),
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

