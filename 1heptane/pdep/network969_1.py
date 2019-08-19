species(
    label = 'O[CH]OO[CH]O(1380)',
    structure = SMILES('O[CH]OO[CH]O'),
    E0 = (-153.587,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3580,3650,1210,1345,900,1100,3000,3050,390,425,1340,1360,335,370,473.778],'cm^-1')),
        HinderedRotor(inertia=(0.115055,'amu*angstrom^2'), symmetry=1, barrier=(18.3272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.79439,'amu*angstrom^2'), symmetry=1, barrier=(41.2566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.5257,'amu*angstrom^2'), symmetry=1, barrier=(12.0869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0758799,'amu*angstrom^2'), symmetry=1, barrier=(12.0869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.5257,'amu*angstrom^2'), symmetry=1, barrier=(12.0869,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0508,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.635543,0.0777041,-0.000119916,8.86333e-08,-2.47381e-11,-18354.5,24.4905], Tmin=(100,'K'), Tmax=(760.343,'K')), NASAPolynomial(coeffs=[14.404,0.0125359,-5.68477e-06,1.04139e-09,-6.974e-14,-20658.2,-39.543], Tmin=(760.343,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-153.587,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(212.019,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsHH) + group(Cs-OsOsHH) + radical(OCJO) + radical(OCJO)"""),
)

species(
    label = 'HOCHO(59)',
    structure = SMILES('O=CO'),
    E0 = (-389.211,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1051.41,1051.94,1052.55,2787.25,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00861396,'amu*angstrom^2'), symmetry=1, barrier=(6.76824,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4005.91,'J/mol'), sigma=(3.626,'angstroms'), dipoleMoment=(1.7,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.89836,-0.00355878,3.55205e-05,-4.385e-08,1.71078e-11,-46778.6,7.34954], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.61383,0.00644964,-2.29083e-06,3.6716e-10,-2.18737e-14,-45330.3,0.847884], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-389.211,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""HOCHO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = 'O=COO[CH]O(8560)',
    structure = SMILES('O=COO[CH]O'),
    E0 = (-392.956,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2782.5,750,1395,475,1775,1000,3615,1277.5,1000,3025,407.5,1350,352.5],'cm^-1')),
        HinderedRotor(inertia=(1.58276,'amu*angstrom^2'), symmetry=1, barrier=(36.3908,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.58521,'amu*angstrom^2'), symmetry=1, barrier=(36.4471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.5834,'amu*angstrom^2'), symmetry=1, barrier=(36.4056,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.58461,'amu*angstrom^2'), symmetry=1, barrier=(36.4332,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (91.0428,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47544,0.0438802,-1.8417e-05,-1.72648e-08,1.19734e-11,-47160.4,22.3051], Tmin=(100,'K'), Tmax=(971.474,'K')), NASAPolynomial(coeffs=[17.3511,0.00586513,-1.95341e-06,4.19608e-10,-3.58278e-14,-51535.6,-60.4695], Tmin=(971.474,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-392.956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(O2s-CsH) + group(Cs-OsOsHH) + group(Cds-OdOsH) + radical(OCJO)"""),
)

species(
    label = '[O]COO[CH]O(8561)',
    structure = SMILES('[O]COO[CH]O'),
    E0 = (-115.149,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,221.008,1295.36],'cm^-1')),
        HinderedRotor(inertia=(0.614634,'amu*angstrom^2'), symmetry=1, barrier=(21.341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.615036,'amu*angstrom^2'), symmetry=1, barrier=(21.3388,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.614871,'amu*angstrom^2'), symmetry=1, barrier=(21.342,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143865,'amu*angstrom^2'), symmetry=1, barrier=(39.4299,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0508,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26234,0.0582904,-6.32735e-05,3.36514e-08,-7.0253e-12,-13749,22.9199], Tmin=(100,'K'), Tmax=(1166.4,'K')), NASAPolynomial(coeffs=[14.091,0.0142958,-6.69549e-06,1.31335e-09,-9.40676e-14,-16741.7,-40.9455], Tmin=(1166.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-115.149,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsHH) + group(Cs-OsOsHH) + radical(OCJO) + radical(OCOJ)"""),
)

species(
    label = '[O][CH]OOCO(8562)',
    structure = SMILES('[O][CH]OOCO'),
    E0 = (-115.149,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,221.008,1295.36],'cm^-1')),
        HinderedRotor(inertia=(0.614634,'amu*angstrom^2'), symmetry=1, barrier=(21.341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.615036,'amu*angstrom^2'), symmetry=1, barrier=(21.3388,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.614871,'amu*angstrom^2'), symmetry=1, barrier=(21.342,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143865,'amu*angstrom^2'), symmetry=1, barrier=(39.4299,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0508,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26234,0.0582904,-6.32735e-05,3.36514e-08,-7.0253e-12,-13749,22.9199], Tmin=(100,'K'), Tmax=(1166.4,'K')), NASAPolynomial(coeffs=[14.091,0.0142958,-6.69549e-06,1.31335e-09,-9.40676e-14,-16741.7,-40.9455], Tmin=(1166.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-115.149,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsHH) + group(Cs-OsOsHH) + radical(OCOJ) + radical(OCJO)"""),
)

species(
    label = '[O][CH]O(171)',
    structure = SMILES('[O][CH]O'),
    E0 = (5.37818,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,1989.68],'cm^-1')),
        HinderedRotor(inertia=(0.0937094,'amu*angstrom^2'), symmetry=1, barrier=(2.15456,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.40204,0.0209571,-4.96202e-05,5.95355e-08,-2.44634e-11,660.882,10.3084], Tmin=(100,'K'), Tmax=(868.739,'K')), NASAPolynomial(coeffs=[-0.343229,0.0179325,-9.40003e-06,1.81361e-09,-1.23835e-13,2076.48,32.2523], Tmin=(868.739,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.37818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsHH) + radical(OCJO) + radical(OCOJ)"""),
)

species(
    label = 'O=COOCO(8563)',
    structure = SMILES('O=COOCO'),
    E0 = (-581.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0508,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53671,0.027387,6.26027e-05,-1.21627e-07,5.3974e-11,-69874.9,22.2953], Tmin=(100,'K'), Tmax=(929.05,'K')), NASAPolynomial(coeffs=[25.3531,-0.00628539,5.77692e-06,-1.06139e-09,6.01961e-14,-77272.4,-106.848], Tmin=(929.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-581.917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(O2s-CsH) + group(Cs-OsOsHH) + group(Cds-OdOsH)"""),
)

species(
    label = 'OC1OOC1O(1381)',
    structure = SMILES('OC1OOC1O'),
    E0 = (-416.173,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0508,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.18325,0.0415348,-3.52002e-05,1.56092e-08,-2.89761e-12,-49990.1,16.5333], Tmin=(100,'K'), Tmax=(1241.62,'K')), NASAPolynomial(coeffs=[8.75232,0.0203719,-9.63326e-06,1.88142e-09,-1.33517e-13,-51621.4,-16.5803], Tmin=(1241.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-416.173,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + ring(12dioxetane)"""),
)

species(
    label = 'HCOH(T)(285)',
    structure = SMILES('[CH]O'),
    E0 = (205.906,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,403.876,3308.82],'cm^-1')),
        HinderedRotor(inertia=(0.0103144,'amu*angstrom^2'), symmetry=1, barrier=(22.7121,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.75938,0.0029613,8.90411e-06,-1.35016e-08,5.39816e-12,24775.6,6.76286], Tmin=(100,'K'), Tmax=(940.429,'K')), NASAPolynomial(coeffs=[5.09112,0.00321239,-9.31686e-07,1.59615e-10,-1.15729e-14,24263.5,-0.971], Tmin=(940.429,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(78.9875,'J/(mol*K)'), label="""HCOH(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[O]O[CH]O(4538)',
    structure = SMILES('[O]O[CH]O'),
    E0 = (9.28844,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,492.5,1135,1000],'cm^-1')),
        HinderedRotor(inertia=(0.156407,'amu*angstrom^2'), symmetry=1, barrier=(3.59612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156804,'amu*angstrom^2'), symmetry=1, barrier=(3.60523,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (62.0248,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.58857,0.0452524,-0.00011197,1.22556e-07,-4.65581e-11,1154.51,14.2974], Tmin=(100,'K'), Tmax=(891.504,'K')), NASAPolynomial(coeffs=[-0.0459,0.0212294,-1.12418e-05,2.13324e-09,-1.41938e-13,3048.61,34.6933], Tmin=(891.504,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(9.28844,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(OCJO)"""),
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
    E0 = (-153.587,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-145.278,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-1.17729,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-50.4228,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (10.7564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-145.219,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-145.303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (215.195,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O[CH]OO[CH]O(1380)'],
    products = ['HOCHO(59)', 'HOCHO(59)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(3)', 'O=COO[CH]O(8560)'],
    products = ['O[CH]OO[CH]O(1380)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(4e+09,'cm^3/(mol*s)'), n=1.39, Ea=(35.8862,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [Od_CO-NdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]COO[CH]O(8561)'],
    products = ['O[CH]OO[CH]O(1380)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.3e+14,'s^-1','+|-',2), n=-0.27, Ea=(113.972,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1800,'K'), comment="""Estimated using template [R2H_S;O_rad_out;Cs_H_out] for rate rule [R2H_S;O_rad_out;Cs_H_out_OOH/H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O][CH]OOCO(8562)'],
    products = ['O[CH]OO[CH]O(1380)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(23000,'s^-1'), n=2.11, Ea=(64.7265,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5Hall;O_rad_out;Cs_H_out_H/NonDeO] for rate rule [R5HJ_1;O_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O][CH]O(171)', '[O][CH]O(171)'],
    products = ['O[CH]OO[CH]O(1380)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.73038e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['O[CH]OO[CH]O(1380)'],
    products = ['O=COOCO(8563)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(6.42e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O[CH]OO[CH]O(1380)'],
    products = ['OC1OOC1O(1381)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/NonDeO;Cpri_rad_out_H/NonDeO]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['HCOH(T)(285)', '[O]O[CH]O(4538)'],
    products = ['O[CH]OO[CH]O(1380)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

network(
    label = '969',
    isomers = [
        'O[CH]OO[CH]O(1380)',
    ],
    reactants = [
        ('HOCHO(59)', 'HOCHO(59)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '969',
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

