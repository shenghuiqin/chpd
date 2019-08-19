species(
    label = 'C=C[CH]CC=[C]O(23340)',
    structure = SMILES('[CH2]C=CCC=[C]O'),
    E0 = (235.175,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,3615,1277.5,1000,3000,3100,440,815,1455,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,312.969,313.311],'cm^-1')),
        HinderedRotor(inertia=(0.248168,'amu*angstrom^2'), symmetry=1, barrier=(17.2615,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.247382,'amu*angstrom^2'), symmetry=1, barrier=(17.2633,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.248086,'amu*angstrom^2'), symmetry=1, barrier=(17.2626,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.247891,'amu*angstrom^2'), symmetry=1, barrier=(17.2646,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.741639,0.0589163,-2.04542e-05,-2.2798e-08,1.5114e-11,28413.9,29.5929], Tmin=(100,'K'), Tmax=(952.355,'K')), NASAPolynomial(coeffs=[17.6548,0.0187554,-5.83105e-06,1.00869e-09,-7.19926e-14,23792.2,-58.5292], Tmin=(952.355,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(235.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(Allyl_P)"""),
)

species(
    label = 'HCCOH(50)',
    structure = SMILES('C#CO'),
    E0 = (80.0402,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3615,1277.5,1000,287.694,287.7,1588.42],'cm^-1')),
        HinderedRotor(inertia=(0.002037,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05541,0.0252003,-3.80822e-05,3.09891e-08,-9.898e-12,9768.72,12.2272], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[6.3751,0.00549429,-1.88137e-06,2.93804e-10,-1.71772e-14,8932.78,-8.24498], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(80.0402,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""HCCOH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = '[CH2][CH]C1CC=C1O(23610)',
    structure = SMILES('[CH2][CH]C1CC=C1O'),
    E0 = (274.165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22747,0.0429655,2.53576e-05,-6.95849e-08,3.12889e-11,33090.6,27.202], Tmin=(100,'K'), Tmax=(953.919,'K')), NASAPolynomial(coeffs=[17.8893,0.0177898,-5.32968e-06,9.74787e-10,-7.45876e-14,27878.5,-63.0536], Tmin=(953.919,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.165,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(Cs_S) + radical(RCCJ)"""),
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
    label = '[CH2]C=CCC=C=O(20121)',
    structure = SMILES('[CH2]C=CCC=C=O'),
    E0 = (117.74,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.601726,'amu*angstrom^2'), symmetry=1, barrier=(13.8349,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0466387,'amu*angstrom^2'), symmetry=1, barrier=(9.36354,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02553,'amu*angstrom^2'), symmetry=1, barrier=(23.5789,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23481,0.059716,-4.88266e-05,2.14387e-08,-3.91704e-12,14261.1,24.0972], Tmin=(100,'K'), Tmax=(1275.05,'K')), NASAPolynomial(coeffs=[11.1768,0.0285271,-1.21357e-05,2.25488e-09,-1.557e-13,11725.8,-26.2829], Tmin=(1275.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.74,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-Cd(CCO)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + radical(Allyl_P)"""),
)

species(
    label = 'C=C=CCC=[C]O(23611)',
    structure = SMILES('C=C=CCC=[C]O'),
    E0 = (260.28,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,1685,370,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.827826,'amu*angstrom^2'), symmetry=1, barrier=(19.0334,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.827866,'amu*angstrom^2'), symmetry=1, barrier=(19.0343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.83732,'amu*angstrom^2'), symmetry=1, barrier=(19.2516,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.778648,0.0610877,-3.8194e-05,-7.79675e-10,6.65864e-12,31429.2,28.233], Tmin=(100,'K'), Tmax=(963.751,'K')), NASAPolynomial(coeffs=[16.58,0.0179254,-5.91127e-06,1.02793e-09,-7.19693e-14,27342.3,-52.8182], Tmin=(963.751,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.28,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C=CCC#CO(23612)',
    structure = SMILES('C=C[CH]CC#CO'),
    E0 = (229.45,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2100,2250,500,550,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,180,180,792.268],'cm^-1')),
        HinderedRotor(inertia=(0.0385837,'amu*angstrom^2'), symmetry=1, barrier=(2.70566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0381413,'amu*angstrom^2'), symmetry=1, barrier=(16.9829,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.199078,'amu*angstrom^2'), symmetry=1, barrier=(88.6222,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59727,'amu*angstrom^2'), symmetry=1, barrier=(36.7243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29754,0.0515312,-2.59394e-05,-9.11102e-10,3.50414e-12,27700.5,24.8163], Tmin=(100,'K'), Tmax=(1061.75,'K')), NASAPolynomial(coeffs=[11.8524,0.0260537,-1.01296e-05,1.83547e-09,-1.26685e-13,24653.9,-30.5296], Tmin=(1061.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(229.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CtH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtOs) + radical(Allyl_S)"""),
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
    label = '[CH2]C=CC[CH]C=O(20123)',
    structure = SMILES('[CH2]C=CCC=C[O]'),
    E0 = (136.893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3000,3100,440,815,1455,1000,288.83,288.868,289.433],'cm^-1')),
        HinderedRotor(inertia=(0.382724,'amu*angstrom^2'), symmetry=1, barrier=(22.4549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.376442,'amu*angstrom^2'), symmetry=1, barrier=(22.4702,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.379495,'amu*angstrom^2'), symmetry=1, barrier=(22.4623,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.949875,0.0501072,1.01307e-05,-5.57741e-08,2.67736e-11,16589.7,27.1272], Tmin=(100,'K'), Tmax=(955.034,'K')), NASAPolynomial(coeffs=[18.4853,0.0184131,-5.66331e-06,1.02492e-09,-7.696e-14,11336.3,-66.6323], Tmin=(955.034,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(136.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = 'C[C]=CCC=[C]O(23613)',
    structure = SMILES('C[C]=CCC=[C]O'),
    E0 = (321.517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.526236,0.0680975,-6.07589e-05,2.83443e-08,-5.25795e-12,38801.7,30.5631], Tmin=(100,'K'), Tmax=(1306.12,'K')), NASAPolynomial(coeffs=[15.5725,0.0220181,-7.83935e-06,1.33311e-09,-8.78251e-14,34871.2,-46.0445], Tmin=(1306.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(321.517,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=CC[C]=CO(23614)',
    structure = SMILES('[CH2]C=CC[C]=CO'),
    E0 = (233.272,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,3615,1277.5,1000,3000,3100,440,815,1455,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.894772,'amu*angstrom^2'), symmetry=1, barrier=(20.5726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.895795,'amu*angstrom^2'), symmetry=1, barrier=(20.5961,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.894823,'amu*angstrom^2'), symmetry=1, barrier=(20.5737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.895314,'amu*angstrom^2'), symmetry=1, barrier=(20.585,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.50828,0.0616464,-1.79746e-05,-3.27952e-08,2.05646e-11,28195.8,27.7757], Tmin=(100,'K'), Tmax=(935.225,'K')), NASAPolynomial(coeffs=[20.5413,0.014238,-3.32393e-06,5.2059e-10,-3.88035e-14,22775,-76.4793], Tmin=(935.225,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(233.272,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = 'CC=[C]CC=[C]O(23615)',
    structure = SMILES('CC=[C]CC=[C]O'),
    E0 = (321.517,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,1670,1700,300,440,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,191.912,193.426],'cm^-1')),
        HinderedRotor(inertia=(0.465503,'amu*angstrom^2'), symmetry=1, barrier=(12.5844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.461016,'amu*angstrom^2'), symmetry=1, barrier=(12.582,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.462048,'amu*angstrom^2'), symmetry=1, barrier=(12.5768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.472012,'amu*angstrom^2'), symmetry=1, barrier=(12.5972,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.526236,0.0680975,-6.07589e-05,2.83443e-08,-5.25795e-12,38801.7,30.5631], Tmin=(100,'K'), Tmax=(1306.12,'K')), NASAPolynomial(coeffs=[15.5725,0.0220181,-7.83935e-06,1.33311e-09,-8.78251e-14,34871.2,-46.0445], Tmin=(1306.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(321.517,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C=C[CH]C=CO(23616)',
    structure = SMILES('[CH2]C=C[CH]C=CO'),
    E0 = (97.1557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.861405,0.0446295,4.46072e-05,-1.07119e-07,4.95686e-11,11820.8,25.7058], Tmin=(100,'K'), Tmax=(918.405,'K')), NASAPolynomial(coeffs=[23.9367,0.00763448,1.30669e-06,-3.95217e-10,2.17645e-14,4904.04,-98.2363], Tmin=(918.405,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(97.1557,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CCJC=C) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=[C]CC=CO(23617)',
    structure = SMILES('[CH2]C=[C]CC=CO'),
    E0 = (233.272,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,3615,1277.5,1000,3000,3100,440,815,1455,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.894772,'amu*angstrom^2'), symmetry=1, barrier=(20.5726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.895795,'amu*angstrom^2'), symmetry=1, barrier=(20.5961,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.894823,'amu*angstrom^2'), symmetry=1, barrier=(20.5737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.895314,'amu*angstrom^2'), symmetry=1, barrier=(20.585,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.50828,0.0616464,-1.79746e-05,-3.27952e-08,2.05646e-11,28195.8,27.7757], Tmin=(100,'K'), Tmax=(935.225,'K')), NASAPolynomial(coeffs=[20.5413,0.014238,-3.32393e-06,5.2059e-10,-3.88035e-14,22775,-76.4793], Tmin=(935.225,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(233.272,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = 'CC=C[CH]C=[C]O(23618)',
    structure = SMILES('CC=C[CH]C=[C]O'),
    E0 = (185.401,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06888,0.0489349,8.84299e-06,-5.42632e-08,2.6872e-11,22418.4,27.8081], Tmin=(100,'K'), Tmax=(932.402,'K')), NASAPolynomial(coeffs=[17.5481,0.0177376,-4.51187e-06,7.18936e-10,-5.18948e-14,17628.4,-59.7484], Tmin=(932.402,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(185.401,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CCJC=C) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C]=CCC=CO(23619)',
    structure = SMILES('[CH2][C]=CCC=CO'),
    E0 = (233.272,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,3615,1277.5,1000,3000,3100,440,815,1455,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.894772,'amu*angstrom^2'), symmetry=1, barrier=(20.5726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.895795,'amu*angstrom^2'), symmetry=1, barrier=(20.5961,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.894823,'amu*angstrom^2'), symmetry=1, barrier=(20.5737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.895314,'amu*angstrom^2'), symmetry=1, barrier=(20.585,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.50828,0.0616464,-1.79746e-05,-3.27952e-08,2.05646e-11,28195.8,27.7757], Tmin=(100,'K'), Tmax=(935.225,'K')), NASAPolynomial(coeffs=[20.5413,0.014238,-3.32393e-06,5.2059e-10,-3.88035e-14,22775,-76.4793], Tmin=(935.225,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(233.272,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = 'CC=CC[C]=[C]O(23620)',
    structure = SMILES('CC=CC[C]=[C]O'),
    E0 = (321.517,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,1670,1700,300,440,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,191.912,193.426],'cm^-1')),
        HinderedRotor(inertia=(0.465503,'amu*angstrom^2'), symmetry=1, barrier=(12.5844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.461016,'amu*angstrom^2'), symmetry=1, barrier=(12.582,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.462048,'amu*angstrom^2'), symmetry=1, barrier=(12.5768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.472012,'amu*angstrom^2'), symmetry=1, barrier=(12.5972,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.526236,0.0680975,-6.07589e-05,2.83443e-08,-5.25795e-12,38801.7,30.5631], Tmin=(100,'K'), Tmax=(1306.12,'K')), NASAPolynomial(coeffs=[15.5725,0.0220181,-7.83935e-06,1.33311e-09,-8.78251e-14,34871.2,-46.0445], Tmin=(1306.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(321.517,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = 'CC=CC[CH][C]=O(20130)',
    structure = SMILES('CC=CC[CH][C]=O'),
    E0 = (168.363,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,1855,455,950,2995,3025,975,1000,1300,1375,400,500,1630,1680,256.368,256.411],'cm^-1')),
        HinderedRotor(inertia=(0.00256398,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169682,'amu*angstrom^2'), symmetry=1, barrier=(7.91637,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00256372,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.310014,'amu*angstrom^2'), symmetry=1, barrier=(14.467,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24409,0.058217,-4.29306e-05,1.72311e-08,-2.90963e-12,20350.4,26.6305], Tmin=(100,'K'), Tmax=(1360.13,'K')), NASAPolynomial(coeffs=[10.6868,0.0304463,-1.23036e-05,2.21896e-09,-1.50257e-13,17781.8,-21.8294], Tmin=(1360.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(168.363,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(CCJCHO)"""),
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
    label = '[CH2]C=[C]O(6303)',
    structure = SMILES('[CH2]C=[C]O'),
    E0 = (224.294,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1277.5,1000,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.206965,'amu*angstrom^2'), symmetry=1, barrier=(4.75853,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.69007,'amu*angstrom^2'), symmetry=1, barrier=(15.8661,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.71481,0.0306825,-3.73315e-05,3.08792e-08,-1.09328e-11,27020.3,15.471], Tmin=(100,'K'), Tmax=(771.374,'K')), NASAPolynomial(coeffs=[4.68419,0.0180544,-8.07735e-06,1.53601e-09,-1.0692e-13,26788.3,6.94699], Tmin=(771.374,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(224.294,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=CJO) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=[C]O(172)',
    structure = SMILES('[CH]=[C]O'),
    E0 = (348.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3120,650,792.5,1650,1058.91,1059.8],'cm^-1')),
        HinderedRotor(inertia=(0.315045,'amu*angstrom^2'), symmetry=1, barrier=(7.24351,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.22339,0.0198737,-3.18466e-05,3.10436e-08,-1.17427e-11,41917.3,13.7606], Tmin=(100,'K'), Tmax=(829.31,'K')), NASAPolynomial(coeffs=[3.78101,0.0111997,-5.33323e-06,1.02849e-09,-7.13893e-14,42030.6,12.4155], Tmin=(829.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(348.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(Cds_P) + radical(C=CJO)"""),
)

species(
    label = 'O[C]=CCC1[CH]C1(23621)',
    structure = SMILES('O[C]=CCC1[CH]C1'),
    E0 = (346.337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22801,0.045524,1.32022e-05,-5.46149e-08,2.5661e-11,41768.6,29.1248], Tmin=(100,'K'), Tmax=(951.063,'K')), NASAPolynomial(coeffs=[16.5183,0.0194338,-5.92568e-06,1.04549e-09,-7.67077e-14,37131.7,-52.9618], Tmin=(951.063,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(346.337,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclopropane) + radical(cyclopropane) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C1[CH]CC=C1O(23550)',
    structure = SMILES('[CH2]C1[CH]CC=C1O'),
    E0 = (167.77,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47162,0.0365978,4.29619e-05,-8.84803e-08,3.86928e-11,20286.6,25.2327], Tmin=(100,'K'), Tmax=(933.596,'K')), NASAPolynomial(coeffs=[17.2368,0.0174643,-4.08035e-06,6.5618e-10,-5.00365e-14,15233.1,-61.0414], Tmin=(933.596,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(167.77,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(cyclopentene-4) + radical(Isobutyl)"""),
)

species(
    label = 'C=C=CCC=CO(23622)',
    structure = SMILES('C=C=CCC=CO'),
    E0 = (20.5358,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.580229,0.0591917,-1.08129e-05,-4.02349e-08,2.31783e-11,2607.74,25.8168], Tmin=(100,'K'), Tmax=(936.018,'K')), NASAPolynomial(coeffs=[20.6413,0.0139334,-3.14137e-06,4.94372e-10,-3.77492e-14,-2920.67,-79.1103], Tmin=(936.018,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.5358,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'CC=CCC=C=O(20133)',
    structure = SMILES('CC=CCC=C=O'),
    E0 = (-33.7596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44166,0.0596874,-4.82584e-05,2.27301e-08,-4.73781e-12,-3971.11,23.0818], Tmin=(100,'K'), Tmax=(1085.87,'K')), NASAPolynomial(coeffs=[7.7278,0.0365311,-1.62705e-05,3.09114e-09,-2.16292e-13,-5336.29,-7.763], Tmin=(1085.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-33.7596,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-Cd(CCO)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH)"""),
)

species(
    label = 'OC1=CCC=CC1(23342)',
    structure = SMILES('OC1=CCC=CC1'),
    E0 = (-120.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50541,0.0318462,6.35903e-05,-1.1223e-07,4.71331e-11,-14410.9,21.664], Tmin=(100,'K'), Tmax=(944.15,'K')), NASAPolynomial(coeffs=[18.9591,0.0161055,-3.87195e-06,6.98663e-10,-5.79402e-14,-20300.9,-75.2751], Tmin=(944.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-120.741,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cyclohexadiene)"""),
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
    label = '[CH]=CCC=[C]O(23623)',
    structure = SMILES('[CH]=CCC=[C]O'),
    E0 = (366.797,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,3615,1277.5,1000,3120,650,792.5,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(0.724801,'amu*angstrom^2'), symmetry=1, barrier=(16.6646,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.724721,'amu*angstrom^2'), symmetry=1, barrier=(16.6628,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.725348,'amu*angstrom^2'), symmetry=1, barrier=(16.6772,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33661,0.0490966,-2.41838e-05,-1.16682e-08,1.04359e-11,44220.1,25.7742], Tmin=(100,'K'), Tmax=(936.975,'K')), NASAPolynomial(coeffs=[15.5183,0.012305,-3.30717e-06,5.31473e-10,-3.75146e-14,40520,-47.2846], Tmin=(936.975,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(366.797,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = 'C=CC=CC=[C]O(23624)',
    structure = SMILES('C=CC=CC=[C]O'),
    E0 = (176.645,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,3615,1277.5,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.02253,'amu*angstrom^2'), symmetry=1, barrier=(23.5101,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02476,'amu*angstrom^2'), symmetry=1, barrier=(23.5613,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02272,'amu*angstrom^2'), symmetry=1, barrier=(23.5142,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.558767,0.0589195,-1.17097e-05,-4.45582e-08,2.67326e-11,21385,25.9711], Tmin=(100,'K'), Tmax=(913.538,'K')), NASAPolynomial(coeffs=[22.6267,0.00681664,7.35319e-07,-2.89898e-10,1.81949e-14,15495.1,-88.667], Tmin=(913.538,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(176.645,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO)"""),
)

species(
    label = 'C=CC[CH]C=[C]O(23568)',
    structure = SMILES('C=CC[CH]C=[C]O'),
    E0 = (235.764,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,355.154,355.155,355.158],'cm^-1')),
        HinderedRotor(inertia=(0.00133645,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.212108,'amu*angstrom^2'), symmetry=1, barrier=(18.9862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.212109,'amu*angstrom^2'), symmetry=1, barrier=(18.9862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.21211,'amu*angstrom^2'), symmetry=1, barrier=(18.9861,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.878558,0.0562531,-1.65872e-05,-2.37346e-08,1.46072e-11,28479.4,28.1256], Tmin=(100,'K'), Tmax=(963.405,'K')), NASAPolynomial(coeffs=[16.5385,0.0206283,-6.88621e-06,1.22209e-09,-8.70956e-14,24098,-53.9204], Tmin=(963.405,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(235.764,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(Allyl_S)"""),
)

species(
    label = 'C=[C]CCC=[C]O(23625)',
    structure = SMILES('C=[C]CCC=[C]O'),
    E0 = (332.493,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,1670,1700,300,440,3615,1277.5,1000,3010,987.5,1337.5,450,1655,332.364,332.397,332.536],'cm^-1')),
        HinderedRotor(inertia=(0.168003,'amu*angstrom^2'), symmetry=1, barrier=(13.1662,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168185,'amu*angstrom^2'), symmetry=1, barrier=(13.1635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167996,'amu*angstrom^2'), symmetry=1, barrier=(13.1655,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168112,'amu*angstrom^2'), symmetry=1, barrier=(13.1645,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.614817,0.066462,-5.16687e-05,1.57396e-08,-5.04226e-14,40118.6,30.6524], Tmin=(100,'K'), Tmax=(1008.6,'K')), NASAPolynomial(coeffs=[15.4593,0.0223555,-8.03201e-06,1.41116e-09,-9.66082e-14,36373.1,-44.814], Tmin=(1008.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(332.493,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = 'C=CCC[C]=[C]O(23626)',
    structure = SMILES('C=CCC[C]=[C]O'),
    E0 = (332.493,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,1670,1700,300,440,3615,1277.5,1000,3010,987.5,1337.5,450,1655,332.364,332.397,332.536],'cm^-1')),
        HinderedRotor(inertia=(0.168003,'amu*angstrom^2'), symmetry=1, barrier=(13.1662,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168185,'amu*angstrom^2'), symmetry=1, barrier=(13.1635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167996,'amu*angstrom^2'), symmetry=1, barrier=(13.1655,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168112,'amu*angstrom^2'), symmetry=1, barrier=(13.1645,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.614817,0.066462,-5.16687e-05,1.57396e-08,-5.04226e-14,40118.6,30.6524], Tmin=(100,'K'), Tmax=(1008.6,'K')), NASAPolynomial(coeffs=[15.4593,0.0223555,-8.03201e-06,1.41116e-09,-9.66082e-14,36373.1,-44.814], Tmin=(1008.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(332.493,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CCCC=[C]O(23627)',
    structure = SMILES('[CH]=CCCC=[C]O'),
    E0 = (341.748,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1685,370,3615,1277.5,1000,3120,650,792.5,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,260.147,260.159],'cm^-1')),
        HinderedRotor(inertia=(0.298971,'amu*angstrom^2'), symmetry=1, barrier=(14.356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.298897,'amu*angstrom^2'), symmetry=1, barrier=(14.3558,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.298892,'amu*angstrom^2'), symmetry=1, barrier=(14.3558,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.29896,'amu*angstrom^2'), symmetry=1, barrier=(14.3561,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.597488,0.0651102,-4.25901e-05,2.65179e-09,5.4949e-12,41234,30.644], Tmin=(100,'K'), Tmax=(970.067,'K')), NASAPolynomial(coeffs=[16.7451,0.0202645,-6.85872e-06,1.19588e-09,-8.30959e-14,37078.3,-52.04], Tmin=(970.067,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(341.748,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = 'C=CCC[CH][C]=O(20138)',
    structure = SMILES('C=CCC[CH][C]=O'),
    E0 = (180.608,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,1855,455,950,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,295.815,299.129,300.345],'cm^-1')),
        HinderedRotor(inertia=(0.17616,'amu*angstrom^2'), symmetry=1, barrier=(10.955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.330855,'amu*angstrom^2'), symmetry=1, barrier=(21.3586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00144199,'amu*angstrom^2'), symmetry=1, barrier=(16.1456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0577846,'amu*angstrom^2'), symmetry=1, barrier=(90.7382,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01996,0.0600114,-4.50256e-05,1.80224e-08,-2.96958e-12,21833.9,27.8781], Tmin=(100,'K'), Tmax=(1418.62,'K')), NASAPolynomial(coeffs=[12.5362,0.0275399,-1.06914e-05,1.88736e-09,-1.26149e-13,18566.5,-31.708], Tmin=(1418.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(180.608,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]=C[CH]CC=CO(23628)',
    structure = SMILES('[CH]C=CCC=CO'),
    E0 = (214.616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.540681,0.0594741,-1.37745e-06,-4.92639e-08,2.57788e-11,25952.1,28.0139], Tmin=(100,'K'), Tmax=(940.315,'K')), NASAPolynomial(coeffs=[19.3656,0.0209084,-6.08037e-06,1.0215e-09,-7.32736e-14,20576.5,-71.4061], Tmin=(940.315,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.616,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'OC1=CC[CH][CH]C1(23581)',
    structure = SMILES('OC1=CC[CH][CH]C1'),
    E0 = (153.089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9188,0.0330257,3.08725e-05,-5.77903e-08,2.29171e-11,18498.6,23.8658], Tmin=(100,'K'), Tmax=(990.371,'K')), NASAPolynomial(coeffs=[10.3334,0.0299413,-1.12581e-05,2.07465e-09,-1.47452e-13,15316.4,-24.2991], Tmin=(990.371,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(153.089,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclohexene) + radical(RCCJCC) + radical(RCCJCC)"""),
)

species(
    label = 'C=CC=CC=CO(23629)',
    structure = SMILES('C=CC=CC=CO'),
    E0 = (-63.0995,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.342679,0.0572292,1.49765e-05,-8.31724e-08,4.29321e-11,-7435.74,23.6184], Tmin=(100,'K'), Tmax=(911.518,'K')), NASAPolynomial(coeffs=[26.7946,0.00264365,3.6095e-06,-8.48075e-10,5.4456e-14,-14812.6,-115.559], Tmin=(911.518,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-63.0995,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CCCC=C=O(20143)',
    structure = SMILES('C=CCCC=C=O'),
    E0 = (-21.4616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06199,0.0573791,-3.83733e-05,1.29098e-08,-1.76247e-12,-2469.58,26.3395], Tmin=(100,'K'), Tmax=(1688.18,'K')), NASAPolynomial(coeffs=[14.3948,0.0257886,-1.03043e-05,1.82539e-09,-1.21016e-13,-6971.26,-44.9654], Tmin=(1688.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-21.4616,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2]C(C=C)C=[C]O(23341)',
    structure = SMILES('[CH2]C(C=C)C=[C]O'),
    E0 = (291.214,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,308.581,311.299],'cm^-1')),
        HinderedRotor(inertia=(0.19307,'amu*angstrom^2'), symmetry=1, barrier=(13.1808,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.192811,'amu*angstrom^2'), symmetry=1, barrier=(13.1621,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.233696,'amu*angstrom^2'), symmetry=1, barrier=(15.8598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.190477,'amu*angstrom^2'), symmetry=1, barrier=(13.1799,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3876.19,'J/mol'), sigma=(6.45323,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=605.45 K, Pc=32.73 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.516779,0.0713614,-6.49459e-05,2.71006e-08,-2.87808e-12,35155.2,29.6822], Tmin=(100,'K'), Tmax=(928.049,'K')), NASAPolynomial(coeffs=[15.0033,0.0222931,-7.24732e-06,1.17648e-09,-7.64266e-14,31890.5,-42.2273], Tmin=(928.049,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(291.214,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(Isobutyl)"""),
)

species(
    label = 'C=CC1CC=C1O(23552)',
    structure = SMILES('C=CC1CC=C1O'),
    E0 = (1.61078,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05915,0.0473195,1.59035e-05,-6.18584e-08,2.91718e-11,315.445,21.5518], Tmin=(100,'K'), Tmax=(947.823,'K')), NASAPolynomial(coeffs=[18.4513,0.0171651,-4.812e-06,8.48424e-10,-6.4447e-14,-4923.95,-71.6703], Tmin=(947.823,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1.61078,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutene)"""),
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
    E0 = (235.175,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (288.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (360.238,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (483.871,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (456.352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (408.77,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (397.954,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (442.906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (467.33,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (483.438,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (384.544,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (627.577,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (352.745,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (307.095,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (382.649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (316.159,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (601.325,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (623.024,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (461.111,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (273.249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (243.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (260.148,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (243.376,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (748.36,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (396.526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (449.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (351.833,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (490.173,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (457.037,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (486.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (436.473,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (417.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (246.423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (298.575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (274.4,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (451.149,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (242.706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C[CH]CC=[C]O(23340)'],
    products = ['HCCOH(50)', 'butadiene13(2459)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C[CH]CC=[C]O(23340)'],
    products = ['[CH2][CH]C1CC=C1O(23610)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(53.5552,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS_D;doublebond_intra;radadd_intra_cdsingle] for rate rule [R5_DS_D;doublebond_intra;radadd_intra_cdsingleNd]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', '[CH2]C=CCC=C=O(20121)'],
    products = ['C=C[CH]CC=[C]O(23340)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.185e+08,'cm^3/(mol*s)'), n=1.63, Ea=(30.7064,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1700,'K'), comment="""Estimated using template [Od_R;HJ] for rate rule [Od_Cdd;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C=C=CCC=[C]O(23611)'],
    products = ['C=C[CH]CC=[C]O(23340)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.42e+08,'cm^3/(mol*s)'), n=1.64, Ea=(11.7989,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2713 used for Ca_Cds-HH;HJ
Exact match found for rate rule [Ca_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', '[CH2]C=CCC#CO(23612)'],
    products = ['C=C[CH]CC=[C]O(23340)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['HCCOH(50)', '[CH2][CH]C=C(2458)'],
    products = ['C=C[CH]CC=[C]O(23340)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.113192,'m^3/(mol*s)'), n=2.418, Ea=(54.0165,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CsJ-CdHH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=C[CH]CC=[C]O(23340)'],
    products = ['[CH2]C=CC[CH]C=O(20123)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=C[CH]CC=[C]O(23340)'],
    products = ['C[C]=CCC=[C]O(23613)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C=CC[C]=CO(23614)'],
    products = ['C=C[CH]CC=[C]O(23340)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.231e+11,'s^-1'), n=0.765, Ea=(234.057,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 82 used for R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CC=[C]CC=[C]O(23615)'],
    products = ['C=C[CH]CC=[C]O(23340)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=C[CH]CC=[C]O(23340)'],
    products = ['[CH2]C=C[CH]C=CO(23616)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.182e+10,'s^-1'), n=0.86, Ea=(149.369,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleNd;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleNd;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C=[C]CC=CO(23617)'],
    products = ['C=C[CH]CC=[C]O(23340)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(9.01194e+11,'s^-1'), n=1.09397, Ea=(394.305,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSD;Cd_rad_out;Cd_H_out_single] for rate rule [R4H_SSD;Cd_rad_out_Cd;Cd_H_out_singleNd]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C[CH]CC=[C]O(23340)'],
    products = ['CC=C[CH]C=[C]O(23618)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(256000,'s^-1'), n=2, Ea=(117.57,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SDS;C_rad_out_2H;Cs_H_out_H/OneDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][C]=CCC=CO(23619)'],
    products = ['C=C[CH]CC=[C]O(23340)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(272058,'s^-1'), n=1.7475, Ea=(73.8225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;Cd_H_out_single] for rate rule [R5H_DSSD;Cd_rad_out;Cd_H_out_singleNd]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CC=CC[C]=[C]O(23620)'],
    products = ['C=C[CH]CC=[C]O(23340)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(770185,'s^-1'), n=1.81245, Ea=(61.132,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5H_RSMS;Cd_rad_out;Cs_H_out_2H] + [R5H_SSMS;Y_rad_out;Cs_H_out_2H] for rate rule [R5H_SSMS;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=C[CH]CC=[C]O(23340)'],
    products = ['CC=CC[CH][C]=O(20130)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(60975.7,'s^-1'), n=1.58648, Ea=(80.9836,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7Hall;C_rad_out_2H;XH_out] for rate rule [R7HJ_5;C_rad_out_2H;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C[CH2](2461)', '[CH2]C=[C]O(6303)'],
    products = ['C=C[CH]CC=[C]O(23340)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.169e+08,'m^3/(mol*s)'), n=-0.455312, Ea=(0.377199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [C_pri_rad;Cd_rad] + [C_rad/H2/Cd;Y_rad] for rate rule [C_rad/H2/Cd;Cd_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=[C]O(172)', '[CH2][CH]C=C(2458)'],
    products = ['C=C[CH]CC=[C]O(23340)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.13324e+07,'m^3/(mol*s)'), n=0.074875, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_rad/H2/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -0.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=C[CH]CC=[C]O(23340)'],
    products = ['O[C]=CCC1[CH]C1(23621)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra_cs2H] for rate rule [R3_D;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C[CH]CC=[C]O(23340)'],
    products = ['[CH2]C1[CH]CC=C1O(23550)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.53e+07,'s^-1'), n=1.05, Ea=(38.0744,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_DS_D;doublebond_intra_pri;radadd_intra_cdsingleNd]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=C[CH]CC=[C]O(23340)'],
    products = ['C=C=CCC=CO(23622)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.21e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=C[CH]CC=[C]O(23340)'],
    products = ['CC=CCC=C=O(20133)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C[CH]CC=[C]O(23340)'],
    products = ['OC1=CCC=CC1(23342)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.53377e+11,'s^-1'), n=0.0685, Ea=(8.20064,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad_out;Cpri_rad_out_2H] for rate rule [R6;CdsingleND_rad_out;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CH2(19)', '[CH]=CCC=[C]O(23623)'],
    products = ['C=C[CH]CC=[C]O(23340)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)', 'C=CC=CC=[C]O(23624)'],
    products = ['C=C[CH]CC=[C]O(23340)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(105.863,'m^3/(mol*s)'), n=1.66278, Ea=(8.08939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-OneDeH_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=[C]O(172)', 'butadiene13(2459)'],
    products = ['C=C[CH]CC=[C]O(23340)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(0.0534234,'m^3/(mol*s)'), n=2.459, Ea=(4.91982,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CdH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=CC[CH]C=[C]O(23568)'],
    products = ['C=C[CH]CC=[C]O(23340)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.169e+11,'s^-1'), n=0.707, Ea=(116.068,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/Cd;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=[C]CCC=[C]O(23625)'],
    products = ['C=C[CH]CC=[C]O(23340)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.66329e+10,'s^-1'), n=0.993, Ea=(157.679,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=CCC[C]=[C]O(23626)'],
    products = ['C=C[CH]CC=[C]O(23340)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.29711e+07,'s^-1'), n=1.52333, Ea=(124.544,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R3H_SS_Cs;Cd_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]=CCCC=[C]O(23627)'],
    products = ['C=C[CH]CC=[C]O(23340)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 194 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=C[CH]CC=[C]O(23340)'],
    products = ['C=CCC[CH][C]=O(20138)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.79021e+08,'s^-1'), n=1.23491, Ea=(201.298,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;C_rad_out_H/Cd;XH_out] for rate rule [R5HJ_3;C_rad_out_H/Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=C[CH]CC=[C]O(23340)'],
    products = ['[CH]=C[CH]CC=CO(23628)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.74e+08,'s^-1'), n=1.713, Ea=(181.895,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleNd;Cd_H_out_singleH] for rate rule [R6HJ_3;Cd_rad_out_singleNd;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=C[CH]CC=[C]O(23340)'],
    products = ['OC1=CC[CH][CH]C1(23581)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.76666e+15,'s^-1'), n=-0.518223, Ea=(11.2479,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6_linear;doublebond_intra_pri;radadd_intra_cdsingleNd] + [R6_linear;doublebond_intra_pri_2H;radadd_intra_cdsingle] for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_cdsingleNd]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=C[CH]CC=[C]O(23340)'],
    products = ['C=CC=CC=CO(23629)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3radExo;Y_rad_NDe;XH_Rrad] for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=C[CH]CC=[C]O(23340)'],
    products = ['C=CCCC=C=O(20143)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(C=C)C=[C]O(23341)'],
    products = ['C=C[CH]CC=[C]O(23340)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=C[CH]CC=[C]O(23340)'],
    products = ['C=CC1CC=C1O(23552)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/OneDe;CdsinglepriND_rad_out]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

network(
    label = '4030',
    isomers = [
        'C=C[CH]CC=[C]O(23340)',
    ],
    reactants = [
        ('HCCOH(50)', 'butadiene13(2459)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4030',
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

