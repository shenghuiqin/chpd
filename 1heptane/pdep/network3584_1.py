species(
    label = '[CH2]C(C=C)C(=C)[O](19794)',
    structure = SMILES('[CH2]C(C=C)C(=C)[O]'),
    E0 = (183.509,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,327.845,329.205,2023.39],'cm^-1')),
        HinderedRotor(inertia=(0.170779,'amu*angstrom^2'), symmetry=1, barrier=(13.0305,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169947,'amu*angstrom^2'), symmetry=1, barrier=(13.0443,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169751,'amu*angstrom^2'), symmetry=1, barrier=(13.012,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.648499,0.0693304,-6.62578e-05,3.47966e-08,-7.35885e-12,22195.4,27.7785], Tmin=(100,'K'), Tmax=(1145.92,'K')), NASAPolynomial(coeffs=[13.1664,0.0256346,-9.06e-06,1.52019e-09,-9.90664e-14,19326.5,-34.3182], Tmin=(1145.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(183.509,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(C=C(C)OJ)"""),
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
    label = '[CH2]C1CC1C(=C)[O](20207)',
    structure = SMILES('[CH2]C1CC1C(=C)[O]'),
    E0 = (209.829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32812,0.0405335,3.38706e-05,-8.18122e-08,3.73566e-11,25349.7,25.8333], Tmin=(100,'K'), Tmax=(921.336,'K')), NASAPolynomial(coeffs=[17.8617,0.0159362,-2.90154e-06,3.80369e-10,-2.85475e-14,20300.5,-63.445], Tmin=(921.336,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(209.829,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(C=C(C)OJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1([O])CC1C=C(20208)',
    structure = SMILES('[CH2]C1([O])CC1C=C'),
    E0 = (341.02,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.828158,0.0559171,-1.09926e-05,-3.08527e-08,1.71848e-11,41141.9,24.5921], Tmin=(100,'K'), Tmax=(969.289,'K')), NASAPolynomial(coeffs=[17.3185,0.0206127,-7.03523e-06,1.28056e-09,-9.29318e-14,36406.8,-62.3853], Tmin=(969.289,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(341.02,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(CJC(C)2O) + radical(CC(C)2OJ)"""),
)

species(
    label = '[CH2]C1OC(=C)C1[CH2](20085)',
    structure = SMILES('[CH2]C1OC(=C)C1[CH2]'),
    E0 = (292.18,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.752835,0.0536933,1.10751e-05,-7.17822e-08,3.83757e-11,35274.9,21.3155], Tmin=(100,'K'), Tmax=(882.568,'K')), NASAPolynomial(coeffs=[21.655,0.00898682,2.03241e-06,-7.25972e-10,5.51671e-14,29637,-87.9526], Tmin=(882.568,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(292.18,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(CJC(C)OC) + radical(Isobutyl)"""),
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
    label = 'C=CC(=C)C(=C)[O](20209)',
    structure = SMILES('C=CC(=C)C(=C)[O]'),
    E0 = (88.7055,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,2950,2980,3010,3040,3070,3100,1330,1380,1430,900,975,1050,1000,1025,1050,1600,1650,1700,341.074,341.076,341.077],'cm^-1')),
        HinderedRotor(inertia=(0.247527,'amu*angstrom^2'), symmetry=1, barrier=(20.4326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.247519,'amu*angstrom^2'), symmetry=1, barrier=(20.4325,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11946,0.0538418,-2.29787e-05,-1.13228e-08,9.15461e-12,10781,22.115], Tmin=(100,'K'), Tmax=(969.684,'K')), NASAPolynomial(coeffs=[14.1031,0.0221045,-7.63969e-06,1.33844e-09,-9.27972e-14,7237.12,-45.4141], Tmin=(969.684,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(88.7055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C2H3(30)',
    structure = SMILES('[CH]=C'),
    E0 = (286.361,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,677.08,1086.68,3788.01],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (27.0452,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36378,0.000265766,2.79621e-05,-3.72987e-08,1.5159e-11,34475,7.9151], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.15027,0.00754021,-2.62998e-06,4.15974e-10,-2.45408e-14,33856.6,1.72812], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(286.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""C2H3""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C=CC(=C)[O](2398)',
    structure = SMILES('C=CC(=C)[O]'),
    E0 = (17.8331,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.914251,'amu*angstrom^2'), symmetry=1, barrier=(21.0204,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98668,0.0359126,-9.09545e-06,-2.04235e-08,1.25939e-11,2225.07,15.9658], Tmin=(100,'K'), Tmax=(920.927,'K')), NASAPolynomial(coeffs=[13.0167,0.00983409,-2.17548e-06,3.0681e-10,-2.11581e-14,-732.214,-41.3653], Tmin=(920.927,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(17.8331,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
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
    label = 'C=C[C](C)C(=C)[O](20210)',
    structure = SMILES('[CH2]C=C(C)C(=C)[O]'),
    E0 = (62.2727,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.753129,0.0628798,-3.80576e-05,6.82091e-10,5.72992e-12,7614.26,24.357], Tmin=(100,'K'), Tmax=(956.437,'K')), NASAPolynomial(coeffs=[14.7682,0.0240139,-8.07407e-06,1.37016e-09,-9.26407e-14,4030.13,-47.3551], Tmin=(956.437,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(62.2727,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=[C]C(C)C(=C)[O](20211)',
    structure = SMILES('C=[C]C(C)C(=C)[O]'),
    E0 = (216.268,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,187.83,187.853,187.854],'cm^-1')),
        HinderedRotor(inertia=(0.00477816,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.406471,'amu*angstrom^2'), symmetry=1, barrier=(10.1774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.406488,'amu*angstrom^2'), symmetry=1, barrier=(10.1773,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.869175,0.0689848,-6.69142e-05,3.62267e-08,-8.08252e-12,26123.7,25.6833], Tmin=(100,'K'), Tmax=(1071.34,'K')), NASAPolynomial(coeffs=[11.3137,0.0299888,-1.23156e-05,2.25161e-09,-1.54398e-13,23885.8,-25.4255], Tmin=(1071.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.268,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2][C](C=C)C(=C)O(20212)',
    structure = SMILES('[CH2]C=C([CH2])C(=C)O'),
    E0 = (75.967,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.415839,0.0639851,-2.0497e-05,-3.20472e-08,2.09296e-11,9279.5,24.3486], Tmin=(100,'K'), Tmax=(924.187,'K')), NASAPolynomial(coeffs=[20.6449,0.0148656,-3.1546e-06,4.41571e-10,-3.13973e-14,3899.04,-80.5301], Tmin=(924.187,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(75.967,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]=C(O)C([CH2])C=C(20213)',
    structure = SMILES('[CH]=C(O)C([CH2])C=C'),
    E0 = (292.8,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,275.637],'cm^-1')),
        HinderedRotor(inertia=(0.275923,'amu*angstrom^2'), symmetry=1, barrier=(14.6373,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.268004,'amu*angstrom^2'), symmetry=1, barrier=(14.6283,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.28803,'amu*angstrom^2'), symmetry=1, barrier=(14.5784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.270258,'amu*angstrom^2'), symmetry=1, barrier=(14.6165,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0218399,0.0790318,-8.31181e-05,4.51014e-08,-9.51583e-12,35366.2,28.9291], Tmin=(100,'K'), Tmax=(1200.89,'K')), NASAPolynomial(coeffs=[17.6546,0.018588,-5.48155e-06,8.15191e-10,-4.93282e-14,31254.6,-58.853], Tmin=(1200.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(292.8,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([C]=C)C(=C)O(20214)',
    structure = SMILES('[CH2]C([C]=C)C(=C)O'),
    E0 = (283.546,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,354.76,354.761],'cm^-1')),
        HinderedRotor(inertia=(0.148912,'amu*angstrom^2'), symmetry=1, barrier=(13.2992,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148911,'amu*angstrom^2'), symmetry=1, barrier=(13.2992,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148912,'amu*angstrom^2'), symmetry=1, barrier=(13.2992,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148912,'amu*angstrom^2'), symmetry=1, barrier=(13.2992,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.234952,0.0780816,-8.41942e-05,4.79074e-08,-1.07255e-11,34242.3,28.2349], Tmin=(100,'K'), Tmax=(1096.66,'K')), NASAPolynomial(coeffs=[15.6371,0.0219023,-7.35155e-06,1.1935e-09,-7.62657e-14,30864.1,-47.4926], Tmin=(1096.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(283.546,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC(C)C(=C)[O](20215)',
    structure = SMILES('[CH]=CC(C)C(=C)[O]'),
    E0 = (225.522,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,250.362,250.616],'cm^-1')),
        HinderedRotor(inertia=(0.23706,'amu*angstrom^2'), symmetry=1, barrier=(10.6093,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.237911,'amu*angstrom^2'), symmetry=1, barrier=(10.6061,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.236865,'amu*angstrom^2'), symmetry=1, barrier=(10.6047,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.678248,0.0697169,-6.52479e-05,3.28424e-08,-6.68524e-12,27246.5,26.2946], Tmin=(100,'K'), Tmax=(1181.37,'K')), NASAPolynomial(coeffs=[13.4569,0.0264496,-1.0311e-05,1.84062e-09,-1.24692e-13,24227.2,-37.485], Tmin=(1181.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(225.522,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C([O])C(C)C=C(20216)',
    structure = SMILES('[CH]=C([O])C(C)C=C'),
    E0 = (225.522,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,250.362,250.616],'cm^-1')),
        HinderedRotor(inertia=(0.23706,'amu*angstrom^2'), symmetry=1, barrier=(10.6093,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.237911,'amu*angstrom^2'), symmetry=1, barrier=(10.6061,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.236865,'amu*angstrom^2'), symmetry=1, barrier=(10.6047,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.678248,0.0697169,-6.52479e-05,3.28424e-08,-6.68524e-12,27246.5,26.2946], Tmin=(100,'K'), Tmax=(1181.37,'K')), NASAPolynomial(coeffs=[13.4569,0.0264496,-1.0311e-05,1.84062e-09,-1.24692e-13,24227.2,-37.485], Tmin=(1181.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(225.522,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=CC([CH2])C(=C)O(20217)',
    structure = SMILES('[CH]=CC([CH2])C(=C)O'),
    E0 = (292.8,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0218399,0.0790318,-8.31181e-05,4.51014e-08,-9.51583e-12,35366.2,28.9291], Tmin=(100,'K'), Tmax=(1200.89,'K')), NASAPolynomial(coeffs=[17.6546,0.018588,-5.48155e-06,8.15191e-10,-4.93282e-14,31254.6,-58.853], Tmin=(1200.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(292.8,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C=C)[C]1CO1(20218)',
    structure = SMILES('[CH2]C(C=C)[C]1CO1'),
    E0 = (330.226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.498366,0.065143,-5.74863e-05,2.91761e-08,-5.76647e-12,39853.3,29.1762], Tmin=(100,'K'), Tmax=(1443.3,'K')), NASAPolynomial(coeffs=[11.7984,0.0243842,-5.3138e-06,5.44931e-10,-2.20623e-14,37574.8,-26.0797], Tmin=(1443.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(330.226,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(C2CsJO) + radical(Isobutyl)"""),
)

species(
    label = 'C=C([O])C1[CH]CC1(20219)',
    structure = SMILES('C=C([O])C1[CH]CC1'),
    E0 = (197.208,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70816,0.0334554,4.23199e-05,-7.84343e-08,3.22604e-11,23816.4,26.0209], Tmin=(100,'K'), Tmax=(965.289,'K')), NASAPolynomial(coeffs=[14.1079,0.024036,-8.25086e-06,1.527e-09,-1.12425e-13,19467.5,-43.489], Tmin=(965.289,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(197.208,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(C=C(C)OJ) + radical(cyclobutane)"""),
)

species(
    label = '[CH2][CH]C1CCC1=O(20120)',
    structure = SMILES('[CH2][CH]C1CCC1=O'),
    E0 = (235.898,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.82988,0.034022,3.1015e-05,-5.91249e-08,2.33943e-11,28462.2,27.1889], Tmin=(100,'K'), Tmax=(1001.14,'K')), NASAPolynomial(coeffs=[11.315,0.0292287,-1.14027e-05,2.15009e-09,-1.54796e-13,24904,-26.8673], Tmin=(1001.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(235.898,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(RCCJ) + radical(CCJCC=O)"""),
)

species(
    label = '[CH2]C1[CH]COC1=C(20220)',
    structure = SMILES('[CH2]C1[CH]COC1=C'),
    E0 = (223.479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80377,0.0199424,0.00010405,-1.61264e-07,6.70887e-11,26983.7,24.5843], Tmin=(100,'K'), Tmax=(921.301,'K')), NASAPolynomial(coeffs=[20.8885,0.010252,6.97364e-07,-2.72883e-10,1.08539e-14,20361.8,-82.7771], Tmin=(921.301,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(223.479,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(Isobutyl) + radical(CCJCO)"""),
)

species(
    label = 'C=CC(=C)C(=C)O(20221)',
    structure = SMILES('C=CC(=C)C(=C)O'),
    E0 = (-49.0993,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.770347,0.0576068,-1.46885e-05,-2.95e-08,1.76938e-11,-5776.76,21.8673], Tmin=(100,'K'), Tmax=(946.952,'K')), NASAPolynomial(coeffs=[17.748,0.0189038,-5.67397e-06,9.67926e-10,-6.91199e-14,-10472.3,-66.9299], Tmin=(946.952,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-49.0993,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CC[CH]C(=C)[O](20222)',
    structure = SMILES('[CH2]C([O])=CCC=C'),
    E0 = (134.887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.874254,0.0594852,-3.2007e-05,-4.44888e-09,7.21713e-12,16344,28.096], Tmin=(100,'K'), Tmax=(967.476,'K')), NASAPolynomial(coeffs=[14.8374,0.022807,-7.78003e-06,1.3482e-09,-9.27239e-14,12657,-43.8984], Tmin=(967.476,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(134.887,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(O)CJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C[CH]CC(=C)[O](19793)',
    structure = SMILES('[CH2]C=CCC(=C)[O]'),
    E0 = (127.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.924674,0.0562641,-1.95367e-05,-1.80395e-08,1.18943e-11,15451.9,27.5065], Tmin=(100,'K'), Tmax=(976.247,'K')), NASAPolynomial(coeffs=[15.5,0.0226346,-7.95269e-06,1.42521e-09,-1.00662e-13,11362.8,-48.8284], Tmin=(976.247,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(127.47,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = 'C=CC1COC1=C(20223)',
    structure = SMILES('C=CC1COC1=C'),
    E0 = (16.4644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14634,0.0453971,2.40144e-05,-7.42331e-08,3.55609e-11,2099.15,17.4914], Tmin=(100,'K'), Tmax=(910.951,'K')), NASAPolynomial(coeffs=[18.358,0.015608,-2.32912e-06,2.22858e-10,-1.53313e-14,-2936.44,-74.3671], Tmin=(910.951,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(16.4644,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(2methyleneoxetane)"""),
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
    label = 'C=C[CH]C(=C)[O](20224)',
    structure = SMILES('[CH2]C=CC(=C)[O]'),
    E0 = (99.8633,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,301.81,301.945],'cm^-1')),
        HinderedRotor(inertia=(0.350061,'amu*angstrom^2'), symmetry=1, barrier=(22.5346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15859,'amu*angstrom^2'), symmetry=1, barrier=(74.5656,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60396,0.0427771,-7.06638e-06,-2.71115e-08,1.55244e-11,12106.2,20.7034], Tmin=(100,'K'), Tmax=(925.164,'K')), NASAPolynomial(coeffs=[13.839,0.0154312,-4.15972e-06,6.48319e-10,-4.42386e-14,8748.73,-43.2822], Tmin=(925.164,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(99.8633,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(C=C(C)OJ)"""),
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
    label = '[CH2]C([C]=C)C=C(20225)',
    structure = SMILES('[CH2]C([C]=C)C=C'),
    E0 = (498.104,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.141251,'amu*angstrom^2'), symmetry=1, barrier=(3.24763,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141273,'amu*angstrom^2'), symmetry=1, barrier=(3.24815,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.816751,'amu*angstrom^2'), symmetry=1, barrier=(18.7787,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.1277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51674,0.0559432,-4.92707e-05,2.65105e-08,-6.18284e-12,59996.4,23.5842], Tmin=(100,'K'), Tmax=(1003.25,'K')), NASAPolynomial(coeffs=[7.6185,0.0316147,-1.28953e-05,2.33814e-09,-1.59194e-13,58772.1,-5.87296], Tmin=(1003.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(498.104,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1CC(=O)C1[CH2](20103)',
    structure = SMILES('[CH2]C1CC(=O)C1[CH2]'),
    E0 = (237.301,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39305,0.0470868,-3.20836e-06,-2.69215e-08,1.31925e-11,28643.6,25.6327], Tmin=(100,'K'), Tmax=(989.216,'K')), NASAPolynomial(coeffs=[11.7968,0.0284314,-1.04226e-05,1.86676e-09,-1.29788e-13,25439.8,-30.2366], Tmin=(989.216,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(237.301,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(Isobutyl) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2][C](C=C)C(C)=O(20226)',
    structure = SMILES('[CH2]C(C=C)=C(C)[O]'),
    E0 = (95.0504,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.542989,0.0649447,-3.5948e-05,-7.64674e-09,1.00431e-11,11566.5,24.6485], Tmin=(100,'K'), Tmax=(948.687,'K')), NASAPolynomial(coeffs=[17.7292,0.0192321,-5.96658e-06,1.00726e-09,-7.00649e-14,7101.82,-63.7043], Tmin=(948.687,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(95.0504,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([C]=C)C(C)=O(20227)',
    structure = SMILES('[CH2]C([C]=C)C(C)=O'),
    E0 = (263.485,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,375,552.5,462.5,1710,284.238,2341.07],'cm^-1')),
        HinderedRotor(inertia=(0.155044,'amu*angstrom^2'), symmetry=1, barrier=(8.96657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154777,'amu*angstrom^2'), symmetry=1, barrier=(8.94628,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222437,'amu*angstrom^2'), symmetry=1, barrier=(12.8741,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158969,'amu*angstrom^2'), symmetry=1, barrier=(8.98521,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.760449,0.0777026,-0.000106617,9.04763e-08,-3.15739e-11,31800.4,27.5938], Tmin=(100,'K'), Tmax=(803.249,'K')), NASAPolynomial(coeffs=[6.72388,0.0391904,-1.82369e-05,3.46058e-09,-2.38983e-13,31126.7,1.90045], Tmin=(803.249,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(263.485,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH]=CC([CH2])C(C)=O(20228)',
    structure = SMILES('[CH]=CC([CH2])C(C)=O'),
    E0 = (272.739,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,375,552.5,462.5,1710,262.594],'cm^-1')),
        HinderedRotor(inertia=(0.141168,'amu*angstrom^2'), symmetry=1, barrier=(6.90277,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141074,'amu*angstrom^2'), symmetry=1, barrier=(6.90243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141134,'amu*angstrom^2'), symmetry=1, barrier=(6.90202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141192,'amu*angstrom^2'), symmetry=1, barrier=(6.90272,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06776,0.0716101,-7.54963e-05,3.774e-08,-1.92934e-12,32901.9,26.4721], Tmin=(100,'K'), Tmax=(598.388,'K')), NASAPolynomial(coeffs=[7.79005,0.0375444,-1.73515e-05,3.31893e-09,-2.31939e-13,31902.8,-4.13324], Tmin=(598.388,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(272.739,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJC(C)C=O) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C]1OCC1C=C(20229)',
    structure = SMILES('[CH2][C]1OCC1C=C'),
    E0 = (330.721,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.655855,0.0652462,-5.85912e-05,3.05225e-08,-6.31634e-12,39904.1,25.7835], Tmin=(100,'K'), Tmax=(1320.74,'K')), NASAPolynomial(coeffs=[11.7339,0.0256369,-6.72531e-06,8.6922e-10,-4.5926e-14,37506.2,-28.7434], Tmin=(1320.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(330.721,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Oxetane) + radical(C2CsJOCs) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]C1[CH]CCC1=O(20132)',
    structure = SMILES('[CH2]C1[CH]CCC1=O'),
    E0 = (166.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62978,0.0377602,2.73279e-05,-5.93555e-08,2.44922e-11,20081.6,24.4922], Tmin=(100,'K'), Tmax=(985.074,'K')), NASAPolynomial(coeffs=[12.6278,0.0278108,-1.03747e-05,1.92937e-09,-1.38899e-13,16230.8,-36.9492], Tmin=(985.074,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(166.152,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclopentanone) + radical(CJC(C)C=O) + radical(CCJCC=O)"""),
)

species(
    label = 'C=CC(=C)C(C)=O(20230)',
    structure = SMILES('C=CC(=C)C(C)=O'),
    E0 = (-65.1872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.938352,0.0563536,-2.31691e-05,-1.02959e-08,7.74855e-12,-7720.49,23.1976], Tmin=(100,'K'), Tmax=(1040.31,'K')), NASAPolynomial(coeffs=[14.9234,0.0248995,-9.99642e-06,1.88496e-09,-1.34469e-13,-11837.9,-50.6293], Tmin=(1040.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-65.1872,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)HHH) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2]C(C=C)C[C]=O(19798)',
    structure = SMILES('[CH2]C(C=C)C[C]=O'),
    E0 = (210.189,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,242.617,3303.35],'cm^-1')),
        HinderedRotor(inertia=(0.0028642,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.285378,'amu*angstrom^2'), symmetry=1, barrier=(11.9226,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.89184,'amu*angstrom^2'), symmetry=1, barrier=(79.0237,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.285411,'amu*angstrom^2'), symmetry=1, barrier=(11.9226,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3598.66,'J/mol'), sigma=(6.12618,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=562.10 K, Pc=35.52 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10296,0.0627311,-5.47256e-05,2.73666e-08,-5.73558e-12,25385.1,29.3539], Tmin=(100,'K'), Tmax=(1127.43,'K')), NASAPolynomial(coeffs=[10.1343,0.0306887,-1.20943e-05,2.15793e-09,-1.45717e-13,23348.6,-15.3005], Tmin=(1127.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(210.189,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O) + radical(Isobutyl)"""),
)

species(
    label = 'C=CC1CCC1=O(19801)',
    structure = SMILES('C=CC1CCC1=O'),
    E0 = (-42.603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.84948,0.0343186,2.85891e-05,-5.46709e-08,2.1263e-11,-5035.21,22.6428], Tmin=(100,'K'), Tmax=(1016.74,'K')), NASAPolynomial(coeffs=[10.6265,0.0309465,-1.24038e-05,2.34864e-09,-1.68421e-13,-8430.48,-27.7668], Tmin=(1016.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-42.603,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutanone)"""),
)

species(
    label = '[CH2][C](C=C)OC=C(20090)',
    structure = SMILES('[CH2]C=C([CH2])OC=C'),
    E0 = (202.413,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,775.311,777.391],'cm^-1')),
        HinderedRotor(inertia=(0.0371832,'amu*angstrom^2'), symmetry=1, barrier=(16.1453,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.701779,'amu*angstrom^2'), symmetry=1, barrier=(16.1353,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.703367,'amu*angstrom^2'), symmetry=1, barrier=(16.1718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.359687,'amu*angstrom^2'), symmetry=1, barrier=(34.711,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.66664,0.0624699,-3.25005e-05,-7.56722e-09,8.8975e-12,24474.3,27.7839], Tmin=(100,'K'), Tmax=(971.847,'K')), NASAPolynomial(coeffs=[16.5687,0.0216335,-7.46292e-06,1.31887e-09,-9.24464e-14,20221.1,-54.4607], Tmin=(971.847,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(202.413,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(O)CJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([C]=O)C=C(20231)',
    structure = SMILES('[CH2]C([C]=O)C=C'),
    E0 = (239.203,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1855,455,950,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,345.456],'cm^-1')),
        HinderedRotor(inertia=(0.115952,'amu*angstrom^2'), symmetry=1, barrier=(9.79792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115863,'amu*angstrom^2'), symmetry=1, barrier=(9.79774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115995,'amu*angstrom^2'), symmetry=1, barrier=(9.79906,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69766,0.0540465,-5.77064e-05,3.67604e-08,-9.98038e-12,28849.4,23.3324], Tmin=(100,'K'), Tmax=(875.26,'K')), NASAPolynomial(coeffs=[7.45677,0.0277279,-1.26042e-05,2.40833e-09,-1.6882e-13,27841.2,-3.68488], Tmin=(875.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(239.203,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CJC(C)C=O) + radical(CC(C)CJ=O)"""),
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
    E0 = (183.509,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (220.746,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (341.02,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (305.552,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (315.221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (317.708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (278.586,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (257.234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (281.976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (368.147,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (340.619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (484.929,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (327.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (269.831,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (269.831,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (294.954,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (434.899,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (414.725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (307.773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (309.273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (239.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (246.909,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (343.444,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (343.444,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (191.793,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (481.426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (741.109,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (242.503,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (300.567,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (307.794,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (305.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (330.721,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (298.15,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (246.909,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (455.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (191.793,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (516.213,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (620.766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(C=C)C(=C)[O](19794)'],
    products = ['CH2CO(28)', 'butadiene13(2459)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(C=C)C(=C)[O](19794)'],
    products = ['[CH2]C1CC1C(=C)[O](20207)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.32e+08,'s^-1'), n=0.97, Ea=(37.2376,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 335 used for R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H
Exact match found for rate rule [R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(C=C)C(=C)[O](19794)'],
    products = ['[CH2]C1([O])CC1C=C(20208)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.24059e+09,'s^-1'), n=0.736667, Ea=(157.512,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 156.1 to 157.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(C=C)C(=C)[O](19794)'],
    products = ['[CH2]C1OC(=C)C1[CH2](20085)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 33 used for R5_SS_D;doublebond_intra_2H_pri;radadd_intra_O
Exact match found for rate rule [R5_SS_D;doublebond_intra_2H_pri;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C=CC(=C)C(=C)[O](20209)'],
    products = ['[CH2]C(C=C)C(=C)[O](19794)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.997e+08,'cm^3/(mol*s)'), n=1.629, Ea=(14.7235,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 64 used for Cds-CdCd_Cds-HH;HJ
Exact match found for rate rule [Cds-CdCd_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C2H3(30)', 'C=CC(=C)[O](2398)'],
    products = ['[CH2]C(C=C)C(=C)[O](19794)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(9050,'cm^3/(mol*s)'), n=2.41, Ea=(13.5143,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 828 used for Cds-CdH_Cds-HH;CdsJ-H
Exact match found for rate rule [Cds-CdH_Cds-HH;CdsJ-H]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=[C][O](173)', 'butadiene13(2459)'],
    products = ['[CH2]C(C=C)C(=C)[O](19794)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.0136943,'m^3/(mol*s)'), n=2.49, Ea=(21.946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CdH_Cds-HH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH2CO(28)', '[CH2][CH]C=C(2458)'],
    products = ['[CH2]C(C=C)C(=C)[O](19794)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7051.36,'m^3/(mol*s)'), n=1.07323, Ea=(43.3386,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cdd_Od;CJ] + [Ck_O;YJ] for rate rule [Ck_O;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(C=C)C(=C)[O](19794)'],
    products = ['C=C[C](C)C(=C)[O](20210)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(58.4615,'s^-1'), n=3.15787, Ea=(98.4673,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=[C]C(C)C(=C)[O](20211)'],
    products = ['[CH2]C(C=C)C(=C)[O](19794)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 204 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(C=C)C(=C)[O](19794)'],
    products = ['[CH2][C](C=C)C(=C)O(20212)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.19923e+08,'s^-1'), n=1.46351, Ea=(157.11,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;XH_out] + [R3H_SS_2Cd;Y_rad_out;XH_out] for rate rule [R3H_SS_2Cd;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=C(O)C([CH2])C=C(20213)'],
    products = ['[CH2]C(C=C)C(=C)[O](19794)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C([C]=C)C(=C)O(20214)'],
    products = ['[CH2]C(C=C)C(=C)[O](19794)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;XH_out] for rate rule [R4H_SS(Cd)S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=CC(C)C(=C)[O](20215)'],
    products = ['[CH2]C(C=C)C(=C)[O](19794)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=C([O])C(C)C=C(20216)'],
    products = ['[CH2]C(C=C)C(=C)[O](19794)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C(C=C)C(=C)[O](19794)'],
    products = ['[CH]=CC([CH2])C(=C)O(20217)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.468e+06,'s^-1','*|/',3), n=1.554, Ea=(111.445,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 290 used for R5H_SSSD;O_rad_out;Cd_H_out_singleH
Exact match found for rate rule [R5H_SSSD;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=[C][O](173)', '[CH2][CH]C=C(2458)'],
    products = ['[CH2]C(C=C)C(=C)[O](19794)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.49215e+07,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(C=C)C(=C)[O](19794)'],
    products = ['[CH2]C(C=C)[C]1CO1(20218)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(C=C)C(=C)[O](19794)'],
    products = ['C=C([O])C1[CH]CC1(20219)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.01e+08,'s^-1'), n=1.02, Ea=(124.265,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 18 used for R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(C=C)C(=C)[O](19794)'],
    products = ['[CH2][CH]C1CCC1=O(20120)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.51071e+08,'s^-1'), n=0.996667, Ea=(125.764,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(C=C)C(=C)[O](19794)'],
    products = ['[CH2]C1[CH]COC1=C(20220)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.07178e+07,'s^-1'), n=1.01592, Ea=(56.0365,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri_2H;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(C=C)C(=C)[O](19794)'],
    products = ['C=CC(=C)C(=C)O(20221)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C(C=C)C(=C)[O](19794)'],
    products = ['C=CC[CH]C(=C)[O](20222)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(C=C)C(=C)[O](19794)'],
    products = ['C=C[CH]CC(=C)[O](19793)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(C=C)C(=C)[O](19794)'],
    products = ['C=CC1COC1=C(20223)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CH2(19)', 'C=C[CH]C(=C)[O](20224)'],
    products = ['[CH2]C(C=C)C(=C)[O](19794)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/TwoDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['O(4)', '[CH2]C([C]=C)C=C(20225)'],
    products = ['[CH2]C(C=C)C(=C)[O](19794)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C(C=C)C(=C)[O](19794)'],
    products = ['[CH2]C1CC(=O)C1[CH2](20103)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.07e+06,'s^-1'), n=1.46, Ea=(58.9944,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 343 used for R5_SS_D;doublebond_intra_2H_pri;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C(C=C)C(=C)[O](19794)'],
    products = ['[CH2][C](C=C)C(C)=O(20226)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.30996,'s^-1'), n=3.5644, Ea=(117.059,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C([C]=C)C(C)=O(20227)'],
    products = ['[CH2]C(C=C)C(=C)[O](19794)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=CC([CH2])C(C)=O(20228)'],
    products = ['[CH2]C(C=C)C(=C)[O](19794)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5H;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C(C=C)C(=C)[O](19794)'],
    products = ['[CH2][C]1OCC1C=C(20229)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(3.27074e+08,'s^-1'), n=0.924088, Ea=(147.213,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_cs2H]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic
Ea raised from 146.2 to 147.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(C=C)C(=C)[O](19794)'],
    products = ['[CH2]C1[CH]CCC1=O(20132)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.12e+09,'s^-1'), n=0.63, Ea=(114.642,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 3 used for R5_SS_D;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C(C=C)C(=C)[O](19794)'],
    products = ['C=CC(=C)C(C)=O(20230)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(C=C)C[C]=O(19798)'],
    products = ['[CH2]C(C=C)C(=C)[O](19794)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C(C=C)C(=C)[O](19794)'],
    products = ['C=CC1CCC1=O(19801)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2][C](C=C)OC=C(20090)'],
    products = ['[CH2]C(C=C)C(=C)[O](19794)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction38',
    reactants = ['CH2(19)', '[CH2]C([C]=O)C=C(20231)'],
    products = ['[CH2]C(C=C)C(=C)[O](19794)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '3584',
    isomers = [
        '[CH2]C(C=C)C(=C)[O](19794)',
    ],
    reactants = [
        ('CH2CO(28)', 'butadiene13(2459)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '3584',
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

