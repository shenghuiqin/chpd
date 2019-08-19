species(
    label = 'C=[C]OO[CH]C=O(1103)',
    structure = SMILES('C=[C]OO[CH]C=O'),
    E0 = (231.868,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2950,3100,1380,975,1025,1650,1685,370,2782.5,750,1395,475,1775,1000,3025,407.5,1350,352.5,234.811],'cm^-1')),
        HinderedRotor(inertia=(1.05724,'amu*angstrom^2'), symmetry=1, barrier=(41.3174,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.355602,'amu*angstrom^2'), symmetry=1, barrier=(13.9074,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00968806,'amu*angstrom^2'), symmetry=1, barrier=(41.3148,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05715,'amu*angstrom^2'), symmetry=1, barrier=(41.3139,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5601,0.0587794,-6.96838e-05,4.85987e-08,-1.44478e-11,27970.6,27.217], Tmin=(100,'K'), Tmax=(801.82,'K')), NASAPolynomial(coeffs=[7.44656,0.0294153,-1.4754e-05,2.93004e-09,-2.09474e-13,27026.6,0.118076], Tmin=(801.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(231.868,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(OCJC=O) + radical(C=CJO)"""),
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
    label = 'OCHCHO(48)',
    structure = SMILES('O=CC=O'),
    E0 = (-225.3,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,180,857.339,857.946,861.419,1717.87],'cm^-1')),
        HinderedRotor(inertia=(0.00964493,'amu*angstrom^2'), symmetry=1, barrier=(5.0669,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3660.03,'J/mol'), sigma=(4.01,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.68412,0.000478013,4.26391e-05,-5.79018e-08,2.31669e-11,-27198.5,4.51187], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[8.72507,0.00633097,-2.35575e-06,3.89783e-10,-2.37487e-14,-29102.4,-20.3904], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-225.3,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""OCHCHO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C=C1OO[CH]C1[O](4722)',
    structure = SMILES('C=C1OO[CH]C1[O]'),
    E0 = (223.757,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.12632,0.0246914,4.25161e-05,-7.47392e-08,3.03808e-11,26994,22.9383], Tmin=(100,'K'), Tmax=(975.306,'K')), NASAPolynomial(coeffs=[14.4515,0.0150467,-5.56076e-06,1.12558e-09,-8.82994e-14,22644.4,-46.1892], Tmin=(975.306,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(223.757,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(CCsJOOC) + radical(CC(C)OJ)"""),
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
    label = 'C=[C]OOC=C=O(4723)',
    structure = SMILES('C=[C]OOC=C=O'),
    E0 = (275.573,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2120,512.5,787.5,1685,370,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,194.675],'cm^-1')),
        HinderedRotor(inertia=(0.491937,'amu*angstrom^2'), symmetry=1, barrier=(12.0952,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.437896,'amu*angstrom^2'), symmetry=1, barrier=(12.0797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.966247,'amu*angstrom^2'), symmetry=1, barrier=(26.397,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5497,0.056545,-6.88164e-05,4.43977e-08,-1.15428e-11,33229.7,26.6256], Tmin=(100,'K'), Tmax=(933.362,'K')), NASAPolynomial(coeffs=[10.214,0.0194126,-9.13971e-06,1.77181e-09,-1.25229e-13,31612.4,-14.5769], Tmin=(933.362,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(275.573,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-(Cdd-O2d)OsH) + group(Cds-CdsHH) + radical(C=CJO)"""),
)

species(
    label = 'C#COO[CH]C=O(4724)',
    structure = SMILES('C#COO[CH]C=O'),
    E0 = (223.64,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,750,770,3400,2100,2175,525,2782.5,750,1395,475,1775,1000,3025,407.5,1350,352.5],'cm^-1')),
        HinderedRotor(inertia=(1.77213,'amu*angstrom^2'), symmetry=1, barrier=(40.7448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.7713,'amu*angstrom^2'), symmetry=1, barrier=(40.7257,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.7724,'amu*angstrom^2'), symmetry=1, barrier=(40.7509,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.77155,'amu*angstrom^2'), symmetry=1, barrier=(40.7313,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34613,0.0632934,-7.83954e-05,5.19206e-08,-1.40908e-11,26989.1,22.465], Tmin=(100,'K'), Tmax=(887.521,'K')), NASAPolynomial(coeffs=[9.97412,0.0244068,-1.26716e-05,2.55062e-09,-1.8377e-13,25457.6,-18.1303], Tmin=(887.521,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(223.64,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCt) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Ct-CtOs) + group(Ct-CtH) + radical(OCJC=O)"""),
)

species(
    label = 'C=[C]OOC[C]=O(4725)',
    structure = SMILES('C=[C]OOC[C]=O'),
    E0 = (248.054,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2950,3100,1380,975,1025,1650,1685,370,1855,455,950,2750,2850,1437.5,1250,1305,750,350,2404.85],'cm^-1')),
        HinderedRotor(inertia=(0.636488,'amu*angstrom^2'), symmetry=1, barrier=(14.6341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.63644,'amu*angstrom^2'), symmetry=1, barrier=(14.633,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.50285,'amu*angstrom^2'), symmetry=1, barrier=(34.5536,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.5029,'amu*angstrom^2'), symmetry=1, barrier=(34.5547,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70309,0.0553406,-5.98686e-05,3.672e-08,-9.59642e-12,29912.6,27.6261], Tmin=(100,'K'), Tmax=(903.17,'K')), NASAPolynomial(coeffs=[7.95643,0.027646,-1.38736e-05,2.7698e-09,-1.99055e-13,28783.1,-1.90591], Tmin=(903.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(248.054,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CsCJ=O) + radical(C=CJO)"""),
)

species(
    label = '[CH]=COO[CH]C=O(4726)',
    structure = SMILES('[CH]=COO[CH]C=O'),
    E0 = (239.22,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,2782.5,750,1395,475,1775,1000,3025,407.5,1350,352.5],'cm^-1')),
        HinderedRotor(inertia=(1.35389,'amu*angstrom^2'), symmetry=1, barrier=(31.1286,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.55003,'amu*angstrom^2'), symmetry=1, barrier=(58.6301,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07673,'amu*angstrom^2'), symmetry=1, barrier=(24.7561,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0777,'amu*angstrom^2'), symmetry=1, barrier=(24.7785,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21264,0.0614842,-6.35255e-05,3.33307e-08,-7.0227e-12,28871.8,25.7253], Tmin=(100,'K'), Tmax=(1140.52,'K')), NASAPolynomial(coeffs=[12.864,0.0206207,-9.78214e-06,1.91607e-09,-1.36656e-13,26214,-32.0179], Tmin=(1140.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(239.22,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(OCJC=O) + radical(Cds_P)"""),
)

species(
    label = 'C=COO[CH][C]=O(4727)',
    structure = SMILES('C=COO[CH][C]=O'),
    E0 = (146.795,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32955,0.0584303,-5.68161e-05,2.79984e-08,-5.56687e-12,17751.7,25.5796], Tmin=(100,'K'), Tmax=(1200.12,'K')), NASAPolynomial(coeffs=[12.5916,0.0208938,-9.90017e-06,1.93664e-09,-1.3788e-13,15048.5,-30.8079], Tmin=(1200.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(146.795,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(OCJC=O) + radical(CsCJ=O)"""),
)

species(
    label = '[CH]=[C]OOCC=O(4728)',
    structure = SMILES('[CH]=[C]OOCC=O'),
    E0 = (340.479,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,1685,370,2782.5,750,1395,475,1775,1000,3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350],'cm^-1')),
        HinderedRotor(inertia=(0.159867,'amu*angstrom^2'), symmetry=1, barrier=(7.49924,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.485958,'amu*angstrom^2'), symmetry=1, barrier=(22.3024,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161813,'amu*angstrom^2'), symmetry=1, barrier=(7.48799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.918761,'amu*angstrom^2'), symmetry=1, barrier=(42.9899,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53546,0.0590153,-6.8827e-05,4.50373e-08,-1.23432e-11,41034.8,27.9517], Tmin=(100,'K'), Tmax=(871.507,'K')), NASAPolynomial(coeffs=[8.54477,0.0268438,-1.34539e-05,2.67854e-09,-1.92004e-13,39813.1,-4.90006], Tmin=(871.507,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(340.479,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CJO)"""),
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
    label = '[O]C=C[O](2537)',
    structure = SMILES('[O]C=C[O]'),
    E0 = (-18.8461,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,714.947,715.113],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.7018,-0.0027881,5.18222e-05,-5.96735e-08,2.03775e-11,-2247.83,13.3237], Tmin=(100,'K'), Tmax=(1035.3,'K')), NASAPolynomial(coeffs=[6.65087,0.0113883,-5.76543e-06,1.26602e-09,-9.87778e-14,-4228.83,-7.62441], Tmin=(1035.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-18.8461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = 'C=[C]OOC1[CH]O1(4729)',
    structure = SMILES('C=[C]OOC1[CH]O1'),
    E0 = (333.027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.816308,0.0625144,-6.82723e-05,3.84748e-08,-8.28898e-12,40175,28.2996], Tmin=(100,'K'), Tmax=(1271.57,'K')), NASAPolynomial(coeffs=[14.3985,0.013451,-2.91891e-06,2.91336e-10,-1.12175e-14,37233.2,-38.4752], Tmin=(1271.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(333.027,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCsJO) + radical(C=CJO)"""),
)

species(
    label = 'C=C1O[CH][CH]OO1(4730)',
    structure = SMILES('[CH2][C]1OC=COO1'),
    E0 = (197.255,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39288,0.0336122,4.7592e-05,-1.09231e-07,5.15856e-11,23840.5,23.3539], Tmin=(100,'K'), Tmax=(901.462,'K')), NASAPolynomial(coeffs=[24.5309,-0.00567827,7.51018e-06,-1.59688e-09,1.06602e-13,17093.7,-100.157], Tmin=(901.462,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(197.255,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(124trioxene) + radical(Cs_P) + radical(CJCOOH)"""),
)

species(
    label = 'C=COOC=C=O(4731)',
    structure = SMILES('C=COOC=C=O'),
    E0 = (35.8286,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.901122,0.0600197,-6.04485e-05,2.98149e-08,-5.69256e-12,4427.63,25.8195], Tmin=(100,'K'), Tmax=(1287.7,'K')), NASAPolynomial(coeffs=[16.4115,0.0118395,-4.32492e-06,7.58598e-10,-5.14319e-14,433.083,-52.9312], Tmin=(1287.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(35.8286,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-(Cdd-O2d)OsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C1OOC1C=O(4732)',
    structure = SMILES('C=C1OOC1C=O'),
    E0 = (-32.4598,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2137,0.0224823,4.63855e-05,-7.61516e-08,2.98806e-11,-3824.79,22.538], Tmin=(100,'K'), Tmax=(996.209,'K')), NASAPolynomial(coeffs=[14.1016,0.0163387,-6.98469e-06,1.46985e-09,-1.14982e-13,-8257.06,-45.1268], Tmin=(996.209,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-32.4598,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Cyclobutane)"""),
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
    label = 'CHCHO(47)',
    structure = SMILES('[CH]=C[O]'),
    E0 = (245.848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.11,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06864,0.0187233,-1.21319e-05,-3.33727e-10,2.32882e-12,29739.4,14.7866], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.96288,0.00799899,-4.30606e-06,1.11076e-09,-1.11415e-13,28725.6,-5.17392], Tmin=(1000,'K'), Tmax=(3000,'K'))], Tmin=(200,'K'), Tmax=(3000,'K'), E0=(245.848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CHCHO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = '[O]O[CH]C=O(3604)',
    structure = SMILES('[O]O[CH]C=O'),
    E0 = (38.8345,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,492.5,1135,1000,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.304071,'amu*angstrom^2'), symmetry=1, barrier=(6.9912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17371,'amu*angstrom^2'), symmetry=1, barrier=(56.7408,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0355,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.48924,0.0359331,-5.02567e-05,4.15061e-08,-1.40856e-11,4722.53,16.7287], Tmin=(100,'K'), Tmax=(782.809,'K')), NASAPolynomial(coeffs=[5.94539,0.0159797,-7.6283e-06,1.46013e-09,-1.01349e-13,4251.69,1.34978], Tmin=(782.809,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(38.8345,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + radical(OCJC=O) + radical(ROOJ)"""),
)

species(
    label = 'C=C1OOC1[CH][O](4733)',
    structure = SMILES('C=C1OOC1[CH][O]'),
    E0 = (294.725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63481,0.0445433,-2.55373e-05,3.63073e-09,8.01796e-13,35538.5,24.8348], Tmin=(100,'K'), Tmax=(1284.88,'K')), NASAPolynomial(coeffs=[13.0685,0.019823,-9.37311e-06,1.83072e-09,-1.29584e-13,31702.7,-36.685], Tmin=(1284.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(294.725,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = 'C=[C]OOC=[C]O(4734)',
    structure = SMILES('C=[C]OOC=[C]O'),
    E0 = (333.841,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2950,3100,1380,975,1025,1650,1670,1700,300,440,3615,1277.5,1000,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(0.783455,'amu*angstrom^2'), symmetry=1, barrier=(18.0132,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.783028,'amu*angstrom^2'), symmetry=1, barrier=(18.0034,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.782722,'amu*angstrom^2'), symmetry=1, barrier=(17.9963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.77982,'amu*angstrom^2'), symmetry=1, barrier=(17.9296,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.801135,0.0647647,-7.08401e-05,3.5538e-08,-6.13768e-12,40272,31.6264], Tmin=(100,'K'), Tmax=(975.375,'K')), NASAPolynomial(coeffs=[16.7052,0.0101872,-3.27716e-06,5.47887e-10,-3.7171e-14,36663.2,-47.3005], Tmin=(975.375,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(333.841,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]OO[C]=CO(4735)',
    structure = SMILES('C=[C]OO[C]=CO'),
    E0 = (333.841,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2950,3100,1380,975,1025,1650,1670,1700,300,440,3615,1277.5,1000,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(0.783455,'amu*angstrom^2'), symmetry=1, barrier=(18.0132,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.783028,'amu*angstrom^2'), symmetry=1, barrier=(18.0034,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.782722,'amu*angstrom^2'), symmetry=1, barrier=(17.9963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.77982,'amu*angstrom^2'), symmetry=1, barrier=(17.9296,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.801135,0.0647647,-7.08401e-05,3.5538e-08,-6.13768e-12,40272,31.6264], Tmin=(100,'K'), Tmax=(975.375,'K')), NASAPolynomial(coeffs=[16.7052,0.0101872,-3.27716e-06,5.47887e-10,-3.7171e-14,36663.2,-47.3005], Tmin=(975.375,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(333.841,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = 'C=COO[C]=C[O](4736)',
    structure = SMILES('C=COO[C]=C[O]'),
    E0 = (235.56,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.01697,'amu*angstrom^2'), symmetry=1, barrier=(23.382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01615,'amu*angstrom^2'), symmetry=1, barrier=(23.3633,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01655,'amu*angstrom^2'), symmetry=1, barrier=(23.3726,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01776,0.0558536,-3.9885e-05,2.06423e-09,5.74134e-12,28447.4,29.1308], Tmin=(100,'K'), Tmax=(966.477,'K')), NASAPolynomial(coeffs=[17.5054,0.00989658,-3.13946e-06,5.71243e-10,-4.27322e-14,24219.8,-55.2341], Tmin=(966.477,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(235.56,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = '[CH]=[C]OOC=CO(4737)',
    structure = SMILES('[CH]=[C]OOC=CO'),
    E0 = (341.193,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,1685,370,3615,1277.5,1000,3120,650,792.5,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680],'cm^-1')),
        HinderedRotor(inertia=(0.993043,'amu*angstrom^2'), symmetry=1, barrier=(22.832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.99184,'amu*angstrom^2'), symmetry=1, barrier=(22.8044,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.990915,'amu*angstrom^2'), symmetry=1, barrier=(22.7831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.992811,'amu*angstrom^2'), symmetry=1, barrier=(22.8267,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.058187,0.0732709,-8.39975e-05,4.42405e-08,-8.60086e-12,41195.9,31.9877], Tmin=(100,'K'), Tmax=(1443.42,'K')), NASAPolynomial(coeffs=[21.7996,0.00188457,1.4251e-06,-4.03756e-10,3.05129e-14,36012.4,-77.5832], Tmin=(1443.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(341.193,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = 'C1=COO1(3609)',
    structure = SMILES('C1=COO1'),
    E0 = (143.913,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31094,-0.00195743,7.27275e-05,-1.01582e-07,4.10226e-11,17349.4,10.0685], Tmin=(100,'K'), Tmax=(927.669,'K')), NASAPolynomial(coeffs=[13.7806,-0.00309387,3.40676e-06,-6.26613e-10,3.47414e-14,13513.3,-49.8617], Tmin=(927.669,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(143.913,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cyclobutene)"""),
)

species(
    label = 'C=C1OC=COO1(4738)',
    structure = SMILES('C=C1OC=COO1'),
    E0 = (1.13085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.11061,0.0250517,4.40865e-05,-7.2269e-08,2.76582e-11,218.472,16.525], Tmin=(100,'K'), Tmax=(1019.75,'K')), NASAPolynomial(coeffs=[13.5705,0.0208295,-9.6139e-06,2.00544e-09,-1.53064e-13,-4236.45,-49.3697], Tmin=(1019.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1.13085,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclohexane)"""),
)

species(
    label = 'C=[C]OC([O])C=O(1102)',
    structure = SMILES('C=[C]OC([O])C=O'),
    E0 = (16.7303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10892,0.0566657,-5.58149e-05,2.71906e-08,-5.16051e-12,2122.1,29.0433], Tmin=(100,'K'), Tmax=(1289.64,'K')), NASAPolynomial(coeffs=[15.2173,0.0129067,-4.91817e-06,8.79998e-10,-6.01366e-14,-1516.84,-42.6099], Tmin=(1289.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(16.7303,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(C=CJO)"""),
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
    E0 = (231.868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (326.061,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (512.683,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (450.542,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (400.548,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (344.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (383.527,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (384.788,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (231.868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (415.323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (290.286,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (271.093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (239.776,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (584.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (440.036,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (306.514,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (231.868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (496.62,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (484.186,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (336.737,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (492.714,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (314.766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (239.399,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (545.668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (792.991,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C]OO[CH]C=O(1103)'],
    products = ['CH2CO(28)', 'OCHCHO(48)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=[C]OO[CH]C=O(1103)'],
    products = ['C=C1OO[CH]C1[O](4722)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.83276e+09,'s^-1'), n=0.690186, Ea=(94.193,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;multiplebond_intra;radadd_intra_cddouble] + [R6;multiplebond_intra;radadd_intra] for rate rule [R6;carbonylbond_intra_H;radadd_intra_cddouble]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'C=[C]OOC=C=O(4723)'],
    products = ['C=[C]OO[CH]C=O(1103)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.3e-15,'cm^3/(molecule*s)'), n=1.43, Ea=(25.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ck_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C#COO[CH]C=O(4724)'],
    products = ['C=[C]OO[CH]C=O(1103)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=[C]OOC[C]=O(4725)'],
    products = ['C=[C]OO[CH]C=O(1103)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.67823e+08,'s^-1'), n=1.48018, Ea=(152.494,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out] for rate rule [R2H_S;CO_rad_out;Cs_H_out_OOH/H]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=COO[CH]C=O(4726)'],
    products = ['C=[C]OO[CH]C=O(1103)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=[C]OO[CH]C=O(1103)'],
    products = ['C=COO[CH][C]=O(4727)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.81182e+08,'s^-1'), n=1.25566, Ea=(151.659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;XH_out] for rate rule [R5HJ_3;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=[C]OOCC=O(4728)'],
    products = ['C=[C]OO[CH]C=O(1103)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/OneDe] for rate rule [R5HJ_1;Cd_rad_out_singleH;Cs_H_out_H/CO]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=[C][O](173)', '[O]C=C[O](2537)'],
    products = ['C=[C]OO[CH]C=O(1103)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.49215e+07,'m^3/(mol*s)'), n=0.027223, Ea=(90.529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 90.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=[C]OO[CH]C=O(1103)'],
    products = ['C=[C]OOC1[CH]O1(4729)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(9.27213e+10,'s^-1'), n=0.543712, Ea=(183.455,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra_cs] for rate rule [R3_CO;carbonyl_intra_H;radadd_intra_csHO]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=[C]OO[CH]C=O(1103)'],
    products = ['C=C1O[CH][CH]OO1(4730)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.19745e+11,'s^-1'), n=0.440371, Ea=(58.418,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6plus;carbonyl_intra_H;radadd_intra] + [R6;multiplebond_intra;radadd_intra_cddouble] + [R6_linear;multiplebond_intra;radadd_intra] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_cddouble]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=[C]OO[CH]C=O(1103)'],
    products = ['C=COOC=C=O(4731)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=[C]OO[CH]C=O(1103)'],
    products = ['C=C1OOC1C=O(4732)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=[C]O[O](110)', 'CHCHO(47)'],
    products = ['C=[C]OO[CH]C=O(1103)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['H2CC(41)', '[O]O[CH]C=O(3604)'],
    products = ['C=[C]OO[CH]C=O(1103)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=[C]OO[CH]C=O(1103)'],
    products = ['C=C1OOC1[CH][O](4733)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(8.48517e+07,'s^-1'), n=1.03851, Ea=(74.6464,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra;radadd_intra] for rate rule [R5_SS_D;doublebond_intra;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CH2CO(28)', '[O]C=C[O](2537)'],
    products = ['C=[C]OO[CH]C=O(1103)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.6e+11,'cm^3/(mol*s)'), n=0, Ea=(311.532,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R_R;O_rad/OneDe] for rate rule [Od_Cdd;O_rad/OneDe]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 310.0 to 311.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=[C]OOC=[C]O(4734)'],
    products = ['C=[C]OO[CH]C=O(1103)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=[C]OO[C]=CO(4735)'],
    products = ['C=[C]OO[CH]C=O(1103)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleNd;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=COO[C]=C[O](4736)'],
    products = ['C=[C]OO[CH]C=O(1103)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cd_H_out_doubleC]
Euclidian distance = 2.82842712475
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=[C]OOC=CO(4737)'],
    products = ['C=[C]OO[CH]C=O(1103)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R7HJ_1;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=[C]OO[CH]C=O(1103)'],
    products = ['C=[C][O](173)', 'C1=COO1(3609)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.11355e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO;Y_rad_intra;OO_intra] for rate rule [R3OO_SD;Y_rad_intra;OO_intra]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=[C]OO[CH]C=O(1103)'],
    products = ['C=C1OC=COO1(4738)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Ypri_rad_out] for rate rule [R6_SSSDS;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=[C]OO[CH]C=O(1103)'],
    products = ['C=[C]OC([O])C=O(1102)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_R]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction25',
    reactants = ['O(4)', '[CH]=COO[C]=C(4571)'],
    products = ['C=[C]OO[CH]C=O(1103)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

network(
    label = '682',
    isomers = [
        'C=[C]OO[CH]C=O(1103)',
    ],
    reactants = [
        ('CH2CO(28)', 'OCHCHO(48)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '682',
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

