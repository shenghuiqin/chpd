species(
    label = 'C=C[C]=CC[C]=O(26510)',
    structure = SMILES('C=C[C]=CC[C]=O'),
    E0 = (317.589,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,1855,455,950,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.955094,'amu*angstrom^2'), symmetry=1, barrier=(21.9595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.955335,'amu*angstrom^2'), symmetry=1, barrier=(21.965,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.955335,'amu*angstrom^2'), symmetry=1, barrier=(21.965,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05581,0.0559359,-3.84769e-05,9.5595e-09,1.66177e-13,38310.5,25.2911], Tmin=(100,'K'), Tmax=(1155.19,'K')), NASAPolynomial(coeffs=[14.6987,0.0209472,-8.95287e-06,1.70189e-09,-1.20199e-13,34341.1,-46.0341], Tmin=(1155.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(317.589,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(CCCJ=O)"""),
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
    label = 'CH2CHCCH(26391)',
    structure = SMILES('C#CC=C'),
    E0 = (274.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,750,770,3400,2100,3010,987.5,1337.5,450,1655,2175,525],'cm^-1')),
        HinderedRotor(inertia=(1.46338,'amu*angstrom^2'), symmetry=1, barrier=(33.6459,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (52.0746,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.18,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.87083,0.0182042,1.06711e-05,-2.72492e-08,1.19478e-11,33023.8,11.2934], Tmin=(100,'K'), Tmax=(955.249,'K')), NASAPolynomial(coeffs=[8.52653,0.0108962,-3.56564e-06,6.31243e-10,-4.51891e-14,31196.2,-19.6435], Tmin=(955.249,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""CH2CHCCH""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C1[C]=CCC1=O(27337)',
    structure = SMILES('[CH2]C1[C]=CCC1=O'),
    E0 = (352.996,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9973,0.031607,2.28017e-05,-4.60556e-08,1.75356e-11,42538.4,21.7913], Tmin=(100,'K'), Tmax=(1059.53,'K')), NASAPolynomial(coeffs=[11.3641,0.0261485,-1.18059e-05,2.35758e-09,-1.72938e-13,38875,-31.8607], Tmin=(1059.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.996,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopentane) + radical(CJC(C)C=O) + radical(cyclopentene-vinyl)"""),
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
    label = 'C=C[C]=CC=C=O(27338)',
    structure = SMILES('C=C[C]=CC=C=O'),
    E0 = (279.181,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2120,512.5,787.5,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.21817,'amu*angstrom^2'), symmetry=1, barrier=(28.0081,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21495,'amu*angstrom^2'), symmetry=1, barrier=(27.934,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30317,0.0602664,-6.51275e-05,3.93821e-08,-9.72249e-12,33674.1,22.2265], Tmin=(100,'K'), Tmax=(978.446,'K')), NASAPolynomial(coeffs=[10.1004,0.0243019,-9.99175e-06,1.81474e-09,-1.23641e-13,31952.6,-20.0231], Tmin=(978.446,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.181,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CCO)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = 'C=CC#CC[C]=O(27339)',
    structure = SMILES('C=CC#CC[C]=O'),
    E0 = (283.153,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2100,2250,500,550,1855,455,950,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,435.32,435.505],'cm^-1')),
        HinderedRotor(inertia=(0.163841,'amu*angstrom^2'), symmetry=1, barrier=(22.0592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163847,'amu*angstrom^2'), symmetry=1, barrier=(22.0553,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163903,'amu*angstrom^2'), symmetry=1, barrier=(22.0572,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7702,0.0362186,9.98354e-06,-3.94156e-08,1.72751e-11,34147,24.8855], Tmin=(100,'K'), Tmax=(1010.58,'K')), NASAPolynomial(coeffs=[13.9451,0.0180219,-7.52553e-06,1.503e-09,-1.12571e-13,30154.7,-41.5571], Tmin=(1010.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(283.153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CtHH) + group(Cds-CdsCtH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(CCCJ=O)"""),
)

species(
    label = 'C=C=C=CC[C]=O(27340)',
    structure = SMILES('C=C=C=CC[C]=O'),
    E0 = (347.979,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,540,563.333,586.667,610,1970,2140,1855,455,950,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.00594,'amu*angstrom^2'), symmetry=1, barrier=(23.1286,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00669,'amu*angstrom^2'), symmetry=1, barrier=(23.1458,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49048,0.0520775,-3.90628e-05,1.37612e-08,-1.93016e-12,41944.7,23.2179], Tmin=(100,'K'), Tmax=(1652.09,'K')), NASAPolynomial(coeffs=[14.8938,0.0196257,-9.59849e-06,1.87148e-09,-1.30967e-13,37516,-48.1743], Tmin=(1652.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(347.979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]=[C]C=C(4699)',
    structure = SMILES('[CH]=C=C[CH2]'),
    E0 = (451.584,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,3000,3100,440,815,1455,1000,180,1024.85,1025.53,1026.61],'cm^-1')),
        HinderedRotor(inertia=(0.00938781,'amu*angstrom^2'), symmetry=1, barrier=(7.01846,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (52.0746,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.76805,0.020302,8.75519e-06,-2.87666e-08,1.37354e-11,54363.7,13.5565], Tmin=(100,'K'), Tmax=(915.031,'K')), NASAPolynomial(coeffs=[9.46747,0.00887314,-1.78262e-06,2.38534e-10,-1.6263e-14,52390.1,-22.2544], Tmin=(915.031,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(451.584,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=C=CJ) + radical(Allyl_P)"""),
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
    label = 'C=C[C]=C[CH]C=O(27341)',
    structure = SMILES('[CH2]C=[C]C=CC=O'),
    E0 = (251.459,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32441,0.0561728,-4.32064e-05,1.72705e-08,-2.8472e-12,30341.8,23.6367], Tmin=(100,'K'), Tmax=(1405.2,'K')), NASAPolynomial(coeffs=[11.8733,0.0261443,-1.11519e-05,2.06282e-09,-1.41578e-13,27377.2,-30.8443], Tmin=(1405.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(251.459,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C=CC=[C]C[C]=O(27342)',
    structure = SMILES('C=CC=[C]C[C]=O'),
    E0 = (356.436,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,1855,455,950,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.739628,'amu*angstrom^2'), symmetry=1, barrier=(17.0055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.739627,'amu*angstrom^2'), symmetry=1, barrier=(17.0055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.739361,'amu*angstrom^2'), symmetry=1, barrier=(16.9994,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12965,0.0536258,-3.27617e-05,3.93017e-09,1.91995e-12,42980.6,25.9048], Tmin=(100,'K'), Tmax=(1149.15,'K')), NASAPolynomial(coeffs=[15.0779,0.0203227,-9.19487e-06,1.80531e-09,-1.29917e-13,38768,-47.7074], Tmin=(1149.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(356.436,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = 'C=[C]C=CC[C]=O(27343)',
    structure = SMILES('C=[C]C=CC[C]=O'),
    E0 = (317.589,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,1855,455,950,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.955094,'amu*angstrom^2'), symmetry=1, barrier=(21.9595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.955335,'amu*angstrom^2'), symmetry=1, barrier=(21.965,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.955335,'amu*angstrom^2'), symmetry=1, barrier=(21.965,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05581,0.0559359,-3.84769e-05,9.5595e-09,1.66177e-13,38310.5,25.2911], Tmin=(100,'K'), Tmax=(1155.19,'K')), NASAPolynomial(coeffs=[14.6987,0.0209472,-8.95287e-06,1.70189e-09,-1.20199e-13,34341.1,-46.0341], Tmin=(1155.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(317.589,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C[C]=[C]CC=O(27344)',
    structure = SMILES('[CH2][CH]C#CCC=O'),
    E0 = (353.523,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2100,2250,500,550,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,1155.21,1155.22],'cm^-1')),
        HinderedRotor(inertia=(0.00507669,'amu*angstrom^2'), symmetry=1, barrier=(4.80708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0191082,'amu*angstrom^2'), symmetry=1, barrier=(18.0957,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.38948,'amu*angstrom^2'), symmetry=1, barrier=(77.9308,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.38947,'amu*angstrom^2'), symmetry=1, barrier=(77.9306,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73336,0.0414384,-6.83762e-06,-1.83801e-08,9.48255e-12,42608,27.5467], Tmin=(100,'K'), Tmax=(1003.22,'K')), NASAPolynomial(coeffs=[10.6684,0.0244102,-9.1836e-06,1.65666e-09,-1.15163e-13,39879.4,-20.2529], Tmin=(1003.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(353.523,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsHH) + group(Cs-(Cds-O2d)CtHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Sec_Propargyl) + radical(RCCJ)"""),
)

species(
    label = 'C=CC=C[CH][C]=O(27345)',
    structure = SMILES('[CH2]C=CC=C[C]=O'),
    E0 = (213.085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00471,0.061723,-5.34234e-05,2.34529e-08,-4.14141e-12,25739.3,24.7144], Tmin=(100,'K'), Tmax=(1346.57,'K')), NASAPolynomial(coeffs=[14.202,0.0225201,-9.75334e-06,1.83241e-09,-1.27382e-13,22185.2,-42.8817], Tmin=(1346.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(213.085,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJ=O) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]=CC=CC[C]=O(27346)',
    structure = SMILES('[CH]=CC=CC[C]=O'),
    E0 = (365.69,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.817449,'amu*angstrom^2'), symmetry=1, barrier=(18.7948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.816641,'amu*angstrom^2'), symmetry=1, barrier=(18.7762,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.81777,'amu*angstrom^2'), symmetry=1, barrier=(18.8021,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11307,0.0523369,-2.4254e-05,-7.88984e-09,6.68817e-12,44095.8,25.8891], Tmin=(100,'K'), Tmax=(1065.98,'K')), NASAPolynomial(coeffs=[15.909,0.0189471,-8.41097e-06,1.6781e-09,-1.23469e-13,39684,-52.3358], Tmin=(1065.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(365.69,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CCCJ=O)"""),
)

species(
    label = 'C=[C][C]=CCC=O(27347)',
    structure = SMILES('[CH2]C#C[CH]CC=O'),
    E0 = (357.13,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2100,2250,500,550,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,686.45,2389.23],'cm^-1')),
        HinderedRotor(inertia=(3.25123,'amu*angstrom^2'), symmetry=1, barrier=(74.7521,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00884174,'amu*angstrom^2'), symmetry=1, barrier=(2.94065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128046,'amu*angstrom^2'), symmetry=1, barrier=(2.94403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.25058,'amu*angstrom^2'), symmetry=1, barrier=(74.7373,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73229,0.0504609,-3.92438e-05,1.76463e-08,-3.44754e-12,43033.9,27.3115], Tmin=(100,'K'), Tmax=(1169.55,'K')), NASAPolynomial(coeffs=[7.89691,0.0293771,-1.2203e-05,2.23246e-09,-1.52726e-13,41591.9,-3.39476], Tmin=(1169.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(357.13,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-CtCsHH) + group(Cs-CtHHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(CCJCC=O) + radical(Propargyl)"""),
)

species(
    label = '[CH]=C[C]=CCC=O(27348)',
    structure = SMILES('[CH]C=C=CCC=O'),
    E0 = (382.099,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000,239.496,239.499,239.499,239.5],'cm^-1')),
        HinderedRotor(inertia=(1.18909,'amu*angstrom^2'), symmetry=1, barrier=(48.3996,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18908,'amu*angstrom^2'), symmetry=1, barrier=(48.3997,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1891,'amu*angstrom^2'), symmetry=1, barrier=(48.3997,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2771,0.0519924,-2.23508e-05,-1.13599e-09,1.97311e-12,46060.2,25.2427], Tmin=(100,'K'), Tmax=(1354.51,'K')), NASAPolynomial(coeffs=[13.4964,0.0303059,-1.42796e-05,2.73927e-09,-1.90599e-13,41429.2,-42.2917], Tmin=(1354.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(382.099,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'O=[C]CC=C1[CH]C1(27349)',
    structure = SMILES('O=[C]CC=C1[CH]C1'),
    E0 = (352.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78483,0.0337222,2.56012e-05,-5.50898e-08,2.20602e-11,42516.4,22.1487], Tmin=(100,'K'), Tmax=(1026.34,'K')), NASAPolynomial(coeffs=[13.9995,0.0220305,-9.79867e-06,1.99806e-09,-1.50228e-13,38117.6,-46.3126], Tmin=(1026.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Methylene_cyclopropane) + radical(CCCJ=O) + radical(Allyl_S)"""),
)

species(
    label = 'O=C1C[CH][C]=CC1(27350)',
    structure = SMILES('O=C1C[CH][C]=CC1'),
    E0 = (303.204,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.52911,0.0156344,7.11472e-05,-9.02514e-08,3.06501e-11,36534.1,17.9166], Tmin=(100,'K'), Tmax=(1076.02,'K')), NASAPolynomial(coeffs=[10.0723,0.0331123,-1.66713e-05,3.47164e-09,-2.59217e-13,32275.7,-31.2723], Tmin=(1076.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.204,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclohexane) + radical(CCJCC=O) + radical(Cds_S)"""),
)

species(
    label = 'C=CC=CC=C=O(27351)',
    structure = SMILES('C=CC=CC=C=O'),
    E0 = (80.186,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11369,0.0568433,-4.67069e-05,1.98587e-08,-3.39265e-12,9753.53,23.3106], Tmin=(100,'K'), Tmax=(1397.17,'K')), NASAPolynomial(coeffs=[13.6724,0.0208885,-8.1057e-06,1.43979e-09,-9.68877e-14,6244.24,-41.4778], Tmin=(1397.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(80.186,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CCO)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C=C=CCC=O(27352)',
    structure = SMILES('C=C=C=CCC=O'),
    E0 = (188.019,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31129,0.0517125,-3.15697e-05,7.39533e-09,-3.65933e-13,22716.3,23.0696], Tmin=(100,'K'), Tmax=(1497.71,'K')), NASAPolynomial(coeffs=[15.945,0.0216448,-1.04848e-05,2.02886e-09,-1.40993e-13,17321.8,-56.816], Tmin=(1497.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.019,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=C[C]=CC(=C)[O](26502)',
    structure = SMILES('C=C[C]=CC(=C)[O]'),
    E0 = (268.601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.118517,0.0727959,-7.66533e-05,4.02971e-08,-8.02712e-12,32455.7,24.9008], Tmin=(100,'K'), Tmax=(1380.29,'K')), NASAPolynomial(coeffs=[18.2599,0.0122707,-2.23671e-06,1.80618e-10,-5.15332e-15,28205.2,-65.7238], Tmin=(1380.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(268.601,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=CC1=CCC1=O(27353)',
    structure = SMILES('C=CC1=CCC1=O'),
    E0 = (100.348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73576,0.0293734,5.00581e-05,-8.70655e-08,3.45668e-11,12168.7,18.1994], Tmin=(100,'K'), Tmax=(1001.03,'K')), NASAPolynomial(coeffs=[17.2246,0.0177641,-7.89189e-06,1.70728e-09,-1.35515e-13,6548.43,-69.1245], Tmin=(1001.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(100.348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutene)"""),
)

species(
    label = 'CO(12)',
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
    label = '[CH2]C=[C]C=C(26613)',
    structure = SMILES('[CH2]C=[C]C=C'),
    E0 = (374.947,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(2.07984,'amu*angstrom^2'), symmetry=1, barrier=(47.8195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.07746,'amu*angstrom^2'), symmetry=1, barrier=(47.7649,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22509,0.0302876,1.05277e-05,-3.68274e-08,1.72727e-11,45167.7,17.3269], Tmin=(100,'K'), Tmax=(923.516,'K')), NASAPolynomial(coeffs=[10.3879,0.0174192,-5.09531e-06,8.16549e-10,-5.51179e-14,42701,-26.5965], Tmin=(923.516,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(374.947,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CC=CCJ)"""),
)

species(
    label = '[C]=O(361)',
    structure = SMILES('[C]=O'),
    E0 = (439.086,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3054.48],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0101,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.08916,0.00200416,-1.61661e-05,2.55058e-08,-1.16424e-11,52802.7,4.52505], Tmin=(100,'K'), Tmax=(856.11,'K')), NASAPolynomial(coeffs=[0.961625,0.00569045,-3.48044e-06,7.19202e-10,-5.08041e-14,53738.7,21.4663], Tmin=(856.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(439.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH2]C=[C]C1CC1=O(27354)',
    structure = SMILES('[CH2]C=[C]C1CC1=O'),
    E0 = (428.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13475,0.0563999,-4.483e-05,1.8671e-08,-3.14749e-12,51638,23.5458], Tmin=(100,'K'), Tmax=(1408.61,'K')), NASAPolynomial(coeffs=[13.1323,0.0223308,-8.55051e-06,1.50067e-09,-1.00101e-13,48258,-38.4458], Tmin=(1408.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(428.44,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(cyclopropanone) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = 'C[C]=C=CC[C]=O(27355)',
    structure = SMILES('CC#C[CH]C[C]=O'),
    E0 = (360.683,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2100,2250,500,550,1855,455,950,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,370.937,1670.8],'cm^-1')),
        HinderedRotor(inertia=(0.0753793,'amu*angstrom^2'), symmetry=1, barrier=(7.37364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0756653,'amu*angstrom^2'), symmetry=1, barrier=(7.36192,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0749113,'amu*angstrom^2'), symmetry=1, barrier=(7.34111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.683237,'amu*angstrom^2'), symmetry=1, barrier=(66.9721,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58771,0.0539532,-4.80413e-05,2.14957e-08,-2.08723e-12,43466.4,27.6048], Tmin=(100,'K'), Tmax=(760.348,'K')), NASAPolynomial(coeffs=[8.3434,0.0264264,-9.54532e-06,1.60346e-09,-1.04041e-13,42207.5,-4.65988], Tmin=(760.348,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(360.683,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-CtCsHH) + group(Cs-CtHHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(CCJCC=O) + radical(CCCJ=O)"""),
)

species(
    label = 'CC=C=[C]C[C]=O(27356)',
    structure = SMILES('C[CH]C#CC[C]=O'),
    E0 = (308.237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64864,0.043381,-1.04734e-05,-1.75186e-08,1.01404e-11,37164.5,26.1892], Tmin=(100,'K'), Tmax=(967.743,'K')), NASAPolynomial(coeffs=[11.4698,0.0222105,-7.76553e-06,1.35592e-09,-9.33084e-14,34354.1,-25.5694], Tmin=(967.743,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.237,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsHH) + group(Cs-(Cds-O2d)CtHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Sec_Propargyl) + radical(CCCJ=O)"""),
)

species(
    label = 'CC=C=C[CH][C]=O(27357)',
    structure = SMILES('CC=[C]C=C[C]=O'),
    E0 = (294.024,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.818069,0.0750475,-9.80814e-05,7.23384e-08,-2.18689e-11,35473,23.4323], Tmin=(100,'K'), Tmax=(802.505,'K')), NASAPolynomial(coeffs=[9.80552,0.0302476,-1.43385e-05,2.76601e-09,-1.93948e-13,34030.6,-17.9489], Tmin=(802.505,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(294.024,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJ=O) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C=C1[CH]CC1=O(27358)',
    structure = SMILES('C=C[C]1[CH]CC1=O'),
    E0 = (277,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.52443,0.0124772,8.53759e-05,-1.16812e-07,4.44728e-11,33386.5,22.8662], Tmin=(100,'K'), Tmax=(966.722,'K')), NASAPolynomial(coeffs=[12.4599,0.0219279,-7.73937e-06,1.50285e-09,-1.15021e-13,29103,-36.9498], Tmin=(966.722,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(277,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutanone) + radical(C=CCJ(C)C=O) + radical(CCJCC=O)"""),
)

species(
    label = 'O=[C]CC1[C]=CC1(27331)',
    structure = SMILES('O=[C]CC1[C]=CC1'),
    E0 = (418.582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61671,0.0440599,-1.43501e-05,-1.06594e-08,6.57939e-12,50436.8,24.7313], Tmin=(100,'K'), Tmax=(1044.25,'K')), NASAPolynomial(coeffs=[11.4231,0.0235167,-9.28967e-06,1.71835e-09,-1.20668e-13,47460.7,-27.4469], Tmin=(1044.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(418.582,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclobutene) + radical(CCCJ=O) + radical(cyclobutene-vinyl)"""),
)

species(
    label = 'CC=C=CC=C=O(27359)',
    structure = SMILES('CC=C=CC=C=O'),
    E0 = (132.967,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54301,0.0584951,-6.24553e-05,4.26923e-08,-1.28723e-11,16076.7,21.8083], Tmin=(100,'K'), Tmax=(782.543,'K')), NASAPolynomial(coeffs=[6.28842,0.0342383,-1.59581e-05,3.07939e-09,-2.16845e-13,15334,0.0781844], Tmin=(782.543,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(132.967,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cd-Cd(CCO)H) + group(Cds-(Cdd-O2d)CsH) + group(Cdd-CdsCds)"""),
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
    label = 'C#C[CH]C[C]=O(27336)',
    structure = SMILES('[CH]=C=CC[C]=O'),
    E0 = (361.877,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,540,610,2055],'cm^-1')),
        HinderedRotor(inertia=(0.942019,'amu*angstrom^2'), symmetry=1, barrier=(21.6589,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.942967,'amu*angstrom^2'), symmetry=1, barrier=(21.6807,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83042,0.0411588,-2.90855e-05,8.37138e-09,-5.06207e-13,43607.3,21.7665], Tmin=(100,'K'), Tmax=(1230.55,'K')), NASAPolynomial(coeffs=[12.1824,0.0151463,-6.68669e-06,1.2802e-09,-9.02369e-14,40481.4,-32.6727], Tmin=(1230.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.877,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(CCCJ=O)"""),
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
    E0 = (317.589,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (393.811,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (501.966,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (507.204,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (576.453,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (436.4,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (444.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (435.27,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (552.895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (532.053,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (495.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (510.667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (588.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (508.789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (533.62,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (611.769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (455.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (378.476,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (406.558,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (342.562,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (563.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (325.874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (317.589,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (814.032,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (428.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (575.151,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (555.709,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (393.399,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (443.825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (436.071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (342.562,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (743.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C[C]=CC[C]=O(26510)'],
    products = ['CH2CO(28)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C[C]=CC[C]=O(26510)'],
    products = ['[CH2]C1[C]=CCC1=O(27337)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.70833e+09,'s^-1'), n=0.576561, Ea=(76.2219,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;doublebond_intra_2H_pri;radadd_intra] for rate rule [R6;doublebond_intra_2H_pri;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'C=C[C]=CC=C=O(27338)'],
    products = ['C=C[C]=CC[C]=O(26510)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.82e-16,'cm^3/(molecule*s)'), n=1.61, Ea=(10.992,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Ck;HJ] for rate rule [Cds-OneDeH_Ck;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C=CC#CC[C]=O(27339)'],
    products = ['C=C[C]=CC[C]=O(26510)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7.66e+08,'cm^3/(mol*s)'), n=1.64, Ea=(12.2591,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2703 used for Ct-Cs_Ct-Cd;HJ
Exact match found for rate rule [Ct-Cs_Ct-Cd;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C=C=C=CC[C]=O(27340)'],
    products = ['C=C[C]=CC[C]=O(26510)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(6.23e+08,'cm^3/(mol*s)'), n=1.429, Ea=(16.6816,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 193 used for Ca_Ca;HJ
Exact match found for rate rule [Ca_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2CO(28)', '[CH]=[C]C=C(4699)'],
    products = ['C=C[C]=CC[C]=O(26510)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.284303,'m^3/(mol*s)'), n=1.93802, Ea=(45.6341,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Ck;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=[C][O](173)', 'CH2CHCCH(26391)'],
    products = ['C=C[C]=CC[C]=O(26510)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.0117197,'m^3/(mol*s)'), n=2.33233, Ea=(9.74423,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-H_Ct-Cd;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=C[C]=CC[C]=O(26510)'],
    products = ['C=C[C]=C[CH]C=O(27341)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.27137e+08,'s^-1'), n=1.53496, Ea=(117.681,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/OneDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=CC=[C]C[C]=O(27342)'],
    products = ['C=C[C]=CC[C]=O(26510)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.75e+11,'s^-1'), n=0.633, Ea=(196.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=[C]C=CC[C]=O(27343)'],
    products = ['C=C[C]=CC[C]=O(26510)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.806e+09,'s^-1'), n=1.172, Ea=(214.463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=C[C]=[C]CC=O(27344)'],
    products = ['C=C[C]=CC[C]=O(26510)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Cd_rad_out;XH_out] for rate rule [R3H_SS_Cs;Cd_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=C[C]=CC[C]=O(26510)'],
    products = ['C=CC=C[CH][C]=O(27345)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.56742e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=CC=CC[C]=O(27346)'],
    products = ['C=C[C]=CC[C]=O(26510)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(191.5,'s^-1'), n=3.05, Ea=(222.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=[C][C]=CCC=O(27347)'],
    products = ['C=C[C]=CC[C]=O(26510)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.81182e+08,'s^-1'), n=1.25566, Ea=(151.659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;XH_out] for rate rule [R5HJ_1;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=C[C]=CCC=O(27348)'],
    products = ['C=C[C]=CC[C]=O(26510)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R6HJ_2;Cd_rad_out_singleH;CO_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=[C][O](173)', '[CH]=[C]C=C(4699)'],
    products = ['C=C[C]=CC[C]=O(26510)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=C[C]=CC[C]=O(26510)'],
    products = ['O=[C]CC=C1[CH]C1(27349)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(5.605e+12,'s^-1'), n=0.275, Ea=(137.649,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 133 used for R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble
Exact match found for rate rule [R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=C[C]=CC[C]=O(26510)'],
    products = ['O=C1C[CH][C]=CC1(27350)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(9.91671e+09,'s^-1'), n=0.30082, Ea=(60.8864,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra_pri_2H;radadd_intra] for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=C[C]=CC[C]=O(26510)'],
    products = ['C=CC=CC=C=O(27351)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C[C]=CC[C]=O(26510)'],
    products = ['C=C=C=CCC=O(27352)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=C[C]=CC[C]=O(26510)'],
    products = ['C=C[C]=CC(=C)[O](26502)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=C[C]=CC[C]=O(26510)'],
    products = ['C=CC1=CCC1=O(27353)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;Y_rad_out;CdsinglepriDe_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CO(12)', '[CH2]C=[C]C=C(26613)'],
    products = ['C=C[C]=CC[C]=O(26510)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2461.18,'m^3/(mol*s)'), n=1.0523, Ea=(61.8613,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [COm;C_pri_rad] for rate rule [COm;C_rad/H2/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_COm
Ea raised from 57.2 to 61.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction24',
    reactants = ['[C]=O(361)', '[CH2]C=[C]C=C(26613)'],
    products = ['C=C[C]=CC[C]=O(26510)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.13464e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=C[C]=CC[C]=O(26510)'],
    products = ['[CH2]C=[C]C1CC1=O(27354)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.64245e+09,'s^-1'), n=0.690807, Ea=(110.85,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 110.8 to 110.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction26',
    reactants = ['C[C]=C=CC[C]=O(27355)'],
    products = ['C=C[C]=CC[C]=O(26510)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(606700,'s^-1'), n=2.347, Ea=(214.468,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 119 used for R2H_S;Cd_rad_out_double;Cs_H_out_2H
Exact match found for rate rule [R2H_S;Cd_rad_out_double;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=C[C]=CC[C]=O(26510)'],
    products = ['CC=C=[C]C[C]=O(27356)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(0.0029655,'s^-1'), n=4.271, Ea=(238.12,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SMM;C_rad_out_2H;Cd_H_out_single] for rate rule [R4H_SMM;C_rad_out_2H;Cd_H_out_singleNd]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=C[C]=CC[C]=O(26510)'],
    products = ['CC=C=C[CH][C]=O(27357)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(4.42353e+06,'s^-1'), n=1.41298, Ea=(75.8094,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;C_rad_out_2H;XH_out] for rate rule [R5H_SMMS;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=C[C]=CC[C]=O(26510)'],
    products = ['[CH2]C=C1[CH]CC1=O(27358)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=C[C]=CC[C]=O(26510)'],
    products = ['O=[C]CC1[C]=CC1(27331)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.53664e+10,'s^-1'), n=0.43543, Ea=(118.482,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra_cs2H] + [R4;doublebond_intra;radadd_intra_cs2H] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=C[C]=CC[C]=O(26510)'],
    products = ['CC=C=CC=C=O(27359)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['CH2(19)', 'C#C[CH]C[C]=O(27336)'],
    products = ['C=C[C]=CC[C]=O(26510)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4478',
    isomers = [
        'C=C[C]=CC[C]=O(26510)',
    ],
    reactants = [
        ('CH2CO(28)', 'CH2CHCCH(26391)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4478',
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

