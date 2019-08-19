species(
    label = 'C=C[C]=CC=[C]O(29679)',
    structure = SMILES('C=C[C]=CC=[C]O'),
    E0 = (375.64,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1670,1700,300,440,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.19572,'amu*angstrom^2'), symmetry=1, barrier=(27.492,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1971,'amu*angstrom^2'), symmetry=1, barrier=(27.5237,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19726,'amu*angstrom^2'), symmetry=1, barrier=(27.5274,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.523776,0.0772253,-8.168e-05,4.12092e-08,-7.64319e-12,45361.1,29.4587], Tmin=(100,'K'), Tmax=(1566.19,'K')), NASAPolynomial(coeffs=[20.8943,0.00669209,1.03514e-06,-4.53709e-10,3.74056e-14,40593.9,-77.2814], Tmin=(1566.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(375.64,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CJO)"""),
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
    label = '[CH2]C1[C]=CC=C1O(30128)',
    structure = SMILES('[CH2]C1[C]=CC=C1O'),
    E0 = (350.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07219,0.0512362,-1.00603e-05,-3.48494e-08,2.06727e-11,42222.5,20.8627], Tmin=(100,'K'), Tmax=(918.367,'K')), NASAPolynomial(coeffs=[18.0459,0.0118138,-2.0331e-06,2.38547e-10,-1.70045e-14,37649.7,-67.503], Tmin=(918.367,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(350.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(Cyclopentadiene) + radical(1,3-cyclopentadiene-vinyl-1) + radical(Isobutyl)"""),
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
    label = 'C=CC#CC=[C]O(30131)',
    structure = SMILES('C=CC#CC=[C]O'),
    E0 = (358.111,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2100,2250,500,550,1685,370,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.32189,'amu*angstrom^2'), symmetry=1, barrier=(30.3928,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32168,'amu*angstrom^2'), symmetry=1, barrier=(30.3879,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.3209,'amu*angstrom^2'), symmetry=1, barrier=(30.3702,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12343,0.0476307,-1.79232e-06,-4.39938e-08,2.36026e-11,43188.7,25.7294], Tmin=(100,'K'), Tmax=(939.52,'K')), NASAPolynomial(coeffs=[19.6907,0.00793289,-1.24104e-06,1.97138e-10,-1.92541e-14,37963.1,-71.9319], Tmin=(939.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(358.111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-CdsCtH) + group(Cds-CdsCtH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-Ct(Cds-Cds)) + radical(C=CJO)"""),
)

species(
    label = 'C=C=C=CC=[C]O(30132)',
    structure = SMILES('C=C=C=CC=[C]O'),
    E0 = (406.03,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,563.333,586.667,610,1970,2140,1685,370,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(1.22376,'amu*angstrom^2'), symmetry=1, barrier=(28.1366,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.22264,'amu*angstrom^2'), symmetry=1, barrier=(28.1109,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0439694,0.0718654,-7.74079e-05,3.97672e-08,-7.63583e-12,48989.6,26.9047], Tmin=(100,'K'), Tmax=(1446.6,'K')), NASAPolynomial(coeffs=[20.2211,0.00687732,-4.85223e-07,-7.72261e-11,9.49411e-15,44114.1,-74.5621], Tmin=(1446.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(406.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=CJO)"""),
)

species(
    label = 'C=C[C]=CC#CO(30133)',
    structure = SMILES('C=C[C]=CC#CO'),
    E0 = (380.009,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2100,2250,500,550,1685,370,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.44998,'amu*angstrom^2'), symmetry=1, barrier=(33.3378,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.44958,'amu*angstrom^2'), symmetry=1, barrier=(33.3287,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.45497,'amu*angstrom^2'), symmetry=1, barrier=(33.4526,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03231,0.0573538,-4.38715e-05,1.07881e-08,1.35379e-12,45818.4,23.1105], Tmin=(100,'K'), Tmax=(1002.91,'K')), NASAPolynomial(coeffs=[15.1539,0.0165547,-6.06805e-06,1.09239e-09,-7.65122e-14,42205.1,-48.9512], Tmin=(1002.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(380.009,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CtH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtOs) + radical(C=CJC=C)"""),
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
    label = 'C=C[C]=C[CH]C=O(27341)',
    structure = SMILES('[CH2]C=[C]C=CC=O'),
    E0 = (251.459,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.48426,'amu*angstrom^2'), symmetry=1, barrier=(34.1261,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48843,'amu*angstrom^2'), symmetry=1, barrier=(34.222,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4825,'amu*angstrom^2'), symmetry=1, barrier=(34.0857,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32441,0.0561728,-4.32064e-05,1.72705e-08,-2.8472e-12,30341.8,23.6367], Tmin=(100,'K'), Tmax=(1405.2,'K')), NASAPolynomial(coeffs=[11.8733,0.0261443,-1.11519e-05,2.06282e-09,-1.41578e-13,27377.2,-30.8443], Tmin=(1405.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(251.459,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C=CC=[C]C=[C]O(30134)',
    structure = SMILES('C=CC=[C]C=[C]O'),
    E0 = (375.64,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1670,1700,300,440,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.19572,'amu*angstrom^2'), symmetry=1, barrier=(27.492,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1971,'amu*angstrom^2'), symmetry=1, barrier=(27.5237,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19726,'amu*angstrom^2'), symmetry=1, barrier=(27.5274,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.523776,0.0772253,-8.168e-05,4.12092e-08,-7.64319e-12,45361.1,29.4587], Tmin=(100,'K'), Tmax=(1566.19,'K')), NASAPolynomial(coeffs=[20.8943,0.00669209,1.03514e-06,-4.53709e-10,3.74056e-14,40593.9,-77.2814], Tmin=(1566.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(375.64,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(C=CJC=C)"""),
)

species(
    label = 'C=[C]C=CC=[C]O(30135)',
    structure = SMILES('C=C=C[CH]C=[C]O'),
    E0 = (362.005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12076,0.0484167,4.43351e-07,-4.68246e-08,2.50554e-11,43656.6,26.6762], Tmin=(100,'K'), Tmax=(922.347,'K')), NASAPolynomial(coeffs=[18.6637,0.011029,-1.67744e-06,1.89215e-10,-1.50551e-14,38774.7,-65.462], Tmin=(922.347,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(362.005,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJC=C) + radical(C=CJO)"""),
)

species(
    label = 'C=C[C]=C[C]=CO(30136)',
    structure = SMILES('[CH2]C=[C]C=C=CO'),
    E0 = (306.733,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.77765,0.076973,-7.72947e-05,3.72114e-08,-6.55991e-12,37087.8,28.8452], Tmin=(100,'K'), Tmax=(1676.06,'K')), NASAPolynomial(coeffs=[20.2867,0.00821761,7.81008e-07,-4.2385e-10,3.51927e-14,32623,-75.9114], Tmin=(1676.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(306.733,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=CJC=C) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=CC=C[C]=[C]O(30137)',
    structure = SMILES('[CH2]C=CC=C=[C]O'),
    E0 = (347.482,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.858064,0.0566006,-2.17625e-05,-2.47674e-08,1.75602e-11,41917.1,26.4848], Tmin=(100,'K'), Tmax=(916.104,'K')), NASAPolynomial(coeffs=[18.7897,0.0113914,-1.912e-06,2.10216e-10,-1.44509e-14,37243.3,-66.031], Tmin=(916.104,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(347.482,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=CC=CCJ) + radical(C=CJO)"""),
)

species(
    label = '[CH]=CC=CC=[C]O(30138)',
    structure = SMILES('[CH]=CC=CC=[C]O'),
    E0 = (423.741,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3120,650,792.5,1650,3615,1277.5,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180],'cm^-1')),
        HinderedRotor(inertia=(1.06023,'amu*angstrom^2'), symmetry=1, barrier=(24.3768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06016,'amu*angstrom^2'), symmetry=1, barrier=(24.3751,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05998,'amu*angstrom^2'), symmetry=1, barrier=(24.371,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.844186,0.0777569,-8.05003e-05,3.88447e-08,-6.82663e-12,51163.5,31.4347], Tmin=(100,'K'), Tmax=(1676.94,'K')), NASAPolynomial(coeffs=[22.0199,0.00446361,1.83624e-06,-5.57768e-10,4.18507e-14,46132.4,-82.8281], Tmin=(1676.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(423.741,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CJO)"""),
)

species(
    label = 'C=C[C]=[C]C=CO(30139)',
    structure = SMILES('C=C[C]=[C]C=CO'),
    E0 = (334.891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.03588,0.0847285,-9.37679e-05,4.86998e-08,-9.17956e-12,40481.9,27.7599], Tmin=(100,'K'), Tmax=(1568.67,'K')), NASAPolynomial(coeffs=[22.4198,0.00416798,3.10764e-06,-9.03567e-10,6.96364e-14,35676,-87.8232], Tmin=(1568.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(334.891,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = 'C=[C][C]=CC=CO(30140)',
    structure = SMILES('[CH2]C#CC=C[CH]O'),
    E0 = (313.067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35673,0.0446777,1.05004e-06,-3.68446e-08,1.79785e-11,37760.4,27.0122], Tmin=(100,'K'), Tmax=(982.142,'K')), NASAPolynomial(coeffs=[15.6974,0.017995,-6.64774e-06,1.26754e-09,-9.40377e-14,33413.5,-49.7037], Tmin=(982.142,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(313.067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(Propargyl) + radical(C=CCJO)"""),
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
    label = '[CH]=C[C]=CC=CO(30141)',
    structure = SMILES('[CH]C=C=CC=CO'),
    E0 = (360.366,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.372816,0.063616,-1.75721e-05,-3.8381e-08,2.40923e-11,43487.6,24.6486], Tmin=(100,'K'), Tmax=(919.096,'K')), NASAPolynomial(coeffs=[22.0509,0.0115717,-1.67082e-06,1.60879e-10,-1.23368e-14,37716.1,-87.8264], Tmin=(919.096,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(360.366,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'O[C]=CC=C1[CH]C1(30142)',
    structure = SMILES('O[C]=CC=C1[CH]C1'),
    E0 = (410.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12042,0.044253,2.0019e-05,-7.25032e-08,3.54031e-11,49526.9,23.0291], Tmin=(100,'K'), Tmax=(921.993,'K')), NASAPolynomial(coeffs=[21.194,0.00667641,6.02457e-07,-2.20098e-10,1.05189e-14,43720.9,-83.5961], Tmin=(921.993,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(410.78,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(Methylene_cyclopropane) + radical(Allyl_S) + radical(C=CJO)"""),
)

species(
    label = 'OC1=CC=[C][CH]C1(30143)',
    structure = SMILES('OC1=CC=[C][CH]C1'),
    E0 = (208.544,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36106,0.0473938,-1.33768e-05,-1.70129e-08,9.74619e-12,25186.3,19.4538], Tmin=(100,'K'), Tmax=(1025.05,'K')), NASAPolynomial(coeffs=[13.9611,0.0211404,-8.49152e-06,1.61864e-09,-1.16823e-13,21399.3,-47.5179], Tmin=(1025.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(208.544,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(Aromatic_pi_S_1_3) + radical(Cds_S)"""),
)

species(
    label = 'C=C=C=CC=CO(30144)',
    structure = SMILES('C=C=C=CC=CO'),
    E0 = (166.286,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.410386,0.0633551,-2.70708e-05,-2.92935e-08,2.14804e-11,20143.4,22.4587], Tmin=(100,'K'), Tmax=(914.767,'K')), NASAPolynomial(coeffs=[23.3448,0.00456596,1.2858e-06,-3.70401e-10,2.35308e-14,14211.3,-95.633], Tmin=(914.767,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(166.286,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds)"""),
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
    label = 'C=CC1=CC=C1O(30107)',
    structure = SMILES('C=CC1=CC=C1O'),
    E0 = (280.504,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25915,0.0458984,2.80428e-06,-4.42562e-08,2.24876e-11,33848.6,18.8814], Tmin=(100,'K'), Tmax=(944.332,'K')), NASAPolynomial(coeffs=[17.2,0.0142359,-3.86284e-06,6.62923e-10,-4.99032e-14,29239,-65.5771], Tmin=(944.332,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(280.504,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(cyclobutadiene_13)"""),
)

species(
    label = '[CH2]C=[C]C1C=C1O(30117)',
    structure = SMILES('[CH2]C=[C]C1C=C1O'),
    E0 = (489.614,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.255008,0.0697028,-6.97843e-05,3.45638e-08,-6.54912e-12,59032.6,24.9472], Tmin=(100,'K'), Tmax=(1403.18,'K')), NASAPolynomial(coeffs=[18.6451,0.0123648,-3.23718e-06,4.50751e-10,-2.66661e-14,54355.4,-68.2803], Tmin=(1403.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(489.614,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopropene) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = 'C[C]=C=CC=[C]O(30145)',
    structure = SMILES('CC#C[CH]C=[C]O'),
    E0 = (441.442,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2100,2250,500,550,2750,2800,2850,1350,1500,750,1050,1375,1000,454.671,454.954],'cm^-1')),
        HinderedRotor(inertia=(0.0771148,'amu*angstrom^2'), symmetry=1, barrier=(11.3204,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0769896,'amu*angstrom^2'), symmetry=1, barrier=(11.3232,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.077012,'amu*angstrom^2'), symmetry=1, barrier=(11.322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.562187,'amu*angstrom^2'), symmetry=1, barrier=(82.5865,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.930003,0.0588747,-5.40531e-05,2.64034e-08,-5.08236e-12,53211.1,29.9754], Tmin=(100,'K'), Tmax=(1334.25,'K')), NASAPolynomial(coeffs=[13.8823,0.0177954,-5.34214e-06,8.01285e-10,-4.8559e-14,49955,-35.497], Tmin=(1334.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(441.442,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtHH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(C=CJO) + radical(Cs_S)"""),
)

species(
    label = 'CC=C=[C]C=[C]O(30146)',
    structure = SMILES('C[CH]C#CC=[C]O'),
    E0 = (383.196,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2100,2250,500,550,2750,2800,2850,1350,1500,750,1050,1375,1000,249.474,249.815],'cm^-1')),
        HinderedRotor(inertia=(0.477069,'amu*angstrom^2'), symmetry=1, barrier=(21.0973,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.82886,'amu*angstrom^2'), symmetry=1, barrier=(80.5944,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.82057,'amu*angstrom^2'), symmetry=1, barrier=(80.6074,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.477129,'amu*angstrom^2'), symmetry=1, barrier=(21.1044,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0514922,0.0669386,-6.35873e-05,3.0144e-08,-5.31259e-12,46252.8,30.8324], Tmin=(100,'K'), Tmax=(1657.88,'K')), NASAPolynomial(coeffs=[16.9689,0.0118732,-1.0991e-06,-7.72182e-11,1.26897e-14,42533.2,-54.0832], Tmin=(1657.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(383.196,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-CdsCtH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(Sec_Propargyl) + radical(C=CJO)"""),
)

species(
    label = 'CC=C=C[C]=[C]O(30147)',
    structure = SMILES('CC=C=C[C]=[C]O'),
    E0 = (428.421,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.10857,'amu*angstrom^2'), symmetry=1, barrier=(25.4882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10831,'amu*angstrom^2'), symmetry=1, barrier=(25.4822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10888,'amu*angstrom^2'), symmetry=1, barrier=(25.4954,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.256659,0.0749862,-8.4744e-05,4.84158e-08,-1.06297e-11,51668.3,26.6773], Tmin=(100,'K'), Tmax=(1188.05,'K')), NASAPolynomial(coeffs=[17.2544,0.0142845,-3.71935e-06,4.89106e-10,-2.67899e-14,47874.5,-57.2246], Tmin=(1188.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(428.421,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=CJC=C) + radical(C=CJO)"""),
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
    label = 'O[C]=CC1[C]=CC1(30148)',
    structure = SMILES('O[C]=CC1[C]=CC1'),
    E0 = (499.607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1943,0.0507818,-1.80095e-05,-1.92641e-08,1.29172e-11,60199.7,24.4707], Tmin=(100,'K'), Tmax=(951.519,'K')), NASAPolynomial(coeffs=[15.7815,0.0159708,-4.92467e-06,8.49311e-10,-6.06093e-14,56223.6,-51.4856], Tmin=(951.519,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(499.607,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(cyclobutene-vinyl) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C=C1[CH]C=C1O(30149)',
    structure = SMILES('[CH2]C=C1[CH]C=C1O'),
    E0 = (189.586,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52392,0.0267887,8.27099e-05,-1.45254e-07,6.34887e-11,22917,19.6793], Tmin=(100,'K'), Tmax=(911.545,'K')), NASAPolynomial(coeffs=[23.6427,0.00203222,4.46739e-06,-1.01297e-09,6.38058e-14,15880.6,-101.46], Tmin=(911.545,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(189.586,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(C=CC=CCJ) + radical(C=CCJC=C)"""),
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
    label = '[CH]=C=CC=[C]O(30150)',
    structure = SMILES('[CH]=C=CC=[C]O'),
    E0 = (419.928,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,540,610,2055,3120,650,792.5,1650,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680],'cm^-1')),
        HinderedRotor(inertia=(1.23817,'amu*angstrom^2'), symmetry=1, barrier=(28.4679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24386,'amu*angstrom^2'), symmetry=1, barrier=(28.5988,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.364388,0.0612203,-6.84904e-05,3.57586e-08,-6.76316e-12,50652.7,25.519], Tmin=(100,'K'), Tmax=(1564.17,'K')), NASAPolynomial(coeffs=[17.4845,0.00229806,2.535e-06,-7.01958e-10,5.34472e-14,47149.3,-58.8135], Tmin=(1564.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(419.928,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CJO)"""),
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
    E0 = (375.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (388.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (521.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (583.406,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (634.504,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (606.911,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (632.242,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (554.309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (538.419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (587.903,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (590.104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (572.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (568.718,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (646.066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (525.986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (684.727,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (602.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (557.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (799.894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (513.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (386.888,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (384.008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (414.865,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (383.924,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (489.614,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (655.91,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (617.671,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (461.462,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (456.624,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (499.607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (513.695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (400.613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (801.491,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C[C]=CC=[C]O(29679)'],
    products = ['HCCOH(50)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C[C]=CC=[C]O(29679)'],
    products = ['[CH2]C1[C]=CC=C1O(30128)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(7.4226e+09,'s^-1'), n=0.3735, Ea=(13.1942,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;doublebond_intra_2H_pri;radadd_intra_cdsingle] for rate rule [R6;doublebond_intra_2H_pri;radadd_intra_cdsingleNd]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'C=C[C]=CC=C=O(27338)'],
    products = ['C=C[C]=CC=[C]O(29679)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.185e+08,'cm^3/(mol*s)'), n=1.63, Ea=(30.7064,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1700,'K'), comment="""Estimated using template [Od_R;HJ] for rate rule [Od_Cdd;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C=CC#CC=[C]O(30131)'],
    products = ['C=C[C]=CC=[C]O(29679)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2371.94,'m^3/(mol*s)'), n=1.49517, Ea=(13.5032,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-De_Ct-De;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C=C=C=CC=[C]O(30132)'],
    products = ['C=C[C]=CC=[C]O(29679)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(6.23e+08,'cm^3/(mol*s)'), n=1.429, Ea=(16.6816,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 193 used for Ca_Ca;HJ
Exact match found for rate rule [Ca_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', 'C=C[C]=CC#CO(30133)'],
    products = ['C=C[C]=CC=[C]O(29679)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=[C]O(172)', 'CH2CHCCH(26391)'],
    products = ['C=C[C]=CC=[C]O(29679)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.0117197,'m^3/(mol*s)'), n=2.33233, Ea=(9.74423,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-H_Ct-Cd;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['HCCOH(50)', '[CH]=[C]C=C(4699)'],
    products = ['C=C[C]=CC=[C]O(29679)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=C[C]=CC=[C]O(29679)'],
    products = ['C=C[C]=C[CH]C=O(27341)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=CC=[C]C=[C]O(30134)'],
    products = ['C=C[C]=CC=[C]O(29679)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.0945e+07,'s^-1'), n=1.951, Ea=(212.263,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_D;Cd_rad_out_singleDe_Cd;Cd_H_out_single] for rate rule [R2H_D;Cd_rad_out_singleDe_Cd;Cd_H_out_singleDe]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=C[C]=CC=[C]O(29679)'],
    products = ['C=[C]C=CC=[C]O(30135)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.806e+09,'s^-1'), n=1.172, Ea=(214.463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=C[C]=CC=[C]O(29679)'],
    products = ['C=C[C]=C[C]=CO(30136)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.75e+11,'s^-1'), n=0.633, Ea=(196.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_singleNd;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C[C]=CC=[C]O(29679)'],
    products = ['C=CC=C[C]=[C]O(30137)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.28371e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=CC=CC=[C]O(30138)'],
    products = ['C=C[C]=CC=[C]O(29679)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(191.5,'s^-1'), n=3.05, Ea=(222.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=C[C]=CC=[C]O(29679)'],
    products = ['C=C[C]=[C]C=CO(30139)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleNd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=C[C]=CC=[C]O(29679)'],
    products = ['C=[C][C]=CC=CO(30140)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.24929e+07,'s^-1'), n=1.719, Ea=(309.087,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_single;Cd_H_out_doubleC] for rate rule [R5HJ_3;Cd_rad_out_singleNd;Cd_H_out_doubleC]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=C[C]=CC=[C]O(29679)'],
    products = ['C=CC=C[CH][C]=O(27345)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.61756e+09,'s^-1'), n=1.29078, Ea=(226.405,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleDe_Cd;XH_out] for rate rule [R5HJ_3;Cd_rad_out_singleDe_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=C[C]=CC=[C]O(29679)'],
    products = ['[CH]=C[C]=CC=CO(30141)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.74e+08,'s^-1'), n=1.713, Ea=(181.895,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleNd;Cd_H_out_singleH] for rate rule [R6HJ_3;Cd_rad_out_singleNd;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=[C]O(172)', '[CH]=[C]C=C(4699)'],
    products = ['C=C[C]=CC=[C]O(29679)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C[C]=CC=[C]O(29679)'],
    products = ['O[C]=CC=C1[CH]C1(30142)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(5.605e+12,'s^-1'), n=0.275, Ea=(137.649,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 133 used for R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble
Exact match found for rate rule [R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=C[C]=CC=[C]O(29679)'],
    products = ['OC1=CC=[C][CH]C1(30143)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.76666e+15,'s^-1'), n=-0.518223, Ea=(11.2479,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6_linear;doublebond_intra_pri;radadd_intra_cdsingleNd] + [R6_linear;doublebond_intra_pri_2H;radadd_intra_cdsingle] for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_cdsingleNd]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=C[C]=CC=[C]O(29679)'],
    products = ['C=C=C=CC=CO(30144)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.21e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C[C]=CC=[C]O(29679)'],
    products = ['C=CC=CC=C=O(27351)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=C[C]=CC=[C]O(29679)'],
    products = ['C=CC1=CC=C1O(30107)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_DSD;CdsingleDe_rad_out;CdsinglepriND_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=C[C]=CC=[C]O(29679)'],
    products = ['[CH2]C=[C]C1C=C1O(30117)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.69117e+12,'s^-1'), n=0.1675, Ea=(113.974,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D_D;doublebond_intra;radadd_intra_cdsingle] for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleNd]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 113.6 to 114.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction26',
    reactants = ['C[C]=C=CC=[C]O(30145)'],
    products = ['C=C[C]=CC=[C]O(29679)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(606700,'s^-1'), n=2.347, Ea=(214.468,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 119 used for R2H_S;Cd_rad_out_double;Cs_H_out_2H
Exact match found for rate rule [R2H_S;Cd_rad_out_double;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CC=C=[C]C=[C]O(30146)'],
    products = ['C=C[C]=CC=[C]O(29679)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.992e-05,'s^-1'), n=4.805, Ea=(234.476,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_MMS;Cd_rad_out_single;Cs_H_out_2H] for rate rule [R4H_MMS;Cd_rad_out_singleDe_Cd;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CC=C=C[C]=[C]O(30147)'],
    products = ['C=C[C]=CC=[C]O(29679)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out;Cs_H_out_2H] for rate rule [R5H_SMMS;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=C[C]=CC=[C]O(29679)'],
    products = ['CC=C=C[CH][C]=O(27357)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(60975.7,'s^-1'), n=1.58648, Ea=(80.9836,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7Hall;C_rad_out_2H;XH_out] for rate rule [R7HJ_5;C_rad_out_2H;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=C[C]=CC=[C]O(29679)'],
    products = ['O[C]=CC1[C]=CC1(30148)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.53664e+10,'s^-1'), n=0.43543, Ea=(123.966,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra_cs2H] + [R4;doublebond_intra;radadd_intra_cs2H] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 121.6 to 124.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=C[C]=CC=[C]O(29679)'],
    products = ['[CH2]C=C1[CH]C=C1O(30149)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D_D;doublebond_intra;radadd_intra_cdsingle] for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleNd]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=C[C]=CC=[C]O(29679)'],
    products = ['CC=C=CC=C=O(27359)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['CH2(19)', '[CH]=C=CC=[C]O(30150)'],
    products = ['C=C[C]=CC=[C]O(29679)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4694',
    isomers = [
        'C=C[C]=CC=[C]O(29679)',
    ],
    reactants = [
        ('HCCOH(50)', 'CH2CHCCH(26391)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4694',
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

