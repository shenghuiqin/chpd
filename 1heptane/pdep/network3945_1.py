species(
    label = '[CH2]C(C=C)C([CH2])OO(22556)',
    structure = SMILES('[CH2]C(C=C)C([CH2])OO'),
    E0 = (246.025,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3615,1310,387.5,850,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,200,800],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.171941,0.0870109,-8.01871e-05,3.94078e-08,-7.82882e-12,29744.2,36.1939], Tmin=(100,'K'), Tmax=(1208.57,'K')), NASAPolynomial(coeffs=[16.1943,0.0328437,-1.29583e-05,2.32342e-09,-1.57711e-13,25788.2,-45.8643], Tmin=(1208.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(246.025,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJCOOH) + radical(Isobutyl)"""),
)

species(
    label = 'CH2CHOOH(64)',
    structure = SMILES('C=COO'),
    E0 = (-53.0705,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.754187,'amu*angstrom^2'), symmetry=1, barrier=(17.3402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.754176,'amu*angstrom^2'), symmetry=1, barrier=(17.34,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (60.052,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3284.22,'J/mol'), sigma=(4.037,'angstroms'), dipoleMoment=(1.3,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.79821,0.0274377,-2.03468e-05,7.62127e-09,-1.12671e-12,-6315.53,9.11829], Tmin=(298,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.82519,0.0274167,-2.04456e-05,7.77399e-09,-1.18661e-12,-6325.27,8.96641], Tmin=(1000,'K'), Tmax=(2000,'K'))], Tmin=(298,'K'), Tmax=(2000,'K'), E0=(-53.0705,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""CH2CHOOH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = '[CH2]C1CC1C([CH2])OO(23266)',
    structure = SMILES('[CH2]C1CC1C([CH2])OO'),
    E0 = (270.728,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.102536,0.0720555,-2.94231e-05,-1.76019e-08,1.32129e-11,32713.6,33.2793], Tmin=(100,'K'), Tmax=(975.756,'K')), NASAPolynomial(coeffs=[18.3739,0.0282701,-9.94696e-06,1.77293e-09,-1.24531e-13,27666.7,-62.0116], Tmin=(975.756,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(270.728,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(Isobutyl) + radical(CJCOOH)"""),
)

species(
    label = '[CH2]C1CC(OO)C1[CH2](23267)',
    structure = SMILES('[CH2]C1CC(OO)C1[CH2]'),
    E0 = (257.429,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.498123,0.0617054,-1.41744e-06,-4.49469e-08,2.27453e-11,31101.5,30.9672], Tmin=(100,'K'), Tmax=(955.462,'K')), NASAPolynomial(coeffs=[16.948,0.0300014,-9.98663e-06,1.73992e-09,-1.21805e-13,26261.7,-56.5212], Tmin=(955.462,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(257.429,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    label = '[CH2]C(OO)C(=C)C=C(23268)',
    structure = SMILES('[CH2]C(OO)C(=C)C=C'),
    E0 = (146.346,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3615,1310,387.5,850,1000,350,440,435,1725,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,305.656,305.673],'cm^-1')),
        HinderedRotor(inertia=(0.304712,'amu*angstrom^2'), symmetry=1, barrier=(20.2,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.304616,'amu*angstrom^2'), symmetry=1, barrier=(20.2,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.304668,'amu*angstrom^2'), symmetry=1, barrier=(20.2003,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.304682,'amu*angstrom^2'), symmetry=1, barrier=(20.2002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.304727,'amu*angstrom^2'), symmetry=1, barrier=(20.2,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.367836,0.084404,-6.74536e-05,2.03425e-08,-7.49753e-14,17768.7,33.2109], Tmin=(100,'K'), Tmax=(1036.48,'K')), NASAPolynomial(coeffs=[20.5576,0.0238862,-9.16061e-06,1.68679e-09,-1.19077e-13,12343.9,-73.7362], Tmin=(1036.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(146.346,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CJCOOH)"""),
)

species(
    label = '[CH2]C(C=C)C(=C)OO(23269)',
    structure = SMILES('[CH2]C(C=C)C(=C)OO'),
    E0 = (170.199,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3615,1310,387.5,850,1000,350,440,435,1725,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,378.643,380.049],'cm^-1')),
        HinderedRotor(inertia=(0.080848,'amu*angstrom^2'), symmetry=1, barrier=(8.23524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0789785,'amu*angstrom^2'), symmetry=1, barrier=(8.23067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0805815,'amu*angstrom^2'), symmetry=1, barrier=(8.24759,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0796693,'amu*angstrom^2'), symmetry=1, barrier=(8.23418,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.710658,'amu*angstrom^2'), symmetry=1, barrier=(71.5977,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.340468,0.0843032,-9.18716e-05,5.81564e-08,-1.53577e-11,20598.7,31.4408], Tmin=(100,'K'), Tmax=(907.99,'K')), NASAPolynomial(coeffs=[10.6289,0.0389783,-1.69932e-05,3.17789e-09,-2.19981e-13,18730.3,-17.2016], Tmin=(907.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(170.199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]OO(104)',
    structure = SMILES('[CH2][CH]OO'),
    E0 = (224.812,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.00920734,'amu*angstrom^2'), symmetry=1, barrier=(3.53679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00921023,'amu*angstrom^2'), symmetry=1, barrier=(3.53685,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43223,'amu*angstrom^2'), symmetry=1, barrier=(32.9297,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.52392,0.0279389,-1.79903e-05,3.58397e-09,4.58838e-13,27095.5,18.6054], Tmin=(100,'K'), Tmax=(1150.47,'K')), NASAPolynomial(coeffs=[9.41961,0.0107482,-4.42273e-06,8.47943e-10,-6.05179e-14,25059.9,-17.5802], Tmin=(1150.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(224.812,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CCsJOOH) + radical(CJCOOH)"""),
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
    label = '[CH2]C(C=C)OO(11483)',
    structure = SMILES('[CH2]C(C=C)OO'),
    E0 = (96.1391,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,289.03],'cm^-1')),
        HinderedRotor(inertia=(0.0306924,'amu*angstrom^2'), symmetry=1, barrier=(1.8193,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.382338,'amu*angstrom^2'), symmetry=1, barrier=(22.665,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.382287,'amu*angstrom^2'), symmetry=1, barrier=(22.665,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.382368,'amu*angstrom^2'), symmetry=1, barrier=(22.6649,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.0972,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08094,0.055951,-4.09109e-05,1.20414e-08,-4.9877e-13,11674.9,26.5848], Tmin=(100,'K'), Tmax=(1131.74,'K')), NASAPolynomial(coeffs=[14.3737,0.0201488,-8.27627e-06,1.54581e-09,-1.08351e-13,7950.13,-42.3535], Tmin=(1131.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(96.1391,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJCOOH)"""),
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
    label = 'HO2(10)',
    structure = SMILES('[O]O'),
    E0 = (2.67648,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1112.81,1388.53,3298.45],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (33.0067,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02956,-0.00263985,1.5223e-05,-1.71671e-08,6.26738e-12,322.677,4.84428], Tmin=(100,'K'), Tmax=(923.913,'K')), NASAPolynomial(coeffs=[4.15133,0.00191146,-4.11274e-07,6.34957e-11,-4.86385e-15,83.4208,3.09341], Tmin=(923.913,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.67648,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HO2""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = '[CH2]C(C=C)C=C(18999)',
    structure = SMILES('[CH2]C(C=C)C=C'),
    E0 = (260.262,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.16486,'amu*angstrom^2'), symmetry=1, barrier=(3.79045,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165408,'amu*angstrom^2'), symmetry=1, barrier=(3.80306,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14666,'amu*angstrom^2'), symmetry=1, barrier=(26.364,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22057,0.0552061,-3.78148e-05,1.41307e-08,-2.20455e-12,31407.2,23.4799], Tmin=(100,'K'), Tmax=(1478.14,'K')), NASAPolynomial(coeffs=[11.1866,0.0282371,-1.04471e-05,1.78735e-09,-1.1691e-13,28461,-28.495], Tmin=(1478.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.262,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(OO)[C](C)C=C(23270)',
    structure = SMILES('[CH2]C=C(C)C([CH2])OO'),
    E0 = (172.558,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3615,1310,387.5,850,1000,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1390,370,380,2900,435,200],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.215559,0.0867672,-7.55495e-05,3.37548e-08,-6.0746e-12,20910.5,34.0795], Tmin=(100,'K'), Tmax=(1323.77,'K')), NASAPolynomial(coeffs=[17.9072,0.0320067,-1.34995e-05,2.506e-09,-1.73191e-13,16112.3,-58.4358], Tmin=(1323.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(172.558,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(CJCOOH)"""),
)

species(
    label = '[CH2]C(C=C)[C](C)OO(23087)',
    structure = SMILES('[CH2]C(C=C)[C](C)OO'),
    E0 = (218.962,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,360,370,350,3000,3100,440,815,1455,1000,3615,1310,387.5,850,1000,200,800],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.365654,0.0827311,-7.72828e-05,4.17769e-08,-9.59874e-12,26463.6,32.5729], Tmin=(100,'K'), Tmax=(1023.94,'K')), NASAPolynomial(coeffs=[10.6545,0.0425376,-1.84017e-05,3.44023e-09,-2.38584e-13,24356.6,-17.3081], Tmin=(1023.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(218.962,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(C2CsJOOH)"""),
)

species(
    label = '[CH2][C](OO)C(C)C=C(23271)',
    structure = SMILES('[CH2][C](OO)C(C)C=C'),
    E0 = (227.842,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.383231,0.0826101,-7.40369e-05,3.69312e-08,-7.80638e-12,27530.7,32.7209], Tmin=(100,'K'), Tmax=(1103.03,'K')), NASAPolynomial(coeffs=[11.5639,0.0420643,-1.88982e-05,3.60515e-09,-2.52966e-13,25064.3,-22.3155], Tmin=(1103.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(227.842,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C2CsJOOH) + radical(CJCOOH)"""),
)

species(
    label = '[CH2]C(OO)C(C)[C]=C(23272)',
    structure = SMILES('[CH2]C(OO)C(C)[C]=C'),
    E0 = (278.785,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,200,800],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0793705,0.086335,-7.98243e-05,3.96891e-08,-8.12825e-12,33671.1,33.9866], Tmin=(100,'K'), Tmax=(1156.78,'K')), NASAPolynomial(coeffs=[14.2501,0.0373341,-1.62844e-05,3.07008e-09,-2.14218e-13,30392.6,-36.4429], Tmin=(1156.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(278.785,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJCOOH) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C](C=C)C(C)OO(23273)',
    structure = SMILES('[CH2]C=C([CH2])C(C)OO'),
    E0 = (110.094,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.377021,0.0841759,-6.28257e-05,1.92942e-08,-1.02841e-12,13409.1,32.9478], Tmin=(100,'K'), Tmax=(1119.65,'K')), NASAPolynomial(coeffs=[19.4382,0.0298599,-1.21298e-05,2.25068e-09,-1.57272e-13,7939.21,-69.5001], Tmin=(1119.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(110.094,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(C=C)C(C)O[O](23274)',
    structure = SMILES('[CH2]C(C=C)C(C)O[O]'),
    E0 = (184.067,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,492.5,1135,1000,408.654,2250.54],'cm^-1')),
        HinderedRotor(inertia=(0.0937967,'amu*angstrom^2'), symmetry=1, barrier=(11.1142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.483396,'amu*angstrom^2'), symmetry=1, barrier=(11.1142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0937917,'amu*angstrom^2'), symmetry=1, barrier=(11.1143,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0937948,'amu*angstrom^2'), symmetry=1, barrier=(11.1142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.2538,'amu*angstrom^2'), symmetry=1, barrier=(74.8113,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.215437,0.0800314,-6.94931e-05,3.34033e-08,-6.62535e-12,22277.1,33.5665], Tmin=(100,'K'), Tmax=(1196.69,'K')), NASAPolynomial(coeffs=[13.4573,0.0357697,-1.40129e-05,2.49575e-09,-1.68448e-13,19107.8,-32.6956], Tmin=(1196.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(184.067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(ROOJ)"""),
)

species(
    label = '[CH2]C([C]=C)C(C)OO(23275)',
    structure = SMILES('[CH2]C([C]=C)C(C)OO'),
    E0 = (269.905,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,200,800],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0910111,0.08613,-8.20186e-05,4.32917e-08,-9.43907e-12,32602.7,33.7326], Tmin=(100,'K'), Tmax=(1092.52,'K')), NASAPolynomial(coeffs=[13.2517,0.037945,-1.58613e-05,2.92148e-09,-2.01125e-13,29727.1,-30.9243], Tmin=(1092.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(269.905,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=CC(C)C([CH2])OO(23276)',
    structure = SMILES('[CH]=CC(C)C([CH2])OO'),
    E0 = (288.039,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,200],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.143278,0.0874206,-7.92923e-05,3.76261e-08,-7.23068e-12,34795.3,34.713], Tmin=(100,'K'), Tmax=(1241.52,'K')), NASAPolynomial(coeffs=[16.6056,0.0334584,-1.40958e-05,2.61733e-09,-1.81156e-13,30636.5,-49.7143], Tmin=(1241.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(288.039,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJCOOH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC([CH2])C(C)OO(23277)',
    structure = SMILES('[CH]=CC([CH2])C(C)OO'),
    E0 = (279.159,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,200],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.102855,0.0868992,-8.04893e-05,4.00867e-08,-8.11617e-12,33725.6,34.3542], Tmin=(100,'K'), Tmax=(1182.63,'K')), NASAPolynomial(coeffs=[15.4384,0.0343331,-1.38154e-05,2.50082e-09,-1.70623e-13,30049.7,-43.23], Tmin=(1182.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C(O[O])C(C)C=C(23278)',
    structure = SMILES('[CH2]C(O[O])C(C)C=C'),
    E0 = (192.948,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,492.5,1135,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0972373,'amu*angstrom^2'), symmetry=1, barrier=(2.23568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0972407,'amu*angstrom^2'), symmetry=1, barrier=(2.23575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.610581,'amu*angstrom^2'), symmetry=1, barrier=(14.0385,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0326652,'amu*angstrom^2'), symmetry=1, barrier=(14.0385,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.610579,'amu*angstrom^2'), symmetry=1, barrier=(14.0384,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.171565,0.0805927,-6.84284e-05,3.10992e-08,-5.79875e-12,23347,33.9376], Tmin=(100,'K'), Tmax=(1267.57,'K')), NASAPolynomial(coeffs=[14.6819,0.0348032,-1.42426e-05,2.6007e-09,-1.78046e-13,19668.4,-39.5068], Tmin=(1267.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.948,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJCOOH) + radical(ROOJ)"""),
)

species(
    label = '[CH2]C(OO)C1[CH]CC1(23279)',
    structure = SMILES('[CH2]C(OO)C1[CH]CC1'),
    E0 = (258.107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.450794,0.065231,-2.12135e-05,-1.50061e-08,9.00092e-12,31181.7,33.5882], Tmin=(100,'K'), Tmax=(1077.23,'K')), NASAPolynomial(coeffs=[15.3972,0.0351018,-1.45867e-05,2.75565e-09,-1.95043e-13,26489.6,-46.4637], Tmin=(1077.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(258.107,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CJCOOH) + radical(cyclobutane)"""),
)

species(
    label = '[CH2]C1[CH]CCC1OO(23280)',
    structure = SMILES('[CH2]C1[CH]CCC1OO'),
    E0 = (174.163,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.970134,0.050528,2.06281e-05,-5.99747e-08,2.58112e-11,21070.4,31.2132], Tmin=(100,'K'), Tmax=(982.668,'K')), NASAPolynomial(coeffs=[14.561,0.0336819,-1.23896e-05,2.2709e-09,-1.61779e-13,16541.7,-43.5698], Tmin=(982.668,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(174.163,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = 'C=CC(=C)C(C)OO(23281)',
    structure = SMILES('C=CC(=C)C(C)OO'),
    E0 = (-67.6166,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.292148,0.0815561,-5.3812e-05,7.08825e-09,4.11796e-12,-7966.7,30.9445], Tmin=(100,'K'), Tmax=(1032.94,'K')), NASAPolynomial(coeffs=[19.8663,0.0277034,-1.07658e-05,1.99635e-09,-1.41391e-13,-13422.7,-73.2136], Tmin=(1032.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-67.6166,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CC(C)C(=C)OO(23282)',
    structure = SMILES('C=CC(C)C(=C)OO'),
    E0 = (-34.8837,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.512497,0.0804338,-7.19639e-05,3.6293e-08,-7.81773e-12,-4073.08,29.038], Tmin=(100,'K'), Tmax=(1078.89,'K')), NASAPolynomial(coeffs=[10.7567,0.0424545,-1.91619e-05,3.66647e-09,-2.57735e-13,-6283.6,-21.1623], Tmin=(1078.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-34.8837,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'OH(5)',
    structure = SMILES('[OH]'),
    E0 = (28.372,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3287.46],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (17.0073,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.4858,0.00133397,-4.70043e-06,5.64379e-09,-2.06318e-12,3411.96,1.99788], Tmin=(100,'K'), Tmax=(1005.25,'K')), NASAPolynomial(coeffs=[2.88225,0.00103869,-2.35652e-07,1.40229e-11,6.34581e-16,3669.56,5.59053], Tmin=(1005.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.372,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = '[CH2]C(C=C)C1CO1(23283)',
    structure = SMILES('[CH2]C(C=C)C1CO1'),
    E0 = (150.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16763,0.0509252,7.06116e-06,-5.5387e-08,3.00123e-11,18168.4,25.6777], Tmin=(100,'K'), Tmax=(857.013,'K')), NASAPolynomial(coeffs=[14.7137,0.0224032,-3.7563e-06,2.76155e-10,-8.06232e-15,14572.2,-45.0191], Tmin=(857.013,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(150.121,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1OCC1C=C(23284)',
    structure = SMILES('[CH2]C1OCC1C=C'),
    E0 = (150.017,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.184,0.0492929,1.1389e-05,-5.752e-08,2.94191e-11,18156.3,23.2606], Tmin=(100,'K'), Tmax=(882.068,'K')), NASAPolynomial(coeffs=[14.8105,0.0231593,-4.81068e-06,5.56184e-10,-3.11827e-14,14365.2,-48.6329], Tmin=(882.068,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(150.017,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Oxetane) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]C(C[CH]C=C)OO(22555)',
    structure = SMILES('[CH2]C=CCC([CH2])OO'),
    E0 = (188.169,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.1103,0.0852949,-7.32149e-05,3.25377e-08,-5.85894e-12,22783.5,34.2796], Tmin=(100,'K'), Tmax=(1317.65,'K')), NASAPolynomial(coeffs=[17.0093,0.0333251,-1.40533e-05,2.60502e-09,-1.79796e-13,18272,-53.0352], Tmin=(1317.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.169,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(CJCOOH)"""),
)

species(
    label = '[CH2]C([CH]CC=C)OO(23285)',
    structure = SMILES('[CH2]C([CH]CC=C)OO'),
    E0 = (249.33,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.195152,0.0857051,-7.94024e-05,3.99938e-08,-8.38724e-12,30122.8,35.6505], Tmin=(100,'K'), Tmax=(1123.02,'K')), NASAPolynomial(coeffs=[13.1576,0.0395353,-1.77346e-05,3.38573e-09,-2.3781e-13,27211.3,-28.39], Tmin=(1123.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(249.33,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJCOOH) + radical(CCJCOOH)"""),
)

species(
    label = '[CH2]C(C=C)C[CH]OO(22554)',
    structure = SMILES('[CH2]C(C=C)C[CH]OO'),
    E0 = (233.288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.190423,0.0856852,-8.21048e-05,4.46282e-08,-1.01545e-11,28193.7,33.4002], Tmin=(100,'K'), Tmax=(1042.05,'K')), NASAPolynomial(coeffs=[11.8991,0.0407397,-1.74058e-05,3.23533e-09,-2.23712e-13,25753.6,-23.5693], Tmin=(1042.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(233.288,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(CCsJOOH)"""),
)

species(
    label = 'C=CC1CCC1OO(22559)',
    structure = SMILES('C=CC1CCC1OO'),
    E0 = (-14.8425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.709704,0.0533388,2.27029e-05,-6.68935e-08,2.88078e-11,-1649.81,28.1274], Tmin=(100,'K'), Tmax=(992.985,'K')), NASAPolynomial(coeffs=[17.4885,0.0310169,-1.19583e-05,2.28626e-09,-1.67668e-13,-7213.75,-63.9398], Tmin=(992.985,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-14.8425,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutane)"""),
)

species(
    label = '[CH2]C(C=C)C([O])CO(23286)',
    structure = SMILES('[CH2]C(C=C)C([O])CO'),
    E0 = (31.0911,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00758983,0.0831263,-7.556e-05,3.73818e-08,-7.51171e-12,3887.06,35.3028], Tmin=(100,'K'), Tmax=(1194.29,'K')), NASAPolynomial(coeffs=[14.9891,0.0329491,-1.25386e-05,2.20259e-09,-1.47657e-13,308.602,-39.6345], Tmin=(1194.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.0911,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C([O])C(C=C)CO(23287)',
    structure = SMILES('[CH2]C([O])C(C=C)CO'),
    E0 = (37.5978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0574624,0.0853816,-7.83996e-05,3.82603e-08,-7.56279e-12,4671.27,34.389], Tmin=(100,'K'), Tmax=(1211.39,'K')), NASAPolynomial(coeffs=[15.8715,0.0327838,-1.327e-05,2.41715e-09,-1.65615e-13,812.052,-45.5137], Tmin=(1211.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(37.5978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CJCO)"""),
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
    label = '[CH2]C([CH]C=C)OO(23288)',
    structure = SMILES('[CH2]C([CH]C=C)OO'),
    E0 = (189.612,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,299.942,689.601],'cm^-1')),
        HinderedRotor(inertia=(0.433255,'amu*angstrom^2'), symmetry=1, barrier=(27.6611,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.375115,'amu*angstrom^2'), symmetry=1, barrier=(23.9445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.43332,'amu*angstrom^2'), symmetry=1, barrier=(27.6603,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.433203,'amu*angstrom^2'), symmetry=1, barrier=(27.6612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.42926,'amu*angstrom^2'), symmetry=1, barrier=(91.246,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.303791,0.0763382,-6.95816e-05,3.22402e-08,-5.98028e-12,22942.2,27.9063], Tmin=(100,'K'), Tmax=(1291.81,'K')), NASAPolynomial(coeffs=[16.5528,0.0260241,-1.11584e-05,2.08945e-09,-1.45277e-13,18744.1,-54.646], Tmin=(1291.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(189.612,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CCJCO) + radical(CJCOOH)"""),
)

species(
    label = '[CH2]C([CH]OO)C=C(21247)',
    structure = SMILES('[CH2]C([CH]OO)C=C'),
    E0 = (257.068,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,240.222,1289.74],'cm^-1')),
        HinderedRotor(inertia=(0.381608,'amu*angstrom^2'), symmetry=1, barrier=(15.6606,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.88075,'amu*angstrom^2'), symmetry=1, barrier=(77.3378,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.379349,'amu*angstrom^2'), symmetry=1, barrier=(15.6225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0157602,'amu*angstrom^2'), symmetry=1, barrier=(77.3484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.67602,'amu*angstrom^2'), symmetry=1, barrier=(77.359,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.80772,0.0712611,-6.9565e-05,3.82961e-08,-8.75836e-12,31032.3,28.9378], Tmin=(100,'K'), Tmax=(1041.91,'K')), NASAPolynomial(coeffs=[10.9254,0.0324185,-1.3645e-05,2.51581e-09,-1.7315e-13,28923.9,-20.2897], Tmin=(1041.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(257.068,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCsJOOH) + radical(Isobutyl)"""),
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
    E0 = (246.025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (283.263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (305.02,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (358.138,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (393.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (343.213,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (396.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (246.025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (301.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (375.9,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (378.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (363.964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (430.664,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (363.964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (324.273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (314.213,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (332.348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (312.199,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (280.895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (499.525,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (370.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (360.667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (309.425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (309.425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (291.669,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (321.881,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (405.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (494.931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (404.652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (254.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (358.575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (362.34,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (571.175,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (638.631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    products = ['CH2CHOOH(64)', 'butadiene13(2459)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    products = ['[CH2]C1CC1C([CH2])OO(23266)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.32e+08,'s^-1'), n=0.97, Ea=(37.2376,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 335 used for R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H
Exact match found for rate rule [R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    products = ['[CH2]C1CC(OO)C1[CH2](23267)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.07e+06,'s^-1'), n=1.46, Ea=(58.9944,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 343 used for R5_SS_D;doublebond_intra_2H_pri;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[CH2]C(OO)C(=C)C=C(23268)'],
    products = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.88727,'m^3/(mol*s)'), n=1.9687, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 18 used for Cds-CdCs_Cds-HH;HJ
Exact match found for rate rule [Cds-CdCs_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -8.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', '[CH2]C(C=C)C(=C)OO(23269)'],
    products = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(170.641,'m^3/(mol*s)'), n=1.56204, Ea=(11.2897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;HJ] for rate rule [Cds-OsCs_Cds;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2][CH]OO(104)', 'butadiene13(2459)'],
    products = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.0136943,'m^3/(mol*s)'), n=2.49, Ea=(21.946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CdH_Cds-HH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C2H3(30)', '[CH2]C(C=C)OO(11483)'],
    products = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6870,'cm^3/(mol*s)'), n=2.41, Ea=(13.7235,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;CdsJ-H] for rate rule [Cds-Cs\O2s/H_Cds-HH;CdsJ-H]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH2CHOOH(64)', '[CH2][CH]C=C(2458)'],
    products = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.616939,'m^3/(mol*s)'), n=2.06429, Ea=(24.382,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cds_Cds;CJ] + [Cds-OsH_Cds;YJ] for rate rule [Cds-OsH_Cds;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 20.0 to 24.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['HO2(10)', '[CH2]C(C=C)C=C(18999)'],
    products = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(21.2,'cm^3/(mol*s)'), n=3.29, Ea=(38.0744,'kJ/mol'), T0=(1,'K'), Tmin=(400,'K'), Tmax=(1100,'K'), comment="""From training reaction 2772 used for Cds-CsH_Cds-HH;OJ-O2s
Exact match found for rate rule [Cds-CsH_Cds-HH;OJ-O2s]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(OO)[C](C)C=C(23270)'],
    products = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.307e+09,'s^-1'), n=1.31, Ea=(203.342,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 163 used for R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    products = ['[CH2]C(C=C)[C](C)OO(23087)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(309968,'s^-1'), n=2.08546, Ea=(132.285,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_2H;Cs_H_out] for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_OOH/Cs]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    products = ['[CH2][C](OO)C(C)C=C(23271)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(OO)C(C)[C]=C(23272)'],
    products = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 204 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    products = ['[CH2][C](C=C)C(C)OO(23273)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(C=C)C(C)O[O](23274)'],
    products = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.18e+08,'s^-1'), n=1.06, Ea=(140.206,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 247 used for R4H_SSS_O(Cs)Cs;O_rad_out;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS_O(Cs)Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C([C]=C)C(C)OO(23275)'],
    products = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=CC(C)C([CH2])OO(23276)'],
    products = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=CC([CH2])C(C)OO(23277)'],
    products = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5H;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(O[O])C(C)C=C(23278)'],
    products = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(5.07e+06,'s^-1'), n=1.55, Ea=(87.9477,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 255 used for R5H_SSSS_OCC;O_rad_out;Cs_H_out_2H
Exact match found for rate rule [R5H_SSSS_OCC;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][CH]OO(104)', '[CH2][CH]C=C(2458)'],
    products = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.49215e+07,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    products = ['[CH2]C(OO)C1[CH]CC1(23279)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.01e+08,'s^-1'), n=1.02, Ea=(124.265,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 18 used for R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    products = ['[CH2]C1[CH]CCC1OO(23280)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.12e+09,'s^-1'), n=0.63, Ea=(114.642,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 3 used for R5_SS_D;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    products = ['C=CC(=C)C(C)OO(23281)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    products = ['C=CC(C)C(=C)OO(23282)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    products = ['OH(5)', '[CH2]C(C=C)C1CO1(23283)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.98e+12,'s^-1','*|/',1.2), n=0, Ea=(45.6433,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 1 used for R2OO_S;C_pri_rad_intra;OOH
Exact match found for rate rule [R2OO_S;C_pri_rad_intra;OOH]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    products = ['OH(5)', '[CH2]C1OCC1C=C(23284)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.47e+11,'s^-1','*|/',1.74), n=0, Ea=(75.8559,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 7 used for R3OO_SS;C_pri_rad_intra;OOH
Exact match found for rate rule [R3OO_SS;C_pri_rad_intra;OOH]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    products = ['[CH2]C(C[CH]C=C)OO(22555)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C([CH]CC=C)OO(23285)'],
    products = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HH)CJ;CsJ;C] for rate rule [cCs(-HH)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    products = ['[CH2]C(C=C)C[CH]OO(22554)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    products = ['C=CC1CCC1OO(22559)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    products = ['[CH2]C(C=C)C([O])CO(23286)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(3.39e+11,'s^-1'), n=0, Ea=(112.55,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 1 used for R2OOH_S;C_rad_out_2H
Exact match found for rate rule [R2OOH_S;C_rad_out_2H]
Euclidian distance = 0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    products = ['[CH2]C([O])C(C=C)CO(23287)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(4.47e+10,'s^-1'), n=0, Ea=(116.315,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 3 used for R3OOH_SS;C_rad_out_2H
Exact match found for rate rule [R3OOH_SS;C_rad_out_2H]
Euclidian distance = 0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['CH2(19)', '[CH2]C([CH]C=C)OO(23288)'],
    products = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeC;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction34',
    reactants = ['CH2(19)', '[CH2]C([CH]OO)C=C(21247)'],
    products = ['[CH2]C(C=C)C([CH2])OO(22556)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/CsO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '3945',
    isomers = [
        '[CH2]C(C=C)C([CH2])OO(22556)',
    ],
    reactants = [
        ('CH2CHOOH(64)', 'butadiene13(2459)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '3945',
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

