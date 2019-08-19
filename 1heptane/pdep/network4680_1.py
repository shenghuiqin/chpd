species(
    label = 'C=C[C]=CO[CH]C=O(29462)',
    structure = SMILES('C=C[C]=COC=C[O]'),
    E0 = (153.251,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,1685,370,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.39521,'amu*angstrom^2'), symmetry=1, barrier=(32.0786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.39529,'amu*angstrom^2'), symmetry=1, barrier=(32.0804,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.39492,'amu*angstrom^2'), symmetry=1, barrier=(32.0719,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0187316,0.0693981,-2.39949e-05,-4.1212e-08,2.77723e-11,18592.2,28.833], Tmin=(100,'K'), Tmax=(906.39,'K')), NASAPolynomial(coeffs=[25.943,0.00424209,2.32726e-06,-6.23601e-10,4.21083e-14,11869.7,-104.849], Tmin=(906.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(153.251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CJC=C)"""),
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
    label = 'C=CC1=COC1[CH][O](29470)',
    structure = SMILES('C=CC1=COC1[CH][O]'),
    E0 = (258.886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.23737,0.0629147,-9.46792e-06,-5.05901e-08,2.87409e-11,31290.5,25.8258], Tmin=(100,'K'), Tmax=(936.932,'K')), NASAPolynomial(coeffs=[25.0816,0.00701031,-2.74497e-07,1.0801e-11,-8.05957e-15,24433.3,-104.165], Tmin=(936.932,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(258.886,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C1[C]=COC=CO1(29471)',
    structure = SMILES('[CH2]C1[C]=COC=CO1'),
    E0 = (246.586,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.484663,0.00875862,0.000238056,-3.84695e-07,1.69475e-10,29848.8,26.5729], Tmin=(100,'K'), Tmax=(894.921,'K')), NASAPolynomial(coeffs=[55.2215,-0.0536209,3.70951e-05,-7.39643e-09,4.96123e-13,12752.7,-272.204], Tmin=(894.921,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(246.586,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(Cds_S) + radical(CJC(C)OC)"""),
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
    label = 'C=C=C=COC=C[O](29472)',
    structure = SMILES('C=C=C=COC=C[O]'),
    E0 = (183.641,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,563.333,586.667,610,1970,2140,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.40993,'amu*angstrom^2'), symmetry=1, barrier=(32.4171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4092,'amu*angstrom^2'), symmetry=1, barrier=(32.4003,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.154362,0.0691187,-3.75343e-05,-1.9276e-08,1.75987e-11,22239.7,27.8316], Tmin=(100,'K'), Tmax=(928.296,'K')), NASAPolynomial(coeffs=[24.0011,0.0061289,-6.33502e-09,-8.17335e-11,2.04265e-15,16099,-94.6693], Tmin=(928.296,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(183.641,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=COJ)"""),
)

species(
    label = 'C=CC#COC=C[O](29473)',
    structure = SMILES('C=CC#COC=C[O]'),
    E0 = (168.263,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2100,2250,500,550,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.42521,'amu*angstrom^2'), symmetry=1, barrier=(32.7683,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.42701,'amu*angstrom^2'), symmetry=1, barrier=(32.8098,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43528,'amu*angstrom^2'), symmetry=1, barrier=(33,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.515957,0.0576134,-5.72427e-06,-5.0099e-08,2.76232e-11,20380.2,27.7983], Tmin=(100,'K'), Tmax=(941.206,'K')), NASAPolynomial(coeffs=[23.7402,0.00663159,-5.22817e-07,8.25848e-11,-1.35552e-14,13894.9,-94.066], Tmin=(941.206,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(168.263,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsCtH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtOs) + radical(C=COJ)"""),
)

species(
    label = 'C=C[C]=COC=C=O(29474)',
    structure = SMILES('C=C[C]=COC=C=O'),
    E0 = (193.264,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2120,512.5,787.5,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.2808,'amu*angstrom^2'), symmetry=1, barrier=(29.448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.28501,'amu*angstrom^2'), symmetry=1, barrier=(29.5449,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.28618,'amu*angstrom^2'), symmetry=1, barrier=(29.5717,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.224701,0.0794111,-8.63791e-05,4.57775e-08,-9.15609e-12,23408,29.0975], Tmin=(100,'K'), Tmax=(1365.32,'K')), NASAPolynomial(coeffs=[20.5632,0.0104354,-1.7299e-06,1.13663e-10,-1.70819e-15,18484,-74.9094], Tmin=(1365.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(193.264,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-(Cdd-O2d)OsH) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
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
    label = 'C=[C]C=COC=C[O](29475)',
    structure = SMILES('C=C=C[CH]OC=C[O]'),
    E0 = (108.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.194799,0.0660784,-2.13154e-05,-3.57636e-08,2.28568e-11,13171.4,28.5084], Tmin=(100,'K'), Tmax=(939.704,'K')), NASAPolynomial(coeffs=[23.8309,0.00960617,-1.62763e-06,2.53586e-10,-2.32831e-14,6780.42,-94.4213], Tmin=(939.704,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(108.24,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=CCJ(O)C)"""),
)

species(
    label = 'C=CC=[C]OC=C[O](29476)',
    structure = SMILES('C=CC=[C]OC=C[O]'),
    E0 = (194,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,1685,370,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.08686,'amu*angstrom^2'), symmetry=1, barrier=(24.9891,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08999,'amu*angstrom^2'), symmetry=1, barrier=(25.061,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09014,'amu*angstrom^2'), symmetry=1, barrier=(25.0645,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.303398,0.0646758,-2.21501e-05,-3.45653e-08,2.28582e-11,23481.3,31.3416], Tmin=(100,'K'), Tmax=(925.94,'K')), NASAPolynomial(coeffs=[23.2781,0.00838809,-5.6173e-07,-7.0766e-14,-3.38936e-15,17385,-87.675], Tmin=(925.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(194,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(C=COJ)"""),
)

species(
    label = 'C=C[C]=COC=[C]O(29477)',
    structure = SMILES('C=C[C]=COC=[C]O'),
    E0 = (251.533,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1670,1700,300,440,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.1804,'amu*angstrom^2'), symmetry=1, barrier=(27.1397,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18085,'amu*angstrom^2'), symmetry=1, barrier=(27.15,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18145,'amu*angstrom^2'), symmetry=1, barrier=(27.1639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1813,'amu*angstrom^2'), symmetry=1, barrier=(27.1603,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.31505,0.0914452,-0.000100857,5.20287e-08,-9.81706e-12,30465.7,35.3404], Tmin=(100,'K'), Tmax=(1539.64,'K')), NASAPolynomial(coeffs=[24.8046,0.00451359,2.41849e-06,-7.35507e-10,5.70228e-14,24683.3,-94.6022], Tmin=(1539.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(251.533,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CJO)"""),
)

species(
    label = '[CH]=CC=COC=C[O](29478)',
    structure = SMILES('[CH]=CC=COC=C[O]'),
    E0 = (201.352,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.22377,'amu*angstrom^2'), symmetry=1, barrier=(28.1369,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.22349,'amu*angstrom^2'), symmetry=1, barrier=(28.1304,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.22351,'amu*angstrom^2'), symmetry=1, barrier=(28.1308,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0192482,0.0664355,-1.184e-05,-5.62057e-08,3.33378e-11,24380.1,29.6367], Tmin=(100,'K'), Tmax=(915.818,'K')), NASAPolynomial(coeffs=[27.6817,0.00139074,3.3412e-06,-7.55703e-10,4.76263e-14,16974.3,-114.156], Tmin=(915.818,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(201.352,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=COJ)"""),
)

species(
    label = 'C=C[C]=CO[C]=CO(29479)',
    structure = SMILES('C=C[C]=CO[C]=CO'),
    E0 = (251.533,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1670,1700,300,440,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.18148,'amu*angstrom^2'), symmetry=1, barrier=(27.1645,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18128,'amu*angstrom^2'), symmetry=1, barrier=(27.1601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18147,'amu*angstrom^2'), symmetry=1, barrier=(27.1644,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18102,'amu*angstrom^2'), symmetry=1, barrier=(27.1541,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.31502,0.0914449,-0.000100856,5.2028e-08,-9.81684e-12,30465.7,35.3403], Tmin=(100,'K'), Tmax=(1539.67,'K')), NASAPolynomial(coeffs=[24.8037,0.00451481,2.41787e-06,-7.35374e-10,5.70126e-14,24683.7,-94.5973], Tmin=(1539.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(251.533,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(C=CJC=C)"""),
)

species(
    label = 'C=CC=CO[C]=C[O](29480)',
    structure = SMILES('C=CC=CO[C]=C[O]'),
    E0 = (194,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,1685,370,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.08686,'amu*angstrom^2'), symmetry=1, barrier=(24.9891,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08999,'amu*angstrom^2'), symmetry=1, barrier=(25.061,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09014,'amu*angstrom^2'), symmetry=1, barrier=(25.0645,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.303398,0.0646758,-2.21501e-05,-3.45653e-08,2.28582e-11,23481.3,31.3416], Tmin=(100,'K'), Tmax=(925.94,'K')), NASAPolynomial(coeffs=[23.2781,0.00838809,-5.6173e-07,-7.0766e-14,-3.38936e-15,17385,-87.675], Tmin=(925.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(194,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = 'C=CC=COC=[C][O](29481)',
    structure = SMILES('[CH2]C=C[CH]OC=C=O'),
    E0 = (123.148,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.383891,0.0683446,-4.62315e-05,3.18078e-09,5.73662e-12,14951.4,28.4506], Tmin=(100,'K'), Tmax=(984.569,'K')), NASAPolynomial(coeffs=[18.7089,0.0180632,-6.4467e-06,1.17283e-09,-8.39397e-14,10171.6,-65.6208], Tmin=(984.569,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(123.148,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)OsH) + radical(Allyl_P) + radical(C=CCJ(O)C)"""),
)

species(
    label = 'C=C[C]=[C]OC=CO(29482)',
    structure = SMILES('C=C[C]=[C]OC=CO'),
    E0 = (251.533,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1670,1700,300,440,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.1804,'amu*angstrom^2'), symmetry=1, barrier=(27.1397,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18085,'amu*angstrom^2'), symmetry=1, barrier=(27.15,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18145,'amu*angstrom^2'), symmetry=1, barrier=(27.1639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1813,'amu*angstrom^2'), symmetry=1, barrier=(27.1603,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.31505,0.0914452,-0.000100857,5.20287e-08,-9.81706e-12,30465.7,35.3404], Tmin=(100,'K'), Tmax=(1539.64,'K')), NASAPolynomial(coeffs=[24.8046,0.00451359,2.41849e-06,-7.35507e-10,5.70228e-14,24683.3,-94.6022], Tmin=(1539.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(251.533,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CJO)"""),
)

species(
    label = 'C=[C][C]=COC=CO(29483)',
    structure = SMILES('[CH2]C#C[CH]OC=CO'),
    E0 = (189.363,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2100,2250,500,550,3615,1277.5,1000,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.43906,0.101776,-0.000115167,5.84802e-08,-1.06166e-11,23041.3,36.346], Tmin=(100,'K'), Tmax=(1652.85,'K')), NASAPolynomial(coeffs=[28.3519,-0.0025386,6.54347e-06,-1.51856e-09,1.08254e-13,16933.1,-115.362], Tmin=(1652.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(189.363,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CtOsHH) + group(Cs-CtHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Propargyl) + radical(CCsJOC(O))"""),
)

species(
    label = '[CH]=C[C]=COC=CO(29484)',
    structure = SMILES('[CH]C=C=COC=CO'),
    E0 = (236.259,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3615,1277.5,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.295213,0.0763391,-3.13053e-05,-3.49349e-08,2.52076e-11,28586.9,30.0911], Tmin=(100,'K'), Tmax=(915.044,'K')), NASAPolynomial(coeffs=[26.0346,0.00932375,-2.68989e-07,-1.21717e-10,7.13133e-15,21755.3,-105.597], Tmin=(915.044,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(236.259,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
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
    label = 'C=C[C]=C[O](26659)',
    structure = SMILES('C=C[C]=C[O]'),
    E0 = (225.587,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.61748,'amu*angstrom^2'), symmetry=1, barrier=(37.1891,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98077,0.0322088,7.53766e-06,-4.54875e-08,2.38742e-11,27216,16.1107], Tmin=(100,'K'), Tmax=(899.949,'K')), NASAPolynomial(coeffs=[16.1069,0.00249504,1.93925e-06,-5.05324e-10,3.47549e-14,23334.1,-57.9915], Tmin=(899.949,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(225.587,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=COJ)"""),
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
    label = '[O]C=COC=C1[CH]C1(29485)',
    structure = SMILES('[O]C=COC=C1[CH]C1'),
    E0 = (188.391,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.653943,0.0482571,3.64964e-05,-1.01447e-07,4.78759e-11,22802.2,26.0292], Tmin=(100,'K'), Tmax=(925.803,'K')), NASAPolynomial(coeffs=[25.9983,0.00410073,2.1644e-06,-4.84809e-10,2.48907e-14,15309,-109.413], Tmin=(925.803,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.391,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Methylene_cyclopropane) + radical(C=COJ) + radical(Allyl_S)"""),
)

species(
    label = 'C=C[C]=COC1[CH]O1(29486)',
    structure = SMILES('C=C[C]=COC1[CH]O1'),
    E0 = (238.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.19084,0.095191,-0.000102764,5.1992e-08,-9.39778e-12,28921.1,33.5598], Tmin=(100,'K'), Tmax=(1689.73,'K')), NASAPolynomial(coeffs=[22.8857,0.00483487,4.95985e-06,-1.36476e-09,1.02556e-13,24871.3,-87.4816], Tmin=(1689.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(238.315,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCsJO) + radical(C=CJC=C)"""),
)

species(
    label = 'C=CC1=CO[CH]C1[O](29487)',
    structure = SMILES('C=CC1=CO[CH]C1[O]'),
    E0 = (173.165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.48717,0.0504683,3.7644e-05,-1.134e-07,5.61295e-11,20978.6,22.2369], Tmin=(100,'K'), Tmax=(898.406,'K')), NASAPolynomial(coeffs=[28.9991,-0.0041251,7.99537e-06,-1.75902e-09,1.19135e-13,12935.7,-128.512], Tmin=(898.406,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(173.165,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(2,3-Dihydrofuran) + radical(CC(C)OJ) + radical(CCsJOC(O))"""),
)

species(
    label = '[C]1[CH]COC=COC=1(29488)',
    structure = SMILES('[C]1[CH]OC=COCC=1'),
    E0 = (159.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70211,0.0156191,0.000129919,-1.94505e-07,7.95299e-11,19262.3,19.9323], Tmin=(100,'K'), Tmax=(933.115,'K')), NASAPolynomial(coeffs=[24.8253,0.00588559,1.87142e-06,-3.57923e-10,8.62455e-15,11055.4,-110.876], Tmin=(933.115,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(159.2,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cyclooctane) + radical(Cds_S) + radical(C=CCJ(O)C)"""),
)

species(
    label = 'C=CC=COC=C=O(29489)',
    structure = SMILES('C=CC=COC=C=O'),
    E0 = (-5.73113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.344402,0.0670731,-3.69964e-05,-1.34702e-08,1.38994e-11,-545.432,27.4592], Tmin=(100,'K'), Tmax=(935.94,'K')), NASAPolynomial(coeffs=[21.0166,0.0122338,-2.81123e-06,4.32984e-10,-3.21206e-14,-5882.7,-78.7443], Tmin=(935.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-5.73113,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-(Cdd-O2d)OsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C=C=COC=CO(29490)',
    structure = SMILES('C=C=C=COC=CO'),
    E0 = (42.1786,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.75063,0.0933101,-9.9579e-05,4.87078e-08,-8.65922e-12,5308.58,33.2855], Tmin=(100,'K'), Tmax=(1650.93,'K')), NASAPolynomial(coeffs=[26.6719,0.00244469,2.96915e-06,-7.74562e-10,5.62426e-14,-1077.86,-109.005], Tmin=(1650.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(42.1786,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=CC1=COC=CO1(29491)',
    structure = SMILES('C=CC1=COC=CO1'),
    E0 = (-49.8393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18692,0.0481653,-7.14611e-07,-3.9062e-08,2.00047e-11,-5880.57,22.3928], Tmin=(100,'K'), Tmax=(955.364,'K')), NASAPolynomial(coeffs=[16.4416,0.0177697,-5.54795e-06,9.8578e-10,-7.21293e-14,-10323,-58.5008], Tmin=(955.364,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-49.8393,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(14dioxin)"""),
)

species(
    label = 'C=C[C]=CC([O])C=O(29458)',
    structure = SMILES('C=C[C]=CC([O])C=O'),
    E0 = (224.163,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.63554,0.0678828,-6.26057e-05,3.01866e-08,-5.80801e-12,27086.8,31.3455], Tmin=(100,'K'), Tmax=(1256.7,'K')), NASAPolynomial(coeffs=[14.8715,0.0225707,-8.52123e-06,1.49538e-09,-1.00369e-13,23508.8,-40.5875], Tmin=(1256.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(224.163,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(C=CJC=C)"""),
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
    label = '[CH]=COC=[C]C=C(27408)',
    structure = SMILES('[CH]=COC=[C]C=C'),
    E0 = (467.677,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.32177,'amu*angstrom^2'), symmetry=1, barrier=(30.39,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.33075,'amu*angstrom^2'), symmetry=1, barrier=(30.5965,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32998,'amu*angstrom^2'), symmetry=1, barrier=(30.5789,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.314286,0.0744218,-7.50011e-05,3.65542e-08,-6.63494e-12,56421.5,29.1713], Tmin=(100,'K'), Tmax=(1570.51,'K')), NASAPolynomial(coeffs=[20.3202,0.00925899,-7.22108e-07,-8.81013e-11,1.16023e-14,51495,-74.7424], Tmin=(1570.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(467.677,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C=[C]C1OC=CO1(29492)',
    structure = SMILES('[CH2]C=[C]C1OC=CO1'),
    E0 = (126.587,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.752715,0.034064,9.91733e-05,-1.82008e-07,8.03561e-11,15376.6,24.0037], Tmin=(100,'K'), Tmax=(919.204,'K')), NASAPolynomial(coeffs=[32.6352,-0.00766434,8.9601e-06,-1.76491e-09,1.07828e-13,5416.9,-149.418], Tmin=(919.204,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(126.587,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cyclopentane) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[O][CH]C1CC=C=CO1(29493)',
    structure = SMILES('[O][CH]C1CC=C=CO1'),
    E0 = (233.385,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.59086,0.129459,-0.000192323,1.44991e-07,-4.28841e-11,28265.1,14.4783], Tmin=(100,'K'), Tmax=(832.77,'K')), NASAPolynomial(coeffs=[18.9415,0.0308354,-1.46794e-05,2.77882e-09,-1.91144e-13,24845.4,-80.821], Tmin=(832.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(233.385,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(six-inringtwodouble-12) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = 'C[C]=C=COC=C[O](29494)',
    structure = SMILES('CC#C[CH]OC=C[O]'),
    E0 = (174.418,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2100,2250,500,550,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.96617,0.0939332,-0.000101528,5.06717e-08,-9.0962e-12,21224.7,35.0043], Tmin=(100,'K'), Tmax=(1672.96,'K')), NASAPolynomial(coeffs=[24.8766,0.0028159,4.32105e-06,-1.13325e-09,8.35645e-14,16013,-97.0432], Tmin=(1672.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(174.418,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CtOsHH) + group(Cs-CtHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(C=COJ) + radical(CCsJOC(O))"""),
)

species(
    label = 'CC=C=[C]OC=C[O](29495)',
    structure = SMILES('C[CH]C#COC=C[O]'),
    E0 = (193.347,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2100,2250,500,550,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.41729,0.0644233,-2.45054e-05,-3.11352e-08,2.21342e-11,23396.8,29.0248], Tmin=(100,'K'), Tmax=(901.306,'K')), NASAPolynomial(coeffs=[21.5054,0.0104306,-5.4617e-07,-1.14298e-10,9.75229e-15,17987.1,-79.4441], Tmin=(901.306,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(193.347,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtOs) + radical(C=COJ) + radical(Sec_Propargyl)"""),
)

species(
    label = 'CC=C=CO[C]=C[O](29496)',
    structure = SMILES('CC=C=CO[C]=C[O]'),
    E0 = (246.781,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,540,610,2055,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.993703,'amu*angstrom^2'), symmetry=1, barrier=(22.8472,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.992653,'amu*angstrom^2'), symmetry=1, barrier=(22.823,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.992396,'amu*angstrom^2'), symmetry=1, barrier=(22.8171,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.295513,0.0715948,-5.67207e-05,1.29812e-08,2.80805e-12,29823.2,31.3991], Tmin=(100,'K'), Tmax=(969.472,'K')), NASAPolynomial(coeffs=[18.8595,0.01677,-5.57646e-06,9.73491e-10,-6.83506e-14,25200.7,-62.8623], Tmin=(969.472,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(246.781,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=CJO) + radical(C=COJ)"""),
)

species(
    label = 'CC=C=COC=[C][O](29497)',
    structure = SMILES('CC=C=CO[CH][C]=O'),
    E0 = (201.055,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,540,610,2055,1855,455,950,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.988811,'amu*angstrom^2'), symmetry=1, barrier=(22.7347,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.990991,'amu*angstrom^2'), symmetry=1, barrier=(22.7848,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.991222,'amu*angstrom^2'), symmetry=1, barrier=(22.7901,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.990199,'amu*angstrom^2'), symmetry=1, barrier=(22.7666,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.747314,0.0882098,-9.50862e-05,4.84505e-08,-9.31735e-12,24366.1,30.7575], Tmin=(100,'K'), Tmax=(1387.25,'K')), NASAPolynomial(coeffs=[24.6251,0.00850139,-1.81763e-06,2.25393e-10,-1.32663e-14,17956.8,-97.6834], Tmin=(1387.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(201.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(CsCJ=O) + radical(CCsJOC(O))"""),
)

species(
    label = '[O]C=COC1[C]=CC1(29498)',
    structure = SMILES('[O]C=COC1[C]=CC1'),
    E0 = (227.881,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.5358,0.0456415,5.52791e-05,-1.29128e-07,5.95499e-11,27560.8,26.9175], Tmin=(100,'K'), Tmax=(925.362,'K')), NASAPolynomial(coeffs=[29.9136,-0.00246227,5.38206e-06,-1.05575e-09,6.06286e-14,18746.3,-130.784], Tmin=(925.362,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(227.881,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(C=COJ) + radical(cyclobutene-vinyl)"""),
)

species(
    label = '[CH2]C=C1[CH]OC=CO1(29499)',
    structure = SMILES('[CH2]C=C1[CH]OC=CO1'),
    E0 = (31.0283,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3514,0.0279037,9.38612e-05,-1.50068e-07,6.0836e-11,3854.73,20.9889], Tmin=(100,'K'), Tmax=(959.689,'K')), NASAPolynomial(coeffs=[22.9047,0.014374,-4.25623e-06,9.40872e-10,-8.42669e-14,-3795.99,-100.413], Tmin=(959.689,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.0283,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cyclohexane) + radical(C=CCJ(O)C) + radical(Allyl_P)"""),
)

species(
    label = '[O]C1[CH]OC=C=CC1(29500)',
    structure = SMILES('[O]C1[CH]OC=C=CC1'),
    E0 = (303.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.435897,0.0592969,-3.67295e-06,-5.14755e-08,2.7584e-11,36591.8,19.828], Tmin=(100,'K'), Tmax=(947.519,'K')), NASAPolynomial(coeffs=[23.0611,0.0111599,-2.46882e-06,4.47164e-10,-3.87801e-14,30177.5,-99.3284], Tmin=(947.519,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(CCsJOC(O)) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC=C=COC=C=O(29501)',
    structure = SMILES('CC=C=COC=C=O'),
    E0 = (47.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.242125,0.0750526,-7.49749e-05,3.79583e-08,-7.53047e-12,5800.54,27.8582], Tmin=(100,'K'), Tmax=(1233.01,'K')), NASAPolynomial(coeffs=[17.4587,0.0192007,-7.02946e-06,1.2217e-09,-8.19535e-14,1554.86,-58.8083], Tmin=(1233.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(47.05,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-(Cdd-O2d)OsH) + group(Cdd-CdsCds)"""),
)

species(
    label = '[CH2]C=C(C=O)C=C[O](29502)',
    structure = SMILES('[CH2]C=C(C=O)C=C[O]'),
    E0 = (-5.33216,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.388752,0.061861,-1.4261e-05,-3.82858e-08,2.22243e-11,-495.354,26.8281], Tmin=(100,'K'), Tmax=(957.945,'K')), NASAPolynomial(coeffs=[22.4207,0.0124481,-3.56796e-06,6.77586e-10,-5.47395e-14,-6670.31,-88.7151], Tmin=(957.945,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-5.33216,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=COJ) + radical(C=CC=CCJ)"""),
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
    label = '[CH]=C=COC=C[O](29503)',
    structure = SMILES('[CH]=C=COC=C[O]'),
    E0 = (197.539,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,610,2055,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.48062,'amu*angstrom^2'), symmetry=1, barrier=(34.0423,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.47852,'amu*angstrom^2'), symmetry=1, barrier=(33.9942,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.78125,0.0547744,-1.51631e-05,-4.17005e-08,2.68438e-11,23889.5,25.3508], Tmin=(100,'K'), Tmax=(896.317,'K')), NASAPolynomial(coeffs=[23.3378,-0.00144465,4.54215e-06,-1.03558e-09,7.14115e-14,18060.6,-90.9615], Tmin=(896.317,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(197.539,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=C=CJ)"""),
)

species(
    label = '[CH2]C1[C]=COC1C=O(29504)',
    structure = SMILES('[CH2]C1[C]=COC1C=O'),
    E0 = (166.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44109,0.0313614,7.01584e-05,-1.31629e-07,5.89695e-11,20082.5,24.782], Tmin=(100,'K'), Tmax=(900.615,'K')), NASAPolynomial(coeffs=[22.2636,0.00502872,3.84408e-06,-9.87759e-10,6.68174e-14,13649.2,-88.3886], Tmin=(900.615,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(166.015,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + ring(2,3-Dihydrofuran) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = 'C=C[C]=COC[C]=O(29505)',
    structure = SMILES('C=C[C]=COC[C]=O'),
    E0 = (153.342,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,1855,455,950,2950,3100,1380,975,1025,1650,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.22447,'amu*angstrom^2'), symmetry=1, barrier=(28.153,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23076,'amu*angstrom^2'), symmetry=1, barrier=(28.2975,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.22956,'amu*angstrom^2'), symmetry=1, barrier=(28.27,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.22681,'amu*angstrom^2'), symmetry=1, barrier=(28.2067,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.06614,0.0855102,-8.67378e-05,4.15346e-08,-7.36388e-12,18647.5,32.0096], Tmin=(100,'K'), Tmax=(1613.28,'K')), NASAPolynomial(coeffs=[24.1306,0.00766407,-6.45389e-08,-1.88482e-10,1.6976e-14,12518.2,-95.4007], Tmin=(1613.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(153.342,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CsCJ=O) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C[C]=[C]OCC=O(29506)',
    structure = SMILES('C=C[C]=[C]OCC=O'),
    E0 = (238.416,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1670,1700,300,440,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.03945,'amu*angstrom^2'), symmetry=1, barrier=(23.8991,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03958,'amu*angstrom^2'), symmetry=1, barrier=(23.9019,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03907,'amu*angstrom^2'), symmetry=1, barrier=(23.8902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03996,'amu*angstrom^2'), symmetry=1, barrier=(23.9106,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.3207,0.0802217,-8.16461e-05,4.07932e-08,-7.76815e-12,28843,31.7677], Tmin=(100,'K'), Tmax=(1406,'K')), NASAPolynomial(coeffs=[21.0899,0.0130146,-3.22973e-06,4.27024e-10,-2.44414e-14,23444.5,-76.6087], Tmin=(1406,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(238.416,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CJO)"""),
)

species(
    label = 'C=[C][C]=COCC=O(29507)',
    structure = SMILES('C=[C][C]=COCC=O'),
    E0 = (197.667,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1670,1700,300,440,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.40159,'amu*angstrom^2'), symmetry=1, barrier=(32.2253,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.39132,'amu*angstrom^2'), symmetry=1, barrier=(31.9893,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.39344,'amu*angstrom^2'), symmetry=1, barrier=(32.0379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.39583,'amu*angstrom^2'), symmetry=1, barrier=(32.0928,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.837548,0.0877692,-9.3842e-05,4.83709e-08,-9.32459e-12,23964,30.0868], Tmin=(100,'K'), Tmax=(1450.97,'K')), NASAPolynomial(coeffs=[23.0317,0.00989445,-8.5512e-07,-8.75415e-11,1.27647e-14,18308.1,-89.574], Tmin=(1450.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(197.667,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=C[C]=COCC=O(29508)',
    structure = SMILES('[CH]C=C=COCC=O'),
    E0 = (223.142,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,540,610,2055,2782.5,750,1395,475,1775,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.065405,0.072645,-3.88139e-05,-1.07342e-08,1.17119e-11,26991.9,28.7895], Tmin=(100,'K'), Tmax=(970.352,'K')), NASAPolynomial(coeffs=[20.4793,0.0203172,-7.11696e-06,1.28611e-09,-9.2433e-14,21531.9,-76.8014], Tmin=(970.352,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(223.142,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=C=C=COCC=O(29509)',
    structure = SMILES('C=C=C=COCC=O'),
    E0 = (29.0615,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.108097,0.0723279,-4.81426e-05,-1.81596e-09,9.1432e-12,3647.39,26.581], Tmin=(100,'K'), Tmax=(965.037,'K')), NASAPolynomial(coeffs=[21.7279,0.0133876,-4.2039e-06,7.65057e-10,-5.741e-14,-1953.63,-84.3525], Tmin=(965.037,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(29.0615,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=CC1=COC1C=O(29510)',
    structure = SMILES('C=CC1=COC1C=O'),
    E0 = (-68.299,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.755529,0.0413126,6.21288e-05,-1.31919e-07,5.93856e-11,-8069.87,23.7633], Tmin=(100,'K'), Tmax=(929.148,'K')), NASAPolynomial(coeffs=[28.2466,0.000154313,3.95795e-06,-7.68814e-10,4.02465e-14,-16510.5,-124.775], Tmin=(929.148,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-68.299,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Cyclobutene)"""),
)

species(
    label = '[CH2]C=[C]C1OC1C=O(29511)',
    structure = SMILES('[CH2]C=[C]C1OC1C=O'),
    E0 = (236.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12257,0.0471712,1.57205e-05,-6.61625e-08,3.31727e-11,28613.2,28.9703], Tmin=(100,'K'), Tmax=(901.336,'K')), NASAPolynomial(coeffs=[18.316,0.0142372,-1.6436e-06,6.27266e-11,-2.34126e-15,23752.2,-61.9641], Tmin=(901.336,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(236.917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Ethylene_oxide) + radical(Cds_S) + radical(Allyl_P)"""),
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
    E0 = (153.251,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (258.886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (262.157,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (412.115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (395.165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (437.551,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (288.395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (367.715,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (390.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (414.312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (423.677,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (401.878,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (588.305,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (186.291,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (315.359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (341.022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (387.779,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (471.435,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (432.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (290.901,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (335.044,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (191.326,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (188.27,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (192.476,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (178.224,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (161.452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (467.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (710.682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (232.595,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (341.531,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (388.886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (427.823,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (389.351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (281.911,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (271.733,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (266.269,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (303.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (178.224,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (467.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (579.102,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (245.171,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (248.041,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (286.916,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (397.826,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (321.304,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (267.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (192.476,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (160.782,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (261.407,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C[C]=CO[CH]C=O(29462)'],
    products = ['OCHCHO(48)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C[C]=CO[CH]C=O(29462)'],
    products = ['C=CC1=COC1[CH][O](29470)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(105.634,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS_D;doublebond_intra;radadd_intra_cdsingle] for rate rule [R5_DS_D;doublebond_intra;radadd_intra_cdsingleDe]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 104.6 to 105.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=C[C]=CO[CH]C=O(29462)'],
    products = ['[CH2]C1[C]=COC=CO1(29471)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.19557e+09,'s^-1'), n=0.501475, Ea=(108.906,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7plus;doublebond_intra_2H_pri;radadd_intra] for rate rule [R8;doublebond_intra_2H_pri;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C=C=C=COC=C[O](29472)'],
    products = ['C=C[C]=CO[CH]C=O(29462)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.23e+08,'cm^3/(mol*s)'), n=1.429, Ea=(16.6816,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 193 used for Ca_Ca;HJ
Exact match found for rate rule [Ca_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C=CC#COC=C[O](29473)'],
    products = ['C=C[C]=CO[CH]C=O(29462)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;HJ] for rate rule [Ct-O_Ct;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', 'C=C[C]=COC=C=O(29474)'],
    products = ['C=C[C]=CO[CH]C=O(29462)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.54607e+06,'m^3/(mol*s)'), n=0.466452, Ea=(32.4943,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C=C[O](2537)', 'CH2CHCCH(26391)'],
    products = ['C=C[C]=CO[CH]C=O(29462)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.04e+06,'m^3/(mol*s)'), n=0, Ea=(33.0536,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;OJ_sec] for rate rule [Ct-H_Ct-Cd;O_rad/OneDe]
Euclidian distance = 2.82842712475
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=C[C]=CO[CH]C=O(29462)'],
    products = ['C=[C]C=COC=C[O](29475)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.806e+09,'s^-1'), n=1.172, Ea=(214.463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=CC=[C]OC=C[O](29476)'],
    products = ['C=C[C]=CO[CH]C=O(29462)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.75e+11,'s^-1'), n=0.633, Ea=(196.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_singleNd;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=C[C]=COC=[C]O(29477)'],
    products = ['C=C[C]=CO[CH]C=O(29462)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=CC=COC=C[O](29478)'],
    products = ['C=C[C]=CO[CH]C=O(29462)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(191.5,'s^-1'), n=3.05, Ea=(222.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=C[C]=CO[C]=CO(29479)'],
    products = ['C=C[C]=CO[CH]C=O(29462)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleNd;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=CC=CO[C]=C[O](29480)'],
    products = ['C=C[C]=CO[CH]C=O(29462)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(9.01194e+11,'s^-1'), n=1.09397, Ea=(394.305,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSD;Cd_rad_out;Cd_H_out_single] for rate rule [R4H_SSD;Cd_rad_out_Cd;Cd_H_out_singleDe]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=C[C]=CO[CH]C=O(29462)'],
    products = ['C=CC=COC=[C][O](29481)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_single;XH_out] for rate rule [R5H_DSSD;Cd_rad_out_singleDe_Cd;XH_out]
Euclidian distance = 3.60555127546
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=C[C]=[C]OC=CO(29482)'],
    products = ['C=C[C]=CO[CH]C=O(29462)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(313415,'s^-1'), n=1.7968, Ea=(63.8264,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5H_RSMS;Cd_rad_out;XH_out] + [R5H_SSMS;Y_rad_out;XH_out] for rate rule [R5H_SSMS;Cd_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=[C][C]=COC=CO(29483)'],
    products = ['C=C[C]=CO[CH]C=O(29462)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.81182e+08,'s^-1'), n=1.25566, Ea=(151.659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;XH_out] for rate rule [R7HJ_1;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C[C]=COC=CO(29484)'],
    products = ['C=C[C]=CO[CH]C=O(29462)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R8Hall;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CHCHO(47)', 'C=C[C]=C[O](26659)'],
    products = ['C=C[C]=CO[CH]C=O(29462)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.63841e+06,'m^3/(mol*s)'), n=0.151, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;O_rad/OneDe] for rate rule [Cd_rad;O_rad/OneDe]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]C=C[O](2537)', '[CH]=[C]C=C(4699)'],
    products = ['C=C[C]=CO[CH]C=O(29462)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.27681e+06,'m^3/(mol*s)'), n=0.151, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [O_rad/OneDe;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C[C]=CO[CH]C=O(29462)'],
    products = ['[O]C=COC=C1[CH]C1(29485)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(5.605e+12,'s^-1'), n=0.275, Ea=(137.649,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 133 used for R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble
Exact match found for rate rule [R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=C[C]=CO[CH]C=O(29462)'],
    products = ['C=C[C]=COC1[CH]O1(29486)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(181.793,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri_HNd_O;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=C[C]=CO[CH]C=O(29462)'],
    products = ['C=CC1=CO[CH]C1[O](29487)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.80657e+08,'s^-1'), n=0.835, Ea=(38.0744,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS_D;doublebond_intra_pri;radadd_intra_cdsingle] for rate rule [R5_DS_D;doublebond_intra_pri;radadd_intra_cdsingleDe]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C[C]=CO[CH]C=O(29462)'],
    products = ['[C]1[CH]COC=COC=1(29488)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.06679e+09,'s^-1'), n=0.473387, Ea=(35.0186,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri_2H;radadd_intra] for rate rule [R8_linear;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=C[C]=CO[CH]C=O(29462)'],
    products = ['C=CC=COC=C=O(29489)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=C[C]=CO[CH]C=O(29462)'],
    products = ['C=C=C=COC=CO(29490)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=C[C]=CO[CH]C=O(29462)'],
    products = ['C=CC1=COC=CO1(29491)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.53377e+11,'s^-1'), n=0.0685, Ea=(8.20064,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad_out;Ypri_rad_out] for rate rule [R6;CdsingleDe_rad_out;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=C[C]=CO[CH]C=O(29462)'],
    products = ['C=C[C]=CC([O])C=O(29458)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction28',
    reactants = ['O(4)', '[CH]=COC=[C]C=C(27408)'],
    products = ['C=C[C]=CO[CH]C=O(29462)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=C[C]=CO[CH]C=O(29462)'],
    products = ['[CH2]C=[C]C1OC=CO1(29492)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(3.49749e+08,'s^-1'), n=0.656505, Ea=(79.3435,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SMS_D;doublebond_intra;radadd_intra] for rate rule [R6_SMS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=C[C]=CO[CH]C=O(29462)'],
    products = ['[O][CH]C1CC=C=CO1(29493)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1e+11,'s^-1'), n=0.21, Ea=(188.28,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R7_SMMS;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C[C]=C=COC=C[O](29494)'],
    products = ['C=C[C]=CO[CH]C=O(29462)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(606700,'s^-1'), n=2.347, Ea=(214.468,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 119 used for R2H_S;Cd_rad_out_double;Cs_H_out_2H
Exact match found for rate rule [R2H_S;Cd_rad_out_double;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['CC=C=[C]OC=C[O](29495)'],
    products = ['C=C[C]=CO[CH]C=O(29462)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(3.992e-05,'s^-1'), n=4.805, Ea=(234.476,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_MMS;Cd_rad_out_single;Cs_H_out_2H] for rate rule [R4H_MMS;Cd_rad_out_singleNd;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['CC=C=CO[C]=C[O](29496)'],
    products = ['C=C[C]=CO[CH]C=O(29462)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(5.59786e+07,'s^-1'), n=1.58088, Ea=(142.57,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_2H] for rate rule [R6H;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['CC=C=COC=[C][O](29497)'],
    products = ['C=C[C]=CO[CH]C=O(29462)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(244756,'s^-1'), n=1.235, Ea=(80.8557,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7H;Y_rad_out;Cs_H_out_2H] for rate rule [R7H;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=C[C]=CO[CH]C=O(29462)'],
    products = ['[O]C=COC1[C]=CC1(29498)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.53664e+10,'s^-1'), n=0.43543, Ea=(118.482,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra_cs2H] + [R4;doublebond_intra;radadd_intra_cs2H] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C=C[C]=CO[CH]C=O(29462)'],
    products = ['[CH2]C=C1[CH]OC=CO1(29499)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(3.68243e+09,'s^-1'), n=0.4695, Ea=(113.018,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SDS_D;doublebond_intra;radadd_intra] for rate rule [R6_SDS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=C[C]=CO[CH]C=O(29462)'],
    products = ['[O]C1[CH]OC=C=CC1(29500)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(163933,'s^-1'), n=1.17455, Ea=(149.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R7_linear;doublebond_intra_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 148.3 to 149.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction38',
    reactants = ['C=C[C]=CO[CH]C=O(29462)'],
    products = ['CC=C=COC=C=O(29501)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C=C[C]=CO[CH]C=O(29462)'],
    products = ['[CH2]C=C(C=O)C=C[O](29502)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction40',
    reactants = ['CH2(19)', '[CH]=C=COC=C[O](29503)'],
    products = ['C=C[C]=CO[CH]C=O(29462)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction41',
    reactants = ['C=C[C]=CO[CH]C=O(29462)'],
    products = ['[CH2]C1[C]=COC1C=O(29504)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(8.3872e+08,'s^-1'), n=0.442111, Ea=(91.9196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_2H_pri;radadd_intra_csHDe] + [R6;doublebond_intra_2H_pri;radadd_intra_cs] for rate rule [R6;doublebond_intra_2H_pri;radadd_intra_csHDe]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction42',
    reactants = ['OCHCHO(48)', '[CH]=[C]C=C(4699)'],
    products = ['C=C[C]=CO[CH]C=O(29462)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.2635,'m^3/(mol*s)'), n=2.48, Ea=(21.7568,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Od_CO-DeH;YJ] for rate rule [Od_CO-DeH;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction43',
    reactants = ['C=C[C]=COC[C]=O(29505)'],
    products = ['C=C[C]=CO[CH]C=O(29462)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(0.00141351,'s^-1'), n=4.515, Ea=(133.574,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/NonDeO] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['C=C[C]=[C]OCC=O(29506)'],
    products = ['C=C[C]=CO[CH]C=O(29462)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(3.32e+07,'s^-1'), n=1.69, Ea=(159.41,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_O;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R3H_SS_O;Cd_rad_out;Cs_H_out_H/CO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['C=[C][C]=COCC=O(29507)'],
    products = ['C=C[C]=CO[CH]C=O(29462)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(2.04875e+10,'s^-1'), n=0.94, Ea=(123.637,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_H/OneDe] for rate rule [R5HJ_1;Cd_rad_out_Cd;Cs_H_out_H/CO]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH]=C[C]=COCC=O(29508)'],
    products = ['C=C[C]=CO[CH]C=O(29462)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/OneDe] for rate rule [R6HJ_2;Cd_rad_out_singleH;Cs_H_out_H/CO]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['C=C[C]=CO[CH]C=O(29462)'],
    products = ['C=C=C=COCC=O(29509)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction48',
    reactants = ['C=C[C]=CO[CH]C=O(29462)'],
    products = ['C=CC1=COC1C=O(29510)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/OneDe;CdsinglepriDe_rad_out]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

reaction(
    label = 'reaction49',
    reactants = ['C=C[C]=CO[CH]C=O(29462)'],
    products = ['[CH2]C=[C]C1OC1C=O(29511)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(108.156,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_csHDe]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

network(
    label = '4680',
    isomers = [
        'C=C[C]=CO[CH]C=O(29462)',
    ],
    reactants = [
        ('OCHCHO(48)', 'CH2CHCCH(26391)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4680',
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

