species(
    label = 'C#C[CH]CO[CH]C=O(29460)',
    structure = SMILES('[CH]=C=CCOC=C[O]'),
    E0 = (151.777,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,610,2055,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.27601,'amu*angstrom^2'), symmetry=1, barrier=(29.3379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.27732,'amu*angstrom^2'), symmetry=1, barrier=(29.3681,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.27675,'amu*angstrom^2'), symmetry=1, barrier=(29.355,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0113538,0.0692281,-2.33584e-05,-4.2182e-08,2.8141e-11,18415.5,28.8598], Tmin=(100,'K'), Tmax=(908.122,'K')), NASAPolynomial(coeffs=[26.2347,0.00364152,2.52032e-06,-6.48694e-10,4.3257e-14,11594.4,-106.459], Tmin=(908.122,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(151.777,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=COJ)"""),
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
    label = 'C#CC1CO[CH][CH]O1(29536)',
    structure = SMILES('C#CC1CO[CH][CH]O1'),
    E0 = (255.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.28772,0.0858222,-9.35724e-05,5.15897e-08,-1.01284e-11,30900.6,27.6842], Tmin=(100,'K'), Tmax=(1590.7,'K')), NASAPolynomial(coeffs=[15.38,0.0131353,3.98912e-06,-1.46065e-09,1.20625e-13,29491.4,-48.2266], Tmin=(1590.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(255.121,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(1,4-Dioxane) + radical(CCsJOCs) + radical(CCsJOCs)"""),
)

species(
    label = 'O=CC1C=[C][CH]CO1(29577)',
    structure = SMILES('O=CC1C=[C][CH]CO1'),
    E0 = (111.187,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28517,0.0430298,2.42152e-05,-6.61951e-08,2.96374e-11,13485.6,21.5188], Tmin=(100,'K'), Tmax=(951.318,'K')), NASAPolynomial(coeffs=[16.3075,0.0212479,-6.68988e-06,1.18865e-09,-8.70377e-14,8754.86,-60.0474], Tmin=(951.318,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(111.187,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(36dihydro2hpyran) + radical(C=CCJCO) + radical(Cds_S)"""),
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
    label = '[CH]=C=CCOC=C=O(29606)',
    structure = SMILES('[CH]=C=CCOC=C=O'),
    E0 = (191.791,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,3120,650,792.5,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.1154,'amu*angstrom^2'), symmetry=1, barrier=(25.6452,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12309,'amu*angstrom^2'), symmetry=1, barrier=(25.822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12084,'amu*angstrom^2'), symmetry=1, barrier=(25.7703,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.250117,0.0794378,-8.63548e-05,4.54903e-08,-9.02801e-12,23232.1,29.1901], Tmin=(100,'K'), Tmax=(1378.71,'K')), NASAPolynomial(coeffs=[20.899,0.00975784,-1.492e-06,7.79359e-11,3.23397e-16,18191.3,-76.766], Tmin=(1378.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.791,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)OsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=C=CCOC=[C]O(29607)',
    structure = SMILES('[CH]=C=CCOC=[C]O'),
    E0 = (250.059,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,2750,2850,1437.5,1250,1305,750,350,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.06739,'amu*angstrom^2'), symmetry=1, barrier=(24.5413,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06669,'amu*angstrom^2'), symmetry=1, barrier=(24.5254,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06646,'amu*angstrom^2'), symmetry=1, barrier=(24.52,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0669,'amu*angstrom^2'), symmetry=1, barrier=(24.5302,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.34644,0.0915336,-0.000101013,5.19318e-08,-9.75407e-12,30290,35.455], Tmin=(100,'K'), Tmax=(1548.15,'K')), NASAPolynomial(coeffs=[25.0761,0.00392251,2.61471e-06,-7.62659e-10,5.84161e-14,24426.8,-96.0805], Tmin=(1548.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.059,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJO) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=C=CCO[C]=CO(29608)',
    structure = SMILES('[CH]=C=CCO[C]=CO'),
    E0 = (250.059,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,2750,2850,1437.5,1250,1305,750,350,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.06739,'amu*angstrom^2'), symmetry=1, barrier=(24.5413,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06669,'amu*angstrom^2'), symmetry=1, barrier=(24.5254,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06646,'amu*angstrom^2'), symmetry=1, barrier=(24.52,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0669,'amu*angstrom^2'), symmetry=1, barrier=(24.5302,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.34644,0.0915336,-0.000101013,5.19318e-08,-9.75407e-12,30290,35.455], Tmin=(100,'K'), Tmax=(1548.15,'K')), NASAPolynomial(coeffs=[25.0761,0.00392251,2.61471e-06,-7.62659e-10,5.84161e-14,24426.8,-96.0805], Tmin=(1548.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.059,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJO) + radical(C=C=CJ)"""),
)

species(
    label = 'C=C=[C]COC=C[O](29609)',
    structure = SMILES('[CH2]C#CCOC=C[O]'),
    E0 = (136.898,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.230311,0.063404,-8.95481e-06,-5.55231e-08,3.25194e-11,16619,28.6815], Tmin=(100,'K'), Tmax=(909.104,'K')), NASAPolynomial(coeffs=[25.5335,0.0043074,2.36552e-06,-6.21165e-10,4.09186e-14,9859.8,-102.852], Tmin=(909.104,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(136.898,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CtOsHH) + group(Cs-CtHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(C=COJ) + radical(Propargyl)"""),
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
    label = '[CH]=C=C[CH]OC=CO(29610)',
    structure = SMILES('[CH]=C=C[CH]OC=CO'),
    E0 = (121.254,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.90862,0.0937154,-9.96879e-05,4.87816e-08,-8.63426e-12,14827.6,35.3044], Tmin=(100,'K'), Tmax=(1678.81,'K')), NASAPolynomial(coeffs=[25.8549,0.00293551,3.42917e-06,-9.05658e-10,6.62488e-14,8976.44,-102.685], Tmin=(1678.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(121.254,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CCJ(O)C)"""),
)

species(
    label = '[CH]=C=[C]COC=CO(29611)',
    structure = SMILES('[CH]=C=[C]COC=CO'),
    E0 = (248.156,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,2750,2850,1437.5,1250,1305,750,350,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.16021,'amu*angstrom^2'), symmetry=1, barrier=(26.6755,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1606,'amu*angstrom^2'), symmetry=1, barrier=(26.6844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15959,'amu*angstrom^2'), symmetry=1, barrier=(26.6613,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15841,'amu*angstrom^2'), symmetry=1, barrier=(26.6342,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.86309,0.0975596,-0.000109905,5.65613e-08,-1.05251e-11,30084.5,34.6578], Tmin=(100,'K'), Tmax=(1584.46,'K')), NASAPolynomial(coeffs=[27.0862,0.000546037,4.59186e-06,-1.14581e-09,8.40802e-14,23914.6,-108.85], Tmin=(1584.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(248.156,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=C=CJ)"""),
)

species(
    label = 'C=C=CCO[C]=C[O](29612)',
    structure = SMILES('C=C=CCO[C]=C[O]'),
    E0 = (237.045,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,540,610,2055,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.04865,'amu*angstrom^2'), symmetry=1, barrier=(24.1106,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05108,'amu*angstrom^2'), symmetry=1, barrier=(24.1663,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04865,'amu*angstrom^2'), symmetry=1, barrier=(24.1106,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.205503,0.0700545,-4.27439e-05,-8.28289e-09,1.20046e-11,28658.8,30.6631], Tmin=(100,'K'), Tmax=(945.198,'K')), NASAPolynomial(coeffs=[21.6075,0.0122963,-3.1571e-06,5.2431e-10,-3.93607e-14,23147.2,-79.1367], Tmin=(945.198,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(237.045,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJO) + radical(C=COJ)"""),
)

species(
    label = 'C=C=CCOC=[C][O](29613)',
    structure = SMILES('C=C=CCOC=[C][O]'),
    E0 = (237.045,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,540,610,2055,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.04865,'amu*angstrom^2'), symmetry=1, barrier=(24.1106,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05108,'amu*angstrom^2'), symmetry=1, barrier=(24.1663,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04865,'amu*angstrom^2'), symmetry=1, barrier=(24.1106,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.205503,0.0700545,-4.27439e-05,-8.28289e-09,1.20046e-11,28658.8,30.6631], Tmin=(100,'K'), Tmax=(945.198,'K')), NASAPolynomial(coeffs=[21.6075,0.0122963,-3.1571e-06,5.2431e-10,-3.93607e-14,23147.2,-79.1367], Tmin=(945.198,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(237.045,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJO) + radical(C=COJ)"""),
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
    label = '[CH]=C=CC[O](26927)',
    structure = SMILES('[CH]=C=CC[O]'),
    E0 = (373.748,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3120,650,792.5,1650,3010,987.5,1337.5,450,1655,540,610,2055,180],'cm^-1')),
        HinderedRotor(inertia=(0.175448,'amu*angstrom^2'), symmetry=1, barrier=(4.03389,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.52908,0.0362944,-4.67924e-05,4.20129e-08,-1.53469e-11,45000.7,18.5473], Tmin=(100,'K'), Tmax=(852.445,'K')), NASAPolynomial(coeffs=[3.23392,0.0243897,-1.07161e-05,1.96772e-09,-1.32886e-13,45192.9,17.0916], Tmin=(852.445,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(373.748,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCOJ) + radical(C=C=CJ)"""),
)

species(
    label = 'C3H2(81)',
    structure = SMILES('[CH]C#C'),
    E0 = (530.398,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.621997,'amu*angstrom^2'), symmetry=1, barrier=(14.3009,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (38.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.6874,0.0353211,-7.31245e-05,7.03801e-08,-2.45313e-11,63833.3,8.93418], Tmin=(100,'K'), Tmax=(908.454,'K')), NASAPolynomial(coeffs=[4.78002,0.0100611,-4.92178e-06,8.868e-10,-5.66859e-14,64115.2,2.68365], Tmin=(908.454,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(530.398,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""HCCCH(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]OC=C[O](2549)',
    structure = SMILES('[CH2]OC=C[O]'),
    E0 = (-23.2817,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.17869,'amu*angstrom^2'), symmetry=1, barrier=(27.1004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17224,'amu*angstrom^2'), symmetry=1, barrier=(26.9521,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.0627,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87374,0.0271007,3.67127e-05,-8.65645e-08,4.08996e-11,-2705.16,18.2416], Tmin=(100,'K'), Tmax=(908.37,'K')), NASAPolynomial(coeffs=[21.3383,-0.00659169,6.44909e-06,-1.31006e-09,8.52988e-14,-8387.52,-85.6062], Tmin=(908.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-23.2817,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-OsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COCJ) + radical(C=COJ)"""),
)

species(
    label = '[CH]=C=CCOC1[CH]O1(29614)',
    structure = SMILES('[CH]=C=CCOC1[CH]O1'),
    E0 = (304.296,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0500545,0.0784315,-8.5878e-05,4.95855e-08,-1.09516e-11,36747.9,30.1437], Tmin=(100,'K'), Tmax=(1256.34,'K')), NASAPolynomial(coeffs=[15.8761,0.0189484,-3.99945e-06,3.74975e-10,-1.25333e-14,33489.1,-46.9628], Tmin=(1256.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(304.296,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(Ethylene_oxide) + radical(C=C=CJ) + radical(CCsJO)"""),
)

species(
    label = '[O]C=COCC1[C]=C1(29615)',
    structure = SMILES('[O]C=COCC1[C]=C1'),
    E0 = (313.207,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.178295,0.075374,-4.27761e-05,-1.94035e-08,1.88784e-11,37835.8,25.9573], Tmin=(100,'K'), Tmax=(923.451,'K')), NASAPolynomial(coeffs=[25.8736,0.00559668,6.08232e-07,-2.19185e-10,1.18949e-14,31187.9,-107.596], Tmin=(923.451,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(313.207,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cyclopropene) + radical(C=COJ) + radical(cyclopropenyl-vinyl)"""),
)

species(
    label = '[CH]=C1[CH]COC=CO1(29552)',
    structure = SMILES('[CH]C1=CCOC=CO1'),
    E0 = (169.444,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06281,-0.00208665,0.000262761,-3.95898e-07,1.69092e-10,20547.9,25.6518], Tmin=(100,'K'), Tmax=(899.423,'K')), NASAPolynomial(coeffs=[48.4911,-0.0375358,2.92291e-05,-5.88279e-09,3.90875e-13,4918.54,-237.592], Tmin=(899.423,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(169.444,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[O]C1[CH]OCC=C=C1(29587)',
    structure = SMILES('[O]C1[CH]OCC=C=C1'),
    E0 = (334.695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.984707,0.056163,-2.7075e-05,-3.59307e-09,4.71129e-12,40371.7,21.9046], Tmin=(100,'K'), Tmax=(1092.06,'K')), NASAPolynomial(coeffs=[14.35,0.0257733,-1.08329e-05,2.05836e-09,-1.46087e-13,36345.6,-48.8206], Tmin=(1092.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(334.695,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(CC(C)OJ) + radical(CCsJOCs)"""),
)

species(
    label = 'C=C=CCOC=C=O(29616)',
    structure = SMILES('C=C=CCOC=C=O'),
    E0 = (37.3138,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.233707,0.0725979,-5.80688e-05,1.33652e-08,2.84959e-12,4632.61,26.8268], Tmin=(100,'K'), Tmax=(971.687,'K')), NASAPolynomial(coeffs=[19.4342,0.0159932,-5.32133e-06,9.37305e-10,-6.64335e-14,-157.92,-70.7032], Tmin=(971.687,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(37.3138,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)OsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C#C[CH]CC([O])C=O(29456)',
    structure = SMILES('C#C[CH]CC([O])C=O'),
    E0 = (239.283,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.723544,0.0675002,-6.3529e-05,3.07116e-08,-5.42362e-12,28901.1,31.6329], Tmin=(100,'K'), Tmax=(975.151,'K')), NASAPolynomial(coeffs=[13.4857,0.0232898,-8.04314e-06,1.33747e-09,-8.72547e-14,26025.2,-31.6001], Tmin=(975.151,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(239.283,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CtCsHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(C=OCOJ)"""),
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
    label = '[CH]=C=CCOC=[CH](27456)',
    structure = SMILES('[CH]=C=CCOC=[CH]'),
    E0 = (466.204,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3115,3125,620,680,785,800,1600,1700,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.17957,'amu*angstrom^2'), symmetry=1, barrier=(27.1205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17736,'amu*angstrom^2'), symmetry=1, barrier=(27.0697,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17998,'amu*angstrom^2'), symmetry=1, barrier=(27.1301,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.347631,0.0745295,-7.52097e-05,3.65084e-08,-6.58786e-12,56246,29.2933], Tmin=(100,'K'), Tmax=(1582.58,'K')), NASAPolynomial(coeffs=[20.604,0.00864841,-5.15284e-07,-1.17639e-10,1.31857e-14,51233.1,-76.2905], Tmin=(1582.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(466.204,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(Cds_P)"""),
)

species(
    label = 'C#CC1COC1[CH][O](29617)',
    structure = SMILES('C#CC1COC1[CH][O]'),
    E0 = (350.633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.623604,0.0692571,-7.16461e-05,4.25148e-08,-9.88535e-12,42297.3,26.4271], Tmin=(100,'K'), Tmax=(1196.62,'K')), NASAPolynomial(coeffs=[11.5286,0.0250483,-6.50669e-06,8.07404e-10,-4.01459e-14,40242.7,-25.8205], Tmin=(1196.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(350.633,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CtCsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Oxetane) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = 'C#CC=COC=C[O](29618)',
    structure = SMILES('C#CC=COC=C[O]'),
    E0 = (131.042,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,750,770,3400,2100,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.43655,'amu*angstrom^2'), symmetry=1, barrier=(33.0292,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43854,'amu*angstrom^2'), symmetry=1, barrier=(33.0749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43442,'amu*angstrom^2'), symmetry=1, barrier=(32.9801,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.23459,0.0625755,-1.09378e-05,-5.2547e-08,3.08262e-11,15915,27.3541], Tmin=(100,'K'), Tmax=(923.222,'K')), NASAPolynomial(coeffs=[26.4392,0.00165515,2.55583e-06,-5.60517e-10,3.2767e-14,8834.2,-109.119], Tmin=(923.222,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(131.042,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsCtH) + group(Cds-CdsOsH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=COJ)"""),
)

species(
    label = 'C#CC[CH]OC=C[O](29529)',
    structure = SMILES('C#CC[CH]OC=C[O]'),
    E0 = (187.415,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,2175,525,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.05974,0.0998693,-0.000112185,5.72268e-08,-1.05444e-11,22787.8,33.9963], Tmin=(100,'K'), Tmax=(1601.27,'K')), NASAPolynomial(coeffs=[28.0803,-0.000242192,4.84661e-06,-1.17803e-09,8.54354e-14,16317.5,-115.665], Tmin=(1601.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(187.415,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CtCsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=COJ) + radical(CCsJOC(O))"""),
)

species(
    label = '[C]#CCCOC=C[O](29619)',
    structure = SMILES('[C]#CCCOC=C[O]'),
    E0 = (330.632,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,2175,525,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.43565,0.0913412,-9.88883e-05,5.02781e-08,-9.32396e-12,39985.8,33.0193], Tmin=(100,'K'), Tmax=(1581.99,'K')), NASAPolynomial(coeffs=[24.3604,0.00526931,2.48979e-06,-7.73733e-10,6.00879e-14,34432.7,-95.0191], Tmin=(1581.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(330.632,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CtCsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=COJ) + radical(Acetyl)"""),
)

species(
    label = 'C#CCCO[C]=C[O](29620)',
    structure = SMILES('C#CCCO[C]=C[O]'),
    E0 = (233.232,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,750,770,3400,2100,1685,370,3010,987.5,1337.5,450,1655,2175,525,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.942329,0.0859359,-9.04333e-05,4.53612e-08,-8.42439e-12,28248.8,33.6352], Tmin=(100,'K'), Tmax=(1539.5,'K')), NASAPolynomial(coeffs=[23.2567,0.0079433,2.87711e-07,-3.03136e-10,2.68613e-14,22589.4,-87.7338], Tmin=(1539.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(233.232,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CtCsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = 'C#CCCOC=[C][O](29621)',
    structure = SMILES('C#CCCOC=[C][O]'),
    E0 = (233.232,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,750,770,3400,2100,1685,370,3010,987.5,1337.5,450,1655,2175,525,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.942329,0.0859359,-9.04333e-05,4.53612e-08,-8.42439e-12,28248.8,33.6352], Tmin=(100,'K'), Tmax=(1539.5,'K')), NASAPolynomial(coeffs=[23.2567,0.0079433,2.87711e-07,-3.03136e-10,2.68613e-14,22589.4,-87.7338], Tmin=(1539.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(233.232,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CtCsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = '[C]#C[CH]COC=CO(29622)',
    structure = SMILES('[C]#C[CH]COC=CO'),
    E0 = (389.071,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2175,525,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,3025,407.5,1350,352.5,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.07844,0.0981094,-0.000111231,5.75477e-08,-1.06685e-11,47044,36.4138], Tmin=(100,'K'), Tmax=(1618.65,'K')), NASAPolynomial(coeffs=[26.2353,-6.20956e-05,5.88012e-06,-1.45077e-09,1.0631e-13,41572.6,-102.407], Tmin=(1618.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(389.071,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CtCsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCJCO) + radical(Acetyl)"""),
)

species(
    label = 'C#CC1CO[CH]C1[O](29623)',
    structure = SMILES('C#CC1CO[CH]C1[O]'),
    E0 = (276.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06062,0.0504285,9.54317e-06,-6.6139e-08,3.64543e-11,33433.1,21.5217], Tmin=(100,'K'), Tmax=(855.502,'K')), NASAPolynomial(coeffs=[18.5896,0.0117459,1.48895e-06,-7.3278e-10,6.15529e-14,28850.3,-69.5653], Tmin=(855.502,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(276.984,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CtCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Tetrahydrofuran) + radical(CC(C)OJ) + radical(CCsJOCs)"""),
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
    label = 'C#CC=COC=CO(29624)',
    structure = SMILES('C#CC=COC=CO'),
    E0 = (-10.4209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.181818,0.0695877,-1.43941e-05,-5.88775e-08,3.57231e-11,-1081.83,27.4395], Tmin=(100,'K'), Tmax=(910.794,'K')), NASAPolynomial(coeffs=[29.7875,-0.00219155,5.27058e-06,-1.13672e-09,7.41345e-14,-9022.95,-127.97], Tmin=(910.794,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-10.4209,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsCtH) + group(Cds-CdsOsH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH)"""),
)

species(
    label = 'C#CCCOC=C=O(29625)',
    structure = SMILES('C#CCCOC=C=O'),
    E0 = (33.5012,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.509165,0.0838231,-8.99596e-05,4.70905e-08,-9.28721e-12,4204.76,28.3373], Tmin=(100,'K'), Tmax=(1395.44,'K')), NASAPolynomial(coeffs=[21.5635,0.0111316,-1.69487e-06,8.4476e-11,8.79417e-16,-1038.24,-82.2189], Tmin=(1395.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(33.5012,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CtCsHH) + group(Cds-(Cdd-O2d)OsH) + group(Ct-CtCs) + group(Ct-CtH)"""),
)

species(
    label = 'C#CC1COC=CO1(29530)',
    structure = SMILES('C#CC1COC=CO1'),
    E0 = (-50.0631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.855707,0.0388017,7.56218e-05,-1.60356e-07,7.64853e-11,-5878.99,16.5362], Tmin=(100,'K'), Tmax=(876.013,'K')), NASAPolynomial(coeffs=[29.8153,-0.00909842,1.32367e-05,-2.98379e-09,2.11528e-13,-14188.7,-137.813], Tmin=(876.013,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-50.0631,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(23dihydro14dioxin)"""),
)

species(
    label = '[CH]=C=CCOC[C]=O(29626)',
    structure = SMILES('[CH]=C=CCOC[C]=O'),
    E0 = (219.323,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1855,455,950,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,540,610,2055,261.083,261.163,261.177],'cm^-1')),
        HinderedRotor(inertia=(0.00247279,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.329529,'amu*angstrom^2'), symmetry=1, barrier=(15.9537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18112,'amu*angstrom^2'), symmetry=1, barrier=(8.74272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.35749,'amu*angstrom^2'), symmetry=1, barrier=(65.5791,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.904605,0.0716473,-7.88592e-05,4.96053e-08,-1.29869e-11,26486.9,29.585], Tmin=(100,'K'), Tmax=(915.481,'K')), NASAPolynomial(coeffs=[9.89181,0.0323794,-1.45188e-05,2.7514e-09,-1.91954e-13,24841.4,-12.9793], Tmin=(915.481,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(219.323,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(CsCJ=O)"""),
)

species(
    label = '[CH]=C=C[CH]OCC=O(29627)',
    structure = SMILES('[CH]=C=C[CH]OCC=O'),
    E0 = (175.592,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,327.593,327.593,327.595],'cm^-1')),
        HinderedRotor(inertia=(0.00565071,'amu*angstrom^2'), symmetry=1, barrier=(53.914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.105043,'amu*angstrom^2'), symmetry=1, barrier=(7.9995,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.707947,'amu*angstrom^2'), symmetry=1, barrier=(53.9139,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.707956,'amu*angstrom^2'), symmetry=1, barrier=(53.914,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.091,0.0685688,-7.28211e-05,4.65915e-08,-1.28114e-11,21219.6,29.1763], Tmin=(100,'K'), Tmax=(861.642,'K')), NASAPolynomial(coeffs=[8.02857,0.0363629,-1.67557e-05,3.21323e-09,-2.25577e-13,20024,-3.2604], Tmin=(861.642,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(175.592,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CCJ(O)C)"""),
)

species(
    label = '[CH]=C=[C]COCC=O(29628)',
    structure = SMILES('[CH]=C=[C]COCC=O'),
    E0 = (302.494,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,540,610,2055,3120,650,792.5,1650,2782.5,750,1395,475,1775,1000,1685,370,250.049,250.053,250.054],'cm^-1')),
        HinderedRotor(inertia=(0.212813,'amu*angstrom^2'), symmetry=1, barrier=(9.44209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.212808,'amu*angstrom^2'), symmetry=1, barrier=(9.44214,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00269621,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49817,'amu*angstrom^2'), symmetry=1, barrier=(66.47,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.884369,0.0747405,-8.77274e-05,5.42422e-08,-1.09799e-11,36488.1,29.4747], Tmin=(100,'K'), Tmax=(642.805,'K')), NASAPolynomial(coeffs=[8.99344,0.0340328,-1.54931e-05,2.92951e-09,-2.03029e-13,35244.1,-7.63064], Tmin=(642.805,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(302.494,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C1[CH]COC1C=O(29603)',
    structure = SMILES('[CH]C1=CCOC1C=O'),
    E0 = (128.719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35074,0.0412145,3.21771e-05,-6.87972e-08,2.82179e-11,15591.8,26.2084], Tmin=(100,'K'), Tmax=(992.542,'K')), NASAPolynomial(coeffs=[14.5644,0.0290231,-1.14515e-05,2.18667e-09,-1.59589e-13,10946.3,-47.6294], Tmin=(992.542,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(128.719,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(25dihydrofuran) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[C]#C[CH]COCC=O(29629)',
    structure = SMILES('[C]#C[CH]COCC=O'),
    E0 = (443.409,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2175,525,2782.5,750,1395,475,1775,1000,3025,407.5,1350,352.5,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.694751,0.075491,-9.25889e-05,6.51897e-08,-1.8703e-11,53446.2,31.1084], Tmin=(100,'K'), Tmax=(848.17,'K')), NASAPolynomial(coeffs=[10.295,0.0302152,-1.25165e-05,2.25112e-09,-1.51288e-13,51817.7,-13.6259], Tmin=(848.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(443.409,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsOsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(CCJCO)"""),
)

species(
    label = 'C#CC=COCC=O(29630)',
    structure = SMILES('C#CC=COCC=O'),
    E0 = (-23.538,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.199404,0.0656585,-2.11351e-05,-3.55563e-08,2.25335e-11,-2677.73,26.0636], Tmin=(100,'K'), Tmax=(946.299,'K')), NASAPolynomial(coeffs=[24.0872,0.009047,-1.71813e-06,3.04254e-10,-2.81728e-14,-9185,-98.3575], Tmin=(946.299,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-23.538,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-CdsCtH) + group(Cds-OdCsH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH)"""),
)

species(
    label = 'C#CC1COC1C=O(29465)',
    structure = SMILES('C#CC1COC1C=O'),
    E0 = (24.2467,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54952,0.0370236,4.07892e-05,-9.16829e-08,4.31045e-11,3020.49,24.493], Tmin=(100,'K'), Tmax=(879.973,'K')), NASAPolynomial(coeffs=[16.7803,0.0144513,-2.72212e-07,-3.1676e-10,2.81283e-14,-1466.65,-57.3045], Tmin=(879.973,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(24.2467,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CtCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Oxetane)"""),
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
    E0 = (151.777,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (273.666,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (277.297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (436.077,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (412.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (400.404,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (283.911,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (196.086,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (240.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (332.937,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (482.478,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (472.927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (432.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (619.596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (508.212,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (333.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (313.959,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (213.948,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (334.695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (176.751,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (465.577,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (709.208,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (350.633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (349.604,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (354.084,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (295.827,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (477.909,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (277.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (272.012,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (513.507,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (308.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (217.198,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (169.559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (191.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (159.309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (271.691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (352.897,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (319.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (346.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (308.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (526.819,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (244.871,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (159.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C#C[CH]CO[CH]C=O(29460)'],
    products = ['OCHCHO(48)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C#C[CH]CO[CH]C=O(29460)'],
    products = ['C#CC1CO[CH][CH]O1(29536)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(9.291e+11,'s^-1'), n=0.234, Ea=(121.888,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_SMSR;doublebond_intra;radadd_intra_O] for rate rule [R7_SMSS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C#C[CH]CO[CH]C=O(29460)'],
    products = ['O=CC1C=[C][CH]CO1(29577)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1e+11,'s^-1'), n=0.21, Ea=(125.52,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;doublebond_intra;radadd_intra_cdsingleH] for rate rule [R7_MMSR;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[CH]=C=CCOC=C=O(29606)'],
    products = ['C#C[CH]CO[CH]C=O(29460)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.54607e+06,'m^3/(mol*s)'), n=0.466452, Ea=(32.4943,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=C=CCOC=[C]O(29607)'],
    products = ['C#C[CH]CO[CH]C=O(29460)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=C=CCO[C]=CO(29608)'],
    products = ['C#C[CH]CO[CH]C=O(29460)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleNd;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C#C[CH]CO[CH]C=O(29460)'],
    products = ['C=C=[C]COC=C[O](29609)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(443913,'s^-1'), n=2.07262, Ea=(132.134,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_singleH;Cd_H_out_singleNd] + [R3H;Cd_rad_out;Cd_H_out_singleNd] + [R3H;Cd_rad_out_singleH;XH_out] for rate rule [R3H;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C#C[CH]CO[CH]C=O(29460)'],
    products = ['C=[C]C=COC=C[O](29475)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R4H_MMS;Cd_rad_out_singleH;Cs_H_out_H/NonDeO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C#C[CH]CO[CH]C=O(29460)'],
    products = ['[CH]=C=C[CH]OC=CO(29610)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(252000,'s^-1'), n=1.85, Ea=(88.2824,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SMSS;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R5H_SMSS;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=C=[C]COC=CO(29611)'],
    products = ['C#C[CH]CO[CH]C=O(29460)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(16140,'s^-1'), n=1.92259, Ea=(84.7802,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_RSSMS;Y_rad_out;XH_out] for rate rule [R6H_RSSMS;Cd_rad_out_double;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=C=CCO[C]=C[O](29612)'],
    products = ['C#C[CH]CO[CH]C=O(29460)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.04e+11,'s^-1'), n=0.627, Ea=(245.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cd_H_out_singleH] for rate rule [R6H;Cd_rad_out_Cd;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=C=CCOC=[C][O](29613)'],
    products = ['C#C[CH]CO[CH]C=O(29460)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.36537e+10,'s^-1'), n=1.06641, Ea=(235.882,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out;Cd_H_out_singleH] for rate rule [R7H;Cd_rad_out;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]C=C[O](2537)', '[CH]=[C]C=C(4699)'],
    products = ['C#C[CH]CO[CH]C=O(29460)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(5.1779e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cs_rad;O_sec_rad] + [C_pri_rad;O_rad] for rate rule [C_rad/H2/Cd;O_rad/OneDe]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['CHCHO(47)', '[CH]=C=CC[O](26927)'],
    products = ['C#C[CH]CO[CH]C=O(29460)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.87866e+07,'m^3/(mol*s)'), n=-0.180109, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_rad/NonDe;Y_rad] for rate rule [O_rad/NonDe;Cd_rad]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['C3H2(81)', '[CH2]OC=C[O](2549)'],
    products = ['C#C[CH]CO[CH]C=O(29460)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.28206e+10,'m^3/(mol*s)'), n=-0.9855, Ea=(1.09644,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_allenic;C_pri_rad] for rate rule [Cd_allenic;C_rad/H2/O]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C#C[CH]CO[CH]C=O(29460)'],
    products = ['[CH]=C=CCOC1[CH]O1(29614)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(181.793,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri_HNd_O;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C#C[CH]CO[CH]C=O(29460)'],
    products = ['[O]C=COCC1[C]=C1(29615)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.2354e+12,'s^-1'), n=-0.1205, Ea=(162.181,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra_cdsingleH] + [R3;doublebond_intra_CdCdd;radadd_intra_cdsingle] for rate rule [R3;doublebond_intra_CdCdd;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C#C[CH]CO[CH]C=O(29460)'],
    products = ['[CH]=C1[CH]COC=CO1(29552)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(8.62196e+06,'s^-1'), n=0.867572, Ea=(62.1704,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;doublebond_intra;radadd_intra] for rate rule [R7_linear;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C#C[CH]CO[CH]C=O(29460)'],
    products = ['[O]C1[CH]OCC=C=C1(29587)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.625e+11,'s^-1'), n=0.16, Ea=(182.918,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;doublebond_intra_pri;radadd_intra_cdsingleH] for rate rule [R7_linear;doublebond_intra_pri;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 180.4 to 182.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction20',
    reactants = ['C#C[CH]CO[CH]C=O(29460)'],
    products = ['C=C=CCOC=C=O(29616)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C#C[CH]CO[CH]C=O(29460)'],
    products = ['C#C[CH]CC([O])C=O(29456)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction22',
    reactants = ['O(4)', '[CH]=C=CCOC=[CH](27456)'],
    products = ['C#C[CH]CO[CH]C=O(29460)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['C#C[CH]CO[CH]C=O(29460)'],
    products = ['C#CC1COC1[CH][O](29617)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.06771e+06,'s^-1'), n=1.35044, Ea=(198.855,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra;radadd_intra_cs] for rate rule [R5_SS_D;doublebond_intra;radadd_intra_csHCt]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 197.8 to 198.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', 'C#CC=COC=C[O](29618)'],
    products = ['C#C[CH]CO[CH]C=O(29460)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.182e+10,'cm^3/(mol*s)'), n=0.859, Ea=(6.76971,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [Cds-OsH_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]C=C[O](2537)', 'CH2CHCCH(26391)'],
    products = ['C#C[CH]CO[CH]C=O(29460)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.00026,'m^3/(mol*s)'), n=3.01, Ea=(98.7424,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-HH_Cds-OneDeH;OJ] for rate rule [Cds-HH_Cds-CtH;O_rad/OneDe]
Euclidian distance = 3.16227766017
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C#CC[CH]OC=C[O](29529)'],
    products = ['C#C[CH]CO[CH]C=O(29460)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.23666e+09,'s^-1'), n=1.04335, Ea=(108.412,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_1H;Cs_H_out_H/OneDe] for rate rule [R2H_S;C_rad_out_H/NonDeO;Cs_H_out_H/Ct]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[C]#CCCOC=C[O](29619)'],
    products = ['C#C[CH]CO[CH]C=O(29460)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_MS;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R3H_TS;Ct_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C#CCCO[C]=C[O](29620)'],
    products = ['C#C[CH]CO[CH]C=O(29460)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out_H/OneDe] for rate rule [R4H_SSS_OCs;Cd_rad_out_Cd;Cs_H_out_H/Ct]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C#CCCOC=[C][O](29621)'],
    products = ['C#C[CH]CO[CH]C=O(29460)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(13813.5,'s^-1'), n=1.88327, Ea=(38.7799,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R5H_DSSS;Cd_rad_out;Cs_H_out_H/Ct]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[C]#C[CH]COC=CO(29622)'],
    products = ['C#C[CH]CO[CH]C=O(29460)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(3728.31,'s^-1'), n=2.31462, Ea=(124.436,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Y_rad_out;O_H_out] + [R8Hall;Y_rad_out;XH_out] for rate rule [R8Hall;Ct_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C#C[CH]CO[CH]C=O(29460)'],
    products = ['C#CC1CO[CH]C1[O](29623)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(6.8435e+15,'s^-1'), n=-1.17677, Ea=(156.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra;radadd_intra_csHDe] for rate rule [R5_SS_D;doublebond_intra_pri;radadd_intra_csHCt]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C#C[CH]CO[CH]C=O(29460)'],
    products = ['[C]1[CH]COC=COC=1(29488)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(8.28217e+09,'s^-1'), n=0.356567, Ea=(65.4205,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;triplebond_intra_H;radadd_intra] for rate rule [R8_linear;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C#C[CH]CO[CH]C=O(29460)'],
    products = ['C#CC=COC=CO(29624)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(2.14e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_De] for rate rule [R5radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C#C[CH]CO[CH]C=O(29460)'],
    products = ['C#CCCOC=C=O(29625)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C#C[CH]CO[CH]C=O(29460)'],
    products = ['C#CC1COC=CO1(29530)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), Tmin=(550,'K'), Tmax=(650,'K'), comment="""Estimated using template [R6_SSSDS;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R6_SSSDS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['OCHCHO(48)', '[CH]=[C]C=C(4699)'],
    products = ['C#C[CH]CO[CH]C=O(29460)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(0.000397355,'m^3/(mol*s)'), n=2.77646, Ea=(45.4073,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_R;CsJ-CdHH] for rate rule [Od_CO-DeH;CsJ-CdHH]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]=C=CCOC[C]=O(29626)'],
    products = ['C#C[CH]CO[CH]C=O(29460)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(0.00141351,'s^-1'), n=4.515, Ea=(133.574,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/NonDeO] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]=C=C[CH]OCC=O(29627)'],
    products = ['C#C[CH]CO[CH]C=O(29460)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(3.8e+10,'s^-1'), n=0.87, Ea=(144.348,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_H/Cd;Cs_H_out_H/OneDe] for rate rule [R3H_SS_O;C_rad_out_H/Cd;Cs_H_out_H/CO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH]=C=[C]COCC=O(29628)'],
    products = ['C#C[CH]CO[CH]C=O(29460)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out_H/OneDe] for rate rule [R4H_SSS;Cd_rad_out_double;Cs_H_out_H/CO]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C#C[CH]CO[CH]C=O(29460)'],
    products = ['[CH]=C1[CH]COC1C=O(29603)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(6.8435e+15,'s^-1'), n=-1.17677, Ea=(156.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra;radadd_intra_csHDe] for rate rule [R5_SS_D;doublebond_intra;radadd_intra_csHCO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[C]#C[CH]COCC=O(29629)'],
    products = ['C#C[CH]CO[CH]C=O(29460)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(49528.1,'s^-1'), n=1.95205, Ea=(83.4098,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R6HJ_2;Ct_rad_out;Cs_H_out_H/CO]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['C#C[CH]CO[CH]C=O(29460)'],
    products = ['C#CC=COCC=O(29630)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.08e+10,'s^-1'), n=-0.305, Ea=(93.094,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3;Y_rad_De;XH_Rrad_De] for rate rule [R3radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction43',
    reactants = ['C#C[CH]CO[CH]C=O(29460)'],
    products = ['C#CC1COC1C=O(29465)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Cpri_rad_out_single] + [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/OneDe;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

network(
    label = '4678',
    isomers = [
        'C#C[CH]CO[CH]C=O(29460)',
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
    label = '4678',
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

