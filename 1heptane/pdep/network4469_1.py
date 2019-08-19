species(
    label = 'C#CC([CH2])C(=C)[O](26501)',
    structure = SMILES('C#CC([CH2])C(=C)[O]'),
    E0 = (325.059,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,2950,3100,1380,975,1025,1650,750,770,3400,2100,350,440,435,1725,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,310.359,311.072],'cm^-1')),
        HinderedRotor(inertia=(0.203268,'amu*angstrom^2'), symmetry=1, barrier=(13.9207,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.202694,'amu*angstrom^2'), symmetry=1, barrier=(13.923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29526,'amu*angstrom^2'), symmetry=1, barrier=(89.006,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.647829,0.0714799,-7.85135e-05,4.17373e-08,-7.43084e-12,39218.5,25.1189], Tmin=(100,'K'), Tmax=(863.516,'K')), NASAPolynomial(coeffs=[14.6554,0.0180653,-5.6555e-06,8.73045e-10,-5.41924e-14,36371.6,-42.8808], Tmin=(863.516,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(325.059,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=C(C)OJ) + radical(Isobutyl)"""),
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
    label = 'C#CC1CC1([CH2])[O](27511)',
    structure = SMILES('C#CC1CC1([CH2])[O]'),
    E0 = (505.97,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.870262,0.058898,-3.3534e-05,-8.15205e-09,1.04708e-11,60976,22.2498], Tmin=(100,'K'), Tmax=(924.484,'K')), NASAPolynomial(coeffs=[16.9972,0.0148699,-3.87541e-06,5.87645e-10,-3.96401e-14,56893.9,-60.2381], Tmin=(924.484,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(505.97,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclopropane) + radical(CJC(C)2O) + radical(CC(C)2OJ)"""),
)

species(
    label = '[CH]=C1CC1C(=C)[O](27559)',
    structure = SMILES('[CH]=C1CC1C(=C)[O]'),
    E0 = (428.989,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.989065,0.057933,-3.95472e-05,4.38819e-09,4.16857e-12,51711.2,22.583], Tmin=(100,'K'), Tmax=(968.104,'K')), NASAPolynomial(coeffs=[15.207,0.0177208,-5.95786e-06,1.03235e-09,-7.14194e-14,48089.9,-50.0351], Tmin=(968.104,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(428.989,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1OC(=C)C1[CH2](27391)',
    structure = SMILES('[CH]=C1OC(=C)C1[CH2]'),
    E0 = (490.336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55253,0.0430631,-2.41741e-06,-2.90529e-08,1.46885e-11,59071.6,23.597], Tmin=(100,'K'), Tmax=(970.853,'K')), NASAPolynomial(coeffs=[13.1729,0.0206809,-7.22622e-06,1.29767e-09,-9.20726e-14,55613.7,-38.3087], Tmin=(970.853,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(490.336,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(Cds_P) + radical(Isobutyl)"""),
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
    label = 'C#CC(=C)C(=C)[O](27560)',
    structure = SMILES('C#CC(=C)C(=C)[O]'),
    E0 = (247.32,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2175,525,750,770,3400,2100,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.34219,'amu*angstrom^2'), symmetry=1, barrier=(30.8595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34875,'amu*angstrom^2'), symmetry=1, barrier=(31.0104,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.246847,0.0673336,-6.77143e-05,3.30092e-08,-6.07169e-12,29893.8,23.7651], Tmin=(100,'K'), Tmax=(1495.74,'K')), NASAPolynomial(coeffs=[19.0492,0.00933342,-1.80905e-06,1.84913e-10,-8.80857e-15,25132.5,-71.6295], Tmin=(1495.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(247.32,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)Ct) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=C(C)OJ)"""),
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
    label = 'C2H(33)',
    structure = SMILES('[C]#C'),
    E0 = (557.301,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (25.0293,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.89868,0.0132988,-2.80733e-05,2.89485e-08,-1.07502e-11,67061.6,6.18548], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.6627,0.00382492,-1.36633e-06,2.13455e-10,-1.23217e-14,67168.4,3.92206], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(557.301,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), label="""C2H""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = 'C#C[C](C)C(=C)[O](27561)',
    structure = SMILES('C#C[C](C)C(=C)[O]'),
    E0 = (272.551,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.777804,0.0646811,-6.04285e-05,2.76299e-08,-4.49833e-12,32901.7,23.0609], Tmin=(100,'K'), Tmax=(1030.44,'K')), NASAPolynomial(coeffs=[14.7834,0.0190693,-6.77741e-06,1.16504e-09,-7.81553e-14,29550.4,-47.1836], Tmin=(1030.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(272.551,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=C(C)OJ) + radical(CCJ(C)CO)"""),
)

species(
    label = 'C#C[C]([CH2])C(=C)O(27562)',
    structure = SMILES('[CH]=C=C([CH2])C(=C)O'),
    E0 = (288.992,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.550383,0.0822665,-9.14748e-05,4.84112e-08,-9.48745e-12,34936.9,25.9359], Tmin=(100,'K'), Tmax=(1448.74,'K')), NASAPolynomial(coeffs=[21.9961,0.00662962,6.9751e-07,-3.813e-10,3.29097e-14,29808.8,-86.3477], Tmin=(1448.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(288.992,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C(O)C([CH2])C#C(27563)',
    structure = SMILES('[CH]=C(O)C([CH2])C#C'),
    E0 = (434.35,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,2175,525,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(0.952197,'amu*angstrom^2'), symmetry=1, barrier=(21.8929,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.954405,'amu*angstrom^2'), symmetry=1, barrier=(21.9436,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.951224,'amu*angstrom^2'), symmetry=1, barrier=(21.8705,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.951369,'amu*angstrom^2'), symmetry=1, barrier=(21.8738,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.207352,0.0789663,-8.75161e-05,4.16385e-08,-5.02731e-12,52381.2,25.6027], Tmin=(100,'K'), Tmax=(867.225,'K')), NASAPolynomial(coeffs=[18.6789,0.011797,-2.52066e-06,2.71892e-10,-1.30063e-14,48499.4,-64.7894], Tmin=(867.225,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(434.35,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])C(C)C#C(27564)',
    structure = SMILES('[CH]=C([O])C(C)C#C'),
    E0 = (367.073,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,750,770,3400,2100,350,440,435,1725,1380,1390,370,380,2900,435,180],'cm^-1')),
        HinderedRotor(inertia=(0.819119,'amu*angstrom^2'), symmetry=1, barrier=(18.8331,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.820063,'amu*angstrom^2'), symmetry=1, barrier=(18.8549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.817634,'amu*angstrom^2'), symmetry=1, barrier=(18.799,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.506351,0.073971,-8.53075e-05,5.06325e-08,-1.17541e-11,44277,24.2439], Tmin=(100,'K'), Tmax=(1059.7,'K')), NASAPolynomial(coeffs=[15.2367,0.01837,-6.60585e-06,1.12137e-09,-7.37719e-14,41155,-47.6757], Tmin=(1059.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(367.073,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_P) + radical(C=C(C)OJ)"""),
)

species(
    label = '[C]#CC(C)C(=C)[O](27565)',
    structure = SMILES('[C]#CC(C)C(=C)[O]'),
    E0 = (457.12,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,2175,525,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,298.229,298.234,298.257],'cm^-1')),
        HinderedRotor(inertia=(0.20144,'amu*angstrom^2'), symmetry=1, barrier=(12.7139,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.263223,'amu*angstrom^2'), symmetry=1, barrier=(16.6133,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.201417,'amu*angstrom^2'), symmetry=1, barrier=(12.7139,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.620011,0.0737876,-9.05111e-05,5.9236e-08,-1.52695e-11,55101.1,24.1754], Tmin=(100,'K'), Tmax=(954.584,'K')), NASAPolynomial(coeffs=[13.2342,0.0209309,-7.45486e-06,1.23147e-09,-7.85944e-14,52692.8,-36.0946], Tmin=(954.584,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(457.12,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(C=C(C)OJ)"""),
)

species(
    label = '[C]#CC([CH2])C(=C)O(27566)',
    structure = SMILES('[C]#CC([CH2])C(=C)O'),
    E0 = (524.398,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,250.47,250.472],'cm^-1')),
        HinderedRotor(inertia=(0.345811,'amu*angstrom^2'), symmetry=1, barrier=(15.3952,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.64352,'amu*angstrom^2'), symmetry=1, barrier=(73.1765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.345901,'amu*angstrom^2'), symmetry=1, barrier=(15.3951,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.345879,'amu*angstrom^2'), symmetry=1, barrier=(15.3953,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0771888,0.0818293,-0.000104219,6.64705e-08,-1.6102e-11,63215.7,26.3978], Tmin=(100,'K'), Tmax=(1119.09,'K')), NASAPolynomial(coeffs=[17.1404,0.013531,-2.87733e-06,2.63037e-10,-7.79849e-15,59854.2,-55.7975], Tmin=(1119.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(524.398,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(Isobutyl)"""),
)

species(
    label = 'C#CC([CH2])[C]1CO1(27567)',
    structure = SMILES('C#CC([CH2])[C]1CO1'),
    E0 = (495.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.695638,0.0664804,-7.51667e-05,4.68783e-08,-1.09356e-11,59680.4,26.2649], Tmin=(100,'K'), Tmax=(1276.26,'K')), NASAPolynomial(coeffs=[11.0019,0.0194763,-2.64219e-06,-3.20084e-11,2.15918e-14,58247.1,-21.2796], Tmin=(1276.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(495.176,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CtCsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Ethylene_oxide) + radical(C2CsJO) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=[C]C1CCC1=O(27385)',
    structure = SMILES('[CH]=[C]C1CCC1=O'),
    E0 = (442.335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70672,0.0427938,-1.37269e-05,-7.46062e-09,4.40941e-12,53289.5,24.0364], Tmin=(100,'K'), Tmax=(1131.95,'K')), NASAPolynomial(coeffs=[10.3156,0.0265823,-1.10744e-05,2.0672e-09,-1.44186e-13,50430.2,-22.5847], Tmin=(1131.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(442.335,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutanone) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'C=C([O])C1[C]=CC1(27552)',
    structure = SMILES('C=C([O])C1[C]=CC1'),
    E0 = (391.901,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3769,0.0481342,-1.71051e-05,-1.44935e-08,9.69487e-12,47237.8,22.3858], Tmin=(100,'K'), Tmax=(980.19,'K')), NASAPolynomial(coeffs=[13.6312,0.0198425,-7.04208e-06,1.26484e-09,-8.91976e-14,43792.3,-41.8102], Tmin=(980.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(391.901,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(cyclobutene-vinyl) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C1[C]=COC1=C(27402)',
    structure = SMILES('[CH2]C1[C]=COC1=C'),
    E0 = (405.081,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9225,0.029736,4.08159e-05,-7.92804e-08,3.4429e-11,48809.6,22.3349], Tmin=(100,'K'), Tmax=(925.74,'K')), NASAPolynomial(coeffs=[14.5353,0.0162086,-3.6518e-06,5.50641e-10,-4.03829e-14,44718.7,-47.0235], Tmin=(925.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(405.081,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = 'C#CC(=C)C(=C)O(27568)',
    structure = SMILES('C#CC(=C)C(=C)O'),
    E0 = (109.515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.464813,0.0645515,-3.7141e-05,-1.32674e-08,1.41409e-11,13311.1,21.4732], Tmin=(100,'K'), Tmax=(929.188,'K')), NASAPolynomial(coeffs=[21.3669,0.00842017,-1.16958e-06,1.28263e-10,-1.11712e-14,7965.45,-85.6949], Tmin=(929.188,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(109.515,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)Ct) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH)"""),
)

species(
    label = 'C#C[CH]CC(=C)[O](26500)',
    structure = SMILES('C#C[CH]CC(=C)[O]'),
    E0 = (299.26,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,750,770,3400,2100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,277.485,277.503,818.472],'cm^-1')),
        HinderedRotor(inertia=(0.361648,'amu*angstrom^2'), symmetry=1, barrier=(19.7611,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00218904,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30596,'amu*angstrom^2'), symmetry=1, barrier=(71.3587,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3806.35,'J/mol'), sigma=(6.32833,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=594.54 K, Pc=34.08 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.571636,0.066372,-6.86045e-05,3.71894e-08,-7.74864e-12,36123.8,26.5803], Tmin=(100,'K'), Tmax=(1311.25,'K')), NASAPolynomial(coeffs=[14.9484,0.0160513,-3.64585e-06,4.03536e-10,-1.83472e-14,32909.2,-44.5561], Tmin=(1311.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(299.26,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=C(C)OJ) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#CC[CH]C(=C)[O](27513)',
    structure = SMILES('C#CCC=C([CH2])[O]'),
    E0 = (304.267,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2175,525,750,770,3400,2100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,423.164,423.178],'cm^-1')),
        HinderedRotor(inertia=(0.0847569,'amu*angstrom^2'), symmetry=1, barrier=(10.7695,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0847552,'amu*angstrom^2'), symmetry=1, barrier=(10.7696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.188257,'amu*angstrom^2'), symmetry=1, barrier=(23.9213,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.688243,0.0650521,-6.51532e-05,3.39492e-08,-6.9157e-12,36720.7,25.9303], Tmin=(100,'K'), Tmax=(1236.69,'K')), NASAPolynomial(coeffs=[15.2323,0.0167811,-5.11388e-06,7.79946e-10,-4.80022e-14,33217.4,-46.946], Tmin=(1236.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(304.267,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=C(O)CJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C#CC1COC1=C(27468)',
    structure = SMILES('C#CC1COC1=C'),
    E0 = (158.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03593,0.0488806,6.93347e-06,-6.08613e-08,3.26885e-11,19126.9,15.2225], Tmin=(100,'K'), Tmax=(893.271,'K')), NASAPolynomial(coeffs=[20.4128,0.00707516,1.63191e-06,-5.55922e-10,4.04733e-14,13871.3,-86.1136], Tmin=(893.271,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(158.015,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(2methyleneoxetane)"""),
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
    label = '[CH]=C=CC(=C)[O](27553)',
    structure = SMILES('[CH]=C=CC(=C)[O]'),
    E0 = (312.889,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,540,610,2055,3120,650,792.5,1650,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.31068,'amu*angstrom^2'), symmetry=1, barrier=(30.1352,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.987821,0.0569827,-6.40071e-05,3.53991e-08,-7.32825e-12,37748.2,21.031], Tmin=(100,'K'), Tmax=(1355.49,'K')), NASAPolynomial(coeffs=[15.0496,0.00757563,-5.78236e-07,-1.02605e-10,1.36376e-14,34662.9,-48.4054], Tmin=(1355.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(312.889,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=C(C)OJ)"""),
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
    label = 'C#CC([CH2])[C]=C(27503)',
    structure = SMILES('C#CC([CH2])[C]=C'),
    E0 = (639.655,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2175,525,2950,3100,1380,975,1025,1650,750,770,3400,2100,3000,3100,440,815,1455,1000,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.564352,'amu*angstrom^2'), symmetry=1, barrier=(12.9756,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2595,0.0612351,-7.30673e-05,4.92288e-08,-1.33647e-11,77030.5,21.837], Tmin=(100,'K'), Tmax=(899.573,'K')), NASAPolynomial(coeffs=[9.94523,0.0226112,-8.6597e-06,1.494e-09,-9.79377e-14,75467.9,-19.1468], Tmin=(899.573,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(639.655,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C1CC(=O)C1[CH2](27361)',
    structure = SMILES('[CH]=C1CC(=O)C1[CH2]'),
    E0 = (421.322,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40917,0.0461271,-1.09873e-05,-1.56193e-08,7.72755e-12,50775.7,22.5945], Tmin=(100,'K'), Tmax=(1126.57,'K')), NASAPolynomial(coeffs=[13.9506,0.0240178,-1.1401e-05,2.29053e-09,-1.66948e-13,46527.2,-45.7196], Tmin=(1126.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(421.322,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(Cds_P) + radical(CJC(C)C=O)"""),
)

species(
    label = 'C#C[C]([CH2])C(C)=O(27569)',
    structure = SMILES('C#CC([CH2])=C(C)[O]'),
    E0 = (268.971,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.670575,0.0654164,-6.3905e-05,3.22189e-08,-6.37296e-12,32476.3,24.059], Tmin=(100,'K'), Tmax=(1238.66,'K')), NASAPolynomial(coeffs=[15.529,0.0174338,-5.79828e-06,9.44685e-10,-6.08089e-14,28795.4,-50.8045], Tmin=(1238.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(268.971,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCtCs) + group(Cds-CdsCsOs) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(Allyl_P) + radical(C=C(C)OJ)"""),
)

species(
    label = '[C]#CC([CH2])C(C)=O(27570)',
    structure = SMILES('[C]#CC([CH2])C(C)=O'),
    E0 = (528.737,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,375,552.5,462.5,1710,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,1762.78],'cm^-1')),
        HinderedRotor(inertia=(0.339127,'amu*angstrom^2'), symmetry=1, barrier=(7.79721,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.339978,'amu*angstrom^2'), symmetry=1, barrier=(7.81676,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.33885,'amu*angstrom^2'), symmetry=1, barrier=(7.79084,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.79376,'amu*angstrom^2'), symmetry=1, barrier=(64.234,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.688509,0.0807269,-0.000127404,1.09611e-07,-3.61131e-11,63704.2,25.7369], Tmin=(100,'K'), Tmax=(906.487,'K')), NASAPolynomial(coeffs=[7.82013,0.0303316,-1.26954e-05,2.21702e-09,-1.42931e-13,63188.8,-3.67983], Tmin=(906.487,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(528.737,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH]=[C]C1COC1=C(27454)',
    structure = SMILES('[CH]=[C]C1COC1=C'),
    E0 = (501.402,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05634,0.0533065,-1.66104e-05,-2.87643e-08,1.92434e-11,60421.5,18.6923], Tmin=(100,'K'), Tmax=(894.715,'K')), NASAPolynomial(coeffs=[17.3987,0.0122945,-1.58474e-06,7.60801e-11,-2.04852e-15,56214.4,-65.5007], Tmin=(894.715,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(501.402,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(Cds_S) + radical(Cds_P)"""),
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
    label = 'C#CC(=C)C(C)=O(27571)',
    structure = SMILES('C#CC(=C)C(C)=O'),
    E0 = (111.599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0334,0.0564496,-3.98538e-05,1.10246e-08,-3.60947e-13,13536.4,21.7571], Tmin=(100,'K'), Tmax=(1176.22,'K')), NASAPolynomial(coeffs=[14.8249,0.0208811,-8.94649e-06,1.698e-09,-1.19614e-13,9508.07,-50.3498], Tmin=(1176.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(111.599,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)HHH) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH)"""),
)

species(
    label = 'C#CC([CH2])C[C]=O(26509)',
    structure = SMILES('C#CC([CH2])C[C]=O'),
    E0 = (375.139,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2175,525,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,216.448],'cm^-1')),
        HinderedRotor(inertia=(0.273079,'amu*angstrom^2'), symmetry=1, barrier=(9.08227,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00207502,'amu*angstrom^2'), symmetry=1, barrier=(9.08177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.14545,'amu*angstrom^2'), symmetry=1, barrier=(71.3339,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.14561,'amu*angstrom^2'), symmetry=1, barrier=(71.3347,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3615.88,'J/mol'), sigma=(6.08218,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=564.79 K, Pc=36.47 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.972278,0.0677876,-8.47319e-05,6.01089e-08,-1.69222e-11,45226.7,27.6291], Tmin=(100,'K'), Tmax=(955.13,'K')), NASAPolynomial(coeffs=[10.0953,0.0244636,-8.6559e-06,1.3992e-09,-8.69669e-14,43717.4,-14.7432], Tmin=(955.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(375.139,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCCJ=O) + radical(Isobutyl)"""),
)

species(
    label = 'C#CC1CCC1=O(26512)',
    structure = SMILES('C#CC1CCC1=O'),
    E0 = (123.347,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83859,0.0370794,9.76956e-06,-3.75634e-08,1.67884e-11,14922.4,20.296], Tmin=(100,'K'), Tmax=(971.264,'K')), NASAPolynomial(coeffs=[11.1401,0.0238287,-8.46275e-06,1.51197e-09,-1.06101e-13,11933.7,-30.3914], Tmin=(971.264,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(123.347,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclobutanone)"""),
)

species(
    label = '[CH]=[C]OC(=C)C=C(27397)',
    structure = SMILES('C#COC([CH2])=C[CH2]'),
    E0 = (411.739,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,750,770,3400,2100,2175,525,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,968.642],'cm^-1')),
        HinderedRotor(inertia=(0.965661,'amu*angstrom^2'), symmetry=1, barrier=(22.2024,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.294189,'amu*angstrom^2'), symmetry=1, barrier=(86.5082,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.966572,'amu*angstrom^2'), symmetry=1, barrier=(22.2234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.964488,'amu*angstrom^2'), symmetry=1, barrier=(22.1755,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.684427,0.0678694,-6.75807e-05,3.44492e-08,-6.93969e-12,49644.4,25.1782], Tmin=(100,'K'), Tmax=(1207.92,'K')), NASAPolynomial(coeffs=[15.2599,0.0196032,-7.64363e-06,1.36927e-09,-9.32382e-14,46123.2,-47.8934], Tmin=(1207.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(411.739,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Ct-CtOs) + group(Ct-CtH) + radical(C=C(O)CJ) + radical(Allyl_P)"""),
)

species(
    label = 'C#CC([CH2])[C]=O(27572)',
    structure = SMILES('C#CC([CH2])[C]=O'),
    E0 = (405.153,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,1855,455,950,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.487557,'amu*angstrom^2'), symmetry=1, barrier=(11.2099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.715122,'amu*angstrom^2'), symmetry=1, barrier=(16.4421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.74859,'amu*angstrom^2'), symmetry=1, barrier=(63.1956,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42283,0.0600401,-8.83487e-05,6.98319e-08,-2.15006e-11,48818.2,21.924], Tmin=(100,'K'), Tmax=(904.058,'K')), NASAPolynomial(coeffs=[9.14647,0.0185846,-7.48393e-06,1.29112e-09,-8.31062e-14,47719.3,-12.9132], Tmin=(904.058,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(405.153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CJC(C)C=O) + radical(CC(C)CJ=O)"""),
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
    E0 = (325.059,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (505.97,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (428.989,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (490.336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (477.103,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (452.061,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (597.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (434.104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (423.526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (482.169,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (626.479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (411.382,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (553.162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (568.696,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (611.769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (556.275,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (450.823,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (455.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (405.081,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (388.459,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (484.994,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (484.994,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (333.343,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (694.452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (882.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (421.322,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (442.118,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (568.636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (501.402,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (383.635,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (388.459,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (620.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (333.343,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (725.539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (786.716,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C#CC([CH2])C(=C)[O](26501)'],
    products = ['CH2CO(28)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C#CC([CH2])C(=C)[O](26501)'],
    products = ['C#CC1CC1([CH2])[O](27511)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.24059e+09,'s^-1'), n=0.736667, Ea=(180.911,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 179.4 to 180.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['C#CC([CH2])C(=C)[O](26501)'],
    products = ['[CH]=C1CC1C(=C)[O](27559)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.881e+08,'s^-1'), n=1.062, Ea=(103.929,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for R4_S_T;triplebond_intra_H;radadd_intra_cs2H
Exact match found for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 102.0 to 103.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['C#CC([CH2])C(=C)[O](26501)'],
    products = ['[CH]=C1OC(=C)C1[CH2](27391)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(165.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_O] for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 161.3 to 165.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C#CC(=C)C(=C)[O](27560)'],
    products = ['C#CC([CH2])C(=C)[O](26501)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.36e+08,'cm^3/(mol*s)'), n=1.64, Ea=(17.9912,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2645 used for Cds-CdCt_Cds-HH;HJ
Exact match found for rate rule [Cds-CdCt_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=[C][O](173)', 'CH2CHCCH(26391)'],
    products = ['C#CC([CH2])C(=C)[O](26501)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.00294841,'m^3/(mol*s)'), n=2.48333, Ea=(17.6885,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CtH_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C2H(33)', 'C=CC(=C)[O](2398)'],
    products = ['C#CC([CH2])C(=C)[O](26501)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00684716,'m^3/(mol*s)'), n=2.49, Ea=(21.946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CdH_Cds-HH;CJ] for rate rule [Cds-CdH_Cds-HH;CtJ_Ct]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH2CO(28)', '[CH]=[C]C=C(4699)'],
    products = ['C#CC([CH2])C(=C)[O](26501)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3525.68,'m^3/(mol*s)'), n=1.07323, Ea=(43.3386,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cdd_Od;CJ] + [Ck_O;YJ] for rate rule [Ck_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C#CC([CH2])C(=C)[O](26501)'],
    products = ['C#C[C](C)C(=C)[O](27561)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(58.4615,'s^-1'), n=3.15787, Ea=(98.4673,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C#CC([CH2])C(=C)[O](26501)'],
    products = ['C#C[C]([CH2])C(=C)O(27562)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.19923e+08,'s^-1'), n=1.46351, Ea=(157.11,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;XH_out] + [R3H_SS_2Cd;Y_rad_out;XH_out] for rate rule [R3H_SS_2Cd;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=C(O)C([CH2])C#C(27563)'],
    products = ['C#CC([CH2])C(=C)[O](26501)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=C([O])C(C)C#C(27564)'],
    products = ['C#CC([CH2])C(=C)[O](26501)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[C]#CC(C)C(=C)[O](27565)'],
    products = ['C#CC([CH2])C(=C)[O](26501)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.39293e+07,'s^-1'), n=1.32074, Ea=(96.0416,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Y_rad_out;Cs_H_out_2H] for rate rule [R4H_TSS;Ct_rad_out;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[C]#CC([CH2])C(=C)O(27566)'],
    products = ['C#CC([CH2])C(=C)[O](26501)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(380071,'s^-1'), n=1.62386, Ea=(44.2978,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;XH_out] for rate rule [R5H_TSSS;Ct_rad_out;O_H_out]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=[C][O](173)', '[CH]=[C]C=C(4699)'],
    products = ['C#CC([CH2])C(=C)[O](26501)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['C#CC([CH2])C(=C)[O](26501)'],
    products = ['C#CC([CH2])[C]1CO1(27567)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C#CC([CH2])C(=C)[O](26501)'],
    products = ['[CH]=[C]C1CCC1=O(27385)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.51071e+08,'s^-1'), n=0.996667, Ea=(125.764,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C#CC([CH2])C(=C)[O](26501)'],
    products = ['C=C([O])C1[C]=CC1(27552)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.27074e+08,'s^-1'), n=0.924088, Ea=(130.68,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C#CC([CH2])C(=C)[O](26501)'],
    products = ['[CH2]C1[C]=COC1=C(27402)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.47e+11,'s^-1'), n=0.15, Ea=(80.0222,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_T;triplebond_intra_H;radadd_intra] for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 74.4 to 80.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction20',
    reactants = ['C#CC([CH2])C(=C)[O](26501)'],
    products = ['C#CC(=C)C(=C)O(27568)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C#CC([CH2])C(=C)[O](26501)'],
    products = ['C#C[CH]CC(=C)[O](26500)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C#CC([CH2])C(=C)[O](26501)'],
    products = ['C#CC[CH]C(=C)[O](27513)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C#CC([CH2])C(=C)[O](26501)'],
    products = ['C#CC1COC1=C(27468)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CH2(19)', '[CH]=C=CC(=C)[O](27553)'],
    products = ['C#CC([CH2])C(=C)[O](26501)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/TwoDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['O(4)', 'C#CC([CH2])[C]=C(27503)'],
    products = ['C#CC([CH2])C(=C)[O](26501)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['C#CC([CH2])C(=C)[O](26501)'],
    products = ['[CH]=C1CC(=O)C1[CH2](27361)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.46159e+06,'s^-1'), n=1.55572, Ea=(96.2631,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 92.7 to 96.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction27',
    reactants = ['C#CC([CH2])C(=C)[O](26501)'],
    products = ['C#C[C]([CH2])C(C)=O(27569)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.30996,'s^-1'), n=3.5644, Ea=(117.059,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[C]#CC([CH2])C(C)=O(27570)'],
    products = ['C#CC([CH2])C(=C)[O](26501)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(263079,'s^-1'), n=1.73643, Ea=(39.8993,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;Cs_H_out_2H] for rate rule [R5H_TSSS;Ct_rad_out;Cs_H_out_2H]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C#CC([CH2])C(=C)[O](26501)'],
    products = ['[CH]=[C]C1COC1=C(27454)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(3.27074e+08,'s^-1'), n=0.924088, Ea=(176.343,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_cs2H]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic
Ea raised from 174.1 to 176.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction30',
    reactants = ['C#CC([CH2])C(=C)[O](26501)'],
    products = ['[CH2]C1[C]=CCC1=O(27337)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(3.47e+11,'s^-1'), n=0.15, Ea=(58.576,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_T;triplebond_intra_H;radadd_intra] for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C#CC([CH2])C(=C)[O](26501)'],
    products = ['C#CC(=C)C(C)=O(27571)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C#CC([CH2])C[C]=O(26509)'],
    products = ['C#CC([CH2])C(=C)[O](26501)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C#CC([CH2])C(=C)[O](26501)'],
    products = ['C#CC1CCC1=O(26512)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]=[C]OC(=C)C=C(27397)'],
    products = ['C#CC([CH2])C(=C)[O](26501)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction35',
    reactants = ['CH2(19)', 'C#CC([CH2])[C]=O(27572)'],
    products = ['C#CC([CH2])C(=C)[O](26501)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4469',
    isomers = [
        'C#CC([CH2])C(=C)[O](26501)',
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
    label = '4469',
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

