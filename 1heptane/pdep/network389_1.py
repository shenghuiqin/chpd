species(
    label = '[O]C(=O)O[CH]C=O(845)',
    structure = SMILES('[O]C=COC([O])=O'),
    E0 = (-355.164,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,322.546,797.152,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.145231,'amu*angstrom^2'), symmetry=1, barrier=(3.33914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145231,'amu*angstrom^2'), symmetry=1, barrier=(3.33914,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.046,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12368,0.0480491,-1.64942e-05,-3.15853e-08,2.04839e-11,-42598.8,22.1204], Tmin=(100,'K'), Tmax=(930.248,'K')), NASAPolynomial(coeffs=[22.127,-0.0032247,3.23384e-06,-6.10314e-10,3.55843e-14,-48195.6,-86.769], Tmin=(930.248,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-355.164,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-O2d)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdOsOs) + radical(OC=OOJ) + radical(C=COJ)"""),
)

species(
    label = 'CO2(13)',
    structure = SMILES('O=C=O'),
    E0 = (-403.131,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([459.166,1086.67,1086.68,2300.05],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0095,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1622.99,'J/mol'), sigma=(3.941,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.2779,0.00275783,7.12787e-06,-1.07855e-08,4.14228e-12,-48475.6,5.97856], Tmin=(100,'K'), Tmax=(988.185,'K')), NASAPolynomial(coeffs=[4.55071,0.00290728,-1.14643e-06,2.25798e-10,-1.69526e-14,-48986,-1.45662], Tmin=(988.185,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-403.131,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), label="""CO2""", comment="""Thermo library: BurkeH2O2"""),
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
    label = '[O][CH]C1OC(=O)O1(6224)',
    structure = SMILES('[O][CH]C1OC(=O)O1'),
    E0 = (-253.386,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.046,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22112,0.0277107,1.31094e-05,-3.59656e-08,1.46549e-11,-30401.5,21.958], Tmin=(100,'K'), Tmax=(1048.71,'K')), NASAPolynomial(coeffs=[12.5214,0.0152869,-7.54427e-06,1.58994e-09,-1.20769e-13,-34039.1,-35.2678], Tmin=(1048.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-253.386,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cds-OdOsOs) + ring(Cyclobutane) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[O]C1([O])OC=CO1(6225)',
    structure = SMILES('[O]C1([O])OC=CO1'),
    E0 = (-194.025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.046,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.50103,0.0194207,3.76585e-05,-6.31307e-08,2.51999e-11,-23269.8,15.0755], Tmin=(100,'K'), Tmax=(980.203,'K')), NASAPolynomial(coeffs=[12.038,0.014044,-5.44211e-06,1.09346e-09,-8.42627e-14,-26750.7,-38.9635], Tmin=(980.203,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-194.025,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cyclopentane) + radical(OCOJ) + radical(OCOJ)"""),
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
    label = '[O]C(=O)OC=C=O(6226)',
    structure = SMILES('[O]C(=O)OC=C=O'),
    E0 = (-315.151,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,3010,987.5,1337.5,450,1655,424.351,425.194,425.277,425.39,425.963,426.106,426.184,426.252],'cm^-1')),
        HinderedRotor(inertia=(0.000930279,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167567,'amu*angstrom^2'), symmetry=1, barrier=(21.506,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.038,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34553,0.052533,-5.9293e-05,2.95347e-08,-5.14681e-12,-37803.2,20.7182], Tmin=(100,'K'), Tmax=(1029.8,'K')), NASAPolynomial(coeffs=[15.9676,0.00434193,-1.63196e-06,3.2048e-10,-2.44249e-14,-41271,-52.4694], Tmin=(1029.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-315.151,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-O2d)H) + group(Cds-(Cdd-O2d)OsH) + group(Cds-OdOsOs) + radical(OC=OOJ)"""),
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
    label = '[O]C(=O)OC=[C]O(6227)',
    structure = SMILES('[O]C(=O)OC=[C]O'),
    E0 = (-256.883,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1277.5,1000,3010,987.5,1337.5,450,1655,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.046,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.162205,0.0656347,-7.73476e-05,4.01666e-08,-7.55636e-12,-30741.5,27.2963], Tmin=(100,'K'), Tmax=(1534.82,'K')), NASAPolynomial(coeffs=[21.4664,-0.00351882,3.55916e-06,-7.62741e-10,5.2989e-14,-35675.6,-79.3805], Tmin=(1534.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-256.883,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-O2d)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdOsOs) + radical(C=CJO) + radical(OC=OOJ)"""),
)

species(
    label = '[O]C(=O)O[C]=CO(6228)',
    structure = SMILES('[O]C(=O)O[C]=CO'),
    E0 = (-256.883,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1277.5,1000,3010,987.5,1337.5,450,1655,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.046,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.162205,0.0656347,-7.73476e-05,4.01666e-08,-7.55636e-12,-30741.5,27.2963], Tmin=(100,'K'), Tmax=(1534.82,'K')), NASAPolynomial(coeffs=[21.4664,-0.00351882,3.55916e-06,-7.62741e-10,5.2989e-14,-35675.6,-79.3805], Tmin=(1534.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-256.883,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-O2d)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdOsOs) + radical(OC=OOJ) + radical(C=CJO)"""),
)

species(
    label = '[O]C=[C]OC(=O)O(6229)',
    structure = SMILES('[O]C=[C]OC(=O)O'),
    E0 = (-359.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.046,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0283,0.0521755,-3.09601e-05,-1.35787e-08,1.30546e-11,-43078.9,25.1151], Tmin=(100,'K'), Tmax=(947.67,'K')), NASAPolynomial(coeffs=[21.2258,-0.000659287,1.3588e-06,-2.19099e-10,8.1404e-15,-48362.7,-78.9205], Tmin=(947.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-359.168,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-O2d)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdOsOs) + radical(C=CJO) + radical(C=COJ)"""),
)

species(
    label = '[O][C]=COC(=O)O(6230)',
    structure = SMILES('O=[C][CH]OC(=O)O'),
    E0 = (-362.992,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.046,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14672,0.0541976,-5.08185e-05,1.81097e-08,-1.18721e-12,-43547.7,25.0556], Tmin=(100,'K'), Tmax=(1056.9,'K')), NASAPolynomial(coeffs=[17.4807,0.0064941,-3.14852e-06,6.77003e-10,-5.26791e-14,-47788.7,-58.38], Tmin=(1056.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-362.992,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-OdOsOs) + radical(CCsJOC(O)) + radical(CsCJ=O)"""),
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
    label = 'CO3t1(671)',
    structure = SMILES('[O]C([O])=O'),
    E0 = (-60.3782,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([180,620.881,889.217,889.555,890.434,1871.84],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.0089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36989,0.0105617,-1.27477e-06,-7.37855e-09,3.82335e-12,-7236.21,10.2409], Tmin=(100,'K'), Tmax=(995.95,'K')), NASAPolynomial(coeffs=[7.17794,0.00288737,-1.19279e-06,2.48518e-10,-1.94626e-14,-8372.65,-10.0125], Tmin=(995.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-60.3782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CO3t1""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[O][C]=O(669)',
    structure = SMILES('[O][C]=O'),
    E0 = (31.9507,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.90478,-0.000175995,8.15126e-06,-1.13656e-08,4.4768e-12,3848.25,8.04855], Tmin=(100,'K'), Tmax=(975.388,'K')), NASAPolynomial(coeffs=[5.59398,-0.00122084,7.11747e-07,-9.7712e-11,3.97995e-15,3238.91,-1.49318], Tmin=(975.388,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.9507,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(OJC=O) + radical((O)CJOH)"""),
)

species(
    label = '[O]C(=O)OC1[CH]O1(6231)',
    structure = SMILES('[O]C(=O)OC1[CH]O1'),
    E0 = (-228.198,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.046,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.650786,0.0524135,-5.24461e-05,2.54278e-08,-4.48534e-12,-27306.9,25.3663], Tmin=(100,'K'), Tmax=(1705.46,'K')), NASAPolynomial(coeffs=[14.045,0.00644531,7.84323e-07,-3.8367e-10,3.17862e-14,-29759,-40.198], Tmin=(1705.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-228.198,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-O2d)) + group(O2s-(Cds-O2d)H) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cds-OdOsOs) + ring(Ethylene_oxide) + radical(CCsJO) + radical(OC=OOJ)"""),
)

species(
    label = '[O]C=CO[C]1OO1(6232)',
    structure = SMILES('[O]C=CO[C]1OO1'),
    E0 = (2.48762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.046,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23,0.0400183,1.80451e-05,-7.79948e-08,4.04747e-11,418.61,22.0097], Tmin=(100,'K'), Tmax=(901.663,'K')), NASAPolynomial(coeffs=[25.3765,-0.0107383,8.71781e-06,-1.77054e-09,1.18078e-13,-6226.95,-104.689], Tmin=(901.663,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.48762,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(Cs-OsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(dioxirane) + radical(Cs_P) + radical(C=COJ)"""),
)

species(
    label = '[O]C1[CH]OC(=O)O1(6220)',
    structure = SMILES('[O]C1[CH]OC(=O)O1'),
    E0 = (-321.524,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.046,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.44604,0.0175122,4.83198e-05,-7.80239e-08,3.11079e-11,-38599.5,21.2631], Tmin=(100,'K'), Tmax=(979.732,'K')), NASAPolynomial(coeffs=[14.4102,0.010203,-4.08451e-06,9.08664e-10,-7.56649e-14,-42937.3,-46.3858], Tmin=(979.732,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-321.524,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cds-OdOsOs) + ring(Cyclopentane) + radical(CCsJOC(O)) + radical(CCOJ)"""),
)

species(
    label = '[O][C]1OC=COO1(6233)',
    structure = SMILES('O=C1O[CH][CH]OO1'),
    E0 = (-201.877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.046,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77726,0.0194871,8.38501e-05,-1.35105e-07,5.47367e-11,-24173.6,18.0667], Tmin=(100,'K'), Tmax=(973.318,'K')), NASAPolynomial(coeffs=[24.3382,0.000470698,-4.25564e-07,4.16255e-10,-5.48353e-14,-32056.4,-108.1], Tmin=(973.318,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-201.877,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-OdOsOs) + ring(Cyclohexanone) + radical(CCsJOOC) + radical(CCsJOC(O))"""),
)

species(
    label = 'O=C=COC(=O)O(6234)',
    structure = SMILES('O=C=COC(=O)O'),
    E0 = (-558.899,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.046,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06472,0.0546174,-4.59084e-05,7.54964e-09,4.13548e-12,-67105.4,21.2497], Tmin=(100,'K'), Tmax=(973.327,'K')), NASAPolynomial(coeffs=[19.0292,0.00307815,-8.29156e-07,1.99569e-10,-1.9407e-14,-71658.3,-70.3564], Tmin=(973.327,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-558.899,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-O2d)H) + group(Cds-(Cdd-O2d)OsH) + group(Cds-OdOsOs)"""),
)

species(
    label = 'O=C1OC=COO1(6235)',
    structure = SMILES('O=C1OC=COO1'),
    E0 = (-469.785,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.046,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.09147,0.00365006,0.000142896,-2.06265e-07,8.28454e-11,-56398.2,14.6947], Tmin=(100,'K'), Tmax=(948.27,'K')), NASAPolynomial(coeffs=[27.6928,-0.00671732,4.87069e-06,-6.62473e-10,1.84492e-14,-65642.9,-130.601], Tmin=(948.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-469.785,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-O2d)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdOsOs) + ring(Cyclohexane)"""),
)

species(
    label = '[O]C(=O)C([O])C=O(844)',
    structure = SMILES('[O]C(=O)C([O])C=O'),
    E0 = (-251.26,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.046,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.50613,0.0286618,-1.42251e-05,2.60564e-09,-1.54191e-13,-30220.2,19.4767], Tmin=(100,'K'), Tmax=(2734.7,'K')), NASAPolynomial(coeffs=[28.4214,-0.0015739,-1.04546e-06,2.22721e-10,-1.2229e-14,-46168.5,-130.034], Tmin=(2734.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-251.26,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-CsCsOsH) + group(Cds-OdCsOs) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCOJ)"""),
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
    label = '[CH]=COC([O])=O(5744)',
    structure = SMILES('[CH]=COC([O])=O'),
    E0 = (-40.7382,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,505.777,505.777,505.777,505.777,505.778,505.778,505.778],'cm^-1')),
        HinderedRotor(inertia=(0.0970506,'amu*angstrom^2'), symmetry=1, barrier=(17.6175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0970506,'amu*angstrom^2'), symmetry=1, barrier=(17.6175,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0462,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67307,0.0428131,-3.21757e-05,9.2322e-10,5.24308e-12,-4808.31,19.282], Tmin=(100,'K'), Tmax=(962.973,'K')), NASAPolynomial(coeffs=[15.9281,0.00309862,-6.84936e-07,1.48425e-10,-1.44561e-14,-8457.79,-53.6464], Tmin=(962.973,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-40.7382,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-O2d)H) + group(Cds-CdsOsH) + group(Cds-OdOsOs) + group(Cds-CdsHH) + radical(Cds_P) + radical(OC=OOJ)"""),
)

species(
    label = '[O]C=CO[C]=O(4793)',
    structure = SMILES('[O]C=CO[C]=O'),
    E0 = (-165.081,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.50427,'amu*angstrom^2'), symmetry=1, barrier=(34.5862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.50379,'amu*angstrom^2'), symmetry=1, barrier=(34.5751,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0462,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4422,0.0438568,-1.98603e-05,-2.20964e-08,1.61314e-11,-19751,18.2045], Tmin=(100,'K'), Tmax=(923.27,'K')), NASAPolynomial(coeffs=[19.5817,-0.00241557,2.81467e-06,-5.59068e-10,3.43845e-14,-24477.9,-75.319], Tmin=(923.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-165.081,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdOsH) + radical((O)CJOC) + radical(C=COJ)"""),
)

species(
    label = '[O]C1([O])OC1C=O(6236)',
    structure = SMILES('[O]C1([O])OC1C=O'),
    E0 = (-112.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.046,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.29669,0.0388241,-4.27971e-05,3.1018e-08,-9.11046e-12,-13511.5,23.4031], Tmin=(100,'K'), Tmax=(994.378,'K')), NASAPolynomial(coeffs=[5.14396,0.0213451,-7.34082e-06,1.15297e-09,-6.99076e-14,-13779.8,11.1807], Tmin=(994.378,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-112.84,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsOsOs) + group(Cds-OdCsH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[O]C(=O)OC[C]=O(6237)',
    structure = SMILES('[O]C(=O)OC[C]=O'),
    E0 = (-313.171,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.046,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81067,0.0423943,-3.5538e-05,1.41659e-08,-2.21887e-12,-37582.2,23.6841], Tmin=(100,'K'), Tmax=(1527.45,'K')), NASAPolynomial(coeffs=[13.4834,0.0118267,-5.52004e-06,1.06441e-09,-7.45479e-14,-41148.1,-37.5744], Tmin=(1527.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-313.171,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-OdOsOs) + radical(CsCJ=O) + radical(OC=OOJ)"""),
)

species(
    label = 'O=CC1OC(=O)O1(6222)',
    structure = SMILES('O=CC1OC(=O)O1'),
    E0 = (-580.571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.046,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.75245,0.00602566,8.46039e-05,-1.16491e-07,4.4614e-11,-69762.4,19.8439], Tmin=(100,'K'), Tmax=(975.766,'K')), NASAPolynomial(coeffs=[15.2357,0.00916883,-3.72586e-06,9.06166e-10,-8.02907e-14,-74784.3,-53.3243], Tmin=(975.766,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-580.571,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-Cs(Cds-O2d)) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + group(Cds-OdOsOs) + ring(Cyclobutane)"""),
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
    E0 = (-355.164,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-233.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-194.025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-70.8647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-355.164,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-94.1036,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-106.537,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-209.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-243.719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (185.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (13.1046,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-173.372,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (2.48762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-297.812,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-201.877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-330.191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-347.633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-41.3644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (202.267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (77.9239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-112.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-171.593,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-179.597,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-347.257,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(=O)O[CH]C=O(845)'],
    products = ['CO2(13)', 'OCHCHO(48)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]C(=O)O[CH]C=O(845)'],
    products = ['[O][CH]C1OC(=O)O1(6224)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5.448e+10,'s^-1','*|/',3), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]C(=O)O[CH]C=O(845)'],
    products = ['[O]C1([O])OC=CO1(6225)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.49749e+08,'s^-1'), n=0.656505, Ea=(161.139,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SMS;multiplebond_intra;radadd_intra] for rate rule [R6_SMS_CO;carbonylbond_intra;radadd_intra_O]
Euclidian distance = 1.73205080757
family: Intra_R_Add_Exocyclic
Ea raised from 157.0 to 161.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[O]C(=O)OC=C=O(6226)'],
    products = ['[O]C(=O)O[CH]C=O(845)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.54607e+06,'m^3/(mol*s)'), n=0.466452, Ea=(32.4943,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CO2(13)', '[O]C=C[O](2537)'],
    products = ['[O]C(=O)O[CH]C=O(845)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.6e+11,'cm^3/(mol*s)'), n=0, Ea=(66.813,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Cd_R;O_rad/OneDe] for rate rule [CO2;O_rad/OneDe]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 65.9 to 66.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]C(=O)OC=[C]O(6227)'],
    products = ['[O]C(=O)O[CH]C=O(845)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C(=O)O[C]=CO(6228)'],
    products = ['[O]C(=O)O[CH]C=O(845)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleNd;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C(=O)O[CH]C=O(845)'],
    products = ['[O]C=[C]OC(=O)O(6229)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(548,'s^-1'), n=3.09, Ea=(145.603,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(=O)O[CH]C=O(845)'],
    products = ['[O][C]=COC(=O)O(6230)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.468e+06,'s^-1','*|/',3), n=1.554, Ea=(111.445,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R5H_SSSD;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CHCHO(47)', 'CO3t1(671)'],
    products = ['[O]C(=O)O[CH]C=O(845)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.91522e+06,'m^3/(mol*s)'), n=0.151, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;O_rad/OneDe] for rate rule [Cd_rad;O_rad/OneDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O][C]=O(669)', '[O]C=C[O](2537)'],
    products = ['[O]C(=O)O[CH]C=O(845)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.27681e+06,'m^3/(mol*s)'), n=0.151, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [O_rad/OneDe;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C(=O)O[CH]C=O(845)'],
    products = ['[O]C(=O)OC1[CH]O1(6231)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(181.793,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri_HNd_O;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]C(=O)O[CH]C=O(845)'],
    products = ['[O]C=CO[C]1OO1(6232)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.55936e+11,'s^-1'), n=0.551275, Ea=(357.652,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra] for rate rule [R3_CO;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic
Ea raised from 356.8 to 357.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]C(=O)O[CH]C=O(845)'],
    products = ['[O]C1[CH]OC(=O)O1(6220)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.33182e+07,'s^-1'), n=1.01661, Ea=(57.3526,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O]C(=O)O[CH]C=O(845)'],
    products = ['[O][C]1OC=COO1(6233)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.78145e+10,'s^-1'), n=0.346137, Ea=(153.288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SMS;multiplebond_intra;radadd_intra] for rate rule [R6_SMS_CO;carbonyl_intra;radadd_intra_O]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 150.1 to 153.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]C(=O)O[CH]C=O(845)'],
    products = ['O=C=COC(=O)O(6234)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]C(=O)O[CH]C=O(845)'],
    products = ['O=C1OC=COO1(6235)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Ypri_rad_out] for rate rule [R6_SSSDS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]C(=O)O[CH]C=O(845)'],
    products = ['[O]C(=O)C([O])C=O(844)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction19',
    reactants = ['O(4)', '[CH]=COC([O])=O(5744)'],
    products = ['[O]C(=O)O[CH]C=O(845)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['O(4)', '[O]C=CO[C]=O(4793)'],
    products = ['[O]C(=O)O[CH]C=O(845)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [CO_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]C(=O)O[CH]C=O(845)'],
    products = ['[O]C1([O])OC1C=O(6236)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.90568e+10,'s^-1'), n=0.237, Ea=(242.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_csHDe] for rate rule [R4_S_CO;carbonylbond_intra;radadd_intra_csHDe]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Exocyclic
Ea raised from 240.3 to 242.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction22',
    reactants = ['OCHCHO(48)', '[O][C]=O(669)'],
    products = ['[O]C(=O)O[CH]C=O(845)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.2635,'m^3/(mol*s)'), n=2.48, Ea=(21.7568,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Od_CO-DeH;YJ] for rate rule [Od_CO-DeH;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O]C(=O)OC[C]=O(6237)'],
    products = ['[O]C(=O)O[CH]C=O(845)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(0.00141351,'s^-1'), n=4.515, Ea=(133.574,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/NonDeO] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]C(=O)O[CH]C=O(845)'],
    products = ['O=CC1OC(=O)O1(6222)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.6e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

network(
    label = '389',
    isomers = [
        '[O]C(=O)O[CH]C=O(845)',
    ],
    reactants = [
        ('CO2(13)', 'OCHCHO(48)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '389',
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

