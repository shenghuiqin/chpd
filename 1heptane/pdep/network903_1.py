species(
    label = '[O]C(O)O[CH]C=O(1317)',
    structure = SMILES('[O]C=COC([O])O'),
    E0 = (-369.572,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,1380,1390,370,380,2900,435,180,180,180,910.925,911.789],'cm^-1')),
        HinderedRotor(inertia=(0.230122,'amu*angstrom^2'), symmetry=1, barrier=(18.5603,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.807226,'amu*angstrom^2'), symmetry=1, barrier=(18.5597,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0314561,'amu*angstrom^2'), symmetry=1, barrier=(18.555,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41963,0.0577655,-5.76161e-05,2.89678e-08,-5.88519e-12,-44357.2,23.5691], Tmin=(100,'K'), Tmax=(1172.93,'K')), NASAPolynomial(coeffs=[12.2793,0.020731,-1.02544e-05,2.04837e-09,-1.47522e-13,-46904.7,-30.5547], Tmin=(1172.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-369.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-OsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(OCOJ)"""),
)

species(
    label = 'HOCHO(59)',
    structure = SMILES('O=CO'),
    E0 = (-389.211,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1051.41,1051.94,1052.55,2787.25,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00861396,'amu*angstrom^2'), symmetry=1, barrier=(6.76824,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4005.91,'J/mol'), sigma=(3.626,'angstroms'), dipoleMoment=(1.7,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.89836,-0.00355878,3.55205e-05,-4.385e-08,1.71078e-11,-46778.6,7.34954], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.61383,0.00644964,-2.29083e-06,3.6716e-10,-2.18737e-14,-45330.3,0.847884], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-389.211,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""HOCHO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = '[O][CH]C1OC(O)O1(7948)',
    structure = SMILES('[O][CH]C1OC(O)O1'),
    E0 = (-222.798,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.70051,0.0353653,-1.85313e-05,3.51106e-09,-2.20109e-13,-26815.7,16.9002], Tmin=(100,'K'), Tmax=(2842.76,'K')), NASAPolynomial(coeffs=[42.2562,-0.0105503,1.29813e-06,-1.07792e-10,7.43537e-15,-52104.6,-215.314], Tmin=(2842.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-222.798,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-OsOsOsH) + ring(Cyclobutane) + radical(CCsJOH) + radical(CCOJ)"""),
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
    label = '[O]C=COC(=O)O(7949)',
    structure = SMILES('[O]C=COC(=O)O'),
    E0 = (-598.913,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.817156,0.0504314,-4.11392e-06,-5.23439e-08,2.92857e-11,-71899.9,22.7444], Tmin=(100,'K'), Tmax=(932.634,'K')), NASAPolynomial(coeffs=[25.3479,-0.00475498,4.18861e-06,-7.66829e-10,4.35371e-14,-78651.1,-105.555], Tmin=(932.634,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-598.913,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-O2d)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdOsOs) + radical(C=COJ)"""),
)

species(
    label = '[O]C(O)OC=C=O(7950)',
    structure = SMILES('[O]C(O)OC=C=O'),
    E0 = (-329.559,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1380,1390,370,380,2900,435,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0963143,'amu*angstrom^2'), symmetry=1, barrier=(2.21446,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.096309,'amu*angstrom^2'), symmetry=1, barrier=(2.21433,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.012028,'amu*angstrom^2'), symmetry=1, barrier=(28.4123,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41931,0.065022,-0.00011112,1.06097e-07,-3.96275e-11,-39551.9,22.9564], Tmin=(100,'K'), Tmax=(806.605,'K')), NASAPolynomial(coeffs=[5.40536,0.029517,-1.58265e-05,3.14674e-09,-2.21482e-13,-39682.9,7.75634], Tmin=(806.605,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-329.559,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cds-(Cdd-O2d)OsH) + radical(OCOJ)"""),
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
    label = '[O]C=COC=O(830)',
    structure = SMILES('[O]C=COC=O'),
    E0 = (-361.529,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.40903,'amu*angstrom^2'), symmetry=1, barrier=(32.3963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41581,'amu*angstrom^2'), symmetry=1, barrier=(32.5523,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.0541,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64366,0.0361947,7.84074e-06,-4.97625e-08,2.51765e-11,-43382.7,18.8534], Tmin=(100,'K'), Tmax=(938.333,'K')), NASAPolynomial(coeffs=[19.2694,0.000634451,1.42126e-06,-2.52749e-10,9.92051e-15,-48432.8,-74.3425], Tmin=(938.333,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-361.529,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdOsH) + radical(C=COJ)"""),
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
    label = '[O]C=CO[C](O)O(7951)',
    structure = SMILES('[O]C=CO[C](O)O'),
    E0 = (-391.723,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0247695,0.0722794,-8.20477e-05,4.23783e-08,-8.09018e-12,-46954.6,28.2746], Tmin=(100,'K'), Tmax=(1457.23,'K')), NASAPolynomial(coeffs=[22.3681,0.00107138,1.27742e-06,-3.2934e-10,2.36733e-14,-52446.7,-84.6409], Tmin=(1457.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-391.723,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-OsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(Cs_P) + radical(C=COJ)"""),
)

species(
    label = '[O]C(O)OC=[C]O(7952)',
    structure = SMILES('[O]C(O)OC=[C]O'),
    E0 = (-271.29,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,340.327,340.332,340.334,4000],'cm^-1')),
        HinderedRotor(inertia=(0.141408,'amu*angstrom^2'), symmetry=1, barrier=(11.6228,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.330212,'amu*angstrom^2'), symmetry=1, barrier=(27.141,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.19754,'amu*angstrom^2'), symmetry=1, barrier=(16.2362,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.330217,'amu*angstrom^2'), symmetry=1, barrier=(27.141,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15347,0.0671193,-8.95593e-05,6.31614e-08,-1.79253e-11,-32530.1,26.2528], Tmin=(100,'K'), Tmax=(857.664,'K')), NASAPolynomial(coeffs=[10.7671,0.0222826,-1.11418e-05,2.20621e-09,-1.57254e-13,-34179.2,-18.6509], Tmin=(857.664,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-271.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-OsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(OCOJ)"""),
)

species(
    label = '[O]C(O)O[C]=CO(7953)',
    structure = SMILES('[O]C(O)O[C]=CO'),
    E0 = (-271.29,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,340.327,340.332,340.334,4000],'cm^-1')),
        HinderedRotor(inertia=(0.141408,'amu*angstrom^2'), symmetry=1, barrier=(11.6228,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.330212,'amu*angstrom^2'), symmetry=1, barrier=(27.141,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.19754,'amu*angstrom^2'), symmetry=1, barrier=(16.2362,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.330217,'amu*angstrom^2'), symmetry=1, barrier=(27.141,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15347,0.0671193,-8.95593e-05,6.31614e-08,-1.79253e-11,-32530.1,26.2528], Tmin=(100,'K'), Tmax=(857.664,'K')), NASAPolynomial(coeffs=[10.7671,0.0222826,-1.11418e-05,2.20621e-09,-1.57254e-13,-34179.2,-18.6509], Tmin=(857.664,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-271.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-OsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(OCOJ) + radical(C=CJO)"""),
)

species(
    label = '[O]C=[C]OC(O)O(7954)',
    structure = SMILES('[O]C=[C]OC(O)O'),
    E0 = (-357.226,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,340.45,340.596,340.628,340.923],'cm^-1')),
        HinderedRotor(inertia=(0.00144959,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182271,'amu*angstrom^2'), symmetry=1, barrier=(15.0058,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.181874,'amu*angstrom^2'), symmetry=1, barrier=(15.0055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182085,'amu*angstrom^2'), symmetry=1, barrier=(15.0074,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.63959,0.0675911,-7.88483e-05,4.3846e-08,-9.31304e-12,-42837.9,28.3252], Tmin=(100,'K'), Tmax=(1167.56,'K')), NASAPolynomial(coeffs=[17.8796,0.00852749,-2.96733e-06,5.18458e-10,-3.56547e-14,-46863.6,-57.5186], Tmin=(1167.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-357.226,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-OsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = '[O][C]=COC(O)O(7955)',
    structure = SMILES('[O][C]=COC(O)O'),
    E0 = (-357.226,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,340.45,340.596,340.628,340.923],'cm^-1')),
        HinderedRotor(inertia=(0.00144959,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182271,'amu*angstrom^2'), symmetry=1, barrier=(15.0058,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.181874,'amu*angstrom^2'), symmetry=1, barrier=(15.0055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182085,'amu*angstrom^2'), symmetry=1, barrier=(15.0074,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.63959,0.0675911,-7.88483e-05,4.3846e-08,-9.31304e-12,-42837.9,28.3252], Tmin=(100,'K'), Tmax=(1167.56,'K')), NASAPolynomial(coeffs=[17.8796,0.00852749,-2.96733e-06,5.18458e-10,-3.56547e-14,-46863.6,-57.5186], Tmin=(1167.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-357.226,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-OsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = '[O][C](O)OC=CO(7956)',
    structure = SMILES('O=C(O)O[CH][CH]O'),
    E0 = (-430.603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.581129,0.0750567,-9.55336e-05,5.88627e-08,-1.40604e-11,-51666.3,26.4534], Tmin=(100,'K'), Tmax=(1030.73,'K')), NASAPolynomial(coeffs=[16.348,0.0138685,-6.48642e-06,1.2669e-09,-9.05038e-14,-54916.5,-50.0897], Tmin=(1030.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-430.603,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(236.962,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-OdOsOs) + radical(CCsJOH) + radical(CCsJOC(O))"""),
)

species(
    label = '[O]C([O])OC=CO(7957)',
    structure = SMILES('[O]C([O])OC=CO'),
    E0 = (-283.637,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,1380,1390,370,380,2900,435,228.369,253.777,263.798,1302.23,1304.15],'cm^-1')),
        HinderedRotor(inertia=(0.149179,'amu*angstrom^2'), symmetry=1, barrier=(7.44887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.726182,'amu*angstrom^2'), symmetry=1, barrier=(32.2301,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.651635,'amu*angstrom^2'), symmetry=1, barrier=(32.2192,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36792,0.0642897,-9.46785e-05,8.60505e-08,-3.27098e-11,-34025,22.1207], Tmin=(100,'K'), Tmax=(728.434,'K')), NASAPolynomial(coeffs=[5.83002,0.0333691,-1.77902e-05,3.58639e-09,-2.56778e-13,-34504.8,3.17624], Tmin=(728.434,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-283.637,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-OsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(OCOJ) + radical(OCOJ)"""),
)

species(
    label = '[O][CH]O(171)',
    structure = SMILES('[O][CH]O'),
    E0 = (5.37818,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,1989.68],'cm^-1')),
        HinderedRotor(inertia=(0.0937094,'amu*angstrom^2'), symmetry=1, barrier=(2.15456,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.40204,0.0209571,-4.96202e-05,5.95355e-08,-2.44634e-11,660.882,10.3084], Tmin=(100,'K'), Tmax=(868.739,'K')), NASAPolynomial(coeffs=[-0.343229,0.0179325,-9.40003e-06,1.81361e-09,-1.23835e-13,2076.48,32.2523], Tmin=(868.739,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.37818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsHH) + radical(OCJO) + radical(OCOJ)"""),
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
    label = '[O]C([O])O(4533)',
    structure = SMILES('[O]C([O])O'),
    E0 = (-145.908,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2358.32,2359.06],'cm^-1')),
        HinderedRotor(inertia=(3.03285e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (62.0248,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.79705,0.0288346,-0.000108603,1.59381e-07,-7.20067e-11,-17564.9,12.6695], Tmin=(100,'K'), Tmax=(855.73,'K')), NASAPolynomial(coeffs=[-15.1242,0.0502973,-2.8813e-05,5.74792e-09,-3.99989e-13,-11874.2,115.335], Tmin=(855.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-145.908,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + radical(OCOJ) + radical(OCOJ)"""),
)

species(
    label = '[O]C(O)OC1[CH]O1(7958)',
    structure = SMILES('[O]C(O)OC1[CH]O1'),
    E0 = (-217.053,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53978,0.0659065,-0.00011602,1.14991e-07,-4.24275e-11,-26028.1,24.5685], Tmin=(100,'K'), Tmax=(881.574,'K')), NASAPolynomial(coeffs=[1.66512,0.0365556,-1.71061e-05,3.15614e-09,-2.10636e-13,-24931.7,30.323], Tmin=(881.574,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-217.053,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-OsOsOsH) + ring(Ethylene_oxide) + radical(CCsJO) + radical(OCOJ)"""),
)

species(
    label = '[O]C1[CH]OC(O)O1(7944)',
    structure = SMILES('[O]C1[CH]OC(O)O1'),
    E0 = (-309.793,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45841,0.053182,-6.89591e-05,4.6233e-08,-1.14209e-11,-37165.5,21.2225], Tmin=(100,'K'), Tmax=(1201.39,'K')), NASAPolynomial(coeffs=[10.4974,0.0107523,-5.83142e-07,-2.55446e-10,3.13433e-14,-38447.2,-20.3394], Tmin=(1201.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-309.793,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-OsOsOsH) + ring(1,3-Dioxolane) + radical(CCOJ) + radical(CCsJOCs)"""),
)

species(
    label = 'O=C(O)OC=CO(7959)',
    structure = SMILES('O=C(O)OC=CO'),
    E0 = (-740.375,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.398235,0.05748,-7.73428e-06,-5.83967e-08,3.40293e-11,-88896.6,22.8384], Tmin=(100,'K'), Tmax=(918.61,'K')), NASAPolynomial(coeffs=[28.6833,-0.00858079,6.89172e-06,-1.34035e-09,8.46863e-14,-96502.5,-124.333], Tmin=(918.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-740.375,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-O2d)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdOsOs)"""),
)

species(
    label = 'O=C=COC(O)O(7960)',
    structure = SMILES('O=C=COC(O)O'),
    E0 = (-556.957,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.939134,0.06689,-8.26498e-05,5.03537e-08,-1.19407e-11,-66875.7,23.519], Tmin=(100,'K'), Tmax=(1036.87,'K')), NASAPolynomial(coeffs=[14.6447,0.0140168,-6.15987e-06,1.17343e-09,-8.27389e-14,-69717.9,-43.0986], Tmin=(1036.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-556.957,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cds-(Cdd-O2d)OsH)"""),
)

species(
    label = 'OC1OC=COO1(7961)',
    structure = SMILES('OC1OC=COO1'),
    E0 = (-365.855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66914,0.0338819,2.69995e-05,-7.35317e-08,3.503e-11,-43902.1,19.0065], Tmin=(100,'K'), Tmax=(913.087,'K')), NASAPolynomial(coeffs=[19.0592,0.00338907,2.03653e-06,-5.05497e-10,3.15999e-14,-48982.4,-73.7385], Tmin=(913.087,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-365.855,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(124trioxene)"""),
)

species(
    label = '[O]C=CO[CH]OO(7962)',
    structure = SMILES('[O]C=CO[CH]OO'),
    E0 = (-146.781,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1310,387.5,850,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.29463,'amu*angstrom^2'), symmetry=1, barrier=(29.766,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29498,'amu*angstrom^2'), symmetry=1, barrier=(29.774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29539,'amu*angstrom^2'), symmetry=1, barrier=(29.7835,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29482,'amu*angstrom^2'), symmetry=1, barrier=(29.7705,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.39354,0.0872848,-0.000101465,5.15637e-08,-9.39557e-12,-17432.4,29.9631], Tmin=(100,'K'), Tmax=(1621.07,'K')), NASAPolynomial(coeffs=[27.1526,-0.00593521,5.87393e-06,-1.24889e-09,8.61358e-14,-23694.1,-112.313], Tmin=(1621.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-146.781,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-OsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(OCJO)"""),
)

species(
    label = '[O]C(O)C([O])C=O(1315)',
    structure = SMILES('[O]C(O)C([O])C=O'),
    E0 = (-263.825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8123,0.0527175,-6.2343e-05,4.46689e-08,-1.37608e-11,-31656.2,29.5754], Tmin=(100,'K'), Tmax=(773.215,'K')), NASAPolynomial(coeffs=[6.66135,0.0276317,-1.36765e-05,2.70743e-09,-1.93253e-13,-32406,7.42883], Tmin=(773.215,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-263.825,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + radical(CCOJ) + radical(C=OCOJ)"""),
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
    label = '[O]C=CO[CH]O(3654)',
    structure = SMILES('[O]C=CO[CH]O'),
    E0 = (-218.286,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.09087,'amu*angstrom^2'), symmetry=1, barrier=(25.0812,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09035,'amu*angstrom^2'), symmetry=1, barrier=(25.0692,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09309,'amu*angstrom^2'), symmetry=1, barrier=(25.1323,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0621,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.703973,0.0538287,-1.33893e-05,-4.98225e-08,3.15716e-11,-26117.4,21.6566], Tmin=(100,'K'), Tmax=(895.327,'K')), NASAPolynomial(coeffs=[26.8921,-0.0113518,8.99636e-06,-1.84754e-09,1.25436e-13,-32883.7,-113.389], Tmin=(895.327,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-218.286,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-OsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(OCJO)"""),
)

species(
    label = '[CH]=COC([O])O(7825)',
    structure = SMILES('[CH]=COC([O])O'),
    E0 = (-55.1457,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1380,1390,370,380,2900,435,180,180,1315.94],'cm^-1')),
        HinderedRotor(inertia=(0.140246,'amu*angstrom^2'), symmetry=1, barrier=(3.22452,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136696,'amu*angstrom^2'), symmetry=1, barrier=(3.14292,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17176,'amu*angstrom^2'), symmetry=1, barrier=(26.941,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0621,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74952,0.0552523,-8.37739e-05,7.72141e-08,-2.92071e-11,-6557.02,21.512], Tmin=(100,'K'), Tmax=(757.148,'K')), NASAPolynomial(coeffs=[5.38531,0.028261,-1.48808e-05,2.97651e-09,-2.11761e-13,-6884.48,6.45624], Tmin=(757.148,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-55.1457,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(OCOJ)"""),
)

species(
    label = '[O]C(O)OC[C]=O(7963)',
    structure = SMILES('[O]C(O)OC[C]=O'),
    E0 = (-302.026,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,3615,1277.5,1000,1380,1390,370,380,2900,435,180,1464.65,1464.89,1465.27],'cm^-1')),
        HinderedRotor(inertia=(0.179587,'amu*angstrom^2'), symmetry=1, barrier=(4.12906,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179457,'amu*angstrom^2'), symmetry=1, barrier=(4.12606,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179626,'amu*angstrom^2'), symmetry=1, barrier=(4.12995,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179626,'amu*angstrom^2'), symmetry=1, barrier=(4.12995,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74129,0.0673101,-0.000140161,1.59619e-07,-6.56519e-11,-36261.1,26.3222], Tmin=(100,'K'), Tmax=(835.571,'K')), NASAPolynomial(coeffs=[-3.72526,0.048856,-2.69256e-05,5.35918e-09,-3.75185e-13,-33789.8,61.0345], Tmin=(835.571,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-302.026,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + radical(CsCJ=O) + radical(OCOJ)"""),
)

species(
    label = '[O][C](O)OCC=O(7964)',
    structure = SMILES('[O][C](O)OCC=O'),
    E0 = (-251.45,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,360,370,350,3615,1277.5,1000,2782.5,750,1395,475,1775,1000,180,1639.44,1639.8,1639.85],'cm^-1')),
        HinderedRotor(inertia=(0.144171,'amu*angstrom^2'), symmetry=1, barrier=(3.31476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144183,'amu*angstrom^2'), symmetry=1, barrier=(3.31506,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144265,'amu*angstrom^2'), symmetry=1, barrier=(3.31693,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144195,'amu*angstrom^2'), symmetry=1, barrier=(3.31532,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72682,0.0676914,-0.000141112,1.60397e-07,-6.58473e-11,-30177.9,26.3834], Tmin=(100,'K'), Tmax=(836.285,'K')), NASAPolynomial(coeffs=[-3.65803,0.0487147,-2.68398e-05,5.33983e-09,-3.73672e-13,-27713,60.7518], Tmin=(836.285,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-251.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-OsOsOsH) + group(Cds-OdCsH) + radical(OCOJ) + radical(Cs_P)"""),
)

species(
    label = '[O]C([O])OCC=O(7965)',
    structure = SMILES('[O]C([O])OCC=O'),
    E0 = (-229.299,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,1512.02,1512.24,1512.62,1512.78,1512.95],'cm^-1')),
        HinderedRotor(inertia=(0.200777,'amu*angstrom^2'), symmetry=1, barrier=(4.61626,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200508,'amu*angstrom^2'), symmetry=1, barrier=(4.61008,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200483,'amu*angstrom^2'), symmetry=1, barrier=(4.6095,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23985,0.0642068,-0.00015568,1.98647e-07,-8.63891e-11,-27539.7,23.6318], Tmin=(100,'K'), Tmax=(839.285,'K')), NASAPolynomial(coeffs=[-12.5246,0.0662166,-3.71019e-05,7.41408e-09,-5.19516e-13,-22653.9,106.618], Tmin=(839.285,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-229.299,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-OsOsOsH) + group(Cds-OdCsH) + radical(OCOJ) + radical(OCOJ)"""),
)

species(
    label = 'OC1O[CH][CH]OO1(7966)',
    structure = SMILES('OC1O[CH][CH]OO1'),
    E0 = (-104.467,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.330749,0.0617989,-6.48273e-05,3.46096e-08,-6.71715e-12,-12416,24.1324], Tmin=(100,'K'), Tmax=(1565.3,'K')), NASAPolynomial(coeffs=[13.4943,0.0105958,1.07164e-06,-6.25764e-10,5.53795e-14,-14385.2,-38.3989], Tmin=(1565.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-104.467,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-OsOsOsH) + ring(124trioxane) + radical(CCsJOCs) + radical(CCsJOOC)"""),
)

species(
    label = 'O=CCOC(=O)O(7967)',
    structure = SMILES('O=CCOC(=O)O'),
    E0 = (-711.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47604,0.0442841,-1.38558e-05,-1.67612e-08,9.79655e-12,-85483.9,23.7217], Tmin=(100,'K'), Tmax=(1031.99,'K')), NASAPolynomial(coeffs=[15.2433,0.0151419,-6.70063e-06,1.35774e-09,-1.01777e-13,-89615.1,-49.3795], Tmin=(1031.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-711.59,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-OdOsOs)"""),
)

species(
    label = 'O=CC1OC(O)O1(1320)',
    structure = SMILES('O=CC1OC(O)O1'),
    E0 = (-549.983,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.24012,0.0347095,-1.13038e-05,-3.81745e-09,1.94216e-12,-66081.1,22.1291], Tmin=(100,'K'), Tmax=(1442.37,'K')), NASAPolynomial(coeffs=[11.395,0.02196,-1.11888e-05,2.2045e-09,-1.54578e-13,-70036.8,-29.9486], Tmin=(1442.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-549.983,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-OsOsOsH) + group(Cds-OdCsH) + ring(Cyclobutane)"""),
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
    E0 = (-369.572,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-222.798,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-340.615,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-85.2722,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-238.432,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-347.937,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-211.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-108.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-120.945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-312.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-312.928,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-281.289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-204.946,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-13.4679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (99.9401,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-187.779,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-309.793,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-344.599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-344.599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-362.041,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-54.0356,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-55.7718,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (24.719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (187.859,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-198.165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-168.452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-92.0403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-127.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-104.467,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-280.603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-361.664,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(O)O[CH]C=O(1317)'],
    products = ['HOCHO(59)', 'OCHCHO(48)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]C(O)O[CH]C=O(1317)'],
    products = ['[O][CH]C1OC(O)O1(7948)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(146.774,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 145.7 to 146.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', '[O]C=COC(=O)O(7949)'],
    products = ['[O]C(O)O[CH]C=O(1317)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.83701,'m^3/(mol*s)'), n=1.71338, Ea=(46.5052,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;HJ] + [CO-NdNd_O;YJ] for rate rule [CO-NdNd_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[O]C(O)OC=C=O(7950)'],
    products = ['[O]C(O)O[CH]C=O(1317)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.54607e+06,'m^3/(mol*s)'), n=0.466452, Ea=(32.4943,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['OH(5)', '[O]C=COC=O(830)'],
    products = ['[O]C(O)O[CH]C=O(1317)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.245e-08,'m^3/(mol*s)'), n=3.486, Ea=(94.7258,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-NdH_O;OJ] for rate rule [CO-NdH_O;OJ_pri]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['HOCHO(59)', '[O]C=C[O](2537)'],
    products = ['[O]C(O)O[CH]C=O(1317)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.6e+11,'cm^3/(mol*s)'), n=0, Ea=(60.12,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO_O;O_rad/OneDe] for rate rule [CO-NdH_O;O_rad/OneDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C(O)O[CH]C=O(1317)'],
    products = ['[O]C=CO[C](O)O(7951)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.43381e+07,'s^-1'), n=1.70481, Ea=(158.036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_NonDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C(O)OC=[C]O(7952)'],
    products = ['[O]C(O)O[CH]C=O(1317)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(O)O[C]=CO(7953)'],
    products = ['[O]C(O)O[CH]C=O(1317)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleNd;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C=[C]OC(O)O(7954)'],
    products = ['[O]C(O)O[CH]C=O(1317)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;XH_out] for rate rule [R4H_SSS_OCs;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O][C]=COC(O)O(7955)'],
    products = ['[O]C(O)O[CH]C=O(1317)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(760143,'s^-1'), n=1.62386, Ea=(44.2978,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;XH_out] for rate rule [R5H_DSSS;Cd_rad_out;O_H_out]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C(O)O[CH]C=O(1317)'],
    products = ['[O][C](O)OC=CO(7956)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(126000,'s^-1'), n=1.85, Ea=(88.2824,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SMSS;Y_rad_out;XH_out] for rate rule [R5H_SMSS;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]C([O])OC=CO(7957)'],
    products = ['[O]C(O)O[CH]C=O(1317)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(293856,'s^-1'), n=1.47286, Ea=(78.6904,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6H;O_rad_out;XH_out] + [R6H_RSSMS;Y_rad_out;XH_out] for rate rule [R6H_RSSMS;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O][CH]O(171)', '[O]C=C[O](2537)'],
    products = ['[O]C(O)O[CH]C=O(1317)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.27681e+06,'m^3/(mol*s)'), n=0.151, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;O_rad/OneDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['CHCHO(47)', '[O]C([O])O(4533)'],
    products = ['[O]C(O)O[CH]C=O(1317)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.75733e+07,'m^3/(mol*s)'), n=-0.180109, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_rad/NonDe;Y_rad] for rate rule [O_rad/NonDe;Cd_rad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -0.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]C(O)O[CH]C=O(1317)'],
    products = ['[O]C(O)OC1[CH]O1(7958)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(181.793,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri_HNd_O;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]C(O)O[CH]C=O(1317)'],
    products = ['[O]C1[CH]OC(O)O1(7944)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.66591e+07,'s^-1'), n=1.01661, Ea=(59.7787,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 57.6 to 59.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]C(O)O[CH]C=O(1317)'],
    products = ['O=C(O)OC=CO(7959)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]C(O)O[CH]C=O(1317)'],
    products = ['O=C=COC(O)O(7960)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]C(O)O[CH]C=O(1317)'],
    products = ['OC1OC=COO1(7961)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Ypri_rad_out] for rate rule [R6_SSSDS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]C=CO[CH]OO(7962)'],
    products = ['[O]C(O)O[CH]C=O(1317)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.263e+10,'s^-1'), n=0, Ea=(92.7451,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnOOH;C_rad_out_1H] for rate rule [ROOH;C_rad_out_H/NonDeO]
Euclidian distance = 1.41421356237
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O]C(O)O[CH]C=O(1317)'],
    products = ['[O]C(O)C([O])C=O(1315)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction23',
    reactants = ['O(4)', '[O]C=CO[CH]O(3654)'],
    products = ['[O]C(O)O[CH]C=O(1317)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/O2;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['O(4)', '[CH]=COC([O])O(7825)'],
    products = ['[O]C(O)O[CH]C=O(1317)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O][CH]O(171)', 'OCHCHO(48)'],
    products = ['[O]C(O)O[CH]C=O(1317)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.2635,'m^3/(mol*s)'), n=2.48, Ea=(21.7568,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Od_CO-DeH;YJ] for rate rule [Od_CO-DeH;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O]C(O)OC[C]=O(7963)'],
    products = ['[O]C(O)O[CH]C=O(1317)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(0.00141351,'s^-1'), n=4.515, Ea=(133.574,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/NonDeO] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O][C](O)OCC=O(7964)'],
    products = ['[O]C(O)O[CH]C=O(1317)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.32e+07,'s^-1'), n=1.69, Ea=(159.41,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_O;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R3H_SS_O;Y_rad_out;Cs_H_out_H/CO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O]C([O])OCC=O(7965)'],
    products = ['[O]C(O)O[CH]C=O(1317)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(726946,'s^-1'), n=1.92229, Ea=(102.145,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R4H_SSS;O_rad_out;Cs_H_out_H/CO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[O]C(O)O[CH]C=O(1317)'],
    products = ['OC1O[CH][CH]OO1(7966)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4.39512e+11,'s^-1'), n=0.277081, Ea=(265.105,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6plus;carbonyl_intra_H;radadd_intra] + [R6_linear;multiplebond_intra;radadd_intra] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic
Ea raised from 263.9 to 265.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction30',
    reactants = ['[O]C(O)O[CH]C=O(1317)'],
    products = ['O=CCOC(=O)O(7967)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[O]C(O)O[CH]C=O(1317)'],
    products = ['O=CC1OC(O)O1(1320)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

network(
    label = '903',
    isomers = [
        '[O]C(O)O[CH]C=O(1317)',
    ],
    reactants = [
        ('HOCHO(59)', 'OCHCHO(48)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '903',
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

