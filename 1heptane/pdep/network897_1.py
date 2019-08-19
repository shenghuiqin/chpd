species(
    label = '[O]C(C=O)O[CH]C=O(1311)',
    structure = SMILES('[O]C=COC([O])C=O'),
    E0 = (-290.344,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,266.809,267.46,267.479,267.595,267.737],'cm^-1')),
        HinderedRotor(inertia=(0.43173,'amu*angstrom^2'), symmetry=1, barrier=(21.9157,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.429668,'amu*angstrom^2'), symmetry=1, barrier=(21.9216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.432166,'amu*angstrom^2'), symmetry=1, barrier=(21.9199,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.386483,0.0624346,-2.59843e-05,-3.11712e-08,2.19752e-11,-34774.4,29.7735], Tmin=(100,'K'), Tmax=(928.8,'K')), NASAPolynomial(coeffs=[24.7494,0.00133247,1.92598e-06,-4.08465e-10,2.24894e-14,-41190.2,-96.1391], Tmin=(928.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-290.344,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=OCOJ) + radical(C=COJ)"""),
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
    label = '[O]C=COC1OC1[O](4864)',
    structure = SMILES('[O]C=COC1OC1[O]'),
    E0 = (-238.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.533745,0.0571357,-1.98385e-06,-6.46682e-08,3.76493e-11,-28485.4,25.9274], Tmin=(100,'K'), Tmax=(883.425,'K')), NASAPolynomial(coeffs=[25.5816,-0.00203942,6.40028e-06,-1.49932e-09,1.06516e-13,-35027.5,-103.788], Tmin=(883.425,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-238.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(C=COJ)"""),
)

species(
    label = '[O][CH]C1OC(C=O)O1(4865)',
    structure = SMILES('[O][CH]C1OC(C=O)O1'),
    E0 = (-159.905,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06464,0.0491964,-3.2789e-05,9.95813e-09,-1.22938e-12,-19169.6,25.0987], Tmin=(100,'K'), Tmax=(1738.34,'K')), NASAPolynomial(coeffs=[11.5351,0.0274042,-1.39844e-05,2.74634e-09,-1.922e-13,-22462.1,-25.8272], Tmin=(1738.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-159.905,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(Cyclobutane) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = '[O]C1OC=COC1[O](4866)',
    structure = SMILES('[O]C1OC=COC1[O]'),
    E0 = (-232.307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35337,0.0420782,1.66652e-05,-6.71144e-08,3.41523e-11,-27829.6,18.2542], Tmin=(100,'K'), Tmax=(893.357,'K')), NASAPolynomial(coeffs=[19.0729,0.00689552,1.59755e-06,-5.42267e-10,3.93316e-14,-32757.6,-75.0961], Tmin=(893.357,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-232.307,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(23dihydro14dioxin) + radical(CCOJ) + radical(CCOJ)"""),
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
    label = '[O]C=COC(=O)C=O(4867)',
    structure = SMILES('[O]C=COC(=O)C=O'),
    E0 = (-455.615,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.463307,0.0667479,-5.67732e-05,1.39469e-08,2.18023e-12,-54660.7,25.5391], Tmin=(100,'K'), Tmax=(1002.62,'K')), NASAPolynomial(coeffs=[20.5862,0.00879001,-3.46048e-06,7.04341e-10,-5.47833e-14,-59817.9,-77.1902], Tmin=(1002.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-455.615,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-O2d)H) + radical(C=COJ)"""),
)

species(
    label = '[O]C(C=O)OC=C=O(4868)',
    structure = SMILES('[O]C(C=O)OC=C=O'),
    E0 = (-250.331,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,375.533,375.533,375.533,375.533],'cm^-1')),
        HinderedRotor(inertia=(0.16349,'amu*angstrom^2'), symmetry=1, barrier=(16.3612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16349,'amu*angstrom^2'), symmetry=1, barrier=(16.3612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16349,'amu*angstrom^2'), symmetry=1, barrier=(16.3612,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.384456,0.0695138,-7.76382e-05,4.10978e-08,-8.26336e-12,-29969.1,29.1772], Tmin=(100,'K'), Tmax=(1274.44,'K')), NASAPolynomial(coeffs=[19.5759,0.00727455,-2.0243e-06,3.09742e-10,-2.01157e-14,-34698,-67.426], Tmin=(1274.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-250.331,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + group(Cds-(Cdd-O2d)OsH) + radical(C=OCOJ)"""),
)

species(
    label = 'HCO(16)',
    structure = SMILES('[CH]=O'),
    E0 = (32.4782,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1131.19,1955.83,1955.83],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.23755,-0.00332075,1.4003e-05,-1.3424e-08,4.37416e-12,3872.41,3.30835], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.92002,0.00252279,-6.71004e-07,1.05616e-10,-7.43798e-15,3653.43,3.58077], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(32.4782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HCO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = '[O]C=CO[C](O)C=O(4869)',
    structure = SMILES('[O]C=CO[C](O)C=O'),
    E0 = (-328.831,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.209619,0.0657735,-3.15002e-05,-2.85594e-08,2.17589e-11,-39396.4,30.3159], Tmin=(100,'K'), Tmax=(927.997,'K')), NASAPolynomial(coeffs=[26.1803,-0.000821388,2.84267e-06,-5.72807e-10,3.3342e-14,-46169.2,-103.558], Tmin=(927.997,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-328.831,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(Cs_P) + radical(C=COJ)"""),
)

species(
    label = '[O]C(C=O)OC=[C]O(4870)',
    structure = SMILES('[O]C(C=O)OC=[C]O'),
    E0 = (-192.062,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1380,1390,370,380,2900,435,290.138,290.138,290.138,290.138],'cm^-1')),
        HinderedRotor(inertia=(0.297991,'amu*angstrom^2'), symmetry=1, barrier=(17.8008,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.297991,'amu*angstrom^2'), symmetry=1, barrier=(17.8008,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.297991,'amu*angstrom^2'), symmetry=1, barrier=(17.8008,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.297992,'amu*angstrom^2'), symmetry=1, barrier=(17.8008,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.704325,0.0815266,-9.20269e-05,4.72155e-08,-8.86242e-12,-22911.5,35.4146], Tmin=(100,'K'), Tmax=(1525.6,'K')), NASAPolynomial(coeffs=[24.4405,0.000447608,2.58845e-06,-6.39877e-10,4.64021e-14,-28820.6,-90.7371], Tmin=(1525.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-192.062,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=OCOJ) + radical(C=CJO)"""),
)

species(
    label = '[O]C=COC(O)[C]=O(4871)',
    structure = SMILES('[O]C=COC(O)[C]=O'),
    E0 = (-379.406,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.223216,0.0654006,-3.05659e-05,-2.93408e-08,2.19717e-11,-45479.7,30.2579], Tmin=(100,'K'), Tmax=(928.674,'K')), NASAPolynomial(coeffs=[26.1247,-0.000699925,2.76846e-06,-5.56197e-10,3.20576e-14,-52250.9,-103.34], Tmin=(928.674,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-379.406,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(CsCJ=O) + radical(C=COJ)"""),
)

species(
    label = '[O]C(C=O)O[C]=CO(4872)',
    structure = SMILES('[O]C(C=O)O[C]=CO'),
    E0 = (-192.062,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1380,1390,370,380,2900,435,290.138,290.138,290.138,290.138],'cm^-1')),
        HinderedRotor(inertia=(0.297991,'amu*angstrom^2'), symmetry=1, barrier=(17.8008,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.297991,'amu*angstrom^2'), symmetry=1, barrier=(17.8008,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.297991,'amu*angstrom^2'), symmetry=1, barrier=(17.8008,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.297992,'amu*angstrom^2'), symmetry=1, barrier=(17.8008,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.704325,0.0815266,-9.20269e-05,4.72155e-08,-8.86242e-12,-22911.5,35.4146], Tmin=(100,'K'), Tmax=(1525.6,'K')), NASAPolynomial(coeffs=[24.4405,0.000447608,2.58845e-06,-6.39877e-10,4.64021e-14,-28820.6,-90.7371], Tmin=(1525.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-192.062,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(C=OCOJ)"""),
)

species(
    label = '[O]C=[C]OC(O)C=O(4873)',
    structure = SMILES('[O]C=[C]OC(O)C=O'),
    E0 = (-294.333,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.474833,0.0766213,-8.10575e-05,3.93477e-08,-7.05343e-12,-35220.2,35.2319], Tmin=(100,'K'), Tmax=(1577.28,'K')), NASAPolynomial(coeffs=[23.5847,0.00311217,7.31678e-07,-2.44068e-10,1.78858e-14,-41255.9,-86.8789], Tmin=(1577.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-294.333,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(C=COJ)"""),
)

species(
    label = '[O][C]=COC(O)C=O(4874)',
    structure = SMILES('[O][C]=COC(O)C=O'),
    E0 = (-294.333,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.474833,0.0766213,-8.10575e-05,3.93477e-08,-7.05343e-12,-35220.2,35.2319], Tmin=(100,'K'), Tmax=(1577.28,'K')), NASAPolynomial(coeffs=[23.5847,0.00311217,7.31678e-07,-2.44068e-10,1.78858e-14,-41255.9,-86.8789], Tmin=(1577.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-294.333,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(C=COJ)"""),
)

species(
    label = '[O][C](C=O)OC=CO(4875)',
    structure = SMILES('[O]C=C([O])OC=CO'),
    E0 = (-321.983,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.67885,0.0909479,-0.000103618,5.22504e-08,-9.4499e-12,-38491.8,33.7705], Tmin=(100,'K'), Tmax=(1640.85,'K')), NASAPolynomial(coeffs=[27.2823,-0.00413632,5.68624e-06,-1.253e-09,8.7396e-14,-44700,-110.249], Tmin=(1640.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-321.983,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C([C]=O)OC=CO(4876)',
    structure = SMILES('[O]C([C]=O)OC=CO'),
    E0 = (-277.136,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,1855,455,950,3615,1277.5,1000,1380,1390,370,380,2900,435,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.95913,'amu*angstrom^2'), symmetry=1, barrier=(22.0523,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.959421,'amu*angstrom^2'), symmetry=1, barrier=(22.059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.959455,'amu*angstrom^2'), symmetry=1, barrier=(22.0598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.960066,'amu*angstrom^2'), symmetry=1, barrier=(22.0738,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.47208,0.0870438,-9.77794e-05,4.86486e-08,-8.69253e-12,-33105.9,35.7392], Tmin=(100,'K'), Tmax=(1658.77,'K')), NASAPolynomial(coeffs=[26.4642,-0.00350636,5.06848e-06,-1.11234e-09,7.70423e-14,-39184.4,-103.561], Tmin=(1658.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-277.136,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=OCOJ) + radical(CsCJ=O)"""),
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
    label = '[O]C([O])C=O(3622)',
    structure = SMILES('[O]C([O])C=O'),
    E0 = (-50.3449,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,1931.96,4000],'cm^-1')),
        HinderedRotor(inertia=(0.430052,'amu*angstrom^2'), symmetry=1, barrier=(9.88775,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0355,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.79889,0.0291308,-4.02657e-05,3.53069e-08,-1.27174e-11,-6014.46,19.0664], Tmin=(100,'K'), Tmax=(803.679,'K')), NASAPolynomial(coeffs=[4.59224,0.0158863,-7.48538e-06,1.42862e-09,-9.89574e-14,-6163.24,11.6742], Tmin=(803.679,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-50.3449,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(C=OCOJ)"""),
)

species(
    label = '[O]C(C=O)OC1[CH]O1(4877)',
    structure = SMILES('[O]C(C=O)OC1[CH]O1'),
    E0 = (-137.825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.890951,0.0658726,-6.69534e-05,3.04097e-08,-3.1276e-12,-16462,29.4038], Tmin=(100,'K'), Tmax=(842.473,'K')), NASAPolynomial(coeffs=[13.9715,0.0174498,-5.09937e-06,7.41008e-10,-4.41224e-14,-19151.6,-34.3418], Tmin=(842.473,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-137.825,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + ring(Ethylene_oxide) + radical(C=OCOJ) + radical(CCsJO)"""),
)

species(
    label = '[O]C=COC1[CH]OO1(4878)',
    structure = SMILES('[O]C=COC1[CH]OO1'),
    E0 = (-24.7609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.380181,0.0629971,-2.64179e-05,-2.81297e-08,1.98611e-11,-2832.48,23.5895], Tmin=(100,'K'), Tmax=(941.922,'K')), NASAPolynomial(coeffs=[23.7741,0.00509243,-1.99818e-07,2.2932e-11,-8.24899e-15,-9077.9,-97.6316], Tmin=(941.922,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-24.7609,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(12dioxetane) + radical(C=COJ) + radical(CCsJOO)"""),
)

species(
    label = '[O]C1[CH]OC(C=O)O1(4879)',
    structure = SMILES('[O]C1[CH]OC(C=O)O1'),
    E0 = (-246.9,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.526907,0.0602904,-6.54641e-05,3.564e-08,-7.05863e-12,-29556.3,26.763], Tmin=(100,'K'), Tmax=(1516.87,'K')), NASAPolynomial(coeffs=[14.173,0.00847466,1.43021e-06,-6.40502e-10,5.4836e-14,-31874.9,-38.7537], Tmin=(1516.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-246.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(1,3-Dioxolane) + radical(CCOJ) + radical(CCsJOCs)"""),
)

species(
    label = '[O]C1[CH]OOC=CO1(4845)',
    structure = SMILES('[O]C1[CH]OOC=CO1'),
    E0 = (32.2129,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0738,0.00764031,0.00019839,-3.15174e-07,1.36856e-10,4033.02,25.569], Tmin=(100,'K'), Tmax=(903.019,'K')), NASAPolynomial(coeffs=[45.1698,-0.0383337,2.66675e-05,-5.24166e-09,3.44985e-13,-10020.3,-216.388], Tmin=(903.019,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(32.2129,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(CCOJ) + radical(CCsJOOC)"""),
)

species(
    label = 'OCHCO(61)',
    structure = SMILES('[O]C=C=O'),
    E0 = (-72.3947,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,3010,987.5,1337.5,450,1655,180],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0281,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3660.03,'J/mol'), sigma=(4.01,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39406,0.0143629,-8.84138e-06,1.60961e-09,3.20389e-13,-8684.14,10.592], Tmin=(298,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.94994,0.010163,-5.5772e-06,1.4572e-09,-1.47435e-13,-9096.51,2.57988], Tmin=(1000,'K'), Tmax=(2900,'K'))], Tmin=(298,'K'), Tmax=(2900,'K'), E0=(-72.3947,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""OCHCO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'O=C[CH]O(684)',
    structure = SMILES('O=C[CH]O'),
    E0 = (-194.023,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(1.63727,'amu*angstrom^2'), symmetry=1, barrier=(37.644,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.64043,'amu*angstrom^2'), symmetry=1, barrier=(37.7167,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (59.044,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.73076,0.0175105,1.90197e-05,-4.1579e-08,1.80494e-11,-23280.5,13.329], Tmin=(100,'K'), Tmax=(963.671,'K')), NASAPolynomial(coeffs=[12.0193,0.00477202,-1.33661e-06,3.02986e-10,-2.77179e-14,-26269.4,-37.3588], Tmin=(963.671,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-194.023,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(OCJC=O)"""),
)

species(
    label = 'O=CC(=O)OC=CO(4880)',
    structure = SMILES('O=CC(=O)OC=CO'),
    E0 = (-597.078,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.695261,0.0823001,-8.91946e-05,4.40654e-08,-8.06216e-12,-71624.8,28.3022], Tmin=(100,'K'), Tmax=(1523.8,'K')), NASAPolynomial(coeffs=[25.3437,0.00245447,7.17581e-07,-2.21191e-10,1.57361e-14,-78226.2,-103.91], Tmin=(1523.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-597.078,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-O2d)H)"""),
)

species(
    label = 'O=C=COC(O)C=O(4881)',
    structure = SMILES('O=C=COC(O)C=O'),
    E0 = (-494.064,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0185268,0.0742758,-7.99333e-05,4.04265e-08,-7.7068e-12,-59265.4,29.848], Tmin=(100,'K'), Tmax=(1401.08,'K')), NASAPolynomial(coeffs=[21.7131,0.00657696,-1.39913e-06,1.76616e-10,-1.07181e-14,-64799.8,-80.3423], Tmin=(1401.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-494.064,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + group(Cds-(Cdd-O2d)OsH)"""),
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
    label = '[O]CO[CH]C=O(971)',
    structure = SMILES('[O]C=COC[O]'),
    E0 = (-179.848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.40766,'amu*angstrom^2'), symmetry=1, barrier=(32.365,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.40681,'amu*angstrom^2'), symmetry=1, barrier=(32.3454,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0621,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4132.15,'J/mol'), sigma=(6.57465,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=645.43 K, Pc=32.99 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28928,0.0348407,4.20423e-05,-1.036e-07,4.89469e-11,-21510,20.2394], Tmin=(100,'K'), Tmax=(916.929,'K')), NASAPolynomial(coeffs=[26.1607,-0.00887106,7.56513e-06,-1.47529e-09,9.27251e-14,-28794.6,-112.445], Tmin=(916.929,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-179.848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-OsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(OCOJ)"""),
)

species(
    label = 'O=CC1OC=COO1(4852)',
    structure = SMILES('O=CC1OC=COO1'),
    E0 = (-302.963,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3226,0.0341215,5.443e-05,-1.15135e-07,5.26334e-11,-36318.5,22.446], Tmin=(100,'K'), Tmax=(912.946,'K')), NASAPolynomial(coeffs=[24.0606,-0.000597569,4.83265e-06,-1.04303e-09,6.5833e-14,-43175,-99.9945], Tmin=(912.946,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-302.963,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + ring(124trioxene)"""),
)

species(
    label = '[O]C(C=O)C([O])C=O(1310)',
    structure = SMILES('[O]C(C=O)C([O])C=O'),
    E0 = (-162.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21395,0.0534636,-4.18756e-05,1.44531e-08,-1.45425e-12,-19401.3,34.1683], Tmin=(100,'K'), Tmax=(1160.74,'K')), NASAPolynomial(coeffs=[14.2824,0.0175632,-7.28639e-06,1.36663e-09,-9.58967e-14,-23050.4,-33.4779], Tmin=(1160.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-162.2,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(C=OCOJ)"""),
)

species(
    label = 'O=C[CH]OO[CH]C=O(1312)',
    structure = SMILES('O=C[CH]OO[CH]C=O'),
    E0 = (-78.8979,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,3000,3050,390,425,1340,1360,335,370,325.144],'cm^-1')),
        HinderedRotor(inertia=(0.00160997,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.647473,'amu*angstrom^2'), symmetry=1, barrier=(48.0889,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.650852,'amu*angstrom^2'), symmetry=1, barrier=(48.0893,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.646056,'amu*angstrom^2'), symmetry=1, barrier=(48.114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.65351,'amu*angstrom^2'), symmetry=1, barrier=(48.0675,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3937.92,'J/mol'), sigma=(6.28783,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=615.09 K, Pc=35.94 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04095,0.071564,-8.42897e-05,5.59919e-08,-1.57016e-11,-9388.28,27.2427], Tmin=(100,'K'), Tmax=(849.401,'K')), NASAPolynomial(coeffs=[9.08222,0.0336952,-1.74139e-05,3.5021e-09,-2.52265e-13,-10754.3,-10.2391], Tmin=(849.401,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-78.8979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(OCJC=O) + radical(OCJC=O)"""),
)

species(
    label = '[O]C=COC([O])=CO(4882)',
    structure = SMILES('[O]C=COC(=O)[CH]O'),
    E0 = (-424.168,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.407602,0.0793464,-8.50069e-05,4.23977e-08,-7.90879e-12,-50841.9,30.261], Tmin=(100,'K'), Tmax=(1468.71,'K')), NASAPolynomial(coeffs=[23.4666,0.00552077,-6.16334e-07,1.00623e-11,1.08124e-15,-56905.1,-90.8625], Tmin=(1468.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-424.168,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(OCJC=O) + radical(C=COJ)"""),
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
    label = '[O]C=CO[CH]C=O(4757)',
    structure = SMILES('[O]C=COC=C[O]'),
    E0 = (-164.846,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.39708,'amu*angstrom^2'), symmetry=1, barrier=(32.1216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.39822,'amu*angstrom^2'), symmetry=1, barrier=(32.1478,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.75379,0.0473164,2.13002e-05,-8.91986e-08,4.5705e-11,-19687,24.5388], Tmin=(100,'K'), Tmax=(907.428,'K')), NASAPolynomial(coeffs=[28.1154,-0.00907459,8.35811e-06,-1.69861e-09,1.11423e-13,-27296.7,-119.376], Tmin=(907.428,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-164.846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[CH]=COC([O])C=O(4829)',
    structure = SMILES('[CH]=COC([O])C=O'),
    E0 = (24.0822,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,344.037,344.52,344.633],'cm^-1')),
        HinderedRotor(inertia=(0.209835,'amu*angstrom^2'), symmetry=1, barrier=(17.6776,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.210651,'amu*angstrom^2'), symmetry=1, barrier=(17.676,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.210336,'amu*angstrom^2'), symmetry=1, barrier=(17.6786,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.933784,0.0572271,-4.17854e-05,1.52604e-09,6.63666e-12,3016.1,26.9423], Tmin=(100,'K'), Tmax=(954.432,'K')), NASAPolynomial(coeffs=[18.5456,0.0076633,-1.99679e-06,3.51162e-10,-2.76209e-14,-1450.12,-62.9885], Tmin=(954.432,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(24.0822,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(Cds_P)"""),
)

species(
    label = '[O]C1OC(C=O)C1[O](4883)',
    structure = SMILES('[O]C1OC(C=O)C1[O]'),
    E0 = (-129.589,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71058,0.0395893,5.95998e-06,-4.02855e-08,2.01204e-11,-15493.5,26.8773], Tmin=(100,'K'), Tmax=(919.065,'K')), NASAPolynomial(coeffs=[13.0128,0.0183569,-5.01642e-06,7.75164e-10,-5.21025e-14,-18751.8,-33.1193], Tmin=(919.065,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-129.589,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + ring(Oxetane) + radical(CC(C)OJ) + radical(CCOJ)"""),
)

species(
    label = '[O]C(C=O)OC[C]=O(4884)',
    structure = SMILES('[O]C(C=O)OC[C]=O'),
    E0 = (-222.798,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,381.644,381.694,381.919,1862.17],'cm^-1')),
        HinderedRotor(inertia=(0.00115764,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0623524,'amu*angstrom^2'), symmetry=1, barrier=(6.4447,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0622649,'amu*angstrom^2'), symmetry=1, barrier=(6.44467,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.187986,'amu*angstrom^2'), symmetry=1, barrier=(19.4351,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35192,0.0638923,-7.76146e-05,5.4931e-08,-1.64658e-11,-26706.1,30.247], Tmin=(100,'K'), Tmax=(797.49,'K')), NASAPolynomial(coeffs=[7.92185,0.0309395,-1.5634e-05,3.1182e-09,-2.23446e-13,-27754,0.0375645], Tmin=(797.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-222.798,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CsCJ=O)"""),
)

species(
    label = '[O][C](C=O)OCC=O(4885)',
    structure = SMILES('[O]C=C([O])OCC=O'),
    E0 = (-335.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.232116,0.0692795,-4.94619e-05,-2.24594e-09,1.02332e-11,-40155.2,26.8834], Tmin=(100,'K'), Tmax=(949.36,'K')), NASAPolynomial(coeffs=[23.1562,0.0056085,-8.69287e-07,1.5272e-10,-1.58644e-14,-45991.2,-90.333], Tmin=(949.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-335.1,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C([C]=O)OCC=O(4886)',
    structure = SMILES('[O]C([C]=O)OCC=O'),
    E0 = (-222.798,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,381.644,381.694,381.919,1862.17],'cm^-1')),
        HinderedRotor(inertia=(0.00115764,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0623524,'amu*angstrom^2'), symmetry=1, barrier=(6.4447,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0622649,'amu*angstrom^2'), symmetry=1, barrier=(6.44467,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.187986,'amu*angstrom^2'), symmetry=1, barrier=(19.4351,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35192,0.0638923,-7.76146e-05,5.4931e-08,-1.64658e-11,-26706.1,30.247], Tmin=(100,'K'), Tmax=(797.49,'K')), NASAPolynomial(coeffs=[7.92185,0.0309395,-1.5634e-05,3.1182e-09,-2.23446e-13,-27754,0.0375645], Tmin=(797.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-222.798,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CsCJ=O) + radical(C=OCOJ)"""),
)

species(
    label = 'O=CC1O[CH][CH]OO1(4842)',
    structure = SMILES('O=CC1O[CH][CH]OO1'),
    E0 = (-41.5739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.958087,0.0725943,-7.2052e-05,3.52809e-08,-6.19296e-12,-4789.9,30.9941], Tmin=(100,'K'), Tmax=(1768.25,'K')), NASAPolynomial(coeffs=[15.0212,0.0108842,1.98154e-06,-8.07013e-10,6.51576e-14,-6444.56,-43.9041], Tmin=(1768.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-41.5739,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(124trioxane) + radical(CCsJOOC) + radical(CCsJOCs)"""),
)

species(
    label = 'O=CCOC(=O)C=O(4887)',
    structure = SMILES('O=CCOC(=O)C=O'),
    E0 = (-568.293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80935,0.0526577,-3.92067e-05,1.43034e-08,-2.18928e-12,-68275.4,24.037], Tmin=(100,'K'), Tmax=(1439.39,'K')), NASAPolynomial(coeffs=[10.3325,0.0289722,-1.45237e-05,2.8712e-09,-2.03662e-13,-70729,-20.1865], Tmin=(1439.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-568.293,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-OdCsH) + group(Cds-O2d(Cds-O2d)H)"""),
)

species(
    label = 'O=CC1OC(C=O)O1(1314)',
    structure = SMILES('O=CC1OC(C=O)O1'),
    E0 = (-487.09,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78978,0.0365152,9.05398e-06,-3.39797e-08,1.37034e-11,-58493.6,24.5316], Tmin=(100,'K'), Tmax=(1082.74,'K')), NASAPolynomial(coeffs=[12.9683,0.0231495,-1.11251e-05,2.27081e-09,-1.67923e-13,-62551.5,-37.8478], Tmin=(1082.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-487.09,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + group(Cds-OdCsH) + ring(Cyclobutane)"""),
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
    E0 = (-290.344,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-188.054,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-159.905,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-168.456,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-211.397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-6.04432,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-235.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-184.026,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-136.436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-29.2832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-215.063,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-41.7171,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-144.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-178.899,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-202.062,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-192.355,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-37.6922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (195.503,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-108.551,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-24.5616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-232.991,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (32.2129,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-84.7143,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-265.371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-265.371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (6.23947,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-282.813,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (23.456,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (234.902,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-280.301,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (78.1585,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (267.087,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (-129.589,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (-89.2239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (-95.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (-120.653,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (-41.5739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (-201.375,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (-282.436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(C=O)O[CH]C=O(1311)'],
    products = ['OCHCHO(48)', 'OCHCHO(48)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]C(C=O)O[CH]C=O(1311)'],
    products = ['[O]C=COC1OC1[O](4864)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(7.785e+11,'s^-1'), n=0.342, Ea=(102.29,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]C(C=O)O[CH]C=O(1311)'],
    products = ['[O][CH]C1OC(C=O)O1(4865)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(130.439,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 129.1 to 130.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]C(C=O)O[CH]C=O(1311)'],
    products = ['[O]C1OC=COC1[O](4866)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(9.291e+11,'s^-1'), n=0.234, Ea=(121.888,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_SMSR;multiplebond_intra;radadd_intra_O] for rate rule [R7_SMSS_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.82842712475
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', '[O]C=COC(=O)C=O(4867)'],
    products = ['[O]C(C=O)O[CH]C=O(1311)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(25.1243,'m^3/(mol*s)'), n=1.86, Ea=(32.426,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO-DeNd_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[O]C(C=O)OC=C=O(4868)'],
    products = ['[O]C(C=O)O[CH]C=O(1311)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.54607e+06,'m^3/(mol*s)'), n=0.466452, Ea=(32.4943,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['HCO(16)', '[O]C=COC=O(830)'],
    products = ['[O]C(C=O)O[CH]C=O(1311)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(5.2e+11,'cm^3/(mol*s)'), n=0, Ea=(93.9308,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO_O;CO_pri_rad] for rate rule [CO-NdH_O;CO_pri_rad]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['OCHCHO(48)', '[O]C=C[O](2537)'],
    products = ['[O]C(C=O)O[CH]C=O(1311)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(5.2e+11,'cm^3/(mol*s)'), n=0, Ea=(60.12,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO_O;O_rad/OneDe] for rate rule [CO-DeH_O;O_rad/OneDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(C=O)O[CH]C=O(1311)'],
    products = ['[O]C=CO[C](O)C=O(4869)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.70223e+09,'s^-1'), n=1.15155, Ea=(153.908,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_OneDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_CO]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C(C=O)OC=[C]O(4870)'],
    products = ['[O]C(C=O)O[CH]C=O(1311)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C(C=O)O[CH]C=O(1311)'],
    products = ['[O]C=COC(O)[C]=O(4871)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;O_rad_out;XH_out] for rate rule [R3H_SS_Cs;O_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C(C=O)O[C]=CO(4872)'],
    products = ['[O]C(C=O)O[CH]C=O(1311)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleNd;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]C(C=O)O[CH]C=O(1311)'],
    products = ['[O]C=[C]OC(O)C=O(4873)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(274,'s^-1'), n=3.09, Ea=(145.603,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]C(C=O)O[CH]C=O(1311)'],
    products = ['[O][C]=COC(O)C=O(4874)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.234e+06,'s^-1','*|/',3), n=1.554, Ea=(111.445,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R5H_SSSD;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O]C(C=O)O[CH]C=O(1311)'],
    products = ['[O][C](C=O)OC=CO(4875)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(126000,'s^-1'), n=1.85, Ea=(88.2824,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SMSS;Y_rad_out;XH_out] for rate rule [R5H_SMSS;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]C([C]=O)OC=CO(4876)'],
    products = ['[O]C(C=O)O[CH]C=O(1311)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(16140,'s^-1'), n=1.92259, Ea=(84.7802,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_RSSMS;Y_rad_out;XH_out] for rate rule [R6H_RSSMS;CO_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]C=C[O](2537)', '[O]C=C[O](2537)'],
    products = ['[O]C(C=O)O[CH]C=O(1311)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.31073e+07,'m^3/(mol*s)'), n=0.151, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [O_rad/OneDe;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 8.0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['CHCHO(47)', '[O]C([O])C=O(3622)'],
    products = ['[O]C(C=O)O[CH]C=O(1311)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.75733e+07,'m^3/(mol*s)'), n=-0.180109, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;O_rad/NonDe] for rate rule [Cd_rad;O_rad/NonDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -0.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]C(C=O)O[CH]C=O(1311)'],
    products = ['[O]C(C=O)OC1[CH]O1(4877)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(181.793,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri_HNd_O;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]C(C=O)O[CH]C=O(1311)'],
    products = ['[O]C=COC1[CH]OO1(4878)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(265.782,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]C(C=O)O[CH]C=O(1311)'],
    products = ['[O]C1[CH]OC(C=O)O1(4879)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.66591e+07,'s^-1'), n=1.01661, Ea=(57.3526,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O]C(C=O)O[CH]C=O(1311)'],
    products = ['[O]C1[CH]OOC=CO1(4845)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(5.649e+12,'s^-1'), n=0.287, Ea=(322.557,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;carbonyl_intra_H;radadd_intra] for rate rule [R7_linear;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 316.4 to 322.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction23',
    reactants = ['OCHCO(61)', 'O=C[CH]O(684)'],
    products = ['[O]C(C=O)O[CH]C=O(1311)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.52781e-06,'m^3/(mol*s)'), n=3.38614, Ea=(181.703,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [doublebond;H_OR] for rate rule [Cdd_Cd;H_OCdpri]
Euclidian distance = 2.2360679775
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]C(C=O)O[CH]C=O(1311)'],
    products = ['O=CC(=O)OC=CO(4880)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]C(C=O)O[CH]C=O(1311)'],
    products = ['O=C=COC(O)C=O(4881)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CO(12)', '[O]CO[CH]C=O(971)'],
    products = ['[O]C(C=O)O[CH]C=O(1311)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(0.118397,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O]C(C=O)O[CH]C=O(1311)'],
    products = ['O=CC1OC=COO1(4852)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Ypri_rad_out] for rate rule [R6_SSSDS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O]C(C=O)O[CH]C=O(1311)'],
    products = ['[O]C(C=O)C([O])C=O(1310)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction23',
    reactants = ['O=C[CH]OO[CH]C=O(1312)'],
    products = ['[O]C(C=O)O[CH]C=O(1311)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(14080,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[O]C=COC([O])=CO(4882)'],
    products = ['[O]C(C=O)O[CH]C=O(1311)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction31',
    reactants = ['O(4)', '[O]C=CO[CH]C=O(4757)'],
    products = ['[O]C(C=O)O[CH]C=O(1311)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(4171.1,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/OneDeO;O_birad]
Euclidian distance = 4.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction32',
    reactants = ['O(4)', '[CH]=COC([O])C=O(4829)'],
    products = ['[O]C(C=O)O[CH]C=O(1311)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction33',
    reactants = ['[O]C(C=O)O[CH]C=O(1311)'],
    products = ['[O]C1OC(C=O)C1[O](4883)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(902977,'s^-1'), n=1.63829, Ea=(160.755,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs] for rate rule [R5_SS_CO;carbonylbond_intra_H;radadd_intra_csHDe]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Exocyclic
Ea raised from 157.3 to 160.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction34',
    reactants = ['[O]C(C=O)OC[C]=O(4884)'],
    products = ['[O]C(C=O)O[CH]C=O(1311)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(0.00141351,'s^-1'), n=4.515, Ea=(133.574,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/NonDeO] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[O]C(C=O)O[CH]C=O(1311)'],
    products = ['[O][C](C=O)OCC=O(4885)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(6.08189e+06,'s^-1'), n=1.81713, Ea=(195.293,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_H/OneDe;XH_out] for rate rule [R3H_SS_O;C_rad_out_H/OneDe;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[O]C([C]=O)OCC=O(4886)'],
    products = ['[O]C(C=O)O[CH]C=O(1311)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(363473,'s^-1'), n=1.92229, Ea=(102.145,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R4H_SSS;CO_rad_out;Cs_H_out_H/CO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[O]C(C=O)O[CH]C=O(1311)'],
    products = ['O=CC1O[CH][CH]OO1(4842)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(4.39512e+11,'s^-1'), n=0.277081, Ea=(248.77,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6plus;carbonyl_intra_H;radadd_intra] + [R6_linear;multiplebond_intra;radadd_intra] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic
Ea raised from 247.3 to 248.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction38',
    reactants = ['[O]C(C=O)O[CH]C=O(1311)'],
    products = ['O=CCOC(=O)C=O(4887)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[O]C(C=O)O[CH]C=O(1311)'],
    products = ['O=CC1OC(C=O)O1(1314)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

network(
    label = '897',
    isomers = [
        '[O]C(C=O)O[CH]C=O(1311)',
    ],
    reactants = [
        ('OCHCHO(48)', 'OCHCHO(48)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '897',
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

