species(
    label = 'C=C([O])C([O])C=O(1100)',
    structure = SMILES('C=C([O])C([O])C=O'),
    E0 = (-103.358,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,494.237,494.455,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0765363,'amu*angstrom^2'), symmetry=1, barrier=(13.2874,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0767246,'amu*angstrom^2'), symmetry=1, barrier=(13.2632,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25999,0.0529281,-5.02374e-05,2.40558e-08,-4.50418e-12,-12326.3,28.4331], Tmin=(100,'K'), Tmax=(1306.51,'K')), NASAPolynomial(coeffs=[14.2309,0.0132161,-4.64403e-06,7.90982e-10,-5.24414e-14,-15715.6,-37.6121], Tmin=(1306.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-103.358,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=OCOJ)"""),
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
    label = '[CH2]C1([O])OC1C=O(4794)',
    structure = SMILES('[CH2]C1([O])OC1C=O'),
    E0 = (16.1413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39341,0.0457112,-6.38261e-06,-3.82768e-08,2.30274e-11,2046.27,23.1785], Tmin=(100,'K'), Tmax=(880.962,'K')), NASAPolynomial(coeffs=[17.01,0.00825638,4.32537e-07,-3.30937e-10,2.72126e-14,-2003.35,-57.5507], Tmin=(880.962,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(16.1413,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Ethylene_oxide) + radical(CJC(O)2C) + radical(CC(C)(O)OJ)"""),
)

species(
    label = 'C=C([O])C1OC1[O](4806)',
    structure = SMILES('C=C([O])C1OC1[O]'),
    E0 = (-51.0445,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.815195,0.0547522,-5.20278e-05,2.57066e-08,-4.7586e-12,-6011.42,26.7039], Tmin=(100,'K'), Tmax=(1578.62,'K')), NASAPolynomial(coeffs=[13.1673,0.0125184,-1.50676e-06,-1.69628e-11,1.00611e-14,-8548.7,-34.2115], Tmin=(1578.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-51.0445,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C1OC([O])C1[O](4807)',
    structure = SMILES('C=C1OC([O])C1[O]'),
    E0 = (27.7908,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6242,0.0362764,1.9567e-05,-6.31528e-08,3.04094e-11,3442.83,19.8457], Tmin=(100,'K'), Tmax=(920.774,'K')), NASAPolynomial(coeffs=[17.9863,0.00670179,1.31466e-07,-1.26102e-10,5.23889e-15,-1329.78,-67.2956], Tmin=(920.774,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(27.7908,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(CC(C)OJ) + radical(CCOJ)"""),
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
    label = 'C=C([O])C(=O)C=O(4808)',
    structure = SMILES('[CH2]C(=O)C(=O)C=O'),
    E0 = (-168.806,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,365,385,505,600,445,480,1700,1720,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.157078,'amu*angstrom^2'), symmetry=1, barrier=(3.61154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156426,'amu*angstrom^2'), symmetry=1, barrier=(3.59654,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15734,'amu*angstrom^2'), symmetry=1, barrier=(3.61757,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28127,0.0658354,-0.000106884,9.30079e-08,-3.16634e-11,-20210.6,24.6996], Tmin=(100,'K'), Tmax=(835.428,'K')), NASAPolynomial(coeffs=[8.22061,0.0225387,-1.10621e-05,2.11265e-09,-1.45037e-13,-21018.6,-5.42747], Tmin=(835.428,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-168.806,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)(Cds-O2d)) + group(Cds-O2d(Cds-O2d)H) + radical(CJCC=O)"""),
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
    label = 'C=C([O])C=O(734)',
    structure = SMILES('[CH2]C(=O)C=O'),
    E0 = (-74.309,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,375,552.5,462.5,1710,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.00741427,'amu*angstrom^2'), symmetry=1, barrier=(9.92444,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34864,'amu*angstrom^2'), symmetry=1, barrier=(31.0078,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.689,0.0301409,-2.3287e-05,9.30131e-09,-1.56936e-12,-8891.4,16.525], Tmin=(100,'K'), Tmax=(1335.27,'K')), NASAPolynomial(coeffs=[7.34755,0.0161853,-7.60949e-06,1.4738e-09,-1.03802e-13,-10135.5,-7.29667], Tmin=(1335.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-74.309,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(CJCC=O)"""),
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
    label = 'C=C([O])[C](O)C=O(4809)',
    structure = SMILES('C=C([O])C(O)=C[O]'),
    E0 = (-263.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.841086,0.0810722,-9.60159e-05,5.07307e-08,-9.62422e-12,-31482.3,26.6298], Tmin=(100,'K'), Tmax=(1558.09,'K')), NASAPolynomial(coeffs=[23.9807,-0.00345792,5.39331e-06,-1.23024e-09,8.8217e-14,-36691.7,-96.0241], Tmin=(1558.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-263.39,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])C(O)[C]=O(4810)',
    structure = SMILES('C=C([O])C(O)[C]=O'),
    E0 = (-187.131,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.899101,0.0594137,-6.37615e-05,3.34435e-08,-6.69548e-12,-22387.5,29.7639], Tmin=(100,'K'), Tmax=(1303.53,'K')), NASAPolynomial(coeffs=[16.4356,0.0088695,-2.29789e-06,3.20682e-10,-1.91494e-14,-26194.2,-48.3743], Tmin=(1303.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-187.131,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(CCCJ=O)"""),
)

species(
    label = 'C=C(O)[C]([O])C=O(4811)',
    structure = SMILES('C=C(O)C([O])=C[O]'),
    E0 = (-263.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.841086,0.0810722,-9.60159e-05,5.07307e-08,-9.62422e-12,-31482.3,26.6298], Tmin=(100,'K'), Tmax=(1558.09,'K')), NASAPolynomial(coeffs=[23.9807,-0.00345792,5.39331e-06,-1.23024e-09,8.8217e-14,-36691.7,-96.0241], Tmin=(1558.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-263.39,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = '[CH]=C(O)C([O])C=O(4812)',
    structure = SMILES('[CH]=C(O)C([O])C=O'),
    E0 = (5.93318,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2782.5,750,1395,475,1775,1000,3120,650,792.5,1650,3615,1277.5,1000,1380,1390,370,380,2900,435,223.04],'cm^-1')),
        HinderedRotor(inertia=(0.00349196,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.384888,'amu*angstrom^2'), symmetry=1, barrier=(13.709,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.379975,'amu*angstrom^2'), symmetry=1, barrier=(13.7092,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.635021,0.0625981,-6.69586e-05,3.41761e-08,-6.59145e-12,844.484,29.5787], Tmin=(100,'K'), Tmax=(1387.56,'K')), NASAPolynomial(coeffs=[18.285,0.00685428,-1.44008e-06,1.71127e-10,-9.56646e-15,-3585.42,-59.6663], Tmin=(1387.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.93318,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(Cds_P)"""),
)

species(
    label = 'C=C(O)C([O])[C]=O(4813)',
    structure = SMILES('C=C(O)C([O])[C]=O'),
    E0 = (-81.2024,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,1855,455,950,3615,1277.5,1000,1380,1390,370,380,2900,435,269.847,270.045],'cm^-1')),
        HinderedRotor(inertia=(0.00231226,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224255,'amu*angstrom^2'), symmetry=1, barrier=(11.6036,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224415,'amu*angstrom^2'), symmetry=1, barrier=(11.6035,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.562669,0.0629253,-6.85499e-05,3.5396e-08,-6.8365e-12,-9631.79,30.2525], Tmin=(100,'K'), Tmax=(1426.16,'K')), NASAPolynomial(coeffs=[18.4643,0.00541821,-3.89763e-07,-5.38074e-11,6.66608e-15,-13995.8,-59.8656], Tmin=(1426.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-81.2024,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]=C([O])C(O)C=O(4814)',
    structure = SMILES('[CH]=C([O])C(O)C=O'),
    E0 = (-99.9951,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2782.5,750,1395,475,1775,1000,3120,650,792.5,1650,3615,1277.5,1000,1380,1390,370,380,2900,435,346.326],'cm^-1')),
        HinderedRotor(inertia=(0.00142746,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133622,'amu*angstrom^2'), symmetry=1, barrier=(11.1977,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131722,'amu*angstrom^2'), symmetry=1, barrier=(11.1852,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.980298,0.058993,-6.18947e-05,3.19394e-08,-6.35835e-12,-11911.6,29.0577], Tmin=(100,'K'), Tmax=(1239.8,'K')), NASAPolynomial(coeffs=[15.906,0.0108364,-3.62957e-06,6.08038e-10,-4.03178e-14,-15612.5,-46.1578], Tmin=(1239.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-99.9951,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=C(C)OJ)"""),
)

species(
    label = '[O]C(C=O)[C]1CO1(4815)',
    structure = SMILES('[O]C(C=O)[C]1CO1'),
    E0 = (43.9461,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31432,0.0521978,-5.27768e-05,2.94181e-08,-6.29821e-12,5388.11,27.3734], Tmin=(100,'K'), Tmax=(1325.27,'K')), NASAPolynomial(coeffs=[10.9087,0.0158398,-3.2497e-06,2.90629e-10,-8.79387e-15,3494.9,-19.1643], Tmin=(1325.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(43.9461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(Ethylene_oxide) + radical(C2CsJO) + radical(C=OCOJ)"""),
)

species(
    label = '[O][CH]C1OCC1=O(4754)',
    structure = SMILES('[O][CH]C1OCC1=O'),
    E0 = (21.6095,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78431,0.0397954,-1.16349e-05,-1.14922e-08,6.25145e-12,2686.32,24.3797], Tmin=(100,'K'), Tmax=(1102.27,'K')), NASAPolynomial(coeffs=[12.2647,0.019817,-9.01513e-06,1.78185e-09,-1.29135e-13,-720.874,-32.1774], Tmin=(1102.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(21.6095,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsCs) + ring(Cyclobutane) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = 'C=C([O])C1[CH]OO1(4816)',
    structure = SMILES('C=C([O])C1[CH]OO1'),
    E0 = (162.225,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37305,0.0520714,-4.56136e-05,2.03273e-08,-3.60694e-12,19610.4,21.8212], Tmin=(100,'K'), Tmax=(1355.29,'K')), NASAPolynomial(coeffs=[13.2605,0.0169868,-6.783e-06,1.22661e-09,-8.3596e-14,16388.2,-39.1431], Tmin=(1355.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(162.225,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(12dioxetane) + radical(CCsJOO) + radical(C=C(C)OJ)"""),
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
    label = 'C=C(O)C(=O)C=O(4817)',
    structure = SMILES('C=C(O)C(=O)C=O'),
    E0 = (-389.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.857558,0.0624127,-6.55473e-05,3.36872e-08,-6.70351e-12,-46681.1,22.3225], Tmin=(100,'K'), Tmax=(1235.87,'K')), NASAPolynomial(coeffs=[16.3629,0.0122283,-4.63723e-06,8.30368e-10,-5.69939e-14,-50513.6,-55.7654], Tmin=(1235.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-389.117,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsHH) + group(Cds-O2d(Cds-O2d)H)"""),
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
    label = 'C=C([O])C[O](939)',
    structure = SMILES('C=C([O])C[O]'),
    E0 = (1.93929,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,350,440,435,1725,180,1871.52,2779.93],'cm^-1')),
        HinderedRotor(inertia=(0.567805,'amu*angstrom^2'), symmetry=1, barrier=(13.055,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.0627,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3931.94,'J/mol'), sigma=(6.33729,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=614.16 K, Pc=35.05 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.88359,0.0288232,-9.37057e-06,-2.85602e-08,3.0654e-11,269.125,17.4204], Tmin=(100,'K'), Tmax=(490.774,'K')), NASAPolynomial(coeffs=[4.65106,0.0216573,-9.59619e-06,1.80409e-09,-1.24992e-13,8.45312,9.2632], Tmin=(490.774,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1.93929,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCOJ) + radical(C=C(C)OJ)"""),
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
    label = 'C=C([O])O[CH]C=O(1101)',
    structure = SMILES('[CH2]C(=O)OC=C[O]'),
    E0 = (-205.81,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3995.21,'J/mol'), sigma=(6.24473,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=624.04 K, Pc=37.23 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.9051,0.0546582,-2.57195e-05,-1.9774e-08,1.52408e-11,-24629.5,25.1425], Tmin=(100,'K'), Tmax=(942.089,'K')), NASAPolynomial(coeffs=[20.3297,0.00555565,-6.73399e-07,1.03189e-10,-1.2058e-14,-29770.4,-75.2712], Tmin=(942.089,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-205.81,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CJCO)"""),
)

species(
    label = 'C=C([O])C([O])=CO(4818)',
    structure = SMILES('C=C([O])C([O])=CO'),
    E0 = (-267.048,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,3615,1277.5,1000,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.973985,'amu*angstrom^2'), symmetry=1, barrier=(22.3938,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.972455,'amu*angstrom^2'), symmetry=1, barrier=(22.3587,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.734618,0.0824598,-0.00010213,5.65128e-08,-1.12235e-11,-31929.2,26.3265], Tmin=(100,'K'), Tmax=(1485.25,'K')), NASAPolynomial(coeffs=[23.6218,-0.00383948,5.93549e-06,-1.37892e-09,1.00726e-13,-36880.7,-93.1262], Tmin=(1485.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-267.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=C(C)OJ)"""),
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
    label = 'C=C([O])[CH]C=O(4597)',
    structure = SMILES('C=C([O])C=C[O]'),
    E0 = (-49.4969,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.20625,'amu*angstrom^2'), symmetry=1, barrier=(27.7341,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3607,0.0446863,-1.00193e-05,-3.58099e-08,2.21363e-11,-5845.51,19.4834], Tmin=(100,'K'), Tmax=(902.98,'K')), NASAPolynomial(coeffs=[19.5497,0.000521001,2.86674e-06,-6.71638e-10,4.54752e-14,-10614.7,-74.6309], Tmin=(902.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-49.4969,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=[C]C([O])C=O(4819)',
    structure = SMILES('C=[C]C([O])C=O'),
    E0 = (211.237,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,898.917,4000],'cm^-1')),
        HinderedRotor(inertia=(0.625925,'amu*angstrom^2'), symmetry=1, barrier=(14.3913,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.823473,'amu*angstrom^2'), symmetry=1, barrier=(18.9333,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23637,0.0384269,-3.00062e-05,1.23045e-08,-2.10242e-12,25469.6,23.8387], Tmin=(100,'K'), Tmax=(1346.51,'K')), NASAPolynomial(coeffs=[8.79016,0.018958,-8.31824e-06,1.56675e-09,-1.08808e-13,23704.6,-9.72941], Tmin=(1346.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.237,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(Cds_S)"""),
)

species(
    label = '[O]C1CC(=O)C1[O](4740)',
    structure = SMILES('[O]C1CC(=O)C1[O]'),
    E0 = (15.8124,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22039,0.0232512,4.24352e-05,-7.42937e-08,3.05082e-11,1980.23,24.2191], Tmin=(100,'K'), Tmax=(963.56,'K')), NASAPolynomial(coeffs=[14.0316,0.0138338,-4.57322e-06,8.97645e-10,-7.07839e-14,-2134.93,-41.8675], Tmin=(963.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(15.8124,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(C=OCOJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC(=O)[C]([O])C=O(4820)',
    structure = SMILES('CC(=O)C([O])=C[O]'),
    E0 = (-156.648,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.969292,0.0564853,-4.02387e-05,4.17309e-10,6.92788e-12,-18722,23.0044], Tmin=(100,'K'), Tmax=(954.032,'K')), NASAPolynomial(coeffs=[18.2999,0.00802911,-2.11142e-06,3.69636e-10,-2.87939e-14,-23130.4,-65.5637], Tmin=(954.032,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-156.648,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)HHH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsOsH) + radical(C=OCOJ) + radical(C=COJ)"""),
)

species(
    label = 'CC(=O)C([O])[C]=O(4821)',
    structure = SMILES('CC(=O)C([O])[C]=O'),
    E0 = (-106.123,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36653,0.0631787,-9.46894e-05,8.05432e-08,-2.72998e-11,-12673.9,25.9653], Tmin=(100,'K'), Tmax=(829.736,'K')), NASAPolynomial(coeffs=[7.49196,0.0255247,-1.19308e-05,2.24841e-09,-1.53885e-13,-13410.7,-0.757812], Tmin=(829.736,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-106.123,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(C=OCCJ=O)"""),
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
    label = '[O]C1[CH]OCC1=O(4764)',
    structure = SMILES('[O]C1[CH]OCC1=O'),
    E0 = (-41.9711,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98534,0.0267379,4.04564e-05,-7.96494e-08,3.45482e-11,-4959.31,23.6588], Tmin=(100,'K'), Tmax=(941.124,'K')), NASAPolynomial(coeffs=[16.751,0.00848885,-1.39683e-06,2.49614e-10,-2.49119e-14,-9709.65,-57.1528], Tmin=(941.124,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-41.9711,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsCs) + ring(Cyclopentane) + radical(CCsJOCs) + radical(C=OCOJ)"""),
)

species(
    label = 'CC(=O)C(=O)C=O(4822)',
    structure = SMILES('CC(=O)C(=O)C=O'),
    E0 = (-380.395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48505,0.0600203,-8.48964e-05,7.09812e-08,-2.43604e-11,-45664.9,22.6429], Tmin=(100,'K'), Tmax=(778.567,'K')), NASAPolynomial(coeffs=[7.15287,0.0267656,-1.286e-05,2.47599e-09,-1.72506e-13,-46422.1,-2.47751], Tmin=(778.567,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-380.395,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)(Cds-O2d)) + group(Cds-O2d(Cds-O2d)H)"""),
)

species(
    label = '[O]C(C=O)C[C]=O(1104)',
    structure = SMILES('[O]C(C=O)C[C]=O'),
    E0 = (-76.0907,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,769.311,4000],'cm^-1')),
        HinderedRotor(inertia=(0.469049,'amu*angstrom^2'), symmetry=1, barrier=(10.7844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.433925,'amu*angstrom^2'), symmetry=1, barrier=(10.7843,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0256679,'amu*angstrom^2'), symmetry=1, barrier=(10.7835,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4220.17,'J/mol'), sigma=(6.6134,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=659.18 K, Pc=33.11 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80208,0.0510391,-5.38874e-05,3.2043e-08,-7.9556e-12,-9074.71,27.9793], Tmin=(100,'K'), Tmax=(959.602,'K')), NASAPolynomial(coeffs=[8.42412,0.0234354,-1.0738e-05,2.0652e-09,-1.45531e-13,-10345.6,-3.695], Tmin=(959.602,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-76.0907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(C=OCOJ)"""),
)

species(
    label = 'O=CC1OCC1=O(1107)',
    structure = SMILES('O=CC1OCC1=O'),
    E0 = (-306.374,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23491,0.0246692,3.43812e-05,-6.0915e-08,2.40891e-11,-36772.2,20.9774], Tmin=(100,'K'), Tmax=(1001.6,'K')), NASAPolynomial(coeffs=[12.5707,0.0180335,-7.55972e-06,1.531e-09,-1.16145e-13,-40580.3,-37.5774], Tmin=(1001.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-306.374,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + ring(Cyclobutane)"""),
)

species(
    label = '[CH2]C(=C[O])OC=O(4823)',
    structure = SMILES('[CH2]C(=C[O])OC=O'),
    E0 = (-244.404,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,351.153,352.686,353.311],'cm^-1')),
        HinderedRotor(inertia=(0.280018,'amu*angstrom^2'), symmetry=1, barrier=(24.655,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.276201,'amu*angstrom^2'), symmetry=1, barrier=(24.6587,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.280075,'amu*angstrom^2'), symmetry=1, barrier=(24.7456,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.968416,0.0556218,-3.51796e-05,-4.66442e-09,8.42684e-12,-29275.9,23.9707], Tmin=(100,'K'), Tmax=(965.842,'K')), NASAPolynomial(coeffs=[18.3927,0.00912987,-2.84267e-06,5.33322e-10,-4.13917e-14,-33839.1,-65.6843], Tmin=(965.842,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-244.404,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-OdOsH) + radical(C=COJ) + radical(C=C(O)CJ)"""),
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
    label = '[O]C([C]=O)C=O(4824)',
    structure = SMILES('[O]C([C]=O)C=O'),
    E0 = (-51.2525,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.157246,'amu*angstrom^2'), symmetry=1, barrier=(3.61539,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157445,'amu*angstrom^2'), symmetry=1, barrier=(3.61997,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0462,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0486,0.0475792,-7.70531e-05,6.75818e-08,-2.30528e-11,-6098.45,21.925], Tmin=(100,'K'), Tmax=(853.02,'K')), NASAPolynomial(coeffs=[6.58301,0.0173887,-8.26551e-06,1.55251e-09,-1.05383e-13,-6547.23,2.67377], Tmin=(853.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-51.2525,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(C=OCCJ=O) + radical(C=OCOJ)"""),
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
    E0 = (-103.358,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (16.1413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-1.06774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (27.7908,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (64.5547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-26.1707,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (52.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-36.3258,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (22.1019,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-28.0775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (53.7522,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (198.062,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (7.81222,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-55.6866,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (141.339,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (127.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (22.8774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (162.424,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (223.757,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-39.958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (188.027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-32.4598,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (107.99,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-103.358,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (193.508,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (454.242,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (15.8124,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (13.7006,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-27.3336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (294.725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-41.9711,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-39.958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (169.51,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (-95.0739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (69.3957,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (330.311,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C([O])C([O])C=O(1100)'],
    products = ['CH2CO(28)', 'OCHCHO(48)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C([O])C([O])C=O(1100)'],
    products = ['[CH2]C1([O])OC1C=O(4794)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.62453e+10,'s^-1'), n=0.5217, Ea=(119.5,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H;radadd_intra] for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 119.5 to 119.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=C([O])C([O])C=O(1100)'],
    products = ['C=C([O])C1OC1[O](4806)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(7.785e+11,'s^-1'), n=0.342, Ea=(102.29,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=C([O])C([O])C=O(1100)'],
    products = ['C=C1OC([O])C1[O](4807)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(131.149,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_O] for rate rule [R5_SS_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 129.7 to 131.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C=C([O])C(=O)C=O(4808)'],
    products = ['C=C([O])C([O])C=O(1100)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(92.1383,'m^3/(mol*s)'), n=1.68375, Ea=(21.5685,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;HJ] for rate rule [CO-DeDe_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=[C][O](173)', 'OCHCHO(48)'],
    products = ['C=C([O])C([O])C=O(1100)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(210.531,'m^3/(mol*s)'), n=1.35673, Ea=(38.9443,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;CJ] + [CO-DeH_O;YJ] for rate rule [CO-DeH_O;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['HCO(16)', 'C=C([O])C=O(734)'],
    products = ['C=C([O])C([O])C=O(1100)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(5.2e+11,'cm^3/(mol*s)'), n=0, Ea=(93.9308,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO_O;CO_pri_rad] for rate rule [CO-CdH_O;CO_pri_rad]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH2CO(28)', '[O]C=C[O](2537)'],
    products = ['C=C([O])C([O])C=O(1100)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7051.36,'m^3/(mol*s)'), n=1.07323, Ea=(43.3386,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cdd_Od;CJ] + [Ck_O;YJ] for rate rule [Ck_O;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=C([O])C([O])C=O(1100)'],
    products = ['C=C([O])[C](O)C=O(4809)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.66008e+10,'s^-1'), n=0.776642, Ea=(125.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;Cs_H_out_noH] + [R2H_S;O_rad_out;Cs_H_out] for rate rule [R2H_S;O_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=C([O])C([O])C=O(1100)'],
    products = ['C=C([O])C(O)[C]=O(4810)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;O_rad_out;XH_out] for rate rule [R3H_SS_Cs;O_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=C([O])C([O])C=O(1100)'],
    products = ['C=C(O)[C]([O])C=O(4811)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.19923e+08,'s^-1'), n=1.46351, Ea=(157.11,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;XH_out] + [R3H_SS_2Cd;Y_rad_out;XH_out] for rate rule [R3H_SS_2Cd;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=C(O)C([O])C=O(4812)'],
    products = ['C=C([O])C([O])C=O(1100)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C(O)C([O])[C]=O(4813)'],
    products = ['C=C([O])C([O])C=O(1100)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(26.449,'s^-1'), n=2.8625, Ea=(89.0146,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;O_H_out] + [R4H_SS(Cd)S;Y_rad_out;XH_out] for rate rule [R4H_SS(Cd)S;CO_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=C([O])C(O)C=O(4814)'],
    products = ['C=C([O])C([O])C=O(1100)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=[C][O](173)', '[O]C=C[O](2537)'],
    products = ['C=C([O])C([O])C=O(1100)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.49215e+07,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=C([O])C([O])C=O(1100)'],
    products = ['[O]C(C=O)[C]1CO1(4815)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=C([O])C([O])C=O(1100)'],
    products = ['[O][CH]C1OCC1=O(4754)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=C([O])C([O])C=O(1100)'],
    products = ['C=C([O])C1[CH]OO1(4816)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(265.782,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=C([O])C([O])C=O(1100)'],
    products = ['C=C1OO[CH]C1[O](4722)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.81184e+09,'s^-1'), n=0.551229, Ea=(327.115,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra] for rate rule [R5_SS_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic
Ea raised from 324.0 to 327.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C([O])C([O])C=O(1100)'],
    products = ['C=C(O)C(=O)C=O(4817)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CO(12)', 'C=C([O])C[O](939)'],
    products = ['C=C([O])C([O])C=O(1100)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.118397,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=C([O])C([O])C=O(1100)'],
    products = ['C=C1OOC1C=O(4732)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(70.8984,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination
Ea raised from 67.4 to 70.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=C([O])O[CH]C=O(1101)'],
    products = ['C=C([O])C([O])C=O(1100)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=C([O])C([O])=CO(4818)'],
    products = ['C=C([O])C([O])C=O(1100)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(163.689,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol
Ea raised from 160.3 to 163.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction25',
    reactants = ['O(4)', 'C=C([O])[CH]C=O(4597)'],
    products = ['C=C([O])C([O])C=O(1100)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/TwoDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['O(4)', 'C=[C]C([O])C=O(4819)'],
    products = ['C=C([O])C([O])C=O(1100)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=C([O])C([O])C=O(1100)'],
    products = ['[O]C1CC(=O)C1[O](4740)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.46159e+06,'s^-1'), n=1.55572, Ea=(119.171,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5_SS_CO;carbonylbond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 115.7 to 119.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=C([O])C([O])C=O(1100)'],
    products = ['CC(=O)[C]([O])C=O(4820)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.30996,'s^-1'), n=3.5644, Ea=(117.059,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=C([O])C([O])C=O(1100)'],
    products = ['CC(=O)C([O])[C]=O(4821)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(463.959,'s^-1'), n=2.50105, Ea=(76.0245,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_2H;XH_out] for rate rule [R4H_SSS;C_rad_out_2H;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=C([O])C([O])C=O(1100)'],
    products = ['C=C1OOC1[CH][O](4733)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(398.083,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 397.3 to 398.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=C([O])C([O])C=O(1100)'],
    products = ['[O]C1[CH]OCC1=O(4764)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.64784e+07,'s^-1'), n=0.990488, Ea=(61.3871,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5_SS_CO;carbonyl_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic
Ea raised from 58.5 to 61.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=C([O])C([O])C=O(1100)'],
    products = ['CC(=O)C(=O)C=O(4822)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]C(C=O)C[C]=O(1104)'],
    products = ['C=C([O])C([O])C=O(1100)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=C([O])C([O])C=O(1100)'],
    products = ['O=CC1OCC1=O(1107)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R4_SSS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C(=C[O])OC=O(4823)'],
    products = ['C=C([O])C([O])C=O(1100)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction36',
    reactants = ['CH2(19)', '[O]C([C]=O)C=O(4824)'],
    products = ['C=C([O])C([O])C=O(1100)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '679',
    isomers = [
        'C=C([O])C([O])C=O(1100)',
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
    label = '679',
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

