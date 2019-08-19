species(
    label = 'C=C[CH]CC([O])C=O(20575)',
    structure = SMILES('[CH2]C=CCC([O])C=O'),
    E0 = (66.2237,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0708243,'amu*angstrom^2'), symmetry=1, barrier=(1.62839,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.961657,'amu*angstrom^2'), symmetry=1, barrier=(22.1104,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0446347,'amu*angstrom^2'), symmetry=1, barrier=(22.0982,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0446231,'amu*angstrom^2'), symmetry=1, barrier=(22.1089,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.690801,0.0620479,-3.07983e-05,-3.37569e-09,5.27838e-12,8093.19,33.916], Tmin=(100,'K'), Tmax=(1065.88,'K')), NASAPolynomial(coeffs=[15.0474,0.0282137,-1.13897e-05,2.12637e-09,-1.49845e-13,3894.17,-41.603], Tmin=(1065.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(66.2237,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(Allyl_P)"""),
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
    label = '[CH2]C=CCC1OC1[O](21182)',
    structure = SMILES('[CH2]C=CCC1OC1[O]'),
    E0 = (117.739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.803828,0.0628873,-3.15422e-05,-5.36306e-09,7.91124e-12,14282.6,28.6078], Tmin=(100,'K'), Tmax=(922.012,'K')), NASAPolynomial(coeffs=[12.6992,0.0304679,-1.00142e-05,1.64072e-09,-1.075e-13,11273.5,-32.2371], Tmin=(922.012,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.739,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Ethylene_oxide) + radical(Allyl_P) + radical(CCOJ)"""),
)

species(
    label = 'C=CC1CC([CH][O])O1(20994)',
    structure = SMILES('C=CC1CC([CH][O])O1'),
    E0 = (180.675,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.896316,0.0602157,-2.44768e-05,-1.3134e-08,1.09283e-11,21849.4,28.4601], Tmin=(100,'K'), Tmax=(913.672,'K')), NASAPolynomial(coeffs=[12.7769,0.0296171,-9.39795e-06,1.51518e-09,-9.88626e-14,18784.6,-32.6754], Tmin=(913.672,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(180.675,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Oxetane) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = '[O]C1CC=CCC1[O](21183)',
    structure = SMILES('[O]C1CC=CCC1[O]'),
    E0 = (94.5118,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22797,0.041158,4.17111e-05,-8.38494e-08,3.47367e-11,11484.7,23.4124], Tmin=(100,'K'), Tmax=(977.818,'K')), NASAPolynomial(coeffs=[16.7574,0.0261177,-9.59622e-06,1.84271e-09,-1.37722e-13,6129.73,-63.0125], Tmin=(977.818,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(94.5118,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclohexene) + radical(CC(C)OJ) + radical(CC(C)OJ)"""),
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
    label = '[CH2]C=CCC(=O)C=O(21184)',
    structure = SMILES('[CH2]C=CCC(=O)C=O'),
    E0 = (-85.657,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,375,552.5,462.5,1710,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0229577,'amu*angstrom^2'), symmetry=1, barrier=(23.7502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0229811,'amu*angstrom^2'), symmetry=1, barrier=(23.7472,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.400482,'amu*angstrom^2'), symmetry=1, barrier=(23.7497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03299,'amu*angstrom^2'), symmetry=1, barrier=(23.7504,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.744489,0.0620956,-3.74344e-05,7.9066e-09,-4.44275e-14,-10177.2,27.9066], Tmin=(100,'K'), Tmax=(1418.6,'K')), NASAPolynomial(coeffs=[17.7458,0.0265175,-1.28842e-05,2.51117e-09,-1.75967e-13,-16244.5,-64.4435], Tmin=(1418.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-85.657,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(Allyl_P)"""),
)

species(
    label = 'C=C=CCC([O])C=O(21185)',
    structure = SMILES('C=C=CCC([O])C=O'),
    E0 = (91.3288,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,362.915,363.249,363.298],'cm^-1')),
        HinderedRotor(inertia=(0.00128078,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127326,'amu*angstrom^2'), symmetry=1, barrier=(11.9158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127307,'amu*angstrom^2'), symmetry=1, barrier=(11.9198,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.654488,0.0649196,-5.01871e-05,1.95959e-08,-3.06843e-12,11111.9,32.8297], Tmin=(100,'K'), Tmax=(1512.52,'K')), NASAPolynomial(coeffs=[16.0969,0.0240805,-9.68567e-06,1.74418e-09,-1.17739e-13,6440.54,-48.0604], Tmin=(1512.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(91.3288,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=OCOJ)"""),
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
    label = 'C=C[CH]CC=O(18209)',
    structure = SMILES('[CH2]C=CCC=O'),
    E0 = (22.3348,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,288.847],'cm^-1')),
        HinderedRotor(inertia=(0.0348698,'amu*angstrom^2'), symmetry=1, barrier=(20.5876,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.897468,'amu*angstrom^2'), symmetry=1, barrier=(20.6346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0347442,'amu*angstrom^2'), symmetry=1, barrier=(20.6385,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90908,0.0343717,1.38992e-05,-3.66954e-08,1.41587e-11,2771.34,21.2657], Tmin=(100,'K'), Tmax=(1086.76,'K')), NASAPolynomial(coeffs=[11.6938,0.0257615,-1.20409e-05,2.42052e-09,-1.77266e-13,-973.64,-34.1992], Tmin=(1086.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(22.3348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=CC[C](O)C=O(21186)',
    structure = SMILES('[CH2]C=CCC(O)=C[O]'),
    E0 = (-77.6652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0435552,0.0687399,-1.2055e-05,-5.15574e-08,2.98477e-11,-9176.75,30.7507], Tmin=(100,'K'), Tmax=(930.982,'K')), NASAPolynomial(coeffs=[25.5672,0.010306,-1.05039e-06,1.01001e-10,-1.23914e-14,-16181.7,-102.986], Tmin=(930.982,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-77.6652,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = 'C[C]=CCC([O])C=O(21187)',
    structure = SMILES('C[C]=CCC([O])C=O'),
    E0 = (152.566,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05077,0.0645996,-4.85567e-05,1.92861e-08,-3.21799e-12,18455.5,32.8106], Tmin=(100,'K'), Tmax=(1365.18,'K')), NASAPolynomial(coeffs=[11.5663,0.0337888,-1.47029e-05,2.75394e-09,-1.90499e-13,15584.4,-21.1938], Tmin=(1365.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(152.566,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Cds_S) + radical(C=OCOJ)"""),
)

species(
    label = '[CH2]C=C[CH]C(O)C=O(21188)',
    structure = SMILES('[CH2]C=C[CH]C(O)C=O'),
    E0 = (-60.5927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.279066,0.068392,-3.29995e-05,-9.73214e-09,9.23018e-12,-7141.91,31.811], Tmin=(100,'K'), Tmax=(1022.25,'K')), NASAPolynomial(coeffs=[18.5585,0.0247524,-9.88428e-06,1.87894e-09,-1.35646e-13,-12336.2,-63.9057], Tmin=(1022.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-60.5927,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(C=CCJCO) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=CCC(O)[C]=O(21189)',
    structure = SMILES('[CH2]C=CCC(O)[C]=O'),
    E0 = (-17.5488,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.433812,0.0673096,-4.00588e-05,5.24818e-10,5.39999e-12,-1972.52,34.8742], Tmin=(100,'K'), Tmax=(1019.46,'K')), NASAPolynomial(coeffs=[17.0737,0.0241917,-9.23863e-06,1.7034e-09,-1.20552e-13,-6517.37,-51.375], Tmin=(1019.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-17.5488,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(Allyl_P)"""),
)

species(
    label = 'CC=[C]CC([O])C=O(21190)',
    structure = SMILES('CC=[C]CC([O])C=O'),
    E0 = (152.566,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,180,621.258,4000],'cm^-1')),
        HinderedRotor(inertia=(0.293069,'amu*angstrom^2'), symmetry=1, barrier=(13.3814,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0952395,'amu*angstrom^2'), symmetry=1, barrier=(2.18974,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0488685,'amu*angstrom^2'), symmetry=1, barrier=(13.3814,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.582005,'amu*angstrom^2'), symmetry=1, barrier=(13.3814,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05077,0.0645996,-4.85567e-05,1.92861e-08,-3.21799e-12,18455.5,32.8106], Tmin=(100,'K'), Tmax=(1365.18,'K')), NASAPolynomial(coeffs=[11.5663,0.0337888,-1.47029e-05,2.75394e-09,-1.90499e-13,15584.4,-21.1938], Tmin=(1365.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(152.566,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Cds_S) + radical(C=OCOJ)"""),
)

species(
    label = '[CH2]C=[C]CC(O)C=O(21191)',
    structure = SMILES('[CH2]C=[C]CC(O)C=O'),
    E0 = (60.3324,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.251368,0.0712984,-5.67674e-05,2.27273e-08,-3.6204e-12,7400.56,35.2013], Tmin=(100,'K'), Tmax=(1500.81,'K')), NASAPolynomial(coeffs=[18.0503,0.0238589,-9.35218e-06,1.66464e-09,-1.11757e-13,2058.16,-57.8941], Tmin=(1500.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(60.3324,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'CC=C[CH]C([O])C=O(21192)',
    structure = SMILES('CC=C[CH]C([O])C=O'),
    E0 = (31.6412,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.25722,0.0707285,-5.34887e-05,2.01671e-08,-3.03152e-12,3949.99,32.4117], Tmin=(100,'K'), Tmax=(1578.95,'K')), NASAPolynomial(coeffs=[18.3641,0.0248576,-9.91124e-06,1.76772e-09,-1.18277e-13,-1767.96,-63.2139], Tmin=(1578.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.6412,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(C=CCJCO)"""),
)

species(
    label = '[CH2][C]=CCC(O)C=O(21193)',
    structure = SMILES('[CH2][C]=CCC(O)C=O'),
    E0 = (60.3324,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.251368,0.0712984,-5.67674e-05,2.27273e-08,-3.6204e-12,7400.56,35.2013], Tmin=(100,'K'), Tmax=(1500.81,'K')), NASAPolynomial(coeffs=[18.0503,0.0238589,-9.35218e-06,1.66464e-09,-1.11757e-13,2058.16,-57.8941], Tmin=(1500.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(60.3324,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'CC=CC[C]([O])C=O(21194)',
    structure = SMILES('CC=CCC([O])=C[O]'),
    E0 = (-91.3596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.302324,0.0675266,-2.92047e-05,-1.94129e-08,1.49226e-11,-10842.4,30.7288], Tmin=(100,'K'), Tmax=(952.576,'K')), NASAPolynomial(coeffs=[19.6721,0.0194868,-5.98917e-06,1.03425e-09,-7.40263e-14,-16043.3,-69.7077], Tmin=(952.576,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-91.3596,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = 'CC=CCC([O])[C]=O(21195)',
    structure = SMILES('CC=CCC([O])[C]=O'),
    E0 = (74.6851,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,1855,455,950,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,294.703,295.429,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00193146,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208809,'amu*angstrom^2'), symmetry=1, barrier=(12.8952,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.207642,'amu*angstrom^2'), symmetry=1, barrier=(12.8989,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208702,'amu*angstrom^2'), symmetry=1, barrier=(12.8941,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.655606,0.0670262,-5.24777e-05,2.13541e-08,-3.52856e-12,9108.23,34.5827], Tmin=(100,'K'), Tmax=(1426.63,'K')), NASAPolynomial(coeffs=[14.7157,0.0276043,-1.10283e-05,1.98462e-09,-1.34287e-13,5096.53,-38.2449], Tmin=(1426.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.6851,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(C=OCOJ)"""),
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
    label = '[CH]=C[CH2](2461)',
    structure = SMILES('[CH]C=C'),
    E0 = (376.654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,229.711,230.18,230.787],'cm^-1')),
        HinderedRotor(inertia=(1.33306,'amu*angstrom^2'), symmetry=1, barrier=(50.5153,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31912,0.00817959,3.34736e-05,-4.36194e-08,1.58213e-11,45331.5,10.6389], Tmin=(100,'K'), Tmax=(983.754,'K')), NASAPolynomial(coeffs=[5.36755,0.0170743,-6.35108e-06,1.1662e-09,-8.2762e-14,44095,-3.44606], Tmin=(983.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(376.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C([O])C=O(2583)',
    structure = SMILES('[CH2]C([O])C=O'),
    E0 = (82.1174,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.347137,'amu*angstrom^2'), symmetry=1, barrier=(18.8168,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0159131,'amu*angstrom^2'), symmetry=1, barrier=(18.8156,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.0627,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.37294,0.0300135,-1.22984e-05,-6.17239e-09,4.64176e-12,9940.08,21.4594], Tmin=(100,'K'), Tmax=(1011.88,'K')), NASAPolynomial(coeffs=[9.84794,0.0127268,-4.85036e-06,8.96671e-10,-6.36509e-14,7799.56,-17.7935], Tmin=(1011.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(82.1174,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CJCO)"""),
)

species(
    label = '[O]C(C=O)CC1[CH]C1(21196)',
    structure = SMILES('[O]C(C=O)CC1[CH]C1'),
    E0 = (178.655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20956,0.0481862,4.48201e-06,-3.7174e-08,1.65691e-11,21599.2,33.3573], Tmin=(100,'K'), Tmax=(1019.54,'K')), NASAPolynomial(coeffs=[13.7,0.02916,-1.16307e-05,2.20194e-09,-1.58095e-13,17494.3,-34.7843], Tmin=(1019.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(178.655,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(cyclopropane) + radical(C=OCOJ)"""),
)

species(
    label = '[CH2]C=CCC1[CH]OO1(21197)',
    structure = SMILES('[CH2]C=CCC1[CH]OO1'),
    E0 = (331.008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.852895,0.0662204,-4.61379e-05,1.63509e-08,-2.39476e-12,39926.4,25.5481], Tmin=(100,'K'), Tmax=(1548.48,'K')), NASAPolynomial(coeffs=[13.5935,0.0333094,-1.42577e-05,2.62569e-09,-1.78871e-13,35980.6,-41.4895], Tmin=(1548.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(331.008,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(12dioxetane) + radical(CCsJOO) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C1[CH]CC(C=O)O1(21021)',
    structure = SMILES('[CH2]C1[CH]CC(C=O)O1'),
    E0 = (58.7932,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33594,0.0387477,5.11508e-05,-1.0291e-07,4.58685e-11,7185.65,26.9224], Tmin=(100,'K'), Tmax=(907.552,'K')), NASAPolynomial(coeffs=[17.5662,0.0193999,-3.12469e-06,3.19221e-10,-2.11744e-14,2090.52,-61.6452], Tmin=(907.552,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(58.7932,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Tetrahydrofuran) + radical(CJC(C)OC) + radical(CCJCO)"""),
)

species(
    label = '[O]C1[CH]OCC=CC1(21050)',
    structure = SMILES('[O]C1[CH]OCC=CC1'),
    E0 = (119.808,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02065,0.0128758,0.000192204,-3.01411e-07,1.29761e-10,14566.3,27.4351], Tmin=(100,'K'), Tmax=(902.467,'K')), NASAPolynomial(coeffs=[40.1501,-0.0220198,1.99394e-05,-4.05624e-09,2.67655e-13,1862.12,-188.582], Tmin=(902.467,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(119.808,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptane) + radical(CC(C)OJ) + radical(CCsJOCs)"""),
)

species(
    label = 'CC=CCC(=O)C=O(21198)',
    structure = SMILES('CC=CCC(=O)C=O'),
    E0 = (-237.156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.908553,0.0624539,-3.77562e-05,9.80779e-09,-9.42815e-13,-28407.2,27.0541], Tmin=(100,'K'), Tmax=(1912.59,'K')), NASAPolynomial(coeffs=[22.9404,0.022242,-1.08191e-05,2.02191e-09,-1.34699e-13,-37907.5,-96.3279], Tmin=(1912.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-237.156,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H)"""),
)

species(
    label = 'C=C=CCC(O)C=O(21199)',
    structure = SMILES('C=C=CCC(O)C=O'),
    E0 = (-152.404,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.484671,0.067127,-4.44551e-05,9.78377e-09,8.61507e-13,-18194.9,32.651], Tmin=(100,'K'), Tmax=(1117.6,'K')), NASAPolynomial(coeffs=[15.9887,0.0269116,-1.09812e-05,2.04567e-09,-1.43227e-13,-22614.3,-48.1386], Tmin=(1117.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-152.404,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
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
    label = 'C=C[CH]CC[O](18939)',
    structure = SMILES('[CH2]C=CCC[O]'),
    E0 = (164.829,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,245.521,2444.63,2444.63],'cm^-1')),
        HinderedRotor(inertia=(0.243124,'amu*angstrom^2'), symmetry=1, barrier=(10.4,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.766284,'amu*angstrom^2'), symmetry=1, barrier=(32.7788,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.766282,'amu*angstrom^2'), symmetry=1, barrier=(32.7788,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3713.27,'J/mol'), sigma=(6.31777,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=580.00 K, Pc=33.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98027,0.0464601,-2.5237e-05,6.38765e-09,-6.46696e-13,19894.2,22.9639], Tmin=(100,'K'), Tmax=(2131.54,'K')), NASAPolynomial(coeffs=[12.6767,0.0263877,-1.11119e-05,1.9699e-09,-1.28561e-13,15334.1,-36.7357], Tmin=(2131.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(164.829,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(CCOJ)"""),
)

species(
    label = 'O=CC1CC=CCO1(20579)',
    structure = SMILES('O=CC1CC=CCO1'),
    E0 = (-236.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69323,0.0295134,7.06349e-05,-1.13654e-07,4.61201e-11,-28387.6,23.4328], Tmin=(100,'K'), Tmax=(945.825,'K')), NASAPolynomial(coeffs=[15.3636,0.0252602,-7.56209e-06,1.33535e-09,-9.89582e-14,-33369.3,-54.422], Tmin=(945.825,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-236.88,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(36dihydro2hpyran)"""),
)

species(
    label = 'C=C[CH]CO[CH]C=O(20577)',
    structure = SMILES('C=C[CH]COC=C[O]'),
    E0 = (-55.6986,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3025,407.5,1350,352.5,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3996.75,'J/mol'), sigma=(6.61607,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=624.28 K, Pc=31.31 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.42488,0.0739037,-1.03914e-05,-6.13178e-08,3.48769e-11,-6518.07,27.9978], Tmin=(100,'K'), Tmax=(931.931,'K')), NASAPolynomial(coeffs=[28.9635,0.00752677,2.54328e-07,-1.21246e-10,7.31996e-16,-14590.8,-125.636], Tmin=(931.931,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-55.6986,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CCJCO) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C=CCC([O])=CO(21200)',
    structure = SMILES('[CH2]C=CCC([O])=CO'),
    E0 = (-81.3231,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3615,1277.5,1000,3000,3100,440,815,1455,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,206.162,206.162,206.162],'cm^-1')),
        HinderedRotor(inertia=(0.623155,'amu*angstrom^2'), symmetry=1, barrier=(18.7949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.623154,'amu*angstrom^2'), symmetry=1, barrier=(18.7949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.623153,'amu*angstrom^2'), symmetry=1, barrier=(18.7949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.623154,'amu*angstrom^2'), symmetry=1, barrier=(18.7949,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0954925,0.0718123,-2.32327e-05,-4.03572e-08,2.64269e-11,-9616.47,31.0284], Tmin=(100,'K'), Tmax=(918.837,'K')), NASAPolynomial(coeffs=[25.1625,0.00984288,-4.0641e-07,-7.99389e-11,3.25029e-15,-16283.7,-99.7118], Tmin=(918.837,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-81.3231,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(C=C(C)OJ)"""),
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
    label = '[CH2]C=CC[CH]C=O(20123)',
    structure = SMILES('[CH2]C=CCC=C[O]'),
    E0 = (136.893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3000,3100,440,815,1455,1000,288.83,288.868,289.433],'cm^-1')),
        HinderedRotor(inertia=(0.382724,'amu*angstrom^2'), symmetry=1, barrier=(22.4549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.376442,'amu*angstrom^2'), symmetry=1, barrier=(22.4702,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.379495,'amu*angstrom^2'), symmetry=1, barrier=(22.4623,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.949875,0.0501072,1.01307e-05,-5.57741e-08,2.67736e-11,16589.7,27.1272], Tmin=(100,'K'), Tmax=(955.034,'K')), NASAPolynomial(coeffs=[18.4853,0.0184131,-5.66331e-06,1.02492e-09,-7.696e-14,11336.3,-66.6323], Tmin=(955.034,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(136.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_P)"""),
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
    label = '[CH]=CCC([O])C=O(21201)',
    structure = SMILES('[CH]=CCC([O])C=O'),
    E0 = (197.846,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,3120,650,792.5,1650,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,320.393,320.394],'cm^-1')),
        HinderedRotor(inertia=(0.00164225,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00164225,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131975,'amu*angstrom^2'), symmetry=1, barrier=(9.61354,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25989,0.052517,-3.54592e-05,8.84922e-09,1.79001e-13,23900.6,30.1913], Tmin=(100,'K'), Tmax=(1130.29,'K')), NASAPolynomial(coeffs=[13.1489,0.0213805,-8.65384e-06,1.60058e-09,-1.1143e-13,20514.3,-31.7133], Tmin=(1130.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(197.846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(Cds_P)"""),
)

species(
    label = 'C=CC1CC([O])C1[O](20971)',
    structure = SMILES('C=CC1CC([O])C1[O]'),
    E0 = (201.904,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02987,0.0454711,3.28333e-05,-7.83599e-08,3.384e-11,24408.2,27.4961], Tmin=(100,'K'), Tmax=(969.237,'K')), NASAPolynomial(coeffs=[18.1802,0.0233047,-8.0946e-06,1.53848e-09,-1.15911e-13,18800.3,-66.4874], Tmin=(969.237,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(201.904,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(CC(C)OJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=CC=CC([O])C=O(21202)',
    structure = SMILES('C=CC=CC([O])C=O'),
    E0 = (25.1675,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,377.958,377.96,377.96],'cm^-1')),
        HinderedRotor(inertia=(0.123558,'amu*angstrom^2'), symmetry=1, barrier=(12.5252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123557,'amu*angstrom^2'), symmetry=1, barrier=(12.5252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123557,'amu*angstrom^2'), symmetry=1, barrier=(12.5252,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.738389,0.0610613,-3.26415e-05,-3.61792e-09,6.24977e-12,3153.7,31.3791], Tmin=(100,'K'), Tmax=(1021.13,'K')), NASAPolynomial(coeffs=[15.8842,0.0233444,-8.9849e-06,1.66462e-09,-1.18092e-13,-1066.27,-47.5249], Tmin=(1021.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(25.1675,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ)"""),
)

species(
    label = 'C=CC[CH]C([O])C=O(20990)',
    structure = SMILES('C=CC[CH]C([O])C=O'),
    E0 = (126.872,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,344.369,344.978,345.183,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00141448,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00144358,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00143331,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.307154,'amu*angstrom^2'), symmetry=1, barrier=(25.4167,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.82474,0.0592781,-2.67447e-05,-6.06316e-09,5.96439e-12,15382.5,36.0893], Tmin=(100,'K'), Tmax=(1063.56,'K')), NASAPolynomial(coeffs=[14.4466,0.0282102,-1.13654e-05,2.12192e-09,-1.49573e-13,11344.6,-35.8289], Tmin=(1063.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(126.872,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCJCO) + radical(C=OCOJ)"""),
)

species(
    label = 'C=[C]CCC([O])C=O(21203)',
    structure = SMILES('C=[C]CCC([O])C=O'),
    E0 = (164.812,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,1685,370,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,180,180,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0299689,'amu*angstrom^2'), symmetry=1, barrier=(15.3807,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0101345,'amu*angstrom^2'), symmetry=1, barrier=(5.19393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0300443,'amu*angstrom^2'), symmetry=1, barrier=(15.382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.66909,'amu*angstrom^2'), symmetry=1, barrier=(15.3837,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.830493,0.06636,-5.05732e-05,2.00142e-08,-3.26159e-12,19938.8,34.0434], Tmin=(100,'K'), Tmax=(1421.41,'K')), NASAPolynomial(coeffs=[13.4688,0.0307942,-1.30407e-05,2.41068e-09,-1.65433e-13,16346,-31.3734], Tmin=(1421.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(164.812,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=OCOJ)"""),
)

species(
    label = 'C=CCC[C]([O])C=O(21204)',
    structure = SMILES('C=CCCC([O])=C[O]'),
    E0 = (-80.3838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.251152,0.0674629,-2.52188e-05,-2.60209e-08,1.78627e-11,-9519.3,31.324], Tmin=(100,'K'), Tmax=(947.043,'K')), NASAPolynomial(coeffs=[20.6608,0.018037,-5.18586e-06,8.82892e-10,-6.41426e-14,-15034.3,-74.7378], Tmin=(947.043,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-80.3838,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = '[CH]=CCCC([O])C=O(21205)',
    structure = SMILES('[CH]=CCCC([O])C=O'),
    E0 = (174.066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.508079,0.0684453,-5.28267e-05,2.07524e-08,-3.28288e-12,21067.8,35.1411], Tmin=(100,'K'), Tmax=(1493.43,'K')), NASAPolynomial(coeffs=[16.1952,0.0264288,-1.06253e-05,1.91367e-09,-1.29265e-13,16382.3,-46.8318], Tmin=(1493.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(174.066,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=OCOJ)"""),
)

species(
    label = 'C=CCCC([O])[C]=O(21206)',
    structure = SMILES('C=CCCC([O])[C]=O'),
    E0 = (86.9305,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,1855,455,950,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,180,180,695.29,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0369054,'amu*angstrom^2'), symmetry=1, barrier=(12.6651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.40367,'amu*angstrom^2'), symmetry=1, barrier=(12.665,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0368921,'amu*angstrom^2'), symmetry=1, barrier=(12.6657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0369224,'amu*angstrom^2'), symmetry=1, barrier=(12.6659,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.419145,0.0689358,-5.48621e-05,2.24048e-08,-3.66406e-12,10592.4,35.8769], Tmin=(100,'K'), Tmax=(1459.47,'K')), NASAPolynomial(coeffs=[16.5393,0.0247547,-9.45386e-06,1.66277e-09,-1.11038e-13,5887.05,-47.988], Tmin=(1459.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(86.9305,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O) + radical(C=OCOJ)"""),
)

species(
    label = '[CH]=C[CH]CC(O)C=O(21207)',
    structure = SMILES('[CH]C=CCC(O)C=O'),
    E0 = (41.6759,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.445437,0.0673999,-3.49538e-05,6.10137e-10,3.55542e-12,5149.48,34.8472], Tmin=(100,'K'), Tmax=(1109.41,'K')), NASAPolynomial(coeffs=[14.7418,0.0338384,-1.38927e-05,2.56635e-09,-1.7822e-13,870.623,-40.5969], Tmin=(1109.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(41.6759,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=CC1CC([O])[CH]O1(21010)',
    structure = SMILES('C=CC1CC([O])[CH]O1'),
    E0 = (107.027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.998962,0.0456139,4.0448e-05,-9.82184e-08,4.59408e-11,12999.6,24.7367], Tmin=(100,'K'), Tmax=(900.342,'K')), NASAPolynomial(coeffs=[19.913,0.0161155,-1.25896e-06,-6.36584e-11,6.3358e-15,7383.51,-76.8014], Tmin=(900.342,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(107.027,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Tetrahydrofuran) + radical(CCsJOCs) + radical(CC(C)OJ)"""),
)

species(
    label = '[O][CH]C1CC=CCO1(21043)',
    structure = SMILES('[O][CH]C1CC=CCO1'),
    E0 = (89.5063,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.963627,0.0591365,-3.11201e-05,4.10459e-09,1.45217e-12,10880.8,24.6807], Tmin=(100,'K'), Tmax=(1114.75,'K')), NASAPolynomial(coeffs=[10.8954,0.034635,-1.31358e-05,2.3109e-09,-1.55363e-13,7974.62,-27.4168], Tmin=(1114.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(89.5063,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(36dihydro2hpyran) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = 'C=CC=CC(O)C=O(21208)',
    structure = SMILES('C=CC=CC(O)C=O'),
    E0 = (-218.566,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.530618,0.0636005,-2.75205e-05,-1.34196e-08,1.04981e-11,-26151.3,31.3441], Tmin=(100,'K'), Tmax=(1001.82,'K')), NASAPolynomial(coeffs=[17.4558,0.0235928,-8.89747e-06,1.65702e-09,-1.19006e-13,-30926,-57.246], Tmin=(1001.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-218.566,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CCCC(=O)C=O(21209)',
    structure = SMILES('C=CCCC(=O)C=O'),
    E0 = (-231.093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20988,0.0672105,-5.4352e-05,2.51081e-08,-5.19421e-12,-27698.5,26.8148], Tmin=(100,'K'), Tmax=(1077.62,'K')), NASAPolynomial(coeffs=[7.84004,0.0426002,-2.00957e-05,3.91569e-09,-2.77748e-13,-29127.5,-5.66758], Tmin=(1077.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-231.093,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-O2d)H) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2]C(C=C)C([O])C=O(20576)',
    structure = SMILES('[CH2]C(C=C)C([O])C=O'),
    E0 = (124.08,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,180,846.518,856.197],'cm^-1')),
        HinderedRotor(inertia=(0.00641395,'amu*angstrom^2'), symmetry=1, barrier=(3.3419,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154816,'amu*angstrom^2'), symmetry=1, barrier=(3.55953,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0237421,'amu*angstrom^2'), symmetry=1, barrier=(12.2235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.529537,'amu*angstrom^2'), symmetry=1, barrier=(12.1751,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4168.81,'J/mol'), sigma=(6.84665,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=651.16 K, Pc=29.47 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.716861,0.0627554,-3.43241e-05,-9.22572e-10,5.18698e-12,15049.9,35.5134], Tmin=(100,'K'), Tmax=(1009.33,'K')), NASAPolynomial(coeffs=[14.4505,0.0273971,-1.01148e-05,1.80447e-09,-1.24549e-13,11306.2,-35.6825], Tmin=(1009.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(124.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(C=OCOJ)"""),
)

species(
    label = 'C=CC1CC(C=O)O1(20580)',
    structure = SMILES('C=CC1CC(C=O)O1'),
    E0 = (-145.711,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49979,0.0320194,7.26372e-05,-1.25507e-07,5.36322e-11,-17413.6,27.668], Tmin=(100,'K'), Tmax=(916.514,'K')), NASAPolynomial(coeffs=[18.2027,0.0186556,-2.92599e-06,3.30337e-10,-2.52754e-14,-22975.7,-65.0985], Tmin=(916.514,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-145.711,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Oxetane)"""),
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
    E0 = (66.2237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (168.514,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (188.267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (207.643,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (158.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (314.92,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (76.9822,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (148.744,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (220.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (273.955,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (208.201,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (141.504,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (314.487,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (211.827,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (334.209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (177.669,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (154.506,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (155.436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (255.868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (458.771,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (292.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (332.01,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (123.576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (133.536,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (91.197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (91.197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (350.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (73.7549,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (258.101,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (66.2237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (379.898,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (579.409,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (201.904,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (238.55,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (82.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (251.974,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (365.008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (259.245,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (334.418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (189.514,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (194.944,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (222.649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (127.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (144.471,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (155.192,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (284.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (74.1315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C[CH]CC([O])C=O(20575)'],
    products = ['OCHCHO(48)', 'butadiene13(2459)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C[CH]CC([O])C=O(20575)'],
    products = ['[CH2]C=CCC1OC1[O](21182)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(7.785e+11,'s^-1'), n=0.342, Ea=(102.29,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=C[CH]CC([O])C=O(20575)'],
    products = ['C=CC1CC([CH][O])O1(20994)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=C[CH]CC([O])C=O(20575)'],
    products = ['[O]C1CC=CCC1[O](21183)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(8.21347e+10,'s^-1'), n=0.43, Ea=(141.419,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_SMSR;multiplebond_intra;radadd_intra_cs2H] for rate rule [R7_SMSS_CO;carbonylbond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.82842712475
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', '[CH2]C=CCC(=O)C=O(21184)'],
    products = ['C=C[CH]CC([O])C=O(20575)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(25.1243,'m^3/(mol*s)'), n=1.86, Ea=(32.426,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO-DeNd_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', 'C=C=CCC([O])C=O(21185)'],
    products = ['C=C[CH]CC([O])C=O(20575)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.42e+08,'cm^3/(mol*s)'), n=1.64, Ea=(11.7989,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2713 used for Ca_Cds-HH;HJ
Exact match found for rate rule [Ca_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['OCHCHO(48)', '[CH2][CH]C=C(2458)'],
    products = ['C=C[CH]CC([O])C=O(20575)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.000217286,'m^3/(mol*s)'), n=3.00879, Ea=(27.5684,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;CsJ-OneDeHH] for rate rule [CO-DeH_O;CsJ-CdHH]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 4.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['HCO(16)', 'C=C[CH]CC=O(18209)'],
    products = ['C=C[CH]CC([O])C=O(20575)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(5.2e+11,'cm^3/(mol*s)'), n=0, Ea=(93.9308,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO_O;CO_pri_rad] for rate rule [CO-CsH_O;CO_pri_rad]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=C[CH]CC([O])C=O(20575)'],
    products = ['[CH2]C=CC[C](O)C=O(21186)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.70223e+09,'s^-1'), n=1.15155, Ea=(153.908,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_OneDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_CO]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=C[CH]CC([O])C=O(20575)'],
    products = ['C[C]=CCC([O])C=O(21187)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=C[CH]CC([O])C=O(20575)'],
    products = ['[CH2]C=C[CH]C(O)C=O(21188)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.07519e+07,'s^-1'), n=1.60667, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_H/Cd] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=C[CH]CC([O])C=O(20575)'],
    products = ['[CH2]C=CCC(O)[C]=O(21189)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;O_rad_out;XH_out] for rate rule [R3H_SS_Cs;O_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CC=[C]CC([O])C=O(21190)'],
    products = ['C=C[CH]CC([O])C=O(20575)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=C[CH]CC([O])C=O(20575)'],
    products = ['[CH2]C=[C]CC(O)C=O(21191)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(274,'s^-1'), n=3.09, Ea=(145.603,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=C[CH]CC([O])C=O(20575)'],
    products = ['CC=C[CH]C([O])C=O(21192)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.15542e+07,'s^-1'), n=1, Ea=(267.985,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SDS;C_rad_out_single;Cs_H_out_H/NonDeC] + [R4H_SDS;C_rad_out_2H;Cs_H_out_1H] for rate rule [R4H_SDS;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=C[CH]CC([O])C=O(20575)'],
    products = ['[CH2][C]=CCC(O)C=O(21193)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.234e+06,'s^-1','*|/',3), n=1.554, Ea=(111.445,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R5H_SSSD;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=C[CH]CC([O])C=O(20575)'],
    products = ['CC=CC[C]([O])C=O(21194)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(126000,'s^-1'), n=1.85, Ea=(88.2824,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5H_SMSS;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CC=CCC([O])[C]=O(21195)'],
    products = ['C=C[CH]CC([O])C=O(20575)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(37753.8,'s^-1'), n=1.925, Ea=(80.7512,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_RSSMS;Y_rad_out;Cs_H_out_2H] for rate rule [R6H_RSSMS;CO_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][CH]C=C(2458)', '[O]C=C[O](2537)'],
    products = ['C=C[CH]CC([O])C=O(20575)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(6.26649e+07,'m^3/(mol*s)'), n=0.074875, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [C_rad/H2/Cd;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: R_Recombination
Ea raised from -0.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C[CH2](2461)', '[CH2]C([O])C=O(2583)'],
    products = ['C=C[CH]CC([O])C=O(20575)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.83556e+08,'m^3/(mol*s)'), n=-0.424942, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;C_rad/H2/Cs] + [Cd_rad;C_pri_rad] for rate rule [Cd_rad;C_rad/H2/Cs]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=C[CH]CC([O])C=O(20575)'],
    products = ['[O]C(C=O)CC1[CH]C1(21196)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra_cs2H] for rate rule [R3_D;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=C[CH]CC([O])C=O(20575)'],
    products = ['[CH2]C=CCC1[CH]OO1(21197)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(265.786,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C[CH]CC([O])C=O(20575)'],
    products = ['[CH2]C1[CH]CC(C=O)O1(21021)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.66591e+07,'s^-1'), n=1.01661, Ea=(57.3526,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=C[CH]CC([O])C=O(20575)'],
    products = ['[O]C1[CH]OCC=CC1(21050)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.97949e+09,'s^-1'), n=0.649948, Ea=(67.3127,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7_linear;multiplebond_intra;radadd_intra_cs2H] + [R7_linear;carbonyl_intra_H;radadd_intra] for rate rule [R7_linear;carbonyl_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=C[CH]CC([O])C=O(20575)'],
    products = ['CC=CCC(=O)C=O(21198)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=C[CH]CC([O])C=O(20575)'],
    products = ['C=C=CCC(O)C=O(21199)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CO(12)', 'C=C[CH]CC[O](18939)'],
    products = ['C=C[CH]CC([O])C=O(20575)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(0.118397,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=C[CH]CC([O])C=O(20575)'],
    products = ['O=CC1CC=CCO1(20579)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), Tmin=(550,'K'), Tmax=(650,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R6_SSSDS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=C[CH]CO[CH]C=O(20577)'],
    products = ['C=C[CH]CC([O])C=O(20575)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C=CCC([O])=CO(21200)'],
    products = ['C=C[CH]CC([O])C=O(20575)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(147.547,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol
Ea raised from 145.6 to 147.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction31',
    reactants = ['O(4)', '[CH2]C=CC[CH]C=O(20123)'],
    products = ['C=C[CH]CC([O])C=O(20575)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/OneDeC;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction32',
    reactants = ['CH2(19)', '[CH]=CCC([O])C=O(21201)'],
    products = ['C=C[CH]CC([O])C=O(20575)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=C[CH]CC([O])C=O(20575)'],
    products = ['C=CC1CC([O])C1[O](20971)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(902977,'s^-1'), n=1.63829, Ea=(135.68,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs] for rate rule [R5_SS_CO;carbonylbond_intra_H;radadd_intra_csHCd]
Euclidian distance = 3.0
family: Intra_R_Add_Exocyclic
Ea raised from 133.9 to 135.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction34',
    reactants = ['H(3)', 'C=CC=CC([O])C=O(21202)'],
    products = ['C=C[CH]CC([O])C=O(20575)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.35e+08,'cm^3/(mol*s)'), n=1.64, Ea=(1.58992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2557 used for Cds-CsH_Cds-CdH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction35',
    reactants = ['butadiene13(2459)', '[O]C=C[O](2537)'],
    products = ['C=C[CH]CC([O])C=O(20575)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(0.106847,'m^3/(mol*s)'), n=2.459, Ea=(4.91982,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CdH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C=CC[CH]C([O])C=O(20990)'],
    products = ['C=C[CH]CC([O])C=O(20575)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.682e+10,'s^-1'), n=0.35, Ea=(125.102,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 160 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=[C]CCC([O])C=O(21203)'],
    products = ['C=C[CH]CC([O])C=O(20575)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.9054e+11,'s^-1'), n=0.853, Ea=(200.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C=C[CH]CC([O])C=O(20575)'],
    products = ['C=CCC[C]([O])C=O(21204)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.32587e+10,'s^-1'), n=0.723333, Ea=(193.022,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_H/Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C=C[CH]CC([O])C=O(20575)'],
    products = ['[CH]=CCCC([O])C=O(21205)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(8.32e+10,'s^-1'), n=0.77, Ea=(268.194,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 195 used for R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleH
Exact match found for rate rule [R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C=CCCC([O])[C]=O(21206)'],
    products = ['C=C[CH]CC([O])C=O(20575)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(657459,'s^-1'), n=1.87237, Ea=(102.584,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R4H_SSS;CO_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C=C[CH]CC([O])C=O(20575)'],
    products = ['[CH]=C[CH]CC(O)C=O(21207)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(3.427,'s^-1'), n=3.311, Ea=(128.721,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;O_rad_out;Cd_H_out_singleH] for rate rule [R6HJ_3;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['C=C[CH]CC([O])C=O(20575)'],
    products = ['C=CC1CC([O])[CH]O1(21010)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(6.8435e+15,'s^-1'), n=-1.17677, Ea=(156.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_csHCd] for rate rule [R5_SS_CO;carbonyl_intra_H;radadd_intra_csHCd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction43',
    reactants = ['C=C[CH]CC([O])C=O(20575)'],
    products = ['[O][CH]C1CC=CCO1(21043)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(9.91671e+09,'s^-1'), n=0.30082, Ea=(60.8864,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra_pri_2H;radadd_intra] for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction44',
    reactants = ['C=C[CH]CC([O])C=O(20575)'],
    products = ['C=CC=CC(O)C=O(21208)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction45',
    reactants = ['C=C[CH]CC([O])C=O(20575)'],
    products = ['C=CCCC(=O)C=O(21209)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(C=C)C([O])C=O(20576)'],
    products = ['C=C[CH]CC([O])C=O(20575)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction47',
    reactants = ['C=C[CH]CC([O])C=O(20575)'],
    products = ['C=CC1CC(C=O)O1(20580)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

network(
    label = '3676',
    isomers = [
        'C=C[CH]CC([O])C=O(20575)',
    ],
    reactants = [
        ('OCHCHO(48)', 'butadiene13(2459)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '3676',
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

