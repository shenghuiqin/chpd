species(
    label = 'C7H7O(470)(469)',
    structure = SMILES('[CH]1C=C2OC2C=CC1'),
    E0 = (169.422,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4069.49,'J/mol'), sigma=(6.651,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=635.64 K, Pc=31.39 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.53274,-0.0059557,0.000184551,-2.36061e-07,8.89148e-11,20464.2,17.6772], Tmin=(100,'K'), Tmax=(962.819,'K')), NASAPolynomial(coeffs=[20.1541,0.0186091,-6.041e-06,1.37592e-09,-1.22272e-13,12539.1,-90.2024], Tmin=(962.819,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(169.422,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s2_3_7_ane) - ring(Cycloheptane) - ring(Ethylene_oxide) + ring(1,4-Cycloheptadiene) + ring(Ethylene_oxide) + radical(Allyl_S)"""),
)

species(
    label = 'C7H7O(473)(472)',
    structure = SMILES('[CH]1CC2C=C3OC3C12'),
    E0 = (288.357,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3934.62,'J/mol'), sigma=(6.54226,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=614.58 K, Pc=31.88 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75729,0.0267855,7.15394e-05,-1.08904e-07,4.16667e-11,34782,17.8466], Tmin=(100,'K'), Tmax=(1000.09,'K')), NASAPolynomial(coeffs=[16.1378,0.0255001,-1.08718e-05,2.25269e-09,-1.73695e-13,29093.5,-65.5918], Tmin=(1000.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(288.357,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_1) + Estimated bicyclic component: polycyclic(s2_3_5_ane) - ring(Cyclopentane) - ring(Ethylene_oxide) + ring(Cyclopentene) + ring(Ethylene_oxide) - ring(Cyclopentene) + radical(cyclobutane)"""),
)

species(
    label = 'C7H7O(476)(475)',
    structure = SMILES('[O]C1C=CCC=CC=1'),
    E0 = (89.4153,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4422.05,'J/mol'), sigma=(6.98257,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=690.71 K, Pc=29.47 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.277261,0.0660949,-2.22778e-05,-3.09166e-08,2.02327e-11,10902.7,21.0868], Tmin=(100,'K'), Tmax=(939.654,'K')), NASAPolynomial(coeffs=[21.5531,0.0149921,-3.7016e-06,6.01876e-10,-4.50583e-14,5161.97,-89.5036], Tmin=(939.654,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(89.4153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C7H7O(477)(476)',
    structure = SMILES('[O]C1C=CC2CC2C=1'),
    E0 = (119.969,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4298.03,'J/mol'), sigma=(6.88159,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=671.34 K, Pc=29.93 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50151,0.0369302,3.97388e-05,-7.87999e-08,3.28524e-11,14535.1,19.764], Tmin=(100,'K'), Tmax=(968.373,'K')), NASAPolynomial(coeffs=[15.5218,0.0234322,-8.15075e-06,1.53248e-09,-1.14138e-13,9737.23,-58.1775], Tmin=(968.373,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(119.969,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_1_3) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C7H7O(478)(477)',
    structure = SMILES('[CH]1CC=CC2=CC1O2'),
    E0 = (113.202,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4069.49,'J/mol'), sigma=(6.651,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=635.64 K, Pc=31.39 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.395696,0.0423229,0.000123745,-2.34618e-07,1.02831e-10,13823.2,-1.94004], Tmin=(100,'K'), Tmax=(943.353,'K')), NASAPolynomial(coeffs=[48.7128,-0.0276317,1.51106e-05,-2.46667e-09,1.31047e-13,-1594.75,-268.607], Tmin=(943.353,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(113.202,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Cyclobutene) + radical(CCJCO)"""),
)

species(
    label = 'C7H7O(481)(480)',
    structure = SMILES('O=C1C=CC2C[CH]C12'),
    E0 = (216.387,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4297.64,'J/mol'), sigma=(6.68421,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=671.28 K, Pc=32.65 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.04598,0.0293115,4.08646e-05,-6.24334e-08,2.21792e-11,26107.3,20.8423], Tmin=(100,'K'), Tmax=(1068.74,'K')), NASAPolynomial(coeffs=[9.73992,0.035063,-1.56965e-05,3.09528e-09,-2.24596e-13,22489.7,-26.0187], Tmin=(1068.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.387,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + polycyclic(s2_4_5_ene_1) + radical(CCJCC=O)"""),
)

species(
    label = 'C7H7O(484)(483)',
    structure = SMILES('[CH]1C2=CC=CCC1O2'),
    E0 = (11.0608,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4069.49,'J/mol'), sigma=(6.651,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=635.64 K, Pc=31.39 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.462793,0.0407404,0.000143636,-2.68545e-07,1.19987e-10,1544.45,-6.65751], Tmin=(100,'K'), Tmax=(920.216,'K')), NASAPolynomial(coeffs=[51.0191,-0.0337956,2.18551e-05,-4.07128e-09,2.53522e-13,-14249.5,-285.082], Tmin=(920.216,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(11.0608,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Oxetane) + radical(C=CCJCO)"""),
)

species(
    label = 'C7H7O(486)(485)',
    structure = SMILES('O=C1[CH]C2CC=CC12'),
    E0 = (172.492,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4297.64,'J/mol'), sigma=(6.68421,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=671.28 K, Pc=32.65 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77487,0.0337575,3.78269e-05,-6.53493e-08,2.46939e-11,20839.2,19.0451], Tmin=(100,'K'), Tmax=(1032.15,'K')), NASAPolynomial(coeffs=[11.7068,0.0322586,-1.37532e-05,2.68875e-09,-1.95968e-13,16818.5,-38.7301], Tmin=(1032.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(172.492,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_1) + radical(CCJC=O)"""),
)

species(
    label = 'C7H7O(495)(494)',
    structure = SMILES('OC1C=C[CH]C=CC=1'),
    E0 = (53.3354,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4422.05,'J/mol'), sigma=(6.98257,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=690.71 K, Pc=29.47 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.243558,0.0574981,2.36137e-05,-9.38989e-08,4.67038e-11,6573.41,18.6852], Tmin=(100,'K'), Tmax=(917.623,'K')), NASAPolynomial(coeffs=[27.4491,0.00461139,2.66105e-06,-6.458e-10,3.87173e-14,-1185.74,-125.3], Tmin=(917.623,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(53.3354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=CCJC=C)"""),
)

species(
    label = 'H(3)(3)',
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
    label = 'C7H6O(492)(491)',
    structure = SMILES('O=C1C=CC=CC=C1'),
    E0 = (69.2232,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4447.44,'J/mol'), sigma=(6.55857,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=694.68 K, Pc=35.77 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8255,-0.0103058,0.000244307,-3.52476e-07,1.47305e-10,8458.46,20.2001], Tmin=(100,'K'), Tmax=(906.109,'K')), NASAPolynomial(coeffs=[40.2319,-0.0280939,2.25315e-05,-4.4704e-09,2.90851e-13,-4731.48,-195.679], Tmin=(906.109,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(69.2232,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cd-Cd(CO)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(Cycloheptane)"""),
)

species(
    label = 'C7H6O(488)(487)',
    structure = SMILES('[O]C1C=C[CH]C=CC=1'),
    E0 = (191.14,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4422.05,'J/mol'), sigma=(6.98257,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=690.71 K, Pc=29.47 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.605953,0.0535798,1.58346e-05,-7.63285e-08,3.83895e-11,23130.6,18.885], Tmin=(100,'K'), Tmax=(919.229,'K')), NASAPolynomial(coeffs=[23.7188,0.00795663,6.12307e-07,-2.55728e-10,1.34208e-14,16559.8,-103.303], Tmin=(919.229,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.14,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=CCJC=C) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C7H6O(487)(486)',
    structure = SMILES('C1=CC=C2OC2C=C1'),
    E0 = (103.464,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4079.36,'J/mol'), sigma=(6.42564,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=637.19 K, Pc=34.89 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.95496,0.0472955,-1.58347e-05,-2.24639e-09,1.11315e-12,12459.4,7.21113], Tmin=(100,'K'), Tmax=(2039.03,'K')), NASAPolynomial(coeffs=[27.6972,0.0218854,-1.41552e-05,2.76709e-09,-1.83559e-13,-2438.3,-141.573], Tmin=(2039.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(103.464,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s2_3_7_ane) - ring(Cycloheptane) - ring(Ethylene_oxide) + ring(1,3,5-Cycloheptatriene) + ring(Ethylene_oxide)"""),
)

species(
    label = 'S(801)(800)',
    structure = SMILES('C1C=CC2C=C(C=1)O2'),
    E0 = (8.87563,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4079.36,'J/mol'), sigma=(6.42564,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=637.19 K, Pc=34.89 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.450588,0.043012,0.000120285,-2.33946e-07,1.03575e-10,1278.23,-6.25478], Tmin=(100,'K'), Tmax=(940.725,'K')), NASAPolynomial(coeffs=[50.1611,-0.0322972,1.73032e-05,-2.8873e-09,1.60273e-13,-14434.1,-280.235], Tmin=(940.725,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(8.87563,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Cyclobutene)"""),
)

species(
    label = 'C7H6O(535)(534)',
    structure = SMILES('O=C1C=CC2C=CC12'),
    E0 = (149.769,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4308.03,'J/mol'), sigma=(6.45303,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=672.90 K, Pc=36.38 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46784,0.0441003,7.28496e-08,-2.73546e-08,1.19664e-11,18114.1,21.1405], Tmin=(100,'K'), Tmax=(1069.85,'K')), NASAPolynomial(coeffs=[12.8397,0.0270689,-1.17812e-05,2.29891e-09,-1.66169e-13,14222.3,-41.3066], Tmin=(1069.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(149.769,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + polycyclic(s2_4_5_diene_1_5)"""),
)

species(
    label = 'O=C1C=C=CCC=C1(2904)',
    structure = SMILES('O=C1C=C=CCC=C1'),
    E0 = (203.723,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96739,0.0370966,6.26737e-06,-2.24968e-08,7.6995e-12,24581.4,18.2703], Tmin=(100,'K'), Tmax=(1234.59,'K')), NASAPolynomial(coeffs=[8.91927,0.0351095,-1.62699e-05,3.14665e-09,-2.21534e-13,21299.8,-23.0721], Tmin=(1234.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(203.723,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)(Cds-Cds)) + group(Cd-Cd(CO)H) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene)"""),
)

species(
    label = '[O]C1C=CC[C]=CC=1(2905)',
    structure = SMILES('[O]C1C=CC[C]=CC=1'),
    E0 = (327.257,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.246637,0.070671,-4.70103e-05,-1.6633e-09,9.24222e-12,39505.9,21.6702], Tmin=(100,'K'), Tmax=(943.141,'K')), NASAPolynomial(coeffs=[20.3521,0.014511,-3.98979e-06,6.53317e-10,-4.65847e-14,34418.8,-81.0139], Tmin=(943.141,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(327.257,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=C(C)OJ) + radical(Cds_S)"""),
)

species(
    label = '[O]C1C=[C]CC=CC=1(2906)',
    structure = SMILES('[O]C1C=[C]CC=CC=1'),
    E0 = (327.257,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.246637,0.070671,-4.70103e-05,-1.6633e-09,9.24222e-12,39505.9,21.6702], Tmin=(100,'K'), Tmax=(943.141,'K')), NASAPolynomial(coeffs=[20.3521,0.014511,-3.98979e-06,6.53317e-10,-4.65847e-14,34418.8,-81.0139], Tmin=(943.141,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(327.257,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=C(C)OJ) + radical(Cds_S)"""),
)

species(
    label = '[O]C1C=CCC=[C]C=1(2907)',
    structure = SMILES('[O]C1C=CCC=[C]C=1'),
    E0 = (288.411,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.547638,0.0813027,-8.10564e-05,3.96784e-08,-7.32737e-12,34867.5,23.6532], Tmin=(100,'K'), Tmax=(1504.38,'K')), NASAPolynomial(coeffs=[21.4039,0.0126043,-2.25678e-06,1.93326e-10,-7.06618e-15,29431.9,-87.3296], Tmin=(1504.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(288.411,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=C(C)OJ) + radical(C=CJC=C)"""),
)

species(
    label = '[O]C1[C]=CCC=CC=1(2908)',
    structure = SMILES('[O]C1[C]=CCC=CC=1'),
    E0 = (288.411,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.547635,0.0813027,-8.10563e-05,3.96783e-08,-7.32735e-12,34867.5,23.6532], Tmin=(100,'K'), Tmax=(1504.39,'K')), NASAPolynomial(coeffs=[21.4039,0.0126044,-2.25685e-06,1.93339e-10,-7.06722e-15,29432,-87.3291], Tmin=(1504.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(288.411,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=CJC=C) + radical(C=C(C)OJ)"""),
)

species(
    label = '[O]C1=[C]C=CCC=C1(2909)',
    structure = SMILES('[O]C1=[C]C=CCC=C1'),
    E0 = (288.411,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.547635,0.0813027,-8.10563e-05,3.96783e-08,-7.32735e-12,34867.5,23.6532], Tmin=(100,'K'), Tmax=(1504.39,'K')), NASAPolynomial(coeffs=[21.4039,0.0126044,-2.25685e-06,1.93339e-10,-7.06722e-15,29432,-87.3291], Tmin=(1504.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(288.411,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=CJC=C) + radical(C=C(C)OJ)"""),
)

species(
    label = '[O]C1C=CC=CC=C1(1734)',
    structure = SMILES('[O]C1C=CC=CC=C1'),
    E0 = (218.031,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.432555,0.0601522,-4.45784e-06,-4.73366e-08,2.49157e-11,26368.1,21.9776], Tmin=(100,'K'), Tmax=(962.06,'K')), NASAPolynomial(coeffs=[21.654,0.0162295,-5.0626e-06,9.56867e-10,-7.44163e-14,20234.2,-90.2402], Tmin=(962.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(218.031,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]C1C=C=CCCC=1(2910)',
    structure = SMILES('[O]C1C=C=CCCC=1'),
    E0 = (219.662,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19479,0.0467979,1.24195e-05,-5.06399e-08,2.29828e-11,26533.6,18.6876], Tmin=(100,'K'), Tmax=(973.799,'K')), NASAPolynomial(coeffs=[15.2508,0.0248269,-8.82955e-06,1.62362e-09,-1.17412e-13,22100.2,-57.4586], Tmin=(973.799,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(219.662,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(C=C(C)OJ)"""),
)

species(
    label = 'OC1[C]=CCC=CC=1(2215)',
    structure = SMILES('OC1[C]=CCC=CC=1'),
    E0 = (150.606,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.35309,0.0902546,-9.0107e-05,4.30324e-08,-7.62278e-12,18330.1,25.0571], Tmin=(100,'K'), Tmax=(1618.78,'K')), NASAPolynomial(coeffs=[24.4866,0.00996358,-4.72605e-07,-1.56214e-10,1.61272e-14,12118.5,-105.397], Tmin=(1618.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(150.606,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=CJC=C)"""),
)

species(
    label = 'OC1=[C]C=CCC=C1(2216)',
    structure = SMILES('OC1=[C]C=CCC=C1'),
    E0 = (150.606,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.35307,0.0902544,-9.01065e-05,4.30319e-08,-7.62265e-12,18330.1,25.057], Tmin=(100,'K'), Tmax=(1618.8,'K')), NASAPolynomial(coeffs=[24.4859,0.00996456,-4.73098e-07,-1.5611e-10,1.61192e-14,12118.9,-105.393], Tmin=(1618.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(150.606,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=CJC=C)"""),
)

species(
    label = 'OC1C=[C]CC=CC=1(2213)',
    structure = SMILES('OC1C=[C]CC=CC=1'),
    E0 = (189.452,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.10964,0.0745168,-3.89797e-05,-1.95504e-08,1.76841e-11,22948.4,21.4484], Tmin=(100,'K'), Tmax=(933.064,'K')), NASAPolynomial(coeffs=[24.0504,0.0112204,-1.97267e-06,2.70734e-10,-2.19104e-14,16686.6,-102.831], Tmin=(933.064,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(189.452,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(Cds_S)"""),
)

species(
    label = 'OC1C=CCC=[C]C=1(2217)',
    structure = SMILES('OC1C=CCC=[C]C=1'),
    E0 = (150.606,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.35309,0.0902546,-9.0107e-05,4.30324e-08,-7.62278e-12,18330.1,25.0571], Tmin=(100,'K'), Tmax=(1618.78,'K')), NASAPolynomial(coeffs=[24.4866,0.00996358,-4.72605e-07,-1.56214e-10,1.61272e-14,12118.5,-105.397], Tmin=(1618.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(150.606,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=CJC=C)"""),
)

species(
    label = 'OC1C=CC[C]=CC=1(2214)',
    structure = SMILES('OC1C=CC[C]=CC=1'),
    E0 = (189.452,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.10964,0.0745168,-3.89797e-05,-1.95504e-08,1.76841e-11,22948.4,21.4484], Tmin=(100,'K'), Tmax=(933.064,'K')), NASAPolynomial(coeffs=[24.0504,0.0112204,-1.97267e-06,2.70734e-10,-2.19104e-14,16686.6,-102.831], Tmin=(933.064,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(189.452,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(Cds_S)"""),
)

species(
    label = 'C1=C[C]2C=CC(C1)O2(1698)',
    structure = SMILES('C1=C[C]2C=CC(C1)O2'),
    E0 = (221.601,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05107,0.0256382,5.82574e-05,-8.36503e-08,3.00196e-11,26737.4,20.6097], Tmin=(100,'K'), Tmax=(1047.26,'K')), NASAPolynomial(coeffs=[11.5923,0.033048,-1.51662e-05,3.07378e-09,-2.27959e-13,22334.2,-37.3432], Tmin=(1047.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(221.601,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s3_5_6_diene_1_5) + radical(C2CsJOCs)"""),
)

species(
    label = '[CH]1CC=CC=C2OC12(2911)',
    structure = SMILES('[CH]1CC=CC=C2OC12'),
    E0 = (215.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.01356,0.0464034,-1.10164e-05,-5.46622e-09,1.74993e-12,25960.8,-0.401604], Tmin=(100,'K'), Tmax=(1969.96,'K')), NASAPolynomial(coeffs=[25.2064,0.0276382,-1.67517e-05,3.25107e-09,-2.16316e-13,12114.4,-135.467], Tmin=(1969.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(215.74,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s2_3_7_ane) - ring(Cycloheptane) - ring(Ethylene_oxide) + ring(1,3-Cycloheptadiene) + ring(Ethylene_oxide) + radical(CCJCO)"""),
)

species(
    label = 'O(4)(4)',
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
    label = '[C]1=CC=CCC=C1(1246)',
    structure = SMILES('[C]1=CC=CCC=C1'),
    E0 = (364.499,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (91.1305,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.898556,0.0536037,-4.67814e-06,-4.06382e-08,2.19826e-11,43964.1,17.7097], Tmin=(100,'K'), Tmax=(939.78,'K')), NASAPolynomial(coeffs=[18.1007,0.0169825,-4.63849e-06,7.70412e-10,-5.5963e-14,39114.8,-72.8103], Tmin=(939.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(364.499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=CJC=C)"""),
)

species(
    label = 'O=C1C=CCC=[C]C1(2176)',
    structure = SMILES('O=C1C=CCC=[C]C1'),
    E0 = (231.55,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05377,-0.0114719,0.000236838,-3.27636e-07,1.32364e-10,27969.5,24.4607], Tmin=(100,'K'), Tmax=(924.411,'K')), NASAPolynomial(coeffs=[34.7942,-0.0132749,1.28068e-05,-2.39298e-09,1.39171e-13,15940.3,-163.243], Tmin=(924.411,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(231.55,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + ring(Cycloheptane) + radical(Cds_S)"""),
)

species(
    label = 'O=C1C=CC[C]=CC1(2177)',
    structure = SMILES('O=C1C=CC[C]=CC1'),
    E0 = (231.55,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05377,-0.0114719,0.000236838,-3.27636e-07,1.32364e-10,27969.5,24.4607], Tmin=(100,'K'), Tmax=(924.411,'K')), NASAPolynomial(coeffs=[34.7942,-0.0132749,1.28068e-05,-2.39298e-09,1.39171e-13,15940.3,-163.243], Tmin=(924.411,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(231.55,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + ring(Cycloheptane) + radical(Cds_S)"""),
)

species(
    label = 'O=C1[C]=CCC=CC1(2179)',
    structure = SMILES('O=C1[C]=CCC=CC1'),
    E0 = (237.67,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93751,-0.0110604,0.000241423,-3.36294e-07,1.36228e-10,28711.7,25.1093], Tmin=(100,'K'), Tmax=(925.135,'K')), NASAPolynomial(coeffs=[36.761,-0.0162394,1.40924e-05,-2.60691e-09,1.51912e-13,16046.7,-173.81], Tmin=(925.135,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(237.67,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + ring(Cycloheptane) + radical(C=CJC=O)"""),
)

species(
    label = 'O=C1C=C[CH]C=CC1(2175)',
    structure = SMILES('O=C1C=C[CH]C=CC1'),
    E0 = (95.4335,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.39604,-0.0283561,0.000298935,-4.01306e-07,1.61079e-10,11594.9,22.4294], Tmin=(100,'K'), Tmax=(919.97,'K')), NASAPolynomial(coeffs=[38.2305,-0.0199491,1.74787e-05,-3.31861e-09,2.00559e-13,-1947.45,-185.23], Tmin=(919.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(95.4335,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + ring(Cycloheptane) + radical(C=CCJC=C)"""),
)

species(
    label = 'O=C1C=[C]CC=CC1(2178)',
    structure = SMILES('O=C1C=[C]CC=CC1'),
    E0 = (231.55,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05377,-0.0114719,0.000236838,-3.27636e-07,1.32364e-10,27969.5,24.4607], Tmin=(100,'K'), Tmax=(924.411,'K')), NASAPolynomial(coeffs=[34.7942,-0.0132749,1.28068e-05,-2.39298e-09,1.39171e-13,15940.3,-163.243], Tmin=(924.411,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(231.55,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + ring(Cycloheptane) + radical(Cds_S)"""),
)

species(
    label = 'O=C1C=CCC2[CH]C12(2912)',
    structure = SMILES('O=C1C=CCC2[CH]C12'),
    E0 = (198.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90315,0.0347505,2.33894e-05,-4.35182e-08,1.54506e-11,23925.7,20.5852], Tmin=(100,'K'), Tmax=(1110.23,'K')), NASAPolynomial(coeffs=[9.49767,0.0356811,-1.60934e-05,3.14392e-09,-2.25395e-13,20495.7,-24.7009], Tmin=(1110.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(198.223,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + polycyclic(s2_3_6_ene_2) + radical(CCJCC=O)"""),
)

species(
    label = 'O=C1C2[CH]CC=CC12(2913)',
    structure = SMILES('O=C1C2[CH]CC=CC12'),
    E0 = (194.691,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69074,0.0401182,1.00331e-05,-3.2753e-08,1.26835e-11,23508,21.1655], Tmin=(100,'K'), Tmax=(1088.66,'K')), NASAPolynomial(coeffs=[10.1734,0.0333583,-1.42822e-05,2.73058e-09,-1.93994e-13,20214.7,-27.1221], Tmin=(1088.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(194.691,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_3_6_ene_1) + radical(CCJCC=O)"""),
)

species(
    label = '[CH]=CCC=CC=C=O(2914)',
    structure = SMILES('[CH]=CCC=CC=C=O'),
    E0 = (322.111,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2120,512.5,787.5,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.56137,'amu*angstrom^2'), symmetry=1, barrier=(12.907,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.56146,'amu*angstrom^2'), symmetry=1, barrier=(12.9091,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.561902,'amu*angstrom^2'), symmetry=1, barrier=(12.9192,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11603,0.0653993,-5.58017e-05,2.61991e-08,-5.21432e-12,38843.1,28.5316], Tmin=(100,'K'), Tmax=(1163.2,'K')), NASAPolynomial(coeffs=[10.2935,0.0338401,-1.51049e-05,2.87471e-09,-2.01373e-13,36708,-17.132], Tmin=(1163.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(322.111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CCO)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'O=[C]C1C=CCC=C1(2915)',
    structure = SMILES('O=[C]C1C=CCC=C1'),
    E0 = (136.627,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,1855,455,950,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4417,0.0447387,3.93264e-06,-3.14139e-08,1.3417e-11,16534.4,25.0785], Tmin=(100,'K'), Tmax=(1050.4,'K')), NASAPolynomial(coeffs=[11.8695,0.0307049,-1.26926e-05,2.4088e-09,-1.71591e-13,12927.3,-32.4849], Tmin=(1050.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(136.627,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(1,4-Cyclohexadiene) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[O]C12C=CCC1C=C2(2429)',
    structure = SMILES('[O]C12C=CCC1C=C2'),
    E0 = (293.871,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78751,0.021698,9.49996e-05,-1.41716e-07,5.58916e-11,35448.5,24.7605], Tmin=(100,'K'), Tmax=(964.465,'K')), NASAPolynomial(coeffs=[18.8474,0.0188705,-6.24635e-06,1.29203e-09,-1.06334e-13,28998.6,-73.3045], Tmin=(964.465,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293.871,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_diene_1_5) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'O=C1C=[C]CCC=C1(2916)',
    structure = SMILES('O=C1C=[C]CCC=C1'),
    E0 = (225.071,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95069,-0.0035602,0.000207024,-2.95402e-07,1.21744e-10,27189,25.3253], Tmin=(100,'K'), Tmax=(913.153,'K')), NASAPolynomial(coeffs=[32.6547,-0.0118606,1.33607e-05,-2.67195e-09,1.6789e-13,16320.1,-148.824], Tmin=(913.153,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(225.071,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)(Cds-Cds)) + ring(Cycloheptane) + radical(Cds_S)"""),
)

species(
    label = 'O=C1[C]=CCCC=C1(2917)',
    structure = SMILES('[O]C1=C=CCCC=C1'),
    E0 = (219.662,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19485,0.0467972,1.2422e-05,-5.06432e-08,2.29843e-11,26533.6,18.6874], Tmin=(100,'K'), Tmax=(973.789,'K')), NASAPolynomial(coeffs=[15.2506,0.0248273,-8.82976e-06,1.62367e-09,-1.17416e-13,22100.3,-57.4574], Tmin=(973.789,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(219.662,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=CC(=O)C=CC=C(2918)',
    structure = SMILES('[CH]=CC(=O)C=CC=C'),
    E0 = (275.713,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3120,650,792.5,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,375,552.5,462.5,1710,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.709665,'amu*angstrom^2'), symmetry=1, barrier=(16.3166,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.709779,'amu*angstrom^2'), symmetry=1, barrier=(16.3192,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.709603,'amu*angstrom^2'), symmetry=1, barrier=(16.3152,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11211,0.0651592,-5.20912e-05,2.17402e-08,-3.80534e-12,33262.9,26.6914], Tmin=(100,'K'), Tmax=(1306.05,'K')), NASAPolynomial(coeffs=[11.5829,0.033091,-1.52614e-05,2.94086e-09,-2.06877e-13,30527.7,-26.6201], Tmin=(1306.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(275.713,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C1C=CC(=O)C=C1(2919)',
    structure = SMILES('[CH2]C1C=CC(=O)C=C1'),
    E0 = (163.305,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51968,0.0527733,-2.97185e-05,8.20686e-09,-9.28549e-13,19731.1,22.5423], Tmin=(100,'K'), Tmax=(1920.95,'K')), NASAPolynomial(coeffs=[12.2463,0.0304374,-1.22773e-05,2.15393e-09,-1.40803e-13,15610.1,-36.2099], Tmin=(1920.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(163.305,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)(Cds-Cds)) + ring(25cyclohexadienone) + radical(Isobutyl)"""),
)

species(
    label = 'C1=CCC=CC2OC=12(4314)',
    structure = SMILES('C1=CCC=CC2OC=12'),
    E0 = (253.589,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.10097,0.0231296,2.68538e-05,-2.82606e-08,6.36557e-12,30473.8,6.65969], Tmin=(100,'K'), Tmax=(1726.32,'K')), NASAPolynomial(coeffs=[11.7044,0.0444115,-2.54381e-05,4.98622e-09,-3.39409e-13,22052.3,-50.9624], Tmin=(1726.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(253.589,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + Estimated bicyclic component: polycyclic(s2_3_7_ane) - ring(Cycloheptane) - ring(Ethylene_oxide) + ring(1_2_cycloheptadiene) + ring(Ethylene_oxide)"""),
)

species(
    label = '[CH]1C=CC2OC2=CC1(4315)',
    structure = SMILES('[CH]1C=CCC=C2OC12'),
    E0 = (103.479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.58824,-0.0172711,0.00023687,-3.03365e-07,1.1593e-10,12540.6,19.5721], Tmin=(100,'K'), Tmax=(950.055,'K')), NASAPolynomial(coeffs=[25.6149,0.00947886,-6.6724e-07,3.64742e-10,-5.72065e-14,2582.75,-119.718], Tmin=(950.055,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(103.479,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s2_3_7_ane) - ring(Cycloheptane) - ring(Ethylene_oxide) + ring(1,4-Cycloheptadiene) + ring(Ethylene_oxide) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = '[C]1CCC=CC2OC=12(4316)',
    structure = SMILES('[C]1CCC=CC2OC=12'),
    E0 = (266.151,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.30123,0.00387118,0.000150805,-1.98301e-07,7.49728e-11,32101.9,20.0886], Tmin=(100,'K'), Tmax=(968.385,'K')), NASAPolynomial(coeffs=[18.9172,0.0206029,-7.33987e-06,1.60102e-09,-1.34766e-13,24881.2,-80.2067], Tmin=(968.385,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(266.151,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s2_3_7_ane) - ring(Cycloheptane) - ring(Ethylene_oxide) + ring(1,4-Cycloheptadiene) + ring(Ethylene_oxide) + radical(Cds_S)"""),
)

species(
    label = '[C]1=CC2OC2=CCC1(4317)',
    structure = SMILES('[C]1=CC2OC2=CCC1'),
    E0 = (266.151,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.30123,0.00387118,0.000150805,-1.98301e-07,7.49728e-11,32101.9,20.0886], Tmin=(100,'K'), Tmax=(968.385,'K')), NASAPolynomial(coeffs=[18.9172,0.0206029,-7.33987e-06,1.60102e-09,-1.34766e-13,24881.2,-80.2067], Tmin=(968.385,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(266.151,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s2_3_7_ane) - ring(Cycloheptane) - ring(Ethylene_oxide) + ring(1,4-Cycloheptadiene) + ring(Ethylene_oxide) + radical(Cds_S)"""),
)

species(
    label = 'C1=C[C]2OC2=CCC1(4318)',
    structure = SMILES('[CH]1C=C2OC2=CCC1'),
    E0 = (181.923,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[5.86643,0.0298449,1.08865e-05,-1.56471e-08,3.19554e-12,21747.6,-13.6112], Tmin=(100,'K'), Tmax=(2143.36,'K')), NASAPolynomial(coeffs=[34.9669,0.0226778,-1.7089e-05,3.31583e-09,-2.13198e-13,-1555.16,-201.45], Tmin=(2143.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(181.923,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s2_3_7_ane) - ring(Cycloheptane) - ring(Ethylene_oxide) + ring(1,3-Cycloheptadiene) + ring(Ethylene_oxide) + radical(Allyl_S)"""),
)

species(
    label = '[C]1=CCCC=C2OC12(4319)',
    structure = SMILES('[C]1=CCCC=C2OC12'),
    E0 = (266.151,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.30123,0.00387118,0.000150805,-1.98301e-07,7.49728e-11,32101.9,20.0886], Tmin=(100,'K'), Tmax=(968.385,'K')), NASAPolynomial(coeffs=[18.9172,0.0206029,-7.33987e-06,1.60102e-09,-1.34766e-13,24881.2,-80.2067], Tmin=(968.385,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(266.151,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s2_3_7_ane) - ring(Cycloheptane) - ring(Ethylene_oxide) + ring(1,4-Cycloheptadiene) + ring(Ethylene_oxide) + radical(Cds_S)"""),
)

species(
    label = '[CH]1C=C2O[C]2C=CC1(4320)',
    structure = SMILES('[CH]1C=C2OC2=C[CH]C1'),
    E0 = (323.036,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[5.8505,0.0268503,1.31594e-05,-1.66842e-08,3.38894e-12,38723.6,-15.3306], Tmin=(100,'K'), Tmax=(2113.07,'K')), NASAPolynomial(coeffs=[32.0477,0.0237271,-1.76093e-05,3.43015e-09,-2.22072e-13,17278.3,-185.863], Tmin=(2113.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(323.036,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s2_3_7_ane) - ring(Cycloheptane) - ring(Ethylene_oxide) + ring(1,3-Cycloheptadiene) + ring(Ethylene_oxide) + radical(Allyl_S) + radical(Allyl_S)"""),
)

species(
    label = '[CH]1[CH]C=C2OC2C=C1(4321)',
    structure = SMILES('[CH]1C=C[CH]C2OC2=C1'),
    E0 = (205.204,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.90372,-0.0296094,0.000274252,-3.47642e-07,1.33506e-10,24769.1,18.1097], Tmin=(100,'K'), Tmax=(941.208,'K')), NASAPolynomial(coeffs=[27.7628,0.00246963,3.63311e-06,-4.89935e-10,1.0471e-15,13989.2,-132.722], Tmin=(941.208,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.204,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s2_3_7_ane) - ring(Cycloheptane) - ring(Ethylene_oxide) + ring(1,4-Cycloheptadiene) + ring(Ethylene_oxide) + radical(C=CCJC=C) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = '[C]1=CC[CH]C=C2OC12(4322)',
    structure = SMILES('[C]1=CC[CH]C=C2OC12'),
    E0 = (407.264,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.50433,-0.0014108,0.000159954,-2.07028e-07,7.80413e-11,49067.3,18.253], Tmin=(100,'K'), Tmax=(966.658,'K')), NASAPolynomial(coeffs=[18.9643,0.0181102,-6.31942e-06,1.42513e-09,-1.23618e-13,41790.8,-81.7762], Tmin=(966.658,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(407.264,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s2_3_7_ane) - ring(Cycloheptane) - ring(Ethylene_oxide) + ring(1,4-Cycloheptadiene) + ring(Ethylene_oxide) + radical(Cds_S) + radical(Allyl_S)"""),
)

species(
    label = '[C]1=CC2OC2=C[CH]C1(4323)',
    structure = SMILES('[C]1=CC2OC2=C[CH]C1'),
    E0 = (407.264,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.50433,-0.0014108,0.000159954,-2.07028e-07,7.80413e-11,49067.3,18.253], Tmin=(100,'K'), Tmax=(966.658,'K')), NASAPolynomial(coeffs=[18.9643,0.0181102,-6.31942e-06,1.42513e-09,-1.23618e-13,41790.8,-81.7762], Tmin=(966.658,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(407.264,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s2_3_7_ane) - ring(Cycloheptane) - ring(Ethylene_oxide) + ring(1,4-Cycloheptadiene) + ring(Ethylene_oxide) + radical(Cds_S) + radical(Allyl_S)"""),
)

species(
    label = '[C]1[CH]CC=CC2OC=12(4324)',
    structure = SMILES('[C]1[CH]CC=CC2OC=12'),
    E0 = (407.264,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.50433,-0.0014108,0.000159954,-2.07028e-07,7.80413e-11,49067.3,18.253], Tmin=(100,'K'), Tmax=(966.658,'K')), NASAPolynomial(coeffs=[18.9643,0.0181102,-6.31942e-06,1.42513e-09,-1.23618e-13,41790.8,-81.7762], Tmin=(966.658,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(407.264,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s2_3_7_ane) - ring(Cycloheptane) - ring(Ethylene_oxide) + ring(1,4-Cycloheptadiene) + ring(Ethylene_oxide) + radical(Cds_S) + radical(Allyl_S)"""),
)

species(
    label = '[CH]1C2CC=CC3OC123(4325)',
    structure = SMILES('[CH]1C2CC=CC3OC123'),
    E0 = (334.022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79349,0.0414287,2.92106e-06,-2.17121e-08,7.28642e-12,40258.5,17.6015], Tmin=(100,'K'), Tmax=(1332.45,'K')), NASAPolynomial(coeffs=[13.0327,0.0339955,-1.8326e-05,3.73572e-09,-2.68297e-13,34928,-48.6103], Tmin=(1332.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(334.022,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_3_6_ene_2) + polycyclic(s2_3_6_ene_1) + polycyclic(s1_3_3_ane) - ring(Cyclopropane) - ring(Ethylene_oxide) - ring(Cyclohexene) + radical(CCJCO)"""),
)

species(
    label = '[CH]1C2OC2=CC2CC12(4326)',
    structure = SMILES('[CH]1C2OC2=CC2CC12'),
    E0 = (307.21,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.84735,0.0202455,9.5743e-05,-1.46864e-07,5.98309e-11,37051.2,20.9235], Tmin=(100,'K'), Tmax=(940.522,'K')), NASAPolynomial(coeffs=[19.9093,0.0124406,-1.87358e-06,3.45844e-10,-3.65275e-14,30601.3,-81.3345], Tmin=(940.522,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(307.21,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + polycyclic(s2_3_6_ene_1) + Estimated bicyclic component: polycyclic(s2_3_6_ane) - ring(Cyclohexane) - ring(Ethylene_oxide) + ring(Cyclohexene) + ring(Ethylene_oxide) - ring(Cyclohexene) + radical(CCJCO)"""),
)

species(
    label = '[CH]=CC1OC1=CC=C(4327)',
    structure = SMILES('[CH]=CC1OC1=CC=C'),
    E0 = (381.082,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,180,180,180,180,1048.36,1161.79,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.184209,'amu*angstrom^2'), symmetry=1, barrier=(4.23533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184209,'amu*angstrom^2'), symmetry=1, barrier=(4.23533,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.406384,0.056071,1.53077e-05,-7.75338e-08,3.84961e-11,45984.1,25.6968], Tmin=(100,'K'), Tmax=(935.778,'K')), NASAPolynomial(coeffs=[25.4868,0.00721563,8.6125e-08,-5.43141e-11,-5.27611e-15,38735.3,-107.288], Tmin=(935.778,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.082,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(Cds_P)"""),
)

species(
    label = 'O(S)(1202)',
    structure = SMILES('O'),
    E0 = (432.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (15.9994,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,51997.4,2.99252], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,51997.4,2.99252], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(432.331,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH2]C1C=CC2OC2=C1(4328)',
    structure = SMILES('[CH2]C1C=CC2OC2=C1'),
    E0 = (328.948,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.908235,0.047432,2.49914e-05,-7.52177e-08,3.43232e-11,39693.3,24.7131], Tmin=(100,'K'), Tmax=(956.886,'K')), NASAPolynomial(coeffs=[20.7602,0.0153255,-4.43648e-06,8.52324e-10,-6.89636e-14,33564.7,-82.3578], Tmin=(956.886,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(328.948,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_3_6_diene_0_3) + radical(Isobutyl)"""),
)

species(
    label = '[C]1C=C2OC2C=CC1(4329)',
    structure = SMILES('[C]1C=C2OC2C=CC1'),
    E0 = (476.524,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2814,0.00233612,0.000156026,-2.078e-07,7.95115e-11,57406.6,19.5912], Tmin=(100,'K'), Tmax=(962.344,'K')), NASAPolynomial(coeffs=[20.9363,0.0151359,-4.73673e-06,1.11637e-09,-1.02369e-13,49633,-91.4262], Tmin=(962.344,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(476.524,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s2_3_7_ane) - ring(Cycloheptane) - ring(Ethylene_oxide) + ring(1,4-Cycloheptadiene) + ring(Ethylene_oxide) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C1=CC2OC=2C=CC1(4330)',
    structure = SMILES('C1=CC2OC=2C=CC1'),
    E0 = (355.837,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.982105,0.103719,-0.000253727,3.20454e-07,-1.40425e-10,42869.3,14.8434], Tmin=(100,'K'), Tmax=(820.154,'K')), NASAPolynomial(coeffs=[-18.2183,0.0969454,-5.76858e-05,1.18177e-08,-8.40998e-13,49396,124.256], Tmin=(820.154,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(355.837,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s2_3_7_ane) - ring(Cycloheptane) - ring(Ethylene_oxide) + ring(1,3,5-Cycloheptatriene) + ring(oxirene)"""),
)

species(
    label = '[C]1=CCC=CC2OC12(1656)',
    structure = SMILES('[C]1=CCC=CC2OC12'),
    E0 = (314.991,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.58226,-8.79434e-05,0.000151908,-1.91673e-07,7.05845e-11,37963.7,21.4901], Tmin=(100,'K'), Tmax=(980.081,'K')), NASAPolynomial(coeffs=[15.9385,0.0256445,-1.02866e-05,2.19264e-09,-1.75995e-13,31491.8,-62.3387], Tmin=(980.081,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(314.991,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s2_3_7_ane) - ring(Cycloheptane) - ring(Ethylene_oxide) + ring(1,4-Cycloheptadiene) + ring(Ethylene_oxide) + radical(Cds_S)"""),
)

species(
    label = '[C]1=CC2OC2C=CC1(1657)',
    structure = SMILES('[C]1=CC2OC2C=CC1'),
    E0 = (314.991,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.58226,-8.79434e-05,0.000151908,-1.91673e-07,7.05845e-11,37963.7,21.4901], Tmin=(100,'K'), Tmax=(980.081,'K')), NASAPolynomial(coeffs=[15.9385,0.0256445,-1.02866e-05,2.19264e-09,-1.75995e-13,31491.8,-62.3387], Tmin=(980.081,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(314.991,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s2_3_7_ane) - ring(Cycloheptane) - ring(Ethylene_oxide) + ring(1,4-Cycloheptadiene) + ring(Ethylene_oxide) + radical(Cds_S)"""),
)

species(
    label = '[CH]1C=CC2OC2C=C1(1655)',
    structure = SMILES('[CH]1C=CC2OC2C=C1'),
    E0 = (178.874,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.92868,-0.0169788,0.000213792,-2.64643e-07,9.8775e-11,21589,18.7483], Tmin=(100,'K'), Tmax=(958.549,'K')), NASAPolynomial(coeffs=[19.1557,0.0193334,-5.82029e-06,1.31491e-09,-1.18539e-13,13699,-83.7795], Tmin=(958.549,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(178.874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s2_3_7_ane) - ring(Cycloheptane) - ring(Ethylene_oxide) + ring(1,4-Cycloheptadiene) + ring(Ethylene_oxide) + radical(C=CCJC=C)"""),
)

species(
    label = '[CH]1C2CC=CC23OC13(4331)',
    structure = SMILES('[CH]1C2CC=CC23OC13'),
    E0 = (330.349,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93923,0.0355039,1.92145e-05,-3.7019e-08,1.21313e-11,39813.7,18.1706], Tmin=(100,'K'), Tmax=(1231.58,'K')), NASAPolynomial(coeffs=[11.2525,0.0360003,-1.88355e-05,3.8475e-09,-2.78863e-13,35188.1,-38.1666], Tmin=(1231.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(330.349,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_1) + polycyclic(s1_3_5_ene_1) + polycyclic(s2_3_4_ane) - ring(Cyclobutane) - ring(Ethylene_oxide) - ring(Cyclopentene) + radical(CCJCO)"""),
)

species(
    label = '[CH]1CC=CC23OC2C13(4332)',
    structure = SMILES('[CH]1CC=CC23OC2C13'),
    E0 = (326.608,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.20906,0.0387265,2.38061e-06,-1.70242e-08,5.12221e-12,39346.4,17.6288], Tmin=(100,'K'), Tmax=(1486.39,'K')), NASAPolynomial(coeffs=[12.6703,0.0353182,-1.91504e-05,3.83233e-09,-2.69421e-13,33503.1,-46.1815], Tmin=(1486.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(326.608,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_3_6_ene_1) + polycyclic(s1_3_6_ene_1) + polycyclic(s2_3_3_ane) - ring(Cyclohexene) - ring(Cyclopropane) - ring(Ethylene_oxide) + radical(Cs_S)"""),
)

species(
    label = '[CH]=CCC=CC1=CO1(4333)',
    structure = SMILES('[CH]=CCC=CC1=CO1'),
    E0 = (496.084,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,1000,180,180,180,180,1600,1608.24,2882.55,3200],'cm^-1')),
        HinderedRotor(inertia=(0.144599,'amu*angstrom^2'), symmetry=1, barrier=(3.32462,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144599,'amu*angstrom^2'), symmetry=1, barrier=(3.32462,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144599,'amu*angstrom^2'), symmetry=1, barrier=(3.32462,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.262864,0.0694394,-4.01777e-05,-7.13945e-09,1.03545e-11,59811.1,28.27], Tmin=(100,'K'), Tmax=(963.135,'K')), NASAPolynomial(coeffs=[19.9022,0.0175821,-5.68053e-06,1.00693e-09,-7.26471e-14,54650.2,-72.894], Tmin=(963.135,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(496.084,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclopropene) + radical(Cds_P)"""),
)

species(
    label = '[CH]1OC12C=CCC=C2(4334)',
    structure = SMILES('[CH]1OC12C=CCC=C2'),
    E0 = (225.758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21686,0.0495054,-4.41333e-06,-2.37458e-08,1.05473e-11,27262.5,20.2246], Tmin=(100,'K'), Tmax=(1097.48,'K')), NASAPolynomial(coeffs=[13.0339,0.0313981,-1.37828e-05,2.67058e-09,-1.9123e-13,23165.4,-44.7339], Tmin=(1097.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(225.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s1_3_6_diene_1_4) + radical(CCsJO)"""),
)

species(
    label = 'C1=C2OC2C2=CCC12(4335)',
    structure = SMILES('C1=C2OC2C2=CCC12'),
    E0 = (232.63,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29272,0.0452975,5.87819e-06,-3.62236e-08,1.50405e-11,28088.4,12.2702], Tmin=(100,'K'), Tmax=(1087.19,'K')), NASAPolynomial(coeffs=[15.231,0.0264137,-1.27673e-05,2.61955e-09,-1.94276e-13,23143,-64.9448], Tmin=(1087.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(232.63,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s2_4_5_ane) - ring(Cyclopentane) - ring(Cyclobutane) + ring(Cyclopentene) + ring(Cyclobutene) + Estimated bicyclic component: polycyclic(s2_3_5_ane) - ring(Cyclopentane) - ring(Ethylene_oxide) + ring(Cyclopentene) + ring(Ethylene_oxide) - ring(Cyclopentene)"""),
)

species(
    label = 'C1=CC2C1C=C1OC12(4336)',
    structure = SMILES('C1=CC2C1C=C1OC12'),
    E0 = (236.007,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27633,0.03175,7.35219e-05,-1.29495e-07,5.41792e-11,28508.7,18.6473], Tmin=(100,'K'), Tmax=(955.989,'K')), NASAPolynomial(coeffs=[23.6064,0.00905632,-1.8629e-06,4.76831e-10,-5.07897e-14,21006.8,-104.984], Tmin=(955.989,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(236.007,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_diene_1_5) + Estimated bicyclic component: polycyclic(s2_3_5_ane) - ring(Cyclopentane) - ring(Ethylene_oxide) + ring(Cyclopentene) + ring(Ethylene_oxide) - ring(Cyclopentene)"""),
)

species(
    label = 'C1=C2OC2[C]2CCC12(4337)',
    structure = SMILES('C1=C2OC2[C]2CCC12'),
    E0 = (253.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86212,0.0195992,0.000100847,-1.44399e-07,5.54016e-11,30541.5,16.0099], Tmin=(100,'K'), Tmax=(981.238,'K')), NASAPolynomial(coeffs=[18.1691,0.0222261,-8.80381e-06,1.86894e-09,-1.49901e-13,24014.6,-79.3041], Tmin=(981.238,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(253.092,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_1) + Estimated bicyclic component: polycyclic(s2_3_5_ane) - ring(Cyclopentane) - ring(Ethylene_oxide) + ring(Cyclopentene) + ring(Ethylene_oxide) - ring(Cyclopentene) + radical(CCJ(C)CO)"""),
)

species(
    label = '[CH]1CC2C1C=C1OC12(4338)',
    structure = SMILES('[CH]1CC2C1C=C1OC12'),
    E0 = (288.357,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75729,0.0267855,7.15394e-05,-1.08904e-07,4.16667e-11,34782,17.8466], Tmin=(100,'K'), Tmax=(1000.09,'K')), NASAPolynomial(coeffs=[16.1378,0.0255001,-1.08718e-05,2.25269e-09,-1.73695e-13,29093.5,-65.5918], Tmin=(1000.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(288.357,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_1) + Estimated bicyclic component: polycyclic(s2_3_5_ane) - ring(Cyclopentane) - ring(Ethylene_oxide) + ring(Cyclopentene) + ring(Ethylene_oxide) - ring(Cyclopentene) + radical(cyclobutane)"""),
)

species(
    label = 'C1[C]2CCC2C2OC=12(4339)',
    structure = SMILES('C1[C]2OC2C2CCC=12'),
    E0 = (352.058,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53827,0.0350725,4.49769e-05,-8.25355e-08,3.32408e-11,42448.4,22.7937], Tmin=(100,'K'), Tmax=(989.507,'K')), NASAPolynomial(coeffs=[15.7606,0.0242964,-9.50521e-06,1.88341e-09,-1.42357e-13,37346.7,-57.2273], Tmin=(989.507,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.058,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_0) + polycyclic(s2_3_5_ene_1) - ring(Cyclopentene) + radical(C2CsJO)"""),
)

species(
    label = 'C1=C2O[C]2C2CCC12(4340)',
    structure = SMILES('C1=C2O[C]2C2CCC12'),
    E0 = (280.623,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59329,0.0322412,5.5328e-05,-9.2272e-08,3.59344e-11,33856.1,17.202], Tmin=(100,'K'), Tmax=(1004.38,'K')), NASAPolynomial(coeffs=[16.0666,0.0255918,-1.08952e-05,2.23202e-09,-1.70345e-13,28376.8,-65.4906], Tmin=(1004.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(280.623,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_1) + Estimated bicyclic component: polycyclic(s2_3_5_ane) - ring(Cyclopentane) - ring(Ethylene_oxide) + ring(Cyclopentene) + ring(Ethylene_oxide) - ring(Cyclopentene) + radical(C2CsJO)"""),
)

species(
    label = '[C]1=C2OC2C2CCC12(4341)',
    structure = SMILES('[C]1=C2OC2C2CCC12'),
    E0 = (358.025,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55493,0.0308673,6.47008e-05,-1.05483e-07,4.14033e-11,43168.9,16.6984], Tmin=(100,'K'), Tmax=(992.902,'K')), NASAPolynomial(coeffs=[17.5203,0.0233091,-9.62919e-06,1.99919e-09,-1.55618e-13,37200.6,-74.3008], Tmin=(992.902,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(358.025,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_1) + Estimated bicyclic component: polycyclic(s2_3_5_ane) - ring(Cyclopentane) - ring(Ethylene_oxide) + ring(Cyclopentene) + ring(Ethylene_oxide) - ring(Cyclopentene) + radical(cyclopentene-vinyl)"""),
)

species(
    label = '[CH]1CC2C=C3OC3[C]12(4342)',
    structure = SMILES('[CH]1CC2C=C3OC3[C]12'),
    E0 = (440.931,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03422,0.0200786,8.3062e-05,-1.18807e-07,4.48279e-11,53123,17.7409], Tmin=(100,'K'), Tmax=(997.623,'K')), NASAPolynomial(coeffs=[15.6353,0.0238547,-1.0289e-05,2.16333e-09,-1.68643e-13,47507.6,-62.387], Tmin=(997.623,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(440.931,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_1) + Estimated bicyclic component: polycyclic(s2_3_5_ane) - ring(Cyclopentane) - ring(Ethylene_oxide) + ring(Cyclopentene) + ring(Ethylene_oxide) - ring(Cyclopentene) + radical(CCJ(C)CO) + radical(cyclobutane)"""),
)

species(
    label = '[CH]1C[C]2C=C3OC3C12(4343)',
    structure = SMILES('[CH]1CC2=C[C]3OC3C12'),
    E0 = (546.6,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73096,0.0338044,3.38956e-05,-6.49651e-08,2.57066e-11,65836.6,24.66], Tmin=(100,'K'), Tmax=(1014.4,'K')), NASAPolynomial(coeffs=[13.8795,0.0248737,-1.0529e-05,2.10557e-09,-1.57205e-13,61366.7,-44.0072], Tmin=(1014.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.6,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_0) + polycyclic(s2_3_5_ene_1) - ring(Cyclopentene) + radical(Cs_S) + radical(C2CsJO)"""),
)

species(
    label = '[CH]1CC2C=C3O[C]3C12(4344)',
    structure = SMILES('[CH]1CC2C=C3O[C]3C12'),
    E0 = (468.462,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7626,0.0327304,3.76292e-05,-6.69884e-08,2.55879e-11,56437.7,18.9444], Tmin=(100,'K'), Tmax=(1036.87,'K')), NASAPolynomial(coeffs=[13.6797,0.0269828,-1.22482e-05,2.496e-09,-1.86616e-13,51804.1,-49.4072], Tmin=(1036.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(468.462,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_1) + Estimated bicyclic component: polycyclic(s2_3_5_ane) - ring(Cyclopentane) - ring(Ethylene_oxide) + ring(Cyclopentene) + ring(Ethylene_oxide) - ring(Cyclopentene) + radical(cyclobutane) + radical(C2CsJO)"""),
)

species(
    label = '[CH]1[CH]C2C1C=C1OC12(4345)',
    structure = SMILES('[CH]1[CH]C2C1C=C1OC12'),
    E0 = (476.195,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92859,0.0272569,5.38747e-05,-8.36179e-08,3.12968e-11,57363.6,19.5815], Tmin=(100,'K'), Tmax=(1027.04,'K')), NASAPolynomial(coeffs=[13.7079,0.0269607,-1.22636e-05,2.52561e-09,-1.90693e-13,52540.1,-49.2642], Tmin=(1027.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(476.195,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_1) + Estimated bicyclic component: polycyclic(s2_3_5_ane) - ring(Cyclopentane) - ring(Ethylene_oxide) + ring(Cyclopentene) + ring(Ethylene_oxide) - ring(Cyclopentene) + radical(cyclobutane) + radical(cyclobutane)"""),
)

species(
    label = '[C]1=C2OC2C2[CH]CC12(4346)',
    structure = SMILES('[C]1=C2OC2C2[CH]CC12'),
    E0 = (545.863,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72473,0.031363,4.69162e-05,-7.99859e-08,3.09172e-11,65750.5,18.4383], Tmin=(100,'K'), Tmax=(1017.64,'K')), NASAPolynomial(coeffs=[15.0584,0.0248206,-1.1049e-05,2.27849e-09,-1.73131e-13,60661.7,-57.7914], Tmin=(1017.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(545.863,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_1) + Estimated bicyclic component: polycyclic(s2_3_5_ane) - ring(Cyclopentane) - ring(Ethylene_oxide) + ring(Cyclopentene) + ring(Ethylene_oxide) - ring(Cyclopentene) + radical(cyclopentene-vinyl) + radical(cyclobutane)"""),
)

species(
    label = '[CH]1C2CC3C2C2OC132(4347)',
    structure = SMILES('[CH]1C2CC3C2C2OC132'),
    E0 = (463.174,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.6615,0.0345435,1.08932e-05,-2.28455e-08,6.23186e-12,55750.3,5.96215], Tmin=(100,'K'), Tmax=(1546.56,'K')), NASAPolynomial(coeffs=[13.9872,0.0367938,-2.18829e-05,4.47082e-09,-3.15605e-13,48474.8,-65.8121], Tmin=(1546.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(463.174,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + polycyclic(s3_4_5_ane) + polycyclic(s2_4_4_ane) + polycyclic(s3_4_5_ane) + polycyclic(s1_3_5_ane) + polycyclic(s2_3_4_ane) - ring(Cyclopentane) - ring(Cyclopentane) - ring(Cyclobutane) - ring(Cyclobutane) - ring(Ethylene_oxide) - ring(Cyclobutane) + radical(CCJCO)"""),
)

species(
    label = 'C1C2C3[C]4OC4C2C13(4348)',
    structure = SMILES('C1C2C3[C]4OC4C2C13'),
    E0 = (465.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.19054,0.0347809,1.46315e-05,-2.77523e-08,8.3211e-12,56108.5,7.41289], Tmin=(100,'K'), Tmax=(1343.16,'K')), NASAPolynomial(coeffs=[9.27373,0.0398064,-2.01504e-05,3.98966e-09,-2.81807e-13,51849.7,-37.6198], Tmin=(1343.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(465.941,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + polycyclic(s3_4_4_ane) + polycyclic(s2_4_5_ane) + polycyclic(s3_4_5_ane) + polycyclic(s2_3_5_ane) - ring(Cyclobutane) - ring(Cyclopentane) - ring(Cyclopentane) - ring(Cyclobutane) + radical(C2CsJO)"""),
)

species(
    label = 'C=CC1[CH]C=C2OC21(4349)',
    structure = SMILES('C=CC1C=C[C]2OC21'),
    E0 = (283.289,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23767,0.049635,-8.35212e-06,-1.93928e-08,9.11965e-12,34180.5,22.3111], Tmin=(100,'K'), Tmax=(1101.46,'K')), NASAPolynomial(coeffs=[13.0102,0.0299992,-1.30923e-05,2.53011e-09,-1.80873e-13,30184.9,-41.9877], Tmin=(1101.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(283.289,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + polycyclic(s2_3_5_ene_1) + radical(C2CsJO)"""),
)

species(
    label = '[CH]1OC1=CC1C=CC1(4350)',
    structure = SMILES('[CH]1OC1=CC1C=CC1'),
    E0 = (315.6,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20825,0.0318786,8.33475e-05,-1.4507e-07,6.15887e-11,38085.6,21.6696], Tmin=(100,'K'), Tmax=(936.721,'K')), NASAPolynomial(coeffs=[24.2694,0.00865267,3.97713e-08,-1.89551e-11,-1.20165e-14,30463.8,-105.701], Tmin=(936.721,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(315.6,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene) + ring(methyleneoxirane) + radical(C=CCJO)"""),
)

species(
    label = 'C1=CC2[CH]CC2C=1(4351)',
    structure = SMILES('C1=CC2[CH]CC2C=1'),
    E0 = (738.532,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (91.1305,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.25149,0.0269774,3.50841e-05,-5.28351e-08,1.83626e-11,88897.4,15.483], Tmin=(100,'K'), Tmax=(1087.05,'K')), NASAPolynomial(coeffs=[8.37229,0.0339492,-1.52347e-05,2.98402e-09,-2.14871e-13,85824.1,-22.5727], Tmin=(1087.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(738.532,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + Estimated bicyclic component: polycyclic(s2_4_5_ane) - ring(Cyclopentane) - ring(Cyclobutane) + ring(1,2-Cyclopentadiene) + ring(Cyclobutane) + radical(cyclobutane)"""),
)

species(
    label = '[CH]1C2C=C3OC3C1C2(4352)',
    structure = SMILES('[CH]1C2C=C3OC3C1C2'),
    E0 = (321.174,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03627,0.0177693,9.56256e-05,-1.408e-07,5.60427e-11,38722.2,20.1551], Tmin=(100,'K'), Tmax=(949.562,'K')), NASAPolynomial(coeffs=[17.6411,0.0165111,-4.23911e-06,8.21338e-10,-6.97751e-14,32851.9,-69.6279], Tmin=(949.562,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(321.174,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + polycyclic(s3_4_6_ene_1) + Estimated bicyclic component: polycyclic(s2_3_6_ane) - ring(Cyclohexane) - ring(Ethylene_oxide) + ring(Cyclohexene) + ring(Ethylene_oxide) - ring(Cyclohexene) + radical(cyclobutane)"""),
)

species(
    label = '[CH2]C1C2C=C3OC3C12(4353)',
    structure = SMILES('[CH2]C1C2C=C3OC3C12'),
    E0 = (293.346,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36301,0.0323427,7.06278e-05,-1.21295e-07,4.9962e-11,35399.5,18.9542], Tmin=(100,'K'), Tmax=(959.572,'K')), NASAPolynomial(coeffs=[20.5899,0.0164513,-4.97593e-06,1.01614e-09,-8.51358e-14,28751.4,-88.4254], Tmin=(959.572,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293.346,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + polycyclic(s2_3_5_ene_1) + Estimated bicyclic component: polycyclic(s2_3_5_ane) - ring(Cyclopentane) - ring(Ethylene_oxide) + ring(Cyclopentene) + ring(Ethylene_oxide) - ring(Cyclopentene) + radical(Isobutyl)"""),
)

species(
    label = 'O=C1C2C3[CH]CC3C12(4354)',
    structure = SMILES('O=C1C2C3[CH]CC3C12'),
    E0 = (342.692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93343,0.0373528,9.52083e-06,-2.59422e-08,8.83474e-12,41297.2,19.1938], Tmin=(100,'K'), Tmax=(1208.72,'K')), NASAPolynomial(coeffs=[8.43762,0.0372127,-1.68427e-05,3.2354e-09,-2.27425e-13,38162.7,-19.8799], Tmin=(1208.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(342.692,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-OdCsCs) + polycyclic(s2_4_4_ane) + polycyclic(s2_3_4_ane) - ring(Cyclobutane) + radical(bicyclo[2.2.0]hexane-secondary)"""),
)

species(
    label = '[C]1CC2C=C3OC3C12(4355)',
    structure = SMILES('[C]1CC2C=C3OC3C12'),
    E0 = (548.734,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53832,0.0293022,6.99827e-05,-1.14986e-07,4.59049e-11,66108.3,16.1891], Tmin=(100,'K'), Tmax=(981.283,'K')), NASAPolynomial(coeffs=[19.4821,0.0179366,-7.07936e-06,1.52693e-09,-1.24237e-13,59612.3,-85.1959], Tmin=(981.283,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(548.734,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_1) + Estimated bicyclic component: polycyclic(s2_3_5_ane) - ring(Cyclopentane) - ring(Ethylene_oxide) + ring(Cyclopentene) + ring(Ethylene_oxide) - ring(Cyclopentene) + radical(CCJ2_triplet)"""),
)

species(
    label = 'OC1=[C]C2CC2C=C1(4356)',
    structure = SMILES('OC1=[C]C2CC2C=C1'),
    E0 = (220.006,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11149,0.0453961,2.28444e-05,-6.7119e-08,3.01361e-11,26581,20.1364], Tmin=(100,'K'), Tmax=(962.411,'K')), NASAPolynomial(coeffs=[18.0041,0.0196844,-6.43492e-06,1.20431e-09,-9.123e-14,21268.7,-71.4195], Tmin=(962.411,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(220.006,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_1_3) + radical(Cds_S)"""),
)

species(
    label = 'OC1[C]=CC2CC2C=1(4357)',
    structure = SMILES('OC1[C]=CC2CC2C=1'),
    E0 = (181.159,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06226,0.0473898,1.83676e-05,-6.33037e-08,2.92557e-11,21909.9,19.4361], Tmin=(100,'K'), Tmax=(951.028,'K')), NASAPolynomial(coeffs=[17.665,0.0202513,-6.16402e-06,1.09479e-09,-8.1051e-14,16821.3,-69.9791], Tmin=(951.028,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(181.159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_1_3) + radical(C=CJC=C)"""),
)

species(
    label = 'OC1C=CC2C[C]2C=1(4358)',
    structure = SMILES('OC1C=CC2C[C]2C=1'),
    E0 = (79.0521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05204,0.0518219,-6.07857e-06,-2.76341e-08,1.34286e-11,9625.18,19.3654], Tmin=(100,'K'), Tmax=(1036.46,'K')), NASAPolynomial(coeffs=[14.9244,0.026793,-1.11147e-05,2.14356e-09,-1.55163e-13,5218.3,-55.4443], Tmin=(1036.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(79.0521,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_1_3) + radical(Aromatic_pi_S_1_3)"""),
)

species(
    label = 'OC1C=[C]C2CC2C=1(4359)',
    structure = SMILES('OC1C=[C]C2CC2C=1'),
    E0 = (220.006,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11149,0.0453961,2.28444e-05,-6.7119e-08,3.01361e-11,26581,20.1364], Tmin=(100,'K'), Tmax=(962.411,'K')), NASAPolynomial(coeffs=[18.0041,0.0196844,-6.43492e-06,1.20431e-09,-9.123e-14,21268.7,-71.4195], Tmin=(962.411,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(220.006,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_1_3) + radical(Cds_S)"""),
)

species(
    label = 'OC1C=C[C]2CC2C=1(4360)',
    structure = SMILES('OC1C=C[C]2CC2C=1'),
    E0 = (79.0521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05204,0.0518219,-6.07857e-06,-2.76341e-08,1.34286e-11,9625.18,19.3654], Tmin=(100,'K'), Tmax=(1036.46,'K')), NASAPolynomial(coeffs=[14.9244,0.026793,-1.11147e-05,2.14356e-09,-1.55163e-13,5218.3,-55.4443], Tmin=(1036.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(79.0521,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_1_3) + radical(Aromatic_pi_S_1_3)"""),
)

species(
    label = 'OC1C=CC2[CH]C2C=1(4361)',
    structure = SMILES('OC1C=CC2[CH]C2C=1'),
    E0 = (208.077,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3095,0.0413747,2.94246e-05,-7.0128e-08,3.01857e-11,25138.7,21.2997], Tmin=(100,'K'), Tmax=(971.275,'K')), NASAPolynomial(coeffs=[16.5975,0.0219133,-7.69816e-06,1.46247e-09,-1.09681e-13,20117.2,-62.5727], Tmin=(971.275,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(208.077,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_1_3) + radical(cyclopropane)"""),
)

species(
    label = '[O]C1C=C[C]2CC2C=1(4362)',
    structure = SMILES('[O]C1C=C[C]2CC2C=1'),
    E0 = (216.857,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38878,0.0481031,-1.4027e-05,-1.06903e-08,5.78564e-12,26183.6,19.6637], Tmin=(100,'K'), Tmax=(1135.44,'K')), NASAPolynomial(coeffs=[11.9411,0.0289374,-1.24988e-05,2.38136e-09,-1.68121e-13,22626.4,-37.6977], Tmin=(1135.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.857,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_1_3) + radical(C=C(C)OJ) + radical(Aromatic_pi_S_1_3)"""),
)

species(
    label = '[O]C1C=CC2C[C]2C=1(4363)',
    structure = SMILES('[O]C1C=CC2C[C]2C=1'),
    E0 = (216.857,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3888,0.0481029,-1.40264e-05,-1.06909e-08,5.78588e-12,26183.6,19.6636], Tmin=(100,'K'), Tmax=(1135.43,'K')), NASAPolynomial(coeffs=[11.9411,0.0289375,-1.24989e-05,2.38138e-09,-1.68122e-13,22626.5,-37.6972], Tmin=(1135.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.857,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_1_3) + radical(Aromatic_pi_S_1_3) + radical(C=C(C)OJ)"""),
)

species(
    label = '[O]C1C=CC2[CH]C2C=1(4364)',
    structure = SMILES('[O]C1C=CC2[CH]C2C=1'),
    E0 = (345.882,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66427,0.0375281,2.14983e-05,-5.25522e-08,2.1968e-11,41696.3,21.5279], Tmin=(100,'K'), Tmax=(994.079,'K')), NASAPolynomial(coeffs=[12.9931,0.0250497,-9.62855e-06,1.82495e-09,-1.32711e-13,37808.1,-41.2874], Tmin=(994.079,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(345.882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_1_3) + radical(cyclopropane) + radical(C=C(C)OJ)"""),
)

species(
    label = '[O]C1=[C]C2CC2C=C1(4365)',
    structure = SMILES('[O]C1=[C]C2CC2C=C1'),
    E0 = (357.811,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46675,0.0415493,1.48879e-05,-4.94522e-08,2.18533e-11,43138.5,20.3626], Tmin=(100,'K'), Tmax=(981.37,'K')), NASAPolynomial(coeffs=[14.369,0.0228708,-8.39337e-06,1.57327e-09,-1.1479e-13,38973.2,-49.9604], Tmin=(981.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(357.811,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_1_3) + radical(C=C(C)OJ) + radical(Cds_S)"""),
)

species(
    label = '[O]C1C=[C]C2CC2C=1(4366)',
    structure = SMILES('[O]C1C=[C]C2CC2C=1'),
    E0 = (357.811,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46675,0.0415493,1.48879e-05,-4.94522e-08,2.18533e-11,43138.5,20.3626], Tmin=(100,'K'), Tmax=(981.37,'K')), NASAPolynomial(coeffs=[14.369,0.0228708,-8.39337e-06,1.57327e-09,-1.1479e-13,38973.2,-49.9604], Tmin=(981.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(357.811,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_1_3) + radical(C=C(C)OJ) + radical(Cds_S)"""),
)

species(
    label = '[O]C1[C]=CC2CC2C=1(4367)',
    structure = SMILES('[O]C1[C]=CC2CC2C=1'),
    E0 = (318.964,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41721,0.0435529,1.03431e-05,-4.54917e-08,2.08816e-11,38467.4,19.663], Tmin=(100,'K'), Tmax=(966.975,'K')), NASAPolynomial(coeffs=[14.0027,0.0234822,-8.14739e-06,1.46951e-09,-1.0508e-13,34537.8,-48.3656], Tmin=(966.975,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.964,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_1_3) + radical(C=CJC=C) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C1=CC2CC2C2O[C]12(4368)',
    structure = SMILES('C1=CC2CC2C2O[C]12'),
    E0 = (299.392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79955,0.034411,3.10313e-05,-5.64688e-08,2.11649e-11,36099.7,18.3436], Tmin=(100,'K'), Tmax=(1051.92,'K')), NASAPolynomial(coeffs=[11.3444,0.032059,-1.40159e-05,2.75522e-09,-2.00476e-13,32213.7,-37.1141], Tmin=(1051.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(299.392,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_3_6_ene_1) + polycyclic(s2_3_6_ene_1) - ring(Cyclohexene) + radical(C2CsJO)"""),
)

species(
    label = '[CH]1C2=CC3CC3C1O2(4369)',
    structure = SMILES('[CH]1C2=CC3CC3C1O2'),
    E0 = (398.791,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22415,0.0301506,9.0736e-05,-1.44044e-07,5.69285e-11,48091,12.7322], Tmin=(100,'K'), Tmax=(985.93,'K')), NASAPolynomial(coeffs=[23.2932,0.0176982,-7.59452e-06,1.74493e-09,-1.46611e-13,39992.8,-112.425], Tmin=(985.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(398.791,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + polycyclic(s3_4_6_ene_0) + polycyclic(s2_3_6_ene_1) - ring(3,4-Dihydro-2H-pyran) + radical(C=CCJCO)"""),
)

species(
    label = 'CH2(S)(21)(22)',
    structure = SMILES('[CH2]'),
    E0 = (418.921,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1358.21,2621.43,3089.55],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.19331,-0.00233105,8.15676e-06,-6.62986e-09,1.93233e-12,50366.2,-0.746734], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.13502,0.00289594,-8.16668e-07,1.13573e-10,-6.36263e-15,50504.1,4.06031], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(418.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[O]c1ccccc1(4370)',
    structure = SMILES('[O]c1ccccc1'),
    E0 = (39.8872,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22108,0.0226145,5.07101e-05,-8.26688e-08,3.32557e-11,4876.28,16.2889], Tmin=(100,'K'), Tmax=(962.489,'K')), NASAPolynomial(coeffs=[13.3453,0.0176701,-5.92751e-06,1.1285e-09,-8.61634e-14,822.531,-46.8882], Tmin=(962.489,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(39.8872,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CbH) + group(Cb-O2s) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + ring(Benzene) + radical(RC=COJ)"""),
)

species(
    label = '[O]C1=CCC2CC2=C1(4371)',
    structure = SMILES('[O]C1=CCC2CC2=C1'),
    E0 = (198.117,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16901,0.044305,2.50237e-05,-6.97978e-08,3.14777e-11,23946.2,23.9812], Tmin=(100,'K'), Tmax=(951.634,'K')), NASAPolynomial(coeffs=[17.7696,0.0191666,-5.71334e-06,1.02646e-09,-7.7377e-14,18765.4,-65.904], Tmin=(951.634,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(198.117,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_0_2) + radical(C=C(C)OJ)"""),
)

species(
    label = '[C]1=CC2CC2C=C1(4372)',
    structure = SMILES('[C]1=CC2CC2C=C1'),
    E0 = (395.052,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (91.1305,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.12445,0.0244179,5.74202e-05,-8.86402e-08,3.46591e-11,47596.5,16.3811], Tmin=(100,'K'), Tmax=(970.105,'K')), NASAPolynomial(coeffs=[12.0678,0.0254255,-9.08952e-06,1.70148e-09,-1.25083e-13,43690.6,-41.4757], Tmin=(970.105,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.052,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_1_3) + radical(C=CJC=C)"""),
)

species(
    label = 'O=C1C=CC2CC2=C1(4373)',
    structure = SMILES('O=C1C=CC2CC2=C1'),
    E0 = (201.688,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49891,0.0469449,-1.55302e-05,-6.26513e-09,3.48984e-12,24353.9,21.7667], Tmin=(100,'K'), Tmax=(1279.82,'K')), NASAPolynomial(coeffs=[12.4105,0.0290741,-1.36099e-05,2.64479e-09,-1.86467e-13,20231.6,-38.761], Tmin=(1279.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(201.688,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)(Cds-Cds)) + polycyclic(s2_3_6_diene_0_3)"""),
)

species(
    label = 'O=C1C=CC2C[C]2C1(4374)',
    structure = SMILES('O=C1C=CC2C[C]2C1'),
    E0 = (179.772,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.21303,0.0397065,-4.01653e-06,-9.22743e-09,2.86111e-12,21685.3,19.0272], Tmin=(100,'K'), Tmax=(1570.7,'K')), NASAPolynomial(coeffs=[10.3319,0.0361185,-1.69084e-05,3.17052e-09,-2.14602e-13,17027,-30.5172], Tmin=(1570.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(179.772,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + polycyclic(s2_3_6_ene_1) + radical(Tertalkyl)"""),
)

species(
    label = 'O=C1C=C[C]2CC2C1(4375)',
    structure = SMILES('O=C1C=C[C]2CC2C1'),
    E0 = (126.33,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07228,0.0312266,3.1274e-05,-4.93262e-08,1.69952e-11,15272.7,16.6212], Tmin=(100,'K'), Tmax=(1104.6,'K')), NASAPolynomial(coeffs=[8.12472,0.0380089,-1.69092e-05,3.27591e-09,-2.33717e-13,12184.7,-21.1057], Tmin=(1104.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(126.33,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + polycyclic(s2_3_6_ene_1) + radical(Allyl_T)"""),
)

species(
    label = 'O=C1C=CC2[CH]C2C1(4376)',
    structure = SMILES('O=C1C=CC2[CH]C2C1'),
    E0 = (220.262,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00393,0.0364731,1.11255e-05,-2.62807e-08,8.57985e-12,26569.1,20.2815], Tmin=(100,'K'), Tmax=(1245.95,'K')), NASAPolynomial(coeffs=[8.30319,0.0383678,-1.77832e-05,3.43499e-09,-2.41397e-13,23282.6,-18.3834], Tmin=(1245.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(220.262,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + polycyclic(s2_3_6_ene_1) + radical(cyclopropane)"""),
)

species(
    label = 'O=C1[C]=CC2CC2C1(4377)',
    structure = SMILES('O=C1[C]=CC2CC2C1'),
    E0 = (238.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64809,0.0415136,6.45271e-06,-2.76648e-08,1.02391e-11,28755.1,19.9082], Tmin=(100,'K'), Tmax=(1164.89,'K')), NASAPolynomial(coeffs=[10.7693,0.0345614,-1.59733e-05,3.12726e-09,-2.2322e-13,24976.7,-32.5853], Tmin=(1164.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(238.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + polycyclic(s2_3_6_ene_1) + radical(C=CJC=O)"""),
)

species(
    label = 'O=C1C=[C]C2CC2C1(4378)',
    structure = SMILES('O=C1C=[C]C2CC2C1'),
    E0 = (232.191,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77344,0.0409159,2.88791e-06,-2.08533e-08,7.38064e-12,28012.7,19.2323], Tmin=(100,'K'), Tmax=(1234.66,'K')), NASAPolynomial(coeffs=[9.50672,0.0364411,-1.66775e-05,3.21122e-09,-2.25629e-13,24534.6,-26.0586], Tmin=(1234.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(232.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + polycyclic(s2_3_6_ene_1) + radical(Cds_S)"""),
)

species(
    label = 'O=C1[CH]C2C3CC3C12(4379)',
    structure = SMILES('[O]C1=CC2C3CC3C12'),
    E0 = (161.761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65806,0.0336628,4.48555e-05,-8.11329e-08,3.28661e-11,19555.7,17.7466], Tmin=(100,'K'), Tmax=(976.859,'K')), NASAPolynomial(coeffs=[14.55,0.0248277,-9.07105e-06,1.73106e-09,-1.28731e-13,14939.8,-54.8818], Tmin=(976.859,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(161.761,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + polycyclic(s2_4_4_ene_1) + polycyclic(s2_3_4_ane) - ring(Cyclobutane) + radical(C=C(C)OJ)"""),
)

species(
    label = 'O=C1C2[CH]C3CC3C12(4380)',
    structure = SMILES('O=C1C2[CH]C3CC3C12'),
    E0 = (264.63,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.84553,0.0357654,2.21674e-05,-4.30527e-08,1.54008e-11,31914.9,19.4453], Tmin=(100,'K'), Tmax=(1111.46,'K')), NASAPolynomial(coeffs=[9.9526,0.03534,-1.60599e-05,3.1497e-09,-2.26314e-13,28336.9,-28.5123], Tmin=(1111.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(264.63,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-OdCsCs) + polycyclic(s2_3_5_ane) + polycyclic(s2_3_5_ane) - ring(Cyclopentane) + radical(CCJCC=O)"""),
)

species(
    label = '[CH]=CC1CC1C=C=O(4381)',
    structure = SMILES('[CH]=CC1CC1C=C=O'),
    E0 = (347.065,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,2120,512.5,787.5,3120,650,792.5,1650,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.96196,0.0521027,-9.16749e-07,-3.9835e-08,1.98663e-11,41864.9,26.8422], Tmin=(100,'K'), Tmax=(975.323,'K')), NASAPolynomial(coeffs=[16.8985,0.0216115,-7.64757e-06,1.42014e-09,-1.03896e-13,37097.8,-58.1459], Tmin=(975.323,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(347.065,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(Cds_P)"""),
)

species(
    label = 'O=C1C=CC2[CH]C1C2(4382)',
    structure = SMILES('O=C1C=CC2[CH]C1C2'),
    E0 = (223.56,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.15511,0.0242874,5.82629e-05,-8.27459e-08,2.98877e-11,26968.6,20.0686], Tmin=(100,'K'), Tmax=(1032.24,'K')), NASAPolynomial(coeffs=[10.517,0.0330884,-1.44019e-05,2.85489e-09,-2.09803e-13,23047.1,-31.1713], Tmin=(1032.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(223.56,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + polycyclic(s3_4_6_ene_1) + radical(CCJCC=O)"""),
)

species(
    label = 'O=[C]C1C=CC2CC12(4383)',
    structure = SMILES('O=[C]C1C=CC2CC12'),
    E0 = (190.593,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67008,0.03538,3.51862e-05,-6.69177e-08,2.66117e-11,23020.9,22.0522], Tmin=(100,'K'), Tmax=(1002.33,'K')), NASAPolynomial(coeffs=[13.1742,0.0279051,-1.11456e-05,2.1545e-09,-1.58046e-13,18784.1,-43.1062], Tmin=(1002.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(190.593,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + polycyclic(s2_3_5_ene_1) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'C1=CC2=CC(=CC1)O2(4384)',
    structure = SMILES('C1=CC2=CC(=CC1)O2'),
    E0 = (38.209,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.221975,0.0491464,7.61791e-05,-1.69875e-07,7.64575e-11,4787.48,-5.2613], Tmin=(100,'K'), Tmax=(950.294,'K')), NASAPolynomial(coeffs=[42.2158,-0.0192766,1.02248e-05,-1.56869e-09,7.48407e-14,-8254.36,-234.018], Tmin=(950.294,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(38.209,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Cyclobutene)"""),
)

species(
    label = 'C1=CC2=C[C](CC1)O2(4385)',
    structure = SMILES('[CH]1C=C2C=C(CC1)O2'),
    E0 = (78.0535,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0979257,0.046548,8.51688e-05,-1.7455e-07,7.64438e-11,9574.89,-4.17276], Tmin=(100,'K'), Tmax=(957.886,'K')), NASAPolynomial(coeffs=[40.2704,-0.012717,6.80481e-06,-8.80495e-10,2.55953e-14,-3173.49,-223.366], Tmin=(957.886,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(78.0535,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Cyclobutene) + radical(Allyl_S)"""),
)

species(
    label = '[C]1=C2C=CCCC1O2(4386)',
    structure = SMILES('[C]1=C2C=CCCC1O2'),
    E0 = (151.141,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.612275,0.0520076,9.07715e-05,-1.96399e-07,8.82283e-11,18389.5,-3.18543], Tmin=(100,'K'), Tmax=(947.056,'K')), NASAPolynomial(coeffs=[46.9636,-0.0238431,1.27812e-05,-2.03007e-09,1.03323e-13,3768.26,-259.741], Tmin=(947.056,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(151.141,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Cyclobutene) + radical(Cds_S)"""),
)

species(
    label = '[C]1=CC2=CC(CC1)O2(4387)',
    structure = SMILES('[C]1=CC2=CC(CC1)O2'),
    E0 = (151.141,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.612275,0.0520076,9.07715e-05,-1.96399e-07,8.82283e-11,18389.5,-3.18543], Tmin=(100,'K'), Tmax=(947.056,'K')), NASAPolynomial(coeffs=[46.9636,-0.0238431,1.27812e-05,-2.03007e-09,1.03323e-13,3768.26,-259.741], Tmin=(947.056,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(151.141,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Cyclobutene) + radical(Cds_S)"""),
)

species(
    label = '[C]1=CCCC2C=C1O2(4388)',
    structure = SMILES('[C]1=CCCC2C=C1O2'),
    E0 = (112.295,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.667032,0.0540649,8.60834e-05,-1.92334e-07,8.72557e-11,13718.6,-3.86586], Tmin=(100,'K'), Tmax=(943.233,'K')), NASAPolynomial(coeffs=[46.6618,-0.023339,1.30882e-05,-2.14808e-09,1.14204e-13,-694.924,-258.51], Tmin=(943.233,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(112.295,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Cyclobutene) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]1CC=CC2=C[C]1O2(4389)',
    structure = SMILES('[CH]1[C]2C=CCC=C1O2'),
    E0 = (182.776,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.238012,0.0233268,0.00017645,-2.95068e-07,1.27439e-10,22173.5,-0.583435], Tmin=(100,'K'), Tmax=(925.017,'K')), NASAPolynomial(coeffs=[49.1785,-0.0347165,2.15164e-05,-3.91011e-09,2.37537e-13,6548.4,-268.396], Tmin=(925.017,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(182.776,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Oxetane) + radical(C=CCJC(O)C=C) + radical(C2CsJOC(O))"""),
)

species(
    label = '[CH]1[CH]C2C=C(C=C1)O2(4390)',
    structure = SMILES('[CH]1C=CC2[CH]C(=C1)O2'),
    E0 = (84.285,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.642671,0.00260305,0.000254377,-3.8812e-07,1.63878e-10,10324.4,-3.86267], Tmin=(100,'K'), Tmax=(919.406,'K')), NASAPolynomial(coeffs=[53.5857,-0.0430331,2.74968e-05,-5.10747e-09,3.17563e-13,-7217.25,-297.287], Tmin=(919.406,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.285,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Oxetane) + radical(C=CCJC=C) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = '[C]1=C2C=CC[CH]C1O2(4391)',
    structure = SMILES('[C]1=C2C=CC[CH]C1O2'),
    E0 = (351.043,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.42404,0.0468706,9.91179e-05,-2.05512e-07,9.19085e-11,42426.4,-1.36469], Tmin=(100,'K'), Tmax=(944.341,'K')), NASAPolynomial(coeffs=[47.5055,-0.0281019,1.4816e-05,-2.41369e-09,1.29392e-13,27664.6,-260.082], Tmin=(944.341,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(351.043,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Cyclobutene) + radical(CCJCO) + radical(Cds_S)"""),
)

species(
    label = '[C]1=CC2=CC([CH]C1)O2(4392)',
    structure = SMILES('[C]1=CC2=CC([CH]C1)O2'),
    E0 = (351.043,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.423738,0.0468669,9.91316e-05,-2.05531e-07,9.19169e-11,42426.4,-1.36577], Tmin=(100,'K'), Tmax=(944.326,'K')), NASAPolynomial(coeffs=[47.5045,-0.0281002,1.4815e-05,-2.41345e-09,1.29373e-13,27665,-260.076], Tmin=(944.326,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(351.043,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Cyclobutene) + radical(Cds_S) + radical(CCJCO)"""),
)

species(
    label = '[C]1=CC[CH]C2C=C1O2(4393)',
    structure = SMILES('[C]1=CC[CH]C2C=C1O2'),
    E0 = (312.197,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.478493,0.0489232,9.44522e-05,-2.01487e-07,9.09584e-11,37755.5,-2.04616], Tmin=(100,'K'), Tmax=(940.569,'K')), NASAPolynomial(coeffs=[47.207,-0.0276033,1.5126e-05,-2.5324e-09,1.40331e-13,23199.9,-258.87], Tmin=(940.569,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(312.197,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Cyclobutene) + radical(C=CJC=C) + radical(CCJCO)"""),
)

species(
    label = '[CH]1C2OC13C=CCC23(4394)',
    structure = SMILES('[CH]1C2OC13C=CCC23'),
    E0 = (399.945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22367,0.0288341,3.9623e-05,-5.55053e-08,1.69895e-11,48174.1,15.0825], Tmin=(100,'K'), Tmax=(1265,'K')), NASAPolynomial(coeffs=[13.1084,0.0381588,-2.33028e-05,4.99248e-09,-3.68826e-13,41920.3,-53.8226], Tmin=(1265,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(399.945,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_1) + polycyclic(s1_4_5_ene_1) + polycyclic(s3_4_4_ane) - ring(Cyclopentene) - ring(Oxetane) - ring(Oxetane) + radical(CCJCO)"""),
)

species(
    label = '[CH]1CC2C3C=C(O3)C12(4395)',
    structure = SMILES('[CH]1CC2C3C=C(O3)C12'),
    E0 = (548.953,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6346,0.0307077,5.79737e-05,-9.00098e-08,3.27071e-11,66127.1,20.6915], Tmin=(100,'K'), Tmax=(1068.91,'K')), NASAPolynomial(coeffs=[16.6959,0.0286986,-1.5479e-05,3.37166e-09,-2.58963e-13,59802.2,-67.4989], Tmin=(1068.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(548.953,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ane) + polycyclic(s3_4_5_ene_3) - ring(Tetrahydrofuran) + radical(cyclobutane)"""),
)

species(
    label = 'C1=C[C]2OC3C(C1)C23(4396)',
    structure = SMILES('C1=C[C]2OC3C(C1)C23'),
    E0 = (318.883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.44769,0.0387116,2.40404e-06,-1.65067e-08,4.72114e-12,38404.3,12.9258], Tmin=(100,'K'), Tmax=(1579.3,'K')), NASAPolynomial(coeffs=[14.8169,0.0349367,-2.01807e-05,4.07403e-09,-2.85459e-13,31061.1,-63.28], Tmin=(1579.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.883,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_6_ene_1) + polycyclic(s2_3_6_ene_2) + polycyclic(s2_3_4_ane) - ring(Cyclopropane) - ring(Oxetane) - ring(Cyclohexene) + radical(C2CsJOCs)"""),
)

species(
    label = '[CH]=CC1=CC(C=C)O1(4397)',
    structure = SMILES('[CH]=CC1=CC(C=C)O1'),
    E0 = (356.516,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3120,650,792.5,1650,2750,3150,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.522384,0.0540467,1.8374e-05,-7.6378e-08,3.66446e-11,43024.7,25.2372], Tmin=(100,'K'), Tmax=(948.07,'K')), NASAPolynomial(coeffs=[23.9836,0.0110689,-2.24118e-06,4.29405e-10,-4.01522e-14,36059,-99.9734], Tmin=(948.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(356.516,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1C=CCC=CO1(4398)',
    structure = SMILES('[CH]=C1C=CCC=CO1'),
    E0 = (299.016,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3120,650,792.5,1650,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29351,-0.00162712,0.000241097,-3.62685e-07,1.54718e-10,36118.4,25.1697], Tmin=(100,'K'), Tmax=(900.219,'K')), NASAPolynomial(coeffs=[44.738,-0.0337936,2.66385e-05,-5.35338e-09,3.54904e-13,21777.9,-216.063], Tmin=(900.219,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(299.016,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptane) + radical(Cds_P)"""),
)

species(
    label = '[CH]1CCC=C2C=C1O2(4399)',
    structure = SMILES('[CH]1C2=CCCC=C1O2'),
    E0 = (52.9745,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.148676,0.0373355,0.000120052,-2.19801e-07,9.63542e-11,6553.85,-3.61504], Tmin=(100,'K'), Tmax=(930.542,'K')), NASAPolynomial(coeffs=[41.323,-0.0178748,1.27431e-05,-2.2837e-09,1.31609e-13,-6381.58,-227.624], Tmin=(930.542,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(52.9745,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Oxetane) + radical(C=CCJCO)"""),
)

species(
    label = '[CH]1C=CC=C2CC1O2(4400)',
    structure = SMILES('[CH]1C=CC2CC(=C1)O2'),
    E0 = (14.1416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.39361,0.0167081,0.000204629,-3.28043e-07,1.40477e-10,1888.99,-3.76733], Tmin=(100,'K'), Tmax=(919.712,'K')), NASAPolynomial(coeffs=[49.4824,-0.0343359,2.29273e-05,-4.2799e-09,2.6555e-13,-14011.2,-273.836], Tmin=(919.712,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(14.1416,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Oxetane) + radical(C=CCJC=C)"""),
)

species(
    label = '[CH]1OC2C=CCC1C=2(4401)',
    structure = SMILES('[CH]1OC2C=CCC1C=2'),
    E0 = (173.987,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2934,0.0412393,2.5944e-05,-6.81678e-08,3.00189e-11,21039.7,17.935], Tmin=(100,'K'), Tmax=(970.221,'K')), NASAPolynomial(coeffs=[18.0464,0.0169341,-5.6853e-06,1.11905e-09,-8.76477e-14,15682,-73.2395], Tmin=(970.221,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(173.987,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Cyclohexane) - ring(Tetrahydrofuran) + ring(1,3-Cyclohexadiene) + ring(2,3-Dihydrofuran) + radical(CCsJOC(O))"""),
)

species(
    label = '[CH2]C1C=CC2=CC1O2(4402)',
    structure = SMILES('[CH2]C1C=CC2=CC1O2'),
    E0 = (490.121,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35121,0.0322954,7.27593e-05,-1.26701e-07,5.3195e-11,59067.1,20.2237], Tmin=(100,'K'), Tmax=(943.833,'K')), NASAPolynomial(coeffs=[21.1196,0.0142286,-2.96201e-06,5.50484e-10,-4.98839e-14,52408.6,-89.51], Tmin=(943.833,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(490.121,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s3_4_6_diene_1_5) + radical(Isobutyl)"""),
)

species(
    label = '[CH]1CC2C=C3OC1C32(4403)',
    structure = SMILES('[CH]1CC2C=C3OC1C32'),
    E0 = (359.279,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77139,0.0327589,4.08273e-05,-6.50789e-08,2.23781e-11,43304.9,12.1027], Tmin=(100,'K'), Tmax=(1139.32,'K')), NASAPolynomial(coeffs=[14.2815,0.0333515,-1.85593e-05,3.96434e-09,-2.9706e-13,37565.2,-62.562], Tmin=(1139.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(359.279,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s2_4_5_ane) - ring(Cyclopentane) - ring(Cyclobutane) + ring(Cyclopentane) + ring(Cyclobutene) + polycyclic(s2_4_5_ane) + Estimated bicyclic component: polycyclic(s2_4_4_ane) - ring(Oxetane) - ring(Cyclobutane) + ring(Oxetane) + ring(Cyclobutene) - ring(Cyclobutene) - ring(Oxetane) - ring(Cyclopentane) + radical(CCJCO)"""),
)

species(
    label = '[C]1CC=CC2=CC1O2(4404)',
    structure = SMILES('[C]1CC=CC2=CC1O2'),
    E0 = (361.515,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.633205,0.0504795,9.60015e-05,-2.05971e-07,9.28308e-11,43694.2,-3.67849], Tmin=(100,'K'), Tmax=(943.914,'K')), NASAPolynomial(coeffs=[49.0179,-0.0293683,1.54174e-05,-2.5224e-09,1.36349e-13,28504.8,-271.159], Tmin=(943.914,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.515,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Cyclobutene) + radical(CCJ2_triplet)"""),
)

species(
    label = '[O]C12C=CC3CC1C32(4405)',
    structure = SMILES('[O]C12C=CC3CC1C32'),
    E0 = (319.827,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86648,0.0345671,2.73588e-05,-4.91808e-08,1.77485e-11,38553.6,16.4885], Tmin=(100,'K'), Tmax=(1085.79,'K')), NASAPolynomial(coeffs=[9.90353,0.0353639,-1.57456e-05,3.07475e-09,-2.21168e-13,35016,-31.2005], Tmin=(1085.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(319.827,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_1) + polycyclic(s2_3_5_ene_1) + polycyclic(s2_3_4_ane) - ring(Cyclopropane) - ring(Cyclopentene) - ring(Cyclobutane) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'O=C1C=CC2CC=C12(4406)',
    structure = SMILES('O=C1C=CC2CC=C12'),
    E0 = (145.757,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.11469,0.0403794,-6.96358e-06,-8.05455e-09,2.76743e-12,17599,13.923], Tmin=(100,'K'), Tmax=(1539.23,'K')), NASAPolynomial(coeffs=[11.7018,0.0325541,-1.59912e-05,3.0684e-09,-2.10662e-13,12623.3,-43.0401], Tmin=(1539.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.757,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)(Cds-Cds)) + group(Cd-Cd(CO)H) + Estimated bicyclic component: polycyclic(s2_4_5_ane) - ring(Cyclopentane) - ring(Cyclobutane) + ring(Cyclopentene) + ring(Cyclobutene)"""),
)

species(
    label = 'O=C1C=CC2CC[C]12(4407)',
    structure = SMILES('O=C1C=CC2CC[C]12'),
    E0 = (173.313,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61983,0.0438349,-2.36796e-06,-1.76952e-08,6.94842e-12,20937.3,17.9768], Tmin=(100,'K'), Tmax=(1183.44,'K')), NASAPolynomial(coeffs=[9.54338,0.0356029,-1.54451e-05,2.916e-09,-2.03531e-13,17763,-27.0723], Tmin=(1183.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(173.313,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + polycyclic(s2_4_5_ene_1) + radical(C2CJCHO)"""),
)

species(
    label = 'O=C1C=CC2[CH]CC12(4408)',
    structure = SMILES('O=C1C=CC2[CH]CC12'),
    E0 = (204.324,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0526,0.0344614,1.7332e-05,-3.25495e-08,1.06294e-11,24651.2,20.6837], Tmin=(100,'K'), Tmax=(1209.55,'K')), NASAPolynomial(coeffs=[8.10687,0.0386926,-1.79918e-05,3.49701e-09,-2.47286e-13,21412.5,-17.0102], Tmin=(1209.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.324,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + polycyclic(s2_4_5_ene_1) + radical(cyclobutane)"""),
)

species(
    label = 'O=C1C=C[C]2CCC21(4409)',
    structure = SMILES('O=C1C=C[C]2CCC21'),
    E0 = (148.466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13086,0.0290747,3.80708e-05,-5.64846e-08,1.94675e-11,17933.7,16.9896], Tmin=(100,'K'), Tmax=(1090.4,'K')), NASAPolynomial(coeffs=[8.19527,0.0379349,-1.69091e-05,3.29218e-09,-2.36024e-13,14761.9,-21.2723], Tmin=(1090.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(148.466,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + polycyclic(s2_4_5_ene_1) + radical(Allyl_T)"""),
)

species(
    label = 'O=C1C=[C]C2CCC12(4410)',
    structure = SMILES('O=C1C=[C]C2CCC12'),
    E0 = (273.992,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8272,0.0388454,9.29242e-06,-2.73655e-08,9.52718e-12,33039,19.6161], Tmin=(100,'K'), Tmax=(1195.96,'K')), NASAPolynomial(coeffs=[9.30841,0.036773,-1.68914e-05,3.27474e-09,-2.31655e-13,29608.3,-24.6767], Tmin=(1195.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(273.992,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + polycyclic(s2_4_5_ene_1) + radical(cyclopentene-vinyl)"""),
)

species(
    label = 'O=C1[C]=CC2CCC12(4411)',
    structure = SMILES('O=C1[C]=CC2CCC12'),
    E0 = (260.446,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70669,0.0393736,1.31527e-05,-3.46221e-08,1.25957e-11,31416.1,20.2756], Tmin=(100,'K'), Tmax=(1140.36,'K')), NASAPolynomial(coeffs=[10.725,0.0346643,-1.6068e-05,3.1647e-09,-2.27211e-13,27608.7,-32.0928], Tmin=(1140.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.446,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + polycyclic(s2_4_5_ene_1) + radical(C=CJC=O)"""),
)

species(
    label = 'O=C1C=C[C]2C[CH]C21(4412)',
    structure = SMILES('O=C1C=C[C]2C[CH]C21'),
    E0 = (348.368,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.32586,0.0239039,4.62972e-05,-6.50565e-08,2.27305e-11,41970.2,18.7832], Tmin=(100,'K'), Tmax=(1062.15,'K')), NASAPolynomial(coeffs=[8.41937,0.0341883,-1.51584e-05,2.97373e-09,-2.15242e-13,38801.2,-19.8063], Tmin=(1062.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(348.368,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + polycyclic(s2_4_5_ene_1) + radical(Allyl_T) + radical(CCJCC=O)"""),
)

species(
    label = 'O=C1C=CC2C[CH][C]12(4413)',
    structure = SMILES('O=C1C=CC2C[CH][C]12'),
    E0 = (373.215,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.82248,0.0386125,5.85935e-06,-2.60172e-08,1.00022e-11,44973.5,19.7404], Tmin=(100,'K'), Tmax=(1116.45,'K')), NASAPolynomial(coeffs=[9.40344,0.0324228,-1.40003e-05,2.66632e-09,-1.8824e-13,41973.8,-23.5217], Tmin=(1116.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(373.215,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + polycyclic(s2_4_5_ene_1) + radical(CCJCC=O) + radical(C2CJCHO)"""),
)

species(
    label = 'O=C1C=CC2[CH][CH]C12(4414)',
    structure = SMILES('O=C1C=CC2[CH][CH]C12'),
    E0 = (404.226,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23348,0.0294983,2.46472e-05,-3.96944e-08,1.31896e-11,48688.3,22.5252], Tmin=(100,'K'), Tmax=(1157.94,'K')), NASAPolynomial(coeffs=[7.95205,0.0355206,-1.65452e-05,3.24584e-09,-2.31809e-13,45635.8,-13.3642], Tmin=(1157.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(404.226,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + polycyclic(s2_4_5_ene_1) + radical(cyclobutane) + radical(CCJCC=O)"""),
)

species(
    label = 'O=C1C=[C]C2C[CH]C12(4415)',
    structure = SMILES('O=C1C=[C]C2C[CH]C12'),
    E0 = (473.894,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.01553,0.0337929,1.69253e-05,-3.49238e-08,1.22618e-11,57075.8,21.431], Tmin=(100,'K'), Tmax=(1142.05,'K')), NASAPolynomial(coeffs=[9.17389,0.0335753,-1.54334e-05,3.02142e-09,-2.16032e-13,53819.9,-21.151], Tmin=(1142.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(473.894,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + polycyclic(s2_4_5_ene_1) + radical(cyclopentene-vinyl) + radical(CCJCC=O)"""),
)

species(
    label = 'O=C1[C]=CC2C[CH]C12(4416)',
    structure = SMILES('O=C1[C]=CC2C[CH]C12'),
    E0 = (460.348,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90107,0.0342303,2.11886e-05,-4.28132e-08,1.56409e-11,55452.7,22.0701], Tmin=(100,'K'), Tmax=(1100.89,'K')), NASAPolynomial(coeffs=[10.781,0.0311789,-1.44582e-05,2.87789e-09,-2.08952e-13,51727.3,-29.6641], Tmin=(1100.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(460.348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + polycyclic(s2_4_5_ene_1) + radical(C=CJC=O) + radical(CCJCC=O)"""),
)

species(
    label = 'C1=CC2CC3O[C]1C23(4417)',
    structure = SMILES('[CH]1C=C2OC3CC1C23'),
    E0 = (289.603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52587,0.0321732,5.73761e-05,-9.24961e-08,3.44293e-11,34939.4,15.1024], Tmin=(100,'K'), Tmax=(1052.2,'K')), NASAPolynomial(coeffs=[17.6891,0.026472,-1.39645e-05,3.05527e-09,-2.36579e-13,28452.2,-78.3615], Tmin=(1052.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(289.603,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_0) + polycyclic(s2_4_5_ene_1) + polycyclic(s2_4_4_ane) - ring(Cyclobutane) - ring(Cyclopentene) - ring(Oxetane) + radical(cyclopentene-allyl)"""),
)

species(
    label = 'O=C1C2[CH]C3CC2C13(4418)',
    structure = SMILES('O=C1C2[CH]C3CC2C13'),
    E0 = (311.82,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.16865,0.0346625,1.24693e-05,-2.54021e-08,7.79098e-12,37573.4,12.4738], Tmin=(100,'K'), Tmax=(1311.45,'K')), NASAPolynomial(coeffs=[8.06415,0.0393249,-1.87629e-05,3.64028e-09,-2.55088e-13,34079.8,-24.9912], Tmin=(1311.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(311.82,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-OdCsCs) + polycyclic(s3_4_5_ane) + polycyclic(s2_4_5_ane) + polycyclic(s2_4_4_ane) - ring(Cyclopentane) - ring(Cyclobutane) - ring(Cyclobutane) + radical(CCJCC=O)"""),
)

species(
    label = 'O=C1[CH]C2C3CC2C13(4419)',
    structure = SMILES('[O]C1=CC2C3CC2C13'),
    E0 = (264.139,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98153,0.0284879,4.89044e-05,-7.60394e-08,2.84686e-11,31855.1,12.4799], Tmin=(100,'K'), Tmax=(1017.36,'K')), NASAPolynomial(coeffs=[11.3359,0.0310761,-1.29546e-05,2.5313e-09,-1.85443e-13,27914.5,-42.8231], Tmin=(1017.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(264.139,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_1) + polycyclic(s2_4_5_ene_1) + polycyclic(s3_4_4_ane) - ring(Cyclobutane) - ring(Cyclobutane) - ring(Cyclopentene) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=CC1[CH]C=CC1=O(4420)',
    structure = SMILES('C=CC1C=CC=C1[O]'),
    E0 = (121.682,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09535,0.0483568,1.00172e-05,-5.18799e-08,2.45662e-11,14753.6,23.2584], Tmin=(100,'K'), Tmax=(958.057,'K')), NASAPolynomial(coeffs=[16.7103,0.0211431,-6.83996e-06,1.2289e-09,-8.97681e-14,10018.6,-60.5025], Tmin=(958.057,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(121.682,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclopentadiene) + radical(C=C(C)OJ)"""),
)

species(
    label = 'O=[C]C=CC1C=CC1(4421)',
    structure = SMILES('O=C=C[CH]C1C=CC1'),
    E0 = (255.078,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,2750,2850,2950,3050,3150,900,950,1000,1050,1100,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19668,0.0484088,2.01991e-06,-3.54973e-08,1.62623e-11,30791.3,27.1057], Tmin=(100,'K'), Tmax=(1015.95,'K')), NASAPolynomial(coeffs=[14.266,0.0270245,-1.08078e-05,2.05608e-09,-1.4836e-13,26583.8,-43.7912], Tmin=(1015.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(255.078,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + ring(Cyclobutene) + radical(CCJC(C)=C=O)"""),
)

species(
    label = '[CH2]C1C2C=CC(=O)C12(4422)',
    structure = SMILES('[CH2]C1C2C=CC(=O)C12'),
    E0 = (209.312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61058,0.0408179,1.24688e-05,-3.80731e-08,1.52129e-11,25270.5,21.9473], Tmin=(100,'K'), Tmax=(1053.69,'K')), NASAPolynomial(coeffs=[11.0134,0.0320574,-1.34026e-05,2.55481e-09,-1.82266e-13,21793.8,-31.0025], Tmin=(1053.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(209.312,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + polycyclic(s2_3_5_ene_1) + radical(Isobutyl)"""),
)

species(
    label = 'OC1C=CC2C[CH]C=12(4423)',
    structure = SMILES('OC1C=CC2C[CH]C=12'),
    E0 = (156.27,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.978528,0.0508988,4.74992e-06,-4.31556e-08,1.98022e-11,18917.4,14.0252], Tmin=(100,'K'), Tmax=(1008.01,'K')), NASAPolynomial(coeffs=[16.5497,0.0249866,-1.00788e-05,1.96089e-09,-1.44459e-13,13955.5,-70.262], Tmin=(1008.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(156.27,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s2_4_5_ane) - ring(Cyclopentane) - ring(Cyclobutane) + ring(Cyclopentadiene) + ring(Cyclobutane) + radical(Allyl_S)"""),
)

species(
    label = 'OC1=C=CC2C[CH]C12(4424)',
    structure = SMILES('OC1=C=CC2C[CH]C12'),
    E0 = (523.973,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23352,0.0460156,1.08811e-05,-4.5031e-08,1.94336e-11,63132,19.1875], Tmin=(100,'K'), Tmax=(1021.56,'K')), NASAPolynomial(coeffs=[14.9315,0.0266765,-1.1081e-05,2.16482e-09,-1.58756e-13,58543.7,-55.9487], Tmin=(1021.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(523.973,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + Estimated bicyclic component: polycyclic(s2_4_5_ane) - ring(Cyclopentane) - ring(Cyclobutane) + ring(1,2-Cyclopentadiene) + ring(Cyclobutane) + radical(cyclobutane)"""),
)

species(
    label = 'O=C1C=CC2C[C]C12(4425)',
    structure = SMILES('O=C1C=CC2C[C]C12'),
    E0 = (464.7,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79884,0.0375163,1.32822e-05,-3.45093e-08,1.27258e-11,55978.8,19.1425], Tmin=(100,'K'), Tmax=(1123.58,'K')), NASAPolynomial(coeffs=[10.5322,0.0325488,-1.49615e-05,2.94183e-09,-2.114e-13,52367.3,-31.3466], Tmin=(1123.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(464.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + polycyclic(s2_4_5_ene_1) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C1=CCC2=CC(=C1)O2(4426)',
    structure = SMILES('C1=CCC2=CC(=C1)O2'),
    E0 = (37.5434,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.162163,0.0449296,9.25295e-05,-1.90177e-07,8.45098e-11,4708,-4.66165], Tmin=(100,'K'), Tmax=(947.216,'K')), NASAPolynomial(coeffs=[43.6748,-0.0221693,1.18906e-05,-1.88106e-09,9.49911e-14,-8891.12,-241.72], Tmin=(947.216,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(37.5434,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Cyclobutene)"""),
)

species(
    label = 'C1=CC[C]2CC(=C1)O2(4427)',
    structure = SMILES('C1=CC[C]2CC(=C1)O2'),
    E0 = (94.3598,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.337429,0.0452031,0.000110796,-2.23682e-07,1.01615e-10,11551.7,-2.8092], Tmin=(100,'K'), Tmax=(921.546,'K')), NASAPolynomial(coeffs=[46.7131,-0.0286076,1.86631e-05,-3.46734e-09,2.14908e-13,-2657.78,-256.002], Tmin=(921.546,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(94.3598,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Oxetane) + radical(C2CsJOC(O))"""),
)

species(
    label = '[C]1C=CCC2CC=1O2(4428)',
    structure = SMILES('[C]1C=CCC2CC=1O2'),
    E0 = (93.1396,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.325897,0.0433615,0.000122219,-2.39492e-07,1.08645e-10,11406.2,-4.75078], Tmin=(100,'K'), Tmax=(915.742,'K')), NASAPolynomial(coeffs=[47.4172,-0.0298033,2.03108e-05,-3.86454e-09,2.45035e-13,-3014.18,-261.874], Tmin=(915.742,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(93.1396,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Oxetane) + radical(C=CJC=C)"""),
)

species(
    label = '[C]1=CC=C2CC(C1)O2(4429)',
    structure = SMILES('[C]1=CC=C2CC(C1)O2'),
    E0 = (131.986,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.273324,0.0413419,0.000126709,-2.43182e-07,1.09394e-10,16077.1,-4.0632], Tmin=(100,'K'), Tmax=(919.516,'K')), NASAPolynomial(coeffs=[47.6804,-0.0302435,1.99678e-05,-3.73815e-09,2.33467e-13,1465.77,-262.886], Tmin=(919.516,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(131.986,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Oxetane) + radical(Cds_S)"""),
)

species(
    label = '[C]1=CCC2CC(=C1)O2(4430)',
    structure = SMILES('[C]1=CCC2CC(=C1)O2'),
    E0 = (93.1396,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.325897,0.0433615,0.000122219,-2.39492e-07,1.08645e-10,11406.2,-4.75078], Tmin=(100,'K'), Tmax=(915.742,'K')), NASAPolynomial(coeffs=[47.4172,-0.0298033,2.03108e-05,-3.86454e-09,2.45035e-13,-3014.18,-261.874], Tmin=(915.742,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(93.1396,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Oxetane) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]1[C]2CC=CC=C1O2(4431)',
    structure = SMILES('[CH]1C=CCC2=C[C]1O2'),
    E0 = (201.602,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.076349,0.0328095,0.000145214,-2.54698e-07,1.09122e-10,24446,0.112354], Tmin=(100,'K'), Tmax=(947.526,'K')), NASAPolynomial(coeffs=[48.8613,-0.0294067,1.5151e-05,-2.37904e-09,1.20181e-13,8691.01,-267.545], Tmin=(947.526,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(201.602,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Cyclobutene) + radical(C=CCJC(O)C=C) + radical(C2CsJOC(O))"""),
)

species(
    label = '[C]1=CC=C2[CH]C(C1)O2(4432)',
    structure = SMILES('[C]1=CC=C2[CH]C(C1)O2'),
    E0 = (248.903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.490145,0.0452759,0.000119054,-2.39501e-07,1.09092e-10,30147.6,-6.08569], Tmin=(100,'K'), Tmax=(919.404,'K')), NASAPolynomial(coeffs=[49.8083,-0.0342598,2.1557e-05,-4.01746e-09,2.51797e-13,15011.3,-276.537], Tmin=(919.404,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(248.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Oxetane) + radical(Cds_S) + radical(C=CCJCO)"""),
)

species(
    label = '[C]1=CCC2[CH]C(=C1)O2(4433)',
    structure = SMILES('[C]1=CCC2[CH]C(=C1)O2'),
    E0 = (210.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.543017,0.0472992,0.00011455,-2.35791e-07,1.08334e-10,25476.6,-6.77221], Tmin=(100,'K'), Tmax=(915.704,'K')), NASAPolynomial(coeffs=[49.5461,-0.0338213,2.1901e-05,-4.14409e-09,2.63385e-13,10531,-275.531], Tmin=(915.704,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(210.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Oxetane) + radical(C=CCJCO) + radical(C=CJC=C)"""),
)

species(
    label = '[C]1C=CCC2[CH]C=1O2(4434)',
    structure = SMILES('[C]1C=CCC2[CH]C=1O2'),
    E0 = (210.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.543017,0.0472992,0.00011455,-2.35791e-07,1.08334e-10,25476.6,-6.77221], Tmin=(100,'K'), Tmax=(915.704,'K')), NASAPolynomial(coeffs=[49.5461,-0.0338213,2.1901e-05,-4.14409e-09,2.63385e-13,10531,-275.531], Tmin=(915.704,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(210.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Oxetane) + radical(C=CCJCO) + radical(C=CJC=C)"""),
)

species(
    label = 'C1=CC2[C]3OC(C1)C32(4435)',
    structure = SMILES('C1=CC2[C]3OC(C1)C32'),
    E0 = (334.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.53565,0.0392229,-9.46481e-08,-1.3898e-08,3.96074e-12,40250.9,12.7758], Tmin=(100,'K'), Tmax=(1636.01,'K')), NASAPolynomial(coeffs=[15.5154,0.0339553,-1.95322e-05,3.91151e-09,-2.71869e-13,32461.9,-67.0589], Tmin=(1636.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(334.28,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_6_ene_2) + polycyclic(s2_3_6_ene_1) + polycyclic(s2_3_4_ane) - ring(Cyclohexene) - ring(Oxetane) - ring(Cyclopropane) + radical(C2CsJOCs)"""),
)

species(
    label = '[CH]1CC2OC3=CC1C32(4436)',
    structure = SMILES('[CH]1CC2OC3=CC1C32'),
    E0 = (353.919,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85561,0.0354875,2.63342e-05,-4.61639e-08,1.50648e-11,42653.1,11.8753], Tmin=(100,'K'), Tmax=(1228.01,'K')), NASAPolynomial(coeffs=[13.339,0.0359172,-2.04058e-05,4.30003e-09,-3.16561e-13,36980,-57.4995], Tmin=(1228.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(353.919,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ane) + Estimated bicyclic component: polycyclic(s2_4_5_ane) - ring(Cyclopentane) - ring(Cyclobutane) + ring(Cyclopentane) + ring(Cyclobutene) + Estimated bicyclic component: polycyclic(s2_4_4_ane) - ring(Oxetane) - ring(Cyclobutane) + ring(Oxetane) + ring(Cyclobutene) - ring(Cyclobutene) - ring(Cyclopentane) - ring(Oxetane) + radical(Cs_S)"""),
)

species(
    label = '[O]C1C=C=CC=CC1(4437)',
    structure = SMILES('[O]C1C=C=CC=CC1'),
    E0 = (353.304,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37045,0.0382118,4.26113e-05,-8.40748e-08,3.49125e-11,42604.9,19.6956], Tmin=(100,'K'), Tmax=(974.12,'K')), NASAPolynomial(coeffs=[16.806,0.0230589,-8.32194e-06,1.60908e-09,-1.21759e-13,37309.5,-66.1124], Tmin=(974.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(353.304,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C=CC=C1C=CO1(4438)',
    structure = SMILES('[CH2]C=CC=C1C=CO1'),
    E0 = (231.055,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2750,3150,900,1100,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.891168,0.04544,4.00777e-05,-9.75749e-08,4.45136e-11,27922.6,24.4762], Tmin=(100,'K'), Tmax=(929.987,'K')), NASAPolynomial(coeffs=[22.1826,0.0125768,-1.61835e-06,2.0285e-10,-2.07167e-14,21423.4,-90.3483], Tmin=(929.987,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(231.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]1C2C=CC=CC1O2(4439)',
    structure = SMILES('[CH]1C2C=CC=CC1O2'),
    E0 = (137.523,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.258109,0.0270945,0.00016074,-2.72393e-07,1.1781e-10,16726.2,-1.51368], Tmin=(100,'K'), Tmax=(925.704,'K')), NASAPolynomial(coeffs=[45.8999,-0.0275247,1.81764e-05,-3.31373e-09,1.99406e-13,2166.14,-251.187], Tmin=(925.704,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(137.523,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Oxetane) + radical(CCJCO)"""),
)

species(
    label = '[CH]1OC2=CC=CCC12(4440)',
    structure = SMILES('[CH]1OC2=CC=CCC12'),
    E0 = (253.652,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30069,0.0348652,6.0823e-05,-1.16842e-07,5.11078e-11,30627.1,18.6057], Tmin=(100,'K'), Tmax=(928.354,'K')), NASAPolynomial(coeffs=[21.7644,0.00962997,-9.36395e-08,-7.05561e-11,-3.29598e-15,24115.5,-93.2054], Tmin=(928.354,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(253.652,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s2_4_6_ane) - ring(Cyclohexane) - ring(Oxetane) + ring(1,3-Cyclohexadiene) + ring(Oxetane) + radical(CCsJOC(O))"""),
)

species(
    label = '[C]1C2=CC=CCC1O2(4441)',
    structure = SMILES('[C]1=C2[CH]C=CCC1O2'),
    E0 = (286.001,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.483125,0.0470495,0.000103491,-2.10115e-07,9.30486e-11,34606.6,-3.05106], Tmin=(100,'K'), Tmax=(948.94,'K')), NASAPolynomial(coeffs=[47.8677,-0.026387,1.34901e-05,-2.10513e-09,1.05354e-13,19560.2,-264.711], Tmin=(948.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(286.001,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Cyclobutene) + radical(C=CCJCO) + radical(Cds_S)"""),
)

species(
    label = 'C1=CCC2C=C(C=1)O2(4442)',
    structure = SMILES('C1=CCC2C=C(C=1)O2'),
    E0 = (128.524,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.355148,0.081737,-7.20494e-05,2.65797e-08,-3.78436e-12,15586.3,-13.6979], Tmin=(100,'K'), Tmax=(1647.92,'K')), NASAPolynomial(coeffs=[26.4848,0.0183128,-1.43185e-05,3.22477e-09,-2.41288e-13,6974.32,-152.81], Tmin=(1647.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(128.524,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(1_2_cycloheptadiene) + ring(Cyclobutene)"""),
)

species(
    label = '[C]1=CCC2C=C(C1)O2(4443)',
    structure = SMILES('[C]1=CCC2C=C(C1)O2'),
    E0 = (169.084,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.266681,0.0431191,0.000111139,-2.13797e-07,9.3357e-11,20536.2,-1.02717], Tmin=(100,'K'), Tmax=(949.595,'K')), NASAPolynomial(coeffs=[45.746,-0.0223811,1.19069e-05,-1.82721e-09,8.71384e-14,6012.03,-251.095], Tmin=(949.595,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(169.084,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Cyclobutene) + radical(Cds_S)"""),
)

species(
    label = '[C]1=C2CC=CCC1O2(4444)',
    structure = SMILES('[C]1=C2CC=CCC1O2'),
    E0 = (169.084,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.266681,0.0431191,0.000111139,-2.13797e-07,9.3357e-11,20536.2,-1.02717], Tmin=(100,'K'), Tmax=(949.595,'K')), NASAPolynomial(coeffs=[45.746,-0.0223811,1.19069e-05,-1.82721e-09,8.71384e-14,6012.03,-251.095], Tmin=(949.595,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(169.084,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Cyclobutene) + radical(Cds_S)"""),
)

species(
    label = '[C]1=CCC2=CC(C1)O2(4445)',
    structure = SMILES('[C]1=CCC2=CC(C1)O2'),
    E0 = (169.084,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.266681,0.0431191,0.000111139,-2.13797e-07,9.3357e-11,20536.2,-1.02717], Tmin=(100,'K'), Tmax=(949.595,'K')), NASAPolynomial(coeffs=[45.746,-0.0223811,1.19069e-05,-1.82721e-09,8.71384e-14,6012.03,-251.095], Tmin=(949.595,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(169.084,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Cyclobutene) + radical(Cds_S)"""),
)

species(
    label = 'C1=CCC2=C[C](C1)O2(4446)',
    structure = SMILES('C1=CCC2=C[C](C1)O2'),
    E0 = (131.458,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.330737,0.0469707,9.53099e-05,-1.94499e-07,8.57138e-11,16010.8,-0.465973], Tmin=(100,'K'), Tmax=(953.57,'K')), NASAPolynomial(coeffs=[44.8182,-0.0208109,1.06394e-05,-1.56508e-09,6.92921e-14,1871.38,-245.128], Tmin=(953.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(131.458,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Cyclobutene) + radical(C2CsJOC(O))"""),
)

species(
    label = '[CH]1C=CCC2=CC1O2(4447)',
    structure = SMILES('[CH]1C=CCC2=CC1O2'),
    E0 = (1.38616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0149111,0.0244223,0.000185635,-3.03073e-07,1.27684e-10,368.266,-1.71317], Tmin=(100,'K'), Tmax=(943.634,'K')), NASAPolynomial(coeffs=[51.0088,-0.0305278,1.67251e-05,-2.69696e-09,1.39914e-13,-16433,-282.801], Tmin=(943.634,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1.38616,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Cyclobutene) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = '[CH]1C2CC3C=C(O3)C12(4448)',
    structure = SMILES('[CH]1C2CC3C=C(O3)C12'),
    E0 = (548.734,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91065,0.0222901,7.97672e-05,-1.1255e-07,4.07675e-11,66093.1,19.6039], Tmin=(100,'K'), Tmax=(1044.76,'K')), NASAPolynomial(coeffs=[16.5741,0.026954,-1.42278e-05,3.13398e-09,-2.44069e-13,59710.6,-67.6624], Tmin=(1044.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(548.734,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + polycyclic(s3_4_6_ene_4) + polycyclic(s2_3_6_ane) - ring(Oxane) + radical(cyclopropane)"""),
)

species(
    label = 'C=CC1[CH]C12C=CO2(4449)',
    structure = SMILES('C=CC1[CH]C12C=CO2'),
    E0 = (357.553,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.994541,0.038129,6.7868e-05,-1.33453e-07,5.93011e-11,43137.8,24.4096], Tmin=(100,'K'), Tmax=(920.931,'K')), NASAPolynomial(coeffs=[25.0998,0.00589767,2.33067e-06,-5.64383e-10,3.0808e-14,35624.9,-106.583], Tmin=(920.931,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(357.553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + polycyclic(s1_3_4_ene) + radical(CCJCO)"""),
)

species(
    label = '[CH]1C2CC3C=CC13O2(4450)',
    structure = SMILES('[CH]1C2CC3C=CC13O2'),
    E0 = (398.546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.1309,0.0321528,3.00528e-05,-4.66707e-08,1.43292e-11,48008.1,16.895], Tmin=(100,'K'), Tmax=(1290.5,'K')), NASAPolynomial(coeffs=[13.658,0.0369129,-2.25426e-05,4.8124e-09,-3.54165e-13,41661.4,-54.7194], Tmin=(1290.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(398.546,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s2_4_5_ane) - ring(Tetrahydrofuran) - ring(Cyclobutane) + ring(Tetrahydrofuran) + ring(Cyclobutene) + polycyclic(s1_4_4_ene_1) + polycyclic(s3_4_5_ane) - ring(Oxetane) - ring(Tetrahydrofuran) - ring(Cyclobutene) + radical(CCJCO)"""),
)

species(
    label = '[CH]1CC=C2C=C(C1)O2(4451)',
    structure = SMILES('[CH]1CC=C2C=C(C1)O2'),
    E0 = (131.399,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0223834,0.051995,5.32577e-05,-1.31379e-07,5.8882e-11,15981.5,-0.772606], Tmin=(100,'K'), Tmax=(966.931,'K')), NASAPolynomial(coeffs=[35.8676,-0.00611213,3.21993e-06,-2.30266e-10,-1.50317e-14,4816.63,-194.558], Tmin=(966.931,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(131.399,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(Cycloheptane) + ring(Cyclobutene) + radical(RCCJCC)"""),
)

species(
    label = 'O=C1C=C2CC=CC12(4452)',
    structure = SMILES('O=C1C=C2CC=CC12'),
    E0 = (150.832,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.01097,0.0429087,-1.25934e-05,-3.99838e-09,1.79745e-12,18213,14.8066], Tmin=(100,'K'), Tmax=(1599.21,'K')), NASAPolynomial(coeffs=[13.1698,0.0301143,-1.47715e-05,2.82038e-09,-1.92526e-13,12710.9,-50.3111], Tmin=(1599.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(150.832,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + Estimated bicyclic component: polycyclic(s2_4_5_ane) - ring(Cyclopentane) - ring(Cyclobutane) + ring(Cyclopentene) + ring(Cyclobutene)"""),
)

species(
    label = 'O=C1C[C]2CC=CC21(4453)',
    structure = SMILES('O=C1C[C]2CC=CC21'),
    E0 = (195.785,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23769,0.0423363,-1.13333e-05,-3.11678e-09,1.28875e-12,23607.8,18.8543], Tmin=(100,'K'), Tmax=(1760.74,'K')), NASAPolynomial(coeffs=[12.876,0.0322069,-1.46634e-05,2.67227e-09,-1.76141e-13,17685.4,-44.6674], Tmin=(1760.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(195.785,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_1) + radical(Tertalkyl)"""),
)

species(
    label = 'O=C1CC2CC=C[C]12(4454)',
    structure = SMILES('O=C1CC2CC=C[C]12'),
    E0 = (130.063,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.14549,0.0243881,6.07139e-05,-8.56527e-08,3.10376e-11,15724,17.4809], Tmin=(100,'K'), Tmax=(1023.89,'K')), NASAPolynomial(coeffs=[10.1562,0.0344564,-1.46343e-05,2.86371e-09,-2.09185e-13,11915.4,-31.9432], Tmin=(1023.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(130.063,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_1) + radical(C=CCJ(C)C=O)"""),
)

species(
    label = 'O=C1CC2[CH]C=CC12(4455)',
    structure = SMILES('O=C1CC2[CH]C=CC12'),
    E0 = (137.668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86595,0.035702,2.21006e-05,-4.269e-08,1.53551e-11,16643.9,17.2318], Tmin=(100,'K'), Tmax=(1099.07,'K')), NASAPolynomial(coeffs=[9.28117,0.0361389,-1.59238e-05,3.07741e-09,-2.19529e-13,13357.5,-26.7781], Tmin=(1099.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(137.668,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_1) + radical(cyclopentene-allyl)"""),
)

species(
    label = 'O=C1CC2CC=[C]C12(4456)',
    structure = SMILES('O=C1CC2CC=[C]C12'),
    E0 = (267.869,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63435,0.0453051,-9.80804e-06,-8.81122e-09,3.69987e-12,32307.7,19.6586], Tmin=(100,'K'), Tmax=(1318.98,'K')), NASAPolynomial(coeffs=[10.1506,0.0348864,-1.5482e-05,2.9132e-09,-2.01058e-13,28720.9,-28.8657], Tmin=(1318.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(267.869,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_1) + radical(cyclopentene-vinyl)"""),
)

species(
    label = 'O=C1CC2C[C]=CC12(4457)',
    structure = SMILES('O=C1CC2C[C]=CC12'),
    E0 = (267.869,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63435,0.0453051,-9.80804e-06,-8.81122e-09,3.69987e-12,32307.7,19.6586], Tmin=(100,'K'), Tmax=(1318.98,'K')), NASAPolynomial(coeffs=[10.1506,0.0348864,-1.5482e-05,2.9132e-09,-2.01058e-13,28720.9,-28.8657], Tmin=(1318.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(267.869,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_1) + radical(cyclopentene-vinyl)"""),
)

species(
    label = 'O=C1[CH][C]2CC=CC21(4458)',
    structure = SMILES('[O]C1=C[C]2CC=CC21'),
    E0 = (315.528,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55239,0.0309534,6.24842e-05,-1.12701e-07,4.81768e-11,38058.4,21.479], Tmin=(100,'K'), Tmax=(933.51,'K')), NASAPolynomial(coeffs=[19.4599,0.0126188,-1.88981e-06,2.8396e-10,-2.76126e-14,32170.6,-77.3107], Tmin=(933.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(315.528,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_diene_1_5) + radical(C=C(C)OJ) + radical(Allyl_T)"""),
)

species(
    label = 'O=C1[CH]C2CC=C[C]12(4459)',
    structure = SMILES('O=C1[CH]C2CC=C[C]12'),
    E0 = (292.192,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.25545,0.0171364,8.59235e-05,-1.17853e-07,4.39655e-11,35224.3,17.4691], Tmin=(100,'K'), Tmax=(990.264,'K')), NASAPolynomial(coeffs=[12.9655,0.0275483,-1.11496e-05,2.23232e-09,-1.69088e-13,30471.4,-47.3837], Tmin=(990.264,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(292.192,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_1) + radical(C=CCJ(C)C=O) + radical(CCJC=O)"""),
)

species(
    label = 'O=C1[CH]C2[CH]C=CC12(4460)',
    structure = SMILES('O=C1[CH]C2[CH]C=CC12'),
    E0 = (299.797,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98359,0.0284176,4.71324e-05,-7.42038e-08,2.77783e-11,36143.7,17.1888], Tmin=(100,'K'), Tmax=(1022.69,'K')), NASAPolynomial(coeffs=[11.6782,0.0298901,-1.28026e-05,2.52903e-09,-1.86143e-13,32100.9,-39.8704], Tmin=(1022.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(299.797,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_1) + radical(cyclopentene-allyl) + radical(CCJC=O)"""),
)

species(
    label = 'O=C1[CH]C2CC=[C]C12(4461)',
    structure = SMILES('O=C1[CH]C2CC=[C]C12'),
    E0 = (429.998,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73385,0.0384031,1.31262e-05,-3.65842e-08,1.41276e-11,51808.1,19.6691], Tmin=(100,'K'), Tmax=(1081.72,'K')), NASAPolynomial(coeffs=[10.87,0.0311893,-1.37148e-05,2.66518e-09,-1.91404e-13,48277,-32.3106], Tmin=(1081.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(429.998,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_1) + radical(CCJC=O) + radical(cyclopentene-vinyl)"""),
)

species(
    label = 'O=C1[CH]C2C[C]=CC12(4462)',
    structure = SMILES('O=C1[CH]C2C[C]=CC12'),
    E0 = (429.998,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73385,0.0384031,1.31262e-05,-3.65842e-08,1.41276e-11,51808.1,19.6691], Tmin=(100,'K'), Tmax=(1081.72,'K')), NASAPolynomial(coeffs=[10.87,0.0311893,-1.37148e-05,2.66518e-09,-1.91404e-13,48277,-32.3106], Tmin=(1081.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(429.998,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_1) + radical(CCJC=O) + radical(cyclopentene-vinyl)"""),
)

species(
    label = 'C1=CC2[C]3OC3C2C1(4463)',
    structure = SMILES('C1=CC2[C]3OC3C2C1'),
    E0 = (321.783,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8952,0.0322454,3.62235e-05,-6.04044e-08,2.21678e-11,38789.4,17.1952], Tmin=(100,'K'), Tmax=(1054.61,'K')), NASAPolynomial(coeffs=[10.7022,0.03344,-1.46858e-05,2.88562e-09,-2.09663e-13,35007.8,-34.8837], Tmin=(1054.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(321.783,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_1) + polycyclic(s2_3_4_ane) - ring(Cyclobutane) + radical(C2CsJO)"""),
)

species(
    label = 'O=C1C2C3[CH]CC2C13(4464)',
    structure = SMILES('O=C1C2C3[CH]CC2C13'),
    E0 = (342.115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.65845,0.0335544,7.51982e-06,-1.67678e-08,4.44137e-12,41191.4,9.31092], Tmin=(100,'K'), Tmax=(1572.67,'K')), NASAPolynomial(coeffs=[9.17099,0.0393888,-1.94087e-05,3.70375e-09,-2.52539e-13,36373.1,-33.8634], Tmin=(1572.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(342.115,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-OdCsCs) + polycyclic(s2_4_5_ane) + polycyclic(s2_4_5_ane) + polycyclic(s3_4_4_ane) - ring(Cyclobutane) - ring(Cyclobutane) - ring(Cyclopentane) + radical(bicyclo[2.1.1]hexane-C2)"""),
)

species(
    label = 'O=C=CC1[CH]C=CC1(2038)',
    structure = SMILES('O=C=CC1C=C[CH]C1'),
    E0 = (134.014,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3528,0.0419954,2.42955e-05,-6.13029e-08,2.61083e-11,16227.9,21.8078], Tmin=(100,'K'), Tmax=(984.87,'K')), NASAPolynomial(coeffs=[15.1443,0.0250521,-9.40455e-06,1.78877e-09,-1.31554e-13,11616.5,-54.1377], Tmin=(984.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(134.014,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCd(CCO)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + ring(Cyclopentene) + radical(cyclopentene-allyl)"""),
)

species(
    label = '[CH2]C=CC1C=CC1=O(4465)',
    structure = SMILES('[CH2]C=CC1C=CC1=O'),
    E0 = (220.524,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2950,3150,900,1000,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09775,0.056451,-3.26444e-05,8.02893e-09,-6.26153e-13,26633.6,25.0553], Tmin=(100,'K'), Tmax=(1613.93,'K')), NASAPolynomial(coeffs=[16.0277,0.0266713,-1.16799e-05,2.14207e-09,-1.43801e-13,20873.7,-57.0337], Tmin=(1613.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(220.524,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + ring(Cyclobutene) + radical(Allyl_P)"""),
)

species(
    label = 'O=C1C2[CH]C1CC=C2(4466)',
    structure = SMILES('O=C1C2[CH]C1CC=C2'),
    E0 = (220.719,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91857,0.031791,3.66932e-05,-6.19577e-08,2.32629e-11,26633.4,20.403], Tmin=(100,'K'), Tmax=(1030.79,'K')), NASAPolynomial(coeffs=[10.6404,0.032092,-1.3434e-05,2.59846e-09,-1.8815e-13,23021.2,-30.7386], Tmin=(1030.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(220.719,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s3_4_6_ene_1) + radical(CCJCC=O)"""),
)

species(
    label = 'O=[C]C1C2C=CCC12(2089)',
    structure = SMILES('O=[C]C1C2C=CCC12'),
    E0 = (191.728,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,1855,455,950,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73379,0.0316856,4.93712e-05,-8.42054e-08,3.33705e-11,23157.3,22.4907], Tmin=(100,'K'), Tmax=(985.853,'K')), NASAPolynomial(coeffs=[14.1768,0.025933,-9.93909e-06,1.92866e-09,-1.4375e-13,18530.1,-48.3873], Tmin=(985.853,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.728,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + polycyclic(s2_3_5_ene_1) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[CH]1C2=CC1CC=CO2(4467)',
    structure = SMILES('[CH]1C2=CC1CC=CO2'),
    E0 = (157.517,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.133602,0.0361433,0.000125204,-2.30939e-07,1.02495e-10,19129.6,-0.231974], Tmin=(100,'K'), Tmax=(921.696,'K')), NASAPolynomial(coeffs=[43.2446,-0.0240735,1.67185e-05,-3.11991e-09,1.92073e-13,5793.32,-233.938], Tmin=(921.696,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.517,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Cyclobutane) + ring(Cycloheptane) + ring(Cyclobutene) + radical(CCJCO)"""),
)

species(
    label = 'OC1=CC2CC=C[C]12(2393)',
    structure = SMILES('OC1=CC2CC=C[C]12'),
    E0 = (198.316,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18391,0.0336335,7.60116e-05,-1.36896e-07,5.87147e-11,23979.4,23.2422], Tmin=(100,'K'), Tmax=(934.673,'K')), NASAPolynomial(coeffs=[23.8968,0.00868226,1.04172e-07,-5.11037e-11,-8.46898e-15,16577.6,-101.683], Tmin=(934.673,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(198.316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_diene_1_5) + radical(CCJ(C)CO)"""),
)

species(
    label = 'O=C1[C]C2CC=CC12(4468)',
    structure = SMILES('[O]C1=[C]C2CC=CC12'),
    E0 = (436.033,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24465,0.0409259,3.22434e-05,-8.05343e-08,3.64148e-11,52559.9,24.1103], Tmin=(100,'K'), Tmax=(940.717,'K')), NASAPolynomial(coeffs=[19.4461,0.0132291,-2.837e-06,4.84828e-10,-4.07636e-14,46936.4,-74.2767], Tmin=(940.717,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(436.033,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_diene_1_5) + radical(cyclobutene-vinyl) + radical(C=C(C)OJ)"""),
)

species(
    label = 'OC1=[C]C2CC=CC12(2395)',
    structure = SMILES('OC1=[C]C2CC=CC12'),
    E0 = (298.228,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.883822,0.0448279,4.00657e-05,-9.81347e-08,4.47271e-11,36002.6,23.9047], Tmin=(100,'K'), Tmax=(934.913,'K')), NASAPolynomial(coeffs=[23.1587,0.00991369,-8.05357e-07,9.87773e-11,-1.57994e-14,29198.4,-96.1739], Tmin=(934.913,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(298.228,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_diene_1_5) + radical(cyclobutene-vinyl)"""),
)

species(
    label = 'OC1=C[C]2CC=CC21(2392)',
    structure = SMILES('OC1=C[C]2CC=CC21'),
    E0 = (177.723,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19062,0.0348657,7.02758e-05,-1.30272e-07,5.64822e-11,21501.2,21.2768], Tmin=(100,'K'), Tmax=(929.668,'K')), NASAPolynomial(coeffs=[23.181,0.00928888,1.5015e-07,-1.04047e-10,-2.48669e-15,14428.9,-99.2564], Tmin=(929.668,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(177.723,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_diene_1_5) + radical(Allyl_T)"""),
)

species(
    label = 'OC1=CC2CC=[C]C12(2396)',
    structure = SMILES('OC1=CC2CC=[C]C12'),
    E0 = (303.249,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.883822,0.0448279,4.00657e-05,-9.81347e-08,4.47271e-11,36606.4,23.9047], Tmin=(100,'K'), Tmax=(934.913,'K')), NASAPolynomial(coeffs=[23.1587,0.00991369,-8.05357e-07,9.87773e-11,-1.57994e-14,29802.3,-96.1739], Tmin=(934.913,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.249,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_diene_1_5) + radical(cyclopentene-vinyl)"""),
)

species(
    label = 'OC1=CC2[CH]C=CC12(2394)',
    structure = SMILES('OC1=CC2[CH]C=CC12'),
    E0 = (173.048,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11337,0.0350185,7.37828e-05,-1.35917e-07,5.87102e-11,20943,21.5007], Tmin=(100,'K'), Tmax=(933.59,'K')), NASAPolynomial(coeffs=[24.4334,0.0078571,5.29128e-07,-1.34644e-10,-2.62154e-15,13418.2,-106.383], Tmin=(933.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(173.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_diene_1_5) + radical(cyclopentene-allyl)"""),
)

species(
    label = 'OC1=CC2C[C]=CC12(2397)',
    structure = SMILES('OC1=CC2C[C]=CC12'),
    E0 = (303.249,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.883822,0.0448279,4.00657e-05,-9.81347e-08,4.47271e-11,36606.4,23.9047], Tmin=(100,'K'), Tmax=(934.913,'K')), NASAPolynomial(coeffs=[23.1587,0.00991369,-8.05357e-07,9.87773e-11,-1.57994e-14,29802.3,-96.1739], Tmin=(934.913,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.249,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_diene_1_5) + radical(cyclopentene-vinyl)"""),
)

species(
    label = '[CH]1C2CC3C=C(O2)C13(4469)',
    structure = SMILES('[CH]1C2CC3C=C(O2)C13'),
    E0 = (289.21,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94201,0.0268066,5.81119e-05,-8.29017e-08,2.85461e-11,34873.4,11.8567], Tmin=(100,'K'), Tmax=(1106.52,'K')), NASAPolynomial(coeffs=[14.046,0.0330646,-1.81689e-05,3.90397e-09,-2.94951e-13,29133,-61.5981], Tmin=(1106.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(289.21,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + polycyclic(s3_5_5_ane) + Estimated bicyclic component: polycyclic(s2_4_5_ane) - ring(Tetrahydrofuran) - ring(Cyclobutane) + ring(Tetrahydrofuran) + ring(Cyclobutene) + Estimated bicyclic component: polycyclic(s2_4_5_ane) - ring(Cyclopentane) - ring(Cyclobutane) + ring(Cyclopentane) + ring(Cyclobutene) - ring(Tetrahydrofuran) - ring(Cyclopentane) - ring(Cyclobutene) + radical(CCJCO)"""),
)

species(
    label = 'C=CC1C2C=CC12[O](4470)',
    structure = SMILES('C=CC1C2C=CC12[O]'),
    E0 = (429.456,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44604,0.0366876,4.40831e-05,-8.4929e-08,3.51623e-11,51761.1,25.0556], Tmin=(100,'K'), Tmax=(972.495,'K')), NASAPolynomial(coeffs=[16.5752,0.0224645,-8.02314e-06,1.54977e-09,-1.17518e-13,46548.4,-59.1829], Tmin=(972.495,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(429.456,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + polycyclic(s2_3_4_ene_1) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[C]1=CC2CC=CC12(1432)',
    structure = SMILES('[C]1=CC2CC=CC12'),
    E0 = (512.787,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (91.1305,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88434,0.0261063,6.25886e-05,-1.02829e-07,4.18762e-11,61768.7,20.256], Tmin=(100,'K'), Tmax=(954.582,'K')), NASAPolynomial(coeffs=[16.0611,0.0180484,-5.43461e-06,1.02663e-09,-8.06988e-14,56722.7,-59.7334], Tmin=(954.582,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(512.787,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_diene_1_5) + radical(cyclobutene-vinyl)"""),
)

species(
    label = 'OC1C=CC=C=CC=1(4471)',
    structure = SMILES('OC1C=CC=C=CC=1'),
    E0 = (181.954,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04777,0.0466037,1.73562e-05,-6.51118e-08,3.0763e-11,22007,17.3286], Tmin=(100,'K'), Tmax=(947.207,'K')), NASAPolynomial(coeffs=[19.5769,0.0141057,-3.629e-06,6.49555e-10,-5.19939e-14,16444.5,-81.8925], Tmin=(947.207,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(181.954,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + ring(1,2,4,6-Cycloheptatetraene)"""),
)

species(
    label = 'OC1C=C=CC=CC=1(4472)',
    structure = SMILES('OC1C=C=CC=CC=1'),
    E0 = (181.954,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04777,0.0466037,1.73562e-05,-6.51118e-08,3.0763e-11,22007,17.3286], Tmin=(100,'K'), Tmax=(947.207,'K')), NASAPolynomial(coeffs=[19.5769,0.0141057,-3.629e-06,6.49555e-10,-5.19939e-14,16444.5,-81.8925], Tmin=(947.207,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(181.954,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + ring(1,2,4,6-Cycloheptatetraene)"""),
)

species(
    label = 'OH(5)(5)',
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
    label = '[C]1=CC=C[CH]C=C1(4473)',
    structure = SMILES('[C]1=C[CH]C=CC=C1'),
    E0 = (466.224,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.1225,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22643,0.0410983,3.34001e-05,-8.60067e-08,4.01216e-11,56192.1,15.5109], Tmin=(100,'K'), Tmax=(918.237,'K')), NASAPolynomial(coeffs=[20.2704,0.00993996,-3.20528e-07,-8.81543e-11,2.59622e-15,50510.9,-86.6323], Tmin=(918.237,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(466.224,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=CJC=C) + radical(C=CCJC=C)"""),
)

species(
    label = 'OC1[C]=C[CH]C=CC=1(4474)',
    structure = SMILES('OC1[C]=C[CH]C=CC=1'),
    E0 = (252.331,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.171491,0.0639513,-5.05339e-06,-6.17698e-08,3.53592e-11,30505.2,19.235], Tmin=(100,'K'), Tmax=(904.941,'K')), NASAPolynomial(coeffs=[25.9676,0.00460168,2.6971e-06,-7.16182e-10,4.83718e-14,23597.8,-115.009], Tmin=(904.941,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(252.331,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=CJC=C) + radical(C=CCJC=C)"""),
)

species(
    label = 'OC1C=C[CH][C]=CC=1(4475)',
    structure = SMILES('OC1C=C[C]=C[CH]C=1'),
    E0 = (252.331,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.171491,0.0639513,-5.05339e-06,-6.17698e-08,3.53592e-11,30505.2,19.235], Tmin=(100,'K'), Tmax=(904.941,'K')), NASAPolynomial(coeffs=[25.9676,0.00460168,2.6971e-06,-7.16182e-10,4.83718e-14,23597.8,-115.009], Tmin=(904.941,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(252.331,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=CCJC=C) + radical(C=CJC=C)"""),
)

species(
    label = 'OC1C=[C][CH]C=CC=1(4476)',
    structure = SMILES('OC1C=[C]C=C[CH]C=1'),
    E0 = (252.331,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.171491,0.0639513,-5.05339e-06,-6.17698e-08,3.53592e-11,30505.2,19.235], Tmin=(100,'K'), Tmax=(904.941,'K')), NASAPolynomial(coeffs=[25.9676,0.00460168,2.6971e-06,-7.16182e-10,4.83718e-14,23597.8,-115.009], Tmin=(904.941,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(252.331,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=CJC=C) + radical(C=CCJC=C)"""),
)

species(
    label = 'OC1=CC=CC2[CH]C12(4477)',
    structure = SMILES('OC1=CC=CC2[CH]C12'),
    E0 = (207.411,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37495,0.0370944,4.59791e-05,-9.06581e-08,3.83153e-11,25059,21.879], Tmin=(100,'K'), Tmax=(960.085,'K')), NASAPolynomial(coeffs=[18.0126,0.0190944,-6.07453e-06,1.15999e-09,-9.03471e-14,19499.2,-70.0275], Tmin=(960.085,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(207.411,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_1_3) + radical(cyclopropane)"""),
)

species(
    label = 'OC1=C[CH]C2C=CC12(2462)',
    structure = SMILES('OC1=C[CH]C2C=CC12'),
    E0 = (173.048,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11337,0.0350185,7.37828e-05,-1.35917e-07,5.87102e-11,20943,21.5007], Tmin=(100,'K'), Tmax=(933.59,'K')), NASAPolynomial(coeffs=[24.4334,0.0078571,5.29128e-07,-1.34644e-10,-2.62154e-15,13418.2,-106.383], Tmin=(933.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(173.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_diene_1_5) + radical(cyclopentene-allyl)"""),
)

species(
    label = 'OC12[CH]C=CC1C=C2(2353)',
    structure = SMILES('OC12[CH]C=CC1C=C2'),
    E0 = (140.57,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28336,0.0278935,9.86174e-05,-1.63423e-07,6.87909e-11,17033.8,24.5339], Tmin=(100,'K'), Tmax=(933.148,'K')), NASAPolynomial(coeffs=[25.1701,0.00699287,1.22008e-06,-2.53926e-10,3.60324e-15,9027.85,-108.064], Tmin=(933.148,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(140.57,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_diene_1_5) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = 'OC1C=[C]C=CCC=1(2296)',
    structure = SMILES('OC1C=[C]C=CCC=1'),
    E0 = (150.606,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.35309,0.0902546,-9.0107e-05,4.30324e-08,-7.62278e-12,18330.1,25.0571], Tmin=(100,'K'), Tmax=(1618.78,'K')), NASAPolynomial(coeffs=[24.4866,0.00996358,-4.72605e-07,-1.56214e-10,1.61272e-14,12118.5,-105.397], Tmin=(1618.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(150.606,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=CJC=C)"""),
)

species(
    label = 'OC1=C=CC=CC=C1(4478)',
    structure = SMILES('OC1=C=CC=CC=C1'),
    E0 = (181.954,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04777,0.0466037,1.73562e-05,-6.51118e-08,3.0763e-11,22007,17.3286], Tmin=(100,'K'), Tmax=(947.207,'K')), NASAPolynomial(coeffs=[19.5769,0.0141057,-3.629e-06,6.49555e-10,-5.19939e-14,16444.5,-81.8925], Tmin=(947.207,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(181.954,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + ring(1,2,4,6-Cycloheptatetraene)"""),
)

species(
    label = 'OC1=[C]CC=CC=C1(2292)',
    structure = SMILES('OC1=[C]CC=CC=C1'),
    E0 = (189.452,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.10964,0.0745168,-3.89797e-05,-1.95504e-08,1.76841e-11,22948.4,21.4484], Tmin=(100,'K'), Tmax=(933.064,'K')), NASAPolynomial(coeffs=[24.0504,0.0112204,-1.97267e-06,2.70734e-10,-2.19104e-14,16686.6,-102.831], Tmin=(933.064,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(189.452,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(Cds_S)"""),
)

species(
    label = 'OC1C=CC=[C]CC=1(2293)',
    structure = SMILES('OC1C=CC=[C]CC=1'),
    E0 = (189.452,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.10964,0.0745168,-3.89797e-05,-1.95504e-08,1.76841e-11,22948.4,21.4484], Tmin=(100,'K'), Tmax=(933.064,'K')), NASAPolynomial(coeffs=[24.0504,0.0112204,-1.97267e-06,2.70734e-10,-2.19104e-14,16686.6,-102.831], Tmin=(933.064,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(189.452,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(Cds_S)"""),
)

species(
    label = 'OC1C=C[C]=CCC=1(2294)',
    structure = SMILES('OC1C=C[C]=CCC=1'),
    E0 = (150.606,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.35308,0.0902545,-9.01067e-05,4.30321e-08,-7.62271e-12,18330.1,25.057], Tmin=(100,'K'), Tmax=(1618.79,'K')), NASAPolynomial(coeffs=[24.4862,0.00996411,-4.72873e-07,-1.56157e-10,1.61228e-14,12118.7,-105.395], Tmin=(1618.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(150.606,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=CJC=C)"""),
)

species(
    label = '[O]C1C=CC=CCC=1(2291)',
    structure = SMILES('[O]C1C=CC=CCC=1'),
    E0 = (89.4153,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.277283,0.0660947,-2.22768e-05,-3.0918e-08,2.02334e-11,10902.7,21.0867], Tmin=(100,'K'), Tmax=(939.65,'K')), NASAPolynomial(coeffs=[21.553,0.0149923,-3.70167e-06,6.01893e-10,-4.50598e-14,5162,-89.5032], Tmin=(939.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(89.4153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=C(C)OJ)"""),
)

species(
    label = 'OC1[C]=CC=CCC=1(2295)',
    structure = SMILES('OC1[C]=CC=CCC=1'),
    E0 = (150.606,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.35309,0.0902546,-9.0107e-05,4.30324e-08,-7.62278e-12,18330.1,25.0571], Tmin=(100,'K'), Tmax=(1618.78,'K')), NASAPolynomial(coeffs=[24.4866,0.00996358,-4.72605e-07,-1.56214e-10,1.61272e-14,12118.5,-105.397], Tmin=(1618.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(150.606,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=CJC=C)"""),
)

species(
    label = 'OC12[CH]C1C=CC=C2(4479)',
    structure = SMILES('OC12[CH]C1C=CC=C2'),
    E0 = (206.084,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51209,0.0331859,5.64529e-05,-1.01911e-07,4.26169e-11,24895.4,23.2741], Tmin=(100,'K'), Tmax=(951.835,'K')), NASAPolynomial(coeffs=[17.9141,0.0183077,-5.27774e-06,9.83015e-10,-7.76073e-14,19324.6,-67.9079], Tmin=(951.835,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(206.084,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_1_3) + radical(CCJCO)"""),
)

species(
    label = 'OC1[CH]C2C=CC2C=1(2316)',
    structure = SMILES('OC1[CH]C2C=CC2C=1'),
    E0 = (245.644,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10176,0.0351226,7.31363e-05,-1.36525e-07,5.94265e-11,29675,24.4523], Tmin=(100,'K'), Tmax=(930.112,'K')), NASAPolynomial(coeffs=[24.9207,0.00610401,1.53575e-06,-3.40541e-10,1.21473e-14,22068.5,-105.806], Tmin=(930.112,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(245.644,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_diene_1_5) + radical(CCJCO)"""),
)

species(
    label = 'OC1=CC=[C]C=CC1(2268)',
    structure = SMILES('OC1=CC=[C]C=CC1'),
    E0 = (149.94,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0922731,0.072192,-2.66325e-05,-3.68551e-08,2.53129e-11,18197.6,21.3247], Tmin=(100,'K'), Tmax=(916.539,'K')), NASAPolynomial(coeffs=[25.2251,0.00880568,1.36001e-08,-1.62582e-10,9.34765e-15,11578.2,-109.404], Tmin=(916.539,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(149.94,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=CJC=C)"""),
)

species(
    label = 'C1=CC=CC=CC=1(4480)',
    structure = SMILES('C1=CC=CC=CC=1'),
    E0 = (395.847,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (90.1225,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.11091,0.0236229,5.64271e-05,-9.0449e-08,3.61545e-11,47693.6,13.5769], Tmin=(100,'K'), Tmax=(964.686,'K')), NASAPolynomial(coeffs=[13.9642,0.0193058,-6.56917e-06,1.25967e-09,-9.63071e-14,43320.6,-53.9941], Tmin=(964.686,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.847,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + ring(1,2,4,6-Cycloheptatetraene)"""),
)

species(
    label = 'OC1=CC=CC=[C]C1(2265)',
    structure = SMILES('OC1=CC=CC=[C]C1'),
    E0 = (188.787,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0472182,0.0702639,-2.24744e-05,-4.01014e-08,2.58684e-11,22868.8,22.0391], Tmin=(100,'K'), Tmax=(928.142,'K')), NASAPolynomial(coeffs=[25.5183,0.00831362,-2.99171e-07,-4.33906e-11,-1.61967e-15,16045.8,-110.584], Tmin=(928.142,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.787,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(Cds_S)"""),
)

species(
    label = 'OC1=[C]C=CC=CC1(2266)',
    structure = SMILES('OC1=[C]C=CC=CC1'),
    E0 = (149.94,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0922731,0.072192,-2.66325e-05,-3.68551e-08,2.53129e-11,18197.6,21.3247], Tmin=(100,'K'), Tmax=(916.539,'K')), NASAPolynomial(coeffs=[25.2251,0.00880568,1.36001e-08,-1.62582e-10,9.34765e-15,11578.2,-109.404], Tmin=(916.539,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(149.94,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=CJC=C)"""),
)

species(
    label = 'OC1=CC=C[C]=CC1(2267)',
    structure = SMILES('OC1=CC=C[C]=CC1'),
    E0 = (149.94,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0922731,0.072192,-2.66325e-05,-3.68551e-08,2.53129e-11,18197.6,21.3247], Tmin=(100,'K'), Tmax=(916.539,'K')), NASAPolynomial(coeffs=[25.2251,0.00880568,1.36001e-08,-1.62582e-10,9.34765e-15,11578.2,-109.404], Tmin=(916.539,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(149.94,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=CJC=C)"""),
)

species(
    label = 'OC1=C[C]=CC=CC1(2269)',
    structure = SMILES('OC1=C[C]=CC=CC1'),
    E0 = (149.94,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0922731,0.072192,-2.66325e-05,-3.68551e-08,2.53129e-11,18197.6,21.3247], Tmin=(100,'K'), Tmax=(916.539,'K')), NASAPolynomial(coeffs=[25.2251,0.00880568,1.36001e-08,-1.62582e-10,9.34765e-15,11578.2,-109.404], Tmin=(916.539,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(149.94,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=CJC=C)"""),
)

species(
    label = 'O[C]1C2C=CC=CC12(4481)',
    structure = SMILES('O[C]1C2C=CC=CC12'),
    E0 = (200.214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31013,0.0442905,1.46292e-05,-5.00994e-08,2.19022e-11,24190.2,20.7971], Tmin=(100,'K'), Tmax=(993.609,'K')), NASAPolynomial(coeffs=[14.8005,0.0251755,-9.64373e-06,1.83439e-09,-1.34015e-13,19772.2,-52.9417], Tmin=(993.609,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(200.214,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_1_3) + radical(C2CsJOH)"""),
)

species(
    label = 'OC1C=C[C]=CC=C1(2246)',
    structure = SMILES('OC1C=C[C]=CC=C1'),
    E0 = (186.666,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.190747,0.0707046,-3.95272e-05,-9.3074e-09,1.14796e-11,22599.7,22.6339], Tmin=(100,'K'), Tmax=(956.488,'K')), NASAPolynomial(coeffs=[20.1867,0.0180949,-5.65793e-06,9.84399e-10,-7.05914e-14,17355.9,-80.3612], Tmin=(956.488,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(186.666,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=CJC=C)"""),
)

species(
    label = 'OC1[C]=CC=CC=C1(2244)',
    structure = SMILES('OC1[C]=CC=CC=C1'),
    E0 = (225.512,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.23542,0.0687635,-3.52268e-05,-1.29103e-08,1.22791e-11,27271,23.3506], Tmin=(100,'K'), Tmax=(972.803,'K')), NASAPolynomial(coeffs=[20.5549,0.017479,-5.90077e-06,1.08733e-09,-8.02257e-14,21790.9,-81.9657], Tmin=(972.803,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(225.512,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(Cds_S)"""),
)

species(
    label = 'OC1C=[C]C=CC=C1(2245)',
    structure = SMILES('OC1C=[C]C=CC=C1'),
    E0 = (186.666,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.190747,0.0707046,-3.95272e-05,-9.3074e-09,1.14796e-11,22599.7,22.6339], Tmin=(100,'K'), Tmax=(956.488,'K')), NASAPolynomial(coeffs=[20.1867,0.0180949,-5.65793e-06,9.84399e-10,-7.05914e-14,17355.9,-80.3612], Tmin=(956.488,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(186.666,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=CJC=C)"""),
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
    label = 'Ar(8)',
    structure = SMILES('[Ar]'),
    E0 = (-6.19426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (39.348,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1134.93,'J/mol'), sigma=(3.33,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,-745,4.3663], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,-745,4.3663], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-6.19426,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ar""", comment="""Thermo library: BurkeH2O2"""),
)

transitionState(
    label = 'TS1',
    E0 = (432.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (291.538,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (402.932,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (539.049,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (539.049,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (504.399,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (504.399,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (500.661,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (263.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (394.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (292.583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (343.684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (233.761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (254.211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (151.414,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (222.492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (303.078,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (220.273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (281.224,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (185.766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (215.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (607.504,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (366.275,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (380.919,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (379.089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (201.289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (275.859,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (334.736,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (220.273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (194.691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (424.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (279.608,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (293.871,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (382.75,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (366.939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (220.273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (119.969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (316.758,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (323.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (325.779,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (481.197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (285.49,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (423.83,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (407.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (600.323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (310.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (534.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (416.996,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (619.056,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (619.056,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (619.056,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (351.426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (300.279,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (310.144,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (422.127,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (796.913,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (488.883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (483.222,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (688.316,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (575.403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (497.573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (454.737,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (441.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (334.975,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS65',
    E0 = (332.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS66',
    E0 = (537.129,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS67',
    E0 = (383.077,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS68',
    E0 = (453.501,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS69',
    E0 = (453.531,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS70',
    E0 = (442.747,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS71',
    E0 = (447.349,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS72',
    E0 = (526.726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS73',
    E0 = (430.724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS74',
    E0 = (402.334,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS75',
    E0 = (652.723,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS76',
    E0 = (758.911,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS77',
    E0 = (680.773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS78',
    E0 = (687.987,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS79',
    E0 = (757.656,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS80',
    E0 = (464.626,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS81',
    E0 = (467.762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS82',
    E0 = (408.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS83',
    E0 = (382.544,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS84',
    E0 = (1170.95,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS85',
    E0 = (467.145,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS86',
    E0 = (481.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS87',
    E0 = (453.281,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS88',
    E0 = (602.157,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS89',
    E0 = (760.526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS90',
    E0 = (370.351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS91',
    E0 = (323.136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS92',
    E0 = (220.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS93',
    E0 = (264.314,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS94',
    E0 = (208.251,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS95',
    E0 = (267.373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS96',
    E0 = (429.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS97',
    E0 = (429.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS98',
    E0 = (557.674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS99',
    E0 = (570.061,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS100',
    E0 = (569.603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS101',
    E0 = (534.953,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS102',
    E0 = (333.631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS103',
    E0 = (398.791,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS104',
    E0 = (307.21,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS105',
    E0 = (458.808,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS106',
    E0 = (383.397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS107',
    E0 = (638.057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS108',
    E0 = (420.513,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS109',
    E0 = (301.402,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS110',
    E0 = (253.105,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS111',
    E0 = (347.874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS112',
    E0 = (379.729,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS113',
    E0 = (276.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS114',
    E0 = (250.826,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS115',
    E0 = (375.193,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS116',
    E0 = (363.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS117',
    E0 = (199.288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS118',
    E0 = (376.322,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS119',
    E0 = (383.495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS120',
    E0 = (333.574,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS121',
    E0 = (259.456,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS122',
    E0 = (241.825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS123',
    E0 = (241.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS124',
    E0 = (210.638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS125',
    E0 = (292.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS126',
    E0 = (292.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS127',
    E0 = (156.604,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS128',
    E0 = (394.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS129',
    E0 = (296.077,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS130',
    E0 = (562.836,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS131',
    E0 = (562.836,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS132',
    E0 = (528.186,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS133',
    E0 = (399.945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS134',
    E0 = (548.953,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS135',
    E0 = (318.883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS136',
    E0 = (398.791,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS137',
    E0 = (397.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS138',
    E0 = (340.061,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS139',
    E0 = (159.154,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS140',
    E0 = (159.154,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS141',
    E0 = (273.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS142',
    E0 = (650.057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS143',
    E0 = (359.279,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS144',
    E0 = (573.307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS145',
    E0 = (320.576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS146',
    E0 = (368.369,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS147',
    E0 = (367.293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS148',
    E0 = (338.142,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS149',
    E0 = (375.379,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS150',
    E0 = (272.955,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS151',
    E0 = (318.3,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS152',
    E0 = (320.047,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS153',
    E0 = (560.679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS154',
    E0 = (585.007,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS155',
    E0 = (616.018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS156',
    E0 = (685.686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS157',
    E0 = (676.337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS158',
    E0 = (341.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS159',
    E0 = (401.403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS160',
    E0 = (395.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS161',
    E0 = (246.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS162',
    E0 = (316.272,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS163',
    E0 = (346.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS164',
    E0 = (461.988,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS165',
    E0 = (535.401,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS166',
    E0 = (295.943,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS167',
    E0 = (667.087,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS168',
    E0 = (676.492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS169',
    E0 = (260.625,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS170',
    E0 = (215.99,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS171',
    E0 = (155.409,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS172',
    E0 = (240.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS173',
    E0 = (176.294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS174',
    E0 = (354.142,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS175',
    E0 = (413.394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS176',
    E0 = (301.841,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS177',
    E0 = (460.695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS178',
    E0 = (426.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS179',
    E0 = (422.306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS180',
    E0 = (334.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS181',
    E0 = (353.919,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS182',
    E0 = (289.603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS183',
    E0 = (499.637,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS184',
    E0 = (272.895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS185',
    E0 = (183.476,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS186',
    E0 = (413.587,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS187',
    E0 = (324.861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS188',
    E0 = (497.793,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS189',
    E0 = (357.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS190',
    E0 = (325.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS191',
    E0 = (318.453,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS192',
    E0 = (318.453,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS193',
    E0 = (237.314,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS194',
    E0 = (429.461,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS195',
    E0 = (548.734,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS196',
    E0 = (504.602,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS197',
    E0 = (220.668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS198',
    E0 = (351.337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS199',
    E0 = (204.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS200',
    E0 = (259.572,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS201',
    E0 = (195.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS202',
    E0 = (398.546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS203',
    E0 = (247.954,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS204',
    E0 = (650.057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS205',
    E0 = (369.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS206',
    E0 = (317.415,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS207',
    E0 = (330.019,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS208',
    E0 = (316.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS209',
    E0 = (312.178,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS210',
    E0 = (312.178,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS211',
    E0 = (527.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS212',
    E0 = (503.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS213',
    E0 = (517.353,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS214',
    E0 = (641.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS215',
    E0 = (641.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS216',
    E0 = (355.947,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS217',
    E0 = (311.82,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS218',
    E0 = (342.318,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS219',
    E0 = (259.178,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS220',
    E0 = (262.364,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS221',
    E0 = (354.626,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS222',
    E0 = (380.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS223',
    E0 = (334.709,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS224',
    E0 = (471.317,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS225',
    E0 = (337.99,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS226',
    E0 = (647.825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS227',
    E0 = (448.573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS228',
    E0 = (266.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS229',
    E0 = (347.558,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS230',
    E0 = (202.085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS231',
    E0 = (336.289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS232',
    E0 = (364.3,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS233',
    E0 = (359.279,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS234',
    E0 = (440.92,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS235',
    E0 = (576.505,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS236',
    E0 = (755.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS237',
    E0 = (410.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS238',
    E0 = (410.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS239',
    E0 = (302.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS240',
    E0 = (302.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS241',
    E0 = (297.883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS242',
    E0 = (297.883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS243',
    E0 = (194.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS244',
    E0 = (498.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS245',
    E0 = (402.932,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS246',
    E0 = (464.581,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS247',
    E0 = (464.123,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS248',
    E0 = (464.123,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS249',
    E0 = (235.339,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS250',
    E0 = (235.339,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS251',
    E0 = (280.384,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS252',
    E0 = (280.384,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS253',
    E0 = (313.246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS254',
    E0 = (193.009,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS255',
    E0 = (410.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS256',
    E0 = (302.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS257',
    E0 = (302.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS258',
    E0 = (297.883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS259',
    E0 = (202.802,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS260',
    E0 = (194.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS261',
    E0 = (235.339,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS262',
    E0 = (280.384,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS263',
    E0 = (280.384,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS264',
    E0 = (313.246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS265',
    E0 = (424.219,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS266',
    E0 = (301.336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS267',
    E0 = (297.217,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS268',
    E0 = (297.217,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS269',
    E0 = (194.249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS270',
    E0 = (235.339,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS271',
    E0 = (313.246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS272',
    E0 = (310.075,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS273',
    E0 = (315.178,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS274',
    E0 = (408.093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS275',
    E0 = (312.604,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['H(3)(3)', 'O=C1C=C=CCC=C1(2904)'],
    products = ['C7H7O(476)(475)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.149e+09,'cm^3/(mol*s)'), n=1.595, Ea=(16.7946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-OneDeH;HJ] for rate rule [Ca_Cds-COH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(3)(3)', 'C7H6O(492)(491)'],
    products = ['C7H7O(476)(475)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(9.22566,'m^3/(mol*s)'), n=2.04274, Ea=(10.5229,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 101 used for Cds-CdH_Cds-CdH;HJ
Exact match found for rate rule [Cds-CdH_Cds-CdH;HJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)(3)', 'C7H6O(488)(487)'],
    products = ['C7H7O(476)(475)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.28e-13,'cm^3/(molecule*s)'), n=0.611, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 53 used for C_rad/H/CdCd;H_rad
Exact match found for rate rule [C_rad/H/CdCd;H_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -1.8 to -1.8 kJ/mol.
Ea raised from -1.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)(3)', '[O]C1C=CC[C]=CC=1(2905)'],
    products = ['C7H7O(476)(475)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)(3)', '[O]C1C=[C]CC=CC=1(2906)'],
    products = ['C7H7O(476)(475)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)(3)', '[O]C1C=CCC=[C]C=1(2907)'],
    products = ['C7H7O(476)(475)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [Cd_rad/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)(3)', '[O]C1[C]=CCC=CC=1(2908)'],
    products = ['C7H7O(476)(475)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [H_rad;Cd_rad/Cd]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)(3)', '[O]C1=[C]C=CCC=C1(2909)'],
    products = ['C7H7O(476)(475)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C1C=CC=CC=C1(1734)'],
    products = ['C7H7O(476)(475)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.25e+09,'s^-1'), n=0.991, Ea=(45.9529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [1_3_pentadiene;CHOC_1;CdHC_2]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C1C=C=CCCC=1(2910)'],
    products = ['C7H7O(476)(475)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.00655e+09,'s^-1'), n=0.985631, Ea=(174.544,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [1_3_unsaturated_pentane_backbone;CH_end;CddC_2] + [1_3_unsaturated_pentane_backbone;CH2(C)_1;unsaturated_end] for rate rule [1_3_unsaturated_pentane_backbone;CH2(C)_1;CddC_2]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction11',
    reactants = ['OC1[C]=CCC=CC=1(2215)'],
    products = ['C7H7O(476)(475)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_2Cd;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['OC1=[C]C=CCC=C1(2216)'],
    products = ['C7H7O(476)(475)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.28371e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleDe_Cd;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['OC1C=[C]CC=CC=1(2213)'],
    products = ['C7H7O(476)(475)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;XH_out] for rate rule [R4H_DSS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['OC1C=CCC=[C]C=1(2217)'],
    products = ['C7H7O(476)(475)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.11e+08,'s^-1'), n=1.1915, Ea=(103.605,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R4H_SDS;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SDS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C7H7O(476)(475)'],
    products = ['C7H7O(495)(494)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(60205.5,'s^-1'), n=1.86417, Ea=(61.9987,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R5H_SS(D)MS;O_rad_out;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 3.31662479036
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['OC1C=CC[C]=CC=1(2214)'],
    products = ['C7H7O(476)(475)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_DSMS;Cd_rad_out_single;XH_out] for rate rule [R5H_DSMS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C7H7O(476)(475)'],
    products = ['C7H7O(470)(469)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(9.4423e+12,'s^-1'), n=-0.141781, Ea=(213.662,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_cyclic;doublebond_intra;radadd_intra] for rate rule [Rn1c7_alpha_short;doublebond_intra_secDe_HCd;radadd_intra_O]
Euclidian distance = 3.74165738677
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C7H7O(476)(475)'],
    products = ['C7H7O(484)(483)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(130.857,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_cyclic;doublebond_intra_pri;radadd_intra] for rate rule [Rn1c7_beta_short;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C7H7O(476)(475)'],
    products = ['C1=C[C]2C=CC(C1)O2(1698)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.87374e+12,'s^-1'), n=-0.166203, Ea=(191.809,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_cyclic;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn1c7_gamma_short;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C7H7O(476)(475)'],
    products = ['C7H7O(478)(477)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.54269e+17,'s^-1'), n=-0.500262, Ea=(96.3509,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_cyclic;doublebond_intra_pri_HCd;radadd_intra] for rate rule [Rn1c7_beta_long;doublebond_intra_pri_HCd;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C7H7O(476)(475)'],
    products = ['[CH]1CC=CC=C2OC12(2911)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.99832e+10,'s^-1'), n=0.37247, Ea=(126.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri;radadd_intra] for rate rule [Rn1c7_alpha_long;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 3.16227766017
family: Intra_R_Add_Endocyclic
Ea raised from 125.9 to 126.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction22',
    reactants = ['O(4)(4)', '[C]1=CC=CCC=C1(1246)'],
    products = ['C7H7O(476)(475)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/OneDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['O=C1C=CCC=[C]C1(2176)'],
    products = ['C7H7O(476)(475)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.82494e+10,'s^-1'), n=0.9, Ea=(134.725,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/OneDe] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/CO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['O=C1C=CC[C]=CC1(2177)'],
    products = ['C7H7O(476)(475)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.182e+10,'s^-1'), n=0.86, Ea=(149.369,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_H/CO]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['O=C1[C]=CCC=CC1(2179)'],
    products = ['C7H7O(476)(475)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.32e+09,'s^-1'), n=0.99, Ea=(141.419,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_Cd;Cs_H_out_1H] for rate rule [R3H_SS;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['O=C1C=C[CH]C=CC1(2175)'],
    products = ['C7H7O(476)(475)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(3.56e+06,'s^-1'), n=1.75, Ea=(105.855,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H;C_rad_out_H/Cd;Cs_H_out_H/OneDe] for rate rule [R4H_SDS;C_rad_out_H/Cd;Cs_H_out_H/CO]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['O=C1C=[C]CC=CC1(2178)'],
    products = ['C7H7O(476)(475)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out_H/Cd] for rate rule [R4H_DSS;Cd_rad_out_Cs;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C7H7O(476)(475)'],
    products = ['O=C1C=CCC2[CH]C12(2912)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(4.00063e+13,'s^-1'), n=-0.283562, Ea=(245.321,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_beta;doublebond_intra_pri_HNd_Cs;radadd_intra_csHDe] for rate rule [Rn0c7_beta_short;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCO]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C7H7O(476)(475)'],
    products = ['C7H7O(486)(485)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(130.857,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_gamma_short;doublebond_intra_pri;radadd_intra_cs] for rate rule [Rn0c7_gamma_short;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCd]
Euclidian distance = 2.82842712475
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C7H7O(476)(475)'],
    products = ['O=C1C2[CH]CC=CC12(2913)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(3.125e+09,'s^-1'), n=0.76, Ea=(105.276,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_cyclic;doublebond_intra_pri_HDe;radadd_intra_csHDe] for rate rule [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HCO;radadd_intra_csHCO]
Euclidian distance = 3.31662479036
family: Intra_R_Add_Endocyclic
Ea raised from 101.0 to 105.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=CCC=CC=C=O(2914)'],
    products = ['C7H7O(476)(475)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.35207e+11,'s^-1'), n=0.0116667, Ea=(102.85,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7;doublebond_intra;radadd_intra_cdsingleH] + [R7_linear;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R7_linear;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction32',
    reactants = ['O=[C]C1C=CCC=C1(2915)'],
    products = ['C7H7O(476)(475)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.59223e+11,'s^-1'), n=0.58375, Ea=(142.981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CJ;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C7H7O(476)(475)'],
    products = ['[O]C12C=CCC1C=C2(2429)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(204.456,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD;multiplebond_intra;radadd_intra_csHNd] for rate rule [R5_SD_CO;carbonylbond_intra_De;radadd_intra_csHNd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 198.4 to 204.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction34',
    reactants = ['O=C1C=[C]CCC=C1(2916)'],
    products = ['C7H7O(476)(475)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(2.66329e+10,'s^-1'), n=0.993, Ea=(157.679,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['O=C1[C]=CCCC=C1(2917)'],
    products = ['C7H7O(476)(475)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_single;Cs_H_out_H/NonDeC] for rate rule [R3H_DS;Cd_rad_out_singleDe;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C7H7O(476)(475)'],
    products = ['C7H7O(481)(480)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(130.857,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_gamma_short;doublebond_intra_pri_HDe;radadd_intra_cs] for rate rule [Rn0c7_gamma_short;doublebond_intra_pri_HCO;radadd_intra_csHCd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C7H7O(476)(475)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2.1e+12,'s^-1'), n=0.14, Ea=(30.5535,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic
Ea raised from 26.3 to 30.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]=CC(=O)C=CC=C(2918)'],
    products = ['C7H7O(476)(475)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(4.625e+11,'s^-1'), n=0.16, Ea=(41.045,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;doublebond_intra_pri;radadd_intra_cdsingleH] for rate rule [R7_linear;doublebond_intra_pri_2H;radadd_intra_cdsingleH]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C1C=CC(=O)C=C1(2919)'],
    products = ['C7H7O(476)(475)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.31121e+11,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction40',
    reactants = ['H(3)(3)', 'C7H6O(487)(486)'],
    products = ['C7H7O(470)(469)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(4.61283,'m^3/(mol*s)'), n=2.04274, Ea=(10.5229,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 101 used for Cds-CdH_Cds-CdH;HJ
Exact match found for rate rule [Cds-CdH_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction41',
    reactants = ['H(3)(3)', 'C1=CCC=CC2OC=12(4314)'],
    products = ['C7H7O(470)(469)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(5.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(15.8155,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2714 used for Ca_Cds-CsH;HJ
Exact match found for rate rule [Ca_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction42',
    reactants = ['C7H7O(470)(469)'],
    products = ['[CH]1C=CC2OC2=CC1(4315)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.169e+11,'s^-1'), n=0.707, Ea=(116.068,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/Cd;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[C]1CCC=CC2OC=12(4316)'],
    products = ['C7H7O(470)(469)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(2.66329e+10,'s^-1'), n=0.993, Ea=(157.679,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[C]1=CC2OC2=CCC1(4317)'],
    products = ['C7H7O(470)(469)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(3.32e+09,'s^-1'), n=0.99, Ea=(141.419,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_1H] for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['C1=C[C]2OC2=CCC1(4318)'],
    products = ['C7H7O(470)(469)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(2e+10,'s^-1'), n=0, Ea=(418.4,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;C_rad_out_single;Cs_H_out_H/NonDeC] for rate rule [R4H_SDS;C_rad_out_OneDe/O;Cs_H_out_H/NonDeC]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[C]1=CCCC=C2OC12(4319)'],
    products = ['C7H7O(470)(469)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out_H/Cd] for rate rule [R4H_DSS;Cd_rad_out_Cs;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['H(3)(3)', '[CH]1C=C2O[C]2C=CC1(4320)'],
    products = ['C7H7O(470)(469)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(7.24e+13,'cm^3/(mol*s)'), n=0.228, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/TwoDe;H_rad] for rate rule [C_rad/TwoDeO;H_rad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -0.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction48',
    reactants = ['H(3)(3)', '[CH]1[CH]C=C2OC2C=C1(4321)'],
    products = ['C7H7O(470)(469)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(4.34601e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction49',
    reactants = ['H(3)(3)', '[C]1=CC[CH]C=C2OC12(4322)'],
    products = ['C7H7O(470)(469)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction50',
    reactants = ['H(3)(3)', '[C]1=CC2OC2=C[CH]C1(4323)'],
    products = ['C7H7O(470)(469)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction51',
    reactants = ['H(3)(3)', '[C]1[CH]CC=CC2OC=12(4324)'],
    products = ['C7H7O(470)(469)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(4.34601e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction52',
    reactants = ['C7H7O(470)(469)'],
    products = ['[CH]1C2CC=CC3OC123(4325)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(2.22857e+12,'s^-1'), n=0, Ea=(182.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_beta_short;doublebond_intra_pri;radadd_intra_cs] for rate rule [Rn0c7_beta_short;doublebond_intra_pri_NdNd;radadd_intra_csHCs]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction53',
    reactants = ['C7H7O(470)(469)'],
    products = ['C7H7O(473)(472)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(130.857,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_gamma_short;doublebond_intra_pri;radadd_intra_cs] for rate rule [Rn0c7_gamma_short;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCd]
Euclidian distance = 2.82842712475
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction54',
    reactants = ['C7H7O(470)(469)'],
    products = ['[CH]1C2OC2=CC2CC12(4326)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(2.1e+12,'s^-1'), n=0.14, Ea=(140.722,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction55',
    reactants = ['[CH]=CC1OC1=CC=C(4327)'],
    products = ['C7H7O(470)(469)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(4.625e+11,'s^-1'), n=0.16, Ea=(41.045,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;doublebond_intra_pri;radadd_intra_cdsingleH] for rate rule [R7_linear;doublebond_intra_pri_2H;radadd_intra_cdsingleH]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction56',
    reactants = ['O(S)(1202)', '[C]1=CC=CCC=C1(1246)'],
    products = ['C7H7O(470)(469)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(5.96556e+06,'m^3/(mol*s)'), n=0, Ea=(0.08368,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [o_atom_singlet;mb_db]
Euclidian distance = 0
family: 1+2_Cycloaddition"""),
)

reaction(
    label = 'reaction57',
    reactants = ['[CH2]C1C=CC2OC2=C1(4328)'],
    products = ['C7H7O(470)(469)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction58',
    reactants = ['C7H7O(470)(469)'],
    products = ['O=C1C2[CH]CC=CC12(2913)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction59',
    reactants = ['H(3)(3)', '[C]1C=C2OC2C=CC1(4329)'],
    products = ['C7H7O(470)(469)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction60',
    reactants = ['H(3)(3)', 'C1=CC2OC=2C=CC1(4330)'],
    products = ['C7H7O(470)(469)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(1.644e+09,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction61',
    reactants = ['[C]1=CCC=CC2OC12(1656)'],
    products = ['C7H7O(470)(469)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_noH] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_NDMustO]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction62',
    reactants = ['[C]1=CC2OC2C=CC1(1657)'],
    products = ['C7H7O(470)(469)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(8.05e+09,'s^-1'), n=0.86, Ea=(139.746,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;Cs_H_out_NonDe] for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_NDMustO]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction63',
    reactants = ['[CH]1C=CC2OC2C=C1(1655)'],
    products = ['C7H7O(470)(469)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(2.66833e+08,'s^-1'), n=0.875, Ea=(262.128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H;C_rad_out_H/Cd;Cs_H_out] + [R4H_SDS;C_rad_out_1H;Cs_H_out] for rate rule [R4H_SDS;C_rad_out_H/Cd;Cs_H_out_NDMustO]
Euclidian distance = 3.60555127546
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction64',
    reactants = ['C7H7O(470)(469)'],
    products = ['[CH]1C2CC=CC23OC13(4331)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(165.553,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_gamma_short;doublebond_intra_pri;radadd_intra_cs] for rate rule [Rn0c7_gamma_short;doublebond_intra_pri_HNd_Cs;radadd_intra_csNdCd]
Euclidian distance = 2.82842712475
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction65',
    reactants = ['C7H7O(470)(469)'],
    products = ['[CH]1CC=CC23OC2C13(4332)'],
    transitionState = 'TS65',
    kinetics = Arrhenius(A=(2.1e+12,'s^-1'), n=0.14, Ea=(162.591,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra_csNdNd]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction66',
    reactants = ['[CH]=CCC=CC1=CO1(4333)'],
    products = ['C7H7O(470)(469)'],
    transitionState = 'TS66',
    kinetics = Arrhenius(A=(4.625e+11,'s^-1'), n=0.16, Ea=(41.045,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_cyclic;doublebond_intra;radadd_intra_cdsingleH] for rate rule [Rn5c3_alpha_short;doublebond_intra_secNd;radadd_intra_cdsingleH]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction67',
    reactants = ['[CH]1OC12C=CCC=C2(4334)'],
    products = ['C7H7O(470)(469)'],
    transitionState = 'TS67',
    kinetics = Arrhenius(A=(2.66e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction68',
    reactants = ['H(3)(3)', 'C1=C2OC2C2=CCC12(4335)'],
    products = ['C7H7O(473)(472)'],
    transitionState = 'TS68',
    kinetics = Arrhenius(A=(7.72e+07,'cm^3/(mol*s)'), n=1.64, Ea=(9.07928,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2568 used for Cds-CsCs_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsCs_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction69',
    reactants = ['H(3)(3)', 'C1=CC2C1C=C1OC12(4336)'],
    products = ['C7H7O(473)(472)'],
    transitionState = 'TS69',
    kinetics = Arrhenius(A=(1.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction70',
    reactants = ['C7H7O(473)(472)'],
    products = ['C1=C2OC2[C]2CCC12(4337)'],
    transitionState = 'TS70',
    kinetics = Arrhenius(A=(7.25e+10,'s^-1'), n=0.6, Ea=(154.39,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_Cs2] for rate rule [R2H_S_cy4;C_rad_out_H/NonDeC;Cs_H_out_Cs2]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction71',
    reactants = ['[CH]1CC2C1C=C1OC12(4338)'],
    products = ['C7H7O(473)(472)'],
    transitionState = 'TS71',
    kinetics = Arrhenius(A=(6.76e+09,'s^-1'), n=0.88, Ea=(158.992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC] for rate rule [R2H_S_cy4;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs/Cs)]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction72',
    reactants = ['C1[C]2CCC2C2OC=12(4339)'],
    products = ['C7H7O(473)(472)'],
    transitionState = 'TS72',
    kinetics = Arrhenius(A=(2.00896e+09,'s^-1'), n=1.09289, Ea=(174.668,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_noH;Cs_H_out_H/NonDeC] for rate rule [R3H_SS_12cy5;C_rad_out_OneDe/Cs;Cs_H_out_H/NonDeC]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction73',
    reactants = ['C7H7O(473)(472)'],
    products = ['C1=C2O[C]2C2CCC12(4340)'],
    transitionState = 'TS73',
    kinetics = Arrhenius(A=(85.4236,'s^-1'), n=3.08967, Ea=(142.368,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_H/NonDeC;Cs_H_out_noH] for rate rule [R3H_SS_12cy4;C_rad_out_H/NonDeC;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction74',
    reactants = ['[C]1=C2OC2C2CCC12(4341)'],
    products = ['C7H7O(473)(472)'],
    transitionState = 'TS74',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out_1H] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction75',
    reactants = ['H(3)(3)', '[CH]1CC2C=C3OC3[C]12(4342)'],
    products = ['C7H7O(473)(472)'],
    transitionState = 'TS75',
    kinetics = Arrhenius(A=(4.34601e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction76',
    reactants = ['H(3)(3)', '[CH]1C[C]2C=C3OC3C12(4343)'],
    products = ['C7H7O(473)(472)'],
    transitionState = 'TS76',
    kinetics = Arrhenius(A=(2.92e+13,'cm^3/(mol*s)'), n=0.18, Ea=(0.518816,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 123 used for H_rad;C_rad/OneDeCs2
Exact match found for rate rule [C_rad/OneDeCs2;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction77',
    reactants = ['H(3)(3)', '[CH]1CC2C=C3O[C]3C12(4344)'],
    products = ['C7H7O(473)(472)'],
    transitionState = 'TS77',
    kinetics = Arrhenius(A=(2.92e+13,'cm^3/(mol*s)'), n=0.18, Ea=(0.518816,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""Estimated using template [C_rad/OneDe;H_rad] for rate rule [C_rad/OneDeO;H_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction78',
    reactants = ['H(3)(3)', '[CH]1[CH]C2C1C=C1OC12(4345)'],
    products = ['C7H7O(473)(472)'],
    transitionState = 'TS78',
    kinetics = Arrhenius(A=(4.34601e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction79',
    reactants = ['H(3)(3)', '[C]1=C2OC2C2[CH]CC12(4346)'],
    products = ['C7H7O(473)(472)'],
    transitionState = 'TS79',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [H_rad;Cd_rad/NonDe]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction80',
    reactants = ['C7H7O(473)(472)'],
    products = ['[CH]1C2CC3C2C2OC132(4347)'],
    transitionState = 'TS80',
    kinetics = Arrhenius(A=(3.24409e+08,'s^-1'), n=0.869444, Ea=(176.269,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R5_cyclic;doublebond_intra_pri_NdNd;radadd_intra_csHCs] for rate rule [Rn1c6_gamma;doublebond_intra_pri_NdNd;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction81',
    reactants = ['C7H7O(473)(472)'],
    products = ['C1C2C3[C]4OC4C2C13(4348)'],
    transitionState = 'TS81',
    kinetics = Arrhenius(A=(6.33004e+14,'s^-1'), n=-0.792922, Ea=(179.405,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Rn0c7_gamma_long_SSS_D;doublebond_intra;radadd_intra_csHCs] for rate rule [Rn0c7_gamma_long_SSS_D;doublebond_intra_secNd_HNd;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction82',
    reactants = ['C=CC1[CH]C=C2OC21(4349)'],
    products = ['C7H7O(473)(472)'],
    transitionState = 'TS82',
    kinetics = Arrhenius(A=(3.72602e+07,'s^-1'), n=1.21458, Ea=(125.171,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs] for rate rule [R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_csHCd]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction83',
    reactants = ['[CH]1OC1=CC1C=CC1(4350)'],
    products = ['C7H7O(473)(472)'],
    transitionState = 'TS83',
    kinetics = Arrhenius(A=(8.83927e+09,'s^-1'), n=0.393333, Ea=(66.944,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;doublebond_intra_pri;radadd_intra_csHNd] for rate rule [R7_linear;doublebond_intra_pri;radadd_intra_csHO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction84',
    reactants = ['O(S)(1202)', 'C1=CC2[CH]CC2C=1(4351)'],
    products = ['C7H7O(473)(472)'],
    transitionState = 'TS84',
    kinetics = Arrhenius(A=(5.96556e+06,'m^3/(mol*s)'), n=0, Ea=(0.08368,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [o_atom_singlet;mb_db]
Euclidian distance = 0
family: 1+2_Cycloaddition"""),
)

reaction(
    label = 'reaction85',
    reactants = ['[CH]1C2OC2=CC2CC12(4326)'],
    products = ['C7H7O(473)(472)'],
    transitionState = 'TS85',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction86',
    reactants = ['[CH]1C2C=C3OC3C1C2(4352)'],
    products = ['C7H7O(473)(472)'],
    transitionState = 'TS86',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction87',
    reactants = ['[CH2]C1C2C=C3OC3C12(4353)'],
    products = ['C7H7O(473)(472)'],
    transitionState = 'TS87',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction88',
    reactants = ['C7H7O(473)(472)'],
    products = ['O=C1C2C3[CH]CC3C12(4354)'],
    transitionState = 'TS88',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond;R_O_C] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond;R_O_C]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction89',
    reactants = ['H(3)(3)', '[C]1CC2C=C3OC3C12(4355)'],
    products = ['C7H7O(473)(472)'],
    transitionState = 'TS89',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction90',
    reactants = ['OC1=[C]C2CC2C=C1(4356)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS90',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;XH_out] for rate rule [R3H_DS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction91',
    reactants = ['OC1[C]=CC2CC2C=1(4357)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS91',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_2Cd;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction92',
    reactants = ['C7H7O(477)(476)'],
    products = ['OC1C=CC2C[C]2C=1(4358)'],
    transitionState = 'TS92',
    kinetics = Arrhenius(A=(1e+11,'s^-1'), n=0, Ea=(100.834,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4H;O_rad_out;Cs_H_out_Cs2] for rate rule [R4H_SDS;O_rad_out;Cs_H_out_Cs2_cy3]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction93',
    reactants = ['OC1C=[C]C2CC2C=1(4359)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS93',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;XH_out] for rate rule [R4H_DSS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction94',
    reactants = ['C7H7O(477)(476)'],
    products = ['OC1C=C[C]2CC2C=1(4360)'],
    transitionState = 'TS94',
    kinetics = Arrhenius(A=(126000,'s^-1'), n=1.85, Ea=(88.2824,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SMSS;Y_rad_out;Cs_H_out] for rate rule [R5H_SMSS;O_rad_out;Cs_H_out_noH]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction95',
    reactants = ['OC1C=CC2[CH]C2C=1(4361)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS95',
    kinetics = Arrhenius(A=(148.499,'s^-1'), n=2.47528, Ea=(59.2966,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5H;C_rad_out_H/NonDeC;XH_out] + [R5H_SSMS;C_rad_out_single;XH_out] for rate rule [R5H_SSMS;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction96',
    reactants = ['H(3)(3)', '[O]C1C=C[C]2CC2C=1(4362)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS96',
    kinetics = Arrhenius(A=(2.92e+13,'cm^3/(mol*s)'), n=0.18, Ea=(0.518816,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 123 used for H_rad;C_rad/OneDeCs2
Exact match found for rate rule [C_rad/OneDeCs2;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction97',
    reactants = ['H(3)(3)', '[O]C1C=CC2C[C]2C=1(4363)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS97',
    kinetics = Arrhenius(A=(2.92e+13,'cm^3/(mol*s)'), n=0.18, Ea=(0.518816,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 123 used for H_rad;C_rad/OneDeCs2
Exact match found for rate rule [C_rad/OneDeCs2;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction98',
    reactants = ['H(3)(3)', '[O]C1C=CC2[CH]C2C=1(4364)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS98',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction99',
    reactants = ['H(3)(3)', '[O]C1=[C]C2CC2C=C1(4365)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS99',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction100',
    reactants = ['H(3)(3)', '[O]C1C=[C]C2CC2C=1(4366)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS100',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction101',
    reactants = ['H(3)(3)', '[O]C1[C]=CC2CC2C=1(4367)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS101',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [H_rad;Cd_rad/Cd]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction102',
    reactants = ['C7H7O(477)(476)'],
    products = ['C1=CC2CC2C2O[C]12(4368)'],
    transitionState = 'TS102',
    kinetics = Arrhenius(A=(9.4423e+12,'s^-1'), n=-0.141781, Ea=(213.662,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_cyclic;doublebond_intra;radadd_intra] for rate rule [Rn1c6_alpha_short;doublebond_intra_secDe_HNd;radadd_intra_O]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction103',
    reactants = ['C7H7O(477)(476)'],
    products = ['[CH]1C2=CC3CC3C1O2(4369)'],
    transitionState = 'TS103',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(278.822,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_cyclic;doublebond_intra_pri;radadd_intra] for rate rule [Rn1c6_beta_short;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic
Ea raised from 278.6 to 278.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction104',
    reactants = ['C7H7O(477)(476)'],
    products = ['[CH]1C2OC2=CC2CC12(4326)'],
    transitionState = 'TS104',
    kinetics = Arrhenius(A=(2.48e+11,'s^-1'), n=0.26, Ea=(187.242,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn1c6_alpha_long;doublebond_intra_pri;radadd_intra] for rate rule [Rn1c6_alpha_long;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 185.0 to 187.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction105',
    reactants = ['CH2(S)(21)(22)', '[O]c1ccccc1(4370)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS105',
    kinetics = Arrhenius(A=(1.54e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [elec_def;mb_db_twocdisub] for rate rule [carbene;mb_db_twocdisub_twoDe]
Euclidian distance = 1.41421356237
family: 1+2_Cycloaddition"""),
)

reaction(
    label = 'reaction106',
    reactants = ['[O]C1=CCC2CC2=C1(4371)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS106',
    kinetics = Arrhenius(A=(9.66e+08,'s^-1'), n=1.162, Ea=(185.28,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_3_unsaturated_pentane_backbone;CH2(C)_1;unsaturated_end] for rate rule [1_3_unsaturated_pentane_backbone;CH2(C)_1;Cd(C)C_2]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction107',
    reactants = ['O(4)(4)', '[C]1=CC2CC2C=C1(4372)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS107',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/OneDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction108',
    reactants = ['H(3)(3)', 'O=C1C=CC2CC2=C1(4373)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS108',
    kinetics = Arrhenius(A=(89.419,'m^3/(mol*s)'), n=1.598, Ea=(7.0333,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsCs_Cds-OneDeH;HJ] for rate rule [Cds-CsCs_Cds-COH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction109',
    reactants = ['O=C1C=CC2C[C]2C1(4374)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS109',
    kinetics = Arrhenius(A=(20108.5,'s^-1'), n=2.606, Ea=(121.63,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_Cs2;Cs_H_out_H/OneDe] for rate rule [R2H_S;C_rad_out_Cs2_cy3;Cs_H_out_H/CO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction110',
    reactants = ['O=C1C=C[C]2CC2C1(4375)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS110',
    kinetics = Arrhenius(A=(1.28e+07,'s^-1'), n=1.56, Ea=(126.775,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_noH;Cs_H_out_H/OneDe] for rate rule [R3H_SS_12cy3;C_rad_out_OneDe/Cs;Cs_H_out_H/CO]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction111',
    reactants = ['O=C1C=CC2[CH]C2C1(4376)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS111',
    kinetics = Arrhenius(A=(6.82e+09,'s^-1'), n=0.73, Ea=(127.612,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_H/NonDeC;Cs_H_out_H/OneDe] for rate rule [R3H_SS_12cy3;C_rad_out_H/NonDeC;Cs_H_out_H/CO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction112',
    reactants = ['O=C1[C]=CC2CC2C1(4377)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS112',
    kinetics = Arrhenius(A=(3.32e+09,'s^-1'), n=0.99, Ea=(141.419,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_Cd;Cs_H_out_H/NonDeC] for rate rule [R3H_SS;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs/Cs)]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction113',
    reactants = ['O=C1C=[C]C2CC2C1(4378)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS113',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out_H/OneDe] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_H/CO]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction114',
    reactants = ['C7H7O(477)(476)'],
    products = ['O=C1[CH]C2C3CC3C12(4379)'],
    transitionState = 'TS114',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(130.857,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_gamma;doublebond_intra_pri;radadd_intra_csHCs] for rate rule [Rn0c6_gamma;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction115',
    reactants = ['C7H7O(477)(476)'],
    products = ['O=C1C2[CH]C3CC3C12(4380)'],
    transitionState = 'TS115',
    kinetics = Arrhenius(A=(1.98e+12,'s^-1'), n=0.5, Ea=(255.224,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c6_beta_long_SS_D;doublebond_intra_pri_HDe;radadd_intra_csHDe] for rate rule [Rn0c6_beta_long_SS_D;doublebond_intra_pri_HCO;radadd_intra_csHCO]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction116',
    reactants = ['[CH]=CC1CC1C=C=O(4381)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS116',
    kinetics = Arrhenius(A=(1.03e+10,'s^-1'), n=0.19, Ea=(16.736,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R6_DSS_D;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction117',
    reactants = ['[CH2]C1C=CC(=O)C=C1(2919)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS117',
    kinetics = Arrhenius(A=(3.53542e+11,'s^-1'), n=0.21, Ea=(35.9824,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7_cyclic;doublebond_intra_pri;radadd_intra_cs2H] + [Rn1c6_alpha_long;doublebond_intra_pri;radadd_intra_cs] for rate rule [Rn1c6_alpha_long;doublebond_intra_pri;radadd_intra_cs2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction118',
    reactants = ['C7H7O(481)(480)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS118',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction119',
    reactants = ['O=C1C=CC2[CH]C1C2(4382)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS119',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction120',
    reactants = ['O=[C]C1C=CC2CC12(4383)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS120',
    kinetics = Arrhenius(A=(1.29612e+11,'s^-1'), n=0.58375, Ea=(142.981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction121',
    reactants = ['H(3)(3)', 'C1=CC2=CC(=CC1)O2(4384)'],
    products = ['C7H7O(478)(477)'],
    transitionState = 'TS121',
    kinetics = Arrhenius(A=(87.5179,'m^3/(mol*s)'), n=1.66467, Ea=(9.4546,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;HJ] for rate rule [Cds-CdOs_Cds;HJ]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction122',
    reactants = ['H(3)(3)', 'S(801)(800)'],
    products = ['C7H7O(478)(477)'],
    transitionState = 'TS122',
    kinetics = Arrhenius(A=(0.0272924,'m^3/(mol*s)'), n=2.81111, Ea=(21.1569,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 26 used for Cds-CdH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CdH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction123',
    reactants = ['C7H7O(478)(477)'],
    products = ['C1=CC2=C[C](CC1)O2(4385)'],
    transitionState = 'TS123',
    kinetics = Arrhenius(A=(2.15233e+09,'s^-1'), n=1.014, Ea=(127.919,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction124',
    reactants = ['C7H7O(484)(483)'],
    products = ['C7H7O(478)(477)'],
    transitionState = 'TS124',
    kinetics = Arrhenius(A=(1.022e+10,'s^-1'), n=1.34, Ea=(199.577,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 9 used for R2H_S;C_rad_out_H/(Cd-Cd-Cd);Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/(Cd-Cd-Cd);Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction125',
    reactants = ['[C]1=C2C=CCCC1O2(4386)'],
    products = ['C7H7O(478)(477)'],
    transitionState = 'TS125',
    kinetics = Arrhenius(A=(3.32e+09,'s^-1'), n=0.99, Ea=(141.419,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_Cd;Cs_H_out_H/NonDeC] for rate rule [R3H_SS_12cy4;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction126',
    reactants = ['[C]1=CC2=CC(CC1)O2(4387)'],
    products = ['C7H7O(478)(477)'],
    transitionState = 'TS126',
    kinetics = Arrhenius(A=(3.32e+09,'s^-1'), n=0.99, Ea=(141.419,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 206 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction127',
    reactants = ['[C]1=CCCC2C=C1O2(4388)'],
    products = ['C7H7O(478)(477)'],
    transitionState = 'TS127',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction128',
    reactants = ['H(3)(3)', '[CH]1CC=CC2=C[C]1O2(4389)'],
    products = ['C7H7O(478)(477)'],
    transitionState = 'TS128',
    kinetics = Arrhenius(A=(4.34601e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction129',
    reactants = ['H(3)(3)', '[CH]1[CH]C2C=C(C=C1)O2(4390)'],
    products = ['C7H7O(478)(477)'],
    transitionState = 'TS129',
    kinetics = Arrhenius(A=(4.34601e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction130',
    reactants = ['H(3)(3)', '[C]1=C2C=CC[CH]C1O2(4391)'],
    products = ['C7H7O(478)(477)'],
    transitionState = 'TS130',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction131',
    reactants = ['H(3)(3)', '[C]1=CC2=CC([CH]C1)O2(4392)'],
    products = ['C7H7O(478)(477)'],
    transitionState = 'TS131',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction132',
    reactants = ['H(3)(3)', '[C]1=CC[CH]C2C=C1O2(4393)'],
    products = ['C7H7O(478)(477)'],
    transitionState = 'TS132',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [H_rad;Cd_rad/Cd]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction133',
    reactants = ['C7H7O(478)(477)'],
    products = ['[CH]1C2OC13C=CCC23(4394)'],
    transitionState = 'TS133',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(286.744,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Rn0c7_gamma_short_SSDS_D;doublebond_intra_pri;radadd_intra_csHCs] for rate rule [Rn0c7_gamma_short_SSDS_D;doublebond_intra_pri_NdCd;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 283.6 to 286.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction134',
    reactants = ['C7H7O(478)(477)'],
    products = ['[CH]1CC2C3C=C(O3)C12(4395)'],
    transitionState = 'TS134',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(435.751,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 116 used for Rn0c7_gamma_short_SSDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCs
Exact match found for rate rule [Rn0c7_gamma_short_SSDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 432.9 to 435.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction135',
    reactants = ['C7H7O(478)(477)'],
    products = ['C1=C[C]2OC3C(C1)C23(4396)'],
    transitionState = 'TS135',
    kinetics = Arrhenius(A=(2.69283e+12,'s^-1'), n=-0.0647836, Ea=(205.681,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_cyclic;doublebond_intra;radadd_intra_csHCs] for rate rule [Rn1c4_alpha_long;doublebond_intra_secDe;radadd_intra_csHCs]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic
Ea raised from 204.4 to 205.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction136',
    reactants = ['C7H7O(478)(477)'],
    products = ['[CH]1C2=CC3CC3C1O2(4369)'],
    transitionState = 'TS136',
    kinetics = Arrhenius(A=(2.1e+12,'s^-1'), n=0.14, Ea=(285.589,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic
Ea raised from 283.3 to 285.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction137',
    reactants = ['[CH]=CC1=CC(C=C)O1(4397)'],
    products = ['C7H7O(478)(477)'],
    transitionState = 'TS137',
    kinetics = Arrhenius(A=(4.625e+11,'s^-1'), n=0.16, Ea=(41.045,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;doublebond_intra_pri;radadd_intra_cdsingleH] for rate rule [R7_linear;doublebond_intra_pri_2H;radadd_intra_cdsingleH]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction138',
    reactants = ['[CH]=C1C=CCC=CO1(4398)'],
    products = ['C7H7O(478)(477)'],
    transitionState = 'TS138',
    kinetics = Arrhenius(A=(4.625e+11,'s^-1'), n=0.16, Ea=(41.045,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_cyclic;doublebond_intra_pri;radadd_intra_cdsingleH] for rate rule [Rn1c7_beta_long;doublebond_intra_pri_HNd_O;radadd_intra_cdsingleH]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction139',
    reactants = ['C7H7O(478)(477)'],
    products = ['[CH]1CCC=C2C=C1O2(4399)'],
    transitionState = 'TS139',
    kinetics = Arrhenius(A=(2.125e+09,'s^-1'), n=0.991, Ea=(45.9529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_3_pentadiene;CH_end;CdHC_2] for rate rule [1_3_pentadiene;CH(CJ)_1;CdHC_2]
Euclidian distance = 1.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction140',
    reactants = ['C7H7O(478)(477)'],
    products = ['[CH]1C=CC=C2CC1O2(4400)'],
    transitionState = 'TS140',
    kinetics = Arrhenius(A=(4.25e+09,'s^-1'), n=0.991, Ea=(45.9529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_3_pentadiene;CH_end;CdHC_2] for rate rule [1_3_pentadiene;CH(CJ)_1;CdHC_2]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction141',
    reactants = ['C7H7O(478)(477)'],
    products = ['[CH]1OC2C=CCC1C=2(4401)'],
    transitionState = 'TS141',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction142',
    reactants = ['[CH2]C1C=CC2=CC1O2(4402)'],
    products = ['C7H7O(478)(477)'],
    transitionState = 'TS142',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction143',
    reactants = ['C7H7O(478)(477)'],
    products = ['[CH]1CC2C=C3OC1C32(4403)'],
    transitionState = 'TS143',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(246.077,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1,3-butadiene_backbone;C=C_1;C=C_2] for rate rule [1,3-butadiene_backbone;CdH(C)_1;CdH(C)_2]
Euclidian distance = 1.41421356237
family: Intra_2+2_cycloaddition_Cd
Ea raised from 243.3 to 246.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction144',
    reactants = ['H(3)(3)', '[C]1CC=CC2=CC1O2(4404)'],
    products = ['C7H7O(478)(477)'],
    transitionState = 'TS144',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction145',
    reactants = ['C7H7O(481)(480)'],
    products = ['[O]C12C=CC3CC1C32(4405)'],
    transitionState = 'TS145',
    kinetics = Arrhenius(A=(6.44049e+09,'s^-1'), n=0.679905, Ea=(104.189,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_csHNd] for rate rule [R4_S_CO;carbonylbond_intra_De;radadd_intra_csHNd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction146',
    reactants = ['H(3)(3)', 'O=C1C=CC2CC=C12(4406)'],
    products = ['C7H7O(481)(480)'],
    transitionState = 'TS146',
    kinetics = Arrhenius(A=(72.3521,'m^3/(mol*s)'), n=1.66655, Ea=(10.8198,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDeCs_Cds;HJ] for rate rule [Cds-COCs_Cds;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction147',
    reactants = ['H(3)(3)', 'C7H6O(535)(534)'],
    products = ['C7H7O(481)(480)'],
    transitionState = 'TS147',
    kinetics = Arrhenius(A=(1.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction148',
    reactants = ['C7H7O(481)(480)'],
    products = ['O=C1C=CC2CC[C]12(4407)'],
    transitionState = 'TS148',
    kinetics = Arrhenius(A=(8.83e+10,'s^-1'), n=0.3, Ea=(121.754,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_OneDe] for rate rule [R2H_S_cy4;C_rad_out_H/NonDeC;Cs_H_out_CO]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction149',
    reactants = ['C7H7O(481)(480)'],
    products = ['O=C1C=CC2[CH]CC12(4408)'],
    transitionState = 'TS149',
    kinetics = Arrhenius(A=(6.76e+09,'s^-1'), n=0.88, Ea=(158.992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC] for rate rule [R2H_S_cy4;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction150',
    reactants = ['C7H7O(481)(480)'],
    products = ['O=C1C=C[C]2CCC21(4409)'],
    transitionState = 'TS150',
    kinetics = Arrhenius(A=(2.76e-23,'s^-1'), n=10.17, Ea=(56.5677,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_H/NonDeC;Cs_H_out_OneDe] for rate rule [R3H_SS_23cy5;C_rad_out_H/NonDeC;Cs_H_out_Cd]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction151',
    reactants = ['O=C1C=[C]C2CCC12(4410)'],
    products = ['C7H7O(481)(480)'],
    transitionState = 'TS151',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out_1H] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction152',
    reactants = ['O=C1[C]=CC2CCC12(4411)'],
    products = ['C7H7O(481)(480)'],
    transitionState = 'TS152',
    kinetics = Arrhenius(A=(236531,'s^-1'), n=1.93, Ea=(59.6011,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;Cs_H_out_H/(NonDeC/Cs)] + [R4H_RSS;Cd_rad_out;Cs_H_out_1H] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction153',
    reactants = ['H(3)(3)', 'O=C1C=C[C]2C[CH]C21(4412)'],
    products = ['C7H7O(481)(480)'],
    transitionState = 'TS153',
    kinetics = Arrhenius(A=(2.92e+13,'cm^3/(mol*s)'), n=0.18, Ea=(0.518816,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 123 used for H_rad;C_rad/OneDeCs2
Exact match found for rate rule [C_rad/OneDeCs2;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction154',
    reactants = ['H(3)(3)', 'O=C1C=CC2C[CH][C]12(4413)'],
    products = ['C7H7O(481)(480)'],
    transitionState = 'TS154',
    kinetics = Arrhenius(A=(4.34601e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction155',
    reactants = ['H(3)(3)', 'O=C1C=CC2[CH][CH]C12(4414)'],
    products = ['C7H7O(481)(480)'],
    transitionState = 'TS155',
    kinetics = Arrhenius(A=(4.34601e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction156',
    reactants = ['H(3)(3)', 'O=C1C=[C]C2C[CH]C12(4415)'],
    products = ['C7H7O(481)(480)'],
    transitionState = 'TS156',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction157',
    reactants = ['H(3)(3)', 'O=C1[C]=CC2C[CH]C12(4416)'],
    products = ['C7H7O(481)(480)'],
    transitionState = 'TS157',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Cd_rad/OneDe]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction158',
    reactants = ['C7H7O(481)(480)'],
    products = ['C1=CC2CC3O[C]1C23(4417)'],
    transitionState = 'TS158',
    kinetics = Arrhenius(A=(4.7342e+08,'s^-1'), n=0.920995, Ea=(125.348,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_csHCs] for rate rule [R4_S_CO;carbonyl_intra_De;radadd_intra_csHCs]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction159',
    reactants = ['C7H7O(481)(480)'],
    products = ['O=C1C2[CH]C3CC2C13(4418)'],
    transitionState = 'TS159',
    kinetics = Arrhenius(A=(1.03212e+12,'s^-1'), n=0.0431209, Ea=(185.016,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_cyclic;doublebond_intra_pri;radadd_intra_csHCs] for rate rule [Rn1c5_beta_long;doublebond_intra_pri_HCO;radadd_intra_csHCs]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction160',
    reactants = ['C7H7O(481)(480)'],
    products = ['O=C1[CH]C2C3CC2C13(4419)'],
    transitionState = 'TS160',
    kinetics = Arrhenius(A=(6.33004e+14,'s^-1'), n=-0.792922, Ea=(179.405,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 112 used for Rn0c7_gamma_long_SSS_D;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs
Exact match found for rate rule [Rn0c7_gamma_long_SSS_D;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction161',
    reactants = ['C=CC1[CH]C=CC1=O(4420)'],
    products = ['C7H7O(481)(480)'],
    transitionState = 'TS161',
    kinetics = Arrhenius(A=(3.72602e+07,'s^-1'), n=1.21458, Ea=(125.171,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs] for rate rule [R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_csHCd]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction162',
    reactants = ['O=[C]C=CC1C=CC1(4421)'],
    products = ['C7H7O(481)(480)'],
    transitionState = 'TS162',
    kinetics = Arrhenius(A=(9.09021e+13,'s^-1'), n=-0.0858201, Ea=(61.1948,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_cyclic;doublebond_intra_pri;radadd_intra] for rate rule [Rn3c4_alpha_long;doublebond_intra_pri;radadd_intra_CO]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction163',
    reactants = ['O=C1C=CC2[CH]C1C2(4382)'],
    products = ['C7H7O(481)(480)'],
    transitionState = 'TS163',
    kinetics = Arrhenius(A=(8.889e+11,'s^-1'), n=0.232, Ea=(122.75,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CO] for rate rule [cCs(-HC)CJ;CsJ-CsH;CO]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction164',
    reactants = ['C7H7O(481)(480)'],
    products = ['[CH2]C1C2C=CC(=O)C12(4422)'],
    transitionState = 'TS164',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HH)CJ;CsJ;C] for rate rule [cCs(-HH)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction165',
    reactants = ['C1=C[C]2C=CC(C1)O2(1698)'],
    products = ['C7H7O(481)(480)'],
    transitionState = 'TS165',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction166',
    reactants = ['OC1C=CC2C[CH]C=12(4423)'],
    products = ['C7H7O(481)(480)'],
    transitionState = 'TS166',
    kinetics = Arrhenius(A=(1290.48,'s^-1'), n=2.90375, Ea=(139.674,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction167',
    reactants = ['OC1=C=CC2C[CH]C12(4424)'],
    products = ['C7H7O(481)(480)'],
    transitionState = 'TS167',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction168',
    reactants = ['H(3)(3)', 'O=C1C=CC2C[C]C12(4425)'],
    products = ['C7H7O(481)(480)'],
    transitionState = 'TS168',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction169',
    reactants = ['H(3)(3)', 'C1=CCC2=CC(=C1)O2(4426)'],
    products = ['C7H7O(484)(483)'],
    transitionState = 'TS169',
    kinetics = Arrhenius(A=(170.641,'m^3/(mol*s)'), n=1.56204, Ea=(11.2897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;HJ] for rate rule [Cds-OsCs_Cds;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction170',
    reactants = ['C1=CC[C]2CC(=C1)O2(4427)'],
    products = ['C7H7O(484)(483)'],
    transitionState = 'TS170',
    kinetics = Arrhenius(A=(20108.5,'s^-1'), n=2.606, Ea=(121.63,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_NonDe;Cs_H_out_H/Cd] for rate rule [R2H_S_cy4;C_rad_out_NDMustO;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction171',
    reactants = ['C7H7O(484)(483)'],
    products = ['[CH]1C=CC=C2CC1O2(4400)'],
    transitionState = 'TS171',
    kinetics = Arrhenius(A=(3.8e+10,'s^-1'), n=0.87, Ea=(144.348,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_H/OneDe;Cs_H_out_H/Cd] for rate rule [R3H_SS_12cy4;C_rad_out_H/(Cd-Cd-Cd);Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction172',
    reactants = ['[C]1C=CCC2CC=1O2(4428)'],
    products = ['C7H7O(484)(483)'],
    transitionState = 'TS172',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_single;Cs_H_out_H/NonDeC] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction173',
    reactants = ['[C]1=CC=C2CC(C1)O2(4429)'],
    products = ['C7H7O(484)(483)'],
    transitionState = 'TS173',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out_H/(Cd-Cd-Cd)] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction174',
    reactants = ['[C]1=CCC2CC(=C1)O2(4430)'],
    products = ['C7H7O(484)(483)'],
    transitionState = 'TS174',
    kinetics = Arrhenius(A=(2.10713e+09,'s^-1'), n=0.59575, Ea=(261.003,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SDS;Y_rad_out;Cs_H_out_H/(NonDeC/Cs)] + [R4H_SDS;Cd_rad_out_Cd;Cs_H_out] for rate rule [R4H_SDS;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction175',
    reactants = ['H(3)(3)', '[CH]1[C]2CC=CC=C1O2(4431)'],
    products = ['C7H7O(484)(483)'],
    transitionState = 'TS175',
    kinetics = Arrhenius(A=(4.34601e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction176',
    reactants = ['H(3)(3)', '[CH]1[CH]C2C=C(C=C1)O2(4390)'],
    products = ['C7H7O(484)(483)'],
    transitionState = 'TS176',
    kinetics = Arrhenius(A=(2.71464e+07,'m^3/(mol*s)'), n=0.107721, Ea=(5.76381,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 36 used for C_rad/H/CdCs;H_rad
Exact match found for rate rule [C_rad/H/CdCs;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction177',
    reactants = ['H(3)(3)', '[C]1=CC=C2[CH]C(C1)O2(4432)'],
    products = ['C7H7O(484)(483)'],
    transitionState = 'TS177',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction178',
    reactants = ['H(3)(3)', '[C]1=CCC2[CH]C(=C1)O2(4433)'],
    products = ['C7H7O(484)(483)'],
    transitionState = 'TS178',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [H_rad;Cd_rad/Cd]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction179',
    reactants = ['H(3)(3)', '[C]1C=CCC2[CH]C=1O2(4434)'],
    products = ['C7H7O(484)(483)'],
    transitionState = 'TS179',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction180',
    reactants = ['C7H7O(484)(483)'],
    products = ['C1=CC2[C]3OC(C1)C32(4435)'],
    transitionState = 'TS180',
    kinetics = Arrhenius(A=(2.22857e+12,'s^-1'), n=0, Ea=(323.219,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_beta_short;doublebond_intra;radadd_intra_cs] for rate rule [Rn0c7_beta_short;doublebond_intra_secNd_HCd;radadd_intra_csHCs]
Euclidian distance = 3.60555127546
family: Intra_R_Add_Endocyclic
Ea raised from 321.9 to 323.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction181',
    reactants = ['C7H7O(484)(483)'],
    products = ['[CH]1CC2OC3=CC1C32(4436)'],
    transitionState = 'TS181',
    kinetics = Arrhenius(A=(3.3472e+15,'s^-1'), n=-0.535753, Ea=(342.858,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_cyclic;doublebond_intra_pri_HCd;radadd_intra_csHCd] for rate rule [Rn1c7_gamma_short;doublebond_intra_pri_HCd;radadd_intra_csHCd]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 340.2 to 342.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction182',
    reactants = ['C7H7O(484)(483)'],
    products = ['C1=CC2CC3O[C]1C23(4417)'],
    transitionState = 'TS182',
    kinetics = Arrhenius(A=(8.03578e+12,'s^-1'), n=-0.324073, Ea=(278.542,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_cyclic;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs] for rate rule [Rn1c7_gamma_short;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 275.7 to 278.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction183',
    reactants = ['[O]C1C=C=CC=CC1(4437)'],
    products = ['C7H7O(484)(483)'],
    transitionState = 'TS183',
    kinetics = Arrhenius(A=(1.44078e+13,'s^-1'), n=-0.00629116, Ea=(146.333,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4_cyclic;multiplebond_intra;radadd_intra_O] + [R4_cyclic;doublebond_intra;radadd_intra] for rate rule [Rn1c7_beta_short;doublebond_intra;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction184',
    reactants = ['[CH2]C=CC=C1C=CO1(4438)'],
    products = ['C7H7O(484)(483)'],
    transitionState = 'TS184',
    kinetics = Arrhenius(A=(1.26e+11,'s^-1'), n=0.16, Ea=(41.84,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_cyclic;doublebond_intra_pri;radadd_intra_cs2H] for rate rule [Rn4c4_beta;doublebond_intra_pri_HNd_O;radadd_intra_cs2H]
Euclidian distance = 2.82842712475
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction185',
    reactants = ['[CH]1C2C=CC=CC1O2(4439)'],
    products = ['C7H7O(484)(483)'],
    transitionState = 'TS185',
    kinetics = Arrhenius(A=(4.25e+09,'s^-1'), n=0.991, Ea=(45.9529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_3_pentadiene;CH_end;CdHC_2] for rate rule [1_3_pentadiene;CH(CJ)_1;CdHC_2]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction186',
    reactants = ['[CH]1OC2=CC=CCC12(4440)'],
    products = ['C7H7O(484)(483)'],
    transitionState = 'TS186',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction187',
    reactants = ['C7H7O(484)(483)'],
    products = ['C7H7O(486)(485)'],
    transitionState = 'TS187',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond;R_O_C] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond;R_O_C]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction188',
    reactants = ['H(3)(3)', '[C]1C2=CC=CCC1O2(4441)'],
    products = ['C7H7O(484)(483)'],
    transitionState = 'TS188',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction189',
    reactants = ['H(3)(3)', 'C1=CCC2C=C(C=1)O2(4442)'],
    products = ['C7H7O(484)(483)'],
    transitionState = 'TS189',
    kinetics = Arrhenius(A=(1.149e+09,'cm^3/(mol*s)'), n=1.595, Ea=(16.7946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-OneDeH;HJ] for rate rule [Ca_Cds-CdH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction190',
    reactants = ['[C]1=CCC2C=C(C1)O2(4443)'],
    products = ['C7H7O(484)(483)'],
    transitionState = 'TS190',
    kinetics = Arrhenius(A=(1.448e+10,'s^-1'), n=0.82, Ea=(156.9,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 154 used for R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction191',
    reactants = ['[C]1=C2CC=CCC1O2(4444)'],
    products = ['C7H7O(484)(483)'],
    transitionState = 'TS191',
    kinetics = Arrhenius(A=(1.182e+10,'s^-1'), n=0.86, Ea=(149.369,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction192',
    reactants = ['[C]1=CCC2=CC(C1)O2(4445)'],
    products = ['C7H7O(484)(483)'],
    transitionState = 'TS192',
    kinetics = Arrhenius(A=(1.182e+10,'s^-1'), n=0.86, Ea=(149.369,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction193',
    reactants = ['C1=CCC2=C[C](C1)O2(4446)'],
    products = ['C7H7O(484)(483)'],
    transitionState = 'TS193',
    kinetics = Arrhenius(A=(7.12e+06,'s^-1'), n=1.75, Ea=(105.855,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SS(Cd)S;C_rad_out_single;Cs_H_out_H/Cd] for rate rule [R4H_SS(Cd)S;C_rad_out_OneDe/Cs;Cs_H_out_H/Cd]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction194',
    reactants = ['C7H7O(484)(483)'],
    products = ['[CH]1C=CCC2=CC1O2(4447)'],
    transitionState = 'TS194',
    kinetics = Arrhenius(A=(2e+10,'s^-1'), n=0, Ea=(418.4,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;C_rad_out_1H;Cs_H_out_H/NonDeC] for rate rule [R4H_SDS;C_rad_out_H/Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction195',
    reactants = ['C7H7O(484)(483)'],
    products = ['[CH]1C2CC3C=C(O3)C12(4448)'],
    transitionState = 'TS195',
    kinetics = Arrhenius(A=(4.00063e+13,'s^-1'), n=-0.283562, Ea=(537.673,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_beta;doublebond_intra_pri_HNd_Cs;radadd_intra_csHDe] for rate rule [Rn0c7_beta_short;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCd]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic
Ea raised from 533.3 to 537.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction196',
    reactants = ['C=CC1[CH]C12C=CO2(4449)'],
    products = ['C7H7O(484)(483)'],
    transitionState = 'TS196',
    kinetics = Arrhenius(A=(6.21184e+10,'s^-1'), n=0.288169, Ea=(147.049,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_5_unsaturated_hexane] for rate rule [1_5_hexadiene]
Euclidian distance = 1.0
family: 6_membered_central_C-C_shift"""),
)

reaction(
    label = 'reaction197',
    reactants = ['H(3)(3)', 'S(801)(800)'],
    products = ['C7H7O(484)(483)'],
    transitionState = 'TS197',
    kinetics = Arrhenius(A=(93.5525,'m^3/(mol*s)'), n=1.3603, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 20 used for Cds-CsH_Cds-(Cd-Cd-Cd)H;HJ
Exact match found for rate rule [Cds-CsH_Cds-(Cd-Cd-Cd)H;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -1.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction198',
    reactants = ['[C]1=CC2=CC(CC1)O2(4387)'],
    products = ['C7H7O(484)(483)'],
    transitionState = 'TS198',
    kinetics = Arrhenius(A=(1.9054e+11,'s^-1'), n=0.853, Ea=(200.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction199',
    reactants = ['C1=CC2=C[C](CC1)O2(4385)'],
    products = ['C7H7O(484)(483)'],
    transitionState = 'TS199',
    kinetics = Arrhenius(A=(1.28e+07,'s^-1'), n=1.56, Ea=(126.775,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_noH;Cs_H_out_H/Cd] for rate rule [R3H_SS_Cs;C_rad_out_OneDe/O;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction200',
    reactants = ['[C]1=CCCC2C=C1O2(4388)'],
    products = ['C7H7O(484)(483)'],
    transitionState = 'TS200',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_single;Cs_H_out_H/NonDeC] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction201',
    reactants = ['[C]1=C2C=CCCC1O2(4386)'],
    products = ['C7H7O(484)(483)'],
    transitionState = 'TS201',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out_H/(Cd-Cd-Cd)] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction202',
    reactants = ['C7H7O(484)(483)'],
    products = ['[CH]1C2CC3C=CC13O2(4450)'],
    transitionState = 'TS202',
    kinetics = Arrhenius(A=(5.65845e+18,'s^-1'), n=-1.57151, Ea=(387.486,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Rn0c7_gamma_long_SDS_D;doublebond_intra_pri;radadd_intra_csHCd] for rate rule [Rn0c7_gamma_long_SDS_D;doublebond_intra_pri_NdCd;radadd_intra_csHCd]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 384.5 to 387.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction203',
    reactants = ['[CH]1CC=C2C=C(C1)O2(4451)'],
    products = ['C7H7O(484)(483)'],
    transitionState = 'TS203',
    kinetics = Arrhenius(A=(3.94565e+09,'s^-1'), n=0.909333, Ea=(116.555,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [1_3_unsaturated_pentane_backbone;CH(CJ)_1;unsaturated_end] + [1_3_pentadiene;CH_end;unsaturated_end] for rate rule [1_3_pentadiene;CH(CJ)_1;unsaturated_end]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction204',
    reactants = ['[CH2]C1C=CC2=CC1O2(4402)'],
    products = ['C7H7O(484)(483)'],
    transitionState = 'TS204',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction205',
    reactants = ['H(3)(3)', 'O=C1C=C2CC=CC12(4452)'],
    products = ['C7H7O(486)(485)'],
    transitionState = 'TS205',
    kinetics = Arrhenius(A=(89.419,'m^3/(mol*s)'), n=1.598, Ea=(7.0333,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsCs_Cds-OneDeH;HJ] for rate rule [Cds-CsCs_Cds-COH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction206',
    reactants = ['O=C1C[C]2CC=CC21(4453)'],
    products = ['C7H7O(486)(485)'],
    transitionState = 'TS206',
    kinetics = Arrhenius(A=(20108.5,'s^-1'), n=2.606, Ea=(121.63,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_Cs2;Cs_H_out_H/OneDe] for rate rule [R2H_S_cy4;C_rad_out_Cs2_cy5;Cs_H_out_H/CO]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction207',
    reactants = ['C7H7O(486)(485)'],
    products = ['O=C1CC2CC=C[C]12(4454)'],
    transitionState = 'TS207',
    kinetics = Arrhenius(A=(7.59089e+08,'s^-1'), n=1.07, Ea=(157.527,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_H/OneDe;Cs_H_out_noH] for rate rule [R3H_SS_23cy5;C_rad_out_H/OneDe;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction208',
    reactants = ['C7H7O(486)(485)'],
    products = ['O=C1CC2[CH]C=CC12(4455)'],
    transitionState = 'TS208',
    kinetics = Arrhenius(A=(3.8e+10,'s^-1'), n=0.87, Ea=(144.348,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_H/OneDe;Cs_H_out_H/Cd] for rate rule [R3H_SS_12cy4;C_rad_out_H/OneDe;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction209',
    reactants = ['O=C1CC2CC=[C]C12(4456)'],
    products = ['C7H7O(486)(485)'],
    transitionState = 'TS209',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out_H/OneDe] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_H/OneDe]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction210',
    reactants = ['O=C1CC2C[C]=CC12(4457)'],
    products = ['C7H7O(486)(485)'],
    transitionState = 'TS210',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out_H/OneDe] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_H/CO]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction211',
    reactants = ['H(3)(3)', 'O=C1[CH][C]2CC=CC21(4458)'],
    products = ['C7H7O(486)(485)'],
    transitionState = 'TS211',
    kinetics = Arrhenius(A=(4.34601e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction212',
    reactants = ['H(3)(3)', 'O=C1[CH]C2CC=C[C]12(4459)'],
    products = ['C7H7O(486)(485)'],
    transitionState = 'TS212',
    kinetics = Arrhenius(A=(3.62e+13,'cm^3/(mol*s)'), n=0.228, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 13 used for C_rad/TwoDeCs;H_rad
Exact match found for rate rule [C_rad/TwoDeCs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction213',
    reactants = ['H(3)(3)', 'O=C1[CH]C2[CH]C=CC12(4460)'],
    products = ['C7H7O(486)(485)'],
    transitionState = 'TS213',
    kinetics = Arrhenius(A=(2.71464e+07,'m^3/(mol*s)'), n=0.107721, Ea=(5.76381,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 36 used for C_rad/H/CdCs;H_rad
Exact match found for rate rule [C_rad/H/CdCs;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction214',
    reactants = ['H(3)(3)', 'O=C1[CH]C2CC=[C]C12(4461)'],
    products = ['C7H7O(486)(485)'],
    transitionState = 'TS214',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction215',
    reactants = ['H(3)(3)', 'O=C1[CH]C2C[C]=CC12(4462)'],
    products = ['C7H7O(486)(485)'],
    transitionState = 'TS215',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction216',
    reactants = ['C7H7O(486)(485)'],
    products = ['C1=CC2[C]3OC3C2C1(4463)'],
    transitionState = 'TS216',
    kinetics = Arrhenius(A=(9.27213e+10,'s^-1'), n=0.543712, Ea=(183.455,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra_cs] for rate rule [R3_CO;carbonyl_intra_Nd;radadd_intra_csHCs]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction217',
    reactants = ['C7H7O(486)(485)'],
    products = ['O=C1C2[CH]C3CC2C13(4418)'],
    transitionState = 'TS217',
    kinetics = Arrhenius(A=(1.57664e+14,'s^-1'), n=-0.00763077, Ea=(139.328,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5_cyclic;doublebond_intra_pri_HNd;radadd_intra_csHCO] + [R5_cyclic;doublebond_intra_pri_HNd_Cs;radadd_intra_csHDe] for rate rule [Rn1c5_beta_long;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCO]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic
Ea raised from 139.2 to 139.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction218',
    reactants = ['C7H7O(486)(485)'],
    products = ['O=C1C2C3[CH]CC2C13(4464)'],
    transitionState = 'TS218',
    kinetics = Arrhenius(A=(3.615e+11,'s^-1'), n=0.385, Ea=(169.826,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_gamma_long_SSS_D;doublebond_intra_pri_HNd;radadd_intra_csHCO] for rate rule [Rn0c7_gamma_long_SSS_D;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction219',
    reactants = ['O=C=CC1[CH]C=CC1(2038)'],
    products = ['C7H7O(486)(485)'],
    transitionState = 'TS219',
    kinetics = Arrhenius(A=(4.94431e+07,'s^-1'), n=1.16299, Ea=(125.164,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_Cs_RR_D;doublebond_intra;radadd_intra_cs] for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_csHCd]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction220',
    reactants = ['[CH2]C=CC1C=CC1=O(4465)'],
    products = ['C7H7O(486)(485)'],
    transitionState = 'TS220',
    kinetics = Arrhenius(A=(1.26e+11,'s^-1'), n=0.16, Ea=(41.84,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_cyclic;doublebond_intra_pri;radadd_intra_cs2H] for rate rule [Rn3c4_alpha_long;doublebond_intra_pri;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction221',
    reactants = ['O=C1C2[CH]CC=CC12(2913)'],
    products = ['C7H7O(486)(485)'],
    transitionState = 'TS221',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction222',
    reactants = ['O=C1C2[CH]C1CC=C2(4466)'],
    products = ['C7H7O(486)(485)'],
    transitionState = 'TS222',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction223',
    reactants = ['O=[C]C1C2C=CCC12(2089)'],
    products = ['C7H7O(486)(485)'],
    transitionState = 'TS223',
    kinetics = Arrhenius(A=(1.29612e+11,'s^-1'), n=0.58375, Ea=(142.981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction224',
    reactants = ['[CH]1C2=CC1CC=CO2(4467)'],
    products = ['C7H7O(486)(485)'],
    transitionState = 'TS224',
    kinetics = Arrhenius(A=(14080,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond;R_O_C] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond;R_O_C]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction225',
    reactants = ['OC1=CC2CC=C[C]12(2393)'],
    products = ['C7H7O(486)(485)'],
    transitionState = 'TS225',
    kinetics = Arrhenius(A=(1290.48,'s^-1'), n=2.90375, Ea=(139.674,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction226',
    reactants = ['H(3)(3)', 'O=C1[C]C2CC=CC12(4468)'],
    products = ['C7H7O(486)(485)'],
    transitionState = 'TS226',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction227',
    reactants = ['OC1=[C]C2CC=CC12(2395)'],
    products = ['C7H7O(486)(485)'],
    transitionState = 'TS227',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;XH_out] for rate rule [R3H_DS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction228',
    reactants = ['OC1=C[C]2CC=CC21(2392)'],
    products = ['C7H7O(486)(485)'],
    transitionState = 'TS228',
    kinetics = Arrhenius(A=(26.449,'s^-1'), n=2.8625, Ea=(89.0146,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;C_rad_out_single;O_H_out] + [R4H_SS(Cd)S;C_rad_out_single;XH_out] for rate rule [R4H_SS(Cd)S;C_rad_out_OneDe/Cs;O_H_out]
Euclidian distance = 3.16227766017
family: intra_H_migration"""),
)

reaction(
    label = 'reaction229',
    reactants = ['OC1=CC2CC=[C]C12(2396)'],
    products = ['C7H7O(486)(485)'],
    transitionState = 'TS229',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;XH_out] for rate rule [R4H_SS(Cd)S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction230',
    reactants = ['OC1=CC2[CH]C=CC12(2394)'],
    products = ['C7H7O(486)(485)'],
    transitionState = 'TS230',
    kinetics = Arrhenius(A=(0.0508,'s^-1'), n=3.24, Ea=(29.037,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_1H;XH_out] for rate rule [R5H_CCC(Cd);C_rad_out_H/Cd;O_H_out]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction231',
    reactants = ['OC1=CC2C[C]=CC12(2397)'],
    products = ['C7H7O(486)(485)'],
    transitionState = 'TS231',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_single;XH_out] for rate rule [R5H_DSSS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 3.74165738677
family: intra_H_migration"""),
)

reaction(
    label = 'reaction232',
    reactants = ['C7H7O(486)(485)'],
    products = ['[CH]1C2CC3C=C(O2)C13(4469)'],
    transitionState = 'TS232',
    kinetics = Arrhenius(A=(7.87374e+12,'s^-1'), n=-0.166203, Ea=(191.809,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_cyclic;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn1c7_gamma_short;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction233',
    reactants = ['C7H7O(486)(485)'],
    products = ['[CH]1CC2C=C3OC1C32(4403)'],
    transitionState = 'TS233',
    kinetics = Arrhenius(A=(1.26e+11,'s^-1'), n=0.16, Ea=(186.787,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn2c5_alpha_long_SS;doublebond_intra_pri;radadd_intra] for rate rule [Rn2c5_alpha_long_SS;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 186.6 to 186.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction234',
    reactants = ['[O]C12C=CCC1C=C2(2429)'],
    products = ['C7H7O(486)(485)'],
    transitionState = 'TS234',
    kinetics = Arrhenius(A=(6.21184e+10,'s^-1'), n=0.288169, Ea=(147.049,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_5_unsaturated_hexane] for rate rule [1_5_hexadiene]
Euclidian distance = 1.0
family: 6_membered_central_C-C_shift"""),
)

reaction(
    label = 'reaction235',
    reactants = ['C=CC1C2C=CC12[O](4470)'],
    products = ['C7H7O(486)(485)'],
    transitionState = 'TS235',
    kinetics = Arrhenius(A=(6.21184e+10,'s^-1'), n=0.288169, Ea=(147.049,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_5_unsaturated_hexane] for rate rule [1_5_hexadiene]
Euclidian distance = 1.0
family: 6_membered_central_C-C_shift"""),
)

reaction(
    label = 'reaction236',
    reactants = ['O(4)(4)', '[C]1=CC2CC=CC12(1432)'],
    products = ['C7H7O(486)(485)'],
    transitionState = 'TS236',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction237',
    reactants = ['H(3)(3)', 'OC1C=CC=C=CC=1(4471)'],
    products = ['C7H7O(495)(494)'],
    transitionState = 'TS237',
    kinetics = Arrhenius(A=(1.149e+09,'cm^3/(mol*s)'), n=1.595, Ea=(16.7946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-OneDeH;HJ] for rate rule [Ca_Cds-CdH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction238',
    reactants = ['H(3)(3)', 'OC1C=C=CC=CC=1(4472)'],
    products = ['C7H7O(495)(494)'],
    transitionState = 'TS238',
    kinetics = Arrhenius(A=(1.149e+09,'cm^3/(mol*s)'), n=1.595, Ea=(16.7946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-OneDeH;HJ] for rate rule [Ca_Cds-CdH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction239',
    reactants = ['OC1C=CC[C]=CC=1(2214)'],
    products = ['C7H7O(495)(494)'],
    transitionState = 'TS239',
    kinetics = Arrhenius(A=(2.3e+10,'s^-1'), n=0.98, Ea=(112.55,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(Cd-Cd-Cd)] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(Cd-Cd-Cd-Cd-Cd)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction240',
    reactants = ['OC1C=[C]CC=CC=1(2213)'],
    products = ['C7H7O(495)(494)'],
    transitionState = 'TS240',
    kinetics = Arrhenius(A=(2.3e+10,'s^-1'), n=0.98, Ea=(112.55,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(Cd-Cd-Cd)] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(Cd-Cd-Cd-Cd-Cd)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction241',
    reactants = ['OC1C=CCC=[C]C=1(2217)'],
    products = ['C7H7O(495)(494)'],
    transitionState = 'TS241',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 3.60555127546
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction242',
    reactants = ['OC1[C]=CCC=CC=1(2215)'],
    products = ['C7H7O(495)(494)'],
    transitionState = 'TS242',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 3.60555127546
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction243',
    reactants = ['OC1=[C]C=CCC=C1(2216)'],
    products = ['C7H7O(495)(494)'],
    transitionState = 'TS243',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H;Cd_rad_out;Cs_H_out_H/(Cd-Cd-Cd)] for rate rule [R4H_SDS;Cd_rad_out_Cd;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 2.82842712475
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction244',
    reactants = ['OH(5)(5)', '[C]1=CC=C[CH]C=C1(4473)'],
    products = ['C7H7O(495)(494)'],
    transitionState = 'TS244',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Cd_rad/Cd] for rate rule [O_pri_rad;Cd_rad/Cd]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction245',
    reactants = ['H(3)(3)', 'C7H6O(488)(487)'],
    products = ['C7H7O(495)(494)'],
    transitionState = 'TS245',
    kinetics = Arrhenius(A=(7e+11,'cm^3/(mol*s)'), n=0.493, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2500,'K'), comment="""From training reaction 25 used for H_rad;O_rad/OneDe
Exact match found for rate rule [H_rad;O_rad/OneDe]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction246',
    reactants = ['H(3)(3)', 'OC1[C]=C[CH]C=CC=1(4474)'],
    products = ['C7H7O(495)(494)'],
    transitionState = 'TS246',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction247',
    reactants = ['H(3)(3)', 'OC1C=C[CH][C]=CC=1(4475)'],
    products = ['C7H7O(495)(494)'],
    transitionState = 'TS247',
    kinetics = Arrhenius(A=(4.34601e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction248',
    reactants = ['H(3)(3)', 'OC1C=[C][CH]C=CC=1(4476)'],
    products = ['C7H7O(495)(494)'],
    transitionState = 'TS248',
    kinetics = Arrhenius(A=(4.34601e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction249',
    reactants = ['C7H7O(495)(494)'],
    products = ['OC1C=CC2[CH]C2C=1(4361)'],
    transitionState = 'TS249',
    kinetics = Arrhenius(A=(4.45714e+12,'s^-1'), n=0, Ea=(182.004,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 127 used for Rn0c7_beta_short_SDSDS_D;doublebond_intra_pri_HCd;radadd_intra_csH(CdCdCd)
Exact match found for rate rule [Rn0c7_beta_short_SDSDS_D;doublebond_intra_pri_HCd;radadd_intra_csH(CdCdCd)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction250',
    reactants = ['C7H7O(495)(494)'],
    products = ['OC1=CC=CC2[CH]C12(4477)'],
    transitionState = 'TS250',
    kinetics = Arrhenius(A=(4.45714e+12,'s^-1'), n=0, Ea=(182.004,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 127 used for Rn0c7_beta_short_SDSDS_D;doublebond_intra_pri_HCd;radadd_intra_csH(CdCdCd)
Exact match found for rate rule [Rn0c7_beta_short_SDSDS_D;doublebond_intra_pri_HCd;radadd_intra_csH(CdCdCd)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction251',
    reactants = ['C7H7O(495)(494)'],
    products = ['OC1=C[CH]C2C=CC12(2462)'],
    transitionState = 'TS251',
    kinetics = Arrhenius(A=(1.13169e+19,'s^-1'), n=-1.57151, Ea=(227.049,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Rn0c7_gamma_long_SDS_D;doublebond_intra;radadd_intra_csHCd] for rate rule [Rn0c7_gamma_long_SDS_D;doublebond_intra_secNd_HCd;radadd_intra_csHCd]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction252',
    reactants = ['C7H7O(495)(494)'],
    products = ['OC12[CH]C=CC1C=C2(2353)'],
    transitionState = 'TS252',
    kinetics = Arrhenius(A=(1.13169e+19,'s^-1'), n=-1.57151, Ea=(227.049,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Rn0c7_gamma_long_SDS_D;doublebond_intra_pri;radadd_intra_csHCd] for rate rule [Rn0c7_gamma_long_SDS_D;doublebond_intra_pri_NdCd;radadd_intra_csHCd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction253',
    reactants = ['C7H7O(495)(494)'],
    products = ['OC1C=[C]C=CCC=1(2296)'],
    transitionState = 'TS253',
    kinetics = Arrhenius(A=(7.85714e+07,'s^-1'), n=1.56, Ea=(259.91,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for 1_3_pentadiene;CH=C_1;CdCJ_2
Exact match found for rate rule [1_3_pentadiene;CH=C_1;CdCJ_2]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction254',
    reactants = ['C7H7O(495)(494)'],
    products = ['O=C1C=C[CH]C=CC1(2175)'],
    transitionState = 'TS254',
    kinetics = Arrhenius(A=(2580.97,'s^-1'), n=2.90375, Ea=(139.674,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond;R_O_H] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond;R_O_H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction255',
    reactants = ['H(3)(3)', 'OC1=C=CC=CC=C1(4478)'],
    products = ['C7H7O(495)(494)'],
    transitionState = 'TS255',
    kinetics = Arrhenius(A=(1.149e+09,'cm^3/(mol*s)'), n=1.595, Ea=(16.7946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-OneDeH;HJ] for rate rule [Ca_Cds-CdH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction256',
    reactants = ['OC1=[C]CC=CC=C1(2292)'],
    products = ['C7H7O(495)(494)'],
    transitionState = 'TS256',
    kinetics = Arrhenius(A=(2.3e+10,'s^-1'), n=0.98, Ea=(112.55,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(Cd-Cd-Cd)] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(Cd-Cd-Cd-Cd-Cd)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction257',
    reactants = ['OC1C=CC=[C]CC=1(2293)'],
    products = ['C7H7O(495)(494)'],
    transitionState = 'TS257',
    kinetics = Arrhenius(A=(2.3e+10,'s^-1'), n=0.98, Ea=(112.55,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(Cd-Cd-Cd)] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(Cd-Cd-Cd-Cd-Cd)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction258',
    reactants = ['OC1C=C[C]=CCC=1(2294)'],
    products = ['C7H7O(495)(494)'],
    transitionState = 'TS258',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 3.60555127546
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction259',
    reactants = ['[O]C1C=CC=CCC=1(2291)'],
    products = ['C7H7O(495)(494)'],
    transitionState = 'TS259',
    kinetics = Arrhenius(A=(2.08e+07,'s^-1'), n=1.61, Ea=(113.386,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;Y_rad_out;Cs_H_out_H/(Cd-Cd-Cd)] for rate rule [R4H_SDS;O_rad_out;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction260',
    reactants = ['OC1[C]=CC=CCC=1(2295)'],
    products = ['C7H7O(495)(494)'],
    transitionState = 'TS260',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H;Cd_rad_out;Cs_H_out_H/(Cd-Cd-Cd)] for rate rule [R4H_SDS;Cd_rad_out_Cd;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 2.82842712475
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction261',
    reactants = ['C7H7O(495)(494)'],
    products = ['OC12[CH]C1C=CC=C2(4479)'],
    transitionState = 'TS261',
    kinetics = Arrhenius(A=(4.45714e+12,'s^-1'), n=0, Ea=(182.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_beta_short_SDSDS_D;doublebond_intra_pri;radadd_intra_csH(CdCdCd)] for rate rule [Rn0c7_beta_short_SDSDS_D;doublebond_intra_pri_NdCd;radadd_intra_csH(CdCdCd)]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction262',
    reactants = ['C7H7O(495)(494)'],
    products = ['OC1[CH]C2C=CC2C=1(2316)'],
    transitionState = 'TS262',
    kinetics = Arrhenius(A=(5.65845e+18,'s^-1'), n=-1.57151, Ea=(227.049,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 131 used for Rn0c7_gamma_long_SDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCd
Exact match found for rate rule [Rn0c7_gamma_long_SDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCd]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction263',
    reactants = ['C7H7O(495)(494)'],
    products = ['OC1=CC2[CH]C=CC12(2394)'],
    transitionState = 'TS263',
    kinetics = Arrhenius(A=(1.13169e+19,'s^-1'), n=-1.57151, Ea=(227.049,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 131 used for Rn0c7_gamma_long_SDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCd
Exact match found for rate rule [Rn0c7_gamma_long_SDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction264',
    reactants = ['C7H7O(495)(494)'],
    products = ['OC1=CC=[C]C=CC1(2268)'],
    transitionState = 'TS264',
    kinetics = Arrhenius(A=(7.85714e+07,'s^-1'), n=1.56, Ea=(259.91,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for 1_3_pentadiene;CH=C_1;CdCJ_2
Exact match found for rate rule [1_3_pentadiene;CH=C_1;CdCJ_2]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction265',
    reactants = ['OH(5)(5)', 'C1=CC=CC=CC=1(4480)'],
    products = ['C7H7O(495)(494)'],
    transitionState = 'TS265',
    kinetics = Arrhenius(A=(41610,'cm^3/(mol*s)'), n=2.487, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds;OJ_pri] for rate rule [Ca_Cds-CdH;OJ_pri]
Euclidian distance = 2.0
family: R_Addition_MultipleBond
Ea raised from -7.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction266',
    reactants = ['OC1=CC=CC=[C]C1(2265)'],
    products = ['C7H7O(495)(494)'],
    transitionState = 'TS266',
    kinetics = Arrhenius(A=(2.3e+10,'s^-1'), n=0.98, Ea=(112.55,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(Cd-Cd-Cd)] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(Cd-Cd-Cd-Cd-Cd)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction267',
    reactants = ['OC1=[C]C=CC=CC1(2266)'],
    products = ['C7H7O(495)(494)'],
    transitionState = 'TS267',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 3.60555127546
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction268',
    reactants = ['OC1=CC=C[C]=CC1(2267)'],
    products = ['C7H7O(495)(494)'],
    transitionState = 'TS268',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 3.60555127546
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction269',
    reactants = ['OC1=C[C]=CC=CC1(2269)'],
    products = ['C7H7O(495)(494)'],
    transitionState = 'TS269',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H;Cd_rad_out;Cs_H_out_H/(Cd-Cd-Cd)] for rate rule [R4H_SDS;Cd_rad_out_Cd;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 2.82842712475
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction270',
    reactants = ['C7H7O(495)(494)'],
    products = ['O[C]1C2C=CC=CC12(4481)'],
    transitionState = 'TS270',
    kinetics = Arrhenius(A=(2.22857e+12,'s^-1'), n=0, Ea=(182.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_beta_short_SDSDS_D;doublebond_intra;radadd_intra_csH(CdCdCd)] for rate rule [Rn0c7_beta_short_SDSDS_D;doublebond_intra_secNd_HCd;radadd_intra_csH(CdCdCd)]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction271',
    reactants = ['C7H7O(495)(494)'],
    products = ['OC1C=C[C]=CC=C1(2246)'],
    transitionState = 'TS271',
    kinetics = Arrhenius(A=(7.85714e+07,'s^-1'), n=1.56, Ea=(259.91,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for 1_3_pentadiene;CH=C_1;CdCJ_2
Exact match found for rate rule [1_3_pentadiene;CH=C_1;CdCJ_2]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction272',
    reactants = ['H(3)(3)', 'C7H6O(492)(491)'],
    products = ['C7H7O(495)(494)'],
    transitionState = 'TS272',
    kinetics = Arrhenius(A=(1.15512,'m^3/(mol*s)'), n=2.412, Ea=(29.0596,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Od_CO;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction273',
    reactants = ['[O]C1C=CC=CC=C1(1734)'],
    products = ['C7H7O(495)(494)'],
    transitionState = 'TS273',
    kinetics = Arrhenius(A=(4.73428e+08,'s^-1'), n=1.24562, Ea=(97.1472,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_CdCd] for rate rule [R2H_S;O_rad_out;Cs_H_out_CdCd]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction274',
    reactants = ['OC1[C]=CC=CC=C1(2244)'],
    products = ['C7H7O(495)(494)'],
    transitionState = 'TS274',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_(CdCdCd)]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction275',
    reactants = ['OC1C=[C]C=CC=C1(2245)'],
    products = ['C7H7O(495)(494)'],
    transitionState = 'TS275',
    kinetics = Arrhenius(A=(4.70399e+08,'s^-1'), n=1.21625, Ea=(125.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_single;Cs_H_out_noH] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_(CdCdCd)]
Euclidian distance = 3.60555127546
family: intra_H_migration"""),
)

network(
    label = '109',
    isomers = [
        'C7H7O(470)(469)',
        'C7H7O(473)(472)',
        'C7H7O(476)(475)',
        'C7H7O(477)(476)',
        'C7H7O(478)(477)',
        'C7H7O(481)(480)',
        'C7H7O(484)(483)',
        'C7H7O(486)(485)',
        'C7H7O(495)(494)',
    ],
    reactants = [
        ('H(3)(3)', 'C7H6O(492)(491)'),
        ('H(3)(3)', 'C7H6O(488)(487)'),
        ('H(3)(3)', 'C7H6O(487)(486)'),
        ('H(3)(3)', 'S(801)(800)'),
        ('H(3)(3)', 'C7H6O(535)(534)'),
    ],
    bathGas = {
        'Ne': 0.333333,
        'N2': 0.333333,
        'Ar(8)': 0.333333,
    },
)

pressureDependence(
    label = '109',
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

