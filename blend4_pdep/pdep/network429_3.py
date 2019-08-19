species(
    label = 'C6H7O(505)(504)',
    structure = SMILES('[O]C1C=CC=CC1'),
    E0 = (145.927,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4135.17,'J/mol'), sigma=(6.71154,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=645.90 K, Pc=31.04 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94943,0.0236292,7.08526e-05,-1.10511e-07,4.40926e-11,17644.2,19.7056], Tmin=(100,'K'), Tmax=(961.861,'K')), NASAPolynomial(coeffs=[15.9872,0.0191697,-6.27642e-06,1.2259e-09,-9.64763e-14,12449.5,-60.438], Tmin=(961.861,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.927,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(CC(C)OJ)"""),
)

species(
    label = 'C6H7O(506)(505)',
    structure = SMILES('[CH2]C=CC=CC=O'),
    E0 = (52.4633,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.09775,'amu*angstrom^2'), symmetry=1, barrier=(25.2395,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09742,'amu*angstrom^2'), symmetry=1, barrier=(25.2318,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09757,'amu*angstrom^2'), symmetry=1, barrier=(25.2352,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3745.2,'J/mol'), sigma=(6.02434,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=584.99 K, Pc=38.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19817,0.051949,-2.18837e-05,-5.89878e-09,4.90491e-12,6418.88,24.4996], Tmin=(100,'K'), Tmax=(1111.19,'K')), NASAPolynomial(coeffs=[13.2992,0.0262162,-1.12126e-05,2.13734e-09,-1.51482e-13,2628.95,-40.109], Tmin=(1111.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(52.4633,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C6H7O(511)(510)',
    structure = SMILES('[CH]1C=CC=COC1'),
    E0 = (52.2376,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3904.52,'J/mol'), sigma=(6.47878,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=609.88 K, Pc=32.58 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29188,-0.0121971,0.000289705,-4.29403e-07,1.82772e-10,6447.86,17.7324], Tmin=(100,'K'), Tmax=(899.402,'K')), NASAPolynomial(coeffs=[51.7645,-0.047906,3.44457e-05,-6.84529e-09,4.54201e-13,-10265.9,-262.861], Tmin=(899.402,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(52.2376,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(C=CCJCO)"""),
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
    label = 'C6H7(464)(463)',
    structure = SMILES('[CH]1C=CC=CC1'),
    E0 = (185.261,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3639.77,'J/mol'), sigma=(6.13475,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=568.52 K, Pc=35.77 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.41084,0.0238203,3.51346e-05,-5.30671e-08,1.90206e-11,22348.5,14.4651], Tmin=(100,'K'), Tmax=(1051.47,'K')), NASAPolynomial(coeffs=[8.26331,0.0294392,-1.26583e-05,2.45535e-09,-1.76993e-13,19576.4,-21.3929], Tmin=(1051.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(185.261,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(Aromatic_pi_S_1_3)"""),
)

species(
    label = 'CO(10)(11)',
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
    label = 'C5H7(211)(210)',
    structure = SMILES('[CH2]C=CC=C'),
    E0 = (175.951,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.52628,'amu*angstrom^2'), symmetry=1, barrier=(35.0921,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.5297,'amu*angstrom^2'), symmetry=1, barrier=(35.1709,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3140.68,'J/mol'), sigma=(5.4037,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=490.57 K, Pc=45.16 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2987,0.0238036,3.93739e-05,-6.93177e-08,2.88419e-11,21235.8,16.7723], Tmin=(100,'K'), Tmax=(942.053,'K')), NASAPolynomial(coeffs=[11.9106,0.0173605,-5.0926e-06,8.77895e-10,-6.40305e-14,17899.7,-37.1203], Tmin=(942.053,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(175.951,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ)"""),
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
    label = 'C6H6O(516)(515)',
    structure = SMILES('C=CC=CC=C=O'),
    E0 = (80.186,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2120,512.5,787.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.925634,'amu*angstrom^2'), symmetry=1, barrier=(21.2821,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.924937,'amu*angstrom^2'), symmetry=1, barrier=(21.2661,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3245.38,'J/mol'), sigma=(5.78871,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=506.92 K, Pc=37.96 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11369,0.0568433,-4.67069e-05,1.98587e-08,-3.39265e-12,9753.53,23.3106], Tmin=(100,'K'), Tmax=(1397.17,'K')), NASAPolynomial(coeffs=[13.6724,0.0208885,-8.1057e-06,1.43979e-09,-9.68877e-14,6244.24,-41.4778], Tmin=(1397.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(80.186,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CCO)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'O=C1C=CC=CC1(4653)',
    structure = SMILES('O=C1C=CC=CC1'),
    E0 = (-20.0282,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.16162,0.027164,3.54173e-05,-6.07398e-08,2.36347e-11,-2330.95,18.4367], Tmin=(100,'K'), Tmax=(999.903,'K')), NASAPolynomial(coeffs=[10.8449,0.0248769,-9.83038e-06,1.88381e-09,-1.37388e-13,-5689.57,-31.5655], Tmin=(999.903,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-20.0282,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-Cds(Cds-Cds)H) + ring(24cyclohexadienone)"""),
)

species(
    label = '[O]C1[CH]C=CC=C1(4660)',
    structure = SMILES('[O]C1[CH]C=CC=C1'),
    E0 = (216.07,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.20833,0.00941517,0.000120939,-1.70937e-07,6.75914e-11,26079.1,19.5746], Tmin=(100,'K'), Tmax=(947.562,'K')), NASAPolynomial(coeffs=[20.0045,0.0106173,-1.78988e-06,4.17835e-10,-4.60759e-14,19280,-83.4042], Tmin=(947.562,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.07,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(C=CCJC(O)C=C) + radical(CC(C)OJ)"""),
)

species(
    label = '[O][C]1C=CC=CC1(4661)',
    structure = SMILES('[O]C1=CC=C[CH]C1'),
    E0 = (108.507,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75406,0.0388812,3.74073e-06,-2.90796e-08,1.26739e-11,13140.3,19.0714], Tmin=(100,'K'), Tmax=(1035.72,'K')), NASAPolynomial(coeffs=[11.5253,0.0248141,-1.01666e-05,1.93753e-09,-1.38982e-13,9846.73,-34.5408], Tmin=(1035.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(108.507,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(Aromatic_pi_S_1_3) + radical(C=C(C)OJ)"""),
)

species(
    label = '[O]C1C=CC=[C]C1(4662)',
    structure = SMILES('[O]C1C=CC=[C]C1'),
    E0 = (383.768,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91824,0.0282077,4.61347e-05,-8.13173e-08,3.31492e-11,46247.4,20.2913], Tmin=(100,'K'), Tmax=(969.756,'K')), NASAPolynomial(coeffs=[14.8083,0.018652,-6.54398e-06,1.27254e-09,-9.76094e-14,41696.7,-52.0734], Tmin=(969.756,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(383.768,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = '[O]C1[C]=CC=CC1(4663)',
    structure = SMILES('[O]C1[C]=CC=CC1'),
    E0 = (383.768,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91824,0.0282077,4.61347e-05,-8.13173e-08,3.31492e-11,46247.4,20.2913], Tmin=(100,'K'), Tmax=(969.756,'K')), NASAPolynomial(coeffs=[14.8083,0.018652,-6.54398e-06,1.27254e-09,-9.76094e-14,41696.7,-52.0734], Tmin=(969.756,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(383.768,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(Cds_S) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]C1C=C[C]=CC1(4664)',
    structure = SMILES('[O]C1C=C[C]=CC1'),
    E0 = (344.922,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86607,0.03024,4.15034e-05,-7.72697e-08,3.21539e-11,41576.5,19.6012], Tmin=(100,'K'), Tmax=(958.733,'K')), NASAPolynomial(coeffs=[14.4673,0.0192213,-6.27408e-06,1.16319e-09,-8.74404e-14,37250.4,-50.6215], Tmin=(958.733,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(344.922,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(C=CJC=C) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]C1C=[C]C=CC1(4665)',
    structure = SMILES('[O]C1C=[C]C=CC1'),
    E0 = (344.922,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86607,0.03024,4.15034e-05,-7.72697e-08,3.21539e-11,41576.5,19.6012], Tmin=(100,'K'), Tmax=(958.733,'K')), NASAPolynomial(coeffs=[14.4673,0.0192213,-6.27408e-06,1.16319e-09,-8.74404e-14,37250.4,-50.6215], Tmin=(958.733,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(344.922,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(C=CJC=C) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]C1=CC=CCC1(4670)',
    structure = SMILES('O=C1C=C[CH]CC1'),
    E0 = (1.79825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.49797,0.0192553,6.14756e-05,-7.69355e-08,2.53392e-11,281.813,15.044], Tmin=(100,'K'), Tmax=(1100.97,'K')), NASAPolynomial(coeffs=[7.67423,0.039335,-1.8861e-05,3.79079e-09,-2.76091e-13,-3214.72,-21.1294], Tmin=(1100.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1.79825,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + ring(Cyclohexane) + radical(cyclohexene-allyl)"""),
)

species(
    label = '[O]C1CC2C=CC12(4671)',
    structure = SMILES('[O]C1CC2C=CC12'),
    E0 = (167.048,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95588,0.0288272,4.51048e-05,-7.58376e-08,2.98392e-11,20179.2,17.8447], Tmin=(100,'K'), Tmax=(987.682,'K')), NASAPolynomial(coeffs=[12.6741,0.0251541,-9.66195e-06,1.86056e-09,-1.37518e-13,16123.9,-43.5427], Tmin=(987.682,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(167.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_4_ene_1) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=CC=CCC=O(4652)',
    structure = SMILES('[CH]=CC=CCC=O'),
    E0 = (205.729,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.878927,'amu*angstrom^2'), symmetry=1, barrier=(20.2083,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.881595,'amu*angstrom^2'), symmetry=1, barrier=(20.2696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.880396,'amu*angstrom^2'), symmetry=1, barrier=(20.242,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17864,0.0494311,-9.2841e-06,-2.22687e-08,1.10615e-11,24855.9,24.8368], Tmin=(100,'K'), Tmax=(1068.93,'K')), NASAPolynomial(coeffs=[15.325,0.0234295,-1.05938e-05,2.12109e-09,-1.55935e-13,20292.8,-51.5523], Tmin=(1068.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.729,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'O[C]1C=CC=CC1(4654)',
    structure = SMILES('OC1=CC=C[CH]C1'),
    E0 = (-29.2978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40527,0.0426762,1.17477e-05,-4.65941e-08,2.07835e-11,-3417.52,18.8205], Tmin=(100,'K'), Tmax=(991.211,'K')), NASAPolynomial(coeffs=[14.9978,0.0218927,-8.35632e-06,1.60277e-09,-1.18211e-13,-7785.73,-55.0776], Tmin=(991.211,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-29.2978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(Aromatic_pi_S_1_3)"""),
)

species(
    label = 'OC1[CH]C=CC=C1(4655)',
    structure = SMILES('OC1[CH]C=CC=C1'),
    E0 = (-14.2907,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.04794,0.0133802,0.000115135,-1.66038e-07,6.60464e-11,-1621.48,20.3421], Tmin=(100,'K'), Tmax=(946.91,'K')), NASAPolynomial(coeffs=[20.0651,0.0124174,-2.37964e-06,5.06204e-10,-5.11314e-14,-8402.56,-83.3865], Tmin=(946.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-14.2907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = 'OC1[C]=CC=CC1(4656)',
    structure = SMILES('OC1[C]=CC=CC1'),
    E0 = (153.408,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7566,0.0321887,4.02698e-05,-7.6331e-08,3.15628e-11,18546.9,21.0633], Tmin=(100,'K'), Tmax=(969.383,'K')), NASAPolynomial(coeffs=[14.8707,0.0204488,-7.1317e-06,1.36041e-09,-1.02622e-13,14013.4,-52.0661], Tmin=(969.383,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(153.408,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(Cds_S)"""),
)

species(
    label = 'OC1C=CC=[C]C1(4657)',
    structure = SMILES('OC1C=CC=[C]C1'),
    E0 = (153.408,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7566,0.0321887,4.02698e-05,-7.6331e-08,3.15628e-11,18546.9,21.0633], Tmin=(100,'K'), Tmax=(969.383,'K')), NASAPolynomial(coeffs=[14.8707,0.0204488,-7.1317e-06,1.36041e-09,-1.02622e-13,14013.4,-52.0661], Tmin=(969.383,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(153.408,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(Cds_S)"""),
)

species(
    label = 'OC1C=[C]C=CC1(4658)',
    structure = SMILES('OC1C=[C]C=CC1'),
    E0 = (114.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70538,0.0342091,3.56826e-05,-7.23447e-08,3.05955e-11,13875.9,20.3699], Tmin=(100,'K'), Tmax=(957.873,'K')), NASAPolynomial(coeffs=[14.5273,0.0210222,-6.86425e-06,1.25165e-09,-9.25028e-14,9568.13,-50.6008], Tmin=(957.873,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(114.561,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(C=CJC=C)"""),
)

species(
    label = 'OC1C=C[C]=CC1(4659)',
    structure = SMILES('OC1C=C[C]=CC1'),
    E0 = (114.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70538,0.0342091,3.56826e-05,-7.23447e-08,3.05955e-11,13875.9,20.3699], Tmin=(100,'K'), Tmax=(957.873,'K')), NASAPolynomial(coeffs=[14.5273,0.0210222,-6.86425e-06,1.25165e-09,-9.25028e-14,9568.13,-50.6008], Tmin=(957.873,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(114.561,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]1C2C=CCC1O2(4666)',
    structure = SMILES('[CH]1C2C=CCC1O2'),
    E0 = (217.99,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.46948,0.00950531,0.000105562,-1.41029e-07,5.31859e-11,26295,18.2578], Tmin=(100,'K'), Tmax=(973.034,'K')), NASAPolynomial(coeffs=[14.4103,0.0218082,-8.04048e-06,1.64461e-09,-1.3006e-13,21065.1,-53.9568], Tmin=(973.034,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(217.99,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s3_4_6_ene_1) + radical(CCJCO)"""),
)

species(
    label = '[CH]1CC2C=CC1O2(4667)',
    structure = SMILES('[CH]1CC2C=CC1O2'),
    E0 = (177.736,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.58418,0.0033332,0.000127697,-1.6681e-07,6.28624e-11,21452.8,18.8459], Tmin=(100,'K'), Tmax=(968.197,'K')), NASAPolynomial(coeffs=[15.9131,0.0191945,-6.76291e-06,1.43882e-09,-1.18953e-13,15547.4,-62.1957], Tmin=(968.197,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(177.736,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s3_5_5_ene_1) + radical(CCJCO)"""),
)

species(
    label = '[CH]1C=CC2CC1O2(4668)',
    structure = SMILES('[CH]1C=CC2CC1O2'),
    E0 = (135.005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06208,0.0185995,8.95049e-05,-1.28207e-07,4.92156e-11,16328.6,13.7301], Tmin=(100,'K'), Tmax=(978.411,'K')), NASAPolynomial(coeffs=[16.0308,0.0219935,-8.45358e-06,1.74134e-09,-1.37178e-13,10699.3,-68.1551], Tmin=(978.411,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(135.005,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s3_4_6_ene_1) + radical(C=CCJCO)"""),
)

species(
    label = '[CH]1C=CCC2OC12(4669)',
    structure = SMILES('[CH]1C=CCC2OC12'),
    E0 = (111.074,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76755,0.0319978,4.33603e-05,-7.50741e-08,2.92382e-11,13454.5,14.2132], Tmin=(100,'K'), Tmax=(1007.51,'K')), NASAPolynomial(coeffs=[13.6129,0.0267965,-1.11692e-05,2.21431e-09,-1.64761e-13,8944.8,-53.5577], Tmin=(1007.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(111.074,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_3_6_ene_2) + radical(C=CCJCO)"""),
)

species(
    label = 'O=C[CH]C1C=CC1(4610)',
    structure = SMILES('[O]C=CC1C=CC1'),
    E0 = (148.839,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43332,0.0373961,3.73022e-05,-8.14725e-08,3.55511e-11,18011,21.4209], Tmin=(100,'K'), Tmax=(949.288,'K')), NASAPolynomial(coeffs=[17.8061,0.0161207,-4.47507e-06,8.15567e-10,-6.41713e-14,12752.6,-68.0401], Tmin=(949.288,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(148.839,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(C=COJ)"""),
)

species(
    label = 'C=C=CC=CC=O(4611)',
    structure = SMILES('C=C=CC=CC=O'),
    E0 = (111.012,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,610,2055,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.07836,'amu*angstrom^2'), symmetry=1, barrier=(24.7936,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.075,'amu*angstrom^2'), symmetry=1, barrier=(24.7163,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14219,0.0583098,-4.64171e-05,1.80867e-08,-2.82455e-12,13458,22.1956], Tmin=(100,'K'), Tmax=(1503.78,'K')), NASAPolynomial(coeffs=[14.912,0.0216827,-9.88212e-06,1.88985e-09,-1.31862e-13,9316.65,-49.8535], Tmin=(1503.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(111.012,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C[C]=CC=CC=O(4612)',
    structure = SMILES('C[C]=CC=CC=O'),
    E0 = (172.249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5022,0.0583251,-4.56188e-05,1.85117e-08,-3.18513e-12,20803.6,22.3142], Tmin=(100,'K'), Tmax=(1307.14,'K')), NASAPolynomial(coeffs=[10.1915,0.0317353,-1.51063e-05,2.94991e-09,-2.08877e-13,18531.9,-21.934], Tmin=(1307.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(172.249,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_S)"""),
)

species(
    label = 'CC=[C]C=CC=O(4613)',
    structure = SMILES('CC=C=CC=C[O]'),
    E0 = (131.144,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,540,610,2055,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.15959,'amu*angstrom^2'), symmetry=1, barrier=(26.6612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15729,'amu*angstrom^2'), symmetry=1, barrier=(26.6083,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.76895,0.0569259,-1.5408e-05,-3.02192e-08,1.83618e-11,15902.2,23.5269], Tmin=(100,'K'), Tmax=(947.117,'K')), NASAPolynomial(coeffs=[18.9353,0.0150296,-4.21071e-06,7.23116e-10,-5.36501e-14,10899.1,-71.3742], Tmin=(947.117,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(131.144,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ)"""),
)

species(
    label = 'CC=C[C]=CC=O(4614)',
    structure = SMILES('CC=C[C]=CC=O'),
    E0 = (133.403,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180],'cm^-1')),
        HinderedRotor(inertia=(0.910243,'amu*angstrom^2'), symmetry=1, barrier=(20.9283,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.910657,'amu*angstrom^2'), symmetry=1, barrier=(20.9378,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.910309,'amu*angstrom^2'), symmetry=1, barrier=(20.9298,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38533,0.0610028,-5.21016e-05,2.46302e-08,-5.00837e-12,16135.9,21.866], Tmin=(100,'K'), Tmax=(1127.17,'K')), NASAPolynomial(coeffs=[9.12659,0.0335312,-1.55432e-05,3.00753e-09,-2.12563e-13,14390.7,-16.4078], Tmin=(1127.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(133.403,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CJC=C)"""),
)

species(
    label = 'CC=CC=[C]C=O(4615)',
    structure = SMILES('CC=CC=C=C[O]'),
    E0 = (131.144,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,540,610,2055,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.15959,'amu*angstrom^2'), symmetry=1, barrier=(26.6612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15729,'amu*angstrom^2'), symmetry=1, barrier=(26.6083,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.76895,0.0569259,-1.5408e-05,-3.02192e-08,1.83618e-11,15902.2,23.5269], Tmin=(100,'K'), Tmax=(947.117,'K')), NASAPolynomial(coeffs=[18.9353,0.0150296,-4.21071e-06,7.23116e-10,-5.36501e-14,10899.1,-71.3742], Tmin=(947.117,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(131.144,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ)"""),
)

species(
    label = 'CC=CC=C[C]=O(4616)',
    structure = SMILES('CC=C[CH]C=C=O'),
    E0 = (67.3711,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,626.033],'cm^-1')),
        HinderedRotor(inertia=(0.784119,'amu*angstrom^2'), symmetry=1, barrier=(18.0284,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.61496,'amu*angstrom^2'), symmetry=1, barrier=(37.1311,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.551652,'amu*angstrom^2'), symmetry=1, barrier=(12.6836,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02962,0.0580588,-4.29012e-05,1.62743e-08,-2.50267e-12,8215.76,23.0833], Tmin=(100,'K'), Tmax=(1525.46,'K')), NASAPolynomial(coeffs=[13.9687,0.0241309,-9.53981e-06,1.69461e-09,-1.13308e-13,4268.13,-44.8043], Tmin=(1525.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(67.3711,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-Cd(CCO)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + radical(C=CCJC=C=O)"""),
)

species(
    label = '[CH]=C[CH2](924)',
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
    label = 'CHCHCHO(1679)',
    structure = SMILES('[CH]=CC=O'),
    E0 = (171.79,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,688.575,691.322,692.403,693.509,695.331],'cm^-1')),
        HinderedRotor(inertia=(0.0100227,'amu*angstrom^2'), symmetry=1, barrier=(3.38495,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (55.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.95468,0.0145151,1.89382e-05,-3.59857e-08,1.47962e-11,20706.8,12.171], Tmin=(100,'K'), Tmax=(981.049,'K')), NASAPolynomial(coeffs=[9.63601,0.00810119,-3.09994e-06,6.30252e-10,-4.9152e-14,18393.6,-25.043], Tmin=(981.049,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(171.79,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CHCHCHO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C=[C]C=CC=O(4617)',
    structure = SMILES('[CH2]C=[C]C=CC=O'),
    E0 = (251.459,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.48426,'amu*angstrom^2'), symmetry=1, barrier=(34.1261,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48843,'amu*angstrom^2'), symmetry=1, barrier=(34.222,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4825,'amu*angstrom^2'), symmetry=1, barrier=(34.0857,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32441,0.0561728,-4.32064e-05,1.72705e-08,-2.8472e-12,30341.8,23.6367], Tmin=(100,'K'), Tmax=(1405.2,'K')), NASAPolynomial(coeffs=[11.8733,0.0261443,-1.11519e-05,2.06282e-09,-1.41578e-13,27377.2,-30.8443], Tmin=(1405.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(251.459,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C=C[C]=CC=O(4618)',
    structure = SMILES('[CH2]C=C[C]=CC=O'),
    E0 = (251.459,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.48426,'amu*angstrom^2'), symmetry=1, barrier=(34.1261,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48843,'amu*angstrom^2'), symmetry=1, barrier=(34.222,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4825,'amu*angstrom^2'), symmetry=1, barrier=(34.0857,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32441,0.0561728,-4.32064e-05,1.72705e-08,-2.8472e-12,30341.8,23.6367], Tmin=(100,'K'), Tmax=(1405.2,'K')), NASAPolynomial(coeffs=[11.8733,0.0261443,-1.11519e-05,2.06282e-09,-1.41578e-13,27377.2,-30.8443], Tmin=(1405.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(251.459,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2][C]=CC=CC=O(4619)',
    structure = SMILES('C=C=C[CH]C=C[O]'),
    E0 = (263.723,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,540,610,2055,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.61406,'amu*angstrom^2'), symmetry=1, barrier=(37.1105,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.61387,'amu*angstrom^2'), symmetry=1, barrier=(37.106,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3233,0.039683,3.0721e-05,-7.93297e-08,3.64771e-11,31832.7,24.2304], Tmin=(100,'K'), Tmax=(931.617,'K')), NASAPolynomial(coeffs=[19.4903,0.0106915,-1.51183e-06,2.05831e-10,-2.00481e-14,26320.9,-73.543], Tmin=(931.617,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(263.723,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=CCJC=C)"""),
)

species(
    label = 'HCO(14)(15)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.23755,-0.00332075,1.4003e-05,-1.3424e-08,4.37416e-12,3872.41,3.30835], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.92002,0.00252279,-6.71004e-07,1.05616e-10,-7.43798e-15,3653.43,3.58077], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(32.4782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HCO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH]=CC=C[CH2](922)',
    structure = SMILES('[CH]=CC=C[CH2]'),
    E0 = (423.048,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(1.52606,'amu*angstrom^2'), symmetry=1, barrier=(35.0872,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.53677,'amu*angstrom^2'), symmetry=1, barrier=(35.3333,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22719,0.0272982,2.2821e-05,-5.20811e-08,2.2993e-11,50955.4,18.1254], Tmin=(100,'K'), Tmax=(937.462,'K')), NASAPolynomial(coeffs=[12.1491,0.0145308,-4.06041e-06,6.79572e-10,-4.92004e-14,47795.9,-36.031], Tmin=(937.462,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(423.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C=CC=[C]C=O(4620)',
    structure = SMILES('[CH2]C=CC=C=C[O]'),
    E0 = (249.2,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,540,610,2055,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.60314,'amu*angstrom^2'), symmetry=1, barrier=(36.8594,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.60245,'amu*angstrom^2'), symmetry=1, barrier=(36.8435,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05776,0.0479039,8.36808e-06,-5.70523e-08,2.88729e-11,30093.3,24.049], Tmin=(100,'K'), Tmax=(928.415,'K')), NASAPolynomial(coeffs=[19.6183,0.0110499,-1.74383e-06,2.26183e-10,-1.93873e-14,24788.9,-74.1231], Tmin=(928.415,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(249.2,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C=CC=C[C]=O(4621)',
    structure = SMILES('[CH2]C=CC=C[C]=O'),
    E0 = (213.085,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(1.17179,'amu*angstrom^2'), symmetry=1, barrier=(26.9417,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.174,'amu*angstrom^2'), symmetry=1, barrier=(26.9927,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1662,'amu*angstrom^2'), symmetry=1, barrier=(26.8131,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00471,0.061723,-5.34234e-05,2.34529e-08,-4.14141e-12,25739.3,24.7144], Tmin=(100,'K'), Tmax=(1346.57,'K')), NASAPolynomial(coeffs=[14.202,0.0225201,-9.75334e-06,1.83241e-09,-1.27382e-13,22185.2,-42.8817], Tmin=(1346.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(213.085,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJ=O) + radical(C=CC=CCJ)"""),
)

species(
    label = 'O=CC=CC1[CH]C1(4622)',
    structure = SMILES('O=CC=CC1[CH]C1'),
    E0 = (216.747,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72157,0.0395919,4.30585e-06,-2.70508e-08,1.08034e-11,26159.5,25.2432], Tmin=(100,'K'), Tmax=(1111.19,'K')), NASAPolynomial(coeffs=[11.3657,0.0286289,-1.29597e-05,2.54509e-09,-1.83264e-13,22549.7,-28.8998], Tmin=(1111.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.747,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclopropane) + radical(cyclopropane)"""),
)

species(
    label = 'O=CC1[CH]C=CC1(4623)',
    structure = SMILES('O=CC1C=C[CH]C1'),
    E0 = (35.8276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.11392,0.0269593,4.38082e-05,-6.98497e-08,2.65915e-11,4389.75,19.1074], Tmin=(100,'K'), Tmax=(1004.52,'K')), NASAPolynomial(coeffs=[10.8226,0.028065,-1.12769e-05,2.17093e-09,-1.58331e-13,834.743,-31.9325], Tmin=(1004.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(35.8276,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclopentene) + radical(cyclopentene-allyl)"""),
)

species(
    label = '[CH2]CC=CC=C=O(4624)',
    structure = SMILES('[CH2]CC=CC=C=O'),
    E0 = (174.989,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2120,512.5,787.5,3000,3100,440,815,1455,1000,180,2791.94],'cm^-1')),
        HinderedRotor(inertia=(0.578556,'amu*angstrom^2'), symmetry=1, barrier=(13.3021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06866,'amu*angstrom^2'), symmetry=1, barrier=(24.5706,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.44898,'amu*angstrom^2'), symmetry=1, barrier=(10.3229,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65847,0.0546242,-4.19395e-05,1.81731e-08,-3.47164e-12,21127.9,25.222], Tmin=(100,'K'), Tmax=(1169.5,'K')), NASAPolynomial(coeffs=[7.77678,0.0336979,-1.50994e-05,2.87312e-09,-2.01002e-13,19696.9,-5.25333], Tmin=(1169.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(174.989,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CCO)H) + group(Cds-(Cdd-O2d)CsH) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=CC=C=CO(4625)',
    structure = SMILES('[CH2]C=CC=C=CO'),
    E0 = (107.737,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,540,610,2055,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.38218,'amu*angstrom^2'), symmetry=1, barrier=(31.7791,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.3813,'amu*angstrom^2'), symmetry=1, barrier=(31.7589,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.38334,'amu*angstrom^2'), symmetry=1, barrier=(31.8056,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.643026,0.0548974,4.97149e-06,-6.3447e-08,3.3789e-11,13096.4,24.1284], Tmin=(100,'K'), Tmax=(912.632,'K')), NASAPolynomial(coeffs=[22.9539,0.00722479,9.58476e-07,-3.47076e-10,2.17362e-14,6937.04,-92.9028], Tmin=(912.632,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(107.737,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=CC=CCJ)"""),
)

species(
    label = 'CH2(17)(18)',
    structure = SMILES('[CH2]'),
    E0 = (381.08,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([971.045,2816.03,3444.23],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.71758,0.00127391,2.17347e-06,-3.48858e-09,1.65209e-12,45872.4,1.75298], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.14632,0.00303671,-9.96474e-07,1.50484e-10,-8.57336e-15,46041.3,4.72342], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(381.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH]=CC=CC=O(4626)',
    structure = SMILES('[CH]=CC=CC=O'),
    E0 = (217.529,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.987604,'amu*angstrom^2'), symmetry=1, barrier=(22.707,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.985774,'amu*angstrom^2'), symmetry=1, barrier=(22.6649,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54247,0.0481101,-3.84326e-05,1.48096e-08,-2.25937e-12,26256,20.3082], Tmin=(100,'K'), Tmax=(1555.32,'K')), NASAPolynomial(coeffs=[14.223,0.0154973,-6.9792e-06,1.3273e-09,-9.21979e-14,22311.6,-46.4688], Tmin=(1555.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(217.529,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH]C=CC=CC=O(4627)',
    structure = SMILES('[CH]C=CC=CC=O'),
    E0 = (305.092,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2782.5,750,1395,475,1775,1000,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.0413,'amu*angstrom^2'), symmetry=1, barrier=(46.9334,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.03446,'amu*angstrom^2'), symmetry=1, barrier=(46.7763,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.03922,'amu*angstrom^2'), symmetry=1, barrier=(46.8857,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.959932,0.0600947,-4.14419e-05,1.38059e-08,-1.84722e-12,36808.9,24.9174], Tmin=(100,'K'), Tmax=(1723.6,'K')), NASAPolynomial(coeffs=[16.0756,0.0250152,-1.0913e-05,1.99763e-09,-1.34479e-13,31598.3,-56.236], Tmin=(1723.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(305.092,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=CC1C=CC1[O](4628)',
    structure = SMILES('C=CC1C=CC1[O]'),
    E0 = (273.723,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52118,0.039214,2.13865e-05,-5.67847e-08,2.45281e-11,33024.2,22.2187], Tmin=(100,'K'), Tmax=(980.966,'K')), NASAPolynomial(coeffs=[14.9453,0.0214197,-7.8957e-06,1.5074e-09,-1.11852e-13,28613,-51.3473], Tmin=(980.966,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(273.723,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=CC=C=CC=O(4629)',
    structure = SMILES('C=CC=C=CC=O'),
    E0 = (111.012,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,610,2055,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.07554,'amu*angstrom^2'), symmetry=1, barrier=(24.7288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07684,'amu*angstrom^2'), symmetry=1, barrier=(24.7587,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14219,0.0583098,-4.6417e-05,1.80867e-08,-2.82455e-12,13458,22.1956], Tmin=(100,'K'), Tmax=(1503.78,'K')), NASAPolynomial(coeffs=[14.912,0.0216827,-9.88211e-06,1.88985e-09,-1.31861e-13,9316.65,-49.8536], Tmin=(1503.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(111.012,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=CC[C]=CC=O(4630)',
    structure = SMILES('C=CC[C]=CC=O'),
    E0 = (203.104,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,260.934,260.951],'cm^-1')),
        HinderedRotor(inertia=(0.261896,'amu*angstrom^2'), symmetry=1, barrier=(12.6533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.261858,'amu*angstrom^2'), symmetry=1, barrier=(12.6532,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.26192,'amu*angstrom^2'), symmetry=1, barrier=(12.6534,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8411,0.0524533,-3.41345e-05,1.05484e-08,-1.34473e-12,24500,24.1019], Tmin=(100,'K'), Tmax=(1683.28,'K')), NASAPolynomial(coeffs=[11.0729,0.0305155,-1.45851e-05,2.80575e-09,-1.9478e-13,21392.1,-25.2433], Tmin=(1683.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(203.104,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]CC=CC=O(4631)',
    structure = SMILES('C=[C]CC=CC=O'),
    E0 = (203.104,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000,260.907,260.938],'cm^-1')),
        HinderedRotor(inertia=(0.261926,'amu*angstrom^2'), symmetry=1, barrier=(12.6548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.261878,'amu*angstrom^2'), symmetry=1, barrier=(12.6521,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.261791,'amu*angstrom^2'), symmetry=1, barrier=(12.6529,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8411,0.0524533,-3.41345e-05,1.05484e-08,-1.34473e-12,24500,24.1019], Tmin=(100,'K'), Tmax=(1683.28,'K')), NASAPolynomial(coeffs=[11.0729,0.0305155,-1.45851e-05,2.80575e-09,-1.9478e-13,21392.1,-25.2433], Tmin=(1683.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(203.104,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_S)"""),
)

species(
    label = 'C=CCC=[C]C=O(4632)',
    structure = SMILES('C=CCC=C=C[O]'),
    E0 = (161.998,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,610,2055,2950,3100,1380,975,1025,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.07773,'amu*angstrom^2'), symmetry=1, barrier=(24.7791,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07891,'amu*angstrom^2'), symmetry=1, barrier=(24.8063,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.990035,0.0522398,-7.46635e-06,-3.39512e-08,1.84062e-11,19604.9,25.7561], Tmin=(100,'K'), Tmax=(962.781,'K')), NASAPolynomial(coeffs=[17.4007,0.0175999,-5.75338e-06,1.0465e-09,-7.71329e-14,14890.4,-60.8666], Tmin=(962.781,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(161.998,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ)"""),
)

species(
    label = '[CH]=CCC=CC=O(4633)',
    structure = SMILES('[CH]=CCC=CC=O'),
    E0 = (212.358,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.635542,'amu*angstrom^2'), symmetry=1, barrier=(14.6124,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.635784,'amu*angstrom^2'), symmetry=1, barrier=(14.6179,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.635937,'amu*angstrom^2'), symmetry=1, barrier=(14.6214,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48094,0.054882,-3.72214e-05,1.20021e-08,-1.56431e-12,25631,25.3434], Tmin=(100,'K'), Tmax=(1723.81,'K')), NASAPolynomial(coeffs=[13.886,0.026097,-1.2174e-05,2.31531e-09,-1.59483e-13,21354.2,-41.2591], Tmin=(1723.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(212.358,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_P)"""),
)

species(
    label = 'C=CCC=C[C]=O(4634)',
    structure = SMILES('C=CC[CH]C=C=O'),
    E0 = (125.841,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,2120,512.5,787.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,297.189,297.198,1594.87],'cm^-1')),
        HinderedRotor(inertia=(0.135226,'amu*angstrom^2'), symmetry=1, barrier=(8.47547,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.341571,'amu*angstrom^2'), symmetry=1, barrier=(21.409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.341574,'amu*angstrom^2'), symmetry=1, barrier=(21.4091,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46375,0.0557229,-4.04047e-05,1.52681e-08,-2.42408e-12,15226.1,25.6131], Tmin=(100,'K'), Tmax=(1420.9,'K')), NASAPolynomial(coeffs=[10.6537,0.0298519,-1.30933e-05,2.45389e-09,-1.69478e-13,12614.5,-21.9515], Tmin=(1420.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(125.841,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + radical(CCJC(C)=C=O)"""),
)

species(
    label = 'C=CC1[CH]C1C=O(4635)',
    structure = SMILES('C=CC1[CH]C1C=O'),
    E0 = (201.884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.82358,0.0315837,3.88578e-05,-7.35423e-08,3.04209e-11,24374.1,25.3845], Tmin=(100,'K'), Tmax=(966.563,'K')), NASAPolynomial(coeffs=[14.1093,0.0210131,-7.23589e-06,1.35665e-09,-1.0109e-13,20117.9,-43.2007], Tmin=(966.563,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(201.884,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(CCJCC=O)"""),
)

species(
    label = 'C=CC1C=C[CH]O1(4636)',
    structure = SMILES('C=CC1[CH]C=CO1'),
    E0 = (8.86229,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00146,0.00619038,0.000155704,-2.26582e-07,9.34673e-11,1173.04,19.1461], Tmin=(100,'K'), Tmax=(915.32,'K')), NASAPolynomial(coeffs=[25.9945,-0.000838938,6.91601e-06,-1.4547e-09,8.87006e-14,-7317.02,-116.868], Tmin=(915.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(8.86229,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(2,3-Dihydrofuran) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = 'C2H3(28)(29)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36378,0.000265766,2.79621e-05,-3.72987e-08,1.5159e-11,34475,7.9151], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.15027,0.00754021,-2.62998e-06,4.15974e-10,-2.45408e-14,33856.6,1.72812], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(286.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""C2H3""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH]C=CC=O(4637)',
    structure = SMILES('[CH]C=CC=O'),
    E0 = (253.32,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.13955,'amu*angstrom^2'), symmetry=1, barrier=(49.1924,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.13866,'amu*angstrom^2'), symmetry=1, barrier=(49.172,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.45042,0.0332273,-1.49834e-05,1.57226e-09,2.41248e-13,30523.8,16.9403], Tmin=(100,'K'), Tmax=(1750.2,'K')), NASAPolynomial(coeffs=[11.9107,0.0188928,-8.94295e-06,1.65009e-09,-1.09645e-13,26096.4,-37.1832], Tmin=(1750.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(253.32,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=C[CH]C1C=CO1(4638)',
    structure = SMILES('C=C[CH]C1C=CO1'),
    E0 = (112.868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44715,0.0181073,0.000131028,-2.03634e-07,8.46999e-11,13702,22.7015], Tmin=(100,'K'), Tmax=(931.042,'K')), NASAPolynomial(coeffs=[28.7014,-0.00139157,5.21199e-06,-9.596e-10,4.77989e-14,4397.17,-129.553], Tmin=(931.042,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(112.868,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = '[CH2]C1C=CC=CO1(4639)',
    structure = SMILES('[CH2]C1C=CC=CO1'),
    E0 = (116.918,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05153,0.0380851,5.85802e-05,-1.20354e-07,5.35544e-11,14192.9,18.7753], Tmin=(100,'K'), Tmax=(930.4,'K')), NASAPolynomial(coeffs=[24.7787,0.00469913,1.77181e-06,-3.75882e-10,1.52093e-14,6807.62,-109.944], Tmin=(930.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(116.918,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(1,3-Cyclohexadiene) + radical(CJC(C)OC)"""),
)

species(
    label = 'C=CC=CC=[C]O(4640)',
    structure = SMILES('C=CC=CC=[C]O'),
    E0 = (176.645,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.02392,'amu*angstrom^2'), symmetry=1, barrier=(23.5419,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02341,'amu*angstrom^2'), symmetry=1, barrier=(23.5302,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02385,'amu*angstrom^2'), symmetry=1, barrier=(23.5404,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.558767,0.0589195,-1.17097e-05,-4.45582e-08,2.67326e-11,21385,25.9711], Tmin=(100,'K'), Tmax=(913.538,'K')), NASAPolynomial(coeffs=[22.6267,0.00681664,7.35319e-07,-2.89898e-10,1.81949e-14,15495.1,-88.667], Tmin=(913.538,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(176.645,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO)"""),
)

species(
    label = 'C=CC=[C]C=CO(4641)',
    structure = SMILES('C=CC=[C]C=CO'),
    E0 = (135.896,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.28045,'amu*angstrom^2'), symmetry=1, barrier=(29.44,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.27915,'amu*angstrom^2'), symmetry=1, barrier=(29.4101,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.27854,'amu*angstrom^2'), symmetry=1, barrier=(29.3961,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.272626,0.0636521,-1.35483e-05,-5.12935e-08,3.17303e-11,16496,23.4682], Tmin=(100,'K'), Tmax=(896.394,'K')), NASAPolynomial(coeffs=[25.3262,0.00261242,3.65767e-06,-9.21272e-10,6.43408e-14,9965.2,-106.035], Tmin=(896.394,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(135.896,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C[C]=CC=CO(4642)',
    structure = SMILES('C=C[C]=CC=CO'),
    E0 = (135.896,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.28045,'amu*angstrom^2'), symmetry=1, barrier=(29.44,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.27915,'amu*angstrom^2'), symmetry=1, barrier=(29.4101,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.27854,'amu*angstrom^2'), symmetry=1, barrier=(29.3961,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.272626,0.0636521,-1.35483e-05,-5.12935e-08,3.17303e-11,16496,23.4682], Tmin=(100,'K'), Tmax=(896.394,'K')), NASAPolynomial(coeffs=[25.3262,0.00261242,3.65767e-06,-9.21272e-10,6.43408e-14,9965.2,-106.035], Tmin=(896.394,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(135.896,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = 'C=[C]C=CC=CO(4643)',
    structure = SMILES('C=C=C[CH]C=CO'),
    E0 = (122.261,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,540,610,2055,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.3242,'amu*angstrom^2'), symmetry=1, barrier=(30.446,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32808,'amu*angstrom^2'), symmetry=1, barrier=(30.5353,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.3272,'amu*angstrom^2'), symmetry=1, barrier=(30.515,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.906533,0.0467041,2.72076e-05,-8.55384e-08,4.12957e-11,14835.8,24.3168], Tmin=(100,'K'), Tmax=(917.215,'K')), NASAPolynomial(coeffs=[22.8223,0.00687185,1.18754e-06,-3.66776e-10,2.10239e-14,8470.73,-92.3023], Tmin=(917.215,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(122.261,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJC=C)"""),
)

species(
    label = '[CH]=CC=CC=CO(4644)',
    structure = SMILES('[CH]=CC=CC=CO'),
    E0 = (183.997,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.14741,'amu*angstrom^2'), symmetry=1, barrier=(26.3811,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14602,'amu*angstrom^2'), symmetry=1, barrier=(26.3493,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14807,'amu*angstrom^2'), symmetry=1, barrier=(26.3964,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.272027,0.0607078,-1.48686e-06,-6.61122e-08,3.71915e-11,22283.9,24.2756], Tmin=(100,'K'), Tmax=(906.758,'K')), NASAPolynomial(coeffs=[27.052,-0.000217445,4.65942e-06,-1.05053e-09,6.96253e-14,15075.4,-115.27], Tmin=(906.758,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(183.997,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'CHCHO(45)(45)',
    structure = SMILES('[CH]=C[O]'),
    E0 = (150.672,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3120,650,792.5,1650],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.11,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.3215,0.0099565,8.23489e-06,-1.92986e-08,8.37464e-12,18150.6,9.47663], Tmin=(100,'K'), Tmax=(969.366,'K')), NASAPolynomial(coeffs=[7.73867,0.0039342,-1.33182e-06,2.69015e-10,-2.15546e-14,16720.8,-14.6539], Tmin=(969.366,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(150.672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""ketene(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CH2CHCHCH(1116)',
    structure = SMILES('[CH]=CC=C'),
    E0 = (346.45,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(1.31937,'amu*angstrom^2'), symmetry=1, barrier=(30.3349,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (53.0825,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.64255,0.0163337,3.86225e-05,-6.71377e-08,2.83603e-11,41729.6,13.282], Tmin=(100,'K'), Tmax=(937.724,'K')), NASAPolynomial(coeffs=[12.9705,0.00669127,-1.0007e-06,1.67602e-10,-1.71436e-14,38279.7,-43.9476], Tmin=(937.724,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(346.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""CH2CHCHCH""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=CC=CC1[CH]O1(4645)',
    structure = SMILES('C=CC=CC1[CH]O1'),
    E0 = (235.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05641,0.0464782,2.363e-05,-8.2232e-08,4.16625e-11,28475.1,23.7946], Tmin=(100,'K'), Tmax=(882.73,'K')), NASAPolynomial(coeffs=[20.8322,0.00788428,2.51806e-06,-8.13822e-10,6.08922e-14,22996,-80.4043], Tmin=(882.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(235.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCsJO)"""),
)

species(
    label = 'C=C=CC=CC[O](4646)',
    structure = SMILES('C=C=CC=CC[O]'),
    E0 = (271.044,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,610,2055,2950,3100,1380,975,1025,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.60286,'amu*angstrom^2'), symmetry=1, barrier=(13.8609,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.603168,'amu*angstrom^2'), symmetry=1, barrier=(13.868,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56707,0.0556824,-4.21885e-05,1.76457e-08,-3.20116e-12,32684.6,24.0088], Tmin=(100,'K'), Tmax=(1237.67,'K')), NASAPolynomial(coeffs=[8.6303,0.0328549,-1.45225e-05,2.7435e-09,-1.91018e-13,30936.2,-11.5733], Tmin=(1237.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(271.044,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCOJ)"""),
)

species(
    label = '[CH]=CC=CC=C(4647)',
    structure = SMILES('[CH]=CC=CC=C'),
    E0 = (392.789,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.1072,'amu*angstrom^2'), symmetry=1, barrier=(25.4567,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10579,'amu*angstrom^2'), symmetry=1, barrier=(25.4242,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31285,0.0449348,2.88003e-06,-4.44384e-08,2.28046e-11,47351.4,20.6785], Tmin=(100,'K'), Tmax=(937.294,'K')), NASAPolynomial(coeffs=[17.2005,0.0128924,-3.06901e-06,4.97417e-10,-3.78225e-14,42802.3,-63.3206], Tmin=(937.294,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(392.789,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C1C=CC1C=O(4648)',
    structure = SMILES('[CH2]C1C=CC1C=O'),
    E0 = (204.323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69057,0.0398055,8.74228e-06,-3.79075e-08,1.70087e-11,24667.3,23.4548], Tmin=(100,'K'), Tmax=(976.874,'K')), NASAPolynomial(coeffs=[11.6085,0.0253105,-9.10149e-06,1.63689e-09,-1.15127e-13,21483.5,-30.5393], Tmin=(976.874,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.323,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclobutene) + radical(Isobutyl)"""),
)

species(
    label = 'C=CC=[C]CC=O(4649)',
    structure = SMILES('C=CC=[C]CC=O'),
    E0 = (196.475,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000,239.581,239.607],'cm^-1')),
        HinderedRotor(inertia=(0.456224,'amu*angstrom^2'), symmetry=1, barrier=(18.574,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.456627,'amu*angstrom^2'), symmetry=1, barrier=(18.5726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.45643,'amu*angstrom^2'), symmetry=1, barrier=(18.5725,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20921,0.0505654,-1.73009e-05,-1.1011e-08,6.50066e-12,23740,24.8016], Tmin=(100,'K'), Tmax=(1125.67,'K')), NASAPolynomial(coeffs=[14.3643,0.0250131,-1.14926e-05,2.27461e-09,-1.64513e-13,19435.6,-46.1855], Tmin=(1125.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(196.475,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C=CC=CC[C]=O(3250)',
    structure = SMILES('C=CC=CC[C]=O'),
    E0 = (118.594,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.783372,'amu*angstrom^2'), symmetry=1, barrier=(18.0113,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.782173,'amu*angstrom^2'), symmetry=1, barrier=(17.9837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.783141,'amu*angstrom^2'), symmetry=1, barrier=(18.006,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19677,0.0487089,-7.28985e-06,-2.55703e-08,1.26837e-11,14375.7,25.1849], Tmin=(100,'K'), Tmax=(1042.88,'K')), NASAPolynomial(coeffs=[15.5436,0.021984,-9.5591e-06,1.90322e-09,-1.40485e-13,9844.15,-52.0117], Tmin=(1042.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(118.594,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O)"""),
)

species(
    label = 'C=C[C]=CCC=O(4650)',
    structure = SMILES('C=C[C]=CCC=O'),
    E0 = (157.629,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.05138,'amu*angstrom^2'), symmetry=1, barrier=(24.1732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05095,'amu*angstrom^2'), symmetry=1, barrier=(24.1635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05009,'amu*angstrom^2'), symmetry=1, barrier=(24.1436,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14085,0.052814,-2.28167e-05,-5.61669e-09,4.8368e-12,19069.7,24.168], Tmin=(100,'K'), Tmax=(1123.45,'K')), NASAPolynomial(coeffs=[13.9442,0.0257035,-1.12871e-05,2.17957e-09,-1.55476e-13,15027,-44.2801], Tmin=(1123.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.629,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = 'C=[C]C=CCC=O(4651)',
    structure = SMILES('C=[C]C=CCC=O'),
    E0 = (157.629,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.05138,'amu*angstrom^2'), symmetry=1, barrier=(24.1732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05095,'amu*angstrom^2'), symmetry=1, barrier=(24.1635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05009,'amu*angstrom^2'), symmetry=1, barrier=(24.1436,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14085,0.052814,-2.28167e-05,-5.61669e-09,4.8368e-12,19069.7,24.168], Tmin=(100,'K'), Tmax=(1123.45,'K')), NASAPolynomial(coeffs=[13.9442,0.0257035,-1.12871e-05,2.17957e-09,-1.55476e-13,15027,-44.2801], Tmin=(1123.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.629,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = 'C1=CC=COC=C1(4672)',
    structure = SMILES('C1=CC=COC=C1'),
    E0 = (62.1925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79695,-0.0230472,0.000304467,-4.41182e-07,1.86994e-10,7626.85,18.0571], Tmin=(100,'K'), Tmax=(895.998,'K')), NASAPolynomial(coeffs=[50.3118,-0.0513285,3.65719e-05,-7.30064e-09,4.87955e-13,-8625.59,-252.851], Tmin=(895.998,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(62.1925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cycloheptane)"""),
)

species(
    label = 'C1=CC=COCC=1(4673)',
    structure = SMILES('C1=CC=COCC=1'),
    E0 = (156.102,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62506,0.0279778,6.84832e-05,-1.1979e-07,5.08622e-11,18882.6,13.1444], Tmin=(100,'K'), Tmax=(935.947,'K')), NASAPolynomial(coeffs=[20.4504,0.00935258,-7.5856e-07,1.1262e-10,-1.79318e-14,12650.5,-90.8983], Tmin=(935.947,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(156.102,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene)"""),
)

species(
    label = '[CH]1CC=CC=CO1(4674)',
    structure = SMILES('[CH]1CC=CC=CO1'),
    E0 = (129.248,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13921,-0.00672645,0.000270561,-4.1076e-07,1.76905e-10,15713.7,20.5556], Tmin=(100,'K'), Tmax=(896.463,'K')), NASAPolynomial(coeffs=[52.1866,-0.0504501,3.57634e-05,-7.1322e-09,4.76592e-13,-834.188,-261.387], Tmin=(896.463,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(129.248,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(CCsJOC(O))"""),
)

species(
    label = '[C]1=CC=COCC1(4675)',
    structure = SMILES('[C]1=CC=COCC1'),
    E0 = (173.163,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48172,-0.0116043,0.000272834,-4.0416e-07,1.72257e-10,20980.5,20.3256], Tmin=(100,'K'), Tmax=(897.455,'K')), NASAPolynomial(coeffs=[48.4407,-0.044379,3.25727e-05,-6.51552e-09,4.34423e-13,5442.96,-240.75], Tmin=(897.455,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(173.163,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(Cds_S)"""),
)

species(
    label = '[C]1=CCCOC=C1(4676)',
    structure = SMILES('[C]1=CCCOC=C1'),
    E0 = (134.316,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42911,-0.00959378,0.000268431,-4.00688e-07,1.71661e-10,16309.6,19.6387], Tmin=(100,'K'), Tmax=(894.435,'K')), NASAPolynomial(coeffs=[48.2145,-0.0440011,3.29514e-05,-6.6503e-09,4.46684e-13,947.3,-239.946], Tmin=(894.435,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(134.316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(C=CJC=C)"""),
)

species(
    label = '[C]1=CC=CCCO1(4677)',
    structure = SMILES('[C]1=CC=CCCO1'),
    E0 = (175.065,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72646,-0.014462,0.000270756,-3.94589e-07,1.66932e-10,21198.1,22.1017], Tmin=(100,'K'), Tmax=(898.49,'K')), NASAPolynomial(coeffs=[45.4653,-0.0397108,2.99788e-05,-6.00694e-09,3.99538e-13,6497.65,-222.299], Tmin=(898.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(175.065,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(C=CJO)"""),
)

species(
    label = '[C]1=COCCC=C1(4678)',
    structure = SMILES('[C]1=COCCC=C1'),
    E0 = (134.316,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42911,-0.00959378,0.000268431,-4.00688e-07,1.71661e-10,16309.6,19.6387], Tmin=(100,'K'), Tmax=(894.435,'K')), NASAPolynomial(coeffs=[48.2145,-0.0440011,3.29514e-05,-6.6503e-09,4.46684e-13,947.3,-239.946], Tmin=(894.435,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(134.316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]1[CH]OC=CC=C1(4679)',
    structure = SMILES('[CH]1C=C[CH]OC=C1'),
    E0 = (172.151,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.41077,-0.0367051,0.000333145,-4.60151e-07,1.90562e-10,20829.6,19.3967], Tmin=(100,'K'), Tmax=(900.294,'K')), NASAPolynomial(coeffs=[46.0759,-0.043464,3.24329e-05,-6.45916e-09,4.27045e-13,5378.99,-228.82], Tmin=(900.294,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(172.151,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(C=CCJC=C) + radical(C=CCJ(O)C)"""),
)

species(
    label = '[C]1[CH]COC=CC=1(4680)',
    structure = SMILES('[C]1[CH]COC=CC=1'),
    E0 = (290.079,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26406,-0.00765961,0.000265138,-4.00418e-07,1.71925e-10,35051,18.3061], Tmin=(100,'K'), Tmax=(897.681,'K')), NASAPolynomial(coeffs=[50.5701,-0.0483978,3.41634e-05,-6.79521e-09,4.52785e-13,18987.9,-254.409], Tmin=(897.681,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(290.079,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(Cds_S) + radical(C=CCJCO)"""),
)

species(
    label = '[C]1=C[CH]COC=C1(4681)',
    structure = SMILES('[C]1=C[CH]COC=C1'),
    E0 = (251.233,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21117,-0.00564528,0.00026072,-3.96923e-07,1.71318e-10,30380.1,17.6201], Tmin=(100,'K'), Tmax=(894.713,'K')), NASAPolynomial(coeffs=[50.3443,-0.0480205,3.45425e-05,-6.93008e-09,4.65054e-13,14492.2,-253.607], Tmin=(894.713,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(251.233,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(C=CJC=C) + radical(C=CCJCO)"""),
)

species(
    label = '[C]1=COC[CH]C=C1(4682)',
    structure = SMILES('[C]1=COC[CH]C=C1'),
    E0 = (251.233,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21117,-0.00564528,0.00026072,-3.96923e-07,1.71318e-10,30380.1,17.6201], Tmin=(100,'K'), Tmax=(894.713,'K')), NASAPolynomial(coeffs=[50.3443,-0.0480205,3.45425e-05,-6.93008e-09,4.65054e-13,14492.2,-253.607], Tmin=(894.713,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(251.233,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(C=CJC=C) + radical(C=CCJCO)"""),
)

species(
    label = '[C]1=CC=C[CH]CO1(4683)',
    structure = SMILES('[C]1=CC=C[CH]CO1'),
    E0 = (291.982,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50936,-0.0105244,0.000263087,-3.90886e-07,1.66618e-10,35268.5,20.0801], Tmin=(100,'K'), Tmax=(898.694,'K')), NASAPolynomial(coeffs=[47.5934,-0.0437273,3.15681e-05,-6.28629e-09,4.17871e-13,20043.2,-235.951], Tmin=(898.694,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(291.982,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(C=CCJCO) + radical(C=CJO)"""),
)

species(
    label = '[CH]1C2C=COCC12(4684)',
    structure = SMILES('[CH]1C2C=COCC12'),
    E0 = (186.666,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90245,0.0215449,8.48438e-05,-1.30547e-07,5.26387e-11,22548.6,16.5638], Tmin=(100,'K'), Tmax=(949.952,'K')), NASAPolynomial(coeffs=[18.0737,0.015657,-4.08272e-06,7.92933e-10,-6.72985e-14,16669.5,-75.3958], Tmin=(949.952,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(186.666,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + polycyclic(s2_3_6_ene_1) + radical(cyclopropane)"""),
)

species(
    label = '[CH]1OCC2C=CC12(4685)',
    structure = SMILES('[CH]1OCC2C=CC12'),
    E0 = (228.894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94762,0.0281676,4.88802e-05,-8.12273e-08,3.21286e-11,27618.6,15.8767], Tmin=(100,'K'), Tmax=(981.015,'K')), NASAPolynomial(coeffs=[13.2509,0.0240568,-9.01868e-06,1.73674e-09,-1.29309e-13,23380.9,-48.7336], Tmin=(981.015,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(228.894,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s2_4_5_ane) - ring(Tetrahydrofuran) - ring(Cyclobutane) + ring(Tetrahydrofuran) + ring(Cyclobutene) + radical(CCsJOCs)"""),
)

species(
    label = '[CH]1C=CC2COC12(4686)',
    structure = SMILES('[CH]1C=CC2COC12'),
    E0 = (132.839,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95946,0.0257675,6.17886e-05,-9.44273e-08,3.61495e-11,16067.3,14.3619], Tmin=(100,'K'), Tmax=(995.324,'K')), NASAPolynomial(coeffs=[13.6538,0.0263228,-1.07117e-05,2.13337e-09,-1.60497e-13,11383.9,-53.834], Tmin=(995.324,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(132.839,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_ene_1) + radical(C=CCJCO)"""),
)

species(
    label = '[C]1C=CC=COC1(4687)',
    structure = SMILES('[C]1=C[CH]C=COC1'),
    E0 = (299.053,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.18214,-0.0299402,0.000314061,-4.42631e-07,1.85223e-10,36099.2,20.4493], Tmin=(100,'K'), Tmax=(897.718,'K')), NASAPolynomial(coeffs=[46.8182,-0.045425,3.34873e-05,-6.69459e-09,4.45648e-13,20694.9,-231.239], Tmin=(897.718,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(299.053,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(C=CCJC=C) + radical(Cds_S)"""),
)

species(
    label = 'C1=CC=CCOC=1(4688)',
    structure = SMILES('C1=CC=CCOC=1'),
    E0 = (156.102,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62513,0.0279769,6.84864e-05,-1.19794e-07,5.08642e-11,18882.6,13.1441], Tmin=(100,'K'), Tmax=(935.94,'K')), NASAPolynomial(coeffs=[20.4502,0.00935295,-7.58778e-07,1.12672e-10,-1.79362e-14,12650.6,-90.8971], Tmin=(935.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(156.102,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene)"""),
)

species(
    label = '[C]1=CCOC=CC1(4689)',
    structure = SMILES('[C]1=CCOC=CC1'),
    E0 = (197.328,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87082,-0.0176256,0.000276614,-3.97994e-07,1.6734e-10,23870.5,21.8955], Tmin=(100,'K'), Tmax=(899.908,'K')), NASAPolynomial(coeffs=[44.5396,-0.0381967,2.90618e-05,-5.81058e-09,3.84975e-13,9344.24,-217.498], Tmin=(899.908,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(197.328,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(Cds_S)"""),
)

species(
    label = '[C]1=COCC=CC1(4690)',
    structure = SMILES('[C]1=COCC=CC1'),
    E0 = (197.328,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87082,-0.0176256,0.000276614,-3.97994e-07,1.6734e-10,23870.5,21.8955], Tmin=(100,'K'), Tmax=(899.908,'K')), NASAPolynomial(coeffs=[44.5396,-0.0381967,2.90618e-05,-5.81058e-09,3.84975e-13,9344.24,-217.498], Tmin=(899.908,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(197.328,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(Cds_S)"""),
)

species(
    label = '[C]1=CCC=COC1(4691)',
    structure = SMILES('[C]1=CCC=COC1'),
    E0 = (197.328,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87082,-0.0176256,0.000276614,-3.97994e-07,1.6734e-10,23870.5,21.8955], Tmin=(100,'K'), Tmax=(899.908,'K')), NASAPolynomial(coeffs=[44.5396,-0.0381967,2.90618e-05,-5.81058e-09,3.84975e-13,9344.24,-217.498], Tmin=(899.908,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(197.328,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(Cds_S)"""),
)

species(
    label = '[C]1=CCC=CCO1(4692)',
    structure = SMILES('[C]1=CCC=CCO1'),
    E0 = (199.231,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.11536,-0.0204811,0.00027453,-3.88418e-07,1.62015e-10,24088.1,23.6723], Tmin=(100,'K'), Tmax=(901.11,'K')), NASAPolynomial(coeffs=[41.5663,-0.0335322,2.647e-05,-5.3025e-09,3.50131e-13,10398.1,-199.059], Tmin=(901.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(199.231,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(C=CJO)"""),
)

species(
    label = '[CH]1C=CCC=CO1(4693)',
    structure = SMILES('[CH]1C=CCC=CO1'),
    E0 = (70.4258,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.09942,-0.0243911,0.000295705,-4.15535e-07,1.72693e-10,8600.93,21.5361], Tmin=(100,'K'), Tmax=(902.662,'K')), NASAPolynomial(coeffs=[43.8008,-0.0362418,2.8011e-05,-5.57598e-09,3.66439e-13,-5973.19,-214.406], Tmin=(902.662,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(70.4258,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(C=CCJ(O)C)"""),
)

species(
    label = '[CH]1C2C=CCOC12(4694)',
    structure = SMILES('[CH]1C2C=CCOC12'),
    E0 = (202.527,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.36635,0.0184317,6.91465e-05,-9.63519e-08,3.59189e-11,24432.7,18.1442], Tmin=(100,'K'), Tmax=(992.501,'K')), NASAPolynomial(coeffs=[11.1146,0.0272596,-1.0823e-05,2.11778e-09,-1.57502e-13,20524.8,-34.9338], Tmin=(992.501,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(202.527,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_3_6_ene_1) + radical(CCJCO)"""),
)

species(
    label = '[C]1=CC=CCOC1(4695)',
    structure = SMILES('[C]1=CC=CCOC1'),
    E0 = (224.193,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.10071,-0.0113427,0.000231419,-3.30915e-07,1.38246e-10,27082.5,21.6444], Tmin=(100,'K'), Tmax=(901.765,'K')), NASAPolynomial(coeffs=[36.2195,-0.0236767,2.07085e-05,-4.19495e-09,2.76687e-13,15277.1,-170.769], Tmin=(901.765,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(224.193,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(Cycloheptane) + radical(Cds_S)"""),
)

species(
    label = '[C]1=CCOCC=C1(4696)',
    structure = SMILES('[C]1=CCOCC=C1'),
    E0 = (185.347,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.04864,-0.00933706,0.000227025,-3.27437e-07,1.37637e-10,22411.6,20.9555], Tmin=(100,'K'), Tmax=(897.923,'K')), NASAPolynomial(coeffs=[35.9852,-0.023285,2.10792e-05,-4.32784e-09,2.88792e-13,10784.9,-169.92], Tmin=(897.923,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(185.347,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(Cycloheptane) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]1COC2C=CC12(4697)',
    structure = SMILES('[CH]1COC2C=CC12'),
    E0 = (243.332,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.43464,0.0138464,8.66627e-05,-1.17231e-07,4.3815e-11,29340.8,17.9895], Tmin=(100,'K'), Tmax=(984.039,'K')), NASAPolynomial(coeffs=[12.5407,0.0249286,-9.74224e-06,1.94913e-09,-1.48721e-13,24826.4,-43.4361], Tmin=(984.039,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.332,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s2_4_5_ane) - ring(Tetrahydrofuran) - ring(Cyclobutane) + ring(Tetrahydrofuran) + ring(Cyclobutene) + radical(CCJCO)"""),
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
    E0 = (223.981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (433.626,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (320.299,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (595.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (595.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (560.911,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (560.911,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (191.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (268.34,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (331.249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (193.883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (308.465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (287.903,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (295.385,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (254.585,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (565.334,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (257.372,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (276.784,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (276.874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (238.798,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (176.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (428.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (248.483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (334.603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (260.195,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (324.222,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (237.008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (164.185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (161.488,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (548.444,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (463.709,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (467.447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (475.516,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (455.526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (465.189,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (424.877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (278.399,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (172.553,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (131.221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (291.544,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (362.039,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (251.604,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (598.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (516.884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (273.723,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (339.599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (360.004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (360.004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (309.275,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (357.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (240.273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (256.433,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (189.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (539.681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (248.483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (184.259,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (324.472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (339.424,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (239.501,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (168.936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (254.156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (335.517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (497.123,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (278.399,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS65',
    E0 = (434.804,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS66',
    E0 = (635.794,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS67',
    E0 = (213.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS68',
    E0 = (331.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS69',
    E0 = (201.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS70',
    E0 = (304.906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS71',
    E0 = (231.585,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS72',
    E0 = (238.77,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS73',
    E0 = (280.754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS74',
    E0 = (383.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS75',
    E0 = (203.682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS76',
    E0 = (330.842,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS77',
    E0 = (304.494,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS78',
    E0 = (219.374,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS79',
    E0 = (237.921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS80',
    E0 = (383.943,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS81',
    E0 = (501.872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS82',
    E0 = (463.483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS83',
    E0 = (467.222,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS84',
    E0 = (503.774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS85',
    E0 = (234.242,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS86',
    E0 = (279.286,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS87',
    E0 = (231.642,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS88',
    E0 = (366.038,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS89',
    E0 = (510.846,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS90',
    E0 = (384.689,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS91',
    E0 = (354.228,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS92',
    E0 = (354.228,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS93',
    E0 = (346.697,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS94',
    E0 = (348.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS95',
    E0 = (185.904,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS96',
    E0 = (297.559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS97',
    E0 = (363.548,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS98',
    E0 = (378.425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS99',
    E0 = (279.286,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction53',
    reactants = ['H(3)(3)', 'O=C1C=CC=CC1(4653)'],
    products = ['C6H7O(505)(504)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(3.97e+07,'cm^3/(mol*s)'), n=1.88, Ea=(32.2168,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2826 used for CO-CdCs_O;HJ
Exact match found for rate rule [CO-CdCs_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction60',
    reactants = ['H(3)(3)', '[O]C1[CH]C=CC=C1(4660)'],
    products = ['C6H7O(505)(504)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5.42928e+07,'m^3/(mol*s)'), n=0.107721, Ea=(5.76381,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 36 used for C_rad/H/CdCs;H_rad
Exact match found for rate rule [C_rad/H/CdCs;H_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction61',
    reactants = ['H(3)(3)', '[O][C]1C=CC=CC1(4661)'],
    products = ['C6H7O(505)(504)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.34601e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction62',
    reactants = ['H(3)(3)', '[O]C1C=CC=[C]C1(4662)'],
    products = ['C6H7O(505)(504)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction63',
    reactants = ['H(3)(3)', '[O]C1[C]=CC=CC1(4663)'],
    products = ['C6H7O(505)(504)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction64',
    reactants = ['H(3)(3)', '[O]C1C=C[C]=CC1(4664)'],
    products = ['C6H7O(505)(504)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [Cd_rad/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction65',
    reactants = ['H(3)(3)', '[O]C1C=[C]C=CC1(4665)'],
    products = ['C6H7O(505)(504)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [H_rad;Cd_rad/Cd]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction70',
    reactants = ['C6H7O(505)(504)'],
    products = ['[O]C1=CC=CCC1(4670)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.125e+09,'s^-1'), n=0.991, Ea=(45.9529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [1_3_pentadiene;CHOC_1;CdHC_2]
Euclidian distance = 0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction71',
    reactants = ['C6H7O(505)(504)'],
    products = ['[O]C1CC2C=CC12(4671)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1,3-butadiene_backbone;C=C_1;C=C_2] for rate rule [1,3-butadiene_backbone;CdH(C)_1;CdH(C)_2]
Euclidian distance = 1.41421356237
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH]=CC=CCC=O(4652)'],
    products = ['C6H7O(505)(504)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+11,'s^-1'), n=0.21, Ea=(125.52,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_DSMS;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R7_DSMS_CO;carbonylbond_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C6H7O(506)(505)'],
    products = ['C6H7O(505)(504)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(8.21347e+10,'s^-1'), n=0.43, Ea=(141.419,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_SMSM;multiplebond_intra;radadd_intra_cs2H] for rate rule [R7_SMSM_CO;carbonylbond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction54',
    reactants = ['C6H7O(505)(504)'],
    products = ['O[C]1C=CC=CC1(4654)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.51445e+09,'s^-1'), n=1.21848, Ea=(162.539,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_(CdCdCd)] for rate rule [R2H_S;O_rad_out;Cs_H_out_(CdCdCd)]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction55',
    reactants = ['C6H7O(505)(504)'],
    products = ['OC1[CH]C=CC=C1(4655)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.07519e+07,'s^-1'), n=1.60667, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_H/Cd] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction56',
    reactants = ['OC1[C]=CC=CC1(4656)'],
    products = ['C6H7O(505)(504)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction57',
    reactants = ['OC1C=CC=[C]C1(4657)'],
    products = ['C6H7O(505)(504)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction58',
    reactants = ['C6H7O(505)(504)'],
    products = ['OC1C=[C]C=CC1(4658)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.71035e+12,'s^-1'), n=1.11009, Ea=(419.408,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSD;Y_rad_out;Cd_H_out_single] for rate rule [R4H_SSD;O_rad_out;Cd_H_out_singleDe]
Euclidian distance = 1.73205080757
family: intra_H_migration"""),
)

reaction(
    label = 'reaction59',
    reactants = ['C6H7O(505)(504)'],
    products = ['OC1C=C[C]=CC1(4659)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.234e+06,'s^-1','*|/',3), n=1.554, Ea=(111.445,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R5H_SSSD;O_rad_out;Cd_H_out_single] for rate rule [R5H_SSSD;O_rad_out;Cd_H_out_singleDe]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction66',
    reactants = ['C6H7O(505)(504)'],
    products = ['[CH]1C2C=CCC1O2(4666)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(130.857,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_cyclic;doublebond_intra_pri_HCd;radadd_intra] for rate rule [Rn1c6_beta_short;doublebond_intra_pri_HCd;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction67',
    reactants = ['C6H7O(505)(504)'],
    products = ['[CH]1CC2C=CC1O2(4667)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.96035e+12,'s^-1'), n=-0.233729, Ea=(130.948,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Rn1c6_gamma;doublebond_intra_pri_HCd;radadd_intra] for rate rule [Rn1c6_gamma;doublebond_intra_pri_HCd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction68',
    reactants = ['C6H7O(505)(504)'],
    products = ['[CH]1C=CC2CC1O2(4668)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.42658e+11,'s^-1'), n=0.0609759, Ea=(92.8718,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_cyclic;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn1c6_beta_long;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction69',
    reactants = ['C6H7O(505)(504)'],
    products = ['[CH]1C=CCC2OC12(4669)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.48e+11,'s^-1'), n=0.26, Ea=(30.1248,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn1c6_alpha_long;doublebond_intra_pri;radadd_intra] for rate rule [Rn1c6_alpha_long;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction72',
    reactants = ['O(4)(4)', 'C6H7(464)(463)'],
    products = ['C6H7O(505)(504)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4171.1,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/OneDeC;O_birad]
Euclidian distance = 4.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction2',
    reactants = ['C6H7O(506)(505)'],
    products = ['O=C[CH]C1C=CC1(4610)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(196.02,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD_D;doublebond_intra;radadd_intra_cs] for rate rule [R5_SD_D;doublebond_intra_HDe_pri;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)(3)', 'C=C=CC=CC=O(4611)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.42e+08,'cm^3/(mol*s)'), n=1.64, Ea=(11.7989,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2713 used for Ca_Cds-HH;HJ
Exact match found for rate rule [Ca_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C6H7O(506)(505)'],
    products = ['C[C]=CC=CC=O(4612)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CC=[C]C=CC=O(4613)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(3.85113e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""From training reaction 288 used for R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CC=C[C]=CC=O(4614)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.33e+08,'s^-1'), n=1.1915, Ea=(103.605,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""From training reaction 284 used for R4H_SDS;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R4H_SDS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CC=CC=[C]C=O(4615)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5H_DSMS;Cd_rad_out_single;Cs_H_out_2H] for rate rule [R5H_DSMS;Cd_rad_out_singleDe;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CC=CC=C[C]=O(4616)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(25020,'s^-1'), n=2.00841, Ea=(94.1172,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6H;Y_rad_out;Cs_H_out_2H] + [R6H_SMSMS;Y_rad_out;Cs_H_out] for rate rule [R6H_SMSMS;CO_rad_out;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=C[CH2](924)', 'CHCHCHO(1679)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(7.23e+13,'cm^3/(mol*s)','+|-',1.2e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""Estimated using an average for rate rule [Cd_pri_rad;Cd_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)(3)', '[CH2]C=[C]C=CC=O(4617)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)(3)', '[CH2]C=C[C]=CC=O(4618)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [Cd_rad/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(3)(3)', '[CH2][C]=CC=CC=O(4619)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(4.34601e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction14',
    reactants = ['HCO(14)(15)', '[CH]=CC=C[CH2](922)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.81e+13,'cm^3/(mol*s)','*|/',3), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 91 used for Cd_pri_rad;CO_pri_rad
Exact match found for rate rule [CO_pri_rad;Cd_pri_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)(3)', '[CH2]C=CC=[C]C=O(4620)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad/OneDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)(3)', '[CH2]C=CC=C[C]=O(4621)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.53107e+07,'m^3/(mol*s)'), n=-0.133333, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;CO_sec_rad] for rate rule [H_rad;CO_rad/OneDe]
Euclidian distance = 1.41421356237
family: R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C6H7O(506)(505)'],
    products = ['O=CC=CC1[CH]C1(4622)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""From training reaction 142 used for R3_D;doublebond_intra_pri_HCd;radadd_intra_cs2H
Exact match found for rate rule [R3_D;doublebond_intra_pri_HCd;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C6H7O(506)(505)'],
    products = ['O=CC1[CH]C=CC1(4623)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.97937e+11,'s^-1'), n=-0.20996, Ea=(120.09,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R5_SD_D;doublebond_intra_pri_HDe;radadd_intra_cs2H] for rate rule [R5_SD_D;doublebond_intra_pri_HCO;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C6H7O(506)(505)'],
    products = ['C6H7O(511)(510)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(5.16404e+09,'s^-1'), n=0.622603, Ea=(78.7572,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7_SDSD_D;multiplebond_intra;radadd_intra_cs2H] + [R7_SDSD_D;carbonyl_intra_H;radadd_intra] for rate rule [R7_SDSD_D;carbonyl_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]CC=CC=C=O(4624)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(3.94565e+09,'s^-1'), n=0.909333, Ea=(116.555,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [1_3_unsaturated_pentane_backbone;CH(CJ)_1;unsaturated_end] + [1_3_pentadiene;CH_end;unsaturated_end] for rate rule [1_3_pentadiene;CH(CJ)_1;unsaturated_end]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction1',
    reactants = ['CO(10)(11)', 'C5H7(211)(210)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(0.236794,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO;R_H] for rate rule [CO;Cd_pri]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 4.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C=CC=C=CO(4625)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CH2(17)(18)', '[CH]=CC=CC=O(4626)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(2.23625e+06,'m^3/(mol*s)'), n=0.36814, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)(3)', '[CH]C=CC=CC=O(4627)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C6H7O(506)(505)'],
    products = ['C=CC1C=CC1[O](4628)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(221.26,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD;multiplebond_intra;radadd_intra_cs] for rate rule [R5_SD_CO;carbonylbond_intra_H;radadd_intra_csHCd]
Euclidian distance = 3.0
family: Intra_R_Add_Exocyclic
Ea raised from 219.7 to 221.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)(3)', 'C=CC=C=CC=O(4629)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(1.149e+09,'cm^3/(mol*s)'), n=1.595, Ea=(16.7946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-OneDeH;HJ] for rate rule [Ca_Cds-CdH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=CC[C]=CC=O(4630)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.448e+10,'s^-1'), n=0.82, Ea=(156.9,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 154 used for R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=[C]CC=CC=O(4631)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1.448e+10,'s^-1'), n=0.82, Ea=(156.9,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 154 used for R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=CCC=[C]C=O(4632)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleDe;Cs_H_out_H/Cd]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]=CCC=CC=O(4633)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=CCC=C[C]=O(4634)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(6.92799e+06,'s^-1'), n=1.7075, Ea=(114.432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R4H_SDS;CO_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C6H7O(506)(505)'],
    products = ['C=CC1[CH]C1C=O(4635)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(1.52971e+10,'s^-1'), n=0.596, Ea=(203.97,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;doublebond_intra_pri_HDe;radadd_intra_csHCd] + [R3_D;doublebond_intra_pri_HDe;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_pri_HCO;radadd_intra_csHCd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C6H7O(506)(505)'],
    products = ['C=CC1C=C[CH]O1(4636)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(8.8675e+11,'s^-1'), n=-0.37996, Ea=(137.453,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD;multiplebond_intra;radadd_intra_csHCd] for rate rule [R5_SD_CO;carbonyl_intra_H;radadd_intra_csHCd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C2H3(28)(29)', '[CH]C=CC=O(4637)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(2.23625e+06,'m^3/(mol*s)'), n=0.36814, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction34',
    reactants = ['C6H7O(506)(505)'],
    products = ['C=C[CH]C1C=CO1(4638)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(196.02,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD_D;doublebond_intra;radadd_intra] for rate rule [R5_SD_D;doublebond_intra_HCd_pri;radadd_intra_O]
Euclidian distance = 3.16227766017
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C6H7O(506)(505)'],
    products = ['[CH2]C1C=CC=CO1(4639)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(3.95e+10,'s^-1'), n=0.53, Ea=(131.796,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_SMSM_D;doublebond_intra_2H_pri;radadd_intra] for rate rule [R7_SMSM_D;doublebond_intra_2H_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction36',
    reactants = ['H(3)(3)', 'C6H6O(516)(515)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(1.54607e+06,'m^3/(mol*s)'), n=0.466452, Ea=(32.4943,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=CC=CC=[C]O(4640)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C=CC=[C]C=CO(4641)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(1.11e+08,'s^-1'), n=1.1915, Ea=(103.605,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R4H_SDS;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SDS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C=C[C]=CC=CO(4642)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_DSMS;Cd_rad_out_single;XH_out] for rate rule [R5H_DSMS;Cd_rad_out_singleDe_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C=[C]C=CC=CO(4643)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(280534,'s^-1'), n=2.03283, Ea=(131.895,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_Cd;XH_out] + [R6H_SMSMS;Y_rad_out;XH_out] for rate rule [R6H_SMSMS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 3.16227766017
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH]=CC=CC=CO(4644)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R7H;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['CHCHO(45)(45)', 'CH2CHCHCH(1116)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(7.23e+13,'cm^3/(mol*s)','+|-',1.2e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""Estimated using an average for rate rule [Cd_rad;Cd_pri_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction43',
    reactants = ['C6H7O(506)(505)'],
    products = ['C=CC=CC1[CH]O1(4645)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri_HCd;radadd_intra] for rate rule [R3_D;doublebond_intra_pri_HCd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction44',
    reactants = ['C=C=CC=CC[O](4646)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS65',
    kinetics = Arrhenius(A=(1.6907e+09,'s^-1'), n=1.02346, Ea=(163.761,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [1_3_unsaturated_pentane_backbone;CH_end;CddC_2] + [1_3_pentadiene;CH_end;unsaturated_end] for rate rule [1_3_pentadiene;CH_end;CddC_2]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction45',
    reactants = ['O(4)(4)', '[CH]=CC=CC=C(4647)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS66',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction46',
    reactants = ['C6H7O(506)(505)'],
    products = ['[CH2]C1C=CC1C=O(4648)'],
    transitionState = 'TS67',
    kinetics = Arrhenius(A=(2.17094e+10,'s^-1'), n=0.2405, Ea=(161.416,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5_SM;doublebond_intra_2H_pri;radadd_intra_cs] + [R5_SD_D;doublebond_intra;radadd_intra_cs] for rate rule [R5_SD_D;doublebond_intra_2H_pri;radadd_intra_csHDe]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction47',
    reactants = ['C=CC=[C]CC=O(4649)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS68',
    kinetics = Arrhenius(A=(1.82494e+10,'s^-1'), n=0.9, Ea=(134.725,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/OneDe] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/CO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['C=CC=CC[C]=O(3250)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS69',
    kinetics = Arrhenius(A=(2.9636e+08,'s^-1'), n=1.49631, Ea=(82.9758,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/(Cd-Cd-Cd)] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['C=C[C]=CCC=O(4650)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS70',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/CO]
Euclidian distance = 2.82842712475
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction50',
    reactants = ['C=[C]C=CCC=O(4651)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS71',
    kinetics = Arrhenius(A=(4.05862e+06,'s^-1'), n=1.71075, Ea=(73.9568,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SDS;Cd_rad_out_Cd;Cs_H_out] + [R4H;Cd_rad_out;Cs_H_out_H/OneDe] for rate rule [R4H_SDS;Cd_rad_out_Cd;Cs_H_out_H/CO]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH]=CC=CCC=O(4652)'],
    products = ['C6H7O(506)(505)'],
    transitionState = 'TS72',
    kinetics = Arrhenius(A=(272000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5H_DSMS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R5H_DSMS;Cd_rad_out_singleH;Cs_H_out_H/CO]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction73',
    reactants = ['H(3)(3)', 'C1=CC=COC=C1(4672)'],
    products = ['C6H7O(511)(510)'],
    transitionState = 'TS73',
    kinetics = Arrhenius(A=(4.364e+10,'cm^3/(mol*s)'), n=0.859, Ea=(6.76971,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [Cds-OsH_Cds;HJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction74',
    reactants = ['H(3)(3)', 'C1=CC=COCC=1(4673)'],
    products = ['C6H7O(511)(510)'],
    transitionState = 'TS74',
    kinetics = Arrhenius(A=(5.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(15.8155,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2714 used for Ca_Cds-CsH;HJ
Exact match found for rate rule [Ca_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction75',
    reactants = ['[CH]1CC=CC=CO1(4674)'],
    products = ['C6H7O(511)(510)'],
    transitionState = 'TS75',
    kinetics = Arrhenius(A=(9.74382e+07,'s^-1'), n=1.61783, Ea=(74.4341,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_1H;Cs_H_out_H/(Cd-Cd-Cd)] for rate rule [R2H_S;C_rad_out_H/NonDeO;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction76',
    reactants = ['[C]1=CC=COCC1(4675)'],
    products = ['C6H7O(511)(510)'],
    transitionState = 'TS76',
    kinetics = Arrhenius(A=(2.66329e+10,'s^-1'), n=0.993, Ea=(157.679,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction77',
    reactants = ['[C]1=CCCOC=C1(4676)'],
    products = ['C6H7O(511)(510)'],
    transitionState = 'TS77',
    kinetics = Arrhenius(A=(6.1583e+09,'s^-1'), n=0.92705, Ea=(170.178,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_DS;Cd_rad_out_single;Cs_H_out_H/NonDeC] + [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction78',
    reactants = ['[C]1=CC=CCCO1(4677)'],
    products = ['C6H7O(511)(510)'],
    transitionState = 'TS78',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out_H/(Cd-Cd-Cd)] for rate rule [R4H_SSS_OCs;Cd_rad_out_Cd;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 2.82842712475
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction79',
    reactants = ['[C]1=COCCC=C1(4678)'],
    products = ['C6H7O(511)(510)'],
    transitionState = 'TS79',
    kinetics = Arrhenius(A=(2.22e+08,'s^-1'), n=1.1915, Ea=(103.605,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R4H_SDS;Cd_rad_out_Cd;Cs_H_out] for rate rule [R4H_SDS;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction80',
    reactants = ['H(3)(3)', '[CH]1[CH]OC=CC=C1(4679)'],
    products = ['C6H7O(511)(510)'],
    transitionState = 'TS80',
    kinetics = Arrhenius(A=(8.69202e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction81',
    reactants = ['H(3)(3)', '[C]1[CH]COC=CC=1(4680)'],
    products = ['C6H7O(511)(510)'],
    transitionState = 'TS81',
    kinetics = Arrhenius(A=(4.34601e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction82',
    reactants = ['H(3)(3)', '[C]1=C[CH]COC=C1(4681)'],
    products = ['C6H7O(511)(510)'],
    transitionState = 'TS82',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction83',
    reactants = ['H(3)(3)', '[C]1=COC[CH]C=C1(4682)'],
    products = ['C6H7O(511)(510)'],
    transitionState = 'TS83',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [Cd_rad/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction84',
    reactants = ['H(3)(3)', '[C]1=CC=C[CH]CO1(4683)'],
    products = ['C6H7O(511)(510)'],
    transitionState = 'TS84',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [H_rad;Cd_rad/NonDe]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction85',
    reactants = ['C6H7O(511)(510)'],
    products = ['[CH]1C2C=COCC12(4684)'],
    transitionState = 'TS85',
    kinetics = Arrhenius(A=(2.22857e+12,'s^-1'), n=0, Ea=(182.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_beta_short;doublebond_intra_pri_HCd;radadd_intra_cs] for rate rule [Rn0c7_beta_short;doublebond_intra_pri_HCd;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction86',
    reactants = ['C6H7O(511)(510)'],
    products = ['[CH]1OCC2C=CC12(4685)'],
    transitionState = 'TS86',
    kinetics = Arrhenius(A=(5.65845e+18,'s^-1'), n=-1.57151, Ea=(227.049,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 131 used for Rn0c7_gamma_long_SDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCd
Exact match found for rate rule [Rn0c7_gamma_long_SDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCd]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction87',
    reactants = ['C6H7O(511)(510)'],
    products = ['[CH]1C=CC2COC12(4686)'],
    transitionState = 'TS87',
    kinetics = Arrhenius(A=(6.33004e+14,'s^-1'), n=-0.792922, Ea=(179.405,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Rn0c7_gamma_long_SSS_D;doublebond_intra_pri_HNd;radadd_intra_csHCs] for rate rule [Rn0c7_gamma_long_SSS_D;doublebond_intra_pri_HNd_O;radadd_intra_csHCs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction88',
    reactants = ['C6H7O(511)(510)'],
    products = ['O=CC1[CH]C=CC1(4623)'],
    transitionState = 'TS88',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_H;R_O_C]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction89',
    reactants = ['H(3)(3)', '[C]1C=CC=COC1(4687)'],
    products = ['C6H7O(511)(510)'],
    transitionState = 'TS89',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction90',
    reactants = ['H(3)(3)', 'C1=CC=CCOC=1(4688)'],
    products = ['C6H7O(511)(510)'],
    transitionState = 'TS90',
    kinetics = Arrhenius(A=(1.149e+09,'cm^3/(mol*s)'), n=1.595, Ea=(16.7946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-OneDeH;HJ] for rate rule [Ca_Cds-CdH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction91',
    reactants = ['[C]1=CCOC=CC1(4689)'],
    products = ['C6H7O(511)(510)'],
    transitionState = 'TS91',
    kinetics = Arrhenius(A=(1.448e+10,'s^-1'), n=0.82, Ea=(156.9,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 154 used for R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction92',
    reactants = ['[C]1=COCC=CC1(4690)'],
    products = ['C6H7O(511)(510)'],
    transitionState = 'TS92',
    kinetics = Arrhenius(A=(1.448e+10,'s^-1'), n=0.82, Ea=(156.9,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 154 used for R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction93',
    reactants = ['[C]1=CCC=COC1(4691)'],
    products = ['C6H7O(511)(510)'],
    transitionState = 'TS93',
    kinetics = Arrhenius(A=(1.182e+10,'s^-1'), n=0.86, Ea=(149.369,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction94',
    reactants = ['[C]1=CCC=CCO1(4692)'],
    products = ['C6H7O(511)(510)'],
    transitionState = 'TS94',
    kinetics = Arrhenius(A=(1.182e+10,'s^-1'), n=0.86, Ea=(149.369,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleNd;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleNd;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction95',
    reactants = ['[CH]1C=CCC=CO1(4693)'],
    products = ['C6H7O(511)(510)'],
    transitionState = 'TS95',
    kinetics = Arrhenius(A=(2.30755e+06,'s^-1'), n=1.805, Ea=(115.478,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;C_rad_out_single;Cs_H_out_H/Cd] for rate rule [R4H_SDS;C_rad_out_H/NonDeO;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction96',
    reactants = ['C6H7O(511)(510)'],
    products = ['[CH]1C2C=CCOC12(4694)'],
    transitionState = 'TS96',
    kinetics = Arrhenius(A=(4.00063e+13,'s^-1'), n=-0.283562, Ea=(245.321,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_beta;doublebond_intra_pri_HNd;radadd_intra_csHDe] for rate rule [Rn0c7_beta_short;doublebond_intra_pri_HNd_O;radadd_intra_csHCd]
Euclidian distance = 1.73205080757
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction97',
    reactants = ['[C]1=CC=CCOC1(4695)'],
    products = ['C6H7O(511)(510)'],
    transitionState = 'TS97',
    kinetics = Arrhenius(A=(1.89098e+10,'s^-1'), n=0.9884, Ea=(139.355,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_1H] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction98',
    reactants = ['[C]1=CCOCC=C1(4696)'],
    products = ['C6H7O(511)(510)'],
    transitionState = 'TS98',
    kinetics = Arrhenius(A=(2.56742e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/NonDeO]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction99',
    reactants = ['C6H7O(511)(510)'],
    products = ['[CH]1COC2C=CC12(4697)'],
    transitionState = 'TS99',
    kinetics = Arrhenius(A=(5.65845e+18,'s^-1'), n=-1.57151, Ea=(227.049,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 131 used for Rn0c7_gamma_long_SDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCd
Exact match found for rate rule [Rn0c7_gamma_long_SDS_D;doublebond_intra_pri_HCd;radadd_intra_csHCd]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

network(
    label = '429',
    isomers = [
        'C6H7O(505)(504)',
        'C6H7O(506)(505)',
        'C6H7O(511)(510)',
    ],
    reactants = [
        ('O(4)(4)', 'C6H7(464)(463)'),
        ('CO(10)(11)', 'C5H7(211)(210)'),
        ('H(3)(3)', 'C6H6O(516)(515)'),
    ],
    bathGas = {
        'Ne': 0.333333,
        'N2': 0.333333,
        'Ar(8)': 0.333333,
    },
)

pressureDependence(
    label = '429',
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

