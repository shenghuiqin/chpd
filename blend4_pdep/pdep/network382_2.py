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
    label = 'C6H7(467)(466)',
    structure = SMILES('[CH]1CC2C=CC12'),
    E0 = (296.996,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3501.78,'J/mol'), sigma=(6.02237,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=546.97 K, Pc=36.38 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.7024,0.0156163,5.63255e-05,-7.46196e-08,2.67302e-11,35778.4,16.0411], Tmin=(100,'K'), Tmax=(1013.01,'K')), NASAPolynomial(coeffs=[7.74818,0.0286307,-1.17185e-05,2.25842e-09,-1.63935e-13,33066.1,-16.7088], Tmin=(1013.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(296.996,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_4_ene_1) + radical(cyclobutane)"""),
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
    label = 'C6H6(468)(467)',
    structure = SMILES('c1ccccc1'),
    E0 = (68.5201,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3319.97,'J/mol'), sigma=(5.949e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.88775,0.00339716,0.0001009,-1.32851e-07,5.08462e-11,8300.31,11.4552], Tmin=(100,'K'), Tmax=(949.238,'K')), NASAPolynomial(coeffs=[12.521,0.0165025,-4.66524e-06,8.85737e-10,-7.161e-14,4052.17,-47.2613], Tmin=(949.238,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(68.5201,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + ring(Benzene)"""),
)

species(
    label = 'C1=CC=CCC=1(6493)',
    structure = SMILES('C1=CC=CCC=1'),
    E0 = (388.503,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.45852,0.0164141,6.44079e-05,-9.56319e-08,3.76994e-11,46797.5,14.9264], Tmin=(100,'K'), Tmax=(957.786,'K')), NASAPolynomial(coeffs=[12.9098,0.016872,-5.38368e-06,1.02574e-09,-7.96299e-14,42772.4,-45.6054], Tmin=(957.786,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(388.503,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + ring(124cyclohexatriene)"""),
)

species(
    label = '[CH]1[CH]C=CC=C1(1038)',
    structure = SMILES('[CH]1[CH]C=CC=C1'),
    E0 = (282.149,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.65505,0.030509,-1.75379e-06,-8.37309e-09,2.50596e-12,33982.1,11.3122], Tmin=(100,'K'), Tmax=(1572.61,'K')), NASAPolynomial(coeffs=[9.21399,0.0284297,-1.37001e-05,2.5963e-09,-1.76603e-13,30113.3,-29.0418], Tmin=(1572.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(282.149,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(Aromatic_pi_S_1_3) + radical(Aromatic_pi_S_1_3)"""),
)

species(
    label = '[C]1=CC=C[CH]C1(1463)',
    structure = SMILES('[C]1=CC=C[CH]C1'),
    E0 = (423.103,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.37414,0.0283839,1.08688e-05,-2.50799e-08,8.87725e-12,50952.1,15.7688], Tmin=(100,'K'), Tmax=(1142.54,'K')), NASAPolynomial(coeffs=[7.65087,0.028021,-1.24316e-05,2.38948e-09,-1.69054e-13,48564.2,-15.5646], Tmin=(1142.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(423.103,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(Cds_S) + radical(Aromatic_pi_S_1_3)"""),
)

species(
    label = '[C]1=CC[CH]C=C1(1599)',
    structure = SMILES('[C]1=CC[CH]C=C1'),
    E0 = (384.256,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.30338,0.0306596,5.26427e-06,-1.95796e-08,7.17224e-12,46281.9,15.1439], Tmin=(100,'K'), Tmax=(1144.3,'K')), NASAPolynomial(coeffs=[7.24701,0.028685,-1.22115e-05,2.29106e-09,-1.59741e-13,44148.4,-13.7513], Tmin=(1144.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(384.256,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(Aromatic_pi_S_1_3) + radical(C=CJC=C)"""),
)

species(
    label = '[C]1=C[CH]CC=C1(1600)',
    structure = SMILES('[C]1=C[CH]CC=C1'),
    E0 = (384.256,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.30338,0.0306596,5.26427e-06,-1.95796e-08,7.17224e-12,46281.9,14.4508], Tmin=(100,'K'), Tmax=(1144.3,'K')), NASAPolynomial(coeffs=[7.24701,0.028685,-1.22115e-05,2.29106e-09,-1.59741e-13,44148.4,-14.4445], Tmin=(1144.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(384.256,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(C=CJC=C) + radical(Aromatic_pi_S_1_3)"""),
)

species(
    label = '[C]1=CC=CCC1(1108)',
    structure = SMILES('[C]1=CC=CCC1'),
    E0 = (326.214,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.45174,0.0176354,6.31059e-05,-9.1148e-08,3.50495e-11,39305,15.9945], Tmin=(100,'K'), Tmax=(970.391,'K')), NASAPolynomial(coeffs=[11.3124,0.0223735,-7.99967e-06,1.52046e-09,-1.13384e-13,35642.6,-36.4967], Tmin=(970.391,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(326.214,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(Cds_S)"""),
)

species(
    label = '[C]1=CCCC=C1(6494)',
    structure = SMILES('[C]1=CCCC=C1'),
    E0 = (287.368,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.39937,0.0196701,5.8465e-05,-8.70866e-08,3.40475e-11,34634.1,15.3051], Tmin=(100,'K'), Tmax=(958.739,'K')), NASAPolynomial(coeffs=[10.9715,0.0229425,-7.72957e-06,1.41106e-09,-1.0321e-13,31196.3,-35.0456], Tmin=(958.739,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(287.368,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3-Cyclohexadiene) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]1C2C=CCC12(2084)',
    structure = SMILES('[CH]1C2C=CCC12'),
    E0 = (371.292,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.82805,0.008202,8.4483e-05,-1.08463e-07,3.97712e-11,44714,16.3084], Tmin=(100,'K'), Tmax=(982.373,'K')), NASAPolynomial(coeffs=[9.73472,0.0249485,-9.59857e-06,1.87676e-09,-1.40694e-13,41192,-27.9089], Tmin=(982.373,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(371.292,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_3_5_ene_1) + radical(cyclopropane)"""),
)

species(
    label = '[CH]1C=CC2CC12(6495)',
    structure = SMILES('[CH]1C=CC2CC12'),
    E0 = (272.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.86057,0.00241223,0.000111571,-1.43076e-07,5.35858e-11,32857.9,12.0439], Tmin=(100,'K'), Tmax=(962.926,'K')), NASAPolynomial(coeffs=[12.3626,0.0207511,-7.05059e-06,1.39673e-09,-1.10012e-13,28347.8,-47.3555], Tmin=(962.926,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(272.685,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_3_5_ene_1) + radical(cyclopentene-allyl)"""),
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
    label = '[CH2]C1C=CC=C1(6496)',
    structure = SMILES('[CH2]C1C=CC=C1'),
    E0 = (293.743,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.09212,0.0280859,3.66007e-05,-6.78951e-08,2.83607e-11,35410.4,16.6699], Tmin=(100,'K'), Tmax=(947.391,'K')), NASAPolynomial(coeffs=[12.1468,0.020431,-6.3738e-06,1.11466e-09,-8.03406e-14,31943.7,-39.5363], Tmin=(947.391,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293.743,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(Cyclopentadiene) + radical(Isobutyl)"""),
)

species(
    label = '[C]1=CCC=CC1(6497)',
    structure = SMILES('[C]1=CCC=CC1'),
    E0 = (331.659,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.47655,0.0176703,6.15711e-05,-8.80657e-08,3.35374e-11,39958.4,18.5952], Tmin=(100,'K'), Tmax=(976.335,'K')), NASAPolynomial(coeffs=[10.7243,0.0236542,-8.73051e-06,1.6643e-09,-1.2325e-13,36452.2,-30.7065], Tmin=(976.335,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(331.659,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cyclohexadiene) + radical(Cds_S)"""),
)

species(
    label = 'C1=CC2CC=C12(6504)',
    structure = SMILES('C1=CC2CC=C12'),
    E0 = (368.085,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13799,0.0288071,2.38854e-05,-5.02784e-08,2.07786e-11,44348.3,10.5553], Tmin=(100,'K'), Tmax=(983.644,'K')), NASAPolynomial(coeffs=[11.6017,0.0198168,-7.38153e-06,1.39574e-09,-1.02202e-13,41059.7,-42.1985], Tmin=(983.644,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(368.085,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s2_4_4_ane) - ring(Cyclobutane) - ring(Cyclobutane) + ring(Cyclobutene) + ring(Cyclobutene)"""),
)

species(
    label = 'C1=CC2C=CC12(1049)',
    structure = SMILES('C1=CC2C=CC12'),
    E0 = (389.58,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.24105,0.0324057,-1.67402e-06,-1.3176e-08,5.20455e-12,46924.1,9.01711], Tmin=(100,'K'), Tmax=(1174.09,'K')), NASAPolynomial(coeffs=[7.9921,0.0263739,-1.12936e-05,2.12407e-09,-1.4809e-13,44638.9,-23.632], Tmin=(1174.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(389.58,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s2_4_4_ane) - ring(Cyclobutane) - ring(Cyclobutane) + ring(Cyclobutene) + ring(Cyclobutene)"""),
)

species(
    label = 'C1=CC2CC[C]12(6505)',
    structure = SMILES('C1=CC2CC[C]12'),
    E0 = (241.139,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.81363,0.00967993,7.97589e-05,-1.03181e-07,3.80339e-11,29059.7,12.2392], Tmin=(100,'K'), Tmax=(975.337,'K')), NASAPolynomial(coeffs=[8.99886,0.0260691,-9.66384e-06,1.83543e-09,-1.35276e-13,25867.1,-27.628], Tmin=(975.337,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(241.139,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_4_ene_1) + radical(Allyl_T)"""),
)

species(
    label = '[C]1=CC2CCC12(6506)',
    structure = SMILES('[C]1=CC2CCC12'),
    E0 = (361.644,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.50087,0.0196907,4.94992e-05,-7.11926e-08,2.64536e-11,43561.3,14.8898], Tmin=(100,'K'), Tmax=(1000.14,'K')), NASAPolynomial(coeffs=[9.11255,0.0264691,-1.04924e-05,2.00873e-09,-1.46169e-13,40577.3,-25.3152], Tmin=(1000.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.644,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_4_ene_1) + radical(cyclobutene-vinyl)"""),
)

species(
    label = '[CH]1C[C]2C=CC12(6507)',
    structure = SMILES('[CH]1C[C]2C=CC12'),
    E0 = (428.977,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.98203,0.010204,6.18138e-05,-7.73758e-08,2.73665e-11,51641.4,13.9834], Tmin=(100,'K'), Tmax=(1004.72,'K')), NASAPolynomial(coeffs=[6.47902,0.027674,-1.11354e-05,2.12659e-09,-1.5375e-13,49354.2,-10.7891], Tmin=(1004.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(428.977,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_4_ene_1) + radical(cyclobutane) + radical(Allyl_T)"""),
)

species(
    label = '[CH]1CC2C=C[C]12(6508)',
    structure = SMILES('[CH]1CC2C=C[C]12'),
    E0 = (428.977,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.98196,0.0102048,6.1811e-05,-7.73721e-08,2.73649e-11,51641.4,13.9836], Tmin=(100,'K'), Tmax=(1004.73,'K')), NASAPolynomial(coeffs=[6.47929,0.0276735,-1.11351e-05,2.12653e-09,-1.53745e-13,49354.1,-10.7907], Tmin=(1004.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(428.977,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_4_ene_1) + radical(Allyl_T) + radical(cyclobutane)"""),
)

species(
    label = '[CH]1[CH]C2C=CC12(1057)',
    structure = SMILES('[CH]1[CH]C2C=CC12'),
    E0 = (484.835,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.87312,0.0160689,3.88556e-05,-4.97889e-08,1.66487e-11,58360,17.0865], Tmin=(100,'K'), Tmax=(1077.19,'K')), NASAPolynomial(coeffs=[5.48658,0.0298224,-1.29623e-05,2.49756e-09,-1.78204e-13,56436,-2.0335], Tmin=(1077.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(484.835,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_4_ene_1) + radical(cyclobutane) + radical(cyclobutane)"""),
)

species(
    label = '[C]1=CC2[CH]CC12(6509)',
    structure = SMILES('[C]1=CC2[CH]CC12'),
    E0 = (549.482,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66614,0.0202224,3.16793e-05,-4.57955e-08,1.60793e-11,66143.2,16.6469], Tmin=(100,'K'), Tmax=(1057.69,'K')), NASAPolynomial(coeffs=[6.7792,0.027773,-1.17969e-05,2.26153e-09,-1.6153e-13,63980.7,-9.53653], Tmin=(1057.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(549.482,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_4_ene_1) + radical(cyclobutane) + radical(cyclobutene-vinyl)"""),
)

species(
    label = '[C]1=CC2C[CH]C12(6510)',
    structure = SMILES('[C]1=CC2C[CH]C12'),
    E0 = (549.482,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66614,0.0202224,3.16793e-05,-4.57955e-08,1.60793e-11,66143.2,16.6469], Tmin=(100,'K'), Tmax=(1057.69,'K')), NASAPolynomial(coeffs=[6.7792,0.027773,-1.17969e-05,2.26153e-09,-1.6153e-13,63980.7,-9.53653], Tmin=(1057.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(549.482,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_4_ene_1) + radical(cyclobutane) + radical(cyclobutene-vinyl)"""),
)

species(
    label = '[CH]1C2C3CC2C13(6511)',
    structure = SMILES('[CH]1C2C3CC2C13'),
    E0 = (582.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.13339,0.00989783,5.74191e-05,-6.40571e-08,1.99976e-11,70048.2,7.54722], Tmin=(100,'K'), Tmax=(1117.69,'K')), NASAPolynomial(coeffs=[3.82914,0.0362815,-1.6739e-05,3.28891e-09,-2.35834e-13,68089.2,-3.95471], Tmin=(1117.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(582.092,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + polycyclic(s2_4_4_ane) + polycyclic(s2_4_4_ane) + polycyclic(s3_4_4_ane) - ring(Cyclobutane) - ring(Cyclobutane) - ring(Cyclobutane) + radical(bicyclo[1.1.1]pentane-C2)"""),
)

species(
    label = '[CH]1C2CC3C1C23(6512)',
    structure = SMILES('[CH]1C2CC3C1C23'),
    E0 = (512.512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.97749,0.0146544,4.47007e-05,-5.17296e-08,1.59488e-11,61684.1,10.6518], Tmin=(100,'K'), Tmax=(1149.15,'K')), NASAPolynomial(coeffs=[4.19744,0.0360225,-1.66265e-05,3.24585e-09,-2.31096e-13,59712.4,-2.76203], Tmin=(1149.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(512.512,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + polycyclic(s2_4_4_ane) + polycyclic(s2_3_4_ane) + polycyclic(s2_3_4_ane) - ring(Cyclobutane) - ring(Cyclopropane) - ring(Cyclobutane) + radical(bicyclo[2.1.1]hexane-C5)"""),
)

species(
    label = 'C=CC1[CH]C=C1(6513)',
    structure = SMILES('C=CC1[CH]C=C1'),
    E0 = (380.712,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26049,0.0233279,4.76149e-05,-7.54783e-08,2.95513e-11,45865.1,16.0766], Tmin=(100,'K'), Tmax=(976.926,'K')), NASAPolynomial(coeffs=[11.4809,0.0226746,-8.34594e-06,1.58296e-09,-1.1677e-13,42293.2,-37.2524], Tmin=(976.926,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(380.712,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(cyclobutene-allyl)"""),
)

species(
    label = '[CH]=CC1C=CC1(6514)',
    structure = SMILES('[CH]=CC1C=CC1'),
    E0 = (463.266,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,2950,3050,3150,900,950,1000,1050,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.97991,0.0321767,2.16535e-05,-4.9167e-08,2.04817e-11,55801.6,18.5934], Tmin=(100,'K'), Tmax=(985.355,'K')), NASAPolynomial(coeffs=[11.6965,0.0222964,-8.31043e-06,1.5549e-09,-1.12621e-13,52451.5,-35.423], Tmin=(985.355,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(463.266,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Cds_P)"""),
)

species(
    label = '[CH]1C2C=CC1C2(2154)',
    structure = SMILES('[CH]1C2C=CC1C2'),
    E0 = (459.155,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.99154,0.00110197,0.000108614,-1.3598e-07,5.00938e-11,55278.9,14.4118], Tmin=(100,'K'), Tmax=(970.474,'K')), NASAPolynomial(coeffs=[10.8477,0.0227184,-8.2563e-06,1.63668e-09,-1.26363e-13,51211.3,-36.3547], Tmin=(970.474,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(459.155,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s3_4_5_ene_1) + radical(bicyclo[2.1.1]hex-2-ene-C5)"""),
)

species(
    label = '[CH2]C1C2C=CC12(6515)',
    structure = SMILES('[CH2]C1C2C=CC12'),
    E0 = (490.239,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.39832,0.0163377,7.36619e-05,-1.08412e-07,4.31163e-11,59037,17.656], Tmin=(100,'K'), Tmax=(944.053,'K')), NASAPolynomial(coeffs=[13.3451,0.0178855,-4.95269e-06,8.8275e-10,-6.81264e-14,54834.3,-45.8377], Tmin=(944.053,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(490.239,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_3_4_ene_1) + radical(Isobutyl)"""),
)

species(
    label = '[C]1CC2C=CC12(6516)',
    structure = SMILES('[C]1CC2C=CC12'),
    E0 = (557.373,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.48659,0.0180988,5.4871e-05,-8.08043e-08,3.09972e-11,67104.5,14.372], Tmin=(100,'K'), Tmax=(979.852,'K')), NASAPolynomial(coeffs=[11.0588,0.0211227,-7.95739e-06,1.53993e-09,-1.15072e-13,63599.6,-36.1223], Tmin=(979.852,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_4_ene_1) + radical(CCJ2_triplet)"""),
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
    E0 = (299.554,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (616.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (493.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (634.895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (600.245,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (596.506,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (483.894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (457.546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (371.292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (316.118,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (339.219,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (400.839,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (453.678,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (488.559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (591.007,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (602.684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (424.915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (503.063,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (641.288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (640.769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (696.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (761.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (761.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (582.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (512.512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (505.884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (464.521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (456.931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (619.09,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (650.174,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (769.165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['H(3)(3)', 'C6H6(468)(467)'],
    products = ['C6H7(464)(463)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(9.221e+08,'cm^3/(mol*s)'), n=1.608, Ea=(4.599,'kcal/mol'), T0=(1,'K'), comment="""Matched reaction 67 C6H6-2 + H <=> C6H7-8 in R_Addition_MultipleBond/training
This reaction matched rate rule [Cds-CdH_Cds-CdH;HJ]
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(3)(3)', 'C1=CC=CCC=1(6493)'],
    products = ['C6H7(464)(463)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(15.8155,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2714 used for Ca_Cds-CsH;HJ
Exact match found for rate rule [Ca_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)(3)', '[CH]1[CH]C=CC=C1(1038)'],
    products = ['C6H7(464)(463)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.60761e+07,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 6.0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)(3)', '[C]1=CC=C[CH]C1(1463)'],
    products = ['C6H7(464)(463)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)(3)', '[C]1=CC[CH]C=C1(1599)'],
    products = ['C6H7(464)(463)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [Cd_rad/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)(3)', '[C]1=C[CH]CC=C1(1600)'],
    products = ['C6H7(464)(463)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[C]1=CC=CCC1(1108)'],
    products = ['C6H7(464)(463)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.66329e+10,'s^-1'), n=0.993, Ea=(157.679,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[C]1=CCCC=C1(6494)'],
    products = ['C6H7(464)(463)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6.1583e+09,'s^-1'), n=0.92705, Ea=(170.178,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_DS;Cd_rad_out_single;Cs_H_out_H/NonDeC] + [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C6H7(464)(463)'],
    products = ['[CH]1C2C=CCC12(2084)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.45714e+12,'s^-1'), n=0, Ea=(186.032,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_beta;doublebond_intra_pri_HCd;radadd_intra_cs] for rate rule [Rn0c6_beta_short;doublebond_intra_pri_HCd;radadd_intra_csHCs]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 184.0 to 186.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['C6H7(464)(463)'],
    products = ['C6H7(467)(466)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.08454e+19,'s^-1'), n=-0.859165, Ea=(130.857,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_gamma;doublebond_intra_pri_HCd;radadd_intra_cs] for rate rule [Rn0c6_gamma;doublebond_intra_pri_HCd;radadd_intra_csHCd]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C6H7(464)(463)'],
    products = ['[CH]1C=CC2CC12(6495)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.307e+12,'s^-1'), n=0.256, Ea=(36.797,'kcal/mol'), T0=(1,'K'), comment="""Matched reaction 83 C6H7-10 <=> C6H7-9 in Intra_R_Add_Endocyclic/training
This reaction matched rate rule [Rn0c6_beta_long_SS_D_HH;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=CC=CC=C(4647)'],
    products = ['C6H7(464)(463)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.57e+10,'s^-1'), n=0.43, Ea=(8.05002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_DSM_D;doublebond_intra_pri_2H;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C1C=CC=C1(6496)'],
    products = ['C6H7(464)(463)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.31121e+11,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[C]1=CCC=CC1(6497)'],
    products = ['C6H7(464)(463)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.448e+10,'s^-1'), n=0.82, Ea=(156.9,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 154 used for R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)(3)', 'C1=CC2CC=C12(6504)'],
    products = ['C6H7(467)(466)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.11e+08,'cm^3/(mol*s)'), n=1.64, Ea=(11.1294,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2594 used for Cds-CdCs_Cds-CsH;HJ
Exact match found for rate rule [Cds-CdCs_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)(3)', 'C1=CC2C=CC12(1049)'],
    products = ['C6H7(467)(466)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1042.23,'m^3/(mol*s)'), n=1.29683, Ea=(1.31187,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 104 used for Cds-CsH_Cds-(CsH-Cds-Cds)_cyc6;HJ
Exact match found for rate rule [Cds-CsH_Cds-(CsH-Cds-Cds)_cyc6;HJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C6H7(467)(466)'],
    products = ['C1=CC2CC[C]12(6505)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.15233e+09,'s^-1'), n=1.014, Ea=(127.919,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_noH] for rate rule [R2H_S_cy4;C_rad_out_H/NonDeC;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[C]1=CC2CCC12(6506)'],
    products = ['C6H7(467)(466)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.32e+09,'s^-1'), n=0.99, Ea=(141.419,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_Cd;Cs_H_out_H/NonDeC] for rate rule [R3H_SS_12cy4;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(3)(3)', '[CH]1C[C]2C=CC12(6507)'],
    products = ['C6H7(467)(466)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.92e+13,'cm^3/(mol*s)'), n=0.18, Ea=(0.518816,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 123 used for H_rad;C_rad/OneDeCs2
Exact match found for rate rule [C_rad/OneDeCs2;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(3)(3)', '[CH]1CC2C=C[C]12(6508)'],
    products = ['C6H7(467)(466)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.34601e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(3)(3)', '[CH]1[CH]C2C=CC12(1057)'],
    products = ['C6H7(467)(466)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(8.69202e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)(3)', '[C]1=CC2[CH]CC12(6509)'],
    products = ['C6H7(467)(466)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [H_rad;Cd_rad/NonDe]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)(3)', '[C]1=CC2C[CH]C12(6510)'],
    products = ['C6H7(467)(466)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C6H7(467)(466)'],
    products = ['[CH]1C2C3CC2C13(6511)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(285.096,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_gamma;doublebond_intra_pri;radadd_intra_csHCs] for rate rule [Rn0c6_gamma;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic
Ea raised from 284.1 to 285.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction25',
    reactants = ['C6H7(467)(466)'],
    products = ['[CH]1C2CC3C1C23(6512)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.03212e+12,'s^-1'), n=0.0431209, Ea=(215.516,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_cyclic;doublebond_intra_pri;radadd_intra_csHCs] for rate rule [Rn1c4_alpha_long;doublebond_intra_pri;radadd_intra_csHCs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 215.2 to 215.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=CC1[CH]C=C1(6513)'],
    products = ['C6H7(467)(466)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(7.45205e+07,'s^-1'), n=1.21458, Ea=(125.171,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs] for rate rule [R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_csHCd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]=CC1C=CC1(6514)'],
    products = ['C6H7(467)(466)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.1e+12,'s^-1'), n=0.14, Ea=(1.2552,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_cyclic;doublebond_intra_pri;radadd_intra_cdsingle] for rate rule [Rn2c4_alpha_long;doublebond_intra_pri;radadd_intra_cdsingleH]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C6H7(467)(466)'],
    products = ['[CH]1C=CC2CC12(6495)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]1C2C=CC1C2(2154)'],
    products = ['C6H7(467)(466)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.31121e+11,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C1C2C=CC12(6515)'],
    products = ['C6H7(467)(466)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.31121e+11,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction31',
    reactants = ['H(3)(3)', '[C]1CC2C=CC12(6516)'],
    products = ['C6H7(467)(466)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '382',
    isomers = [
        'C6H7(464)(463)',
        'C6H7(467)(466)',
    ],
    reactants = [
        ('H(3)(3)', 'C6H6(468)(467)'),
    ],
    bathGas = {
        'Ne': 0.333333,
        'N2': 0.333333,
        'Ar(8)': 0.333333,
    },
)

pressureDependence(
    label = '382',
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

